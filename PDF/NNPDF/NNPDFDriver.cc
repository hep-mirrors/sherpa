/***
 *
 * NNPDF C++ Driver
 *
 * Stefano Carrazza for the NNPDF Collaboration
 * email: stefano.carrazza@mi.infn.it
 *
 * October 2014
 *
 * Usage:
 *
 *  NNPDFDriver *pdf = new NNPDFDriver("gridname.LHgrid");
 *
 *  pdf->initPDF(0); // select replica [0,fMem]
 *
 *  or 
 * 
 *  NNPDFDriver *pdf = new NNPDFDriver("gridname.LHgrid", 0);
 *
 *  then
 *
 *  pdf->xfx(x,Q,fl); // -> returns double
 *
 *  // with fl = [-6,7], LHAPDF format
 *
 */

#include "NNPDFDriver.h"
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cmath>

#define XMINGRID 1e-9
#define NNDriverVersion "1.0.5"

using namespace std;

void split(vector<string>& results, string const& input)
{
  stringstream strstr(input);
  std::istream_iterator<string> it(strstr);
  std::istream_iterator<string> end;
  results.assign(it, end);
  return;
}

/**
 * @brief Default constructor
 * @param gridfilename the string containing the LHgrid file name and location.
 */
NNPDFDriver::NNPDFDriver(string const& gridfilename, int const& rep):
  fNFL(13),
  fNX(100),
  fMem(0),
  fRep(0),
  fAlphas(0),
  fXGrid(NULL),
  fLogXGrid(NULL),
  fHasPhoton(false),
  fSingleMem(true),
  //HS
  fxmin(0),
  fxmax(0),
  fQmin(0),
  fQmax(0),
  fMz(0),  
  forder(0),
  fNFlavors(5),
  // Masses
  fMDown(0),
  fMUp(0),  
  fMStrange(0),
  fMCharm(0),
  fMBottom(0), 
  fMTop(0)
  
{
  // Logo --- off
  //cout << " ****************************************" << endl;
  //cout << "      NNPDFDriver version " << NNDriverVersion << endl;
  //cout << "  Grid: " << gridfilename << endl;
  //cout << " ****************************************" << endl;


  // Read PDFs from file
  readPDFSet(gridfilename,rep);
}

/**
 * @brief NNPDFDriver::~NNPDFDriver the destructor
 */
NNPDFDriver::~NNPDFDriver()
{  
  for (size_t s = 0; s < fPDFGrid.size(); s++)
    for (int imem = 0; imem <= fMem; imem++)
      {
	for (int i = 0; i < fNFL; i++)
	  {
	    for (int j = 0; j < fNX; j++)
	      if (fPDFGrid[s][imem][i][j]) delete[] fPDFGrid[s][imem][i][j];
	    if (fPDFGrid[s][imem][i]) delete[] fPDFGrid[s][imem][i];
	  }
	if (fPDFGrid[s][imem]) delete[] fPDFGrid[s][imem];
      }
  fPDFGrid.clear();

  if (fXGrid) delete[] fXGrid;
  if (fLogXGrid) delete[] fLogXGrid;

  for (size_t t = 0; t < fQ2Grid.size(); t++)
    if (fQ2Grid[t]) delete[] fQ2Grid[t];
  fQ2Grid.clear();

  for (size_t t = 0; t < fLogQ2Grid.size(); t++)
    if (fLogQ2Grid[t]) delete[] fLogQ2Grid[t];  
  fLogQ2Grid.clear();

  fNQ2.clear();
}

/**
 * @brief NNPDFDriver::initPDF
 * @param irep
 */
void NNPDFDriver::initPDF(int irep)
{
  if (fSingleMem)
    {
      cout << "Error: initPDF not available due to the constructor" << endl;
      exit(-1);
    }
  else
    {      
      if (irep > fMem || irep < 0)
	{
	  cout << "Error: replica out of range [0," << fMem << "]" << endl;
	  exit(-1);
	}
      else
	fRep = irep;
    }
}

inline bool hasKey(string const &t, string const &key) {
  if (t.find(key) != string::npos) return true;
  else return false;
}

inline double readDouble(string const &t, string const &key) {
  vector<string> splitstring;
  split(splitstring,t);
  return atof(splitstring[1].c_str());
}

inline int readInt(string const &t, string const &key) {
  vector<string> splitstring;
  split(splitstring,t);
  return atoi(splitstring[1].c_str());
}

/**
 * @brief NNPDFDriver::readPDFSet read the LHgrid file into an array
 * @param grid the LHgrid filename
 */
void NNPDFDriver::readPDFSet(string const& grid, int const& rep)
{
  fstream f;
  stringstream file("");
  int firstindex = (int) grid.find_last_of("/") + 1;
  int lastindex  = (int) grid.length() - firstindex;
  
  string name = grid.substr(firstindex, lastindex);      
  file << grid << "/" << name << ".info";
  f.open(file.str().c_str(), ios::in);
  
  if (f.fail())
    {
      cout << "In NNPDFDriver::readPDFSet fLHAPDF6 --- Error: cannot open file " << file.str() << endl;
      exit(-1);
    }

  // Read in Meta data
  string tmp;
  vector<string> splitstring;	      
  for (;;)
    {
      getline(f,tmp);

      // Be a bit quieter
      //if (tmp.find("SetDesc:") != string::npos)
        //cout << tmp << endl;
      
      if (hasKey(tmp, "NumMembers:")) fMem = readInt(tmp, "NumMembers:") -1;
      if (tmp.find("Flavors: [") != string::npos) 	  
        if (tmp.find("22") != string::npos) { fHasPhoton = true;  fNFL++; }
     
      

      if (hasKey(tmp, "AlphaS_MZ:")) fAlphas = readDouble(tmp, "AlphaS_MZ:");

      if (hasKey(tmp, "XMin:")) {
        fxmin = readDouble(tmp, "XMin:");
        if ( fabs(fxmin - XMINGRID) > 1e-10) {
          cout << "Problem with XMINGRID" << endl;
          exit(-1);
        }
      }
      
      if (hasKey(tmp, "NumFlavors:")) fNFlavors = readInt(tmp, "NumFlavors:");


      if (hasKey(tmp, "XMax:")) fxmax = readDouble(tmp, "XMax:");
      if (hasKey(tmp, "QMin:")) fQmin = readDouble(tmp, "QMin:");
      if (hasKey(tmp, "QMax:")) fQmax = readDouble(tmp, "QMax:");
      if (hasKey(tmp, "MZ:") and !hasKey(tmp, "AlphaS_MZ:")  ) { // Not that great workaround for apparently
        fMz = readDouble(tmp, "MZ:");                            // slightly imperfect hadKey method
                                                                 // needs imrpovement in the future
      }
      if (hasKey(tmp, "AlphaS_OrderQCD:")) forder = readInt(tmp, "AlphaS_OrderQCD:");

      if (hasKey(tmp, "MDown"))    fMDown    = readDouble(tmp, "MDown");
      if (hasKey(tmp, "MUp"))      fMUp      = readDouble(tmp, "MUp");
      if (hasKey(tmp, "MStrange")) fMStrange = readDouble(tmp, "MStrange");
      if (hasKey(tmp, "MCharm"))   fMCharm   = readDouble(tmp, "MCharm");
      if (hasKey(tmp, "MBottom"))  fMBottom  = readDouble(tmp, "MBottom");
      if (hasKey(tmp, "MTop"))     fMTop     = readDouble(tmp, "MTop");

      
      if (f.eof()) break;
    }
  f.close();

  //// single member switcher
  fMem=0;

  // Read in PDF stuff
  stringstream gfile("");
  if (rep < 10)
    gfile << grid << "/" << name << "_000" << rep << ".dat";
  else if (rep < 100)
    gfile << grid << "/" << name << "_00" << rep << ".dat";
  else if (rep < 1000)
    gfile << grid << "/" << name << "_0" << rep << ".dat";
  else
    gfile << grid << "/" << name << "_" << rep << ".dat";

  // Some status/debug message
  // cout <<"NNPDFDriver 1.0.5, opening file " << gfile.str().c_str() << " to read PDF data " << endl;
  // cout << "The NNPDF collaboration:   arXiv:1410.8849 [hep-ph]" << endl;
  fstream g;
  g.open(gfile.str().c_str(), ios::in);
  
  getline(g, tmp);
  getline(g, tmp);
  getline(g, tmp);

  int sub = 0;
  for (;;)
    {	
      if (sub == 0)
        {
          // reading Xgrid
          getline(g, tmp);
          split(splitstring,tmp);
      
          fNX = splitstring.size();
          fXGrid = new double[fNX];
          fLogXGrid = new double[fNX];
          for (int ix = 0; ix < fNX; ix++)
            {
              fXGrid[ix] = atof(splitstring[ix].c_str());
              fLogXGrid[ix] = log(fXGrid[ix]);
            }
        }

      getline(g, tmp);
      split(splitstring, tmp);
      
      fNQ2.push_back(splitstring.size());	    
      fQ2Grid.push_back(new double[fNQ2[sub]]);
      fLogQ2Grid.push_back(new double[fNQ2[sub]]);

      for (int iq = 0; iq < fNQ2[sub]; iq++)
        {
          fQ2Grid[sub][iq] = pow(atof(splitstring[iq].c_str()), 2.0);
          fLogQ2Grid[sub][iq] = log(fQ2Grid[sub][iq]);	       
        }

      // skip flavor line
      getline(g, tmp);

      // This is nasty --- a vector is created that keeps track of the
      // order in which parton flavours are stored in the file
      vector<int> fls;
      
      split(splitstring,tmp);
      for (int i = 0; i < splitstring.size(); i++)
        {	      
          if (atoi(splitstring[i].c_str()) == 21)
            fls.push_back(6);
          else
            fls.push_back(atoi(splitstring[i].c_str())+6);
        }

      // building PDFgrid
      fPDFGrid.push_back(new double***[fMem+1]);     
      for (int imem = 0; imem <= fMem; imem++)
        {
          fPDFGrid[sub][imem] = new double**[fNFL];
          for (int i = 0; i < fNFL; i++)
            {
              fPDFGrid[sub][imem][i] = new double*[fNX];
              for (int j = 0; j < fNX; j++)
                {
                  fPDFGrid[sub][imem][i][j] = new double[fNQ2[sub]];
                  for (int z = 0; z < fNQ2[sub]; z++)
                    fPDFGrid[sub][imem][i][j][z] = 0.0;
                }
            }
        }
      
      // read PDF grid points	      
      for (int imem = 0; imem <= fMem; imem++)
        for (int ix = 0; ix < fNX; ix++)
          for (int iq = 0; iq < fNQ2[sub]; iq++)
            for (int fl = 0; fl < fls.size(); fl++) {
              // The filling of the grid here is done in the same order as they appear in the file
              // fls[fl] has values in [1, 11] the last iteration fills the gluon
              g >> fPDFGrid[sub][imem][fls[fl]][ix][iq];
            }

      getline(g, tmp);
      getline(g, tmp);

      if (tmp.find("---") != string::npos)
        {
          sub++;
          getline(g, tmp);
          if (g.eof()) break;
          continue;
        }
    }
  
  g.close();

}

/**
 * @brief NNPDFDriver::xfx
 * @param x
 * @param Q
 * @param id
 * @return
 */
double NNPDFDriver::xfx(double const&X, double const& Q2_glob, int const& ID)
{
  double res= 0; // Initialise return value
  double Q2 = Q2_glob;
  double x  = X;

  int id    = ID; // 0 is anti top, 6 is gluon, 12 is top
  //int id    = ID+6;
  int sub   = 0;

  for (size_t i = 0; i < fNQ2.size(); i++)
      if ( Q2 >= fQ2Grid[i][0]) sub = i;

  // check bounds
  if (x < XMINGRID || x < fXGrid[0] || x > fXGrid[fNX-1]) {
    cout << "Parton interpolation: x out of range -- freezed" << endl;  
    if (x < fXGrid[0])  x = fXGrid[0];
    if (x < XMINGRID)   x = XMINGRID;
     if (x > fXGrid[fNX-1]) x = fXGrid[fNX-1];
  }
  if (Q2 < fQ2Grid[sub][0] || Q2 > fQ2Grid[sub][fNQ2[sub]-1]) {
    cout << "Parton interpolation: Q2 out of range -- freezed" << endl;
    cout << Q2 << "\t" << fQ2Grid[sub][0] << "\t" << sub << endl;
    if (Q2 < fQ2Grid[sub][0]) Q2 = fQ2Grid[sub][0];
    if (Q2 > fQ2Grid[sub][fNQ2[sub]-1]) Q2 = fQ2Grid[sub][fNQ2[sub]-1];
  }

  // find nearest points in the grid
  int minx = 0;
  int maxx = fNX;

  while (maxx-minx > 1)
    {
      int midx = (minx+maxx)/2;
      if (x < fXGrid[midx]) 
	maxx = midx;
      else
	minx = midx;      
    }
  int ix = minx;

  int minq = 0;
  int maxq = fNQ2[sub];
  while (maxq-minq > 1)
    {
      int midq = (minq+maxq)/2;
      if (Q2 < fQ2Grid[sub][midq])
	maxq = midq;
      else
	minq = midq;      
    }
  int iq2 = minq;

  // Assign grid for interpolation. M,N -> order of polyN interpolation
  int   ix1a[fM], ix2a[fN];
  double x1a[fM], x2a[fN];
  double ya[fM][fN];

  for (int i = 0; i < fM; i++)
    {
      // one of the following if statements must overwrite,
      // this initialization, otherwise we will exit below
      ix1a[i] = -1; 
      if (ix+1 >= fM/2 && ix+1 <= (fNX-fM/2))
        ix1a[i] = ix+1 - fM/2 + i;
      if (ix+1 < fM/2)
        ix1a[i] = i;
      if (ix+1 > (fNX-fM/2))
        ix1a[i] = (fNX-fM) + i;
	
      // Check grids
      if (ix1a[i] < 0 || ix1a[i] >= fNX) {
        cout << "Error in grids! i, ixia[i] = " 
             << i << "\t" << ix1a[i] << endl;
        exit(-1);
    	}
    }

  for (int j = 0; j < fN; j++)
    {
      if (iq2+1 >= fN/2 && iq2+1 <= (fNQ2[sub]-fN/2)) ix2a[j] = iq2+1 - fN/2 + j;
      if (iq2+1 < fN/2) ix2a[j] = j;
      if (iq2+1 > (fNQ2[sub]-fN/2)) ix2a[j] = (fNQ2[sub]-fN) + j;

      // Check grids
      if (ix2a[j] < 0 || ix2a[j] >= fNQ2[sub])
	{
	  cout << "Error in grids! j, ix2a[j] = "
	       << j << "\t" << ix2a[j] << endl;
	  exit(-1);
	}
    }

  // define points where to evaluate interpolation
  // choose between linear or logarithmic (x,Q2) interpolation
  const double xch = 1e-1;

  double x1;
  if (x < xch) x1 = log(x);
  else x1 = x;
  double x2 = log(Q2);
  
  if (id < 0 || id > fNFL-1)
    {
      cout << "Error: flavor out of range:" << id << endl;
      exit(-1);
    }
  else
    {
      // Choose betwen linear or logarithmic (x,Q2) interpolation
      for (int i = 0; i < fM; i++)
	{
	  if (x < xch)
	    x1a[i] = fLogXGrid[ix1a[i]];
	  else
	    x1a[i] = fXGrid[ix1a[i]];
	  
	  for (int j = 0; j < fN; j++)
	    {
	      x2a[j] = fLogQ2Grid[sub][ix2a[j]];
	      ya[i][j] = fPDFGrid[sub][fRep][id][ix1a[i]][ix2a[j]];
	    }
	}
      
      // 2D polynomial interpolation
      double y = 0, dy = 0;
      lh_polin2(x1a,x2a,ya,x1,x2,y,dy);
      res = y;
    }

  return res;
}

void NNPDFDriver::lh_polin2(double x1a[], double x2a[],
			    double ya[][fN],
			    double x1, double x2, 
			    double& y, double& dy)
{
  double yntmp[fN];
  double ymtmp[fM];

  for (int j = 0; j < fM; j++)
    {
      for(int k = 0; k < fN; k++)
	yntmp[k] = ya[j][k];
      
      lh_polint(x2a,yntmp,fN,x2,ymtmp[j],dy);
    }
  lh_polint(x1a,ymtmp,fM,x1,y,dy);
}

void NNPDFDriver::lh_polint(double xa[], double ya[], int n, double x,
			    double& y, double& dy)
{
  int ns = 0;  
  double dif = abs(x-xa[0]);
  double c[fM > fN ? fM : fN];
  double d[fM > fN ? fM : fN];
  
  for (int i = 0; i < n; i++)
    {
      double dift = abs(x-xa[i]);
      if (dift < dif)
	{
	  ns = i;
	  dif = dift;
	}
      c[i] = ya[i];
      d[i] = ya[i];
    }
  y = ya[ns];
  ns--;
  for (int m = 1; m < n; m++)
    {
      for (int i = 0; i < n-m; i++)
	{
	  double ho = xa[i]-x;
	  double hp = xa[i+m]-x;
	  double w = c[i+1]-d[i];
	  double den = ho-hp;
	  if (den == 0)
	    {
	      cout << "failure in polint" << endl;
	      exit(-1);	       
	    }
	  den = w/den;
	  d[i] = hp*den;
	  c[i] = ho*den;
	}
      if (2*(ns+1) < n-m)
	dy = c[ns+1];
      else {
	dy = d[ns];
	ns--;
      }
      y+=dy;
    }
}
