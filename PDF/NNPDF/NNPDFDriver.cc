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
 *    or
 *
 *  pdf->xfx(x,Q); // -> returns double*
 *
 *  // with double* index from [0,13], always in the LHAPDF format
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
#define NNDriverVersion "1.0.4"

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
  fMem(1),
  fRep(0),
  fAlphas(0),
  fXGrid(NULL),
  fLogXGrid(NULL),
  fHasPhoton(false),
  fSingleMem(false),
  fLHAPDF6(false),
  //HS
  fxmin(0),
  fxmax(0),
  fQmin(0),
  fQmax(0),
  fMz(0),  
  forder(0),
  // Masses
  fMDown(0),
  fMUp(0),  
  fMStrange(0),
  fMCharm(0),
  fMBottom(0), 
  fMTop(0)
  
{
  // Logo
  cout << " ****************************************" << endl;
  cout << " HS: TODO: find a way to load PDF set only once" << endl;
  cout << "      NNPDFDriver version " << NNDriverVersion << endl;
  cout << "  Grid: " << gridfilename << endl;
  cout << " ****************************************" << endl;

  // Check the grid version
  if (gridfilename.find(".LHgrid") == string::npos) fLHAPDF6 = true;

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
  if (fLHAPDF6)
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
	  //cout << "In LHAPDF6 --- Error: cannot open file " << grid << endl;
	  cout << "In LHAPDF6 --- Error: cannot open file " << file.str() << endl;
	  exit(-1);
	}

      string tmp;
      vector<string> splitstring;	      
      for (;;)
	{
	  getline(f,tmp);

	  if (tmp.find("SetDesc:") != string::npos)
	    cout << tmp << endl;
	  
          if (tmp.find("Flavors: [") != string::npos) 	  
	    if (tmp.find("22") != string::npos) { fHasPhoton = true;  fNFL++; }
	 
	  //if (tmp.find("NumMembers:") != string::npos)
	    //{
	      //split(splitstring,tmp);
	      //fMem = atof(splitstring[1].c_str())-1;
	    //}
	  //if (tmp.find("AlphaS_MZ:") != string::npos)
	    //{
	      //split(splitstring,tmp);
              //fAlphas = atof(splitstring[1].c_str());
	    //}
	  //if (tmp.find("XMin:") != string::npos)
	    //{
	      //split(splitstring,tmp);
	      //if ( fabs(atof(splitstring[1].c_str()) - XMINGRID) > 1e-10)
		//{
		  //cout << "Problem with XMINGRID" << endl;
		  //exit(-1);
		//}
              //// HS
              //fxmin = atof(splitstring[1].c_str());
		  // HS
          //if (tmp.find("XMax:") != string::npos)
            //{
	      //split(splitstring,tmp);
              //fxmax = atof(splitstring[1].c_str());
            //}
	  //if (tmp.find("QMin:") != string::npos)
	    //{
	      //split(splitstring,tmp);
              //// HS
              //fQmin = atof(splitstring[1].c_str());
	    //}
	  //// HS
          //if (tmp.find("QMax:") != string::npos)
            //{
	      //split(splitstring,tmp);
              //fQmax = atof(splitstring[1].c_str());
            //}    //}
          // HS
          //if (tmp.find("MZ:") != string::npos)
            //{
	      //split(splitstring,tmp);
              //fMz = atof(splitstring[1].c_str());
            //}
          //// HS
          //if (tmp.find("AlphaS_OrderQCD:") != string::npos)
            //{
	      //split(splitstring,tmp);
              //forder = atoi(splitstring[1].c_str());
            //}
	  
          
          if (hasKey(tmp, "NumMembers:")) fMem = readInt(tmp, "NumMembers:") -1;

          if (hasKey(tmp, "AlphaS_MZ:")) fAlphas = readDouble(tmp, "AlphaS_MZ:");

          if (hasKey(tmp, "XMin:")) {
            fxmin = readDouble(tmp, "XMin:");
	    if ( fabs(fxmin - XMINGRID) > 1e-10) {
              cout << "Problem with XMINGRID" << endl;
              exit(-1);
            }
          }
          if (hasKey(tmp, "XMax:")) fxmax = readDouble(tmp, "XMax:");
          if (hasKey(tmp, "QMin:")) fQmin = readDouble(tmp, "QMin:");
          if (hasKey(tmp, "QMax:")) fQmax = readDouble(tmp, "QMax:");
          if (hasKey(tmp, "MZ:")) fMz = readDouble(tmp, "MZ:");
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

      // single member switcher
      if (rep >= 0) { fMem = 0; fSingleMem = true; }

      if (fSingleMem)
	{
          //cout << "SingleMem" << endl;
	  stringstream file("");
	  if (rep < 10)
	    file << grid << "/" << name << "_000" << rep << ".dat";
	  else if (rep < 100)
	    file << grid << "/" << name << "_00" << rep << ".dat";
	  else if (rep < 1000)
	    file << grid << "/" << name << "_0" << rep << ".dat";
	  else
	    file << grid << "/" << name << "_" << rep << ".dat";

	  f.open(file.str().c_str(), ios::in);
          cout <<"Opening file " << file.str().c_str() << " to read PDF data " << endl;
	  
	  getline(f, tmp);
	  getline(f, tmp);
	  getline(f, tmp);

	  int sub = 0;
	  for (;;)
	    {	      
	      if (sub == 0)
		{
		  // reading Xgrid
		  getline(f, tmp);
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

	      getline(f, tmp);
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
	      getline(f, tmp);
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
		    for (int fl = 0; fl < fls.size(); fl++) 
		      f >> fPDFGrid[sub][imem][fls[fl]][ix][iq];	

	      getline(f, tmp);
	      getline(f, tmp);

	      if (tmp.find("---") != string::npos)
		{
		  sub++;
		  getline(f, tmp);
		  if (f.eof()) break;
		  continue;
		}
	    }
	  
	  f.close();
	}
      else
	{
	  for (int imem = 0; imem <= fMem; imem++)
	    {
	      stringstream file("");
	      if (imem < 10)
		file << grid << "/" << name << "_000" << imem << ".dat";
	      else if (imem < 100)
		file << grid << "/" << name << "_00" << imem << ".dat";
	      else if (imem < 1000)
		file << grid << "/" << name << "_0" << imem << ".dat";
	      else
		file << grid << "/" << name << "_" << imem << ".dat";

	      f.open(file.str().c_str(), ios::in);
	  
	      getline(f, tmp);
	      getline(f, tmp);
	      getline(f, tmp);
	      
	      int sub = 0;
	      for (;;)
		{	      
		  if (sub == 0 && imem == 0)
		    {
		      // reading Xgrid
		      getline(f, tmp);
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
		  else if (sub == 0) getline(f,tmp);
		  
		  if (imem == 0)
		    {
		      getline(f, tmp);
		      split(splitstring, tmp);
		      
		      fNQ2.push_back(splitstring.size());	    
		      fQ2Grid.push_back(new double[fNQ2[sub]]);
		      fLogQ2Grid.push_back(new double[fNQ2[sub]]);
		  
		      for (int iq = 0; iq < fNQ2[sub]; iq++)
			{
			  fQ2Grid[sub][iq] = pow(atof(splitstring[iq].c_str()), 2.0);
			  fLogQ2Grid[sub][iq] = log(fQ2Grid[sub][iq]);	       
			}

		      fPDFGrid.push_back(new double***[fMem+1]);     
		    }
		  else getline(f, tmp);

		  // skip flavor line
		  getline(f, tmp);
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
		  
		  // read PDF grid points	      
		  for (int ix = 0; ix < fNX; ix++)
		    for (int iq = 0; iq < fNQ2[sub]; iq++)
		      for (int fl = 0; fl < fls.size(); fl++) 
			f >> fPDFGrid[sub][imem][fls[fl]][ix][iq];	
		  
		  getline(f, tmp);
		  getline(f, tmp);
		  
		  if (tmp.find("---") != string::npos)
		    {
		      sub++;
		      getline(f, tmp);
		      if (f.eof()) break;
		      continue;
		    }
		}
	      
	      f.close();	      
	    }
	}

    }
  else
    {
      // reading Q2 and setting subgrids to 0 because it is a LHAPDF5 grid
      int sub = 0;
      
      fstream f;
      f.open(grid.c_str(),ios::in);
      if (f.fail())
	{
	  cout << "Error: cannot open file " << grid << endl;
	  exit(-1);
	}
      
      // Read header
      string tmp;
      for (;;) {
	getline(f,tmp);
	if (tmp.find("Alphas") != string::npos) break;
	cout << tmp << endl;
      }
      cout << endl;
      
      // Extract alphas and total number of members
      for (;;) {
	getline(f,tmp);
	if (tmp.find("Parameterlist:") != string::npos){
	  f >> tmp >> fMem >> tmp >> tmp >> fAlphas;
	  break;
        }
      }
      
      // Reading grid: removing header
      for (;;)
	{
	  getline(f,tmp);
	  if (tmp.find("NNPDF20intqed") != string::npos)
	    {
	      fHasPhoton = true;
	      fNFL++;
	      getline(f,tmp);
	      break;
	    }
	  else if (tmp.find("NNPDF20int") != string::npos)
	    {
	      getline(f,tmp);
	      break;
	    }
	}

      // getting nx
      f >> fNX;
      fXGrid = new double[fNX];
      for (int ix = 0; ix < fNX; ix++)
	f >> fXGrid[ix];
      
      fLogXGrid = new double[fNX];
      for (int ix = 0; ix < fNX; ix++)
	fLogXGrid[ix] = log(fXGrid[ix]);
      
      fNQ2.push_back(50);
      f >> fNQ2[sub];
      
      f >> tmp; // skip first value
      fQ2Grid.push_back(new double[fNQ2[sub]]);
      for (int iq = 0; iq < fNQ2[sub]; iq++)
	f >> fQ2Grid[sub][iq];
      
      fLogQ2Grid.push_back(new double[fNQ2[sub]]);
      for (int iq = 0; iq < fNQ2[sub]; iq++)
	fLogQ2Grid[sub][iq] = log(fQ2Grid[sub][iq]);
      
      // Prepare grid array
      if (rep >= 0) { fMem = 0; fSingleMem = true; }
      
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
      
      if (fSingleMem)
	{
	  // skip replicas before reading
	  for (int imem = 0; imem < rep; imem++)
	    for (int ix = 0; ix < fNX; ix++)
	      for (int iq = 0; iq < fNQ2[sub]; iq++)
		for (int fl = 0; fl < fNFL; fl++) 
		  f >> tmp;	      
	}
      
      // read PDF grid points
      f >> tmp; // remove replica number
      for (int imem = 0; imem <= fMem; imem++)
	for (int ix = 0; ix < fNX; ix++)
	  for (int iq = 0; iq < fNQ2[sub]; iq++)
	    for (int fl = 0; fl < fNFL; fl++) 
	      f >> fPDFGrid[sub][imem][fl][ix][iq];	
      
      f.close();
    }
}

/**
 * @brief NNPDFDriver::xfx
 * @param x
 * @param Q
 * @return
 */
double NNPDFDriver::xfx(double const& x, double const& Q)
//double* NNPDFDriver::xfx(double const& x, double const& Q)
{
  return xfxevolve(x,Q);
}

/**
 * @brief NNPDFDriver::xfx
 * @param x
 * @param Q
 * @param id
 * @return
 */
double NNPDFDriver::xfx(double const& x, double const& Q, int const& id)
{
  return xfxevolve(x,Q,id+6);
  //return xfxevolve(x,Q,id+6)[0];
}

/**
 * @brief NNPDFDriver::xfxevolve
 * @param x
 * @param Q
 * @param id
 * @return
 */
//double* NNPDFDriver::xfxevolve(double x, double Q, int id)
double NNPDFDriver::xfxevolve(double x, double Q, int id)
{
  double Q2 = Q*Q;
  int sub = 0;
  if (fLHAPDF6)
    {
      for (size_t i = 0; i < fNQ2.size(); i++)
	if ( Q2 >= fQ2Grid[i][0]) sub = i;
    }

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
      if (ix+1 >= fM/2 && ix+1 <= (fNX-fM/2)) ix1a[i] = ix+1 - fM/2 + i;
      if (ix+1 < fM/2) ix1a[i] = i;
      if (ix+1 > (fNX-fM/2)) ix1a[i] = (fNX-fM) + i;
	
      // Check grids
      if (ix1a[i] < 0 || ix1a[i] >= fNX)
	{
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

  //double xpdf[fNFL];// = double[fNFL];
  //double *xpdf = new double[fNFL];
  //for (int i = 0; i < fNFL; i++) xpdf[i] = 0.0;
  
  //if (id < 0 || id > fNFL-1)
    //{
      //for (int ipdf = 0; ipdf < fNFL; ipdf++)
	//{
	  //// Choose betwen linear or logarithmic (x,Q2) interpolation
	  //for (int i = 0; i < fM; i++)
	    //{
	      //if (x < xch)
		//x1a[i] = fLogXGrid[ix1a[i]];
	      //else
		//x1a[i] = fXGrid[ix1a[i]];
	      
	      //for (int j = 0; j < fN; j++)
		//{
		  //x2a[j] = fLogQ2Grid[sub][ix2a[j]];
		  //ya[i][j] = fPDFGrid[sub][fRep][ipdf][ix1a[i]][ix2a[j]];
		//}
	    //}
	  
	  //// 2D polynomial interpolation
	  //double y = 0, dy = 0;
	  //lh_polin2(x1a,x2a,ya,x1,x2,y,dy);
	  
	  //xpdf[ipdf] = y;
	//}
    //}
  //else
    //{
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
      
      //xpdf[0] = y;
    //}

  return y;
  //return xpdf;
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
