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

#pragma once

#include <iostream>
#include <vector>
#include <string>
using std::string;
using std::vector;

class NNPDFDriver {

 private:

  // Interpolation order
  static const int fM = 4;
  static const int fN = 4;

  int fNFL;           //! Total flavour number
  int fNX;            //! Total number of x points in the grid
  vector<int> fNQ2;   //! Total number of Q2 points in the grid (subgrids)
  int fMem;           //! Total number of Members
  int fRep;           //! Select the current replica
  double fAlphas;     //! AlphaS value
  double *fXGrid;     //! x grid
  double *fLogXGrid;  //! x grid
  vector<double*>     fQ2Grid;    //! q2 grid
  vector<double*>     fLogQ2Grid; //! q2 grid
  vector<double****>  fPDFGrid;   //! PDF grid
  bool fHasPhoton;    //! bool with photon information
  bool fSingleMem;    //! bool which determines the constructor
  bool fLHAPDF6;      //! bool which determines the grid version
  // Metainfo read from .info
  double fxmin; // xmin accessor
  double fxmax; // xmax accessor
  double fQmin; // Qmin accessor
  double fQmax; // Qmax accessor
  double fMz; // well Z mass
  int forder; // order alphas  QCD
  // Masses
  double fMDown;
  double fMUp;  
  double fMStrange;
  double fMCharm;
  double fMBottom; 
  double fMTop;


 public:
  /// The constructor
  NNPDFDriver(string const& gridfilename = "", int const& rep = -1);
  /// The destructor
  ~NNPDFDriver();

  //! Init PDF member irep = [0, fMem]
  void initPDF(int irep);

  //! returns the x*pdf of flavour id = [-6,7], LHA order.
  double xfx(double const& x, double const& Q, int const& id);

  //! returns the x*pdf array with all LHA order.
  double xfx(double const& x, double const& Q);
  //double* xfx(double const& x, double const& Q);

  //! Get NFL method, returns total number of flavours
  int GetNFL() { return fNFL; }

  //! Get AlphaS method, returns the alphas at Mz
  double GetAlphaSMz() { return fAlphas; }

  // Get xmin
  double GetXMin() {return fxmin; }

  // Get xmax
  double GetXMax() {return fxmax; }
  
  // Get Qmin
  double GetQMin() {return fQmin; }

  // Get Qmax
  double GetQMax() {return fQmax; }
  
  // Get MZ
  double GetMz() { return fMz; }

  // Get order alphas
  int GetOrderAlphaS() { return forder; }


  //! Returns true if the set contains the photon PDF
  bool hasPhoton() { return fHasPhoton; }

  //! Return the total number of replicas
  int GetMembers() { return fMem+1; }

 private:
  /// Reads the PDF from file
  void readPDFSet(string const&, int const&);
  /// Apply the interpolation algorithm
  double xfxevolve(double,double,int = -10);
  //double* xfxevolve(double,double,int = -10);
  /// Performs the 2D polynomial interpolation
  void lh_polin2(double[],double[],double[][fN],
		 double,double,double&,double&);
  /// Performs the 1D polynomial interpolation
  void lh_polint(double[],double[],int,double,double&,double&);
};
