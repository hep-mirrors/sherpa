#include "Vector.H"

using namespace ATOOLS;

namespace HADRONS {

// elementary matrices

double G( int mu, int nu ); 
double E( short a, short b, short c, short d ); 
double Delta( int mu, int nu ); 

// Traces with 2 gammas
double Trace( int mu, int nu ); 
double Trace( Vec4D a, int B ); 
double Trace( Vec4D a, Vec4D b ); 

// Traces with 4 gammas
double Trace( int mu, int nu, int rho, int sigma ); 
double Trace( Vec4D a, int B, int C, int D ); 
double Trace( Vec4D a, Vec4D b, int C, int D ); 
double Trace( Vec4D a, int C, Vec4D b, int D ); 
double Trace( Vec4D a, Vec4D b, Vec4D c, int D ); 
double Trace( Vec4D a, Vec4D b, Vec4D c, Vec4D d ); 

// Traces with 4 gammas and gamma5
Complex Trace5( int mu, int nu, int rho, int sigma ); 
Complex Trace5( Vec4D a, int B, int C, int D ); 
Complex Trace5( Vec4D a, Vec4D b, int C, int D ); 
Complex Trace5( Vec4D a, int C, Vec4D b, int D ); 
Complex Trace5( Vec4D a, Vec4D b, Vec4D c, int D ); 
Complex Trace5( Vec4D a, Vec4D b, Vec4D c, Vec4D d ); 

// Traces with 6 gammas
double Trace( int mu, int nu, int rho, int sigma, int alpha, int beta ); 
double Trace( Vec4D a, int mu, int nu, Vec4D b, int rho, int sigma ); 

// Traces with 6 gammas and gamma5
Complex Trace5( int mu, int nu, int rho, int sigma, int alpha, int beta ); 
Complex Trace5( Vec4D a, int mu, int nu, Vec4D b, int rho, int sigma ); 

// Propagator structure function
double Zeta( Vec4D a, Vec4D b, Vec4D c, double MP2 ); 
double Zeta( int mu, int nu, Vec4D c, double MP2 ); 

/*!
  \file Traces.H
  \brief Tools to calculates traces of gamma matrices

  This file can be found in directory \c ME_Library.
*/  
} // end of namespace
