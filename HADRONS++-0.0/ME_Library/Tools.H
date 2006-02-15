#ifndef HadTools_h
#define HadTools_h
#include "MathTools.H"
#include "Vector.H"

namespace HADRONS
{  
  class Tools {
  public:
    static double   Lambda(double,double,double);
    static Complex  BreitWigner( double s, double Mass2, double MassWidth );
    static Complex  BreitWigner( double s, double Mass2, 
				 double Width, double ms );
    static Complex  BreitWigner( double s, double Mass2, 
				 double Width, double ms1, double ms2 );
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms );
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms1, double ms2 );
    static ATOOLS::Vec4D RealBosonPolarizationVector(
        ATOOLS::Vec4D p, int lambda, double M2, bool & iszero );
    static void ComplexBosonPolarizationVector(
        ATOOLS::Vec4D p, int lambda, double M2, ATOOLS::Vec4D* eps);
    static void  ComplexBosonPolarizationVector(
        ATOOLS::Vec4D p, int lambda, ATOOLS::Vec4D* eps);

    static ATOOLS::Vec4D Cross( ATOOLS::Vec4D a, ATOOLS::Vec4D b, ATOOLS::Vec4D c );
    static double Epsilon( short a, short b, short c, short d );
    
    
    
    /*!
    \fn Tools::RealBosonPolarizationVector( ATOOLS::Vec4D p, int lambda, double M2, bool & iszero )
    \brief Calculates real linear polarization vectors of a massive spin-1 boson.

    Calculations according to Diploma Tanju Gleisberg, p. 35
    Input: 
      p: momentum of particle
      M2: (offshell) mass
      lambda 0 -> eps^s, 1 -> eps^l, 2 -> eps^t1, 3 -> eps^t2
      iszero: indicates whether the returned pol. vector is zero.
    Output:
      real polarization vector as Vec4D
    */
    
    /*!
    \fn Tools::ComplexBosonPolarizationVector( ATOOLS::Vec4D p, int lambda, double M2 )
    \brief Calculates complex linear polarization vectors of a massive spin-1 boson.

    Calculations according to Diploma Tanju Gleisberg, p. 31
    Input:
      p: momentum of particle
      M2: mass if offshell
      lambda 0 -> eps^0, 1 -> eps^+, 2 -> eps^-, 3 -> eps^s
    Output:
      real part of polarization vector as Vec4D[0]
      imaginary part of polarization vector as Vec4D[1]
    */
  }; 

}

#endif
