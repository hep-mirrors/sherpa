#ifndef HadTools_h
#define HadTools_h
#include "MathTools.H"
#include "Vector.H"
#include "ComplexVector.H"
#include "Matrix.H"
#include "Flavour.H"
#include "Message.H"

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif


namespace HADRONS
{

  struct Flavour_Info {
    int nout;
    ATOOLS::Flavour* flavs;
    int* indices;
  };

  typedef std::map<std::string,double> String_Map;
  
  struct GeneralModel: public String_Map {
  public:
    inline double operator()(const std::string &tag,const double &def) const
    {
      String_Map::const_iterator fit(find(tag));
//       if( fit==end() ) {
//         std::cout<<METHOD<<" Warning: Did not find a value for "<<tag<<" in "
//           <<"decay channel file. Will use default value of "<<def<<std::endl;
//         return def;
//       }
//       return fit->second;
      return fit!=end()?fit->second:def;
    }
  };

  class Tools {
  public:
    static double   Lambda(double a,double b,double c);
    static Complex  BreitWigner( double s, double Mass2, double MassWidth );
    static Complex  BreitWignerFix( double s, double Mass2, double MassWidth );
    static ATOOLS::Vec4D RealBosonPolarizationVector(
        ATOOLS::Vec4D p, int lambda, double M2, bool & iszero );
    static ATOOLS::Vec4D* ComplexBosonPolarizationVector(ATOOLS::Vec4D p, int lambda, double M2);
    static ATOOLS::Vec4D*  ComplexBosonPolarizationVector(ATOOLS::Vec4D p, int lambda);

    static ATOOLS::ComplexVec4D ComplexBosonPolarizationVectorC(ATOOLS::Vec4D p, int lambda, double M2);
    static ATOOLS::ComplexVec4D ComplexBosonPolarizationVectorC(ATOOLS::Vec4D p, int lambda);
    static ATOOLS::CMatrix ComplexSpin2BosonPolarizationVectorC(ATOOLS::Vec4D p, int lambda, double M2 );
    static ATOOLS::CMatrix ComplexSpin2BosonPolarizationVectorC(ATOOLS::Vec4D p, int lambda );
    static ATOOLS::Vec4D Cross( ATOOLS::Vec4D a, ATOOLS::Vec4D b, ATOOLS::Vec4D c );
//     static double Epsilon( short a, short b, short c, short d );
    
    
    
  }; 
    /*!
      \file Tools.H
      \brief Declares the class HADRONS::Tools
      */
    /*!
      \class Tools
      \brief A small collection of frequently used functions for ME's

      This class is completely static so that it does not need a constructor.
      Its methods can be called directly.
      */

    /*!
    \fn Tools::RealBosonPolarizationVector( ATOOLS::Vec4D p, int lambda, double M2, bool & iszero )
    \brief Calculates real linear polarization vectors of a massive spin-1 boson.

    Calculations according to the Diploma thesis by Tanju Gleisberg, p. 35
    Input: 
      - <var>p</var>: momentum of particle
      - <var>M2</var>: (offshell) mass
      - <var>lambda</var>
        - \f$0 \to \epsilon^s\f$
        - \f$1 \to \epsilon^l\f$
        - \f$2 \to \epsilon^{t_1}\f$
        - \f$3 \to \epsilon^{t_2}\f$
        .
      - <var>iszero</var>: indicates whether the returned polarisation vector is zero.
      .
    Output:
      - real polarization vector as  Vec4D
      .
    */
    
    /*!
    \fn Tools::ComplexBosonPolarizationVector( ATOOLS::Vec4D p, int lambda, double M2, ATOOLS::Vec4D * eps )
    \brief Calculates complex linear polarization vectors of a massive spin-1 boson.

    Calculations according to the Diploma thesis by Tanju Gleisberg, p. 35
    Input: 
      - <var>p</var>: momentum of particle
      - <var>M2</var>: (offshell) mass
      - <var>lambda</var>
        - \f$0 \to \epsilon^+\f$
        - \f$1 \to \epsilon^0\f$
        - \f$2 \to \epsilon^-\f$
        - \f$3 \to \epsilon^s\f$
        .
      .
    Output:
      - real part of polarization vector as Vec4D[0]
      - imaginary part of polarization vector as Vec4D[1]
      .
    */

    /*!
    \fn Tools::ComplexBosonPolarizationVector( ATOOLS::Vec4D p, int lambda, ATOOLS::Vec4D * eps )
    \brief Calculates complex linear polarization vectors of a massive spin-1 boson onshell.
      \sa Tools::ComplexBosonPolarizationVector( ATOOLS::Vec4D p, int lambda, double M2, ATOOLS::Vec4D * eps )
    */
  /*!
    \fn Tools::Lambda( double a, double b, double c )
    \brief Returns the value of the triangulation function
    */
  /*!
    \fn Tools::Cross( Vec4D a, Vec4D b, Vec4D c )
    \brief Returns 4-vector cross product

    This method returns the 4-vector
    \f[ p^\mu = \epsilon^{\mu \alpha \beta \gamma} a_\alpha b_\beta c_\gamma, \f]
    where \f$a,b,c\f$ are 4-vectors and \f$\epsilon\f$ the totally antisymmetric \f$\epsilon\f$-tensor.
    */

}

#endif
