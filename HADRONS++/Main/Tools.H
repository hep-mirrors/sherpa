#ifndef HADRONS_Main_Tools_H
#define HADRONS_Main_Tools_H

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

#include <map>
#include <iostream>
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"

namespace HADRONS
{

  struct Flavour_Info {
    int n;
    ATOOLS::Flavour* flavs;
    int* indices;
  };

  typedef std::map<std::string,double> String_Map;

  struct GeneralModel: public String_Map {
    std::map<kf_code, kf_code> m_aliases;
    
  public:
    inline double operator()(const std::string &tag,const double &def) const 
    {
      String_Map::const_iterator fit(find(tag));
      return fit!=end()?fit->second:def;
    }
  };

  class Tools {
  public:
    static double   Lambda(double a,double b,double c);
    static Complex  BreitWigner( double s, double Mass2, double MassWidth );
    static Complex  BreitWignerFix( double s, double Mass2, double MassWidth );
    
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms );
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms1, double ms2 );
    static const int k0;
    static const double Vud;
    static const double Vus;
    static const double Vub;
    static const double Vcd;
    static const double Vcs;
    static const double Vcb;
    static const double Vtd;
    static const double Vts;
    static const double Vtb;
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
    \fn Tools::Lambda( double a, double b, double c )
    \brief Returns the value of the triangulation function
    */
}

#endif
