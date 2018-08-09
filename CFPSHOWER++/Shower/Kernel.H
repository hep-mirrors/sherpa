#ifndef CFPSHOWER_Shower_Kernel_H
#define CFPSHOWER_Shower_Kernel_H

#include "CFPSHOWER++/Shower/Splitting.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include <iostream>

namespace CFPSHOWER {
  class Kernel {
  private:
    kernel_type::code      m_type;
    ATOOLS::Flavour_Vector m_flavs;
    bool                   m_swapped;
    
    SF_Base    * p_sfunction; 
    Gauge_Base * p_gauge;

    Weight * MakeWeight(const Splitting & split,const double & overfac=1.);
  public:
    Kernel(const Kernel_Info & info);
    inline ~Kernel() {
      delete p_sfunction;
      delete p_gauge;
    }
    
    double Integral(Splitting & split,const ATOOLS::Mass_Selector * msel);
    bool   Generate(Splitting & split,const ATOOLS::Mass_Selector * msel,
		    const double & overfac=1.);
    bool   FillOffsprings(Splitting & split);

    inline void SetSF(SF_Base * sfunction)       { p_sfunction = sfunction; }
    inline void SetGauge(Gauge_Base * gaugepart) { p_gauge     = gaugepart; }
    
    inline SF_Base    * GetSF()    const { return p_sfunction; } 
    inline Gauge_Base * GetGauge() const { return p_gauge; }
    inline const ATOOLS::Flavour_Vector & GetFlavs() const { return m_flavs; }
  };
  
  std::ostream &operator<<(std::ostream &s,const Kernel & kernel);

  typedef ATOOLS::Getter_Function<Kernel,Kernel_Info,std::less<std::string> > Kernel_Getter;
  typedef std::vector<Kernel *> Kernel_Vector;
}

#endif
