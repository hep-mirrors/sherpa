#ifndef METOOLS_HadronCurrents_Resonance_Base_H
#define METOOLS_HadronCurrents_Resonance_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Lookup_Tables.H"
#include <set>
#include <string>
#include <vector>

namespace METOOLS {
  class Propagator_Base;
  
  class Partial_Width_Base {
  protected:
    ATOOLS::Flavour              m_inflav;
    std::vector<ATOOLS::Flavour> m_outflavs;
    std::vector<double>          m_decmasses, m_decmasses2;
    double                       m_mass, m_mass2, m_mmin, m_smin;
    double                       m_BR, m_width, m_width2, m_partialwidth;
    double                       m_prefactor;
    
    inline const double Flux(const double & s) const {
      return 1./(2.*sqrt(s));
    }
    virtual const double Calculate(const double & s) = 0;
    virtual void FixPrefactor();
  public:
    Partial_Width_Base(const ATOOLS::Flavour & inflav,
		       const std::vector<ATOOLS::Flavour> & outflavs,
		       const double & BR);
    virtual void   Init(const std::vector<Propagator_Base *> & props,
			const ATOOLS::axis & m_Qrange) = 0;
    virtual const double operator()(const double & s);
    virtual const double Gamma() const { return m_partialwidth; }
    virtual const double BR()    const { return m_BR; }
    virtual void SetPrefactor(const double & pref=1.) { m_prefactor = pref; }
    virtual void Output();
  };


  class Total_Width_Base {
  protected:
    ATOOLS::Flavour m_inflav;
    double          m_mass, m_mass2, m_width, m_width2;
    std::set<Partial_Width_Base *> m_channels;
  public:
    Total_Width_Base(const ATOOLS::Flavour & inflav);
    ~Total_Width_Base();
    
    const double operator()(const double & s);
    Partial_Width_Base * SelectChannel(const double & s=-1);

    void PrintBRs() const;
    void OutputLineshape(const double & mmin,const double &mmax,
			 const size_t & steps=100);
    inline ATOOLS::Flavour & Flav() { return m_inflav; }
  };
}



#endif

