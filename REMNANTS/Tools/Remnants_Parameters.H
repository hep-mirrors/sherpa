#ifndef REMNANTS_Tools_Remnants_Parameters_H
#define REMNANTS_Tools_Remnants_Parameters_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>

namespace REMNANTS {
  struct pkform {
    enum code {
      none           = 0,
      gauss          = 1,
      gauss_limited  = 2,
      dipole         = 3,
      dipole_limited = 4,
      undefined      = 99
    };
  };
  struct pkrecoil {
    enum code {
      democratic     = 1,
      beam_vs_shower = 2,
      undefined      = 99
    };
  };

  struct pkparams {
    pkform::code   m_form;
    pkrecoil::code m_recoil;
    std::map<std::string, double> m_params;

    pkparams() {}
    pkparams(const pkparams & pkp) {
      m_form   = pkp.m_form;
      m_recoil = pkp.m_recoil;
      for (std::map<std::string, double>::const_iterator it=pkp.m_params.begin();
	   it!=pkp.m_params.end();it++)
	m_params[it->first] = it->second;
    }
    void Output();
  };
  
  class Remnants_Parameters {
  private:
    std::map<long unsigned int, pkparams * > m_params;
    
    void Initialize();
    pkform::code   SelectForm(const std::string & form=std::string("dipole_limited"));
    pkrecoil::code SelectRecoil(const std::string & form=std::string("beam_vs_shower"));
  public:
    Remnants_Parameters();
    ~Remnants_Parameters();

    void ReadParameters();
    const pkform::code   & GetForm(const ATOOLS::Flavour & beamflav);
    const pkrecoil::code & GetRecoil(const ATOOLS::Flavour & beamflav);
    const double         & operator()(const ATOOLS::Flavour & beamflav,const std::string & tag);
  };

  extern Remnants_Parameters * rempars;
}

#endif
