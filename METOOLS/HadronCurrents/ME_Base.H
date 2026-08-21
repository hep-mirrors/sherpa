#ifndef METOOLS_HadronCurrents_ME_Base_H
#define METOOLS_HadronCurrents_ME_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <string>
#include <utility>
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#define DEFINE_ME_GETTER(CLASS,TAG)				       \
  DECLARE_GETTER(CLASS,TAG,METOOLS::ME_Base,METOOLS::ME_Parameters);    \
  ME_Base* ATOOLS::Getter<METOOLS::ME_Base,METOOLS::ME_Parameters,CLASS>:: \
  operator()(const METOOLS::ME_Parameters &parameters) const	      \
  { return new CLASS(parameters.flavs, parameters.indices,TAG); }

namespace METOOLS {
  class ME_Base : public Spin_Amplitudes {
  protected:
    std::string m_name;
    const ATOOLS::Flavour_Vector& m_flavs;
    std::vector<int>              m_pi;
    std::vector<double>           m_masses, m_masses2;

    virtual double lambdaNorm(const double,const double,const double);
  public:
    ME_Base(const ATOOLS::Flavour_Vector& flavs,
	    const std::vector<int>& decayindices,
	    const std::string& name);
    virtual ~ME_Base() {}
    virtual void Calculate(const ATOOLS::Vec4D_Vector& momenta,bool anti)=0;
    virtual void SetModelParameters(ATOOLS::Scoped_Settings& s);
    virtual void SetModelParameters(METOOLS::GeneralModel _md ) {};
    virtual bool SetColorFlow(std::vector<ATOOLS::Particle*> outparts,
			      int n_q, int n_g, bool anti) {
      return false;
    }
    virtual std::string Name() { return m_name; }
  };

  typedef ATOOLS::Getter_Function<ME_Base,METOOLS::ME_Parameters>
  ME_Getter_Function;
}


#endif
