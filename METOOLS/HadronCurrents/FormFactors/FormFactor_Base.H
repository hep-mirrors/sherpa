#ifndef METOOLS_Current_Library_FormFactor_Base_H
#define METOOLS_Current_Library_FormFactor_Base_H

#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>

namespace METOOLS {
  enum class ff_model {
    none      = 0,
    KS        = 1,
    RChiPT    = 2,
    unknown   = 999
  };

  class GeneralModel;
  
  struct FF_Parameters {
    ff_model                       m_ffmodel;
    const ATOOLS::Flavour_Vector & m_flavs;
    std::vector<double>            m_masses, m_masses2;
    std::vector<int>               m_pi;
    std::string                    m_name;
    GeneralModel                 * p_model;

    FF_Parameters(const ff_model & ffmodel,
		  const ATOOLS::Flavour_Vector & flavs,
		  const std::vector<int> & decayindices,
		  const std::string & name,
		  GeneralModel * model=NULL) :
      m_ffmodel(ffmodel), m_flavs(flavs), m_pi(decayindices), m_name(name),
      p_model(model)
    {
      m_masses.resize(m_flavs.size());
      m_masses2.resize(m_flavs.size());      
      for (size_t i=0;i<m_flavs.size();i++) {
	m_masses[i]  = m_flavs[i].Mass(true);
	m_masses2[i] = m_masses[i]*m_masses[i]; 
      }
    }
  };
  
  class FormFactor_Base {
  protected:
    ff_model               m_ffmodel;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<double>    m_masses, m_masses2;
    std::vector<int>       m_pi;
    std::string            m_name;
    GeneralModel         * p_model;
  public:
    FormFactor_Base(const FF_Parameters & params);
    inline virtual ~FormFactor_Base() {};
    virtual Complex operator()(const ATOOLS::Vec4D_Vector& moms) = 0;
    inline void        SetName(std::string name) { m_name = name; }
    inline std::string Name() const              { return m_name; }
  };
}

#define DECLARE_FF_GETTER(CLASS,TAG)				\
  DECLARE_GETTER(CLASS,TAG,METOOLS::FormFactor_Base,METOOLS::FF_Parameters); \
  void ATOOLS::Getter<METOOLS::FormFactor_Base,METOOLS::FF_Parameters,CLASS>:: \
  PrintInfo(std::ostream &str,const size_t width) const { str<<"FF::"<<TAG; }

typedef ATOOLS::Getter_Function
<METOOLS::FormFactor_Base,METOOLS::FF_Parameters> FF_Getter;

#endif
