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
#include <type_traits>

namespace METOOLS {
  // Form-factor IDs are local to each family.  Universal values are:
  // 0 -> exactly zero, 1 -> exactly one, 99 -> invalid (always throws).
  // Families assign their physical models from 10 onwards.
  static constexpr int ff_zero = 0;
  static constexpr int ff_one = 1;
  static constexpr int ff_invalid = 99;

  class GeneralModel;
  
  struct FF_Parameters {
    int                            m_ffmodel;
    const ATOOLS::Flavour_Vector & m_flavs;
    std::vector<double>            m_masses, m_masses2;
    std::vector<int>               m_pi;
    std::map<std::string,double>   m_parameterMap;
    std::string                    m_name;
    GeneralModel                 * p_model;

    FF_Parameters(const int ffmodel,
		  const ATOOLS::Flavour_Vector & flavs,
		  const std::vector<int> & decayindices,
		  std::map<std::string,double> & parameterMap,
		  const std::string & name,
		  GeneralModel * model=NULL) :
      m_ffmodel(ffmodel), m_flavs(flavs), m_pi(decayindices),
      m_parameterMap(parameterMap), m_name(name),
      p_model(model)
    {
      // Validate reserved IDs before any getter can decline construction.
      if (m_ffmodel==ff_invalid || m_ffmodel<0 ||
          (m_ffmodel>1 && m_ffmodel<10))
        THROW(fatal_error,"Invalid FORM_FACTOR="+std::to_string(m_ffmodel)+
              "; 0=zero, 1=one, 99=invalid, family models start at 10.");
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
    int                    m_ffmodel;
    template<typename Enum>
    Enum ModelAs() const {
      static_assert(std::is_enum<Enum>::value, "Form-factor model must be an enum");
      return static_cast<Enum>(m_ffmodel);
    }
    static void CheckReservedModel(int model);

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
