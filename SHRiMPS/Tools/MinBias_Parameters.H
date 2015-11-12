#ifndef SHRIMPS_Tools_MinBias_Parameters_H
#define SHRIMPS_Tools_MinBias_Parameters_H

#include "SHRiMPS/Tools/Parameter_Structures.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <map>
#include <string>

namespace SHRIMPS {
  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class MinBias_Parameters {
  private:
    FormFactor_Parameters m_ff_params;
    Eikonal_Parameters    m_eik_params;
    run_mode::code        m_runmode;
    weight_mode::code     m_weightmode;
    ff_form::code         m_ffform;
    absorption::code      m_absorp;
    ladder_weight::code   m_ladderweight;
    ktform::code          m_ktform;
    ordering::code        m_ordering;
    resc_ktmin::code      m_resc_ktmin;
    resc_nosing::code     m_resc_nosing;
    resc_over_sing::code  m_resc_over_sing;
    resc_mode::code       m_rescmode;
    reconn_mode::code     m_reconnmode;
    StringDoubleMap       m_params;
    MODEL::asform::code   m_as_form;

    void FillFFParams(ATOOLS::Data_Reader * dr,const int test);
    void FillEikonalParams(ATOOLS::Data_Reader * dr,const int test);
  public:
    MinBias_Parameters();
    ~MinBias_Parameters();

    void   Init(ATOOLS::Data_Reader * dr);
    double operator()(std::string keyword);
    void   Set(std::string keyword, double value) {
      if (m_params.find(keyword)!=m_params.end()) m_params[keyword] = value;
      else {
	std::cerr<<"Error in MinBias_Parameters::Set: '"
		 <<keyword<<"' not found.\n";
      }
    }

    run_mode::code       RunMode()      const { return m_runmode; }
    weight_mode::code    WeightMode()   const { return m_weightmode; }
    ff_form::code        FF_Form()      const { return m_ffform; }  
    absorption::code     Absorption()   const { return m_absorp; }  
    ladder_weight::code  LadderWeight() const { return m_ladderweight; }
    ktform::code         KTForm()       const { return m_ktform; }
    ordering::code       Ordering()     const { return m_ordering; }
    resc_mode::code      RescMode()     const { return m_rescmode; }
    resc_ktmin::code     RescKTMin()    const { return m_resc_ktmin; }
    resc_nosing::code    RescNoSing()   const { return m_resc_nosing; }
    resc_over_sing::code RescOverSing() const { return m_resc_over_sing; }
    reconn_mode::code    ReconnMode()   const { return m_reconnmode; }
    MODEL::asform::code  AsForm()       const { return m_as_form; }

    const FormFactor_Parameters & GetFFParameters() const { 
      return m_ff_params; 
    }
    const Eikonal_Parameters & GetEikonalParameters() const { 
      return m_eik_params; 
    }

  };
  
  extern MinBias_Parameters MBpars;
}

#endif
