#ifndef SHRIMPS_Main_Single_Diffractive_Event_Generator_H
#define SHRIMPS_Main_Single_Diffractive_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SHRiMPS/Cross_Sections/Sigma_SD.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include <vector>
#include <map>

namespace SHRIMPS {
  class Single_Diffractive_Event_Generator : public Event_Generator_Base{
  private:
    Sigma_SD        * p_sigma;
    Beam_Remnant_Handler * p_beams;
    ATOOLS::Flavour   m_beam1, m_beam2, m_out1, m_out2;                 
    ATOOLS::Vec4D     m_p1, m_p2, m_p1out, m_p2out;
    double            m_pl12, m_pl22, m_sign1;
    bool              m_needsboost;
    double            m_accu;
    int               m_test;
 
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    bool   FixKinematics();
  public:
    Single_Diffractive_Event_Generator();
    Single_Diffractive_Event_Generator(Sigma_SD * sigma,Beam_Remnant_Handler * beams,
				       const int & test=0);
    ~Single_Diffractive_Event_Generator();
				       
    bool SingleDiffractiveEvent(ATOOLS::Blob_List * blobs,
				const double & xsec=0.);
    const double & XSec() const { return p_sigma->XSec(); }
    Omega_ik * GetEikonal() const { return NULL; }
    const double Smin()   const { return 0.; }
  };
}
#endif
