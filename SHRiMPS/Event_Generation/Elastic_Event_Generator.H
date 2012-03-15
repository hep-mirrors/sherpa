#ifndef SHRIMPS_Main_Elastic_Event_Generator_H
#define SHRIMPS_Main_Elastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <vector>

namespace SHRIMPS {
  class Elastic_Event_Generator : public Event_Generator_Base {
  private:
    Sigma_Elastic        * p_sigma;
    Beam_Remnant_Handler * p_beams;
    ATOOLS::Flavour        m_beam1, m_beam2;                 
    ATOOLS::Vec4D          m_p1, m_p2, m_p1out, m_p2out;
    double                 m_pl12, m_pl22, m_sign1;
    bool                   m_needsboost;
    double                 m_accu;
    int                    m_test;

    void   FixKinematics();
  public:
    Elastic_Event_Generator() {}
    Elastic_Event_Generator(Sigma_Elastic * sigma,Beam_Remnant_Handler * beams,
			    const int & test=0);

    bool ElasticEvent(ATOOLS::Blob_List * blobs,const double & xsec=0.);
    const double & XSec() const { return p_sigma->XSec(); }
    Omega_ik * GetEikonal() const { return NULL; }
    const double Smin()   const { return 0.; }
  };
}
#endif
