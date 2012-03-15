#ifndef SHRIMPS_Main_Quasi_Elastic_Event_Generator_H
#define SHRIMPS_Main_Quasi_Elastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Event_Generation/Elastic_Event_Generator.H"
#include "SHRiMPS/Event_Generation/Single_Diffractive_Event_Generator.H"
#include "SHRiMPS/Event_Generation/Double_Diffractive_Event_Generator.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <vector>

namespace SHRIMPS {
  class Quasi_Elastic_Event_Generator : public Event_Generator_Base , public ATOOLS::Function_Base {
  private:
    Elastic_Event_Generator            * p_elastic;
    Single_Diffractive_Event_Generator * p_singdiff;
    Double_Diffractive_Event_Generator * p_doubdiff;
    double                               m_el, m_sd, m_dd, m_stot;
  public:
    Quasi_Elastic_Event_Generator() {}
    Quasi_Elastic_Event_Generator(Elastic_Event_Generator * elastic,
				  Single_Diffractive_Event_Generator * sdiff,
				  Double_Diffractive_Event_Generator * ddiff);

    bool QuasiElasticEvent(ATOOLS::Blob_List * blobs,const double & xsec=0.);
    const double & XSec() const { return m_stot; }
  };
}
#endif
