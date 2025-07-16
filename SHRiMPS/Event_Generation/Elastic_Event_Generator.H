#ifndef SHRIMPS_Main_Elastic_Event_Generator_H
#define SHRIMPS_Main_Elastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include <vector>
#include <map>

namespace SHRIMPS {
  class Elastic_Event_Generator : public Event_Generator_Base {
  private:
    Sigma_Elastic * p_sigma;
    ATOOLS::Flavour m_beam1, m_beam2;                 
    ATOOLS::Vec4D   m_p1, m_p2, m_p1out, m_p2out;
    ATOOLS::Blob  * p_blob; 
    double          m_E12, m_pl12, m_pl1, m_E22, m_pl22, m_pl2, m_abs_t;
    int             m_sign1;
    bool            m_ana;

    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    void FixKinematics();
    void FillBlob();
  public:
    Elastic_Event_Generator(Sigma_Elastic * sigma,const int & test=0);
    ~Elastic_Event_Generator();

    int InitEvent(ATOOLS::Blob_List * blobs); 
    ATOOLS::Blob * GenerateEvent() { return NULL; }
  };
}
#endif
