#ifndef SHRIMPS_Main_Inelastic_Event_Generator_H
#define SHRIMPS_Main_Inelastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Event_Generation/Ladder_Generator.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace SHRIMPS {
  class Inelastic_Event_Generator : public Event_Generator_Base {
  private:
    Ladder_Generator m_laddergenerator;

    std::map<Omega_ik *,double>                 m_xsecs; 
    std::map<Omega_ik *,std::vector<double> * > m_Bgrids; 
    Omega_ik * p_eikonal;
    double     m_sigma,m_B;
    int        m_Nladders, m_Nprim, m_Ngen;
    bool       m_first, m_done;

    void Initialise();

    void InitInelasticEvent();
    void ResetEvent();
    bool SelectEikonal();
    bool SelectB();

    int  AddScatter(ATOOLS::Blob_List * blobs);
    ATOOLS::Blob * CreateBlob(ATOOLS::Blob_List * blobs);

    void TestSelectB(const std::string & dirname);
    void TestNumberOfLadders(const std::string & dirname);
    void TestInitialLadders(const std::string & dirname);
  public:
    Inelastic_Event_Generator();
    ~Inelastic_Event_Generator();

    void Reset();
    int  GenerateEvent(ATOOLS::Blob_List * blobs,const bool & isUE=false);

    inline double XSec(Omega_ik * eikonal) {
      if (eikonal!=NULL) return m_xsecs[eikonal];
      return m_sigma;
    }

    void Test(const std::string & dirname);
  };
}
#endif
