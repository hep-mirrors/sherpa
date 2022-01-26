#ifndef SHRIMPS_Main_Soft_Diffractive_Event_Generator_H
#define SHRIMPS_Main_Soft_Diffractive_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Cross_Sections/Sigma_SD.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include <vector>
#include <map>
#include <cmath>

namespace SHRIMPS {
  class Soft_Diffractive_Event_Generator : public Event_Generator_Base{
  private:
    Sigma_SD      * p_sigma;
    ATOOLS::Flavour m_beam[2], m_out[4];
    std::map<ATOOLS::Flavour,double> m_hadronmaps[2];
    ATOOLS::Vec4D   m_p[2], m_pout[4];
    double          m_E[2], m_abs_t;
    double          m_rate[3], m_sigma;
    int             m_sign1, m_mode;
    int m_contMassRange0 = 0;
    int m_contMassRange1 = 0;
    double m_Prob1440 = 0.4;
    double m_Prob1710 = 0.2;
    int m_contMassIndex;
    
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    void InitialiseHadronMaps();
    void SelectModeAndFSHadrons();
    ATOOLS::Flavour SelectFromMap(const size_t & hadron);
    void FixKinematics();
    std::vector<ATOOLS::Vec4D> SplitIntoQandQQ(ATOOLS::Vec4D pmu, double Mqq, double Mq);
    std::vector<double> ComputePxPyPz(double p, int sign1, int mode);
    ATOOLS::Vec4D Get4Vector(double M2[], double Etot);
    double drawMass();
    void FillBlob(ATOOLS::Blob * blob);
  public:
    Soft_Diffractive_Event_Generator(Sigma_SD * sigma,const int & test=0);
    ~Soft_Diffractive_Event_Generator();
				       
    int GenerateEvent(ATOOLS::Blob_List * blobs,const bool & flag=false);
  };
}
#endif
