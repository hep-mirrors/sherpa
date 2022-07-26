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
    std::map<std::pair<ATOOLS::Flavour,ATOOLS::Flavour>,double> m_hadronmaps[2];
    ATOOLS::Blob  * p_blob;
    ATOOLS::Vec4D   m_p[2], m_pout[4];
    double          m_E[2], m_abs_t;
    double          m_rate[3], m_sigma, m_massexp;
    int             m_sign1, m_mode;
    double          m_Prob1440,  m_Prob1710;
    double          m_minmass2[2], m_maxmass2[2], m_expargUp[2], m_expargLow[2];
    bool            m_contMassRange[2];
    
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    void   InitialiseHadronMaps();
    void   SelectMode();
    void   SelectFS();
    void   SelectFlavours(size_t beam);
    void   FixKinematics();
    double SelectMass2(size_t beam);
    void   FixBinarySystem(ATOOLS::Vec4D (& moms)[2]);
    void   SplitQandQQ(size_t beam,ATOOLS::Vec4D & mom);
    void   FillBlob();
    //std::vector<ATOOLS::Vec4D> SplitIntoQandQQ(ATOOLS::Vec4D pmu, double Mqq, double Mq);
    //std::vector<double> ComputePxPyPz(double p, int sign1, int mode);
    //ATOOLS::Vec4D Get4Vector(double M2[], double Etot);
  public:
    Soft_Diffractive_Event_Generator(Sigma_SD * sigma,const int & test=0);
    ~Soft_Diffractive_Event_Generator();
				       
    int InitEvent(ATOOLS::Blob_List * blobs);
    ATOOLS::Blob * GenerateEvent() { return NULL; }
  };
}
#endif
