#ifndef SHRIMPS_Main_Inelastic_Event_Generator_H
#define SHRIMPS_Main_Inelastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Ladders/Primary_Ladders.H"
#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "SHRiMPS/Event_Generation/Collinear_Emission_Generator.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace SHRIMPS {
  class Cluster_Algorithm;

  class Inelastic_Event_Generator : public Event_Generator_Base {
  private:
    Sigma_Inelastic * p_sigma;
    Primary_Ladders   m_primaries;
    bool m_mustinit;
    int  m_test;
    
    std::map<Omega_ik *,double>                 m_xsecs; 
    std::map<Omega_ik *,std::vector<double> * > m_Bgrids; 
    double m_sigma, m_B;
    int    m_Nladders, m_Nprim, m_Ngen;

    Cluster_Algorithm * p_cluster;

    Collinear_Emission_Generator * p_collemgen;
        
    bool InitEvent(ATOOLS::Blob_List * blobs);
    bool SelectEikonal();
    bool SelectB();
    bool MakePrimaryScatterBlobs(ATOOLS::Blob_List * blobs);
  public:
    Inelastic_Event_Generator(Sigma_Inelastic * sigma,const int & test);
    ~Inelastic_Event_Generator();

    void Reset();
    void Initialise(Remnant_Handler * remnants,Cluster_Algorithm * cluster);
    int  GenerateEvent(ATOOLS::Blob_List * blobs,const bool & flag=false);


    inline Colour_Generator * GetColourGenerator()  {
      return m_primaries.GetColourGenerator();
    }
    inline void   SetCluster(Cluster_Algorithm * cluster) {
    }
    inline double XSec(Omega_ik * eikonal=NULL) {
      if (eikonal!=NULL) return m_xsecs[eikonal];
      return m_sigma;
    }
    inline const double B()      const { return m_B; }
    inline const double Ymax()   const { return m_primaries.Ymax(); }
    inline const double KT2min() const { return m_primaries.KT2min(); }
    void Test(const std::string & dirname);
  };
}
#endif
