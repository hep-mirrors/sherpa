#ifndef SHRIMPS_Main_Inelastic_Event_Generator_H
#define SHRIMPS_Main_Inelastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Cross_Sections/Sigma_D.H"
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
    Sigma_Elastic * p_sigma_el;
    Sigma_D * p_sigma_sd;
    Primary_Ladders   m_primaries;
    bool m_mustinit;
    int  m_test;

    double m_abs_t;
    
    std::map<Omega_ik *,double>                 m_xsecs; 
    std::map<Omega_ik *,std::vector<double> * > m_Bgrids; 
    double m_sigma, m_B;
    int    m_Nladders, m_Nprim, m_Ngen;

    Cluster_Algorithm * p_cluster;

    Collinear_Emission_Generator * p_collemgen;
        
    bool SelectEikonal();
    bool SelectB();
    ATOOLS::Blob * MakePrimaryScatterBlob();
  public:
    Inelastic_Event_Generator(Sigma_Inelastic * sigma,const int & test,Sigma_Elastic * sigma_el,Sigma_D * sigma_sd);
    ~Inelastic_Event_Generator();

    void Reset();
    void Initialise(Remnant_Handler * remnants,Cluster_Algorithm * cluster);
    int InitEvent(ATOOLS::Blob_List * blobs);
    ATOOLS::Blob * GenerateEvent();


    inline Colour_Generator * GetColourGenerator()  {
      return m_primaries.GetColourGenerator();
    }
    inline void   SetCluster(Cluster_Algorithm * cluster) {
    }
    inline void SetMaxEnergies(const double & E1,const double & E2) {
      m_primaries.SetMaxEnergies(E1,E2);
    }
    inline double XSec(Omega_ik * eikonal=NULL) {
      if (eikonal!=NULL) return m_xsecs[eikonal];
      return m_sigma;
    }
    inline const double B()      const { return m_B; }
    inline const double Yhat()   const { return m_primaries.Yhat(); }
    inline const double Ymax()   const { return m_primaries.Ymax(); }
    inline const double KT2min() const { return m_primaries.KT2min(); }
    void Test(const std::string & dirname);
  };
}
#endif
