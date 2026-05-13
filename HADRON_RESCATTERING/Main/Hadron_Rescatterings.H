#ifndef HADRON_RESCATTERING_Main_Hadron_Rescatterings_H
#define HADRON_RESCATTERING_Main_Hadron_Rescatterings_H

#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"
#include "HADRON_RESCATTERING/XSecs/MesonMeson.H"
#include "HADRON_RESCATTERING/XSecs/BaryonMeson.H"
#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "HADRON_RESCATTERING/XSecs/ScatteringBase.H"

#include "ATOOLS/Phys/Hadron_Collision.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace HADRON_RESCATTERING {
  class Hadron_Rescatterings {
  private:
    double            m_pref, m_S;
    bool              m_on;
    long unsigned int m_nfails;

    ATOOLS::Flavour  m_flav1, m_flav2;
    ScatteringBase * p_interaction, * p_mm, * p_bm, * p_bb;

    bool   SelectInteraction(ATOOLS::Particle * inpart1,ATOOLS::Particle * inpart2);
    bool   WillRescatter(const double & s,const double & b2);

    ATOOLS::Blob * ElasticRescatter(ATOOLS::Hadron_Collision * candidate);
  public:
    Hadron_Rescatterings(const bool & on);
    ~Hadron_Rescatterings();
    
    ATOOLS::Blob * operator()(ATOOLS::Hadron_Collision * candidate);
    void Initialize();
    void Reset();
  };
}

#endif
