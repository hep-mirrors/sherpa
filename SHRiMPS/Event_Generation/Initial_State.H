#ifndef SHRIMPS_Main_Initial_State_H
#define SHRIMPS_Main_Initial_State_H

#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Event_Generation/Parton_Luminosity.H"
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "ATOOLS/Phys/Blob_List.H"

#include <list>

namespace SHRIMPS {
  class Initial_State {
  private:
    Parton_Luminosity * p_lumi;
    Omega_ik          * p_eikonal;
    Eikonal_Weights   * p_eikonalweights;
    double              m_B,m_b1,m_b2;
    ATOOLS::Vec4D       m_pos;
    double              m_ycms,m_shat,m_weight;
    bool                m_isrescatter;
  public:
    Initial_State(Parton_Luminosity * lumi=NULL) : p_lumi(lumi) {}
    inline ~Initial_State() {}

    void InitNewCollision(Omega_ik * eikonal,const double & B,
			  Ladder * ladder=NULL);
    bool ProvideIS(Ladder * ladder);
    void DefineIS(Ladder_Particle *& lpart1,Ladder_Particle *& lpart2,
		  bool isresc=false);


    inline const ATOOLS::Vec4D & Pos()  const { return m_pos; }
    inline const double & B1()          const { return m_b1; }
    inline const double & B2()          const { return m_b2; }
    inline const double & B()           const { return m_B; }
    inline const double & Ycms()        const { return m_ycms; }
    inline const double & Weight()      const { return m_weight; }
    inline const double & Shat()        const { return m_shat; }
    inline const double & Smin()        const { return p_lumi->Smin(); }
    inline const bool   & IsRescatter() const { return m_isrescatter; }
  };
}
#endif
