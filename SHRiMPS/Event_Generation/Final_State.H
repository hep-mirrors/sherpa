#ifndef SHRIMPS_Main_Final_State_H
#define SHRIMPS_Main_Final_State_H

#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Eikonals/Eikonal_Weights.H"
#include "MODEL/Main/Strong_Coupling.H"


namespace SHRIMPS {
  class Final_State :  public ATOOLS::Function_Base {
  private:
    MODEL::Strong_Coupling * p_alphaS;
    Eikonal_Weights          m_weights;
    Eikonal_Contributor    * p_omegaik, * p_omegaki;
    double   m_origY, m_Ymax, m_DeltaY;
    double   m_b1, m_b2;
    Ladder * p_ladder;
    double   m_availableE, m_E,m_pz;
    
    void   AddInitialRapiditiesToLadder();
    void   AddPropagators();
    void   ConstructKinematics();
    void   ConstructPropagator(const double & y1,const double & y2,
			       T_Prop* prop,const int & sign);
    colour_type::code PropColour(const double & y1,const double & y2);
    double SelectQT2(const double & qt2max,const double & expo);
    double QT2(const double & qt2max,const double & y1,const double & y2,
	       const colour_type::code & col);
  public:
    Final_State();
    ~Final_State();
    void   SetEikonal(Omega_ik * eikonal);
    void   SetLadder(Ladder * ladder);
    void   SetImpactParameters(const double & b1, const double & b2);
    void   SetAvailableEnergy(const double & E=-1.);
    void   FillPrimaryLadder();

    const double Pplus()  const { return (m_E+m_pz)/2.; }
    const double Pminus() const { return (m_E-m_pz)/2.; }

    void Test(const std::string & dirname);
  };
}
#endif
