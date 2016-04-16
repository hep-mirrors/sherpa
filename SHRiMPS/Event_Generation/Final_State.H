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
    size_t   m_Nladders;
    double   m_E[2], m_pp, m_pm;
    ATOOLS::Vec4D m_totkt;
    
    void   AddInitialRapiditiesToLadder();
    void   AddPropagators();
    void   SanitizePropagators();
    void   ConstructKinematics();
    void   ConstructPropagator(const double & y1,const double & y2,
			       TPropList::iterator & prop);
    void   ConstructEmission(const LadderMap::iterator & emit,
			     const ATOOLS::Vec4D & ktvec);
    colour_type::code PropColour(const double & y1,const double & y2);
    void   CheckNextPropagator(TPropList::iterator & prop,
			       LadderMap::iterator & lit2);
    double SelectQT2(const double & qt2max,const double & expo);
    double QT2(double qt2max,const double & y1,const double & y2,
	       TPropList::iterator & prop);
  public:
    Final_State();
    ~Final_State();
    void   SetEikonal(Omega_ik * eikonal);
    void   SetNLadders(const size_t & N); 
    void   SetLadder(Ladder * ladder);
    void   SetImpactParameters(const double & b1, const double & b2);
    void   SetAvailableEnergies(const double E[2]);
    void   FillPrimaryLadder();

    const double ASratio(const double & that) const {
      return (*p_alphaS)(that)/p_alphaS->MaxValue();
    }
    const double Pplus()  const { return m_pp; }
    const double Pminus() const { return m_pm; }

    void Test(const std::string & dirname);
  };
}
#endif
