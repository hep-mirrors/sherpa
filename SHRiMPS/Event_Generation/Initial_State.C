#include "SHRiMPS/Event_Generation/Initial_State.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;

void Initial_State::
InitNewCollision(Omega_ik * eikonal,const double & B,Ladder * ladder) {
  p_eikonal        = eikonal;
  p_eikonalweights = eikonal->GetWeights();
  m_B = B;
}

bool Initial_State::ProvideIS(Ladder * ladder) {
 return true; 
}

void Initial_State::
DefineIS(Ladder_Particle *& lpart1,Ladder_Particle *& lpart2,bool isresc) {
  if (isresc) {
    m_pos = (lpart1->m_pos+lpart2->m_pos)/2.;
    m_b1  = sqrt(m_pos[1]*m_pos[1]+m_pos[2]*m_pos[2]);
    m_b2  = sqrt((m_B-m_pos[1])*(m_B-m_pos[1])+m_pos[2]*m_pos[2]);
  }
  else m_pos = p_eikonal->SelectB1B2(m_b1,m_b2,m_B);
  m_ycms   = (lpart1->m_mom + lpart2->m_mom).Y();
  m_shat   = (lpart1->m_mom+lpart2->m_mom).Abs2();
  double intercept(p_eikonalweights->EffectiveIntercept(m_b1,m_b2,0.));
  m_weight = p_lumi->Weight(m_shat,intercept);
}



