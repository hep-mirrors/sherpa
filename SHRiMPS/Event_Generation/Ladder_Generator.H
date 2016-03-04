#ifndef SHRIMPS_Event_Generation_Ladder_Generator_H
#define SHRIMPS_Event_Generation_Ladder_Generator_H

#include "SHRiMPS/Event_Generation/Final_State.H"
#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"

namespace SHRIMPS {
  class Ladder_Generator {
  private:
    Final_State      m_FS;
    Colour_Generator m_Colours;
    double           m_Y, m_Ymax;
    Ladder         * p_ladder;

    Omega_ik * p_eikonal;
    double     m_B, m_E[2], m_mandels, m_b1, m_b2;

    void Reset();
    void InitLadder(ATOOLS::Blob * blob);
    void AddInParticles();
    void FillBlob(ATOOLS::Blob * blob);
    bool CheckTotalMomentum();
    
    void TestPositioning(const std::string & dirname);
  public:
    Ladder_Generator();
    ~Ladder_Generator();
    
    void InitCollision(Omega_ik * eikonal,const double & B);
    bool MakePrimaryLadder(ATOOLS::Blob * blob);

    inline const double LeftoverEnergy(const int & beam) const {
      if (beam==0 || beam==1) return m_E[beam];
      return -1.;
    }
    
    void Test(const std::string & dirname);
  };
}
#endif
