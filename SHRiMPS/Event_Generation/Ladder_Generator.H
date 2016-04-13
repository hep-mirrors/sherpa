#ifndef SHRIMPS_Event_Generation_Ladder_Generator_H
#define SHRIMPS_Event_Generation_Ladder_Generator_H

#include "SHRiMPS/Event_Generation/Final_State.H"
#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator {
  private:
    Final_State      m_FS;
    Colour_Generator m_Colours;
    double           m_Y, m_Ymax;
    Ladder         * p_ladder;

    Omega_ik * p_eikonal;
    double     m_B, m_b1, m_b2, m_Ecms, m_E[2];
    size_t     m_Nladders;
    
    void Reset();
    void InitLadder(ATOOLS::Blob * blob);
    void AddInParticles();
    bool BreakPrimaryLadderGenerationLoop(const bool & isfirst,bool & sucess);
    bool CheckTotalMomentum();
    bool AcceptLadderForHardness();
    void UpdateInitialEnergies();
    void FillBlob(ATOOLS::Blob * blob);
    
    void TestPositioning(const std::string & dirname);

    std::map<std::string,ATOOLS::Histogram * > m_histos;
    void AnalyseLadder();
  public:
    Ladder_Generator();
    ~Ladder_Generator();
    
    void InitCollision(Omega_ik * eikonal,const double & B,const size_t & N);
    bool MakePrimaryLadder(ATOOLS::Blob * blob,bool isfirst);

    inline const double LeftoverEnergy(const int & beam) const {
      if (beam==0 || beam==1) return m_E[beam];
      return -1.;
    }
    inline const double TMax() const {
      return p_ladder->That();
    }
    
    void Test(const std::string & dirname);
  };
}
#endif
