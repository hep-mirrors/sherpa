#ifndef SHRIMPS_Ladders_Ladder_Generator_Base_H
#define SHRIMPS_Ladders_Ladder_Generator_Base_H

#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "SHRiMPS/Eikonals/Rapidity_Density.H"
#include "SHRiMPS/Ladders/MEs.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Beam_Remnants/Remnant_Handler.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_Base {
  protected:
    Sigma_Partonic           m_partonic;
    double                   m_Ymax, m_deltaY;
    double                   m_beamY = m_Ymax + m_deltaY;
    double                   m_Yrange;
    double                   m_qt2min, m_qt2minFF, m_kt2min, m_kt2minShower;
    MODEL::Strong_Coupling * p_alphaS;
    Rapidity_Density         m_density;
    MEs                      m_me;
    Omega_ik               * p_eikonal;
    Ladder                 * p_ladder;
    LadderMap              * p_emissions;
    TPropList              * p_props;
    size_t                   m_dir;
    double                   m_B, m_b1, m_b2, m_E[2], m_weight;
    double                   m_shat, m_qt2, m_ylimits[2], m_y[2][2];
    double                   m_qt2prev[2];
    ATOOLS::Vec4D            m_eqt, m_qt, m_qtprev[2];

    bool                     m_test;
    std::map<std::string,ATOOLS::Histogram * > m_histos;


    virtual void   CalculateWeight() = 0;
    virtual void   InitLadder(const ATOOLS::Vec4D & pos);
    virtual void   ResetFSFlavours();
    virtual void   ConstructISKinematics();
    virtual void   ConstructSimpleLadder();
    virtual double RescaleLadder(Ladder * ladder,const ATOOLS::Vec4D & P_in);
    virtual double TotalReggeWeight(Ladder * ladder);
    virtual double AlphaSWeight(const double & kt2);
    virtual double ReggeWeight(const double & qt2, const double & y1,
			       const double y2);
    virtual double LDCWeight(const double & qt2, const double & qt2prev);
    virtual void   MakeTransverseUnitVector();
    inline virtual const double AlphaS(const double q2) const {
      return (*p_alphaS)(q2);
    }
    inline virtual const double AlphaSMax() const {
      return p_alphaS->MaxValue();
    }
 public:
    Ladder_Generator_Base();
    virtual ~Ladder_Generator_Base();

    virtual void Initialise(Remnant_Handler * remnants);
    virtual void Test();
    virtual Ladder * operator()(const ATOOLS::Vec4D & pos) = 0;
    virtual void QuarkReplace();
    inline virtual void FixLadderType() {
      //p_ladder->FixType(-m_Ymax + m_deltaY, m_Ymax - m_deltaY);
      m_Yrange = 1.;
      //m_Yrange = 0.1*m_beamY;
      //m_Yrange = 0.01*m_beamY;
      msg_Out() << "----> in Ladders, m_beamY = " << m_beamY << " and m_Yrange = " << m_Yrange << std::endl;
      p_ladder->FixType(-m_beamY+m_Yrange, m_beamY-m_Yrange);
    }
    inline virtual void InitCollision(Omega_ik * eikonal,const double & B) {
      p_eikonal  = eikonal; m_B = B; m_density.SetEikonal(p_eikonal);
    }
    inline virtual void SetImpactParameters(const double & b1,
					    const double & b2) {
      m_b1 = b1; m_b2 = b2;
    }
    inline virtual void SetMaximalScale(const double & E1,const double & E2) {
      m_E[0] = E1; m_E[1] = E2;
    }
    inline virtual const double Ymax()         const { return m_Ymax; }
    inline virtual const double KT2min()       const { return m_kt2minShower; }
    inline virtual const double Weight()       const { return m_weight; } 
  };
}
#endif
