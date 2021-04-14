#ifndef SHRIMPS_Ladders_Ladder_Generator_QT_H
#define SHRIMPS_Ladders_Ladder_Generator_QT_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_QT : public Ladder_Generator_Base {
  private:
    Sigma_Partonic  m_partonic;
    double          m_Ebeam[2];
    double          m_S, m_maxsigma, m_qt2minFF, m_shatmin, m_shat, m_seff;
    bool            m_fixflavour;
    ATOOLS::Flavour m_flavs[2];
    ATOOLS::Vec4D   m_qT, m_q[2], m_k[2];
    
    bool   FixInitialPartons();
    bool   MakeTrialLadder();
    bool   TrialEmission(size_t dir);
    bool   LastEmissions();
    void   AddEmission(size_t dir,TPropList::iterator & pit);
    double AbsorptionWeight(const ATOOLS::Vec4D & k,const double & y);
    bool   FixSimpleKinematics();
    void   SelectPropagatorColours();
    void   CalculateWeight();
    
    inline double QT2Min(size_t dir=2) {
      double qt2min = m_qt2min;
      for (size_t beam=0;beam<2;beam++) {
	if ((dir==2 || dir==beam) && (ATOOLS::dabs(m_y[beam][1])>m_Ymax)) qt2min = m_qt2minFF;
      }
      return qt2min;
    }
    inline double QT2Max() {
      return m_seff/ATOOLS::sqr(cosh(ATOOLS::Max(ATOOLS::dabs(m_y[0][1]),
						 ATOOLS::dabs(m_y[1][1]))));
    }
    inline ATOOLS::Vec4D MakeFSMomentum(size_t dir) {
      double kt = sqrt((m_q[dir]+(dir==0 ? -m_qT : m_qT)).PPerp2());
      return ( kt * ATOOLS::Vec4D(cosh(m_y[dir][1]),0.,0.,sinh(m_y[dir][1])) +
	       (m_q[dir] + (dir==0 ? -m_qT : m_qT)).Perp());
    }
    inline ATOOLS::Vec4D MakePropMomentum(const double & qt2min,const double & qt2max,
					  Form_Factor * ff=NULL) {
      m_qt2 = ATOOLS::Max(0.,(ff ?
			      ff->SelectQT2(qt2max,qt2min) :
			      qt2min * pow(qt2max/qt2min, ATOOLS::ran->Get())));
      return sqrt(m_qt2) * m_eqt;
    }
  public:
    Ladder_Generator_QT();
    ~Ladder_Generator_QT() {}

    void Initialise() { m_partonic.Initialise(); }
    Ladder * operator()(const ATOOLS::Vec4D & pos);
    const  double Weight() const { return m_weight; }    
    inline void SetRemnantHandler(Remnant_Handler * remnants) {
      for (size_t beam=0;beam<2;beam++) {
	p_pdf[beam] = remnants->GetPDF(beam);
	m_partonic.SetPDF(beam,p_pdf[beam]);
      }
    }
  };
}
#endif
