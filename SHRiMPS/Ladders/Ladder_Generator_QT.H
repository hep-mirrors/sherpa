#ifndef SHRIMPS_Ladders_Ladder_Generator_QT_H
#define SHRIMPS_Ladders_Ladder_Generator_QT_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_QT : public Ladder_Generator_Base {
  private:
    double           m_sigmahat, m_seff;
    ATOOLS::Vec4D    m_q[2], m_k[2], m_qini[2];
    ATOOLS::Flavour  m_flavs[2];
    ATOOLS::Poincare m_cms, m_zax;
    double           m_shat, m_yhat;
    
    bool     FixInitialPartons();
    Ladder * InitializeRescatterLadder(Ladder_Particle * parts[2]);
    bool     MakeTrialLadder(Ladder * ladder);
    Ladder * MakeRescatterLadder(LadderMap::iterator lit[2],TPropList::iterator & pit1);
    void     MergeLadders(Ladder * ladder,Ladder * rescatter,
			  LadderMap::iterator lit[2],TPropList::iterator pit);
    bool     TrialEmission(Ladder * ladder,size_t dir,const double & yhat=0);
    bool     LastEmissions(Ladder * ladder);
    void     AddEmission(Ladder * ladder,size_t dir,TPropList::iterator & pit);
    double   AbsorptionWeight(const ATOOLS::Vec4D & k,const double & y);
    double   MEWeight();
    bool     FixSimpleKinematics();
    void     SelectPropagatorColours(Ladder * ladder);
    void     CalculateWeight();
    
    double   QT2Min(size_t dir=2);
    double   QT2Max();
    ATOOLS::Vec4D MakeFSMomentum(size_t dir);
    ATOOLS::Vec4D MakePropMomentum(const double & qt2min,const double & qt2max,
				   Form_Factor * ff=NULL,const bool isrescatter=false);
  public:
    Ladder_Generator_QT();
    ~Ladder_Generator_QT() {}

    Ladder * operator()(const ATOOLS::Vec4D & pos);
    void   AddRescatters(Ladder * ladder);
  };
}
#endif
