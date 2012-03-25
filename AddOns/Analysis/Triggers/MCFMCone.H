#ifndef ANALYSIS_Triggers_MCFMCone_H
#define ANALYSIS_Triggers_MCFMCone_H

#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Histogram.H"
#include <list>
#include <vector>

namespace ATOOLS {class ConeMCFM;}
namespace ANALYSIS {
  class Primitive_Analysis;

  class MCFMCone : public Jet_Algorithm_Base {

    
    //   temporary variable to store the jet momenta
    ATOOLS::Particle_List       * p_jets;
    std::vector<double> * p_kts;

    double m_rsep;

    ATOOLS::ConeMCFM* p_ConeCluster;

    void AddToKtlist(double );
    void AddToJetlist(const ATOOLS::Vec4D &);

  public:
    MCFMCone(ATOOLS::Particle_Qualifier_Base * const qualifier,double f=.5);
    ~MCFMCone();

    bool ConstructJets(const ATOOLS::Particle_List * pl = 0,ATOOLS::Particle_List * jets = 0,
		       std::vector<double> * kt2 = 0,double R=.7);
    int    NumberOfJets();

  };



  inline int    MCFMCone::NumberOfJets() { return p_jets->size(); }

  //
}

#endif
