#ifndef ANALYSIS_Triggers_MySISCone_H
#define ANALYSIS_Triggers_MySISCone_H

#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Histogram.H"
#include <list>
#include <vector>

namespace siscone {class Csiscone;}
namespace ANALYSIS {
  class Primitive_Analysis;

  class SISCone : public Jet_Algorithm_Base {

    
    //   temporary variable to store the jet momenta
    ATOOLS::Particle_List       *p_jets;
    const ATOOLS::Particle_List *p_orig;
    std::vector<double> * p_kts;

    double m_f;

    siscone::Csiscone* p_ConeCluster;

    ATOOLS::Flavour GetBFlavour(const ATOOLS::Vec4D & mom);
    void AddToKtlist(double );
    void AddToJetlist(const ATOOLS::Vec4D &);

  public:
    SISCone(ATOOLS::Particle_Qualifier_Base * const qualifier,double f=.5);
    ~SISCone();

    bool ConstructJets(const ATOOLS::Particle_List * pl = 0,ATOOLS::Particle_List * jets = 0,
		       std::vector<double> * kt2 = 0,double R=.7);
    int    NumberOfJets();

  };



  inline int    SISCone::NumberOfJets() { return p_jets->size(); }

  //
}

#endif
