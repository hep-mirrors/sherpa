#ifndef AHADIC_Tools_Soft_Cluster_Handler_H
#define AHADIC_Tools_Soft_Cluster_Handler_H

#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Single_Transitions.H"
#include "AHADIC++/Tools/Double_Transitions.H"
#include "AHADIC++/Tools/Cluster.H"

namespace AHADIC {
  class Soft_Cluster_Handler {
  private:
    std::list<Proto_Particle *> * p_hadrons; 
    Constituents                * p_constituents;
    Single_Transitions          * p_singletransitions;
    Double_Transitions          * p_doubletransitions;
    KT_Selector                   m_ktselector;
    double m_light;
    double m_piphoton_threshold, m_dipion_threshold, m_open_threshold;

    Cluster    * p_cluster;
    double       m_mass2, m_mass;
    Flavour_Pair m_flavs;

    ATOOLS::Vec4D   m_moms[2];
    ATOOLS::Flavour m_hads[2];
    
    bool   FillFlavours(Cluster * cluster);
    bool   CheckOutsideRange();
    bool   Decay();
    double PhaseSpace(const double & m2,const double & m3);
    double DecayWeight();
    double RadiationWeight();
    double Annihilation();
    bool   FixKinematics();
    bool   DiQuarkToQuarks(const ATOOLS::Flavour & di,
			   ATOOLS::Flavour & q1,ATOOLS::Flavour & q2);
    double DefineHadronsInAnnihilation(const Flavour_Pair & one,
				       const Flavour_Pair & two);

 public:
    Soft_Cluster_Handler(std::list<Proto_Particle *> * hadrons);
    ~Soft_Cluster_Handler();

    void Init();
    bool MustPromptDecay(Cluster * cluster);
    bool Treat(Cluster * cluster,bool force=false);
    bool TreatTwoGluons(Cluster * cluster);
    bool RadiativeDecay(Cluster * cluster);

    double MinSingleMass(const ATOOLS::Flavour & fl1,
			 const ATOOLS::Flavour & fl2);
    ATOOLS::Flavour LowestTransition(const ATOOLS::Flavour & fl1,
				     const ATOOLS::Flavour & fl2);
    double MinDoubleMass(const ATOOLS::Flavour & fl1,
			 const ATOOLS::Flavour & fl2);

    void PrintHadrons();
    std::list<Proto_Particle *> * GetHadrons() const { return p_hadrons; } 
  };
}
#endif
