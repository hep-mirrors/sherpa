#ifndef AHADIC_Formations_Singlet_Checker_H
#define AHADIC_Formations_Singlet_Checker_H

#include "AHADIC++/Formation/Trivial_Splitter.H"
#include "AHADIC++/Tools/Singlet_Tools.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include <list>
#include <map>

namespace AHADIC {
  class Singlet_Checker : public Singlet_Tools {
  private:
    std::list<Singlet *>                    * p_singlets;
    std::list<std::list<Singlet *>::iterator> m_badones;
    std::map<Singlet *,ATOOLS::Flavour>       m_transitions;
    std::list<Proto_Particle *>             * p_hadrons;
    
    Soft_Cluster_Handler * p_softclusters;
    Trivial_Splitter       m_splitter;
    
    Proto_Particle * p_part1, * p_part2;
    double           m_minSmass2;
    bool             m_isring;

    ATOOLS::Vec4D    m_singletmom;
    double           m_targetmass;
    Proto_Particle * p_recoiler;

    bool      CheckSinglet();
    bool      FusePartonsInLowMassSinglet();
    bool      DealWithProblematicSinglets();
    void      SortProblematicSinglets();
    bool      FindOtherSingletToTransit();
    bool      FindRecoilerForTransit();
    bool      TestSingletForRecoiler();
    bool      TestForOtherRecoilers();
    bool      TestForRecoilerBeams();
    bool      TestRecoiler(Proto_Particle * part);
    bool      TransitProblematicSinglets();
    bool      TransitProblematicSingletWithRecoiler();
    void      ForcedDecays();
    bool      ForcedDecayOfTwoPartonSinglet();
    bool      ExtractAndCheckFlavours();
    bool      TwoGluonSingletToHadrons();
    bool      TwoQuarkSingletToHadrons();
 public:
    Singlet_Checker(std::list<Singlet *> * p_singlets,
		    Soft_Cluster_Handler * softclusters);
    ~Singlet_Checker();
    void Init();
    
    bool operator()();
  };
}

#endif
