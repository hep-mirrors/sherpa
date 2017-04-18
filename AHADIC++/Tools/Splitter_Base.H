#ifndef AHADIC_Tools_Splitter_Base_H
#define AHADIC_Tools_Splitter_Base_H

#include "AHADIC++/Tools/Flavour_Selector.H"
#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Z_Selector.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Math/Poincare.H"

namespace AHADIC {
  class Splitter_Base {
  protected:
    std::list<Cluster *> * p_cluster_list; 
    Soft_Cluster_Handler * p_softclusters;

    Single_Transitions   * p_singletransitions;
    Double_Transitions   * p_doubletransitions;
    Constituents         * p_constituents;
    Flavour_Selector       m_flavourselector;
    KT_Selector            m_ktselector;
    Z_Selector             m_zselector;
    
    Proto_Particle * p_part1, * p_part2, * p_part3;
    bool   m_barrd;
    double m_minmass;
    double m_Q2, m_E, m_mass1, m_mass2, m_mass3, m_Emax;
    double m_m12, m_m22, m_m32, m_alpha, m_beta;

    Flavour_Pair    m_flavs1, m_flavs2;
    double          m_popped_mass, m_popped_mass2;
    ATOOLS::Flavour m_newflav1, m_newflav2;
    double          m_minQ_1, m_minQ_2;
    double          m_ktmax, m_kt, m_kt2, m_z1, m_z2, m_phi;
    ATOOLS::Vec4D   m_ktvec;
    
    ATOOLS::Poincare m_boost, m_rotat;

    virtual bool InitSplitting(Proto_Particle * part1,Proto_Particle * part2,
			       Proto_Particle * part3);
    virtual void FillMasses();
    virtual void ConstructLightCone(const double & kt2=0.);
    virtual void ConstructPoincare();

    virtual bool MakeSplitting();
    virtual void PopFlavours();
    virtual void DetermineMinimalMasses();
    virtual bool MakeKinematics();
    virtual void MakeTransverseMomentum();
    virtual bool MakeLongitudinalMomenta() = 0;
    virtual void FillParticlesInLists()    = 0;
    virtual bool CheckKinematics()         = 0;
  public:
    Splitter_Base(std::list<Cluster *> * cluster_list,
		  Soft_Cluster_Handler * softclusters);
    virtual void Init();
    
    virtual bool operator()(Proto_Particle * part1,Proto_Particle * part2,
			    Proto_Particle * part3=0);
  };
  
}

#endif
