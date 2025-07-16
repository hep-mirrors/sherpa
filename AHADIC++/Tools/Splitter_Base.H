#ifndef AHADIC_Tools_Splitter_Base_H
#define AHADIC_Tools_Splitter_Base_H

#include "AHADIC++/Tools/Flavour_Selector.H"
#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Z_Selector.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

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

    bool   m_ktorder;
    size_t m_attempts;
    
    Proto_Particle * p_part[3];
    ATOOLS::Vec4D    m_Qvec;
    bool   m_barrd;
    double m_minmass;
    double m_Q2, m_Q, m_E, m_mass[3], m_m2[3], m_Emax;
    double m_lc[2];

    Flavour_Pair    m_flavs1, m_flavs2;
    double          m_popped_mass, m_popped_mass2, m_msum[2], m_msum2[2];
    ATOOLS::Flavour m_newflav[2];
    double          m_minQ[2], m_minQ2[2], m_mdec[2], m_mdec2[2];
    
    double          m_zmin[2], m_zmax[2], m_z[2];
    double          m_ktmax, m_kt2max, m_ktfac, m_kt, m_kt2, m_phi;
    ATOOLS::Vec4D   m_ktvec;
    
    ATOOLS::Poincare m_boost, m_rotat;

    bool m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

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
    virtual bool FillParticlesInLists()    = 0;
    virtual bool CheckKinematics()         = 0;
  public:
    Splitter_Base(std::list<Cluster *> * cluster_list,
		  Soft_Cluster_Handler * softclusters);
    ~Splitter_Base();
    virtual void Init(const bool & isgluon);
    
    virtual bool operator()(Proto_Particle * part1,Proto_Particle * part2,
			    Proto_Particle * part3=0);
    virtual double WeightFunction(const double & z,
				  const double & zmin,const double & zmax,
				  const unsigned int & cnt)=0;
  };
  
}

#endif
