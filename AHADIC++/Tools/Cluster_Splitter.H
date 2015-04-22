#ifndef AHADIC_Tools_Cluster_Splitter_H
#define AHADIC_Tools_Cluster_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"
#include "AHADIC++/Tools/Cluster.H"

namespace AHADIC {
  class Cluster_Splitter : public Splitter_Base {
  private:
    size_t m_nmax;
    double m_etax, m_etax_lead, m_etay, m_etay_lead;
    bool   m_added;

    size_t m_pairs;
    double m_alphanew,m_betanew;
    ATOOLS::Vec4D m_cms;

    std::list<PoppedPair *>::iterator m_popsplit, m_popspect;
    size_t m_popspliti, m_popspecti;
    Proto_Particle * p_trip, * p_anti;

    long int m_fourmom_err, m_kt2_err, m_LC_err, m_flav_err;

    size_t SelectNumberOfPairs(const size_t & nmax);
    bool   SelectSplitter(Proto_Particle * part1,Proto_Particle * part2);
    std::pair<double,double> FixExponents();
    bool   ConstructSystem(Cluster * cluster);
    void   ConstructKinematics(const double & etay,const double & etax);
    bool   AcceptSystem(const double & pt2max);
    bool   PoppedMassPossible(const double & m2);
    void   MakeKinematics();
    void   MakePairKinematics(PoppedPair * pop,
			      ATOOLS::Vec4D & test,ATOOLS::Vec4D & test1);
    void   MakeSplitterAndSpectatorMoms(ATOOLS::Vec4D & test,
					ATOOLS::Vec4D & test1);
    void   MakeClusters(Cluster * cluster);
    void   SelectPartners();
    bool   Rearrange();
    void   MakeSplitterAndSpectatorClusters(Cluster * cluster);
    void   MakeOtherClusters(Cluster * cluster);
    bool   EnforceMomentum(Cluster * cluster);
  public:
    Cluster_Splitter();
    ~Cluster_Splitter();
    
    bool operator()(Cluster * cluster);
  };
}

#endif
