#ifndef AHADIC_Formation_Gluon_Splitter_H
#define AHADIC_Formation_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"


namespace AHADIC {
  class Gluon_Splitter : public Splitter_Base {
  private:
    int    m_mode;
    double m_alpha, m_arg, m_x, m_y, m_lastmass;
    bool   m_lastB, m_lastC;
    
    long int  m_kin_fails{ 0 };
    
    bool      MakeLongitudinalMomenta();
    void      CalculateLimits();
    bool      CalculateXY();

    bool      FillParticlesInLists();
    bool      CheckKinematics();
    bool      CheckConstituentKinematics(const ATOOLS::Vec4D & newmom11,
					 const ATOOLS::Vec4D & newmom12);
    Cluster * MakeCluster();
    void      ReplaceClusterWithHadron(const ATOOLS::Flavour & fl,ATOOLS::Vec4D & mom);
    void      UpdateSpectator(const ATOOLS::Vec4D & clumom);
  public:
    Gluon_Splitter(std::list<Cluster *> * cluster_list,
		   Soft_Cluster_Handler * softclusters) :
      Splitter_Base(cluster_list,softclusters) {}
    ~Gluon_Splitter();

    void   Init(const bool & isgluon=true);
    double WeightFunction(const double & z,
			  const double & zmin=0.,const double & zmax=1.,
			  const unsigned int & cnt=0);
    inline void GetLast(double & lastmass, bool & isB, bool & isC) {
      lastmass = m_lastmass; isB = m_lastB; isC = m_lastC;
    }
  };
}

#endif
