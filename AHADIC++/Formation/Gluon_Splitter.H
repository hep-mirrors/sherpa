#ifndef AHADIC_Formation_Gluon_Splitter_H
#define AHADIC_Formation_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"


namespace AHADIC {
  class Gluon_Splitter : public Splitter_Base {
  private:
    int    m_mode;
    double m_alpha, m_arg, m_x, m_y;
    
    long int  m_kin_fails{ 0 };
    
    bool      MakeLongitudinalMomenta();
    void      CalculateLimits();
    bool      CalculateXY();

    bool      FillParticlesInLists();
    bool      CheckKinematics();
    Cluster * MakeCluster();
    void      UpdateSpectator();
  public:
    Gluon_Splitter(std::list<Cluster *> * cluster_list,
		   Soft_Cluster_Handler * softclusters) :
      Splitter_Base(cluster_list,softclusters) {}
    ~Gluon_Splitter();

    void   Init(const bool & isgluon=true);
    double WeightFunction(const double & z,
			  const double & zmin=0.,const double & zmax=1.,
			  const unsigned int & cnt=0);
  };
}

#endif
