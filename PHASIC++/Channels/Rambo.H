#ifndef PHASIC_Channels_Rambo_H
#define PHASIC_Channels_Rambo_H

#include "PHASIC++/Channels/Single_Channel.H"

namespace PHASIC {
  class Rambo : public Single_Channel {
    //int      nin,nout;
    //double * ms;
    bool     massflag;
    bool     m_use_psmass;
    double   Z_N, * xm2, * p2, * E;
    short int itmax;
    double accu;
  public:
    Rambo(int,int,const ATOOLS::Flavour *,bool use_psmass=false);
    ~Rambo();
    void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);  
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *);
    void GenerateWeight(ATOOLS::Vec4D * p) { GenerateWeight(p,NULL); }  
    void GeneratePoint(ATOOLS::Vec4D * p)  { GeneratePoint(p,NULL); }
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    void MassiveWeight(ATOOLS::Vec4D *,double);
    void MassivePoint(ATOOLS::Vec4D *,double);
    std::string Name() { return std::string("Rambo"); }
    void UsePSMass(bool use_psmass) { m_use_psmass = use_psmass; }
    int CountResonances(ATOOLS::Flavour *&) { return 0; }
  };
  //! Flag to see wether massive or massless points are to be generated.
  //! Constants needed for the weight and the phase space distribution.
  //! Maximal number of iterations to generate the massive weight
  //! Accuracy when generating the massive weight.
  /*!
    Generates the weight. In case of massless Rambo this is given by
    a constant, since Rambo uses uniform distribution. For massive
    points, this constant is corrected via MassiveWeight.
    /sa MassiveWeight
  */
  /*!
    To generate a phase space point according to a uniform distribution.
    In case massive vectors are needed this point gets modified by 
    MassivePoint.
    /sa MassivePoint
  */
  /*!
    This method is only a dummy method meant to comply with the multi-channel.
    In fact it calls the GeneratePoint method without random numbers.
  */
  /*!
    To generate a massive weight. The idea is to invert the generation of
    massive vectors, i.e. we start with massive vectors and force them
    on a zero-mass shell by rescaling and boosting all momenta by common 
    factors, such that both four momentum conservation and mass shell 
    conditions are fulfilled.
  */
  /*!
    To generate a massive point. The idea is to start from a massless point
    given by ordinary RAMBO::GeneratePoint and then to force them
    on their mass shell by rescaling and boosting all momenta by common 
    factors, such that both four momentum conservation and mass shell 
    conditions are fulfilled.
  */
}

#endif
