#ifndef PHASIC_Channels_RamboKK_H
#define PHASIC_Channels_RamboKK_H

#include "PHASIC++/Channels/Single_Channel.H"

namespace PHASIC {
  class RamboKK : public Single_Channel {
    //int      nin,nout;
    //double * ms;
    //! Flag to see wether massive or massless points are to be generated.
    bool     massflag;
    //! Constants needed for the weight and the phase space distribution.
    double   Z_N, * xm2, * p2, * E;
    //! Maximal number of iterations to generate the massive weight
    short int itmax;
    //! Accuracy when generating the massive weight.
    double accu;
    //! Constants needed to generate the mass of the KK particle 
    int     kkp,ed;
    double  r2,maxm2,maxn,mpss,gn,prevET,gam;
    /*!
      Determines the mass for the KK particle. This has to be done every time
      before a phase space point is generated.

      The distribution of the mass states is
      \f[
      m^{2}_{\vec{n}}=\frac{4\pi^{2}\vec{n}^{2}}{R^{2}},
      \f]
      where \f$\vec{n}\f$ is a $n$-dimensional vector of integers and with $R$ 
      contributes to the size of the extra dimensions, is given by
      \f[
      R^{n}=\frac{(4\pi)^{n/2}\Gamma(n/2)}{2M_s^{n+2}}.
      \f]
    */
    void Set_KKmass();
  public:
    RamboKK(int,int,const ATOOLS::Flavour *);
    ~RamboKK();
    /*!
      Generates the weight. In case of massless Rambo this is given by
      a constant, since Rambo uses uniform distribution. For massive
      points, this constant is corrected via MassiveWeight.
      /sa MassiveWeight

      In attendance of a KK particle the weight has to be multiplied by an 
      extra volume factor of
      \f[
      \frac{1}{n}\frac{\hat{m}^{n}}{M_s^{n+2}G_N},
      \f]
      due to the multiplicity of the mass states
    */
    void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);  
    /*!
      To generate a phase space point according to a uniform distribution.
      In case massive vectors are needed this point gets modified by 
      MassivePoint.
      /sa MassivePoint
    */
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *);
    /*!
      This method is only a dummy method meant to comply with the multi-channel.
      In fact it calls the GeneratePoint method without random numbers.
    */
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    /*!
      To generate a massive weight. The idea is to invert the generation of
      massive vectors, i.e. we start with massive vectors and force them
      on a zero-mass shell by rescaling and boosting all momenta by common 
      factors, such that both four momentum conservation and mass shell 
      conditions are fulfilled.
    */
    void MassiveWeight(ATOOLS::Vec4D *,double);
    /*!
      To generate a massive point. The idea is to start from a massless point
      given by ordinary RAMBO::GeneratePoint and then to force them
      on their mass shell by rescaling and boosting all momenta by common 
      factors, such that both four momentum conservation and mass shell 
      conditions are fulfilled.
    */
    void MassivePoint(ATOOLS::Vec4D *,double);
    std::string Name() { return std::string("RamboKK"); }
    int CountResonances(ATOOLS::Flavour *&) { return 0; }
  };
}

#endif
