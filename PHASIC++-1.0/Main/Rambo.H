#ifndef Rambo_h
#define Rambo_h

#include "Single_Channel.H"

namespace PHASIC {
  class Rambo : public Single_Channel {
    int      nin,nout;
    double * ms;
    //! Flag to see wether massive or massless points are to be generated.
    bool     massflag;
    //! Constants needed for the weight and the phase space distribution.
    double   Z_N, * xm2, * p2, * E;
    //! Maximal number of iterations to generate the massive weight
    short int itmax;
    //! Accuracy when generating the massive weight.
    double accu;
  public:
    Rambo(int,int,APHYTOOLS::Flavour *);
    ~Rambo();
    /*!
      Generates the weight. In case of massless Rambo this is given by
      a constant, since Rambo uses uniform distribution. For massive
      points, this constant is corrected via MassiveWeight.
      /sa MassiveWeight
    */
    void GenerateWeight(AMATOOLS::vec4d *,APHYTOOLS::Cut_Data *);  
    /*!
      To generate a phase space point according to a uniform distribution.
      In case massive vectors are needed this point gets modified by 
      MassivePoint.
      /sa MassivePoint
    */
    void GeneratePoint(AMATOOLS::vec4d *,APHYTOOLS::Cut_Data *);
    /*!
      This method is only a dummy method meant to comply with the multi-channel.
      In fact it calls the GeneratePoint method without random numbers.
    */
    void GeneratePoint(AMATOOLS::vec4d *,APHYTOOLS::Cut_Data *,double *);
    /*!
      To generate a massive weight. The idea is to invert the generation of
      massive vectors, i.e. we start with massive vectors and force them
      on a zero-mass shell by rescaling and boosting all momenta by common 
      factors, such that both four momentum conservation and mass shell 
      conditions are fulfilled.
    */
    void MassiveWeight(AMATOOLS::vec4d *,double);
    /*!
      To generate a massive point. The idea is to start from a massless point
      given by ordinary RAMBO::GeneratePoint and then to force them
      on their mass shell by rescaling and boosting all momenta by common 
      factors, such that both four momentum conservation and mass shell 
      conditions are fulfilled.
    */
    void MassivePoint(AMATOOLS::vec4d *,double);
    std::string Name() { return std::string("Rambo"); }
    int CountResonances(APHYTOOLS::Flavour *&) { return 0; }
  };
}

#endif
