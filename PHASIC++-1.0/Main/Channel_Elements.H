#ifndef Channel_Elements_h
#define Channel_Elements_h

#include "Vector.H"
#include "Channel_Basics.H"

namespace PHASIC {
  class  Channel_Elements;
  extern Channel_Elements CE;                                                  

  class Channel_Elements {
    Channel_Basics CB;
  public:
    /*!
      Weights and momenta for isotropic two-body decays.
      The idea for the gbeneration of the momenta is to boost the
      decaying vector into its rest frame, select the orientation of
      the two decay products with given masses 4-pi isotropically and 
      to boost back. The weight is just a constant.
    */
    double Isotropic2Weight(const AMATOOLS::vec4d&,const AMATOOLS::vec4d&);
    void Isotropic2Momenta(AMATOOLS::vec4d,double&,double&,
			   AMATOOLS::vec4d&,AMATOOLS::vec4d&,double,double);
    /*!
      Weights and momenta for anisotropic two-body decays as they occur
      in splittings involving one outgoing massless vector boson.
      The idea here is to use a peaked distribution for the polar angle
      of one decay product (not the "vector boson") along the original direction
      of the decaying particle. 
    */
    double Anisotropic2Weight(double,double,double,
			      const AMATOOLS::vec4d&,const AMATOOLS::vec4d&);
    void Anisotropic2Momenta(AMATOOLS::vec4d,double,double,double,double,double,
			     AMATOOLS::vec4d&,AMATOOLS::vec4d&,double,double);
    /*!
      Weights and momenta for a Tchannel process. Again, there is a strong peak in
      the forward direction in order to minimize the mass of the propagator.
    */
    double TChannelWeight(const AMATOOLS::vec4d&,const AMATOOLS::vec4d&,
			  const AMATOOLS::vec4d&,const AMATOOLS::vec4d&,
			  double,double,double,double,double,int);
    int TChannelMomenta(AMATOOLS::vec4d,AMATOOLS::vec4d,AMATOOLS::vec4d&,AMATOOLS::vec4d&,
			double&,double&,double,double,double,double,double,int,double,double);
    /* Decay p -> q + p1, q is space-like with energy Eq given from outside
       cos(pq) is constrianed by ctmin and ctmax. */
    double BremsstrahlungWeight(double,double,double,
				const AMATOOLS::vec4d&,const AMATOOLS::vec4d&);
    void BremsstrahlungMomenta(AMATOOLS::vec4d&,const double,const double,const double,
			       const double,const double,const double,
			       AMATOOLS::vec4d &, AMATOOLS::vec4d &,const double,const double);
    /*!
      Scale and weight for a massless propagator. It is distributed 
      according to a simple pole structure peaking at 0 (or the minimal
      scale that is kinematically allowed).
    */
    double MasslessPropWeight(double,double,double,const double&);
    double MasslessPropMomenta(double,double,double,double);
    /*!
      Scale and weight for a massive propagator. It is distributed 
      according to a physical Breit-Wigner distribution with given
      mass and width.
    */
    double MassivePropWeight(double,double,int,double,double,double);
    double MassivePropMomenta(double,double,int,double,double,double);
    /*!
      Scale and weight for a typical electron initial state radiation 
      event. It is distributed according to a simple pole structure peaking 
      at the maximal scale that is kinematically allowed.
    */
    double LLPropWeight(double,double,double,double,const double&);
    double LLPropMomenta(double,double,double,double,double);
    double ThresMomenta(double,double,double,double,double);
    double ThresWeight(double,double,double,double,double);
    /*
      double Isotropic3Weight(const AMATOOLS::vec4d&,const AMATOOLS::vec4d&,const AMATOOLS::vec4d&);
      void Isotropic3Momenta(const AMATOOLS::vec4d&,const double,
      AMATOOLS::vec4d&,AMATOOLS::vec4d&,AMATOOLS::vec4d&,
      const double,const double,const double);
    */
    double DiceYUniform(double, double *, double *, int, double);
    double WeightYUniform(double, double *, double *, int, double);
    double DiceYCentral(double, double *, double *, int, double);
    double WeightYCentral(double, double *, double *, int, double);
    double DiceYForward(double, double *, double *, double, int, double);
    double WeightYForward(double, double *, double *, double, int, double);
    double DiceYBackward(double, double *, double *, double, int, double);
    double WeightYBackward(double, double *, double *, double, int, double);
  };
}
#endif
