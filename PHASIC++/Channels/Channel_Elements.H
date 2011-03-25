#ifndef PHASIC_Channels_Channel_Elements_H
#define PHASIC_Channels_Channel_Elements_H

#include "ATOOLS/Math/Vector.H"
#include "PHASIC++/Channels/Channel_Basics.H"
#include "ATOOLS/Org/Info_Key.H"

namespace PHASIC {

  class Channel_Elements {
  private:
    
    Channel_Basics CB;
    void QcdAntenna(ATOOLS::Vec4D* &,int,double, double);
    void BasicAntenna(ATOOLS::Vec4D,ATOOLS::Vec4D,ATOOLS::Vec4D &, double*, double);
    void PermP(int,int* &);
    void Polytope(int,double* &);
    void CheckMasses(const double & s1,ATOOLS::Vec4D & p1,
		     const double & s2,ATOOLS::Vec4D & p2) const;
  public:

    double Isotropic2Weight(const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,double=-1.0,double=1.0);
    double Isotropic2Weight(ATOOLS::Vec4D&,ATOOLS::Vec4D&,double&,double&,double=-1.0,double=1.0);
    void Isotropic2Momenta(ATOOLS::Vec4D,double,double,
			   ATOOLS::Vec4D&,ATOOLS::Vec4D&,double,double,double=-1.0,double=1.0);
    void Isotropic2Grid(ATOOLS::Vec4D&,ATOOLS::Vec4D&,double&,double&,double=-1.0,double=1.0);

    double Anisotropic2Weight(double,double,double,
			      const ATOOLS::Vec4D&,const ATOOLS::Vec4D&);
    void Anisotropic2Momenta(ATOOLS::Vec4D,double,double,double,double,double,
			     ATOOLS::Vec4D&,ATOOLS::Vec4D&,double,double);

    double TChannelWeight(const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,
			  const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,
			  double,double,double,double,double,int);
    double TChannelWeight(const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,
			  const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,
			  double,double,double,double,double,int,double&,double&);
    int TChannelMomenta(ATOOLS::Vec4D,ATOOLS::Vec4D,ATOOLS::Vec4D&,ATOOLS::Vec4D&,
			double,double,double,double,double,double,double,int,double,double);
    void TChannelGrid(const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,
		      const ATOOLS::Vec4D&,const ATOOLS::Vec4D&,  
		      double,double,double,double,double&,double&);
    double BremsstrahlungWeight(double,double,double,
				const ATOOLS::Vec4D&,const ATOOLS::Vec4D&);
    void BremsstrahlungMomenta(ATOOLS::Vec4D&,const double,const double,const double,
			       const double,const double,const double,
			       ATOOLS::Vec4D &, ATOOLS::Vec4D &,const double,const double);

    double MasslessPropWeight(double,double,double,double);
    double MasslessPropWeight(double,double,double,double,double&);
    double MasslessPropMomenta(double,double,double,double);

    double AntennaWeight(double amin,double amax,const double a,double &ran);
    double AntennaMomenta(double amin,double amax,double ran);

    double MassivePropWeight(double,double,int,double,double,double);
    double MassivePropWeight(double,double,int,double,double,double,double&);
    double MassivePropMomenta(double,double,int,double,double,double);

    double QCDAPWeight(ATOOLS::Vec4D*,int,double);
    void QCDAPMomenta(ATOOLS::Vec4D*,ATOOLS::Vec4D,int,double);

    double LLPropWeight(double,double,double,double,double);
    double LLPropWeight(double,double,double,double,double,double&);
    double LLPropMomenta(double,double,double,double,double);

    double ThresholdMomenta(double,double,double,double);
    double ThresholdWeight(double,double,double,double);
    double ThresholdMomenta(double,double,double,double,double);
    double ThresholdWeight(double,double,double,double,double);
    double ThresholdWeight(double,double,double,double,double,double&);

    double DiceYUniform(const double tau,const ATOOLS::Double_Container &xinfo,
			const ATOOLS::Double_Container &yinfo,const double ran,const int mode) const;
    double WeightYUniform(const double tau,const ATOOLS::Double_Container &xinfo,
			  const ATOOLS::Double_Container &yinfo,double&,const int mode) const;
    double DiceYCentral(const double tau,const ATOOLS::Double_Container &xinfo,
			const ATOOLS::Double_Container &yinfo,const double ran,const int mode) const;
    double WeightYCentral(const double tau,const ATOOLS::Double_Container &xinfo,
			  const ATOOLS::Double_Container &yinfo,double&,const int mode) const;
    double DiceYForward(const double yexponent,const double tau,
			const ATOOLS::Double_Container &xinfo,
			const ATOOLS::Double_Container &yinfo,const double ran,const int mode) const;
    double WeightYForward(const double yexponent,const double tau,
			  const ATOOLS::Double_Container &xinfo,
			  const ATOOLS::Double_Container &yinfo,double&,const int mode) const;
    double DiceYBackward(const double yexponent,const double tau,
			 const ATOOLS::Double_Container &xinfo,
			 const ATOOLS::Double_Container &yinfo,const double ran,const int mode) const;
    double WeightYBackward(const double yexponent,const double tau,
			   const ATOOLS::Double_Container &xinfo,
			   const ATOOLS::Double_Container &yinfo,double&,const int mode) const;

  };// end of class Channel_Elements

  extern Channel_Elements CE;                                                  

  /*!
    Weights and momenta for isotropic two-body decays.
    The idea for the gbeneration of the momenta is to boost the
    decaying vector into its rest frame, select the orientation of
    the two decay products with given masses 4-pi isotropically and 
    to boost back. The weight is just a constant.
  */

  /*!
    Weights and momenta for anisotropic two-body decays as they occur
    in splittings involving one outgoing massless vector boson.
    The idea here is to use a peaked distribution for the polar angle
    of one decay product (not the "vector boson") along the original direction
    of the decaying particle. 
  */

  /*!
    Weights and momenta for a Tchannel process. Again, there is a strong peak in
    the forward direction in order to minimize the mass of the propagator.
  */

  /*! 
     Decay p -> q + p1, q is space-like with energy Eq given from outside
     cos(pq) is constrained by ctmin and ctmax. 
  */

  /*!
    Scale and weight for a massless propagator. It is distributed 
    according to a simple pole structure peaking at 0 (or the minimal
    scale that is kinematically allowed).
  */

  /*!
    Scale and weight for a massive propagator. It is distributed 
    according to a physical Breit-Wigner distribution with given
    mass and width.
  */

  /*!
    Scale and weight for a typical electron initial state radiation 
    event. It is distributed according to a simple pole structure peaking 
    at the maximal scale that is kinematically allowed.
  */

}// end of namespace PHASIC

#endif
