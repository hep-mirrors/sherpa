#ifndef Helicity_h
#define Helicity_h

#include "Flavour.H"
//#include "Pol_Info.H"

namespace AMEGIC {

  class Pol_Info;
  //! A list of signs, i.e. helicities.(ready)
  struct signlist {
    //! The list og signs
    int *s;
    //! Switching certain combinations on/off.
    int on;
    //! The multiplicity of the helicity combination, used when reducing the list to independent helicities. 
    int    Mult;
    double polfactor;
    int    partner;
  };
  
  //! Determines all helicity combinations. (ready)
  class Helicity {
    //! The list of helicity combinations.
    signlist* Slist;    
    //! The number of helicities
    int nsign;
    //! List of helicity combinations, neccessary to calculate physical fermion helicities
    signlist* FSlist;
    int fnsign;

    int fermion_hels;

  public:
    //! type of polarisation: 'c'=circular, 'l'=linear, 'h'=helicity state
    char *p_type;
    //! angle for linear polarisation
    double *angle;
    //! The constructor determines all helicity combinations
    /*!
      The helicity combinations are determined using a loop over loop technique.
      Having at hand all possible helicity combinations, the second helicity for 
      scalar particles (which have no helicity at all) and massive vector bosons
      (where the sum over the polarizations is carried out differently to the 
      massless vector bosons).
    */
    Helicity(int,int,APHYTOOLS::Flavour*,Pol_Info*);
    //! Returns the maximum number of helicities.
    int Max_Hel() {return nsign;}
    //! Returns the maximum number of helicities for calculation of physical fermion helicities.
    int Max_FHel() {return fnsign;}
    //! Returns a certain helicity combination.
    int* operator[] (int i) {return FSlist[i].s;}
    //! Switches a certain helicity combination off.
    void switch_off(int i) {Slist[i].on = 0;}
    void SetPartner(int i,int j) {Slist[j].partner = i;}
    inline const int&  Partner(const int i) const {return Slist[i].partner;}
    //! Returns the status of a certain helicity combination.
    int  On(int i) {return Slist[i].on;}
    //! Increases the multiplicity of a combination by one.
    void Inc_Mult(int i) {(Slist[i].Mult)++;}  
    //! Returns the multiplicity of a certain helicity combination.
    int  Multiplicity(int i) {return Slist[i].Mult;}  
    //! Polarization factors
    double PolarizationFactor(int i) {return Slist[i].polfactor;}  
    //! Returns true when physical fermion helicities should be calculated
    int iscalcferm() {return fermion_hels;}
  };
}
#endif
