#ifndef Helicity_h
#define Helicity_h

#include "Flavour.H"

namespace AMEGIC {

  class Pol_Info;

  //! A list of signs, i.e. helicities.(ready)
  struct Sign_List {
    //! The list og signs
    int *s;
    //! Switching certain combinations on/off.
    int on;
    //! The multiplicity of the helicity combination, used when reducing the list to independent helicities. 
    int    multi;
    double polfactor;
    int    partner;

    Sign_List(): s(NULL),on(1),multi(1),polfactor(1.),partner(-1) {}
    ~Sign_List() {
      if (s) delete [] s;
    }
  };
  
  //! Determines all helicity combinations. (ready)
  class Helicity {
    //! The list of helicity combinations.
    Sign_List* p_slist;    
    //! The number of helicities
    int m_nsign;
    //! List of helicity combinations, neccessary to calculate physical fermion helicities
    Sign_List* p_fslist;
    int m_fnsign;

    int m_fermion_hels;

    //! type of polarisation: 'c'=circular, 'l'=linear, 'h'=helicity state and '+'/'-'  for polarised particles
    char   * p_pol_types;
    //! angle for linear polarisation
    double * p_angles;
  public:
    //! The constructor determines all helicity combinations
    /*!
      The helicity combinations are determined using a loop over loop technique.
      Having at hand all possible helicity combinations, the second helicity for 
      scalar particles (which have no helicity at all) and massive vector bosons
      (where the sum over the polarizations is carried out differently to the 
      massless vector bosons).
    */
    Helicity(int,int,ATOOLS::Flavour*,Pol_Info*);
    //! Returns the maximum number of helicities.
    int MaxHel() {return m_nsign;}
    //! Returns the maximum number of helicities for calculation of physical fermion helicities.
    int MaxFHel() {return m_fnsign;}
    //! Returns a certain helicity combination.
    int* operator[] (int i) {return p_fslist[i].s;}
    //! Switches a certain helicity combination off.
    void SwitchOff(int i) {p_slist[i].on = 0;}
    void SetPartner(int i,int j) {p_slist[j].partner = i;}
    inline const int&  Partner(const int i) const {return p_slist[i].partner;}
    //! Returns the status of a certain helicity combination.
    int  On(int i) {return p_slist[i].on;}
    //! Increases the multiplicity of a combination by one.
    void IncMultiplicity(int i) {(p_slist[i].multi)++;}  
    //! Returns the multiplicity of a certain helicity combination.
    int  Multiplicity(int i) {return p_slist[i].multi;}  
    //! Polarization factors
    double PolarizationFactor(int i) {return p_slist[i].polfactor;}  
    //! Returns true when physical fermion helicities should be calculated
    int IsCalcFerm() {return m_fermion_hels;}

    char   * PolTypes()  { return p_pol_types; }
    double * PolAngles() { return p_angles; }
    bool IsContrib(int,int*,int);

    ~Helicity();
  };
}
#endif
