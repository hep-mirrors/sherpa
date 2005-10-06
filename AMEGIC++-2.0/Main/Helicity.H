#ifndef Helicity_h
#define Helicity_h

#include "Flavour.H"
#include <iostream>
#include "MyComplex.H"
#include "Basic_Sfuncs.H"
#include "Matrix.H"

namespace AMEGIC {

  class Pol_Info;

  //! A list of signs, i.e. helicities.(ready)
  struct Sign_List {
    //! The list of signs
    int *s;
    //! Switching certain combinations on/off.
    int on;
    //! The multiplicity of the helicity combination, used when reducing the list to independent helicities. 
    int    multi;
    double polfactor;
    int    partner;

    Sign_List(): s(0),on(1),multi(1),polfactor(1.),partner(-1) {}
    ~Sign_List() {
      //      if (s) delete [] s;
    }
  };
  
 
class Helicity {
    std::vector<ATOOLS::Flavour> m_flavours; // list of the external particles
    std::vector<int> m_nPols;  // the number of possible polarisation states for each flavour

 
    Sign_List* p_slist;    
    size_t m_nsign;               // the number of polarisation combinations
  //    size_t m_nPhysHel;            // the relevant number of helicity combinations


  //    Sign_List* p_fslist;
  //    int m_fnsign;

    //! type of polarisation: 'c'=circular, 'l'=linear, 'h'=helicity state and '+'/'-'  for polarised particles
    char   * p_pol_types;
    //! angle for linear polarisation
    double * p_angles;

    bool m_allowTrafo;
    bool m_needsTrafo;            // is a trafo to the helicity base nessecary
    std::vector<int> m_trafoList; // list of the flavours that are polarized spinors
    std::vector<ATOOLS::CMatrix> m_trafoMatrices; // trafo matrices for pol->hel

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
    size_t MaxHel() {return m_nsign;}
    //! Returns the maximum number of helicities for calculation of physical fermion helicities.
  //    int MaxFHel() {return m_fnsign;}
    //! Returns a certain helicity combination.
  int* operator[] (int i) {return p_slist[i].s;} // was fslist!!
    //! Switches a certain helicity combination off.
    void SwitchOff(int i) { if (m_trafoList.size()==0) p_slist[i].on = 0;}
    int  GetPartnerPol(const int heli,const int flav, int& lambda);
    void SetPartner(int i,int j) {if (m_trafoList.size()==0) p_slist[j].partner = i;}
    inline const int&  Partner(const int i) const {return p_slist[i].partner;}
    //! Returns the status of a certain helicity combination.
    int  On(int i) {return p_slist[i].on;}
    //! Increases the multiplicity of a combination by one.
    void IncMultiplicity(int i) {if (m_trafoList.size()==0) (p_slist[i].multi)++;}  
    //! Returns the multiplicity of a certain helicity combination.
    int  Multiplicity(int i) {return p_slist[i].multi;}  
    //! Polarization factors
    double PolarizationFactor(int i) {return p_slist[i].polfactor;}  
    //! Returns true when physical fermion helicities should be calculated
    //  int IsCalcFerm() {return m_fermion_hels;}

    char   * PolTypes()  { return p_pol_types; }
    double * PolAngles() { return p_angles; }
    bool IsContrib(int,int*,int);
    int  Compare(Helicity*,int);

    void SpinorTransformation(std::vector<Complex>& helAmpls);
    void InitializeSpinorTransformation(Basic_Sfuncs * BS);

    inline void ForceNoTransformation() { m_allowTrafo = false; }
    inline void AllowTransformation() { m_allowTrafo = true; }
    inline bool UseTransformation() { return  m_needsTrafo && m_allowTrafo; }
    ~Helicity();

  };
}
#endif
