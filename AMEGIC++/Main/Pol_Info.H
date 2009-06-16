#ifndef AMEGIC_Main_Pol_Info_H
#define AMEGIC_Main_Pol_Info_H

#define Explicit_Pols

#include "ATOOLS/Phys/Flavour.H"
//syntax for polarised particles: (in processes.dat)
//
//    degree of polarisation (0.=unpolarized,1.=fully polarized) (only 1 or 0 for outgoing)
//                        | l for linear
//                        | | polarization angle
//                        | | |
//   Bosons:   linear:   (.8l45)
//                        |
//                        | +/-/0(for massive bosons)
//                        |  |
//             circular: (.8c+)
//                        |
//                        | +/-
//                        |  |
//   Fermion helicities: (.8h+) (up to now only ok for massless)

namespace AMEGIC {

  class mt {
  public:
    enum momtype{  mom    = 10,  //external momentum
		   prop   = 11,  //propagator
		   cmprop = 12,  //center of momentum propagator
		   p_p    = 1,   //polarisation +1 (circular)
		   p_m    =-1,   //             -1 (circular)
		   p_l    = 0,   //             longitudinal
		   p_s    = 2,   //             scalar
		   p_si   = 7,   //             massless scalar
		   p_l0   = 3,   //             linear with direction
		   p_l1   = 4,   //             linear, perpendicular to p_l0
		   p_lh   = 5,   //             horizontaly (linear)
		   p_lp   = 6,    //             perpendicular (linear)
		   p_none = 8,
		   p_spec = 9,
		   
		   p_t1   = 101,  //dummies for graviton polarization tensors
		   p_t2   = 102,
		   p_t3   = 103,
		   p_t4   = 104,
		   p_t5   = 105
    };
  };
  /*!
    \class mt
    \brief Types for explicit 4-vectors.

    This class specifies the types for explicit 4-vectors in AMEGIC::Momfunc. 
    They are used to label the arguments in AMEGIC::Basic_Func and derived classes. 

    The following types are defined:
    - mom:         momentum of an external particle  
    - prop:        momentum for a propagator
    - cmprop:      momentum for a propagator, where the spacial part is identical 0
    - p_p/p_m:     circular polarization vectors (+/-)
    - p_l:         longitudinal polarization
    - p_s/p_si:    scalar polarization for massive/massless off-shell vector boson
    - p_l0/p_l1:   linear polarization vectors with a direction (defined in AMEGIC::Pol_Info::angle)
    - p_lh/p_lp:   linear polarization vectors with fixed direction
    - p_t1-p_t5:   dummies for polarization tensors of external spin-2 particles, will be replaced in AMEGIC::Single_Amplitude_Base::Single_ZvalueTensor(Zfunc* z,std::vector<int>*,std::vector<int>*,int) by combinations of polarization vectors

    For the calculation of these vectors see AMEGIC::Basic_Sfuncs::Calc_Momlist().
  */
    
  class Pol_Info{
  public:
    //Polarisation type: 'c'=circular 'l'=linear 'h'=helicity 't'=tensor
    char pol_type;
    //direction for linear polarisation
    double angle;
    int      num;
    int    * type;
    double * factor;
    void Init(int);
    friend std::ostream & operator<<(std::ostream &, Pol_Info &);
  public:
    Pol_Info();
    Pol_Info(const ATOOLS::Flavour &);
    Pol_Info(const Pol_Info &);
    Pol_Info& operator=(const Pol_Info &);
   ~Pol_Info();

    void SetPol(char);
    char GetPol();
    char GetPolarisationType() { return pol_type; }
    int  DoFNumber() { return num; }
    int    Polarisation(const int i) { return type[i]; }
    double Degree(const int i)       { return factor[i]; }
  };

  std::ostream & operator<<(std::ostream &, Pol_Info &);

  /*!
    \class Pol_Info
    \brief Container for everything about the polarization state of external particles.

    This class is a container to transfer the polarization state information from 
    AMEGIC::Amegic::ReadProcesses(std::string) to the construction of the helicity states
    in AMEGIC::Helicity.
  */
  /*!
    \var char Pol_Info::pol_type
    Type of the polarization state.
    Possible values:
    - c for circular (vector boson)
    - l for linear (vector boson)
    - h for fermion helicities
    - t for tensor (spin-2 particles)    
  */
  /*!
    \var double Pol_Info::angle
    angle for linear polarization vectors of type mt::p_l0 and mt::p_l1
  */
  /*!
    \var int Pol_Info::num
    number of polarization vectors (tensors) to calculate a given polarization state
    Possible values:
    - 1 for totally polarized beams/final states
    - 2 for partially polarized/unpolarized massless vector bosons
    - 3 for partially polarized/unpolarized massive vector bosons
    - 5 for massive spin-2 tensors 
  */
  /*!
    \var int* Pol_Info::type
    Array of num AMEGIC::mt::momtype 
  */
  /*!
    \var double * Pol_Info::factor
    The weight of each polarization vector for the right partial polarization
  */  


  class Tensor_Struc{                //defines the spin2 tensors
  public:
    void GetPolCombos(int, std::vector<std::vector<int> >*, std::vector<int>*);
    double GetTfactor(int);
  };
  /*!
    \class Tensor_Struc
    \brief Defines tensors for external spin-2 particles 

    This class defines the 5 polarization tensors for a massive external spin-2 particle
    in terms of (outer) products of two polarization vectors.

    Polarization tensors are defined as follows:
    - AMEGIC::mt::p_t1 \f$\epsilon_{\mu}^{+}\epsilon_{\nu}^{+}\f$
    - AMEGIC::mt::p_t2 \f$\sqrt{2}\,\epsilon_{\mu}^{+}\epsilon_{\nu}^{0}\f$
    - AMEGIC::mt::p_t3 \f$\sqrt{\frac{2}{3}}\left(\epsilon_{\mu}^{+}\epsilon_{\nu}^{-}-\epsilon_{\mu}^{0}\epsilon_{\nu}^{0}\right)\f$
    - AMEGIC::mt::p_t4 \f$\sqrt{2}\,\epsilon_{\mu}^{-}\epsilon_{\nu}^{0}\f$
    - AMEGIC::mt::p_t5 \f$\epsilon_{\mu}^{-}\epsilon_{\nu}^{-}\f$
  */
  /*!
    \fn void Tensor_Struc::GetPolCombos(int, vector<vector<int> >*, vector<int>*)
    Gives a list of polarization vectors and relative signs to define a polarization tensor.
  */
  /*!
    \fn double Tensor_Struc::GetTfactor(int) 
    Returns a normalization factor for each polarization tensor.
  */
}
#endif

