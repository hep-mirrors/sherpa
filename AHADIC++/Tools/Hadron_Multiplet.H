#ifndef AHADIC_Tools_Hadron_Multiplet_H
#define AHADIC_Tools_Hadron_Multiplet_H

#include "AHADIC++/Tools/Hadron_Wave_Function.H"
#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>


namespace AHADIC {

  class Hadron_Multiplet {
  private:
    std::string               m_name;
    std::set<ATOOLS::Flavour> m_elements;
    double                    m_spinweight,m_extraweight,m_weight;
  public:
    Hadron_Multiplet() : 
      m_name(std::string("")),
      m_spinweight(1.),m_extraweight(1.),m_weight(1.) {}
    void SetName(std::string name)                 { m_name        = name;  }
    void AddToElements(ATOOLS::Flavour fl)         { m_elements.insert(fl); }
    void SetSpinWeight(double wt)                  { m_spinweight  = wt;    }
    void SetExtraWeight(double wt)                 { m_extraweight = wt;    }
    void SetWeight(double wt=0.)                   { 
      if (wt==0.) m_weight = m_spinweight * m_extraweight;
             else m_weight = wt;
    }
    std::string                 Name()        { return m_name;            }
    std::set<ATOOLS::Flavour> * GetElements() { return &m_elements;       }
    int                         Size()        { return m_elements.size(); }
    double                      SpinWeight()  { return m_spinweight;      }
    double                      ExtraWeight() { return m_extraweight;     }
    double                      Weight()      { return m_weight;          }
  };

  typedef std::map<int,Hadron_Multiplet *>  Hadron_Multiplet_Map;
  typedef Hadron_Multiplet_Map::iterator    Hadron_Multiplet_Miter;

  class All_Hadron_Multiplets {
  private:
    double                  m_singletsuppression, m_octetsuppression;
    Hadron_WF_Map         * p_wavefunctions;
    Hadron_Multiplet_Map  * p_multiplets;

    void ConstructWaveFunctions();
    Hadron_Wave_Function * ConstructMesonWaveFunction(const int,const int,const int,
						      const int,const int,const int);
    Hadron_Wave_Function * ConstructBaryonWaveFunction(int,int,int,int,int);
    void ConstructAntiWaveFunctions();
    void CreateMultiplets();
    void AddMultipletWeights();

    void LookUpAngles(const int angular,const int spin,double & costh,double & sinth);
    bool LookUpWaveFunction(const int kfcode,const int fl1,const int fl2,
			    const Hadron_Wave_Function * wf);
  public:
    All_Hadron_Multiplets();
    ~All_Hadron_Multiplets();

    Hadron_WF_Map         * GetWaveFunctions()    { return p_wavefunctions; }
    Hadron_Multiplet_Map  * GetMultiplets()       { return p_multiplets; }
    Hadron_Multiplet      * GetSU3Pseudoscalars() { 
      if (p_multiplets->find(1)!=p_multiplets->end()) 
	return p_multiplets->find(1)->second; 
      return NULL;
    }

    Hadron_Wave_Function * GetWaveFunction(ATOOLS::Flavour);
    void                   PrintWaveFunctions();
    void                   PrintMultiplets();
  };

  /*!
    \file Hadron_Multiplet.H
    \brief Contains the classes AHADIC::Hadron_Multiplet and AHADIC::All_Hadron_Multiplets.
  */

  /*!
    \class Hadron_Multiplet
    \brief This structure serves for arranging various hadrons into multiplets. It consists mainly
    of a set of the (hadron) flavour, ordered in mass, a name for the multiplet, and weights
    common to the multiplet, such as the spin- or multiplet-weight. The former is given by the
    particles spin \f$s\f$, i.e. by \f$(2s+1)\f$, yielding, e.g., 1 for pseudoscalars and 3 for vector
    mesons.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Hadron_Multiplet
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var std::string Hadron_Multiplet::m_name
    \brief The name of the multiplet.
  */
  /*!
    \var std::set<ATOOLS::Flavour> Hadron_Multiplet::m_elements
    \brief A list of the memebers of the multiplets (like pseudoscalar mesons etc.)
  */
  /*!
    \var double Hadron_Multiplet::m_spinweight
    \brief The spinweight, given by \f$(2s+1)\$ for a multiplet with given pin\f$s\f$.
  */
  /*!
    \var double Hadron_Multiplet::m_extraweight
    \brief Some extra weight to enhance or diminish te probability for transitions into
    that multiplet.  These weights are external parameters of the model.
  */
  /*!
    \var double Hadron_Multiplet::m_weight
    \brief The product of spin and extra weight, a part of the probability that a given colour
    anticolour combination produces a given hadron.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member functions for Hadron_Multiplet
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Hadron_Multiplet::Hadron_Multiplet()
    \brief Default constructor, all weights equal one and the name is void.
  */
  /*!
    \fn void Hadron_Multiplet::SetName(std::string name)
    \brief Sets the name for the multiplet.
  */
  /*!
    \fn void Hadron_Multiplet::AddToElements(ATOOLS::Flavour fl)
    \brief Adds a given (hadron-) flavour to the set of elements, Hadron_Multiplet::m_elements.
  */
  /*!
    \fn void Hadron_Multiplet::SetSpinWeight(double wt)
    \brief Sets the spin weight for a multiplet.
  */
  /*!
    \fn void Hadron_Multiplet::SetExtraWeight(double wt)
    \brief Sets the extra weight for a multiplet.
  */
  /*!
    \fn void Hadron_Multiplet::SetWeight(double wt=0.)
    \brief Sets the weight for a multiplet or calculates it from extra and spin weight.
  */
  /*!
    \fn std::string Hadron_Multiplet::Name()
    \brief The multiplet's name.
  */
  /*!
    \fn std::set<ATOOLS::Flavour> * Hadron_Multiplet::GetElements()
    \brief Returns a pointer to the list of elements, Hadron_Multiplet::m_elements.
  */
  /*!
    \fn int Hadron_Multiplet::Size()
    \brief Number of elements in that multiplet.
  */
  /*!
    \fn double Hadron_Multiplet::SpinWeight()
    \brief Simple access to Hadron_Multiplet::m_spinweight.
  */
  /*!
    \fn double Hadron_Multiplet::ExtraWeight()
    \brief Simple access to Hadron_Multiplet::m_extraweight.
  */
  /*!
    \fn double Hadron_Multiplet::Weight()
    \brief Simple access to Hadron_Multiplet::m_weight.
  */

  /*!
    \class All_Hadron_Multiplets
    \brief This class sets up and manages the various hadron multiplets. 

    The full work of this class is steered by the constructor, which, first, iterates over all hadrons
    that are switched on, and constructs their wave functions. These wavefunctions are
    stored in the All_Hadron_Multiplets::p_wavefunctions and sorted into different \f$SU(6)\f$ multiplets, 
    stored in All_Hadron_Multiplets::p_multiplets. This class is initialised by the object
    Hadronisation_Parameters and made available to the rest of the code through a corresponding Get-method, 
    AHADIC::Hadronisation_Parameters::GetMultiplets().
  */


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for All_Hadron_Multiplets
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Hadron_WF_Map * All_Hadron_Multiplets::p_wavefunctions
    \brief The list of all wavefunctions (realised through AHADIC::Hadron_Wave_Function)
    that have been constructed by this class.
  */
  /*!
    \var Hadron_Multiplet_Map * All_Hadron_Multiplets::p_multiplets
    \brief The list of all hadron multiplets, they are sorted by a number given
    by (+/-) for baryons/antibaryons (always + for mesons), and 2s+1+100*L+1000*R, and
    eventually another 10000 for the "wrong" wave function parity.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for All_Hadron_Multiplets
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn All_Hadron_Multiplets::All_Hadron_Multiplets()
    \brief In the constructor all wave functions are constructed through calls to
    All_Hadron_Multiplets::ConstructWaveFunctions() and 
    All_Hadron_Multiplets::ConstructAntiWaveFunctions(); then the hadrons are arranged in multiplets, 
    by All_Hadron_Multiplets::CreateMultiplets(), before specific multiplet weights are added through 
    All_Hadron_Multiplets::AddMultipletWeights().
  */
  /*!
    \fn All_Hadron_Multiplets::~All_Hadron_Multiplets()
    \brief Deletes both the Hadron_WF_Map, All_Hadron_Multiplets::p_wavefunctions, and 
    the Hadron_Multiplet_Map, All_Hadron_Multiplets::p_multiplets.
  */
  /*!
    \fn Hadron_Multiplet_Map * All_Hadron_Multiplets::GetMultiplets()
    \brief Returns the constructed hadron multiplets, All_Hadron_Multiplets::p_mulitplets, which are 
    arranged in AHADIC::Hadron_Multiplet_Map (a map of ints and hadron multiplets).
  */
  /*!
    \fn Hadron_Multiplet_Map * All_Hadron_Multiplets::GetSU3Pseudoscalars()
    \brief Returns the pseudoscalar hadron multiplet, if existent.
  */
  /*!
    \fn void All_Hadron_Multiplets::PrintMultiplets()
    \brief Prints all hadron multiplets.
  */
  /*!
    \fn Hadron_WF_Map * All_Hadron_Multiplets::GetWaveFunctions()
    \brief Returns all hadron wave functions, which are arranged in AHADIC::Hadron_WF_Map
    (a map of (hadron) flavours) and the corresponding pointers to their AHADIC::Hadron_Wave_Function.
  */
  /*!
    \fn Hadron_Wave_Function * All_Hadron_Multiplets::GetWaveFunction(ATOOLS::Flavour)
    \brief Returns the pointer to a specific AHADIC::Hadron_Wave_Function, corresponding to
    the respective hadon flavour given in the argument.
  */
  /*!
    \fn void All_Hadron_Multiplets::PrintWaveFunctions()
    \brief Prints all hadron wave functions, arranged in multiplets.
  */
  /*!
    \fn void All_Hadron_Multiplets::ConstructWaveFunctions()
    \brief Constructs all wave functions through corresponding calls to specialised methods.

    In so doing, this method iterates over all flavours in the particle list, kf_table and 
    analyses the respective codes. From the codes individual quark flavours can be extracted for 
    hadrons; depending on whether there are two or three flavours specified, either 
    All_Hadron_Multiplets::ConstructMesonWaveFunction(const int,const int,const int,const int,const int,const int) 
    or All_Hadron_Multiplets::ConstructBaryonWaveFunction(int,int,int,int,int), respectively, is 
    called.  The resulting hadron wave function is added to the list of wavefunctions, realised through 
    the object All_Hadron_Multiplets::p_wavefunctions (Hadron_WF_Map).
  */
  /*!
    \fn Hadron_Wave_Function * All_Hadron_Multiplets::ConstructMesonWaveFunction(int,int,int,int,int,int)
    \brief Here the meson wave functions are constructed from the spin, the angular momentum
    and the two quark flavour indices \f$q_1\f$ and \f$\bar q_2\f$. 

    There are a number of cases to consider: 
    - The simplemost one is, when both quark flavour indices are different. Then the wave function 
    is just given by
    \f[
    \Psi_{q_1\bar q_2} = \left|q_1\bar q_2\right\rangle.
    \f]
    - If both indices are identical things are a bit more complicated since in this case there are
      differnt components. Things become even more messy when singlet-octet mixing, parametrised through
      an angle \f$\theta\f$, is included. 
      - However, for \f$q_1 = \bar q_2 = 1\f$ (d-quarks) the wave function reads
      \f[
      \Psi_{11} = \frac{1}{\sqrt{2}}\left(\left|u\bar u\right\rangle-\left|d\bar d\right\rangle\right).
      \f]
      - Furthermore, for \f$q_1 = \bar q_2 = 2\f$ (u-quarks) the wave function reads
      \f[
      \Psi_{22} = 
      \left(\frac{\cos\theta}{\sqrt{3}}+\frac{\sin\theta}{\sqrt{6}}\right)
      \left(\left|u\bar u\right\rangle+\left|d\bar d\right\rangle\right)+
      \left(\frac{\cos\theta}{\sqrt{3}}-\frac{2\sin\theta}{\sqrt{6}}\right)
      \left(\left|s\bar s\right\rangle\right).
      \f]
      - Finally, if \f$q_1 = \bar q_2 = 3\f$ (s-quarks) the wave function reads
      \f[
      \Psi_{33} = 
      \left(\frac{2\cos\theta}{\sqrt{6}}+\frac{\sin\theta}{\sqrt{3}}\right)
      \left(\left|s\bar s\right\rangle\right)+
      \left(\frac{\cos\theta}{\sqrt{6}}-\frac{\sin\theta}{\sqrt{3}}\right)
      \left(\left|u\bar u\right\rangle+\left|d\bar d\right\rangle\right).
      \f] 
      - For the heavy quarks the wave function just reads
      \f[
      \Psi_{QQ} = \left|Q\bar Q\right\rangle.
      \f]
  */
  /*!
    \fn Hadron_Wave_Function * All_Hadron_Multiplets::ConstructBaryonWaveFunction(int,int,int,int,int)
    \brief Here the baryon wave functions are constructed from the spin, the angular momentum
    and the three quark flavour indices \f$q_1\f$, \f$q_2\f$, and \f$q_3\f$. The problem here is to
    map them onto \f$SU(6)\f$ wave functions made of quark-diquark pairs. A solution can be found 
    in D.B.Lichtenberg, W.Namgung, J.G.Wills, and E.Predazzi, Z.Phys.C19 (1983) 19.

    In the following \f$d^{S,V}_{ij}\f$ stands for scalar or vector diquarks made of the quarks i and j. 
    For the octet baryons, there are various cases to be considered:
    - Proton-like baryons: Two identical quarks, then the wave function reads
      \f[
      \Psi_{221} = 
      \frac{1}{\sqrt{3}}\left|q_1d^S_{22}\right\rangle+
      \frac{1}{\sqrt{6}}\left|q_2d^V_{12}\right\rangle+
      \frac{1}{\sqrt{2}}\left|q_2d^S_{12}\right\rangle.
      \f]
    - \f$\Lambda\f$-like baryons: Three different quarks, then the wave function reads
      \f[
      \Psi_{312} =
      \frac{1}{\sqrt{12}}\left|q_3d^S_{12}\right\rangle+
      \frac{1}{\sqrt{12}}\left|q_1d^S_{23}\right\rangle+
      \frac{1}{\sqrt{12}}\left|q_2d^S_{13}\right\rangle+
      \frac{1}{\sqrt{4}}\left|q_3d^V_{12}\right\rangle+
      \frac{1}{\sqrt{4}}\left|q_1d^V_{23}\right\rangle+
      \frac{1}{\sqrt{4}}\left|q_2d^V_{13}\right\rangle.
      \f]
    - \f$\Sigma\f$-like baryons: Three different quarks, then the wave function reads
      \f[
      \Psi_{321} =
      \frac{1}{\sqrt{12}}\left|q_3d^V_{12}\right\rangle+
      \frac{1}{\sqrt{12}}\left|q_1d^V_{23}\right\rangle+
      \frac{1}{\sqrt{12}}\left|q_2d^V_{13}\right\rangle+
      \frac{1}{\sqrt{4}}\left|q_3d^S_{12}\right\rangle+
      \frac{1}{\sqrt{4}}\left|q_1d^S_{23}\right\rangle+
      \frac{1}{\sqrt{4}}\left|q_2d^S_{13}\right\rangle.
      \f]

    For the decuplet baryons, there are again various cases to be considered:
    - \f$\Delta^{++,-}\f$-like baryons: Three identical quarks with wave function
      \f[
      \Psi_{111} = \left|q_1d^V_{11}\right\rangle.
      \f]
    - \f$\Delta^{+,0}\f$-like baryons: Two identical quarks with wave function
      \f[
      \Psi_{221} = 
      \frac{1}{\sqrt{3}}\left|q_1d^V_{22}\right\rangle+
      \frac{2}{\sqrt{3}}\left|q_2d^V_{12}\right\rangle.
      \f]
    - \f$\Sigma\f$(1385)-like baryons:
      \f[
      \Psi_{321} = 
      \frac{1}{\sqrt{3}}\left|q_1d^V_{23}\right\rangle+\mbox{\rm cycl.}.
      \f]

    \todo Here, more work has to be done!!!   
    For baryons containing one heavy quark, there are the following wavefunctions
    - \f$\Sigma_c\f$-like baryons:
      \f[
      \Psi_{Qqq}^{1/2,3/2} = 
      \left|Qd^{S,V}_{qq}\right\rangle,
      \f]
      where the diquark spin depends on whether the baryons has spin 1/2 (S) or 3/2 (V).
  */
  /*!
    \fn void All_Hadron_Multiplets::ConstructAntiWaveFunctions()
    \brief Here, wavefunctions for the anti-particles are constructed, where needed. This
    method fully relies on Hadron_Wave_Function::Anti().
  */
  /*!
    \fn void All_Hadron_Multiplets::CreateMultiplets()
    \brief This method takes the wave functions constructed so far and, according to the spin
    of the respective hadrons, places them in different multiplets. 
    
    These multiplets basically are \f$SU(3)_F\f$ or \f$SU(3)_F\otimes SU(2)_S\f$ multiplets,
    and thus constrained on the light quarks. Hadrons containing a heavy quark are treated separately
    such that "reduced" multiplets are constructed, which are \f$SU(3)_F\f$ symmetric only in
    the light quark sector.
  */
  /*!
    \fn void All_Hadron_Multiplets::AddMultipletWeights()
    \brief Weights are added to the multiplets, in particular spin weights taking care of the
    spin d.o.f. in each multiplet. The addition of further, multiplet-specific weights that allow
    to steer abundances beyond phasespace, spin, etc., is envisioned but needs to be done.
    
    \todo Add multiplet-specific weights!!!
  */
  /*!
    \fn void All_Hadron_Multiplets::LookUpAngles(const int angular,const int spin,double & costh,double & sinth)
    \brief This method is used to steer mixing of flavour eigenfunctions to hadron wavefunctions.

    Examples for this include the flavour content of the \f$\eta\f$ particles in the pseudoscalar multiplet
    or \f$\omega-\phi\f$-mixing in the vector octet. The angles have to be given explicitly in the
    corresponding Cluster.dat file, otherwise they are set to the following default values:
    - \f$\theta_{0^-} = -0.3491\f$
    - \f$\theta_{1^+} = 0.6155\f$.

    \todo Further mixing angles for higher multiplets to follow.
  */
  /*!
    \fn bool All_Hadron_Multiplets::LookUpWaveFunction(const int,const int,const int,const Hadron_Wave_Function *)
    \brief Not realised yet.
  */
}

#endif
