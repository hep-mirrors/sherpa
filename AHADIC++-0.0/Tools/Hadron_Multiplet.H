#ifndef Hadron_Multiplet_H
#define Hadron_Multiplet_H

#include "Hadron_Wave_Function.H"
#include "Flavour.H"
#include <map>
#include <string>


namespace AHADIC {

  class Hadron_Multiplet {
  private:
    std::string        m_name;
    ATOOLS::FlavourSet m_elements;
    double             m_spinweight,m_extraweight,m_weight;
  public:
    Hadron_Multiplet() : 
      m_name(std::string("")),
      m_spinweight(1.),m_extraweight(1.),m_weight(1.) {}
    void SetName(std::string name)                 { m_name        = name;  }
    void AddToElements(ATOOLS::Flavour & fl)       { m_elements.insert(fl); }
    void AddToElements(const ATOOLS::Flavour & fl) { m_elements.insert(fl); }
    void SetSpinWeight(double wt)                  { m_spinweight  = wt;    }
    void SetExtraWeight(double wt)                 { m_extraweight = wt;    }
    void SetWeight(double wt=0.)                   { 
      if (wt==0.) m_weight = m_spinweight * m_extraweight;
             else m_weight = wt;
    }
    std::string          Name()        { return m_name;            }
    ATOOLS::FlavourSet * GetElements() { return &m_elements;       }
    int                  Size()        { return m_elements.size(); }
    double               SpinWeight()  { return m_spinweight;      }
    double               ExtraWeight() { return m_extraweight;     }
    double               Weight()      { return m_weight;          }
  };

  typedef std::map<int,Hadron_Multiplet *>  Hadron_Multiplet_Map;
  typedef Hadron_Multiplet_Map::iterator    Hadron_Multiplet_Miter;

  class All_Hadron_Multiplets {
  private:
    Hadron_WF_Map         * p_wavefunctions;
    Hadron_Multiplet_Map  * p_multiplets;

    void ConstructWaveFunctions();
    Hadron_Wave_Function * ConstructMesonWaveFunction(int,int,int,int);
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

    Hadron_WF_Map        * GetWaveFunctions() { return p_wavefunctions; }

    Hadron_Wave_Function * GetWaveFunction(ATOOLS::Flavour);
    void                   PrintWaveFunctions();
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
  /*!
    \class All_Hadron_Multiplets
    \brief This class sets up and manages the various hadron multiplets.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for All_Hadron_Multiplets
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Hadron_WF_Map * All_Hadron_Multiplets::p_wavefunctions
    \brief
  */
  /*!
    \var Hadron_Multiplet_Map * All_Hadron_Multiplets::p_multiplets
    \brief
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for All_Hadron_Multiplets
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn All_Hadron_Multiplets::All_Hadron_Multiplets()
    \brief In the constructor all wave functions are constructed through calls to
    ConstructWaveFunctions() and ConstructAntiWaveFunctions(); then the hadrons are
    arranged in multiplets, by CreateMultiplets(), before specific multiplet weights are
    added through AddMultipletWeights().
  */
  /*!
    \fn All_Hadron_Multiplets::~All_Hadron_Multiplets()
    \brief Deletes both the Hadron_WF_Map, p_wavefunctions, and the Hadron_Multiplet_Map,
    p_multiplets.
  */
  /*!
    \fn Hadron_WF_Map * All_Hadron_Multiplets::GetWaveFunctions()
    \brief Returns all hadron wave funxctions, which are arranged in a map of (hadron) flavours
    and the corresponding pointers to their AHADIC::Hadron_Wave_Function.
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

    In so doing, this method iterates over all flavours in the particle list, ATOOLS::kf_table and 
    analyses the respective codes. From the codes individual quark flavours can be extracted for 
    hadrons; depending on whether there are two or three flavours specified, either 
    ConstructMesonWaveFunction(int,int,int,int) or ConstructBaryonWaveFunction(int,int,int,int,int), 
    respectively, is called. The resulting hadron wave function is added to the list of wavefunctions, 
    realised through the object p_wavefunctions (Hadron_WF_Map).
  */
  /*!
    \fn Hadron_Wave_Function * All_Hadron_Multiplets::ConstructMesonWaveFunction(int,int,int,int)
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
    \brief 
  */
  /*!
    \fn void All_Hadron_Multiplets::ConstructAntiWaveFunctions()
    \brief
  */
  /*!
    \fn void All_Hadron_Multiplets::CreateMultiplets()
    \brief
  */
  /*!
    \fn void All_Hadron_Multiplets::AddMultipletWeights()
    \brief
  */
  /*!
    \fn void All_Hadron_Multiplets::LookUpAngles(const int angular,const int spin,double & costh,double & sinth)
    \brief
  */
  /*!
    \fn bool All_Hadron_Multiplets::LookUpWaveFunction(const int,const int,const int,const Hadron_Wave_Function *)
    \brief
  */
}

#endif
