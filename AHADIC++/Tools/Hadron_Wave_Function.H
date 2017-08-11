#ifndef AHADIC_Tools_Hadron_Wave_Function_H
#define AHADIC_Tools_Hadron_Wave_Function_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>


namespace AHADIC {
  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour> Flavour_Pair;
  typedef std::map<Flavour_Pair *,double>            WFcomponent;
  typedef WFcomponent::iterator                      WFcompiter;


  class Hadron_Wave_Function {
  private:
    ATOOLS::Flavour m_hadron;
    int             m_spin, m_kfcode;
    double          m_mpletwt;
    WFcomponent     m_waves;
    bool            m_barrable;
  public:
    Hadron_Wave_Function();
    Hadron_Wave_Function(const ATOOLS::Flavour &);    
    ~Hadron_Wave_Function();

    void AddToWaves(Flavour_Pair *,double);
    Hadron_Wave_Function * Anti();
    double WaveWeight(ATOOLS::Flavour,ATOOLS::Flavour);

    double MultipletWeight()                 { return m_mpletwt;    }
    void SetMultipletWeight(const double wt) { m_mpletwt = wt;      }
    void SetFlavour(ATOOLS::Flavour & flav)  { m_hadron  = flav;    }
    ATOOLS::Flavour GetFlavour() const       { return m_hadron;     }
    int  KfCode() const                      { return m_kfcode;     }
    void SetKfCode(const int _kfcode)        { m_kfcode  = _kfcode; }
    int  Spin() const                        { return m_spin;       }
    void SetSpin(const int _spin)            { m_spin    = _spin;   }
    WFcomponent * GetWaves()                 { return &m_waves;     }
    bool ExistAnti()                         { return m_barrable;   }

    friend std::ostream & 
    operator<<(std::ostream & s, Hadron_Wave_Function & wf);
  };

  typedef std::map<ATOOLS::Flavour,Hadron_Wave_Function *> Hadron_WF_Map; 
  typedef Hadron_WF_Map::iterator                          Hadron_WF_Miter;

  /*!
    \file Hadron_Wave_Function.H
    \brief Contains the class AHADIC::Hadron_Wave_Function
  */

  /*!
    \class Hadron_Wave_Function
    \brief A class that describes the wave function of one hadron plus a few operations possible
    with such a wave function. 

    The wave functions are repesented as maps of flavour pairs, Flavour_Pair = 
    pair<ATOOLS::Flavour,ATOOLS::Flavour>), and doubles, which equal their respective squared wave 
    function contribution. Additional information consists of a specific multiplet weight with 
    which relative occurence of hadron species, sorted by their \f$ SU(3)_F\f$ /
    \f$ SU(6)_{F\otimes S}\f$ multiplet for mesons/baryons, can be steered.  Also, information about 
    the hadron, whether barring it affects the wave function, etc.. The Hadron_Wave_Function objects 
    are created by the AHADIC::All_Hadron_Multiplets, in particular through the methods called from 
    All_Hadron_Multiplets::ConstructWaveFunctions().

    Since this is mainly a container class, there are only a few methods with actual actions.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Hadron_Wave_Function
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var ATOOLS::Flavour Hadron_Wave_Function::m_hadron
    \brief The hadron for which the wave function is set up.
  */
  /*!
    \var int Hadron_Wave_Function::m_spin
    \brief The hadron's spin.
  */
  /*!
    \var int Hadron_Wave_Function::m_kfcode
    \brief The hadron's kfcode.
  */
  /*!
    \var double Hadron_Wave_Function::m_mpletwt
    \brief An extra weight which can be set per mutliplet from the outside. Will, e.g, be used to
    suppress resonances with a relative angular momentum.
  */
  /*!
    \var WFcomponent Hadron_Wave_Function::m_waves
    \brief A map of triplet-antitriplet pairs and double numbers, representing different
    individual contributions to the total wave function and their weight (amplitude squared).
    This is the central piece of the class, most of the methods are just operating on this map.
  */
  /*!
    \var bool Hadron_Wave_Function::m_barrable
    \brief Indicates whether the wave function barred is different.

    By default, this switch is "false", unless any of the components of the wave
    function has equal triplet and barred antitriplet flavours. Then the particle equals the
    antiparticle and does not need to be treated further to account for anti-particles. Examples
    include \f$\pi_0\f$ and the \f$\eta\f$-particles.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Hadron_Wave_Function
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Hadron_Wave_Function::Hadron_Wave_Function()
    \brief Constructor of a void hadron wave function with no specified flavour, no weights, etc..
  */
  /*!
    \fn Hadron_Wave_Function::Hadron_Wave_Function(const ATOOLS::Flavour &)
    \brief Constructor of a wave function for a specific flavour.
  */
  /*!
    \fn Hadron_Wave_Function::~Hadron_Wave_Function()
    \brief Deletes the WFcomponent, i.e the m_waves.
  */
  /*!
    \fn void Hadron_Wave_Function::AddToWaves(Flavour_Pair *,double)
    \brief Adds a certain contribution, i.e. a triplet-antitriplet pair, given by the flavour pair
    and its squared amplitude, i.e. weight, to the total wave function.
  */
  /*!
    \fn Hadron_Wave_Function * Hadron_Wave_Function::Anti()
    \brief If possible, i.e. if m_barrable=true, the Hadron_Wave_Function of the antiparticle
    is constructed by simply barring all ingredients of the particle's wavefunction and the 
    pointer to the new Hadron_Wave_Function is returned. 
  */
  /*!
    \fn double Hadron_Wave_Function::WaveWeight(ATOOLS::Flavour,ATOOLS::Flavour)
    \brief Returns the weight for a triplet-antitriplet contribution to the total wave function.
    If this does not contribute, the method yields zero.
  */
}

#endif
