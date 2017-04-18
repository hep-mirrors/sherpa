#ifndef AHADIC_Tools_Wave_Function_H
#define AHADIC_Tools_Wave_Function_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>


namespace AHADIC {
  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour> Flavour_Pair;
  typedef std::map<Flavour_Pair *,double>            WaveComponents;

  // This is more or less a container class:
  // For every hadron it contains the probabilities related to a decomposition
  // into it constituent quark-antiquark or quark-diquark pairs, the
  // WaveComponents.
  // This information is hard-coded in the Multiplet_Constructor class - not
  // very elegant.
  // The only real functionality it adds is to allow "barring", i.e.
  // constructing the wave function for the anti-hadron.

  class Wave_Function {
  private:
    ATOOLS::Flavour m_hadron;
    int             m_kfcode, m_spin2;
    double          m_mpletwt, m_extrawt;
    WaveComponents  m_waves;
    bool            m_barrable;
  public:
    Wave_Function();
    Wave_Function(const ATOOLS::Flavour &);    
    ~Wave_Function();

    Wave_Function *  GetAnti();
    ATOOLS::Flavour  GetFlavour() const        { return m_hadron; }
    WaveComponents * GetWaves()                { return &m_waves; }
    const size_t     size() const              { return m_waves.size(); }
    
    void   AddToWaves(Flavour_Pair *,double);
    double WaveWeight(ATOOLS::Flavour,ATOOLS::Flavour);

    void   SetMultipletWeight(const double & wt) { m_mpletwt = wt; }
    void   SetExtraWeight(const double & wt)     { m_extrawt = wt; }
    void   SetFlavour(ATOOLS::Flavour & flav)    { m_hadron  = flav; }
    void   SetKfCode(const int & _kfcode)        { m_kfcode  = _kfcode; }
    void   SetSpin(const double & spin)          { m_spin2   = spin; }

    const double & MultipletWeight() const     { return m_mpletwt; }
    const double & ExtraWeight() const         { return m_extrawt; }
    const int    & Spin() const                { return m_spin2; }
    const double   SpinWeight() const          { return double(m_spin2); }
    const int    & KfCode() const              { return m_kfcode; }
    const bool   & ExistAnti() const           { return m_barrable; }

    friend std::ostream & operator<<(std::ostream & s, Wave_Function & wf);
  };

  typedef std::map<ATOOLS::Flavour,Wave_Function *> Wave_Functions; 
}

#endif
