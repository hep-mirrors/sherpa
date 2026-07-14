#ifndef AHADIC_Tools_Wave_Function_H
#define AHADIC_Tools_Wave_Function_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>
#include <vector>


namespace AHADIC {
  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour> Flavour_Pair;
  // component amplitudes are per-variation vectors for the reweighting of
  // the mixing angles; angle-independent components hold a single entry
  typedef std::map<Flavour_Pair *,std::vector<double> > WaveComponents;

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
    std::vector<double> m_mpletwts;
    std::vector<double> m_extrawts;
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
    void   AddToWaves(Flavour_Pair *,const std::vector<double> &);
    double WaveWeight(ATOOLS::Flavour,ATOOLS::Flavour);

    void   SetMultipletWeights(const std::vector<double> & wts) { m_mpletwts = wts; }
    void   SetExtraWeights(const std::vector<double> & wts) { m_extrawts = wts; }
    void   SetFlavour(ATOOLS::Flavour & flav)    { m_hadron  = flav; }
    void   SetKfCode(const int & _kfcode)        { m_kfcode  = _kfcode; }
    void   SetSpin(const double & spin)          { m_spin2   = spin; }

    const std::vector<double> & MultipletWeights() const { return m_mpletwts; }
    const std::vector<double> & ExtraWeights() const { return m_extrawts; }
    const int    & Spin() const                { return m_spin2; }
    const double   SpinWeight() const          { return double(m_spin2); }
    const int    & KfCode() const              { return m_kfcode; }
    const bool   & ExistAnti() const           { return m_barrable; }

    friend std::ostream & operator<<(std::ostream & s, Wave_Function & wf);
  };

  typedef std::map<ATOOLS::Flavour,Wave_Function *> Wave_Functions; 
}

#endif
