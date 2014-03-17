#ifndef MCATNLO_Tools_Singlet_H
#define MCATNLO_Tools_Singlet_H

#include "MCATNLO/Tools/Parton.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS { class Mass_Selector; }

namespace PDF { class Shower_Base; }

namespace PHASIC { class Jet_Finder; }

namespace MCATNLO {

  class Singlet;
  class Sudakov;

  typedef std::vector<Singlet*> All_Singlets;
  typedef All_Singlets::iterator ASiter;

  class Singlet : public Parton_List {
    PHASIC::Jet_Finder *p_jf;
    PDF::Shower_Base *p_shower;
    All_Singlets *p_all;
    ATOOLS::Mass_Selector *p_ms;
    void *p_procs;
    double m_mur2;
    Parton_List m_dels;
  public :
    Singlet() : 
      p_jf(NULL), p_all(NULL), p_procs(NULL) {}
    ~Singlet();

    int       SplitParton(Parton *, Parton *, Parton *);
    bool      ArrangeColours(Parton *, Parton *, Parton *);

    void        BoostAllFS(Parton *l,Parton *r,Parton *s);
    void        BoostBackAllFS(Parton *l,Parton *r,Parton *s);

    bool JetVeto(Sudakov *const sud) const;
      
    friend std::ostream& operator<<(std::ostream &,Singlet &);

    inline void SetShower(PDF::Shower_Base *const s) { p_shower=s; }

    inline void SetJF(PHASIC::Jet_Finder *const jf) { p_jf=jf; }
    inline void SetAll(All_Singlets *const all) { p_all=all; }

    inline PHASIC::Jet_Finder *JF() const { return p_jf; }
    inline All_Singlets *All() const { return p_all; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

    inline void SetProcs(void *const procs) { p_procs=procs; }

    inline double MuR2() const { return m_mur2; }

    inline void SetMuR2(const double &mu2) { m_mur2=mu2; }

    template <class Type> inline Type *Procs() const 
    { return static_cast<Type*>(p_procs); }

  };

  std::ostream& operator<<(std::ostream &,All_Singlets &); 
}

#endif
