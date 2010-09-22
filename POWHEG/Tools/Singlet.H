#ifndef POWHEG_Tools_Singlet_H
#define POWHEG_Tools_Singlet_H

#include "POWHEG/Tools/Parton.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS { class Mass_Selector; }

namespace PHASIC { class Jet_Finder; }

namespace POWHEG {

  class Singlet;
  class Sudakov;

  typedef std::vector<Singlet*> All_Singlets;
  typedef All_Singlets::iterator ASiter;

  class Singlet : public Parton_List {
    PHASIC::Jet_Finder *p_jf;
    All_Singlets *p_all;
    ATOOLS::Mass_Selector *p_ms;
    ATOOLS::RB_Map *p_rbmap;
    void *p_procs;
    bool ArrangeColours(Parton *, Parton *, Parton *);
  public :
    Singlet() : 
      p_jf(NULL), p_all(NULL), p_rbmap(NULL), p_procs(NULL) {}
    ~Singlet();

    int       SplitParton(Parton *, Parton *, Parton *);

    void        BoostAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			   const ATOOLS::Flavour &mo,const int mode);
    void        BoostBackAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			       const ATOOLS::Flavour &mo,const int mode);

    bool JetVeto(Sudakov *const sud,
		 PHASIC::Jet_Finder *const jf,const double &q2,
		 const ATOOLS::Flavour &fj,const ATOOLS::Vec4D &pj) const;
      
    friend std::ostream& operator<<(std::ostream &,Singlet &);

    inline void SetJF(PHASIC::Jet_Finder *const jf) { p_jf=jf; }
    inline void SetAll(All_Singlets *const all) { p_all=all; }

    inline PHASIC::Jet_Finder *JF() const { return p_jf; }
    inline All_Singlets *All() const { return p_all; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

    inline void SetRBMap(ATOOLS::RB_Map *const rbmap) { p_rbmap=rbmap; }

    inline void SetProcs(void *const procs) { p_procs=procs; }

    inline ATOOLS::RB_Map *RBMap() const { return p_rbmap; }

    template <class Type> inline Type *Procs() const 
    { return static_cast<Type*>(p_procs); }

  };

  std::ostream& operator<<(std::ostream &,All_Singlets &); 
}

#endif
