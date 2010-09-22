#ifndef CSSHOWER_Tools_Singlet_H
#define CSSHOWER_Tools_Singlet_H

#include "CSSHOWER++/Tools/Parton.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS { class Mass_Selector; }

namespace PHASIC { class Jet_Finder; }

namespace CSSHOWER {

  class Singlet;
  class Sudakov;

  typedef std::vector<Singlet*> All_Singlets;
  typedef All_Singlets::iterator ASiter;

  class Singlet : public Parton_List {
    Parton *p_split, *p_left, *p_right, *p_spec;
    PHASIC::Jet_Finder *p_jf;
    All_Singlets *p_all;
    ATOOLS::Mass_Selector *p_ms;
    int m_nlo;
    void SetColours(Singlet *const sing,Parton *const np,
		      const int oldc[2],const int newc[2]);
    void SetColours(Parton *const pa,Parton *const p,
		    const int oldc[2],const int newc[2]);
  public :
    Singlet() : 
      p_split(NULL), p_left(NULL), p_right(NULL), p_spec(NULL),
      p_jf(NULL), p_all(NULL), m_nlo(0) {}
    ~Singlet();

    void AddParton(Parton *const p);
    void RemoveParton(Parton *const p,const int mode=0);

    bool RearrangeColours(Parton *, Parton *, Parton *);
    bool ArrangeColours(Parton *, Parton *, Parton *);
    int  SplitParton(Parton *, Parton *, Parton *);

    void      ExtractPartons(ATOOLS::Blob *,ATOOLS::Mass_Selector *const ms);

    void        UpdateDaughters();
    void        BoostAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			   const ATOOLS::Flavour &mo,const int mode);
    void        BoostBackAllFS(Parton *l,Parton *r,Parton *s,Parton *f,
			       const ATOOLS::Flavour &mo,const int mode);

    bool JetVeto(Sudakov *const sud,
		 PHASIC::Jet_Finder *const jf,const double &q2,
		 const ATOOLS::Flavour &fj,const ATOOLS::Vec4D &pj) const;
      
    friend std::ostream& operator<<(std::ostream &,Singlet &);

    inline void SetSplit(Parton *const part) { p_split=part; }
    inline void SetLeft(Parton *const part)  { p_left=part;  }
    inline void SetRight(Parton *const part) { p_right=part; }
    inline void SetSpec(Parton *const part)  { p_spec=part;  }
 
    inline Parton *GetSplit() const { return p_split; }
    inline Parton *GetLeft() const  { return p_left;  }
    inline Parton *GetRight() const { return p_right; }
    inline Parton *GetSpec() const  { return p_spec;  }

    inline void SetJF(PHASIC::Jet_Finder *const jf) { p_jf=jf; }
    inline void SetAll(All_Singlets *const all) { p_all=all; }

    inline PHASIC::Jet_Finder *JF() const { return p_jf; }
    inline All_Singlets *All() const { return p_all; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

    inline void SetNLO(const int nlo) { m_nlo=nlo; }

    inline int NLO() const { return m_nlo; }

  };

  std::ostream& operator<<(std::ostream &,All_Singlets &); 
}

#endif
