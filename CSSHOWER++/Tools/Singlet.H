#ifndef CSSHOWER_Tools_Singlet_H
#define CSSHOWER_Tools_Singlet_H

#include "CSSHOWER++/Tools/Parton.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace ATOOLS { class Mass_Selector; }

namespace PHASIC { class Jet_Finder; }

namespace CSSHOWER {

  class Singlet;
  class Sudakov;
  class Shower;

  typedef std::vector<Singlet*> All_Singlets;
  typedef All_Singlets::iterator ASiter;

  class Singlet : public Parton_List {
    Parton *p_split, *p_left, *p_right, *p_spec;
    PHASIC::Jet_Finder *p_jf;
    All_Singlets *p_all;
    ATOOLS::Mass_Selector *p_ms;
    Shower *p_shower;
    double m_kt2_next, m_lkf, m_mur2;
    int m_nlo, m_nskip, m_nme, m_nmax;
    ATOOLS::DecayInfo_Vector m_decs;
  public :
    Singlet() : 
      p_split(NULL), p_left(NULL), p_right(NULL), p_spec(NULL),
      p_jf(NULL), p_all(NULL),
      m_kt2_next(0.0), m_lkf(0.0), m_mur2(0.0),
      m_nlo(0), m_nskip(0), m_nme(0), m_nmax(0) {}
    ~Singlet();

    void AddParton(Parton *const p);
    void RemoveParton(Parton *const p,const int mode=0);

    bool RearrangeColours(Parton *, Parton *, Parton *);
    bool ArrangeColours(Parton *, Parton *, Parton *);
    int  SplitParton(Parton *, Parton *, Parton *);

    void      ExtractPartons(ATOOLS::Blob *,ATOOLS::Mass_Selector *const ms);

    void        BoostAllFS(Parton *l,Parton *r,Parton *s,bool onlyFS=false);
    void        BoostBackAllFS(Parton *l,Parton *r,Parton *s,bool onlyFS=false);

    Parton *IdParton(const size_t &id) const;

    double JetVeto(Sudakov *const sud) const;

    void Reduce();
      
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
    inline void SetShower(Shower *const shower) { p_shower=shower; }

    inline void SetDecays(const ATOOLS::DecayInfo_Vector &d) { m_decs=d; }

    inline void SetNLO(const int nlo) { m_nlo=nlo; }
    inline void SetNSkip(const int n) { m_nskip=n; }
    inline void SetNME(const int nme) { m_nme=nme; }
    inline void SetNMax(const int nmax) { m_nmax=nmax; }
    inline void SetLKF(const double &lkf) { m_lkf=lkf; } 

    inline int NLO() const { return m_nlo; }
    inline int NSkip() const { return m_nskip; }
    inline int NME() const { return m_nme; }
    inline int NMax() const { return m_nmax; }
    inline double LKF() const { return m_lkf; }

    inline void SetKtNext(const double &kt2) { m_kt2_next=kt2; } 

    inline double KtNext() const { return m_kt2_next; }

    inline void SetMuR2(const double &mur2) { m_mur2=mur2; } 

    inline double MuR2() const { return m_mur2; }

  };

  std::ostream& operator<<(std::ostream &,Singlet &);
  std::ostream& operator<<(std::ostream &,All_Singlets &); 
}

#endif
