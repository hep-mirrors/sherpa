#ifndef Hadron_Part_H
#define Hadron_Part_H

#include "Transitions.H"
#include "Cluster.H"


namespace AHADIC {
  struct chtrans {
    enum code {
      HH_only      = 0,
      CH_incl      = 1,
      CH_plus_Anni = 2
    };
  };

  struct hadsel {
    enum code {
      keep    = 0,
      newpair = 1
    };
  };

  class Hadron_Part {
  protected:
    chtrans::code m_cht;
    hadsel::code  m_hadsel;
    
    All_Single_Transitions * p_stransitions;
    All_Double_Transitions * p_dtransitions;
  public:
    Hadron_Part();
    ~Hadron_Part() {}
    virtual void RedoDecay(Cluster *,Part_List *,int,
			   ATOOLS::Flavour &,ATOOLS::Flavour &)=0;
    virtual bool ForcedDecay(Cluster *,Part_List *)=0;
  };

  class Isotropic : public Hadron_Part {
  private:
    void TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    void CheckDecayKinematics(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    bool CheckDecayKinematics(Cluster *,Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
  public:
    Isotropic();
    ~Isotropic() {}
    void RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &);
    bool ForcedDecay(Cluster *,Part_List *);
  };

  class Retain : public Hadron_Part {
  public:
    Retain() {}
    ~Retain() {}
    void RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &) {}
    bool ForcedDecay(Cluster *,Part_List *) {}
  };
}
#endif
