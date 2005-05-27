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
  public:
    Hadron_Part(All_Single_Transitions * trans) : p_stransitions(trans) {}
    ~Hadron_Part() {}
    virtual void RedoDecay(Cluster *,Part_List *,int,
			   ATOOLS::Flavour &,ATOOLS::Flavour &) {}
  };

  class Isotropic : public Hadron_Part {
  private:
    void TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    void CheckDecayKinematics(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    void CheckDecayKinematics(Cluster *,Cluster *,ATOOLS::Flavour &);
  public:
    Isotropic(All_Single_Transitions *);
    ~Isotropic() {}
    void RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &);
  };

  class Retain : public Hadron_Part {
  public:
    Retain(All_Single_Transitions * trans) : Hadron_Part(trans) {}
    ~Retain() {}
    void RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &) {}
  };
}
#endif
