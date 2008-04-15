#ifndef Flavour_H
#define Flavour_H

#define kf_code long unsigned int
#define kf_none 0

#include "Flavour_Tags.H"

#include <string> 
#include <vector>
#include <set>
#include <map>

namespace ATOOLS {

  class Flavour;

  typedef std::vector<Flavour*> PFlavour_Vector;

  struct Particle_Info {
  public:

    kf_code m_kfc;
    double  m_mass, m_yuk, m_width, m_dg, m_dm, m_qoverp2;
    int     m_icharge, m_isoweak, m_strong;
    int     m_spin, m_stable, m_masssign, m_dummy, m_majorana;
    bool    m_on, m_massive, m_hadron;

    std::string m_idname, m_texname;

    PFlavour_Vector m_content;

  public:

    // default constructor
    inline Particle_Info(): 
      m_kfc(kf_none), m_mass(0.0), m_yuk(0.0), m_width(0.0), m_dg(0.0),
      m_dm(0.0), m_qoverp2(1.0), m_icharge(0), 
      m_isoweak(0), m_strong(0), m_spin(0), m_stable(1), m_masssign(0), 
      m_dummy(1), m_majorana(0), m_on(0), m_massive(0), m_hadron(0) {}
    Particle_Info(const Particle_Info &info);
    Particle_Info(const kf_code &kfc,const double &mass,const double &width,
		  const int icharge,const int isoweak,const int strong,
		  const int spin,const int majorana,const bool on,
		  const bool stable,bool massive,const std::string &idname,
		  const std::string &texname,const bool dummy=0);
    Particle_Info(const kf_code &kfc,const double &mass,const double &width, 
		  const int icharge,const int isoweak,const int spin,
		  const bool on,const bool stable,const std::string &idname,
		  const std::string &texname);

    ~Particle_Info();

    // member functions
    bool Includes(const Flavour &fl) const;
    
    void Add(const Flavour &fl);
    void Clear();

    Flavour operator[](const size_t &i) const;

    inline size_t Size() const  { return m_content.size();   }
    inline bool   Group() const { return m_content.size()>1; }

  };// end of class Particle_Info

  typedef std::map<kf_code,Particle_Info*> KFCode_ParticleInfo_Map;

  class KF_Table: public KFCode_ParticleInfo_Map {
  public:

    ~KF_Table();

    kf_code KFFromIDName(const std::string &idname) const;
    kf_code KFFromTexName(const std::string &texname) const;

  };// end of class KF_Table

  extern KF_Table s_kftable;

  class Flavour {
  protected:

    kf_code m_kfc;
    int     m_anti;

   friend std::ostream &operator<<(std::ostream &os, const Flavour &fl);

  public:

    inline Flavour(const kf_code &kfc=kf_none,
		   const bool &anti=0): m_kfc(kfc), m_anti(0) 
    { if (anti && s_kftable[kfc]->m_majorana==0) m_anti=anti; }
    inline Flavour(const Flavour &fl): m_kfc(fl.m_kfc), m_anti(fl.m_anti) {}

    // member functions
    int  Ctq() const;
    void FromCtq(const int code);

    int HepEvt();
    void FromHepEvt(long int code);
    void FromHadron(long int code);

    std::string IDName() const;
    std::string ShellName() const;
    std::string TexName() const;
    std::string RootName() const;

    bool IsDiQuark() const;
    bool IsBaryon() const;
    bool IsB_Hadron() const;
    bool IsC_Hadron() const;

    double DiceLifeTime() const;

    // inline functions
    inline Flavour Bar() const { return Flavour(m_kfc,!m_anti); }

    inline kf_code Kfcode() const { return m_kfc; }

    inline const size_t Size() const { return s_kftable[m_kfc]->Size(); }

    inline bool Includes(const Flavour &fl) const 
    {
      if (s_kftable[m_kfc]->Group()) return s_kftable[m_kfc]->Includes(fl);
      return m_kfc==fl.m_kfc && m_anti==fl.m_anti;
    }

    inline Flavour operator[](const size_t &i) const  
    { 
      if (!s_kftable[m_kfc]->Group()) return *this;
      return m_anti?(*s_kftable[m_kfc])[i].Bar():(*s_kftable[m_kfc])[i];
    }

    inline operator long int() const 
    { return m_anti?-(long int)m_kfc:(long int)m_kfc; }

    inline Flavour &operator=(const Flavour& fl) 
    { if (this!=&fl) { m_kfc=fl.m_kfc; m_anti=fl.m_anti; } return *this; }

    inline bool IsAnti() const { return m_anti; }

    inline void SetIntCharge(const int icharge) const 
    { s_kftable[m_kfc]->m_icharge=icharge; }

    inline int    IntCharge() const 
    { int iq(s_kftable[m_kfc]->m_icharge); return m_anti?-iq:iq;     }
    inline double Charge() const    
    { double c(s_kftable[m_kfc]->m_icharge/3.0); return m_anti?-c:c; }

    inline double IsoWeak() const 
    { double c(s_kftable[m_kfc]->m_isoweak/2.0); return m_anti?-c:c; }

    inline int  StrongCharge() const 
    { int c(s_kftable[m_kfc]->m_strong); return m_anti?-c:c; }
    inline bool Strong() const
    { return s_kftable[m_kfc]->m_strong!=0; }

    inline int IntSpin() const { return s_kftable[m_kfc]->m_spin;     }
    inline double Spin() const { return s_kftable[m_kfc]->m_spin/2.0; }

    inline bool SelfAnti() const { return s_kftable[m_kfc]->m_majorana!=0; }
    inline bool Majorana() const { return s_kftable[m_kfc]->m_majorana==1; }

    inline void SetOn(const bool on) const { s_kftable[m_kfc]->m_on=on; }

    inline bool IsOn() const { return s_kftable[m_kfc]->m_on; }

    inline void SetStable(const int stable) const 
    { s_kftable[m_kfc]->m_stable=stable; }

    inline int  Stable() const   { return s_kftable[m_kfc]->m_stable;   }
    inline bool IsStable() const { return s_kftable[m_kfc]->m_stable&1; }

    inline void SetMassOn(const bool on) const    
    { s_kftable[m_kfc]->m_massive=on; }

    inline void SetMass(const double &mass) const 
    { s_kftable[m_kfc]->m_mass=mass;  }
    inline void SetYuk(const double &yuk) const   
    { s_kftable[m_kfc]->m_yuk=yuk;    }

    inline void SetMassSign(const int ms) const 
    { s_kftable[m_kfc]->m_masssign=ms; }

    inline bool IsMassive() const 
    { return s_kftable[m_kfc]->m_mass?s_kftable[m_kfc]->m_massive:0; }

    inline double Mass() const    
    { return s_kftable[m_kfc]->m_massive?s_kftable[m_kfc]->m_mass:0.0; }
    inline double SelMass() const 
    { return s_kftable[m_kfc]->m_massive&&!IsKK()?s_kftable[m_kfc]->m_mass:0.0; }
    inline double PSMass() const  
    { return s_kftable[m_kfc]->m_mass; }
    inline double Yuk() const     
    { return s_kftable[m_kfc]->m_yuk; }
    inline double DeltaGamma() const
    { return s_kftable[m_kfc]->m_dg; }
    inline void SetDeltaGamma(double dgamma) const
    {  s_kftable[m_kfc]->m_dg = dgamma; }
    inline double DeltaM() const
    {  return s_kftable[m_kfc]->m_dm; }
    inline void SetDeltaM(double dm) const
    {  s_kftable[m_kfc]->m_dm = dm; }
    inline double QOverP2() const
    {  return s_kftable[m_kfc]->m_qoverp2; }
    inline void SetQOverP2(double qoverp2) const
    {  s_kftable[m_kfc]->m_qoverp2 = qoverp2; }
    inline int MassSign() const { return s_kftable[m_kfc]->m_masssign; }

    inline void SetWidth(const double &width) const 
    { s_kftable[m_kfc]->m_width=width; }

    inline double Width() const 
    { return s_kftable[m_kfc]->m_width; }

    inline bool IsHadron() const { return s_kftable[m_kfc]->m_hadron; }

    inline bool IsFermion() const { return IntSpin()==1;   }
    inline bool IsBoson() const   { return IntSpin()%2==0; }
    inline bool IsScalar() const  { return IntSpin()==0;   }
    inline bool IsVector() const  { return IntSpin()==2;   }
    inline bool IsRaritaSchwinger() const { return IntSpin()==3; }
    inline bool IsTensor() const  { return IntSpin()==4;   }

    inline int LeptonFamily() const 
    { if (IsLepton()) return (m_kfc-9)/2; return 0; }
    inline int QuarkFamily() const 
    { if (IsQuark()) return (m_kfc+1)/2; return 0; }

    inline int LeptonNumber() 
    { if (IsLepton()||IsSlepton()||IsSneutrino()) return m_anti?-1:1; return 0; }
    inline double BaryonNumber() 
    { if (IsQuark()||IsSquark()) return m_anti?-1./3.:1./3.; return 0.; }

    inline bool IsPhoton() const { return m_kfc==kf_photon;   }
    inline bool IsLepton() const { return m_kfc>10&&m_kfc<19; }

    inline bool IsQuark() const { return m_kfc<10;                           }
    inline bool IsGluon() const { return m_kfc==kf_gluon||m_kfc==kf_shgluon; }
    inline bool IsJet() const   { return m_kfc==kf_jet;                      }

    inline bool IsChargino() const 
    { return m_kfc==kf_Chargino1||m_kfc==kf_Chargino2; }
    inline bool IsNeutralino() const 
    { return m_kfc==kf_Neutralino1||m_kfc==kf_Neutralino2||
	m_kfc==kf_Neutralino3||m_kfc==kf_Neutralino4; }
    inline bool IsSlepton() const 
    { return (m_kfc>1000010&&m_kfc<1000017)||(m_kfc>2000010&&m_kfc<2000017); }
    inline bool IsSneutrino() const 
    { return m_kfc>1000010&&m_kfc<1000017&&m_kfc%2==0; }

    inline bool IsSquark() const 
    { return Strong()&&IntSpin()==0&&!Majorana(); }
    inline bool IsGluino() const 
    { return m_kfc==kf_Gluino; }

    inline bool IsIno() const    
    { return IsGluino()||IsNeutralino()||IsChargino(); }

    inline bool IsUptype() const   { return s_kftable[m_kfc]->m_isoweak==1;  }
    inline bool IsDowntype() const { return s_kftable[m_kfc]->m_isoweak==-1; }

    inline const bool IsSusy() const 
    { return 1000000<m_kfc&&m_kfc<3000000; }

    inline const bool IsKK() const 
    { if (m_kfc==kf_graviton || m_kfc==kf_gscalar) return 1;
      return 5000000<m_kfc && m_kfc<7000000; }
    inline const int KKGeneration() const 
    { if (!IsKK()) return 0; return (m_kfc-1000000*(m_kfc/1000000))/100000; }
    inline bool Is5VDummy() const 
    { return m_kfc==kf_shgluon; }   

  };// end of class Flavour

  class Flavour_Sorting_Criterion {
  public :
    bool operator()(const Flavour &fl1,const Flavour &fl2) 
    {
      kf_code kf1(fl1.Kfcode()), kf2(fl2.Kfcode());
      if (kf1>kf2) return true;
      if (kf1<kf2) return false;
      return !(fl1.IsAnti()&&!fl2.IsAnti());
    }
  };// end of class Flavour_Sorting_Criterion

  typedef std::set<Flavour,Flavour_Sorting_Criterion> FlavourSet;
  typedef FlavourSet::iterator FlSetIter;
  typedef FlavourSet::const_iterator FlSetConstIter;

}// end of namespace ATOOLS

#endif
