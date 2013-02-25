#ifndef ATOOLS_Phys_Flavour_H
#define ATOOLS_Phys_Flavour_H

#define kf_code long unsigned int
#define kf_none 0

#include "ATOOLS/Phys/Flavour_Tags.H"
#include "ATOOLS/Math/MathTools.H"

#include <string> 
#include <vector>
#include <set>
#include <iostream>
#include <map>

namespace ATOOLS {

  class Flavour;
  class Mass_Selector;
  
  void OutputHadrons(std::ostream&);
  void OutputParticles(std::ostream&);
  void OutputContainers(std::ostream&);
  
  typedef std::vector<Flavour*> PFlavour_Vector;

  struct Particle_Info {
  public:

    kf_code m_kfc;
    double  m_mass, m_hmass, m_yuk, m_width, m_dg, m_dm, m_qoverp2;
    int     m_icharge, m_isoweak, m_strong, m_resummed, m_priority;
    int     m_spin, m_stable, m_masssign, m_dummy, m_majorana, m_formfactor;
    bool    m_on, m_massive, m_hadron, m_isgroup;

    std::string m_idname, m_texname;

    PFlavour_Vector m_content;

  public:

    // default constructor
    inline Particle_Info(): 
      m_kfc(kf_none), m_mass(0.0), m_hmass(0.0), m_yuk(-1.0), m_width(0.0),
      m_dg(0.0), m_dm(0.0), m_qoverp2(1.0), m_icharge(0), m_isoweak(0),
      m_strong(0), m_resummed(0), m_priority(0), m_spin(0), m_stable(1),
      m_masssign(0), m_dummy(1), m_majorana(0), m_formfactor(0), m_on(0),
      m_massive(0), m_hadron(0), m_isgroup(0) {}
    Particle_Info(const Particle_Info &info);
    Particle_Info(const kf_code &kfc,const double &mass,const double &width,
		  const int icharge,const int isoweak,const int strong,
		  const int spin,const int majorana,const bool on,
		  const int stable,bool massive,const std::string &idname,
		  const std::string &texname,const bool dummy=0, 
		  const bool isgroup=false);
    Particle_Info(const kf_code &kfc,const double &mass,const double &width, 
		  const int icharge,const int isoweak,const int spin,
		  const bool on,const int stable,const std::string &idname,
		  const std::string &texname);
    Particle_Info(const kf_code &kfc,const double &mass,const int icharge,
                  const int spin,const int formfactor,const std::string &idname,
                  const std::string &texname);

    ~Particle_Info();

    // member functions
    bool Includes(const Flavour &fl) const;
    
    void Add(const Flavour &fl);
    void Clear();

    Flavour operator[](const size_t &i) const;
    
    inline size_t Size() const  { return m_content.size();   }
    inline bool   Group() const { return m_isgroup||m_content.size()>1; }
    inline void   SetIsGroup(bool isgroup) { m_isgroup = isgroup; };

    void SetResummed();

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

    Particle_Info *p_info;

    int m_anti;

    friend std::ostream &operator<<(std::ostream &os, const Flavour &fl);

  public:

    inline Flavour(Particle_Info &info,const bool &anti=0): 
      p_info(&info), m_anti(0)
    { if (anti && p_info->m_majorana==0) m_anti=anti; }
    inline Flavour(const kf_code &kfc=kf_none,const bool &anti=0): 
      p_info(NULL), m_anti(0)
    { KFCode_ParticleInfo_Map::iterator it(s_kftable.find(kfc));
      if (it!=s_kftable.end()) p_info=it->second; else return;
      if (anti && p_info->m_majorana==0) m_anti=anti; }
    inline Flavour(const Flavour &fl): 
      p_info(fl.p_info), m_anti(fl.m_anti) {}

    // member functions
    int  Ctq() const;
    void FromCtq(const int code);

    int HepEvt() const;
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

    double GenerateLifeTime() const;
    double RelBWMass(const double& min, const double& max,
                     double peak=-1.0, double width=-1.0) const;

    // inline functions
    inline Flavour Bar() const { return Flavour(*p_info,!m_anti); }

    inline kf_code Kfcode() const { return p_info->m_kfc; }

    inline size_t Size() const { return p_info->Size(); }
    inline size_t IsGroup() const { return p_info->Group(); }

    inline bool Includes(const Flavour &fl) const 
    {
      if (p_info->Group()) return p_info->Includes(fl);
      return p_info==fl.p_info && m_anti==fl.m_anti;
    }

    inline Flavour operator[](const size_t &i) const  
    { 
      if (!p_info->Group()) return *this;
      return m_anti?(*p_info)[i].Bar():(*p_info)[i];
    }

    inline bool operator==(const Flavour &fl)
    { return p_info==fl.p_info && m_anti==fl.m_anti; }

    inline operator long int() const 
    { return m_anti?-(long int)Kfcode():(long int)Kfcode(); }

    inline Flavour &operator=(const Flavour& fl) 
    { if (this!=&fl) { p_info=fl.p_info; m_anti=fl.m_anti; } return *this; }

    inline bool IsAnti() const { return m_anti; }

    inline void SetIntCharge(const int icharge) const 
    { p_info->m_icharge=icharge; }

    inline int    IntCharge() const
    { int iq(p_info->m_icharge); return m_anti?-iq:iq;     }
    inline double Charge() const
    { double c(p_info->m_icharge/3.0); return m_anti?-c:c; }
    inline bool   Electromagnetic() const
    { return p_info->m_icharge!=0&&!IsDiQuark(); }

    inline double IsoWeak() const 
    { double c(p_info->m_isoweak/2.0); return m_anti?-c:c; }
    inline bool Weak() const 
    { return p_info->m_isoweak!=0&&!IsDiQuark(); }

    inline int  StrongCharge() const 
    { int c(p_info->m_strong); return m_anti?-c:c; }
    inline bool Strong() const
    { return p_info->m_strong!=0&&!IsDiQuark(); }

    inline bool Resummed() const
    { return p_info->m_resummed; }

    inline int Priority() const
    { return p_info->m_priority; }

    inline int IntSpin() const { return p_info->m_spin;     }
    inline double Spin() const { return p_info->m_spin/2.0; }

    inline bool SelfAnti() const { return p_info->m_majorana!=0; }
    inline bool Majorana() const { return p_info->m_majorana==1; }

    inline int FormFactor() const {return p_info->m_formfactor; }

    inline bool IsIon() const { return p_info->m_kfc>1000000000; }

    inline int GetAtomicNumber() const { return (p_info->m_kfc/10)%1000; }

    inline void SetOn(const bool on) const { p_info->m_on=on; }

    inline bool IsOn() const { return p_info->m_on; }

    inline void SetStable(const int stable) const 
    { p_info->m_stable=stable; }

    inline int  Stable() const   { return p_info->m_stable;   }
    bool IsStable() const;

    inline void SetMassOn(const bool on) const    
    { p_info->m_massive=on; }

    inline void SetMass(const double &mass) const 
    { p_info->m_mass=mass;  }
    inline void SetHadMass(const double &hmass) const 
    { p_info->m_hmass=hmass;  }

    inline void SetMassSign(const int ms) const 
    { p_info->m_masssign=ms; }

    inline bool IsMassive() const 
    { return p_info->m_mass?p_info->m_massive:0; }

    inline double Mass(const bool set=0) const    
    { return set||p_info->m_massive?p_info->m_mass:0.0; }
    inline double SelMass() const 
    { return p_info->m_massive&&!IsKK()?p_info->m_mass:0.0; }
    inline double HadMass() const  
    { return p_info->m_hmass; }
    inline double Yuk() const     
    { return p_info->m_yuk>=0.0 ? p_info->m_yuk : (p_info->m_massive ? p_info->m_mass : 0.0); }
    inline double DeltaGamma() const
    { return p_info->m_dg; }
    inline void SetDeltaGamma(double dgamma) const
    {  p_info->m_dg = dgamma; }
    inline double DeltaM() const
    {  return p_info->m_dm; }
    inline void SetDeltaM(double dm) const
    {  p_info->m_dm = dm; }
    inline double QOverP2() const
    {  return p_info->m_qoverp2; }
    inline void SetQOverP2(double qoverp2) const
    {  p_info->m_qoverp2 = qoverp2; }
    inline int MassSign() const { return p_info->m_masssign; }

    inline void SetWidth(const double &width) const 
    { p_info->m_width=width; }

    inline double Width() const 
    { return p_info->m_width; }

    inline bool IsHadron() const { return p_info->m_hadron; }

    inline bool IsFermion() const { return IntSpin()==1;   }
    inline bool IsBoson() const   { return IntSpin()%2==0; }
    inline bool IsScalar() const  { return IntSpin()==0;   }
    inline bool IsVector() const  { return IntSpin()==2;   }
    inline bool IsRaritaSchwinger() const { return IntSpin()==3; }
    inline bool IsTensor() const  { return IntSpin()==4;   }

    inline int LeptonFamily() const 
    { if (IsLepton()) return (Kfcode()-9)/2; return 0; }
    inline int QuarkFamily() const 
    { if (IsQuark()) return (Kfcode()+1)/2; return 0; }

    inline int LeptonNumber() 
    { if (IsLepton()||IsSlepton()||IsSneutrino()) return m_anti?-1:1; return 0; }
    inline double BaryonNumber() 
    { if (IsQuark()||IsSquark()) return m_anti?-1./3.:1./3.; return 0.; 
    if (abs(StrongCharge())==3) return 1./double(StrongCharge()); }

    inline bool IsPhoton() const { return Kfcode()==kf_photon;   }
    inline bool IsLepton() const { return Kfcode()>10&&Kfcode()<19; }

    inline bool IsQuark() const { return Kfcode()<10; }
    inline bool IsGluon() const 
    { return Kfcode()==kf_gluon||Kfcode()==kf_shgluon; }
    inline bool IsJet() const   { return Kfcode()==kf_jet; }

    inline bool IsChargino() const 
    { return (Kfcode()==kf_Chargino1||Kfcode()==kf_Chargino2) && IntSpin()==1; }
    inline bool IsNeutralino() const 
    { return (Kfcode()==kf_Neutralino1||Kfcode()==kf_Neutralino2||
	      Kfcode()==kf_Neutralino3||Kfcode()==kf_Neutralino4) && IntSpin()==1; }
    inline bool IsSlepton() const 
    { return ((Kfcode()>1000010&&Kfcode()<1000017)||
	(Kfcode()>2000010&&Kfcode()<2000017)) && IntSpin()==0; }
    inline bool IsSneutrino() const 
    { return Kfcode()>1000010&&Kfcode()<1000017&&Kfcode()%2==0&&IntSpin()==0; }

    inline bool IsSquark() const 
    { return Strong()&&(StrongCharge()==3 || StrongCharge()==-3)&&IntSpin()==0&&!Majorana(); }
    inline bool IsGluino() const 
    { return Kfcode()==kf_Gluino; }

    inline bool IsIno() const    
    { return IsGluino()||IsNeutralino()||IsChargino(); }

    inline bool IsUptype() const   { return p_info->m_isoweak==1;  }
    inline bool IsDowntype() const { return p_info->m_isoweak==-1; }

    inline bool IsSusy() const 
    { return 1000000<Kfcode()&&Kfcode()<3000000; }

    inline bool IsKK() const 
    { if (Kfcode()==kf_graviton || Kfcode()==kf_gscalar) return 1;
      return 0; }
    inline int KKGeneration() const 
    { if (!IsKK()) return 0; 
      return (Kfcode()-1000000*(Kfcode()/1000000))/100000; }
    inline bool Is5VDummy() const 
    { return Kfcode()==kf_shgluon; }   

    inline bool IsDummy() const 
    { return p_info->m_dummy; }   

    inline bool operator<(const Flavour &f) const 
    {
      if (Kfcode()<f.Kfcode()) return true;
      if (Kfcode()>f.Kfcode()) return false;
      return m_anti<f.m_anti;
    }
    
    static kf_code PdgToSherpa(const unsigned long& pdg);

  };// end of class Flavour

  typedef std::vector<Flavour> Flavour_Vector;
  typedef std::map<Flavour,Flavour> Flavour_Map;
  typedef std::set<Flavour> Flavour_Set;

  class Mass_Selector {
  public:

    virtual ~Mass_Selector();

    virtual double Mass(const Flavour &fl) const = 0;

    inline double Mass2(const Flavour &fl) const
    { double m(Mass(fl)); return m*m; }

  };// end of class Mass_Selector

}// end of namespace ATOOLS

#endif
