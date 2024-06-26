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

  typedef std::vector<Flavour*> PFlavour_Vector;

  struct Particle_Info {
  public:
    // masses, yukawa's etc. in GeV, radius in fm
    kf_code m_kfc;
    double  m_mass, m_hmass, m_radius, m_yuk, m_width, m_dg, m_dm, m_qoverp2;
    int     m_icharge, m_strong, m_resummed, m_priority;
    int     m_spin, m_stable, m_masssign, m_dummy, m_majorana, m_formfactor;
    bool    m_on, m_massive, m_hadron, m_isgroup;

    std::string m_idname, m_antiname, m_texname, m_antitexname;

    PFlavour_Vector m_content;

  public:

    // default constructor
    inline Particle_Info():
      m_kfc(kf_none), m_mass(0.0), m_radius(0.0), m_hmass(0.0), m_yuk(-1.0), m_width(0.0),
      m_dg(0.0), m_dm(0.0), m_qoverp2(1.0), m_icharge(0),
      m_strong(0), m_resummed(0), m_priority(0), m_spin(0), m_stable(1),
      m_masssign(0), m_dummy(1), m_majorana(0), m_formfactor(0), m_on(0),
      m_massive(0), m_hadron(0), m_isgroup(0) {}
    Particle_Info(const Particle_Info &info);
    Particle_Info(const kf_code &kfc,const double &mass,const double &radius,const double &width,
      const int icharge,const int strong,
      const int spin,const int majorana,const bool on,
      const int stable,bool massive,const std::string &idname,
      const std::string &antiname, const std::string &texname,
      const std::string &antitexname ,const bool dummy=0,
      const bool isgroup=false);
    Particle_Info(const kf_code &kfc,const double &mass,const double &radius,const double &width,
      const int icharge,const int spin,
      const bool on,const int stable,const std::string &idname,
      const std::string &texname);
    Particle_Info(const kf_code &kfc,const double &mass,const double &radius,const int icharge,
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

  class Flavour {

  private:

    void InitializeParticleInfo(const kf_code);

  protected:

    Particle_Info *p_info;

    int m_anti;

    friend std::ostream &operator<<(std::ostream &os, const Flavour &fl);

  public:

    Flavour(Particle_Info&, bool anti=0);
    Flavour(long int kfc=kf_none);
    Flavour(kf_code kfc, bool anti);
    Flavour(const Flavour& fl);

    // member functions
    std::string IDName() const;
    std::string ShellName() const;
    std::string LegacyShellName() const;
    std::string TexName() const;
    std::string RootName() const;

    bool IsDiQuark() const;
    bool IsBaryon() const;
    bool IsMeson() const;
    bool IsNucleon() const;
    bool IsB_Hadron() const;
    bool IsC_Hadron() const;

    double GenerateLifeTime() const;
    double RelBWMass(const double& min, const double& max,
                     double peak=-1.0, double width=-1.0) const;




    static double ISSymmetryFactor(const std::vector<Flavour>& flavs);
    static double FSSymmetryFactor(const std::vector<Flavour>& flavs);

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

    inline bool operator==(const Flavour &fl) const
    { return p_info==fl.p_info && m_anti==fl.m_anti; }

    inline operator long int() const
    { return m_anti?-(long int)Kfcode():(long int)Kfcode(); }

    inline Flavour &operator=(const Flavour& fl)
    { if (this!=&fl) { p_info=fl.p_info; m_anti=fl.m_anti; } return *this; }

    inline bool IsAnti() const { return m_anti; }

    inline int    IntCharge() const
    { int iq(p_info->m_icharge); return m_anti?-iq:iq;     }
    inline double Charge() const
    { double c(p_info->m_icharge/3.0); return m_anti?-c:c; }
    inline bool IsQED() const
    { return (IsPhoton()||IntCharge());}

    inline double IsoWeak() const
    { double c((p_info->m_kfc<7 || (p_info->m_kfc>10 && p_info->m_kfc<17))
	       ?((p_info->m_kfc%2)?-0.5:0.5):0); return m_anti?-c:c; }

    Flavour IsoWeakPartner() const;
    Flavour GoldstoneBosonPartner() const;

    inline int  StrongCharge() const
    { int c(p_info->m_strong); return m_anti?-c:c; }
    inline bool Strong() const
    { return p_info->m_strong!=0&&!IsDiQuark(); }
    inline bool IsQCD() const
    { return p_info->m_strong; }

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
    { p_info->m_mass=dabs(mass);
      p_info->m_masssign = mass < 0.0 ? -1 : 1; }

    inline void SetHadMass(const double &hmass) const
    { p_info->m_hmass=hmass;  }

    inline bool IsMassive() const
    { return p_info->m_mass?p_info->m_massive:0; }

    inline void SetRadius(const double & radius)
    { p_info->m_radius=radius; }

    inline double Mass(const bool set=0) const
    { return set||p_info->m_massive?p_info->m_mass:0.0; }
    inline double SelMass() const
    { return p_info->m_massive&&!IsKK()?p_info->m_mass:0.0; }
    inline double Radius() const
    { return p_info->m_radius ;}
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

    inline bool IsPhoton() const { return Kfcode()==kf_photon;   }
    inline bool IsLepton() const { return Kfcode()>10&&Kfcode()<19; }
    inline bool IsChargedLepton() const { return IsLepton()&&IsDowntype(); }
    inline bool IsNeutrino()      const { return IsLepton()&&IsUptype(); }

    inline bool IsQuark() const { return Kfcode()<10; }
    inline bool IsGluon() const
    { return Kfcode()==kf_gluon||Kfcode()==kf_shgluon; }
    inline bool IsJet() const   { return Kfcode()==kf_jet; }

    inline bool IsUptype() const   { return IsAnti()?IsoWeak()<0:IsoWeak()>0; }
    inline bool IsDowntype() const { return IsAnti()?IsoWeak()>0:IsoWeak()<0; }

    inline bool IsKK() const
    { if (Kfcode()==kf_graviton || Kfcode()==kf_gscalar) return 1;
      return 0; }
    inline int KKGeneration() const
    { if (!IsKK()) return 0;
      return (Kfcode()-1000000*(Kfcode()/1000000))/100000; }

    inline bool IsDummy() const
    { return p_info->m_dummy; }

    inline bool operator<(const Flavour &f) const
    {
      if (Kfcode()<f.Kfcode()) return true;
      if (Kfcode()>f.Kfcode()) return false;
      return m_anti<f.m_anti;
    }

  };// end of class Flavour

  std::ostream &operator<<(std::ostream &os, const Flavour &fl);

  typedef std::vector<Flavour> Flavour_Vector;
  typedef std::map<Flavour,Flavour> Flavour_Map;
  typedef std::set<Flavour> Flavour_Set;

  class Mass_Selector {
  public:
    virtual ~Mass_Selector();
    virtual double Mass(const Flavour &fl)  const = 0;
    inline  double Mass2(const Flavour &fl) const {
      double m(Mass(fl)); return m*m;
    }
  };// end of class Mass_Selector

}// end of namespace ATOOLS

#endif
