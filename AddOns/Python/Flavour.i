//%module Flavour
%include "std_string.i"
%{

#define kf_code long unsigned int
#define kf_none 0

#include <ATOOLS/Phys/Flavour_Tags.H>
#include <ATOOLS/Math/MathTools.H>
#include <ATOOLS/Phys/Flavour.H>
#include <ATOOLS/Org/MyStrStream.H>
#include <string> 
#include <vector>
#include <set>
#include <iostream>
#include <map>

using namespace ATOOLS;
#define kf_code long unsigned int

%}


// also tell SWIG about the definition of kf_code
#define kf_code long unsigned int

namespace ATOOLS {

  class Flavour;
  class Mass_Selector;
  
  class Flavour {
  public:

    inline Flavour(const kf_code &kfc=kf_none,const bool &anti=0): 
      p_info(NULL), m_anti(0)
    { KFCode_ParticleInfo_Map::iterator it(s_kftable.find(kfc));
      if (it!=s_kftable.end()) p_info=it->second; else return;
      if (anti && p_info->m_majorana==0) m_anti=anti; }

    // member functions
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

    inline bool operator==(const Flavour &fl)
    { return p_info==fl.p_info && m_anti==fl.m_anti; }

    inline bool IsAnti() const { return m_anti; }

    inline int    IntCharge() const
    { int iq(p_info->m_icharge); return m_anti?-iq:iq;     }
    inline double Charge() const
    { double c(p_info->m_icharge/3.0); return m_anti?-c:c; }

    inline int  StrongCharge() const 
    { int c(p_info->m_strong); return m_anti?-c:c; }
    inline bool Strong() const
    { return p_info->m_strong!=0&&!IsDiQuark(); }

    inline bool Resummed() const
    { return p_info->m_resummed; }

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

    inline bool IsPhoton() const { return Kfcode()==kf_photon;   }
    inline bool IsLepton() const { return Kfcode()>10&&Kfcode()<19; }

    inline bool IsQuark() const { return Kfcode()<10; }
    inline bool IsGluon() const 
    { return Kfcode()==kf_gluon||Kfcode()==kf_shgluon; }
    inline bool IsJet() const   { return Kfcode()==kf_jet; }

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
    
    %extend {
      std::string __str__() {
	MyStrStream conv;
	conv<<*self;
	return conv.str();
      };
    };

  };// end of class Flavour

 }// end of namespace ATOOLS

