#ifndef CFPSHOWER_Tools_Splitting_H
#define CFPSHOWER_Tools_Splitting_H

#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "CFPSHOWER++/Tools/Weight.H"
#include "CFPSHOWER++/Tools/Parton.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS {
  class Mass_Selector;
}
  
namespace CFPSHOWER {
  class Kernel;
  
  class Splitting {
    static size_t    s_cnt;
    Parton         * p_splitter, * p_spectator;
    Kernel         * p_kernel;
    Weight         * p_weight;
    double           m_t, m_t0, m_Q2, m_z, m_phi, m_sijk, m_x, m_y, m_eta;
    double           m_mij2, m_mi2, m_mj2, m_mk2;
    ATOOLS::Vec4D    m_moms[3], m_specmom;
    ATOOLS::Poincare_Sequence m_seq;
    Color            m_cols[3];
    Parton *         p_outs[3];
    int              m_kinscheme, m_clustered;
    PHASIC::Kin_Args m_ff;
    void SetMasses(const ATOOLS::Mass_Selector * ms);
  public:
    Splitting(Parton * splitter=NULL,Parton * spectator=NULL,
	      const double  & t=0., const double  & t0=0.);
    inline ~Splitting() {
      if (p_weight) delete p_weight;
      s_cnt--;
    }
    void InitKinematics(const ATOOLS::Mass_Selector * ms);
    
    inline void SetKinScheme(const int & kin)   { m_kinscheme = kin; }
    inline const int & KinScheme() const        { return m_kinscheme; }
    inline void SetClustered(const int & clu)   { m_clustered = clu; }
    inline const int & Clustered() const        { return m_clustered; }
    inline void SetKernel(Kernel * kernel)      { p_kernel    = kernel; }
    inline void SetWeight(Weight * weight)      { p_weight    = weight; }

    inline void SetT(const double & t)          { m_t = t; }
    inline void SetT0(const double & t0)        { m_t0 = t0; }
    inline void SetZ(const double & z=-1.)      {
      if (z>=0.) m_z = z;
      else m_z = ATOOLS::ran->Get();
    }
    inline void SetX(const double & x)           { m_x   = x; }
    inline void SetY(const double & y)           { m_y   = y; }
    inline void SetQ2() {
      if (p_splitter!=NULL && p_spectator!=NULL) {
	m_sijk = ((p_splitter->Beam()>0?1.:-1.)*p_splitter->Mom() +
		  (p_spectator->Beam()>0?1.:-1.)*p_spectator->Mom()).Abs2();
	m_Q2 = ATOOLS::dabs(m_sijk);
      }
    }
    inline void SetQ2(const double & Q2)         { m_Q2  = Q2; }
    inline void SetSijk(const double & sijk)     { m_sijk = sijk; }
    inline void SetEta(const double & eta)       { m_eta = eta; }
    inline void Setmij2(const double & mij2)     { m_mij2 = mij2; }
    inline void Setmi2(const double & mi2)       { m_mi2 = mi2; }
    inline void Setmj2(const double & mj2)       { m_mj2 = mj2; }
    inline void Setmk2(const double & mk2)       { m_mk2 = mk2; }
    inline void Setphi(const double & phi=-1.) {
      if (phi>=0.) m_phi = phi;
      else m_phi = 2.*M_PI*ATOOLS::ran->Get();
    }
    inline void SetSequence(const ATOOLS::Poincare_Sequence & seq) { m_seq = seq; }
    
    inline void SetSplitter(Parton * parton)               { p_splitter  = parton; }
    inline void SetSpectator(Parton * parton)              { p_spectator = parton; }
    inline void UpdateSpectatorMomentum()                  { p_spectator->SetMom(m_specmom); }
    inline void SetParton(size_t i,Parton * const parton)  { p_outs[i] = parton; }
    inline void SetSpecMom(const ATOOLS::Vec4D & mom)      { m_specmom = mom; }
    inline void SetMom(size_t i,const ATOOLS::Vec4D & mom) { m_moms[i] = mom; }
    inline void SetCol(size_t i,const Color & col)         { m_cols[i] = col; }

    inline const double & T()    const { return m_t; }
    inline const double & T0()   const { return m_t0; }
    inline const double & Z()    const { return m_z; }
    inline const double & X()    const { return m_x; }
    inline const double & Y()    const { return m_y; }
    inline const double & Q2()   const { return m_Q2; }
    inline const double & Eta()  const { return m_eta; }
    inline const double & sijk() const { return m_sijk; }
    inline const double & mij2() const { return m_mij2; }
    inline const double & mi2()  const { return m_mi2; }
    inline const double & mj2()  const { return m_mj2; }
    inline const double & mk2()  const { return m_mk2; }
    inline const double & phi()  const { return m_phi; }
    inline const ATOOLS::Vec4D & Mom(size_t i) const { return m_moms[i]; }
    inline const Color           Col(size_t i) const { return m_cols[i]; }
    
    inline const ATOOLS::Poincare_Sequence & GetSequence() const { return m_seq; }   
    inline Parton           * GetSplitter()       const { return p_splitter; }
    inline Parton           * GetSpectator()      const { return p_spectator; }
    inline Parton           * GetParton(size_t i) const { return p_outs[i]; }
    inline Kernel           * GetKernel()         const { return p_kernel; }
    inline Weight           * GetWeight()         const { return p_weight; }
    inline PHASIC::Kin_Args * GetKinematics()           { return &m_ff; }

    inline static size_t Count() { return s_cnt; }
  };

  std::ostream &operator<<(std::ostream &s,const Splitting & split);
}// end of namespace CFPSHOWER

#endif
