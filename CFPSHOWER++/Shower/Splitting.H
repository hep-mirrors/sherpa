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
    static size_t s_cnt;
    Parton      * p_splitter, * p_spectator;
    Kernel      * p_kernel;

    double        m_msplit2, m_mspect2, m_m2[3];
    double        m_t, m_tstart, m_t0, m_Q2, m_s;
    double        m_x, m_y, m_eta, m_phi, m_z[3], m_sij[3];
    
    Weight      * p_weight;
    Color         m_cols[3];
    ATOOLS::Vec4D m_moms[3], m_specmom;
    Parton *      p_outs[3];
    
    int           m_kinscheme, m_clustered;
    PHASIC::Kin_Args          m_kinargs;
    ATOOLS::Poincare_Sequence m_seq;
    
    void SetMasses(const ATOOLS::Mass_Selector * ms);
  public:
    Splitting(Parton * splitter=NULL,Parton * spectator=NULL,
	      const double  & t=0., const double  & t0=0.);
    inline ~Splitting() {
      if (p_weight) delete p_weight;
      s_cnt--;
    }
    void InitKinematics(const ATOOLS::Mass_Selector * ms);
    void UpdateSpectatorMomentum();
    
    inline void SetSplitter(Parton * parton)    { p_splitter  = parton; }
    inline void SetSpectator(Parton * parton)   { p_spectator = parton; }
    inline void SetKernel(Kernel * kernel)      { p_kernel    = kernel; }
    inline void Set_t(const double & t)         { m_t      = t; }
    inline void Set_tStart(const double & t)    { m_tstart = t; }
    inline void Set_t0(const double & t0)       { m_t0     = t0; }
    inline void Set_Q2(const double & Q2)       { m_Q2     = Q2; }
    inline void Set_Q2()                        {
      if (p_splitter!=NULL && p_spectator!=NULL) {
	m_s = ((p_splitter->Beam()>0?1.:-1.)*p_splitter->Mom() +
	       (p_spectator->Beam()>0?1.:-1.)*p_spectator->Mom()).Abs2();
	m_Q2 = ATOOLS::dabs(m_s);
      }
    }
    inline void Set_s(const double & s)         { m_s      = s; }
    inline void Set_x(const double & x)         { m_x      = x; }
    inline void Set_y(const double & y)         { m_y      = y; }
    inline void Set_eta(const double & eta)     { m_eta    = eta; }
    inline void Set_phi(const double & phi=-1.) {
      m_phi    = (phi<0)?2.*M_PI*ATOOLS::ran->Get():phi;
    }
    inline void Set_z(const size_t & pos,const double & z) { m_z[pos] = z; }
    inline void Set_sij(size_t & pos,const double & sij)   { m_sij[pos] = sij; }
    
    inline void SetWeight(Weight * weight)                 { p_weight    = weight; }
    inline void SetCol(size_t i,const Color & col)         { m_cols[i] = col; }
    inline void SetMom(size_t i,const ATOOLS::Vec4D & mom) { m_moms[i] = mom; }
    inline void SetSpecMom(const ATOOLS::Vec4D & mom)      { m_specmom = mom; }
    inline void SetParton(size_t i,Parton * const parton)  { p_outs[i] = parton; }
    inline void SetKinScheme(const int & kin)   { m_kinscheme = kin; }
    inline void SetClustered(const int & clu)   { m_clustered = clu; }
    inline void SetSequence(const ATOOLS::Poincare_Sequence & seq) { m_seq = seq; }
    
    inline Parton       * GetSplitter()           const { return p_splitter; }
    inline Parton       * GetSpectator()          const { return p_spectator; }
    inline Kernel       * GetKernel()             const { return p_kernel; }
    inline const double & msplit2()               const { return m_msplit2; }
    inline const double & mspect2()               const { return m_mspect2; }
    inline const double & m2(const size_t & pos)  const { return m_m2[pos]; }
    inline const double & t()                     const { return m_t; }
    inline const double & tStart()                const { return m_tstart; }
    inline const double & t0()                    const { return m_t0; }
    inline const double & Q2()                    const { return m_Q2; }
    inline const double & s()                     const { return m_s; }
    inline const double & y()                     const { return m_y; }
    inline const double & x()                     const { return m_x; }
    inline const double & eta()                   const { return m_eta; }
    inline const double & phi()                   const { return m_phi; }
    inline const double & z(const size_t & pos)   const { return m_z[pos]; }
    inline const double & sij(const size_t & pos) const { return m_sij[pos]; }
    inline const Color           Col(size_t i)    const { return m_cols[i]; }
    inline const ATOOLS::Vec4D & Mom(size_t i)    const { return m_moms[i]; }
    

    inline Weight           * GetWeight()         const { return p_weight; }
    inline Parton           * GetParton(size_t i) const { return p_outs[i]; }
    inline const int        & KinScheme()         const { return m_kinscheme; }
    inline const int        & Clustered()         const { return m_clustered; }
    inline PHASIC::Kin_Args * const GetKinematics()     { return &m_kinargs; }
    inline const ATOOLS::Poincare_Sequence & GetSequence() const { return m_seq; }   

    inline static size_t Count() { return s_cnt; }
  };

  std::ostream &operator<<(std::ostream &s,const Splitting & split);
}// end of namespace CFPSHOWER

#endif
