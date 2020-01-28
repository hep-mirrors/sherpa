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

  struct splitting_mode {
    enum code {
      diff = 0,
      coll = 1,
      soft = 2
    };
  };
  std::ostream &operator<<(std::ostream &s,const splitting_mode::code & mode);

  
  class Splitting {
    static size_t s_cnt;
    Parton * p_splitter, * p_spectator;
    Kernel * p_kernel;
    splitting_mode::code m_mode;

    double   m_msplit2, m_mspect2;
    double   m_Q2, m_Q2red, m_eta;
    double   m_tstart, m_tcut;
    double   m_t, m_t2, m_z, m_z2, m_phi, m_phi2, m_y;
    bool     m_ismassive, m_isend;
    
    Weight  * p_weight;
    size_t                m_nflavs;
    std::vector<Parton *> m_outs;
    std::vector<double>   m_m2s;
    int       m_kinscheme, m_clustered;

    PHASIC::Kin_Args          m_kinargs;
    ATOOLS::Poincare_Sequence m_seq;
  public:
    Splitting(Parton * splitter=NULL,Parton * spectator=NULL,
	      const double  & t=0., const double  & tcut=0.);
    ~Splitting();
    void InitSplitting(const ATOOLS::Mass_Selector * ms);
    void UpdateSpectatorMomentum();
    
    inline void SetSplitter(Parton * parton)       { p_splitter  = parton; }
    inline void SetSpectator(Parton * parton)      { p_spectator = parton; }
    inline void SetKernel(Kernel * kernel)         { p_kernel    = kernel; }
    inline void SetMode(splitting_mode::code mode) { m_mode   = mode; }
    inline void Set_eta(const double & eta)        { m_eta    = eta; }
    inline void Set_tstart(const double & tstart)  { m_tstart = tstart; }
    inline void Set_tcut(const double & tcut)      { m_tcut   = tcut; }
    inline void Set_t(const double & t)            { m_t      = t; }
    inline void Set_t2(const double & t2)          { m_t2     = t2; }
    inline void Set_z(const double & z)            { m_z      = z; }
    inline void Set_z2(const double & z2)          { m_z2     = z2; }
    inline void Set_phi(const double & phi=-1.) {
      m_phi = (phi<0)?2.*M_PI*ATOOLS::ran->Get():phi;
    }
    inline void Set_phi2(const double & phi2=-1.) {
      m_phi2 = (phi2<0)?2.*M_PI*ATOOLS::ran->Get():phi2;
    }
    inline void Set_y(const double & y)           { m_y      = y; }
    inline void SetEndPoint(const double & isend) { m_isend  = isend; }
    
    inline void SetWeight(Weight * weight)                 { p_weight = weight; }
    inline void SetParton(size_t i,Parton * const parton)  { m_outs[i] = parton; }
    inline void SetKinScheme(const int & kin)   { m_kinscheme = kin; }
    inline void SetClustered(const int & clu)   { m_clustered = clu; }
    inline void SetSequence(const ATOOLS::Poincare_Sequence & seq) { m_seq = seq; }
    
    inline Parton       * GetSplitter()           const { return p_splitter; }
    inline Parton       * GetSpectator()          const { return p_spectator; }
    inline Kernel       * GetKernel()             const { return p_kernel; }
    inline const splitting_mode::code & Mode()    const { return m_mode; } 
    inline const bool HasMode(splitting_mode::code mode) const { return m_mode & mode; }
    inline const bool   & IsMassive()             const { return m_ismassive; }
    inline const double & msplit2()               const { return m_msplit2; }
    inline const double & mspect2()               const { return m_mspect2; }
    inline const double & m2(const size_t & pos)  const { return m_m2s[pos]; }
    inline const double & Q2()                    const { return m_Q2; }
    inline const double & Q2red()                 const { return m_Q2red; }
    inline const double & eta()                   const { return m_eta; }
    inline const double & tstart()                const { return m_tstart; }
    inline const double & tcut()                  const { return m_tcut; }
    inline const double & t()                     const { return m_t; }
    inline const double & t2()                    const { return m_t2; }
    inline const double & z()                     const { return m_z; }
    inline const double & z2()                    const { return m_z2; }
    inline const double & phi()                   const { return m_phi; }
    inline const double & phi2()                  const { return m_phi2; }
    inline const double & y()                     const { return m_y; }
    inline const bool   & IsEndPoint()            const { return m_isend; }

    inline Weight           * GetWeight()         const { return p_weight; }
    inline const size_t       NPartons()          const { return m_nflavs; }
    inline Parton           * GetParton(size_t i) const { return m_outs[i]; }
    inline const int        & KinScheme()         const { return m_kinscheme; }
    inline const int        & Clustered()         const { return m_clustered; }
    inline PHASIC::Kin_Args * const GetKinematics()     { return &m_kinargs; }
    inline const ATOOLS::Poincare_Sequence & GetSequence() const { return m_seq; }   

    inline static size_t Count() { return s_cnt; }
  };

  std::ostream &operator<<(std::ostream &s,const Splitting & split);
}// end of namespace CFPSHOWER

#endif
