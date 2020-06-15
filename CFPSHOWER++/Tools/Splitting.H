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
    double   m_tstart, m_tcut, m_zmin, m_zmax;
    double   m_t[2], m_z[2], m_ztilde[3], m_phi[2], m_y;
    bool     m_ismassive, m_isend;
    
    Weight      * p_weight;
    size_t        m_nflavs;
    Parton      * m_outs[3];
    double        m_m2s[3];
    ATOOLS::Vec4D m_moms[3], m_specmom; 
    int           m_kinscheme, m_clustered;

    ATOOLS::Poincare_Sequence m_seq;
    PHASIC::Kin_Args          m_kinargs;
  public:
    Splitting(Parton * splitter=NULL,Parton * spectator=NULL,
	      const double  & t=0., const double  & tcut=0.);
    ~Splitting();
    void InitSplitting(const ATOOLS::Mass_Selector * ms);
    bool InitLimits();
    
    inline void SetSplitter(Parton * parton)       { p_splitter  = parton; }
    inline void SetSpectator(Parton * parton)      { p_spectator = parton; }
    inline void SetKernel(Kernel * kernel)         { p_kernel    = kernel; }
    inline void SetMode(splitting_mode::code mode) { m_mode   = mode; }
    inline void Set_eta(const double & eta)        { m_eta    = eta; }
    inline void Set_tstart(const double & tstart)  { m_tstart = tstart; }
    inline void Set_tcut(const double & tcut)      { m_tcut   = tcut; }
    inline void Set_t(const size_t & pos,const double & t)         { m_t[pos]      = t; }
    inline void Set_sai(const double & sai)                        { m_t[1]        = sai; }
    inline void Set_zmin(const double & zmin=0.)                   { m_zmin        = zmin; }
    inline void Set_zmax(const double & zmax=1.)                   { m_zmax        = zmax; }
    inline void Set_z(const size_t & pos,const double & z=0.)      { m_z[pos]      = z; }
    inline void Set_ztilde(const size_t & pos,const double & z=0.) { m_ztilde[pos] = z; }
    inline void Set_phi(const size_t & pos,const double & phi=-1.) {
      m_phi[pos] = (phi<0)?2.*M_PI*ATOOLS::ran->Get():phi;
    }
    inline void Set_y(const double & y)           { m_y      = y; }
    inline void SetEndPoint(const double & isend) { m_isend  = isend; }
    
    inline void SetWeight(Weight * weight)  {
      if (p_weight!=NULL) delete p_weight;
      p_weight  = weight;
    }
    inline void SetParton(size_t i,Parton * const parton)        { m_outs[i] = parton; }
    inline void SetMomentum(size_t i,ATOOLS::Vec4D const & mom)  { m_moms[i] = mom; }
    inline void SetSpectatorMomentum(ATOOLS::Vec4D const & mom)  { m_specmom = mom; }
    inline void SetKinScheme(const int & kin)   { m_kinscheme = kin; }
    inline void SetClustered(const int & clu)   { m_clustered = clu; }
    inline void SetSequence(const ATOOLS::Poincare_Sequence & seq) { m_seq = seq; }
    
    inline Parton       * GetSplitter()             const { return p_splitter; }
    inline Parton       * GetSpectator()            const { return p_spectator; }
    inline Kernel       * GetKernel()               const { return p_kernel; }
    inline const splitting_mode::code & Mode()      const { return m_mode; } 
    inline const bool HasMode(splitting_mode::code mode) const { return m_mode & mode; }
    inline const bool   & IsMassive()                const { return m_ismassive; }
    inline const double & msplit2()                  const { return m_msplit2; }
    inline const double & mspect2()                  const { return m_mspect2; }
    inline const double & m2(const size_t & pos)     const { return m_m2s[pos]; }
    inline const double & Q2()                       const { return m_Q2; }
    inline const double & Q2red()                    const { return m_Q2red; }
    inline const double & eta()                      const { return m_eta; }
    inline const double & tstart()                   const { return m_tstart; }
    inline const double & tcut()                     const { return m_tcut; }
    inline const double & t(const size_t & pos)      const { return m_t[pos]; }
    inline const double & zmin()                     const { return m_zmin; }
    inline const double & zmax()                     const { return m_zmax; }
    inline const double & z(const size_t & pos)      const { return m_z[pos]; }
    inline const double & ztilde(const size_t & pos) const { return m_ztilde[pos]; }
    inline const double & phi(const size_t & pos)    const { return m_phi[pos]; }
    inline const double & y()                        const { return m_y; }
    inline const double & sai()                      const { return m_t[1]; }
    inline const bool   & IsEndPoint()               const { return m_isend; }

    inline Weight              * GetWeight()         const { return p_weight; }
    inline const size_t          NPartons()          const { return m_nflavs; }
    inline Parton              * GetParton(size_t i) const { return m_outs[i]; }
    inline const ATOOLS::Vec4D & Momentum(size_t i)  const { return m_moms[i]; }
    inline const ATOOLS::Vec4D & SpectatorMomentum() const { return m_specmom; }
    inline const int           & KinScheme()         const { return m_kinscheme; }
    inline const int           & Clustered()         const { return m_clustered; }
    inline const ATOOLS::Poincare_Sequence & GetSequence()   const { return m_seq; }   
    inline const PHASIC::Kin_Args          * GetKinematics() const { return &m_kinargs; }
    
    inline static size_t Count() { return s_cnt; }
  };

  std::ostream &operator<<(std::ostream &s,const Splitting & split);
}// end of namespace CFPSHOWER

#endif
