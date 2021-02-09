#ifndef CFPSHOWER_Tools_Splitting_H
#define CFPSHOWER_Tools_Splitting_H

#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "CFPSHOWER++/Tools/Weight.H"
#include "CFPSHOWER++/Tools/Parton.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"

namespace CFPSHOWER {
  class Kernel;
  
  class Splitting {
    static size_t s_cnt;
    Parton * p_splitter, * p_spectator;
    Kernel * p_kernel;
    Weight * p_weight;

    double m_Q2, m_tstart, m_tcut;
    double m_t, m_y, m_eta, m_phi, m_zmin, m_zmax, m_z, m_kt2;
    double m_momscale;
    
    ATOOLS::Vec4D              m_cmstot, m_systot;
    std::vector<ATOOLS::Vec4D> m_moms;
    std::vector<Parton *>      m_partons;
  public:
    Splitting(Parton * splitter=NULL,Parton * spectator=NULL,
	      const double  & t=0., const double  & tcut=0.);
    ~Splitting();
    bool InitLimits();
    
    inline void SetSplitter(Parton * parton)  { p_splitter  = parton; }
    inline void SetSpectator(Parton * parton) { p_spectator = parton; }
    inline void SetKernel(Kernel * kernel)    { p_kernel    = kernel; }
    inline void SetWeight(Weight * weight)    {
      if (p_weight!=NULL) delete p_weight;
      p_weight  = weight;
    }
    inline void SetCMSTot(ATOOLS::Vec4D & cms) { m_cmstot = cms; }
    inline void SetSysTot(ATOOLS::Vec4D & sys) { m_systot = sys; }
    inline void SetMom(const size_t & pos,const ATOOLS::Vec4D & p) { m_moms[pos] = p; }
    inline void SetMomScale(const double & u)  { m_momscale = u; }
    inline void AddParton(Parton * parton)     { m_partons.push_back(parton); }
    inline void Clear()                        { m_partons.clear(); }

    inline Kernel      * GetKernel()               const { return p_kernel; }
    inline Parton      * GetSplitter()             const { return p_splitter; }
    inline Parton      * GetSpectator()            const { return p_spectator; }
    inline Weight      * GetWeight()               const { return p_weight; }
    inline Parton      * GetParton(const size_t i) const { return m_partons[i]; }
    inline ATOOLS::Vec4D GetCMSTot()               const { return m_cmstot; }
    inline ATOOLS::Vec4D GetSysTot()               const { return m_systot; }
    inline ATOOLS::Vec4D Mom(const size_t pos)     const { return m_moms[pos]; }
    
    inline void SetTstart(const double & tstart) { m_tstart = tstart; }
    inline void SetTcut(const double & tcut)     { m_tcut   = tcut; }
    inline void SetT(const double & t)           { m_t      = t; }
    inline void SetPhi(const double & phi)       { m_phi    = phi; }
    inline void SetZ(const double & z)           { m_z      = z; }
    inline void SetY(const double & y)           { m_y      = y; }
    inline void SetEta(const double & eta)       { m_eta      = eta; }
    inline const double & T()        const { return m_t; }
    inline const double & Phi()      const { return m_phi; }
    inline const double & Tstart()   const { return m_tstart; }
    inline const double & Tcut()     const { return m_tcut; }
    inline const double & Zmin()     const { return m_zmin; }
    inline const double & Zmax()     const { return m_zmax; }
    inline const double & Z()        const { return m_z; }
    inline const double & Y()        const { return m_y; }
    inline const double & Eta()      const { return m_eta; }
    inline const double & Q2()       const { return m_Q2; }
    inline const size_t   NPartons() const { return m_partons.size(); }
    inline const double & MomScale() const { return m_momscale; } 
    
    inline static size_t Count() { return s_cnt; }
  };

  std::ostream &operator<<(std::ostream &s,Splitting & split);
}// end of namespace CFPSHOWER

#endif
