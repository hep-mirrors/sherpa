#ifndef CFPSHOWER_Tools_Splitting_H
#define CFPSHOWER_Tools_Splitting_H

#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "CFPSHOWER++/Tools/Weight.H"
#include "CFPSHOWER++/Tools/Parton.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include <algorithm>

namespace CFPSHOWER {
  class Kernel;
  
  class Splitting {
    static size_t s_cnt;
    Parton * p_splitter, * p_spectator;
    Kernel * p_kernel;
    Weight * p_weight;

    double m_Q, m_Q2, m_tstart, m_tcut;
    double m_t, m_y, m_eta, m_phi, m_zmin, m_zmax, m_z, m_kt2;
    double m_IntPhi;
    bool   m_ismassive;
    size_t m_nout;
    double        m_msplit, m_msplit2, m_mspect, m_mspect2, m_masssum;
    ATOOLS::Vec4D m_psplit, m_pspect, m_pdipole, m_systot, m_kinspect;
    std::vector<double>        m_masses, m_masses2;
    std::vector<ATOOLS::Vec4D> m_moms, m_allmoms;
    std::vector<Parton *>      m_partons;
    //std::map<Parton *,int>     m_updatemap;
  public:
    Splitting(Parton * splitter=NULL,Parton * spectator=NULL,
	      const double  & t=0., const double  & tcut=0.);
    ~Splitting();
    void Init(Kernel * kernel,const ATOOLS::Mass_Selector * msel);
    
    inline void SetSplitter(Parton * parton)  { p_splitter  = parton; }
    inline void SetSpectator(Parton * parton) { p_spectator = parton; }
    inline void SetWeight(Weight * weight)    {
      if (p_weight!=NULL) delete p_weight;
      p_weight  = weight;
    }
    inline void SetCMSTot(ATOOLS::Vec4D cms)   { m_pdipole  = cms; }
    inline void SetSysTot(ATOOLS::Vec4D sys)   { m_systot   = sys; }
    inline void SetKinSpect(ATOOLS::Vec4D kin) { m_kinspect = kin; }
    inline void AddParton(Parton * part)       { m_partons.push_back(part); }
    //inline void AddUpdate(Parton * part,int pos) { m_updatemap[part] = pos; }
    inline void Clear()                          { m_partons.clear(); }
    inline void SetMom(const size_t & pos,
		       const ATOOLS::Vec4D & p)  { m_moms[pos] = p; }
    inline void SetAllMoms(std::vector<ATOOLS::Vec4D> & moms) {
      m_allmoms.clear();
      std::copy(moms.begin(), moms.end(), std::back_inserter(m_allmoms));
    }

    inline Kernel      * GetKernel()               const { return p_kernel; }
    inline Parton      * GetSplitter()             const { return p_splitter; }
    inline Parton      * GetSpectator()            const { return p_spectator; }
    inline Weight      * GetWeight()               const { return p_weight; }
    inline Parton      * GetParton(const size_t i) const { return m_partons[i]; }
    //inline int UpdatePos(Parton * part) {
    //  if (m_updatemap.find(part)==m_updatemap.end()) return -1;
    //  return m_updatemap[part];
    //}
    inline ATOOLS::Vec4D GetCMSTot()               const { return m_pdipole; }
    inline ATOOLS::Vec4D GetSysTot()               const { return m_systot; }
    inline ATOOLS::Vec4D GetKinSpect()             const { return m_kinspect; }
    inline ATOOLS::Vec4D Mom(const size_t pos)     const { return m_moms[pos]; }
    inline size_t        NumberP()                 const { return m_allmoms.size(); }
    inline ATOOLS::Vec4D P(const size_t pos)       const { return m_allmoms[pos]; }
    inline double        Mass(const size_t pos)    const { return m_masses[pos]; }
    inline double        Mass2(const size_t pos)   const { return m_masses2[pos]; }
    inline double        MassSum()                 const { return m_masssum; }
    inline double        Mass2Splitter()           const { return m_msplit2;} 
    inline double        Mass2Spectator()          const { return m_mspect2;} 
    
    inline void SetTstart(const double & tstart) { m_tstart = tstart; }
    inline void SetTcut(const double & tcut)     { m_tcut   = tcut; }
    inline void SetT(const double & t)           { m_t      = t; }
    inline void SetPhi(const double & phi)       { m_phi    = phi; }
    inline void SetZmin(const double & zmin)     { m_zmin   = zmin; }
    inline void SetZmax(const double & zmax)     { m_zmax   = zmax; }
    inline void SetZ(const double & z)           { m_z      = z; }
    inline void SetY(const double & y)           { m_y      = y; }
    inline void SetKT2(const double & kt2)       { m_kt2    = kt2; }
    inline void SetEta(const double & eta)       { m_eta    = eta; }
    inline void SetIntPhi(const double & IntPhi) { m_IntPhi = IntPhi; }
    inline const double & T()        const { return m_t; }
    inline const double & Phi()      const { return m_phi; }
    inline const double & Tstart()   const { return m_tstart; }
    inline const double & Tcut()     const { return m_tcut; }
    inline const double & Zmin()     const { return m_zmin; }
    inline const double & Zmax()     const { return m_zmax; }
    inline const double & Z()        const { return m_z; }
    inline const double & Y()        const { return m_y; }
    inline const double & Eta()      const { return m_eta; }
    inline const double & KT2()      const { return m_kt2; }
    inline const double & Q2()       const { return m_Q2; }
    inline const double & IntPhi()   const { return m_IntPhi; }
    inline const size_t   NPartons() const { return m_partons.size(); }
    inline const size_t   NOut()     const { return m_nout; }
    inline static size_t Count()           { return s_cnt; }
  };

  std::ostream &operator<<(std::ostream &s,Splitting & split);
}// end of namespace CFPSHOWER

#endif
