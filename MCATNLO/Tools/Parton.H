#ifndef MCATNLO_Tools_Parton_H
#define MCATNLO_Tools_Parton_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Org/Message.H"
#include <list>

namespace MCATNLO {

  class Singlet;
  class Splitting_Function_Base;

  struct pst {
    enum code {
      IS   = -1,
      FS   = 1,
      none = 0
    };
  };

  class Parton;

  struct Color_Info {
    int m_i[2], m_j[2], m_k[2], m_new;
    inline Color_Info(): m_new(0)
    { m_i[0]=m_i[1]=m_j[0]=m_j[1]=m_k[0]=m_k[1]=0; }
    inline Color_Info(ATOOLS::ColorID ci,ATOOLS::ColorID cj,
		      ATOOLS::ColorID ck,const int nc=-1): m_new(nc)
    {
      m_i[0]=ci.m_i; m_i[1]=ci.m_j;
      m_j[0]=cj.m_i; m_j[1]=cj.m_j;
      m_k[0]=ck.m_i; m_k[1]=ck.m_j;
    }
  };// end of struct Color_Info

  std::ostream &operator<<(std::ostream &str,const Color_Info &ci);

  class Parton {
  private:
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_mom;
    ATOOLS::Flow    m_flow,m_meflow;
    pst::code       m_pst;
    int             m_beam,m_kin,m_col;
    double          m_kt_start, m_kt_veto, m_kt_test;
    double          m_kt_max, m_phi;
    double          m_z_test, m_y_test, m_xBj;
    Parton        * p_spect;
    Singlet       * p_sing;
    size_t          m_id, m_idx;
    ATOOLS::Poincare_Sequence m_lt;

    std::vector<Parton*> m_specs;
    Splitting_Function_Base *p_sf;

    Color_Info m_ci;

  public:
    inline Parton();
    inline Parton(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,pst::code=pst::none);
    inline Parton(const ATOOLS::Particle *,pst::code=pst::none);

    inline std::vector<Parton*> &Specs() { return m_specs; }

    inline void SetLT(const ATOOLS::Poincare_Sequence &lt) { m_lt=lt; }
    inline const ATOOLS::Poincare_Sequence &LT() const { return m_lt; }

    inline ATOOLS::Flavour   const GetFlavour()  const;
    inline ATOOLS::Vec4D     const Momentum()    const;
    inline ATOOLS::Vec4D           Momentum();
    inline pst::code               GetType()     const;
    inline int                     GetFlow(const int) const;
    inline int                     GetMEFlow(const int) const;
    inline int                     GetRFlow(const int) const;
    inline int                     Beam()        const;
    inline int                     Kin()         const { return m_kin; }
    inline int                     Col()         const { return m_col; }
    inline Color_Info             &Color() { return m_ci; }
    inline size_t                  Id()          const;
    inline size_t                  Idx()         const { return m_idx; }
    inline Splitting_Function_Base *SF()         const { return p_sf; }
    inline double                  KtMax()       const;
    inline double                  KtStart()     const;
    inline double                  KtVeto()      const;
    inline double                  KtTest()      const;
    inline double                  ZTest()       const;
    inline double                  YTest()       const;
    inline double                  Xbj()         const;
    inline double                  Phi()         const;
    inline Parton *                GetSpect();
    inline Singlet*                GetSing();
    inline void SetFlavour(const ATOOLS::Flavour &);
    inline void SetMomentum(const ATOOLS::Vec4D &);
    inline void SetStart(const double);
    inline void SetKtMax(const double);
    inline void SetVeto(const double);
    inline void SetKtTest(const double);
    inline void SetZTest(const double);
    inline void SetYTest(const double);
    inline void SetTest(const double,const double,const double,const double);
    inline void SetXbj(const double);
    inline void SetPhi(const double);
    inline void SetSpect(Parton *);
    inline void SetSing(Singlet *);
    inline void SetFlow(int code_index,int code=0);
    inline void SetMEFlow(int code_index,int code=0);
    inline void SetBeam(int);
    inline void SetKin(int kin) { m_kin=kin; }
    inline void SetCol(int col) { m_col=col; }
    inline void SetId(size_t);
    inline void SetIdx(size_t idx) { m_idx=idx; }
    inline void SetSF(Splitting_Function_Base *const sf) { p_sf=sf; }
    friend std::ostream& operator<<(std::ostream &,const Parton &);
  };

  Parton::Parton() : 
    m_flav(ATOOLS::Flavour(kf_none)),
    m_mom(ATOOLS::Vec4D(0.,0.,0.,0.)), 
    m_pst(pst::none), m_beam(0), m_kin(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_kt_max(0.), m_phi(2.0*M_PI),
    m_z_test(1.), m_y_test(1.), m_xBj(1.),
    p_spect(NULL), p_sing(NULL), m_id(0), m_idx(0)
  { 
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Flavour & flav,const ATOOLS::Vec4D & mom,pst::code pst) : 
    m_flav(flav), m_mom(mom), m_pst(pst), m_beam(0), m_kin(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_kt_max(0.), m_phi(2.0*M_PI),
    m_z_test(1.), m_y_test(1.), m_xBj(1.),
    p_spect(NULL), p_sing(NULL), m_id(0), m_idx(0)
  {
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Particle * part,pst::code pst) : 
    m_flav(part->Flav()), m_mom(part->Momentum()), 
    m_pst(pst), m_beam(0), m_kin(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_kt_max(0.), m_phi(0.),
    m_z_test(1.), m_y_test(1.), m_xBj(1.),
    p_spect(NULL), p_sing(NULL), m_id(0), m_idx(0)
  {
    if (m_pst==pst::FS) {
      m_flow.SetCode(1,part->GetFlow(1));
      m_flow.SetCode(2,part->GetFlow(2));
    }
    else {
      m_flow.SetCode(1,part->GetFlow(2));
      m_flow.SetCode(2,part->GetFlow(1));
    }
  }
  ATOOLS::Flavour   const Parton::GetFlavour()  const { return m_flav; }
  ATOOLS::Vec4D     const Parton::Momentum()    const { return m_mom; }
  ATOOLS::Vec4D           Parton::Momentum()          { return m_mom; }
  pst::code               Parton::GetType()     const { return m_pst; }
  int                     Parton::GetFlow(const int index) const {
    return m_flow.Code(index);
  }
  int                     Parton::GetMEFlow(const int index) const {
    return m_meflow.Code(index);
  }
  int                     Parton::GetRFlow(const int index) const {
    if (m_meflow.Code(index)>0) return m_meflow.Code(index); 
    return m_flow.Code(index);
  }
  double                  Parton::KtStart()     const { return m_kt_start; }
  double                  Parton::KtMax()       const { return m_kt_max; }
  double                  Parton::KtVeto()      const { return m_kt_veto; }
  double                  Parton::KtTest()      const { return m_kt_test; }
  double                  Parton::ZTest()       const { return m_z_test; }
  double                  Parton::YTest()       const { return m_y_test; }
  double                  Parton::Xbj()         const { return m_xBj; }
  double                  Parton::Phi()         const { return m_phi; }
  int                     Parton::Beam()        const { return m_beam; }
  size_t                  Parton::Id()          const { return m_id; }
  Parton *                Parton::GetSpect()          { return p_spect; }
  Singlet*                Parton::GetSing()           { return p_sing; }
  void Parton::SetFlavour(const ATOOLS::Flavour & fl) { m_flav      = fl; }
  void Parton::SetMomentum(const ATOOLS::Vec4D & mom) { m_mom      = mom; }
  void Parton::SetStart(const double kt)              { m_kt_start = kt; }
  void Parton::SetKtMax(const double kt)              { m_kt_max   = kt; }
  void Parton::SetVeto(const double kt)               { m_kt_veto  = kt; }
  void Parton::SetKtTest(const double kt)             { m_kt_test  = kt; }
  void Parton::SetZTest(const double z)               { m_z_test   = z; }
  void Parton::SetYTest(const double y)               { m_y_test   = y; }  
  void Parton::SetTest(const double kt,
		       const double z, 
		       const double y,
		       const double phi)              { m_kt_test  = kt; 
                                                        m_z_test   = z; 
                                                        m_y_test   = y; 
                                                        m_phi      = phi; }
  void Parton::SetXbj(const double x)                 { m_xBj      = x; }
  void Parton::SetPhi(const double phi)               { m_phi      = phi; }
  void Parton::SetSpect(Parton * part)                { p_spect    = part; }  
  void Parton::SetSing(Singlet* sing)                 { p_sing    = sing; }  
  void Parton::SetFlow(const int index, const int code) {
    if ((!m_flav.IsDiQuark()) && (!m_flav.Strong())) return;
    m_flow.SetCode(index,code);
  }
  void Parton::SetMEFlow(const int index, const int code) {
    if ((!m_flav.IsDiQuark()) && (!m_flav.Strong())) return;
    m_meflow.SetCode(index,code);
  }
  void Parton::SetBeam(const int beam)                { m_beam = beam;}          
  void Parton::SetId(const size_t id)                 { m_id = id;}          
  typedef std::list<Parton *>   Parton_List;
  typedef Parton_List::iterator PLiter;

}

#endif
