#ifndef CSSHOWER_Tools_Parton_H
#define CSSHOWER_Tools_Parton_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include <list>

namespace CSSHOWER {

  class Singlet;

  struct pst {
    enum code {
      IS   = -1,
      FS   = 1,
      none = 0
    };
  };

  class Parton {
  private:
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_mom;
    ATOOLS::Flow    m_flow,m_meflow;
    pst::code       m_pst;
    int             m_beam,m_stat,m_kin;
    double          m_kt_start, m_kt_veto, m_kt_test;
    double          m_kt_max, m_kt_prev, m_kt_next, m_phi;
    double          m_z_test, m_y_test, m_weight, m_xBj;
    Parton        * p_left, * p_right;
    Parton        * p_spect, * p_prev, * p_next;
    Singlet       * p_sing;
    size_t          m_id;

    std::vector<std::pair<double,double> > m_weights;
  public:
    inline Parton();
    inline Parton(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,pst::code=pst::none);
    inline Parton(const ATOOLS::Particle *,pst::code=pst::none);

    void DeleteAll();

    Parton *FollowUp();
    void    UpdateDaughters();
    void    UpdateNewDaughters();
    void    UpdateColours();

    inline std::vector<std::pair<double,double> > &Weights() { return m_weights; }

    double Weight(const double &scale);

    inline ATOOLS::Flavour   const GetFlavour()  const;
    inline ATOOLS::Vec4D     const Momentum()    const;
    inline ATOOLS::Vec4D           Momentum();
    inline pst::code               GetType()     const;
    inline int                     GetFlow(const int) const;
    inline int                     GetMEFlow(const int) const;
    inline int                     GetRFlow(const int) const;
    inline int               const Beam()        const;
    inline int               const Stat()        const;
    inline int               const Kin()         const { return m_kin; }
    inline size_t            const Id()          const;
    inline double            const KtMax()       const;
    inline double            const KtPrev()      const;
    inline double            const KtNext()      const;
    inline double            const KtStart()     const;
    inline double            const KtVeto()      const;
    inline double            const KtTest()      const;
    inline double            const ZTest()       const;
    inline double            const YTest()       const;
    inline double            const Weight()        const;
    inline double            const Xbj()         const;
    inline double            const Phi()         const;
    inline Parton *                GetLeft();
    inline Parton *                GetRight();
    inline Parton *                GetSpect();
    inline Parton *                GetPrev();
    inline Parton *                GetNext();
    inline Singlet*                GetSing();
    inline void SetFlavour(const ATOOLS::Flavour &);
    inline void SetMomentum(const ATOOLS::Vec4D &);
    inline void SetStart(const double);
    inline void SetKtMax(const double);
    inline void SetKtPrev(const double);
    inline void SetKtNext(const double);
    inline void SetVeto(const double);
    inline void SetKtTest(const double);
    inline void SetZTest(const double);
    inline void SetYTest(const double);
    inline void SetWeight(const double);
    inline void SetTest(const double,const double,const double,const double);
    inline void SetXbj(const double);
    inline void SetPhi(const double);
    inline void SetLeftOf(Parton *);
    inline void SetRightOf(Parton *);
    inline void SetLeft(Parton *);
    inline void SetRight(Parton *);
    inline void SetSpect(Parton *);
    inline void SetPrev(Parton *);
    inline void SetNext(Parton *);
    inline void SetSing(Singlet *);
    inline void SetFlow(int code_index,int code=0);
    inline void SetMEFlow(int code_index,int code=0);
    inline void SetRFlow();
    inline void SetBeam(int);
    inline void SetStat(int);
    inline void SetKin(int kin) { m_kin=kin; }
    inline void SetId(size_t);
    friend std::ostream& operator<<(std::ostream &,const Parton &);
  };

  Parton::Parton() : 
    m_flav(ATOOLS::Flavour(kf_none)),
    m_mom(ATOOLS::Vec4D(0.,0.,0.,0.)), 
    m_pst(pst::none), m_beam(0), m_stat(0), m_kin(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_kt_max(0.), m_kt_prev(0.), m_kt_next(0.), m_phi(2.0*M_PI),
    m_z_test(1.), m_y_test(1.), m_weight(1.0), m_xBj(1.), p_left(NULL), p_right(NULL),
    p_spect(NULL), p_prev(NULL), p_next(NULL), p_sing(NULL), m_id(0)
  { 
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Flavour & flav,const ATOOLS::Vec4D & mom,pst::code pst) : 
    m_flav(flav), m_mom(mom), m_pst(pst), m_beam(0), m_stat(0), m_kin(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_kt_max(0.), m_kt_prev(0.), m_kt_next(0.), m_phi(2.0*M_PI),
    m_z_test(1.), m_y_test(1.), m_weight(1.0), m_xBj(1.), p_left(NULL), p_right(NULL),
    p_spect(NULL), p_prev(NULL), p_next(NULL), p_sing(NULL), m_id(0)
  {
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Particle * part,pst::code pst) : 
    m_flav(part->Flav()), m_mom(part->Momentum()), 
    m_pst(pst), m_beam(0), m_stat(0), m_kin(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_kt_max(0.), m_kt_prev(0.), m_kt_next(0.), m_phi(0.),
    m_z_test(1.), m_y_test(1.), m_weight(1.0), m_xBj(1.), p_left(NULL), p_right(NULL),
    p_spect(NULL), p_prev(NULL), p_next(NULL), p_sing(NULL), m_id(0)
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
  double            const Parton::KtStart()     const { return m_kt_start; }
  double            const Parton::KtMax()       const { return m_kt_max; }
  double            const Parton::KtPrev()      const { return m_kt_prev; }
  double            const Parton::KtNext()      const { return m_kt_next; }
  double            const Parton::KtVeto()      const { return m_kt_veto; }
  double            const Parton::KtTest()      const { return m_kt_test; }
  double            const Parton::ZTest()       const { return m_z_test; }
  double            const Parton::YTest()       const { return m_y_test; }
  double            const Parton::Weight()      const { return m_weight; }
  double            const Parton::Xbj()         const { return m_xBj; }
  double            const Parton::Phi()         const { return m_phi; }
  int               const Parton::Beam()        const { return m_beam; }
  int               const Parton::Stat()        const { return m_stat; }
  size_t            const Parton::Id()          const { return m_id; }
  Parton *                Parton::GetLeft()           { return p_left;  }
  Parton *                Parton::GetRight()          { return p_right; }
  Parton *                Parton::GetSpect()          { return p_spect; }
  Parton *                Parton::GetPrev()           { return p_prev; }
  Parton *                Parton::GetNext()           { return p_next; }
  Singlet*                Parton::GetSing()           { return p_sing; }
  void Parton::SetFlavour(const ATOOLS::Flavour & fl) { m_flav      = fl; }
  void Parton::SetMomentum(const ATOOLS::Vec4D & mom) { m_mom      = mom; }
  void Parton::SetStart(const double kt)              { m_kt_start = kt; }
  void Parton::SetKtMax(const double kt)              { m_kt_max   = kt; }
  void Parton::SetKtPrev(const double kt)             { m_kt_prev  = kt; }
  void Parton::SetKtNext(const double kt)             { m_kt_next  = kt; }
  void Parton::SetVeto(const double kt)               { m_kt_veto  = kt; }
  void Parton::SetKtTest(const double kt)             { m_kt_test  = kt; }
  void Parton::SetZTest(const double z)               { m_z_test   = z; }
  void Parton::SetYTest(const double y)               { m_y_test   = y; }  
  void Parton::SetWeight(const double weight)         { m_weight   = weight; }  
  void Parton::SetTest(const double kt,
		       const double z, 
		       const double y,
		       const double phi)              { m_kt_test  = kt; 
                                                        m_z_test   = z; 
                                                        m_y_test   = y; 
                                                        m_phi      = phi; }
  void Parton::SetXbj(const double x)                 { m_xBj      = x; }
  void Parton::SetPhi(const double phi)               { m_phi      = phi; }
  void Parton::SetLeftOf(Parton * part)
  { part->SetLeft(p_left); if (p_left) p_left->SetRight(part); }  
  void Parton::SetRightOf(Parton * part)
  { part->SetRight(p_right); if (p_right) p_right->SetLeft(part); }  
  void Parton::SetLeft(Parton * part)                 { p_left    = part; }  
  void Parton::SetRight(Parton * part)                { p_right   = part; }  
  void Parton::SetSpect(Parton * part)                { p_spect    = part; }  
  void Parton::SetPrev(Parton * part)                 { p_prev    = part; }  
  void Parton::SetNext(Parton * part)                 { p_next    = part; }  
  void Parton::SetSing(Singlet* sing)                 { p_sing    = sing; }  
  void Parton::SetFlow(const int index, const int code) {
    if ((!m_flav.IsDiQuark()) && (!m_flav.Strong())) return;
    m_flow.SetCode(index,code);
  }
  void Parton::SetMEFlow(const int index, const int code) {
    if ((!m_flav.IsDiQuark()) && (!m_flav.Strong())) return;
    m_meflow.SetCode(index,code);
  }
  void Parton::SetRFlow() {
    m_meflow.SetCode(1,m_flow.Code(1));
    m_meflow.SetCode(2,m_flow.Code(2));
  }
  void Parton::SetBeam(const int beam)                { m_beam = beam;}          
  void Parton::SetStat(const int stat)                { m_stat = stat;}          
  void Parton::SetId(const size_t id)                 { m_id = id;}          
  typedef std::list<Parton *>   Parton_List;
  typedef Parton_List::iterator PLiter;

  inline size_t IdCount(const size_t &id)
  {
    size_t ic(id), cn(0);
    for (size_t i(0);ic>0;++i) {
      size_t c(1<<i);
      if (ic&c) { ++cn;	ic-=c; }
    }
    return cn;
  }

}

#endif
