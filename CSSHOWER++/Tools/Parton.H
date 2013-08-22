#ifndef CSSHOWER_Tools_Parton_H
#define CSSHOWER_Tools_Parton_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
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
    ATOOLS::Vec4D   m_mom, m_oldmom, m_fixspec;
    ATOOLS::Flow    m_flow,m_meflow;
    pst::code       m_pst;
    int             m_beam,m_stat,m_kin,m_kscheme,m_col;
    double          m_kt_start, m_kt_veto, m_kt_test, m_kt_max;
    double          m_z_test, m_y_test, m_phi, m_weight, m_xBj, m_t;
    Parton        * p_left, * p_right;
    Parton        * p_spect, * p_prev, * p_next;
    Singlet       * p_sing;
    size_t          m_id;
    bool            m_connected;
    ATOOLS::Poincare_Sequence m_lt;
    bool            m_invs;
    Parton        * p_invleft, * p_invright;

    std::vector<std::pair<double,double> > m_weights;

  public:
    inline Parton();
    inline Parton(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,pst::code=pst::none);
    inline Parton(const ATOOLS::Particle *,pst::code=pst::none);

    void DeleteAll();

    Parton *FollowUp();
    bool    Splits();
    void    UpdateDaughters();
    void    UpdateNewDaughters(Parton *ref);
    void    UpdateColours();

    inline void SetKScheme(const int kscheme) { m_kscheme=kscheme; }
    inline int KScheme() const { return m_kscheme; }

    inline void SetLT(const ATOOLS::Poincare_Sequence &lt) { m_lt=lt; }
    inline const ATOOLS::Poincare_Sequence &LT() const { return m_lt; }

    inline std::vector<std::pair<double,double> > &Weights() { return m_weights; }

    double Weight(const double &scale);

    inline ATOOLS::Flavour   const GetFlavour()  const;
    inline ATOOLS::Vec4D     const Momentum()    const;
    inline ATOOLS::Vec4D     const OldMomentum() const { return m_oldmom; }
    inline ATOOLS::Vec4D     const FixSpec()     const { return m_fixspec; }
    inline ATOOLS::Vec4D           Momentum();
    inline pst::code               GetType()     const;
    inline int                     GetFlow(const int) const;
    inline int                     GetMEFlow(const int) const;
    inline int                     GetRFlow(const int) const;
    inline int                     Beam()        const;
    inline int                     Stat()        const;
    inline int                     Kin()         const { return m_kin; }
    inline int                     Col()         const { return m_col; }
    inline size_t                  Id()          const;
    inline double                  KtMax()       const;
    inline double                  KtStart()     const;
    inline double                  KtVeto()      const;
    inline double                  KtTest()      const;
    inline double                  ZTest()       const;
    inline double                  YTest()       const;
    inline double                  Weight()        const;
    inline double                  Xbj()         const;
    inline double                  Mass2()       const { return m_t; }
    inline double                  Phi()         const;
    inline Parton *                GetLeft();
    inline Parton *                GetRight();
    inline Parton *                GetSpect();
    inline Parton *                GetPrev();
    inline Parton *                GetNext();
    inline Parton *                GetInvLeft();
    inline Parton *                GetInvRight();
    inline Singlet*                GetSing();
    inline bool                    Connected()   const { return m_connected; }
    inline bool                    InvColours()  const { return m_invs; }

    inline void SetInvColours(const bool tag) { m_invs = tag; }
    inline void SetFlavour(const ATOOLS::Flavour &);
    inline void SetMomentum(const ATOOLS::Vec4D &);
    inline void SetOldMomentum(const ATOOLS::Vec4D &p) { m_oldmom=p; }
    inline void SetFixSpec(const ATOOLS::Vec4D &p) { m_fixspec=p; }
    inline void SetStart(const double);
    inline void SetKtMax(const double);
    inline void SetVeto(const double);
    inline void SetKtTest(const double);
    inline void SetZTest(const double);
    inline void SetYTest(const double);
    inline void SetWeight(const double);
    inline void SetTest(const double,const double,const double,const double);
    inline void SetXbj(const double);
    inline void SetMass2(const double t) { m_t=t; }
    inline void SetPhi(const double);
    void SetLeftOf(Parton *);
    void SetRightOf(Parton *);
    inline void SetLeft(Parton *);
    inline void SetRight(Parton *);
    inline void SetInvLeft(Parton *);
    inline void SetInvRight(Parton *);
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
    inline void SetCol(int col) { m_col=col; }
    inline void SetId(size_t);
    inline void SetConnected(const bool & c) { m_connected = c; }
    friend std::ostream& operator<<(std::ostream &,const Parton &);
  };

  Parton::Parton() : 
    m_flav(ATOOLS::Flavour(kf_none)),
    m_mom(ATOOLS::Vec4D(0.,0.,0.,0.)), 
    m_pst(pst::none), m_beam(0), m_stat(0), m_kin(0), m_kscheme(0), m_col(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.), m_kt_max(0.),
    m_z_test(1.), m_y_test(1.), m_phi(2.0*M_PI), m_weight(1.0), m_xBj(1.),
    p_left(NULL), p_right(NULL), p_spect(NULL), p_prev(NULL), p_next(NULL), 
    p_sing(NULL), 
    m_id(0), m_connected(true),
    m_invs(false), p_invleft(NULL), p_invright(NULL)
  { 
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Flavour & flav,const ATOOLS::Vec4D & mom,pst::code pst) : 
    m_flav(flav), m_mom(mom), m_pst(pst), m_beam(0), m_stat(0), m_kin(0), m_kscheme(0), m_col(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.), m_kt_max(0.),
    m_z_test(1.), m_y_test(1.), m_phi(2.0*M_PI), m_weight(1.0), m_xBj(1.),
    p_left(NULL), p_right(NULL), p_spect(NULL), p_prev(NULL), p_next(NULL), 
    p_sing(NULL), 
    m_id(0), m_connected(true),
    m_invs(false), p_invleft(NULL), p_invright(NULL)

  {
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Particle * part,pst::code pst) : 
    m_flav(part->Flav()), m_mom(part->Momentum()), 
    m_pst(pst), m_beam(0), m_stat(0), m_kin(0), m_kscheme(0), m_col(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.), m_kt_max(0.),
    m_z_test(1.), m_y_test(1.), m_phi(2.0*M_PI), m_weight(1.0), m_xBj(1.),
    p_left(NULL), p_right(NULL), p_spect(NULL), p_prev(NULL), p_next(NULL), 
    p_sing(NULL), 
    m_id(0), m_connected(true),
    m_invs(false), p_invleft(NULL), p_invright(NULL)

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
  double                  Parton::Weight()      const { return m_weight; }
  double                  Parton::Xbj()         const { return m_xBj; }
  double                  Parton::Phi()         const { return m_phi; }
  int                     Parton::Beam()        const { return m_beam; }
  int                     Parton::Stat()        const { return m_stat; }
  size_t                  Parton::Id()          const { return m_id; }
  Parton *                Parton::GetLeft()           { return p_left;  }
  Parton *                Parton::GetRight()          { return p_right; }
  Parton *                Parton::GetInvLeft()        { return p_invleft;  }
  Parton *                Parton::GetInvRight()       { return p_invright; }
  Parton *                Parton::GetSpect()          { return p_spect; }
  Parton *                Parton::GetPrev()           { return p_prev; }
  Parton *                Parton::GetNext()           { return p_next; }
  Singlet*                Parton::GetSing()           { return p_sing; }
  void Parton::SetFlavour(const ATOOLS::Flavour & fl) { m_flav      = fl; }
  void Parton::SetMomentum(const ATOOLS::Vec4D & mom) { m_mom      = mom; }
  void Parton::SetStart(const double kt)              { m_kt_start = kt; }
  void Parton::SetKtMax(const double kt)              { m_kt_max   = kt; }
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
  void Parton::SetLeft(Parton * part)  { 
    p_left = part;  
    //if (m_invs) p_invright = part;
  }  
  void Parton::SetRight(Parton * part) { 
    p_right = part; 
    //if (m_invs) p_invleft = part;
  }  
  void Parton::SetInvLeft(Parton * part)              { p_invleft  = part; }  
  void Parton::SetInvRight(Parton * part)             { p_invright = part; }  
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

}

#endif
