#ifndef CSSHOWER_Tools_Parton_H
#define CSSHOWER_Tools_Parton_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"

namespace CSSHOWER {

  class Singlet;
  class Splitting_Function_Base;

  struct pst {
    enum code {
      IS   = -1,
      FS   = 1,
      none = 0
    };
  };

  struct Sudakov_Reweighting_Info {
    bool accepted;
    double scale;
    double accwgt;
    double lastj;
    double lastcpl;
    Splitting_Function_Base* sf;
    double x, y, z;
    ATOOLS::Flavour flspec;
  };
  typedef std::vector<Sudakov_Reweighting_Info> Sudakov_Reweighting_Infos;

  class Parton {
  private:
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_mom;
    ATOOLS::Flow    m_flow;
    pst::code       m_pst;
    bool            m_forcedsplitting;
    int             m_beam,m_stat,m_kin,m_kscheme,m_col,m_nforced;
    double          m_kt_start, m_kt_veto, m_kt_test;
    double          m_z_test, m_y_test, m_phi, m_weight, m_xBj, m_t;
    double          m_kt_soft[2];
    bool            m_fromdec;
    Parton        * p_left, * p_right;
    Parton        * p_spect, * p_prev, * p_next;
    Singlet       * p_sing;
    size_t          m_id;
    ATOOLS::Poincare_Sequence m_lt;

    std::vector<std::pair<double,double> > m_weights;
    Sudakov_Reweighting_Infos m_sudakovrewinfos;

  public:
    inline Parton();
    inline Parton(const ATOOLS::Flavour &,const ATOOLS::Vec4D &,pst::code=pst::none);
    inline Parton(const ATOOLS::Particle *,pst::code=pst::none);

    void DeleteAll();

    Parton *FollowUp();
    void    UpdateDaughters();
    void    UpdateColours(int newr,int newa);

    inline void SetKScheme(const int kscheme) { m_kscheme=kscheme; }
    inline int KScheme() const { return m_kscheme; }

    inline void SetLT(const ATOOLS::Poincare_Sequence &lt) { m_lt=lt; }
    inline const ATOOLS::Poincare_Sequence &LT() const { return m_lt; }

    inline std::vector<std::pair<double,double> > &Weights() { return m_weights; }
    Sudakov_Reweighting_Infos& SudakovReweightingInfos() { return m_sudakovrewinfos; }

    double Weight(const double &scale);

    inline ATOOLS::Flavour   const GetFlavour()   const;
    inline ATOOLS::Vec4D     const Momentum()     const;
    inline ATOOLS::Vec4D           Momentum();
    inline pst::code               GetType()      const;
    inline int                     GetFlow(const int) const;
    inline int                     Beam()         const;
    inline int                     Stat()         const;
    inline int                     Kin()          const { return m_kin; }
    inline int                     Col()          const { return m_col; }
    inline size_t                  Id()           const;
    inline double                  KtStart()      const;
    inline double                  KtSoft(int i)  const;
    inline double                  KtVeto()       const;
    inline double                  KtTest()       const;
    inline double                  ZTest()        const;
    inline double                  YTest()        const;
    inline double                  Weight()       const;
    inline double                  Xbj()          const;
    inline double                  Mass2()        const { return m_t; }
    inline double                  Phi()          const;
    inline bool                    ForcedSplitting() const;
    inline size_t                  ForcedTrials() const { return m_nforced; }
    inline Parton *                GetLeft();
    inline Parton *                GetRight();
    inline Parton *                GetSpect();
    inline Parton *                GetPrev();
    inline Parton *                GetNext();
    inline Singlet*                GetSing();

    inline void SetFlavour(const ATOOLS::Flavour &);
    inline void SetMomentum(const ATOOLS::Vec4D &);
    inline void SetStart(const double);
    inline void SetSoft(int i,const double);
    inline void SetVeto(const double);
    inline void SetKtTest(const double);
    inline void SetKtStart(const double);
    inline void SetZTest(const double);
    inline void SetYTest(const double);
    inline void SetWeight(const double);
    inline void SetTest(const double,const double,const double,const double);
    inline void SetXbj(const double);
    inline void SetMass2(const double t) { m_t=t; }
    inline void SetPhi(const double);
    inline void SetForcedSplitting(const bool);
    inline void IncForcedTrials() { m_nforced++; }
    void SetLeftOf(Parton *);
    void SetRightOf(Parton *);
    inline void SetLeft(Parton *);
    inline void SetRight(Parton *);
    inline void SetSpect(Parton *);
    inline void SetPrev(Parton *);
    inline void SetNext(Parton *);
    inline void SetSing(Singlet *);
    inline void SetFlow(int code_index,int code=0);
    inline void SetBeam(int);
    inline void SetStat(int);
    inline void SetKin(int kin) { m_kin=kin; }
    inline void SetCol(int col) { m_col=col; }
    inline void SetId(size_t);
    inline void SetFromDec(const bool &val){m_fromdec=val;}
    inline bool FromDec() const {return m_fromdec;}
    friend std::ostream& operator<<(std::ostream &,const Parton &);
  };

  Parton::Parton() : 
    m_flav(ATOOLS::Flavour(kf_none)),
    m_mom(ATOOLS::Vec4D(0.,0.,0.,0.)), 
    m_pst(pst::none), m_forcedsplitting(false),
    m_beam(0), m_stat(0), m_kin(0), m_kscheme(0), m_col(0), m_nforced(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_z_test(1.), m_y_test(1.), m_phi(2.0*M_PI), m_weight(1.0), m_xBj(1.),
    p_left(NULL), p_right(NULL), p_spect(NULL), p_prev(NULL), p_next(NULL), p_sing(NULL), 
    m_id(0),m_fromdec(false)
  { 
    m_kt_soft[0]=m_kt_soft[1]=std::numeric_limits<double>::max();
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Flavour & flav,const ATOOLS::Vec4D & mom,pst::code pst) : 
    m_flav(flav), m_mom(mom), m_pst(pst), m_forcedsplitting(false),
    m_beam(0), m_stat(0), m_kin(0), m_kscheme(0), m_col(0), m_nforced(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_z_test(1.), m_y_test(1.), m_phi(2.0*M_PI), m_weight(1.0), m_xBj(1.),
    p_left(NULL), p_right(NULL), p_spect(NULL), p_prev(NULL), p_next(NULL), p_sing(NULL), 
    m_id(0), m_fromdec(false)
  {
    m_kt_soft[0]=m_kt_soft[1]=std::numeric_limits<double>::max();
    m_flow.SetCode(1,0);
    m_flow.SetCode(2,0);
  }

  Parton::Parton(const ATOOLS::Particle * part,pst::code pst) : 
    m_flav(part->Flav()), m_mom(part->Momentum()), 
    m_pst(pst), m_forcedsplitting(false),
    m_beam(0), m_stat(0), m_kin(0), m_kscheme(0), m_col(0), m_nforced(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_z_test(1.), m_y_test(1.), m_phi(2.0*M_PI), m_weight(1.0), m_xBj(1.),
    p_left(NULL), p_right(NULL), p_spect(NULL), p_prev(NULL), p_next(NULL), p_sing(NULL), 
    m_id(0), m_fromdec(part->FromDec())
  {
    m_kt_soft[0]=m_kt_soft[1]=std::numeric_limits<double>::max();
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
  double                  Parton::KtStart()     const { return m_kt_start; }
  double                  Parton::KtSoft(int i) const { return m_kt_soft[i]; }
  double                  Parton::KtVeto()      const { return m_kt_veto; }
  double                  Parton::KtTest()      const { return m_kt_test; }
  double                  Parton::ZTest()       const { return m_z_test; }
  double                  Parton::YTest()       const { return m_y_test; }
  double                  Parton::Weight()      const { return m_weight; }
  double                  Parton::Xbj()         const { return m_xBj; }
  double                  Parton::Phi()         const { return m_phi; }
  bool                    Parton::ForcedSplitting() const { return m_forcedsplitting; }
  int                     Parton::Beam()        const { return m_beam; }
  int                     Parton::Stat()        const { return m_stat; }
  size_t                  Parton::Id()          const { return m_id; }
  Parton *                Parton::GetLeft()           { return p_left;  }
  Parton *                Parton::GetRight()          { return p_right; }
  Parton *                Parton::GetSpect()          { return p_spect; }
  Parton *                Parton::GetPrev()           { return p_prev; }
  Parton *                Parton::GetNext()           { return p_next; }
  Singlet*                Parton::GetSing()           { return p_sing; }
  void Parton::SetFlavour(const ATOOLS::Flavour & fl) { m_flav      = fl; }
  void Parton::SetMomentum(const ATOOLS::Vec4D & mom) { m_mom      = mom; }
  void Parton::SetStart(const double kt)              { m_kt_start = kt; }
  void Parton::SetSoft(int i,const double kt)         { m_kt_soft[i] = kt; }
  void Parton::SetVeto(const double kt)               { m_kt_veto  = kt; }
  void Parton::SetKtTest(const double kt)             { m_kt_test  = kt; }
  void Parton::SetKtStart(const double kt)            { m_kt_start = kt; }
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
  void Parton::SetForcedSplitting(const bool fs)      { m_forcedsplitting = fs; }
  void Parton::SetLeft(Parton * part)                 { p_left    = part; }  
  void Parton::SetRight(Parton * part)                { p_right   = part; }  
  void Parton::SetSpect(Parton * part)                { p_spect    = part; }  
  void Parton::SetPrev(Parton * part)                 { p_prev    = part; }  
  void Parton::SetNext(Parton * part)                 { p_next    = part; }  
  void Parton::SetSing(Singlet* sing)                 { p_sing    = sing; }  
  void Parton::SetFlow(const int index, const int code) {
    m_flow.SetCode(index,code);
  }
  void Parton::SetBeam(const int beam)                { m_beam = beam;}          
  void Parton::SetStat(const int stat)                { m_stat = stat;}          
  void Parton::SetId(const size_t id)                 { m_id = id;}          
  typedef std::list<Parton *>   Parton_List;
  typedef Parton_List::iterator PLiter;

}

#endif
