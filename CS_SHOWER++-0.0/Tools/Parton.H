#ifndef Parton_h
#define Parton_h

#include "Flavour.H"
#include "Particle.H"
#include "Vector.H"
#include "Message.H"
#include <list>

namespace CS_SHOWER {
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
    pst::code       m_pst;
    int             m_cpartner;
    double          m_kt_start, m_kt_veto, m_kt_test;
    double          m_z_test, m_y_test, m_xBj;
    Parton        * p_spect;
  public:
    inline Parton();
    inline Parton(ATOOLS::Flavour &,ATOOLS::Vec4D &,pst::code=pst::none);
    inline Parton(ATOOLS::Particle * const,pst::code=pst::none);

    inline ATOOLS::Flavour   const GetFlavour()  const;
    inline int               const GetCPartner() const;
    inline ATOOLS::Vec4D     const Momentum()    const;
    inline ATOOLS::Vec4D           Momentum();
    inline pst::code               GetType()     const;
    inline double            const KtStart()     const;
    inline double            const KtVeto()      const;
    inline double            const KtTest()      const;
    inline double            const ZTest()       const;
    inline double            const YTest()       const;
    inline double            const Xbj()         const;
    inline Parton *          const GetSpect()    const;
    inline void SetMomentum(const ATOOLS::Vec4D &);
    inline void SetStart(const double);
    inline void SetVeto(const double);
    inline void SetKtTest(const double);
    inline void SetZTest(const double);
    inline void SetYTest(const double);
    inline void SetTest(const double,const double,const double);
    inline void SetXbj(const double);
    inline void SetSpect(Parton *);
    friend std::ostream& operator<<(std::ostream &,const Parton &);
  };

  Parton::Parton() : 
    m_flav(ATOOLS::Flavour(ATOOLS::kf::none)),
    m_mom(ATOOLS::Vec4D(0.,0.,0.,0.)), 
    m_pst(pst::none), m_cpartner(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_z_test(1.), m_y_test(1.), m_xBj(1.), p_spect(NULL)
  { }

  Parton::Parton(ATOOLS::Flavour & flav,ATOOLS::Vec4D & mom,pst::code pst) : 
    m_flav(flav), m_mom(mom), m_pst(pst), m_cpartner(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_z_test(1.), m_y_test(1.), m_xBj(1.), p_spect(NULL)
  {
    if (m_flav.IsQuark()) {
      if (m_flav.IsAnti()) m_cpartner = -1;
                      else m_cpartner =  1;
    }
  }

  Parton::Parton(ATOOLS::Particle * const part,pst::code pst) : 
    m_flav(part->Flav()), m_mom(part->Momentum()), 
    m_pst(pst), m_cpartner(0),
    m_kt_start(0.), m_kt_veto(0.), m_kt_test(0.),
    m_z_test(1.), m_y_test(1.), m_xBj(1.), p_spect(NULL)
  {
    if (m_flav.IsQuark()) {
      if (m_flav.IsAnti()) m_cpartner = -1;
                      else m_cpartner =  1;
    }
  }


  ATOOLS::Flavour   const Parton::GetFlavour()  const { return m_flav; }
  ATOOLS::Vec4D     const Parton::Momentum()    const { return m_mom; }
  ATOOLS::Vec4D           Parton::Momentum()          { return m_mom; }
  pst::code               Parton::GetType()     const { return m_pst; }
  int               const Parton::GetCPartner() const { return m_cpartner; } 
  double            const Parton::KtStart()     const { return m_kt_start; }
  double            const Parton::KtVeto()      const { return m_kt_veto; }
  double            const Parton::KtTest()      const { return m_kt_test; }
  double            const Parton::ZTest()       const { return m_z_test; }
  double            const Parton::YTest()       const { return m_y_test; }
  double            const Parton::Xbj()         const { return m_xBj; }
  Parton *          const Parton::GetSpect()    const { return p_spect; }
  void Parton::SetMomentum(const ATOOLS::Vec4D & mom) { m_mom      = mom; }
  void Parton::SetStart(const double kt)              { m_kt_start = kt; }
  void Parton::SetVeto(const double kt)               { m_kt_veto  = kt; }
  void Parton::SetKtTest(const double kt)             { m_kt_test  = kt; }
  void Parton::SetZTest(const double z)               { m_z_test   = z; }
  void Parton::SetYTest(const double y)               { m_y_test   = y; }  
  void Parton::SetTest(const double kt,
		       const double z, 
		       const double y)                { m_kt_test  = kt; 
                                                        m_z_test   = z; 
                                                        m_y_test   = y; }
  void Parton::SetXbj(const double x)                 { m_xBj      = x; }
  void Parton::SetSpect(Parton * part)                { p_spect    = part; }  
  typedef std::list<Parton *>   Parton_List;
  typedef Parton_List::iterator PLiter;
}

#endif
