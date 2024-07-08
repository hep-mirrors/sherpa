#ifndef AMISIC_Perturbative_MI_Process_H
#define AMISIC_Perturbative_MI_Process_H
/*!
  \file MI_Process.H
  \brief Declares the class MI_Process
*/

#include "EXTRA_XS/Main/ME2_Base.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"

namespace AMISIC { class MI_Integrator; }

namespace AMISIC {
  class XS_Base {
  protected:
    std::string m_name;
    double      m_Ms, m_Mt, m_Mu, m_lastxs, m_shat, m_that, m_uhat;
    std::vector<std::vector<int> > m_colours;
    std::vector<double>            m_masses, m_masses2;
  public:
    XS_Base();
    explicit XS_Base(const std::vector<double> & masses);
    virtual ~XS_Base();

    virtual void Calc(const double & s,const double & t,const double & u)=0;
    virtual bool SetColours(const ATOOLS::Flavour_Vector & flavs)=0;
    virtual double operator()() { return m_lastxs; }
    virtual const int & Colour(const size_t & i,const size_t & j) const {
      return m_colours[i][j];
    }
    virtual const std::string & Name() const { return m_name; }
  };

#define DECLARE_XSBASE_GETTER(NAME,TAG)				     \
  DECLARE_GETTER(NAME,TAG,AMISIC::XS_Base,ATOOLS::Flavour_Vector);   \
  void ATOOLS::Getter<AMISIC::XS_Base,ATOOLS::Flavour_Vector,NAME>:: \
  PrintInfo(std::ostream &str,const size_t width) const		     \
  {								     \
    str<<#TAG;							     \
  }

  class MI_Process {
  private:
    std::string                  m_name;
    ATOOLS::Momenta_Stretcher    m_stretcher;
    std::vector<ATOOLS::Flavour> m_flavs;
    std::vector<double>          m_masses, m_masses2;
    std::vector<ATOOLS::Vec4D>   m_momenta;
    XS_Base*                     p_me2;
    double                       m_Emin, m_sumInMasses, m_sumOutMasses;
    bool                         m_masslessIS;

    const bool AllowedRemnants(std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> remnants)
            const;
    void MasslessISKinematics(MI_Integrator* integrator);
  public:
    explicit MI_Process(const std::vector<ATOOLS::Flavour>& flavs);
    ~MI_Process();
    ///////////////////////////////////////////////////////////////////////////
    // Creating four-vectors for a legal phase space point, checking for its
    // viability re remnants, masses, PDF limits, etc..
///////////////////////////////////////////////////////////////////////////
    bool MakeKinematics(
            MI_Integrator*                                         integrator,
            std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2>  remnants);
    ///////////////////////////////////////////////////////////////////////////
    // Making sure that the c.m. energy of the scatter is larger than the
    // IS or FS sum of masses.
    ///////////////////////////////////////////////////////////////////////////
    inline const bool AllowedKinematics(const double & Ehat) const {
      return (m_sumInMasses<Ehat && m_sumOutMasses<Ehat);
    }

    ATOOLS::Particle * GetParticle(const size_t & i);
    inline void SetME2(XS_Base * me2) { p_me2 = me2; }
    inline const ATOOLS::Flavour& Flav(const size_t& i) const
    {
      return m_flavs[i];
    }
    inline const ATOOLS::Vec4D& Momentum(size_t i) const
    {
      return m_momenta[i];
    }
    inline bool SetColours() const { return p_me2->SetColours(m_flavs); }
    inline const std::string& Name() const { return m_name; }
    inline double             operator()() const { return (*p_me2)(); }
  };
}

#endif
