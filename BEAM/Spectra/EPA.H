#ifndef BEAM_Main_EPA_H
#define BEAM_Main_EPA_H

#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Spectra/EPA_FF.H"

namespace ATOOLS { class Settings; }

namespace BEAM {
  class EPA : public Beam_Base {
  private:
    EPA_ff_type   m_type;
    EPA_FF_Base * p_ff;
    double m_mass, m_charge, m_gamma;
    double m_aqed, m_pref;
    double m_x, m_q2;
    double m_q2min, m_q2max, m_theta_max, m_pt2min, m_pt2max;
    double m_Rmin, m_Rmax, m_Rcut, m_xmin, m_xmax, m_bmin, m_bmax;

    double m_Lambda2, m_mu, m_Q02, m_WSd;
    
    size_t m_plotting;
    size_t m_approx;
    size_t m_nxbins, m_nbbins;
    bool   m_analytic;

    void   Initialise();
    void   RegisterDefaults() const;
    void   InitFormFactor(ATOOLS::Settings &s);
    double ExtractParameter(ATOOLS::Settings &s,const std::string & tag);
    void   WriteDebugFiles(ATOOLS::Settings &s);

    void   Tests();

    // This methids can be deleted once we are happy everything works
    double phi(double x, double qq);
    double OldWeight(double x, double q2);
    void   selfTest(std::string filename);
  public:
    EPA(const ATOOLS::Flavour beam,const double energy,const double pol,
	const int dir);
    ~EPA() override = default;
    Beam_Base * Copy() { return new EPA(*this); }


    bool CalculateWeight(double x, double q2) {
      m_x = x; m_q2 = q2;
      m_weight = (*this)(m_x); return true;
    }
    inline const double operator()(const double & x) {
      return ( (x>=m_xmin && x<=m_xmax) ?
	       ATOOLS::Max(0., m_pref/x * p_ff->N(x)) : 0. ); 
    }
    inline double Xmax() const override { return m_xmax; }
    inline double Xmin() const override { return m_xmin; }

    void FixPosition();
    void SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i=0);
  };// end of class EPA

}// end of namespace BEAM

#endif
