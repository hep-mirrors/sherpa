#ifndef BEAM_Main_EPA_H
#define BEAM_Main_EPA_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Function_Base.H"

namespace BEAM {

  class EPA : public Beam_Base {
  private:

    class CosInt : public ATOOLS::Function_Base {
    public:
      double GetValue(double T) { return -cos(T)/T; }
      double operator()(double T) { return GetValue(T); }
      double GetCosInt(double X);
    };

    double m_mass, m_charge, m_gamma, m_minR, m_maxR;
    double m_q2Max, m_pt_min, m_aqed, m_xmin, m_theta_max;

    bool m_lo_epa;

    void RegisterDefaults();

    double phi(double x, double qq);

    void selfTest(std::string filename);

    int m_formfactor;

  public:

    EPA(const ATOOLS::Flavour beam,const double energy,const double pol,
	const int dir);

    ~EPA();

    Beam_Base *Copy();     
    bool CalculateWeight(double, double);
    void FixPosition();
    void SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i=0);
  };// end of class EPA

}// end of namespace BEAM

#endif
