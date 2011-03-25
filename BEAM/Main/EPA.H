#ifndef BEAM_Main_EPA_H
#define BEAM_Main_EPA_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Function_Base.H"

namespace ATOOLS { class Data_Reader; }

namespace BEAM {

  class EPA : public Beam_Base {
  private:

    class CosInt : public ATOOLS::Function_Base {
    public:
      double GetValue(double T) { return -cos(T)/T; }
      double operator()(double T) { return GetValue(T); }
      double GetCosInt(double X);
    };

    double m_mass, m_charge, m_q2Max, m_pt_min, m_aqed;

    double phi(double x, double qq);

    void selfTest(std::string filename);

    int m_debug;
    int m_formfactor;

  public:

    EPA(const ATOOLS::Flavour beam,const double mass,
	const double charge,const double energy,const double pol,
	const int dir,ATOOLS::Data_Reader *const read);

    ~EPA();

    Beam_Base *Copy();     

    ATOOLS::Vec4D   OutMomentum();
    ATOOLS::Flavour Remnant();

    bool   CalculateWeight(double, double);
    double Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none));
    bool   On() { return 1; }

  };// end of class EPA

}// end of namespace BEAM

#endif
