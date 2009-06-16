#ifndef BEAM_Main_EPA_H
#define BEAM_Main_EPA_H

#include "BEAM/Main/Beam_Base.H"

namespace ATOOLS { class Data_Reader; }

namespace BEAM {

  class EPA : public Beam_Base {
  private:

    double m_mass, m_charge, m_q2Max, m_pt_min, m_aqed;

    double phi(double x, double qq);

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
