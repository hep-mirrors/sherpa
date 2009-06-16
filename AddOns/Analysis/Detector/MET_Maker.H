#ifndef Analysis_Detector_MET_Maker_H
#define Analysis_Detector_MET_Maker_H

#include "AddOns/Analysis/Detector/Object_Definition_Base.H"
#include "ATOOLS/Phys/Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class MET_Maker : public Object_Definition_Base {
  private:
    int    m_corrmode;
    double m_cutHad, m_cutEM;
    double m_deviation;
    double m_Rinfty_ECal, m_ampl_ECal, m_slope_ECal;
    double m_Rinfty_HCal, m_ampl_HCal, m_slope_HCal;

    double ReconstructedE_ECal(const double);
    double ReconstructedE_HCal(const double);
  public:
    MET_Maker(Primitive_Analysis *,const std::string="ET_UP");
    ~MET_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector);
    void SetCuts(const double cutHad,const double cutEM);
    void SetECorrection(const int,
			const double,const double=0.,const double=0.,
			const double=0.,const double=0.,const double=0.);
  };
}

#endif
