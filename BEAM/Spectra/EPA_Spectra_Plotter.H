#ifndef BEAM_Spectra_EPA_Spectra_Plotter_H
#define BEAM_Spectra_EPA_Spectra_Plotter_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace BEAM {
  class EPA;
  enum class EPA_ff_type;
  
  class EPA_Spectra_Plotter {
  private:
    EPA        * p_epa;
    ATOOLS::Flavour m_beam;
    double       m_R, m_thetamax, m_q2max;
    double       m_bmin, m_bmax,  m_bsize;
    double       m_xmin, m_xmax,  m_xsize;
    size_t       m_nxbins, m_nbbins, m_nxaxis, m_nbaxis;
    double       m_minx, m_maxx, m_minb, m_maxb;
    std::string  m_dir;
    ATOOLS::axis * p_xaxis, * p_baxis;
    std::vector<double> m_xvals;
    std::map<std::string, ATOOLS::Histogram * > m_histograms;

    void   Init();
    void   TestBesselFunctions();
    void   TestPoint();
    void   TestDipole();
    void   TestGauss();
    void   TestWoodSaxon();
    void   IonFFs();
    void   FillAnalytic(const enum EPA_ff_type & type,const double & Q2max=-1.);
    void   FillNumerical(const enum EPA_ff_type & type,const double & Q2max=-1.);
  public:
    EPA_Spectra_Plotter(EPA * epa,const std::string & dir=std::string("EPA_Spectra")) :
      p_epa(epa), m_dir(dir) {}
    ~EPA_Spectra_Plotter();
    void operator()(const size_t disc);
  };
}

#endif
