#ifndef Analysis_Detector_Simple_Cone_H
#define Analysis_Detector_Simple_Cone_H

#include "AddOns/Analysis/Detector/ElMag_Calorimeter.H"
#include "AddOns/Analysis/Detector/Had_Calorimeter.H"
#include "AddOns/Analysis/Detector/Muon_Chambers.H"
#include "AddOns/Analysis/Detector/Reconstructed_Object.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Histogram.H"
#include <vector>

namespace ANALYSIS {
  // a simplied version of the UA1  jet algorithm with jetradius and 
  //   minimum scalar transverse energy
  //        (RJET=1., EJCUT=5. FOR UA1)
  //  getjet.f obtained from MLM Feb 2004
  //    "Simple simulation obtained from Frank Paige 23 March 1988"
  //    (RESEM=.15, RESHAD=.35 FOR URANIUM CALORIMETER)
  //    (RESEM=.15, RESHAD=.70 FOR IRON CALORIMETER)
  //    (RESEM=.11, RESHAD=.70 FOR CDF)

  //

  class Simple_Cone {
    double          m_Etcut, m_Etmin, m_dR, m_dR2;
    ATOOLS::Flavour m_flav;

    Had_Calorimeter   * p_HCal;
    ElMag_Calorimeter * p_ECal;
    Muon_Chambers     * p_MC;

    void CalcJets(ObjectList *);
  public:
    Simple_Cone(const double Etcut,const double Etmin,const double sep=1.);
    ~Simple_Cone();

    bool ConstructJets(ObjectList *);
    void FillShape(const int jetno,ATOOLS::Histogram * histo,
		   const double weight=1.,const double ncount=1);

    double Radius()       { return m_dR; }
    double Et_cut()       { return m_Etcut; }

    void   SetCalorimeters(Had_Calorimeter * hcal,ElMag_Calorimeter * ecal) {
      p_HCal = hcal; p_ECal = ecal;
    }
    void   SetMuonChambers(Muon_Chambers * mc) { p_MC = mc; }
  };
}

#endif
