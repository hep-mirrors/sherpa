#ifndef Calorimeter_Cone_H
#define Calorimeter_Cone_H

#include "Primitive_Calorimeter.H"
#include "Kt_Algorithm.H"
#include "Particle_List.H"
#include "Vector.H"

namespace ANALYSIS {
  struct Jet_Data {
    int           i,j;
    ATOOLS::Vec4D mom;
    double        et;
    Jet_Data(int _i, int _j, const ATOOLS::Vec4D & _mom, double _et) :
      i(_i),j(_j),mom(_mom),et(_et) {}
  };

  class Calorimeter_Cone : public ATOOLS::Jet_Algorithm_Base {
    double m_dR, m_dR2, m_Etcut, m_Etstop, m_maxetajet, m_minetajet;
    int    m_neta, m_nphi, m_dneta, m_dnphi;
    double m_mineta, m_maxeta, m_delta_eta, m_delta_phi;
    int    m_etamode;

    int    ** p_jetno;

    Primitive_Calorimeter * p_calorimeter;
    
    std::vector<Jet_Data> m_jets;

    void CalcJets();
  public:
    Calorimeter_Cone(const double dR, const double Et_cut, 
		     Primitive_Calorimeter * const);
    ~Calorimeter_Cone();

    bool ConstructJets(ATOOLS::Particle_List *,std::vector<double> *);
    bool ConstructJets(const ATOOLS::Particle_List * pl,ATOOLS::Particle_List * jets,
		       std::vector<double> * kt2,double sep);
    void SetEtaRangeForJets(const double mineta,const double maxeta,const int mode);
  };

  inline bool Calorimeter_Cone::ConstructJets(const ATOOLS::Particle_List * pl,
					      ATOOLS::Particle_List * jets,
					      std::vector<double> * kt2,double sep) 
  {
    return ConstructJets(jets,kt2);
  }

  inline void Calorimeter_Cone::SetEtaRangeForJets(const double mineta,
						   const double maxeta,const int mode) 
  {
    if (mineta<maxeta) { m_minetajet = mineta; m_maxetajet = maxeta; }
                  else { m_minetajet = maxeta; m_maxetajet = mineta; }
    m_etamode = mode;
  }

  // a simplied version of the UA1  jet algorithm with jetradius and 
  //   minimum scalar transverse energy
  //        (RJET=1., EJCUT=5. FOR UA1)
  //  getjet.f obtained from MLM Feb 2004
  //    "Calorimeter simulation obtained from Frank Paige 23 March 1988"
  //    (RESEM=.15, RESHAD=.35 FOR URANIUM CALORIMETER)
  //    (RESEM=.15, RESHAD=.70 FOR IRON CALORIMETER)
  //    (RESEM=.11, RESHAD=.70 FOR CDF)

  //
}

#endif
