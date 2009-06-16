#ifndef Analysis_Detector_Photon_Maker_H
#define Analysis_Detector_Photon_Maker_H

#include "AddOns/Analysis/Detector/Object_Definition_Base.H"
#include "ATOOLS/Phys/Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Photon_Maker : public Object_Definition_Base {
    bool   m_elements;
    int    m_dim;
    double m_Estart, m_R2track, m_trackcut;
    double m_R2hadiso, m_relhad, m_minhad;
    double m_R2EMiso, m_totEM, m_minEM;
    double m_inv;

    void BuildMatchedClusters();
    void IsolateClusters();
    void CorrectEnergies();
  public:
    Photon_Maker(Primitive_Analysis *,const std::string="ET_UP");
    ~Photon_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector);

    void SetClustering(const double E,const int dim);
    void SetTracking(const double R,const double cut);
    void SetHadIso(const double R,const double rel,const double cut);
    void SetEMIso(const double R,const double rel,const double cut);

    void SetECorrection(const double inv);
  };
}

#endif
