#ifndef Analysis_Detector_Muon_Maker_H
#define Analysis_Detector_Muon_Maker_H

#include "AddOns/Analysis/Detector/Object_Definition_Base.H"
#include "ATOOLS/Phys/Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Muon_Maker : public Object_Definition_Base {
    int    m_dim;
    double m_Estart, m_R2track, m_cut_Ebyp;
    double m_R2hadiso, m_relhad, m_minhad;
    double m_R2EMiso, m_totEM, m_minEM;
    double m_inv;

    std::vector<Cell *>   * p_cluster;

    void GetTracks();
    void IsolateTracks();
    void CorrectEnergies();
  public:
    Muon_Maker(Primitive_Analysis *,const std::string="ET_UP");
    ~Muon_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector);

    void SetTracking(const double E);
    void SetHadIso(const double R,const double rel,const double cut);
    void SetEMIso(const double R,const double rel,const double cut);

    void SetECorrection(const double inv);
  };
}

#endif
