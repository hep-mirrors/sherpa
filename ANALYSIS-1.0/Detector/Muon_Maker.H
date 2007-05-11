#ifndef Muon_Maker_h
#define Muon_Maker_h

#include "Object_Definition_Base.H"
#include "Particle.H"

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
    virtual ~Muon_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *);

    void SetTracking(const double E);
    void SetHadIso(const double R,const double rel,const double cut);
    void SetEMIso(const double R,const double rel,const double cut);

    void SetECorrection(const double inv);
  };
}

#endif
