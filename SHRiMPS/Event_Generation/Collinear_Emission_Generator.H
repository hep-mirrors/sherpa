#ifndef SHRIMPS_Main_Collinear_Emission_Generator_H
#define SHRIMPS_Main_Collinear_Emission_Generator_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Strong_Coupling.H"
#include <list>

namespace SHRIMPS {
  class Collinear_Emission_Generator {
  private:
    ATOOLS::Particle_Vector  m_inparts, m_outparts, m_beamparts, m_gluons;
    MODEL::Strong_Coupling * p_alphaS;
    double                   m_kt2min, m_kt2max, m_sladder;
    std::map<int,double>     m_kt2starts;

    void AddEmissions();
    void GetKt2(bool isgluon, double & kt2, double & z);
    double GetZ(bool isgluon);
    ATOOLS::Particle * FindSpectator(ATOOLS::Particle * part);
    bool FixKinematics(double kt2, double z, ATOOLS::Vec4D & split, ATOOLS::Vec4D & gluon, ATOOLS::Vec4D & spec);
    void FillBlob(ATOOLS::Blob * blob);
    void CleanUp();
  public:
    Collinear_Emission_Generator();
    ~Collinear_Emission_Generator();

    int GenerateEmissions(ATOOLS::Blob_List * blobs);
  };
}
#endif
