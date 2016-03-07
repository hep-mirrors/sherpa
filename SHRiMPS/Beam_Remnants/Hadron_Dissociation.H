#ifndef SHRIMPS_Beam_Remnants_Hadron_Dissociation_H
#define SHRIMPS_Beam_Remnants_Hadron_Dissociation_H

#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <set>

namespace BEAM   { class Beam_Base; }

namespace SHRIMPS {

  class Hadron_Dissociation {
  private:
    Continued_PDF   * p_pdf;
    ATOOLS::Flavour m_beamflav;
    ATOOLS::Vec4D   m_beamvec, m_outmom;
    int    m_dir;
    double m_xmin;

    ATOOLS::Blob  * p_blob;
    std::set<int>   m_cols[2];
    ATOOLS::Flavour m_quark, m_diquark;

    void MakeBeamBlob();
    void AddBeamParticle();
    void AddPartonsFromCollision(ATOOLS::Blob_List * blobs);
    void HarvestBlob(ATOOLS::Blob * blob);
    void IdentifyAndFillSoftBlob(ATOOLS::Blob_List * blobs);
    void AddSpectatorPartons(ATOOLS::Blob * softblob);
    void FixConstituentFlavours();
    void CalculateParallelMomenta(ATOOLS::Vec4D & qmom,ATOOLS::Vec4D & dimom);
  public:
    Hadron_Dissociation() {};
    Hadron_Dissociation(BEAM::Beam_Base * beambase,Continued_PDF * pdf);
    ~Hadron_Dissociation() {}

    void Reset();
    bool FillBeamBlob(ATOOLS::Blob_List * blobs);

    inline ATOOLS::Blob * GetBeamBlob() const { return p_blob; }
    inline Continued_PDF * GetPDF() const { return p_pdf; }
  };
}
#endif
