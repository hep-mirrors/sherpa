#ifndef SHRIMPS_Beam_Remnants_Beam_Remnant_Handler_H
#define SHRIMPS_Beam_Remnants_Beam_Remnant_Handler_H

#include "SHRiMPS/Beam_Remnants/Hadron_Dissociation.H"
#include "SHRiMPS/Eikonals/Omega_ik.H"

namespace BEAM   { class Beam_Spectra_Handler; }

namespace SHRIMPS {

  class Beam_Remnant_Handler {
  private:
    std::vector<Hadron_Dissociation *> m_hadrons;
    Omega_ik     * p_eikonal; 
    ATOOLS::Blob * p_softblob;

    std::map<ATOOLS::Particle *,ATOOLS::Vec4D> m_qtmap;

    void ConstructAndAddSoftBlob(ATOOLS::Blob_List * blobs);
    void AddBeamBlobs(ATOOLS::Blob_List * blobs);
    void AddTransverseMomentaToSpectators(ATOOLS::Blob_List * blobs);
  public:
    Beam_Remnant_Handler() {};
    Beam_Remnant_Handler(BEAM::Beam_Spectra_Handler * beamspectra,
			 std::vector<Continued_PDF> & pdfs);
    ~Beam_Remnant_Handler();

    void Reset();
    bool InitialiseCollision();
    void SetEikonal(Omega_ik * eikonal);
    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List * blobs);

    void SetBeamBlob(ATOOLS::Blob *const beamblob,const int & beam);
    inline Continued_PDF * GetPDF(const int & beam) {
      if (beam!=0 && beam!=1) return NULL;
      return m_hadrons[beam]->GetPDF();
    }
  };
}
#endif
