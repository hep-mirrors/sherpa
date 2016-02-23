#ifndef SHRIMPS_Beam_Remnants_Beam_Remnant_Handler_H
#define SHRIMPS_Beam_Remnants_Beam_Remnant_Handler_H

#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <vector>

namespace BEAM   { class Beam_Spectra_Handler; }

namespace SHRIMPS {

  class Beam_Remnant_Handler {
  private:
    ATOOLS::Blob                 * p_blob;
    std::vector<Continued_PDF>   * p_pdfs;
    std::vector<BEAM::Beam_Base *> m_beams;
    std::vector<ATOOLS::Vec4D>     m_beamvecs, m_outmoms;

    void AddBeamBlobs(ATOOLS::Blob_List * blobs);
    void UnsetNeedsBeamsStatus(ATOOLS::Blob_List * blobs);
public:
    Beam_Remnant_Handler() {};
    Beam_Remnant_Handler(BEAM::Beam_Spectra_Handler * beamspectra,
			 std::vector<Continued_PDF> & pdfs);
    ~Beam_Remnant_Handler() {}

    void Reset();
    bool InitialiseCollision();
    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List * blobs);

    inline void SetBeamBlob(ATOOLS::Blob *const beamblob,const int & beam) {} 
    inline Continued_PDF * GetPDF(const int & beam) {
      if (beam!=0 && beam!=1) return NULL;
      return &(*p_pdfs)[beam];
    }
  };
}
#endif
