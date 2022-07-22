#ifndef SHRIMPS_Beam_Remnants_Beam_Remnant_Handler_H
#define SHRIMPS_Beam_Remnants_Beam_Remnant_Handler_H

#include "SHRiMPS/Beam_Remnants/Hadron_Dissociation.H"
#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "SHRiMPS/Eikonals/Omega_ik.H"

namespace BEAM{ class Beam_Spectra_Handler; }
namespace PDF { class ISR_Handler; }

namespace SHRIMPS {

  class Remnant_Handler {
  private:
    Colour_Generator * p_colourgenerator; 
    Form_Factor      * p_ff[2];
    std::vector<Hadron_Dissociation *> m_hadrons;
    
    void InitialiseCollision(ATOOLS::Blob_List * blobs);
  public:
    Remnant_Handler(PDF::ISR_Handler *const isr);
    ~Remnant_Handler();

    void Reset();
    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List * blobs,
					     const double & B=0.);
    void SetBeamBlob(ATOOLS::Blob *const beamblob,const int & beam);

    inline void SetFormFactors(Form_Factor * ff1,Form_Factor * ff2) {
      p_ff[0] = ff1; p_ff[1] = ff2;
    }
    inline void SetColourGenerator(Colour_Generator * colourgenerator) {
      p_colourgenerator = colourgenerator;
    }
    inline Hadron_Dissociation * GetHadronDissociation(const size_t & beam) {
      return m_hadrons[beam];
    }
    inline Continued_PDF * GetPDF(const size_t & beam) {
      return m_hadrons[beam]->GetPDF();
    }
  };
}
#endif
