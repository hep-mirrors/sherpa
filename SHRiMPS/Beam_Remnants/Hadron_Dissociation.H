#ifndef SHRIMPS_Beam_Remnants_Hadron_Dissociation_H
#define SHRIMPS_Beam_Remnants_Hadron_Dissociation_H

#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <set>
#include <map>

namespace BEAM   { class Beam_Base; }

namespace SHRIMPS {

  class Hadron_Dissociation {
  private:
    Continued_PDF * p_pdf;
    Form_Factor   * p_ff;
    ATOOLS::Flavour m_beamflav;
    ATOOLS::Vec4D   m_beamvec, m_outmom;
    int    m_dir;
    double m_xmin, m_QT2max;

    ATOOLS::Blob  * p_blob;
    std::set<int>   m_cols[2];
    ATOOLS::Flavour m_quark, m_diquark;
    std::map<ATOOLS::Particle *,ATOOLS::Vec4D> * m_qtmap;
    
    void SpecifyBeamBlob();
    void AddPartonsFromCollision(ATOOLS::Blob_List * blobs);
    void HarvestBlob(ATOOLS::Blob * blob);
    void IdentifyAndFillSoftBlob(ATOOLS::Blob_List * blobs);
    void AddSpectatorPartons(ATOOLS::Blob * softblob);
    void FixConstituentFlavours();
    void CalculateParallelMomenta(ATOOLS::Vec4D & qmom,ATOOLS::Vec4D & dimom);
    void SelectTrialTransverseMomenta();
  public:
    Hadron_Dissociation() {};
    Hadron_Dissociation(BEAM::Beam_Base * beambase,Continued_PDF * pdf);
    ~Hadron_Dissociation() {}

    void Reset();
    bool FillBeamBlob(ATOOLS::Blob_List * blobs);
    void SetBeamBlob(ATOOLS::Blob * blob=NULL) { p_blob = blob; }
    
    inline void SetFormFactor(Form_Factor * ff) { p_ff = ff; }
    inline void SetQTMap(std::map<ATOOLS::Particle *,ATOOLS::Vec4D> * qtmap) {
      m_qtmap = qtmap;
    }
    inline ATOOLS::Blob  * GetBeamBlob() const { return p_blob; }
    inline Continued_PDF * GetPDF()      const { return p_pdf; }
  };
}
#endif
