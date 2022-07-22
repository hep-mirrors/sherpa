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
    int             m_beam, m_dir;
    double          m_xmin, m_QT2max, m_expo;

    ATOOLS::Blob                             * p_blob, * p_softblob;
    std::list<ATOOLS::Flavour>                 m_beamflavs;
    std::set<int>                              m_beamcols[2];
    ATOOLS::Flavour                            m_quark, m_diquark;
    std::map<ATOOLS::Particle *,ATOOLS::Vec4D> m_qtmap;

    bool CheckResiduals();
    void AddBeamBlob(ATOOLS::Blob_List * blobs,const double & B);
    void HarvestCollision(ATOOLS::Blob_List * blobs);
    bool CompensateFlavours();
    bool AddFlavourCompensator(const ATOOLS::Flavour & flav);
    bool AddColourCompensator(const ATOOLS::Flavour & flav);
    bool SelectCompensatorMomentum(ATOOLS::Particle * part);
    bool SelectCompensatorTripletColours(ATOOLS::Particle * part);
    bool CompensateColours();
    void CleanColours();
    void AddSpectatorPartons();
    void FixConstituentFlavours();
    void CalculateParallelMomenta(ATOOLS::Vec4D & qmom,ATOOLS::Vec4D & dimom);
    void SelectTrialTransverseMomenta();
  public:
    Hadron_Dissociation() {};
    Hadron_Dissociation(const int & beam,BEAM::Beam_Base * beambase,Continued_PDF * pdf);
    ~Hadron_Dissociation();

    void Reset();
    bool FillBeamBlob(ATOOLS::Blob_List * blobs, const double & B=0.);
    inline void SetSoftBlob(ATOOLS::Blob * blob=NULL) { p_softblob = blob; }
    inline void SetFormFactor(Form_Factor * ff)       { p_ff = ff; }
    inline Continued_PDF * GetPDF() const { return p_pdf; }
    inline ATOOLS::Blob  * GetSoftBlob()  { return p_softblob; }
   };
}
#endif
