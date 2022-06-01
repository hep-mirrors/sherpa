#ifndef SHERPA_Tools_HepMC2_Interface_H
#define SHERPA_Tools_HepMC2_Interface_H
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__HEPMC2

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "ATOOLS/Phys/Variations.H"

#include <map>
#include <set>

namespace HepMC {
  class GenEvent;
  class GenVertex;
  class GenParticle;
}

namespace ATOOLS {
  class Blob_List;
  class Particle;
  class ME_Weight_Info;
  class NLO_subevtlist;
  class PDF_Info;
}

namespace SHERPA {

  class EventInfo {
  private:
    ATOOLS::Blob * p_sp;
    bool m_usenamedweights, m_extendedweights;
    std::vector<ATOOLS::Variations_Source> m_variationsources;
    std::vector<double> m_orders;
    double m_wgt, m_mewgt, m_wgtnorm, m_ntrials, m_pswgt, m_pwgt;
    double m_mur2, m_muf12, m_muf22;
    double m_alphas, m_alpha;
    bool m_userhook;
    double m_userweight;
    ATOOLS::nlo_type::code m_type;
    ATOOLS::ME_Weight_Info * p_wgtinfo;
    ATOOLS::PDF_Info * p_pdfinfo;
    ATOOLS::NLO_subevtlist * p_subevtlist;
    ATOOLS::Weights_Map m_wgtmap;

    void ReadIn(ATOOLS::Blob_Data_Base * &db, std::string name, bool abort);

  public:
    EventInfo(ATOOLS::Blob * sp, const double &wgt,
              bool namedweights,
              bool extendedweights,
              bool includemeonlyweights);
    EventInfo(const EventInfo &evtinfo);
    ~EventInfo() {}
    bool WriteTo(HepMC::GenEvent &evt, const int& idx=-1);
    inline const ATOOLS::NLO_subevtlist * SubEvtList() { return p_subevtlist; }
    inline void SetWeight(const double &wgt)          { m_wgt=wgt; }
    inline void SetMEWeight(const double &mewgt)      { m_mewgt=mewgt; }
    inline void SetPartonicWeight(const double& pwgt) { m_pwgt=pwgt; }
    inline void SetMuR2(const double &mur2)           { m_mur2=mur2; }
    inline void SetMuF12(const double &muf12)         { m_muf12=muf12; }
    inline void SetMuF22(const double &muf22)         { m_muf22=muf22; }
    void SetAlphaS();
    void SetAlpha();
  };

  class HepMC2_Interface {
  private :
    bool m_usenamedweights, m_extendedweights, m_includemeonlyweights;
    bool m_hepmctree;

    std::map<ATOOLS::Blob*,HepMC::GenVertex *> m_blob2genvertex;
    std::map<ATOOLS::Particle*,HepMC::GenParticle *> m_particle2genparticle;
    std::set<ATOOLS::btp::code> m_ignoreblobs;

    HepMC::GenEvent *p_event;
    std::vector<HepMC::GenEvent*> m_subeventlist;

    bool Sherpa2ShortHepMC(const ATOOLS::Vec4D&,
                           const ATOOLS::Flavour&,
                           bool incoming,
                           HepMC::GenParticle*&);

  public :

    // constructor
    HepMC2_Interface();

    // destructor
    ~HepMC2_Interface();

    // member functions
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs, HepMC::GenEvent& event);
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs);
    bool Sherpa2HepMC(ATOOLS::Blob *,     HepMC::GenVertex *& ); 
    bool Sherpa2HepMC(ATOOLS::Particle *, HepMC::GenParticle *&);

    bool Sherpa2ShortHepMC(ATOOLS::Blob_List* const, HepMC::GenEvent&);
    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs);
    bool SubEvtList2ShortHepMC(EventInfo &evtinfo);

    void AddCrossSection(HepMC::GenEvent&, const Uncertain<double>& xs);

    // inline functions
    inline HepMC::GenEvent* GenEvent() { return p_event; }
    // For RivetInterface to enable ME/shower only analysis of fully hadronised event
    inline void Ignore(ATOOLS::btp::code type) { m_ignoreblobs.insert(type); }

    inline HepMC::GenEvent* GenSubEvent(size_t i)
    { if (i<m_subeventlist.size()) return m_subeventlist[i]; return NULL; }
    inline std::vector<HepMC::GenEvent*> GenSubEventList()
    { return m_subeventlist; }

    inline void SetHepMCNamedWeights(bool x)   { m_usenamedweights=x; }
    inline bool HepMCNamedWeights()            { return m_usenamedweights; }
    inline void SetHepMCExtendedWeights(bool x){ m_extendedweights=x; }
    inline bool HepMCExtendedWeights()         { return m_extendedweights; }
    inline void SetHepMCTreeLike(bool x)       { m_hepmctree=x; }
    inline bool HepMCTreeLike()                { return m_hepmctree; }
    inline void SetHepMCIncludeMEOnlyVariations(bool x)
    { m_includemeonlyweights=x; }

    // For RivetInterface to check if a weight container key could be a variation weight name
    bool StartsLikeVariationName(const std::string&);

    void DeleteGenSubEventList();

  };// end of class HepMC2_Interface

}// end of namespace SHERPA

#endif
#endif
