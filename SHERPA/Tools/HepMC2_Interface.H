#ifndef SHERPA_Tools_HepMC2_Interface_H
#define SHERPA_Tools_HepMC2_Interface_H
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__HEPMC2

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHASIC++/Process/Subprocess_Info.H"


#include <map>
#include <set>

namespace HepMC {
  class GenEvent;
  class GenVertex;
  class GenParticle;
  class WeightContainer;
}

namespace ATOOLS {
  class Blob_List;
  class Particle;
  class ME_Weight_Info;
  class NLO_subevtlist;
  class PDF_Info;
}

namespace SHERPA {

  class Scale_Variation;
  typedef std::map<std::string, Scale_Variation *> NamedScaleVariationMap;

  class EventInfo {
  private:
    size_t m_usenamedweights;
    ATOOLS::Blob * p_sp;
    int    m_oqcd, m_oew;
    double m_wgt, m_mewgt, m_wgtnorm, m_ntrials, m_mur2, m_muf12, m_muf22;
    double m_alphas, m_alpha;
    PHASIC::nlo_type::code m_type;
    ATOOLS::ME_Weight_Info * p_wgtinfo;
    ATOOLS::PDF_Info * p_pdfinfo;
    ATOOLS::NLO_subevtlist * p_subevtlist;
    NamedScaleVariationMap * p_nsvmap;

    bool ReadIn(ATOOLS::Blob_Data_Base * &db, std::string name, bool abort);

  public:
    EventInfo(ATOOLS::Blob * sp, const double &wgt, bool namedweights);
    EventInfo(const EventInfo &evtinfo);
    ~EventInfo() {}
    bool WriteTo(HepMC::GenEvent &evt, const int& idx=-1);
    inline const ATOOLS::NLO_subevtlist * SubEvtList() { return p_subevtlist; }
    inline void SetWeight(const double &wgt)     { m_wgt=wgt; }
    inline void SetMEWeight(const double &mewgt) { m_mewgt=mewgt; }
    inline void SetMuR2(const double &mur2)      { m_mur2=mur2; }
    inline void SetMuF12(const double &muf12)    { m_muf12=muf12; }
    inline void SetMuF22(const double &muf22)    { m_muf22=muf22; }
  };

  class HepMC2_Interface {
  private :
    bool m_usenamedweights;

    std::map<ATOOLS::Blob*,HepMC::GenVertex *> m_blob2genvertex;
    std::map<ATOOLS::Particle*,HepMC::GenParticle *> m_particle2genparticle;
    std::set<ATOOLS::btp::code> m_ignoreblobs;

    HepMC::GenEvent *p_event;
    std::vector<HepMC::GenEvent*> m_subeventlist;

  public :

    // constructor
    HepMC2_Interface();

    // destructor
    ~HepMC2_Interface();

    // member functions
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs, HepMC::GenEvent& event,
                      double weight=1.0);
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs, double weight=1.0);
    bool Sherpa2HepMC(ATOOLS::Blob *,HepMC::GenVertex *&,
		      const std::map<size_t,size_t> &decids);
    bool Sherpa2HepMC(ATOOLS::Particle *,HepMC::GenParticle *&,
		      const std::map<size_t,size_t> &decids);

    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                           HepMC::GenEvent& event,
                           double weight=1.0);
    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs, double weight=1.0);
    bool SubEvtList2ShortHepMC(EventInfo &evtinfo);

    // inline functions
    inline HepMC::GenEvent* GenEvent() { return p_event; }
    inline void Ignore(ATOOLS::btp::code type) { m_ignoreblobs.insert(type); }

    inline HepMC::GenEvent* GenSubEvent(size_t i)
    { if (i<m_subeventlist.size()) return m_subeventlist[i]; return NULL; }
    inline std::vector<HepMC::GenEvent*> GenSubEventList()
    { return m_subeventlist; }

    inline void SetHepMCNamedWeights(bool x) { m_usenamedweights=x; }
    inline bool HepMCNamedWeights()          { return m_usenamedweights; }

    void DeleteGenSubEventList();

  };// end of class HepMC2_Interface

}// end of namespace SHERPA

#endif
#endif
