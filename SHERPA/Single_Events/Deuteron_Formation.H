#ifndef SHERPA_Single_Events_Deuteron_Formation_H
#define SHERPA_Single_Events_Deuteron_Formation_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Math/Vector.H"
#include <map>
#include <string>
#include <vector>

namespace ATOOLS {
  class Blob_List;
  class Histogram;
  class Particle;
}
namespace HADRON_RESCATTERING {
  class Deuteron_Coalescence;
}

namespace SHERPA {
  // (Anti)deuteron formation following arXiv:1504.07242: after all hadron
  // decays, every unique pair of final state (anti)nucleons not stemming
  // from weak decays forms a deuteron in channel i with P_i = sigma_i/sigma_0.
  // Besides the formed deuterons, "expected" spectra are filled with weight
  // P_i from every pair, sampled NSAMP times, for low-noise spectra.
  class Deuteron_Formation : public Event_Phase_Handler {
  private :
    typedef std::map<std::string,ATOOLS::Histogram *>  Histogram_Map;
    typedef std::map<std::string,double>               Count_Map;
    typedef std::map<std::string,std::vector<double> > Vector_Map;

    bool        m_on, m_done;
    size_t      m_nsamp;
    double      m_nevents, m_widthcut;
    std::string m_outdir;

    HADRON_RESCATTERING::Deuteron_Coalescence * p_deuterons;

    // all histograms, written to <OUTPUT>/<name>.dat in Finish
    Histogram_Map m_histos;
    // numbers counted in the current event (deuterons per channel and
    // momentum window, nucleons, pairs, ...), keyed by name; reset for every
    // event and added to the run sums at its end
    Count_Map     m_event_counts;
    // run sums of these per-event numbers and of their squares, giving the
    // per-event means and their errors in summary.dat
    Count_Map     m_count_sums, m_count_square_sums;
    // pair counts in (k, |p1+p2|), cf. Fig. 6 of the paper
    Vector_Map    m_pairmaps;

    // ALICE multiplicity classes (arXiv:1902.09290).  For the current event,
    // per species (d, dbar, p, pbar): number at |y|<0.5, i.e. dN/dy, and the
    // invariant yield 1/(2 pi pT) dN/dpT dy at |y|<0.5 in pT bins
    Count_Map     m_event_midy_count;
    Vector_Map    m_event_midy_spectrum;
    // run sums over INEL>0 events, in rows of the V0M charged multiplicity nV0:
    // number of events, charged particles at |eta|<0.5, dN/dy per species,
    // and the spectra (and their squares, for errors) per species
    std::vector<double> m_mult_events, m_mult_nch_sum;
    Vector_Map    m_mult_midy_count_sum, m_mult_spectrum_sum, m_mult_spectrum_square_sum;

    bool FromWeakDecay(ATOOLS::Particle * part) const;
    std::string ProtonOrigin(ATOOLS::Particle * part) const;
    void InitHistograms();
    // steps of Treat
    bool EventComplete(ATOOLS::Blob_List * bloblist) const;
    void ResetEvent();
    void CountEventActivity(ATOOLS::Blob_List * bloblist,size_t & nV0,
                            size_t & nch_mid,bool & inel0) const;
    void CollectNucleons(ATOOLS::Blob_List * bloblist,
                         std::vector<ATOOLS::Particle *> & nucleons);
    void FillProton(ATOOLS::Particle * part);
    void FillMesons(ATOOLS::Blob_List * bloblist);
    void FillPartners(const std::vector<ATOOLS::Particle *> & nucleons,
                      ATOOLS::Blob_List * bloblist);
    void FillExpected(const std::vector<ATOOLS::Particle *> & nucleons);
    bool FormDeuterons(const std::vector<ATOOLS::Particle *> & nucleons,
                       ATOOLS::Blob_List * bloblist);
    void FillDeuteron(const std::string & tag,const bool & anti,
                      const ATOOLS::Vec4D & mom,const double & weight);
    void AddEventToRunSums(const size_t & nV0,const size_t & nch_mid,
                           const bool & inel0);
    // steps of Finish
    void WriteHistograms() const;
    void WritePairMaps() const;
    void WriteMultiplicityClasses() const;
    void WriteSummary() const;
  public :
    Deuteron_Formation();
    ~Deuteron_Formation();

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}

#endif
