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
    bool        m_on, m_done;
    size_t      m_nsamp;
    double      m_nevents, m_widthcut;
    std::string m_outdir;

    HADRON_RESCATTERING::Deuteron_Coalescence * p_deuterons;

    std::map<std::string,ATOOLS::Histogram *> m_histos;
    // per-run sums of per-event yields and of their squares (for errors)
    std::map<std::string,double>              m_counts, m_counts2, m_evcounts;
    // pair counts in (k, |p1+p2|), cf. Fig. 6 of the paper
    std::map<std::string,std::vector<double> > m_pairmaps;
    // ALICE multiplicity classes (arXiv:1902.09290): per value of the V0M
    // charged multiplicity, sums over INEL>0 events of the mid-rapidity
    // invariant yields (and of their squares) of d, dbar, p, pbar
    std::vector<double> m_mult_events, m_mult_nch;
    std::map<std::string,std::vector<double> > m_mult_yield, m_mult_yield2, m_mult_dndy;
    std::map<std::string,std::vector<double> > m_ev_yield;
    std::map<std::string,double>               m_ev_dndy;

    bool FromWeakDecay(ATOOLS::Particle * part) const;
    void FillDeuteron(const std::string & tag,const bool & anti,
		      const ATOOLS::Vec4D & mom,const double & weight);
    void InitHistograms();
  public :
    Deuteron_Formation();
    ~Deuteron_Formation();

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}

#endif
