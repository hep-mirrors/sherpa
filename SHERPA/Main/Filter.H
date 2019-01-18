#ifndef SHERPA_Filter_H
#define SHERPA_Filter_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob_List.H"
#include <map>

namespace SHERPA {
  struct FilterCriterion {
    ATOOLS::Flavour m_flav;
    double m_etamin, m_etamax, m_pTmin, m_pTmax;
    size_t m_Nmin, m_Nmax;
  };


  class Filter {
  private:
    std::map<ATOOLS::Flavour,FilterCriterion *> m_filters;
    std::list<ATOOLS::Particle * >              m_particles;
    std::map<ATOOLS::Flavour, int>              m_accepted;
    bool m_on;

    void Add(const std::vector<std::string> & entry);
    void Reset();
    void HarvestActiveParticles(ATOOLS::Blob_List * blobs);
    void FilterAccepted();
    bool Check();
  public:
    Filter();
    ~Filter();
    bool Init(const std::string & path, const std::string & file);
    static void ShowSyntax(int mode);

    bool operator()(ATOOLS::Blob_List * blobs);
  };
}

#endif
