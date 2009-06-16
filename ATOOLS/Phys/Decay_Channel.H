#ifndef ATOOLS_Phys_Decay_Channel_H
#define ATOOLS_Phys_Decay_Channel_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <set>

namespace HELICITIES {
  class Spin_Amplitudes;
}

namespace PHASIC {
  class Multi_Channel;
}

namespace ATOOLS {
  class Decay_Channel {
  protected :
    double      m_width, m_deltawidth, m_minmass, m_max;
    Flavour     m_flin;
    FlavourSet  m_flouts;
    std::vector<HELICITIES::Spin_Amplitudes*> m_diagrams;
    PHASIC::Multi_Channel* p_channels;
  public :
    Decay_Channel() {}
    Decay_Channel(const Flavour &);
    Decay_Channel(const Decay_Channel &);
    virtual ~Decay_Channel();
    inline void AddDecayProduct(const Flavour & flout) {
      m_flouts.insert(flout); 
      m_minmass += flout.HadMass();
    }
    inline void SetWidth(const double & width) { m_width =width; }
    inline void SetDeltaWidth(const double & delta) { m_deltawidth=delta; }

    std::string Name() const;
    inline const double& Width() const { return m_width; }
    inline const double& DeltaWidth() const { return m_deltawidth; }
    inline const double& MinimalMass() const { return m_minmass; }
    double DiceMass(const double& min, const double& max) const;

    inline int NOut() const { return m_flouts.size(); }
    inline const Flavour& GetDecaying() const { return m_flin; }      
    inline const FlavourSet& GetDecayProducts() const { return m_flouts; }
    inline const Flavour& GetDecayProduct(size_t i) const {
      if (i<0 || i>m_flouts.size())
        THROW(fatal_error,ToString(i)+" out of bounds.");
      for (FlSetConstIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) {
        if (i==0) return (*flit);
        --i;
      }
      return *m_flouts.begin();
    }
    inline const double& Max() const { return m_max; }
    inline void SetMax(const double& max) { m_max=max; }
    inline std::vector<HELICITIES::Spin_Amplitudes*>& GetDiagrams() {
      return m_diagrams;
    }
    inline void AddDiagram(HELICITIES::Spin_Amplitudes* amp) {
      m_diagrams.push_back(amp);
    }
    inline PHASIC::Multi_Channel* Channels() const { return p_channels; }
    inline void SetChannels(PHASIC::Multi_Channel* chans) { p_channels=chans; }
    
    void Output() const;

    friend std::ostream &operator<<(std::ostream &os, const Decay_Channel &dt);
  };
}
#endif
