#ifndef REMNANTS_Tools_Colour_Generator_H
#define REMNANTS_Tools_Colour_Generator_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"

namespace REMNANTS {
  class Remnant_Handler;
  class Remnant_Base;

  class Colour_Generator {
    class MassCrit {
    public:
      inline bool operator()(const std::pair<ATOOLS::Particle *, ATOOLS::Particle *> & s1,
              const std::pair<ATOOLS::Particle*, ATOOLS::Particle*>& s2) const
      {
        return ((s1.first->Momentum() + s1.second->Momentum()).Y() <
                (s2.first->Momentum() + s2.second->Momentum()).Y());
      }
    };

    typedef std::set<std::pair<ATOOLS::Particle*, ATOOLS::Particle*>, MassCrit>
            Singlets;

  private:
    std::list<int>     m_cols[2][2];
    std::set<int>      m_vetoed[2][2];
    std::array<std::shared_ptr<Remnant_Base>, 2> p_remnants;
    std::array<ATOOLS::Particle*, 2> p_inparts;

    bool   TChannelColourFlows();
    bool   SChannelColourFlows();
    bool   ConstrainedColourFlows(const size_t& tbeam);
    bool   ConstrainedGGFlows(const size_t& tbeam);
    bool   ConstrainedGQFlows(const size_t& tbeam);
    bool   ConstrainedQGFlows(const size_t& tbeam);
    bool   ConstrainedQQFlows(const size_t& tbeam);
    void   AssignColours(const size_t& beam, ATOOLS::Particle* trip,
                         ATOOLS::Particle* anti);
    size_t AvailableColours(const size_t& beam);
    int    DefineColourDonor(const size_t& beam);
    void   ReplaceBoth(const int& beam, const size_t& index);
    void Replace(const int& beam, const size_t& index, ATOOLS::Particle* part);
    void ReplaceInFS(const int& oldcol, const int& newcol, const size_t& index,
                     ATOOLS::Blob* blob);
    void ReplaceInIS(const int& oldcol, const int& newcol, const size_t& index,
                     ATOOLS::Blob* blob) const;

  public:
    Colour_Generator();
    ~Colour_Generator() = default;

    void Initialize(Remnant_Handler * remnants);
    bool ConnectColours(ATOOLS::Blob *const showerblob);
    void AddColour(const size_t& beam, const size_t& pos,
                   ATOOLS::Particle* part);
    int  NextColour(const size_t& beam, const size_t& pos);

    inline bool HasColours(const size_t & beam) const {
      return ((m_cols[beam][0].size()+m_cols[beam][1].size()) > 0);
    }

    inline const std::list<int>& Colours(const size_t& beam,
                                         const size_t& pos) const
    {
      return m_cols[beam][pos];
    }

    inline void Reset() {
      for (size_t beam=0;beam<2;beam++) {
	for (size_t pos=0;pos<2;pos++) {
	  m_cols[beam][pos].clear();
	  m_vetoed[beam][pos].clear();
	}
      }
    }
    inline void ResetFlags() {
      for (size_t beam=0;beam<2;beam++) {
	for (size_t pos=0;pos<2;pos++) m_vetoed[beam][pos].clear();
      }
    }
    void Output() const;
  };
}

#endif
