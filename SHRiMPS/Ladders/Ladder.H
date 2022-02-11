#ifndef SHRIMPS_Ladders_Ladder_H
#define SHRIMPS_Ladders_Ladder_H

#include "SHRiMPS/Ladders/Ladder_Particle.H"
#include "SHRiMPS/Ladders/T_Prop.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include <map>
#include <list>

namespace SHRIMPS {
  class Ladder {
  private:
    LadderMap       m_emissions;
    TPropList       m_tprops;
    ATOOLS::Vec4D   m_position;
    Ladder_Particle m_inpart[2];
    bool            m_isrescatter;
  public:
    Ladder(const ATOOLS::Vec4D & position=ATOOLS::Vec4D(0.,0.,0.,0.),
	   const bool & rescatter=false);
    ~Ladder();

    Ladder_Particle *
    AddRapidity(const double y,
		const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_gluon),
		const ATOOLS::Vec4D & mom=ATOOLS::Vec4D(0.,0.,0.,0.));
    void DeleteRapidity(LadderMap::iterator yiter);
    void AddPropagator(T_Prop prop);
    void DeletePropagator(TPropList::iterator piter);
    void UpdatePropagatorKinematics();
    void Reset(const bool & all=true);
    void ResetFS();
    
    void OutputRapidities();    
    bool ExtractHardest(TPropList::iterator & winner,const double & qt2min=0.);
    void HardestIncomingMomenta(const TPropList::iterator & winner,
				ATOOLS::Vec4D & q0,ATOOLS::Vec4D & q1);

    inline LadderMap * GetEmissions()                   { return &m_emissions; }
    inline TPropList * GetProps()                       { return &m_tprops; }
    inline const ATOOLS::Vec4D & Position()       const { return m_position; }
    inline Ladder_Particle * InPart(const size_t & pos) { return pos<2?&m_inpart[pos]:NULL; }
    inline const bool IsRescatter()               const { return m_isrescatter; }
    inline Ladder_Particle * FSPart(const double y) {
      return (m_emissions.find(y)==m_emissions.end() ?
	      NULL : &m_emissions.find(y)->second);
    }
    inline ATOOLS::Vec4D FSMomentum() {
      ATOOLS::Vec4D sum(0.,0.,0.,0.);
      for (LadderMap::iterator lit=m_emissions.begin();lit!=m_emissions.end();lit++)
	sum += lit->second.Momentum();
      return sum;
    }
    bool CheckFourMomentum();
    friend std::ostream & operator<<(std::ostream & s, const Ladder &);
  };

  std::ostream & operator<<(std::ostream & s, Ladder &);
}
#endif
