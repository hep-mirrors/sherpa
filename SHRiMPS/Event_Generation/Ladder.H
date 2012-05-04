#ifndef SHRIMPS_Main_Ladder_H
#define SHRIMPS_Main_Ladder_H

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include <map>
#include <list>

namespace SHRIMPS {
  struct colour_type {
    enum code {
      none,
      singlet,
      triplet,
      octet
    };
  };

  std::ostream & operator<<(std::ostream & s, const colour_type::code & colour);

  struct Ladder_Particle {
    ATOOLS::Particle * p_part;
    ATOOLS::Flavour    m_flav;
    ATOOLS::Vec4D      m_mom, m_pos;
    ATOOLS::Flow       m_flow;
    bool               m_marked,m_beam,m_IS;

    Ladder_Particle() :
      p_part(NULL), m_flav(ATOOLS::Flavour(kf_none)), 
      m_mom(ATOOLS::Vec4D(0.,0.,0.,0.)), m_pos(ATOOLS::Vec4D(0.,0.,0.,0.)),
      m_flow(ATOOLS::Flow(NULL)),m_marked(false),m_beam(false),m_IS(false)
    {}

    Ladder_Particle(const Ladder_Particle & part) :
      p_part(NULL), m_flav(part.m_flav), m_mom(part.m_mom), m_pos(part.m_pos),
      m_flow(ATOOLS::Flow(NULL)),
      m_marked(part.m_marked),m_beam(part.m_beam),m_IS(part.m_IS)
    {}

    Ladder_Particle(const ATOOLS::Particle * part);

    Ladder_Particle(const ATOOLS::Flavour & flav,
		    const ATOOLS::Vec4D & mom=ATOOLS::Vec4D(0.,0.,0.,0.),
		    const ATOOLS::Vec4D & pos=ATOOLS::Vec4D(0.,0.,0.,0.)) :
      p_part(NULL), m_flav(flav), m_mom(mom), m_pos(pos),
      m_flow(ATOOLS::Flow(NULL)),m_marked(false),m_beam(false),m_IS(false)
    {}

    ~Ladder_Particle() { }

    inline void SetFlow(const unsigned int & pos,const int & code=-1) {
      if(code==-1) m_flow.SetCode(pos);
      else m_flow.SetCode(pos,code);
      if (p_part) p_part->SetFlow(pos,code);
    }
    inline int GetFlow(const unsigned int & pos,const bool & test=true) const {
      unsigned int code(m_flow.Code(pos));
      if (test && p_part && p_part->GetFlow(pos)!=code) {
	msg_Error()<<"Error in "<<METHOD<<": "<<std::endl
		   <<"   colours do not coincide "
		   <<"("<<code<<" vs. "<<p_part->GetFlow(pos)<<") "
		   <<"for pos = "<<pos<<"."<<std::endl;
	exit(1);
      }
      return code;
    }

    inline void BoostBack(const ATOOLS::Poincare & booster) {
      booster.BoostBack(m_mom);
    }

    inline void Boost(const ATOOLS::Poincare & booster) {
      booster.Boost(m_mom);
    }

    ATOOLS::Particle * GetParticle();
  };
  std::ostream & operator<<(std::ostream & s, const Ladder_Particle & part);
  typedef std::map<double,Ladder_Particle,std::less_equal<double> > LadderMap;


  struct T_Prop {
    T_Prop(const colour_type::code & col,
	   const ATOOLS::Vec4D & q,const double & q02=0.) : 
      m_col(col), m_q(q), m_q2(ATOOLS::dabs(m_q.Abs2())), 
      m_qt2(m_q.PPerp2()), m_q02(q02) {}
    T_Prop(const colour_type::code & col) : 
      m_col(col), m_q(ATOOLS::Vec4D(0.,0.,0.,0.)), m_q2(0.), 
      m_qt2(0.), m_q02(0.) {}
    T_Prop(const ATOOLS::Vec4D & q,const double & q02=0.) : 
      m_col(colour_type::octet), m_q(q), m_q2(ATOOLS::dabs(m_q.Abs2())), 
      m_qt2(m_q.PPerp2()), m_q02(q02) {}
    colour_type::code m_col;
    ATOOLS::Vec4D     m_q;
    double            m_q2, m_qt2, m_q02;
  };
  typedef std::list<T_Prop> TPropList;
  std::ostream & operator<<(std::ostream & s, const T_Prop & tprop);
  std::ostream & operator<<(std::ostream & s, const TPropList & props);

  class Ladder {
  private:
    ATOOLS::Vec4D       m_position;
    double              m_weight, m_Y, m_Ymin, m_Ymax, m_deltaY;
    LadderMap           m_emissions;
    TPropList           m_tprops;
    bool                m_diffractive, m_rescatter, m_harddiffractive;
    double              m_maxkt2, m_minkt2;
    double              m_shat, m_that, m_uhat, m_mu2, m_Yhat, m_DeltaYhat;
    colour_type::code   m_propcol;
    Ladder_Particle   * p_inpart1, * p_inpart2, * p_part1, * p_part2;
    T_Prop            * p_hardestprop;    
    bool                m_enforceup;

    bool FixFirstColours(LadderMap::iterator & liter,
			 int & col1,int & col2,const size_t & fix,
			 TPropList::iterator & citer);
    bool FixIntermediateColours(LadderMap::iterator & liter,
				int & col1,int & col2,size_t & fix,
				TPropList::iterator & citer);
    bool FixLastColours(LadderMap::iterator & liter,
			const int & col1,const int & col2,const size_t & fix,
			TPropList::iterator & citer);
    bool MoreSinglets(const TPropList::iterator & citer);
    void UpdatePropagators();
  public:
    Ladder(const ATOOLS::Vec4D & position);
    ~Ladder();
    // Colour on the ladder
    bool GenerateColourIndices(size_t & fix);
    bool CanReplaceColour(const size_t & oldc,const size_t & newc,
			  const size_t & pos,const bool & inclIS=true);

    // ME reweighting
    inline const double            & Shat()       const { return m_shat; }
    inline const double            & That()       const { return m_that; }
    inline const double            & Uhat()       const { return m_uhat; }
    inline const double            & Yhat()       const { return m_Yhat; }
    inline const double            & DeltaYhat()  const { return m_DeltaYhat; }
    inline const colour_type::code & Propcolour() const { return m_propcol; }
    const bool & IsHardDiffractive() const { return m_harddiffractive; } 
    bool ExtractHardest(); 
    bool ReconstructMEFlavours(ATOOLS::Flavour & i1,ATOOLS::Flavour & i2,
			       ATOOLS::Flavour & o1,ATOOLS::Flavour & o2);
    inline Ladder_Particle * GetMEFSParticle(const int & i) const {
      if (i==1) return p_part1;
      return p_part2;
    }
    inline void SetMEFSParticles(Ladder_Particle * part1,
				 Ladder_Particle * part2) {
      p_part1 = part1; p_part2 = part2;
    }

    // Simple service functions
    bool CheckFourMomentum();
    inline void Reset(const bool & all=true) { 
      m_emissions.clear(); 
      m_tprops.clear();
      if (all) {
	m_deltaY  = m_Ymax = m_Ymin = m_Y = 0.;
	m_maxkt2  = m_minkt2 = m_mu2 = 0.;
	m_diffractive = m_rescatter = m_harddiffractive = false;
	m_weight  = 1.;
	p_inpart1 = NULL;
	p_inpart2 = NULL;
      }
    }
    inline void ResetFS() {
      m_emissions.clear(); 
      m_tprops.clear();
      m_maxkt2 = m_mu2 = 0.;
      m_diffractive = m_harddiffractive = false;
      m_weight = 1.;

    }
    // Global properties of ladders
    const double & Weight()        const { return m_weight; }
    const double & MaxKT2()        const { return m_maxkt2; }
    const double & MinKT2()        const { return m_minkt2; }
    const double & Mu2()           const { return m_mu2; }
    const double & Y()             const { return m_Y; }
    const bool   & IsDiffractive() const { return m_diffractive; } 
    const bool   & IsRescatter()   const { return m_rescatter; } 
    double MRKweight() const;
    inline Ladder_Particle * GetIn(const int & pos) {
      if (pos==0) return p_inpart1;
      return p_inpart2;
    }
    inline Ladder_Particle * GetIn1() const { return p_inpart1; }
    inline Ladder_Particle * GetIn2() const { return p_inpart2; }
    inline void GetInFlavours(ATOOLS::Flavour & flav1,
			      ATOOLS::Flavour & flav2) const {
      flav1 = p_inpart1->m_flav; flav2 = p_inpart2->m_flav;
    }
    inline void GetInFlows(ATOOLS::Flow & flow1,
			      ATOOLS::Flow & flow2) const {
      flow1 = p_inpart1->m_flow; flow2 = p_inpart2->m_flow;
    }
    inline const ATOOLS::Vec4D & Position() const { return m_position; }
    inline void SetWeight(const double & weight)  { m_weight = weight; }
    inline void SetMaxKT2(const double & maxkt2)  { 
      if (maxkt2>m_maxkt2) m_maxkt2 = maxkt2; 
    }
    inline void SetMinKT2(const double & minkt2)  { m_minkt2 = minkt2; }
    inline void SetMu2(const double & mu2)        { m_mu2    = mu2; }
    inline void SetDiffractive(const bool & diff) { m_diffractive = diff; }
    inline void SetRescatter(const bool & resc)   { m_rescatter = resc; }

    inline void SetInParticles(Ladder_Particle * part1,
			       Ladder_Particle * part2) {
      if (part1 && part2) {
	part1->m_IS= part2->m_IS = true;
	if (part1->m_mom.Y()<part2->m_mom.Y()) {
	  p_inpart1 = part1; 
	  p_inpart2 = part2;
	}
	else {
	  p_inpart1 = part2; 
	  p_inpart2 = part1;
	}
	p_inpart1->m_pos = m_position; 
	p_inpart2->m_pos = m_position;
      }
      else { 
	p_inpart1 = NULL; 
	p_inpart2 = NULL; 
      }
    }

    inline T_Prop * GetHardestProp()  { return p_hardestprop; }
    // Organizing the emissions and propagators
    bool SwapColourIndices();
    inline size_t Size()  const { return m_emissions.size(); }
    inline LadderMap * GetEmissions() { return &m_emissions; }
    inline TPropList * GetProps()     { return &m_tprops; }
    inline void MakeDiffractive() {
      if (IsDiffractive()) return;
      TPropList::iterator tit(m_tprops.begin()), softest(m_tprops.begin());
      while (tit!=m_tprops.end()) {
	if (ATOOLS::dabs(tit->m_q2)<ATOOLS::dabs(softest->m_q2)) softest=tit;
	tit++;
      }
      softest->m_col = colour_type::singlet;
      m_diffractive  = true;
    }
    inline void AddRapidity(const double & y) {
      Ladder_Particle part;
      m_emissions[y] = part;
      part.m_pos = m_position;
    }
    inline void DeleteRapidity(const double & y) {
      if (m_emissions.find(y)!=m_emissions.end()) m_emissions.erase(y);
    }
    inline void DeleteRapidity(LadderMap::iterator yiter) {
      if (yiter!=m_emissions.end()) m_emissions.erase(yiter);
    }

    inline void AddParticle(const double & y,Ladder_Particle part) {
      part.m_pos = m_position;
      m_emissions[y] = part;
    }
    inline void SetParticle(LadderMap::iterator yiter,
			    const Ladder_Particle part) {
      yiter->second = part;
    }
    LadderMap::iterator GetEmissionsBegin() { 
      return m_emissions.begin(); 
    }
    LadderMap::iterator GetEmissionsEnd() { 
      return m_emissions.end(); 
    }
    LadderMap::reverse_iterator GetEmissionsRBegin() { 
      return m_emissions.rbegin(); 
    }
    LadderMap::reverse_iterator GetEmissionsREnd() { 
      return m_emissions.rend(); 
    }
    TPropList::iterator GetPropsBegin() { 
      return m_tprops.begin(); 
    }
    TPropList::iterator GetPropsEnd() { 
      return m_tprops.end(); 
    }
    TPropList::reverse_iterator GetPropsRBegin() { 
      return m_tprops.rbegin(); 
    }
    TPropList::reverse_iterator GetPropsREnd() { 
      return m_tprops.rend(); 
    }

    friend std::ostream & operator<<(std::ostream & s, const Ladder &);
  };

  std::ostream & operator<<(std::ostream & s, const Ladder &);
}
#endif
