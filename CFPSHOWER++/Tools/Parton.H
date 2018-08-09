#ifndef CFPShower_Tools_Parton_H
#define CFPShower_Tools_Parton_H

#include "CFPSHOWER++/Tools/Color.H"
#include "CFPSHOWER++/Tools/Weight.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <iostream>
#include <list>
#include <algorithm>

namespace CFPSHOWER {
  class Splitting;

  struct BranchingWeight {
    double m_t, m_weight;
    inline BranchingWeight(const double & t=-1.0,const double & weight=1.) :
      m_t(t), m_weight(weight) {}
  };
  typedef std::vector<BranchingWeight>            Weight_Vector;
  typedef std::map<const Parton *, Weight_Vector> PWV_Map;

  
  class Parton {
  private:
    static size_t s_cid, s_cnt;
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_mom, m_sp;
    Color           m_color;
    size_t          m_id, m_beam;
    bool            m_on;
    
    Parton        * p_in;
    std::vector<Parton *> m_offsprings;
    std::list<Parton *>   m_spectators;
    PWV_Map         m_specweights;
  public:
    Parton(const ATOOLS::Flavour & m_flav,const ATOOLS::Vec4D & mom,
	   const Color & color=Color(),const size_t & beam=0,const size_t & id=0);
    ~Parton();

    inline const bool              On()                  const { return m_on; }
    inline const ATOOLS::Flavour & Flav()                const { return m_flav; }
    inline const ATOOLS::Vec4D   & Mom()                 const { return m_mom; }
    inline const size_t          & Beam()                const { return m_beam; }
    inline const size_t          & Id()                  const { return m_id; }
    inline Parton                * Out(const size_t & i) const { return m_offsprings[i]; }
    inline Parton                * In()                  const { return p_in; }
    inline Color                 & GetColor()                  { return m_color; }
    inline const PWV_Map         & GetWeights()                { return m_specweights; }
    inline const std::list<Parton *>   * GetSpectators() const { return &m_spectators; }
    inline const std::vector<Parton *> * GetOffsprings() const { return &m_offsprings; }
    
    inline bool FindSpectator(Parton * parton)    {
      return std::find(m_spectators.begin(),m_spectators.end(),parton)!=m_spectators.end();
    }
    inline void AddSpectator(Parton * parton)     { m_spectators.push_back(parton); }
    inline void RemoveSpectator(Parton * parton)  { m_spectators.remove(parton); }
    inline void AddOffspring(Parton * offspring)  { m_offsprings.push_back(offspring); }
    inline void SetMom(const ATOOLS::Vec4D &mom)  { m_mom=mom; }
    inline void SetFlav(const ATOOLS::Flavour &f) { m_flav=f; }
    inline void SetBeam(const int & beam)         { m_beam=beam; }
    inline void SetId(const size_t &id)           { m_id=id; }
    inline void SetColor(const Color & color)     { m_color = color; }
    inline void SetOut(int i,Parton *const p)     { m_offsprings[i]=p; }
    inline void SetIn(Parton *const p)            { p_in=p; }
    inline void SwitchOff()                       { m_on = false; }       
    inline void SwitchOn()                        { m_on = true; }
    inline void ClearWeights()                    { m_specweights.clear(); }

    const double GetXB() const;
    
    double GetWeight(const double & t) const;
    void   AddWeight(const Splitting & split,const bool & accept);
    
    inline static void Reset(const size_t & id=0) { s_cid = id; }
    inline static size_t Count()                  { return s_cnt; }

  };

  std::ostream &operator<<(std::ostream & s,Parton & part);
  typedef std::list<Parton*> Parton_List;
  typedef std::vector<Parton*> Parton_Vector;
}

#endif
