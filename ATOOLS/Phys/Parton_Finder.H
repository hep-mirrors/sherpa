#ifndef ATOOLS_Phys_Parton_Finder_H
#define ATOOLS_Phys_Parton_Finder_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob.H"

namespace ATOOLS {

  class Parton_Tester {
  public:
    
    // destructor
    virtual ~Parton_Tester();

    // member functions
    virtual void Turn();
    virtual bool Test(const Particle *parton) const;

  };// end of class Parton_Tester

  class Parton_Finder {
  private:
    
    const Particle *m_start, *m_end;

    std::set<btp::code> m_excludeblobs;
    std::set<kf_code>  m_excludeflavours;

    Parton_Tester *p_criterion;
    Particle_List  m_track;

    bool m_forward;

    void Turn();
    bool Test(const Particle *cur);

    const Particle *FindConstConnectedForward(const Particle *start);
    const Particle *FindConstConnectedBackward(const Particle *start);

    const Particle *FindConstConnected(const Particle *start,
				       bool forward=true);

  public:

    // constructor
    Parton_Finder(Parton_Tester &criterion);

    // member funnctions
    void Clear();

    Particle *FindConnected(const Particle *start,
			    bool forward=true);

    const Particle *FindConstConnected(); 

    // inline functions
    inline void SetStart(Particle *const start) { m_start=start; }
    inline void SetEnd(Particle *const end)     { m_end=end;     }

    inline void Exclude(const btp::code code) 
    { m_excludeblobs.insert(code); }
    inline void Exclude(const kf_code code)  
    { m_excludeflavours.insert(code); }

    inline Particle *Start() const { return (Particle*)m_start; }
    inline Particle *End() const   { return (Particle*)m_end;   }

    inline const Particle *ConstStart() const { return m_start; }
    inline const Particle *ConstEnd() const   { return m_end;   }

    inline const Particle_List &Track() { return m_track; }

  };// end of class Parton_Finder

}

#endif
