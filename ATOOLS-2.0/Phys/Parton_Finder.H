#ifndef Parton_Finder_H
#define Parton_Finder_H

#include "Blob.H"

namespace ATOOLS {

  struct pfc {
    
    enum code{
      color = 1
    };

  };// end of struct pfc

  class Parton_Finder {
  private:
    
    const Particle *m_start, *m_end;

    std::set<btp::code> m_excludeblobs;
    std::set<kf::code>  m_excludeflavours;

    std::set<const Particle*> m_tested;

    std::pair<unsigned int,unsigned int> m_color;

    pfc::code m_criterion;

    void Turn();
    bool Test(const Particle *cur);

    const Particle *FindConstConnectedForward(const Particle *start);
    const Particle *FindConstConnectedBackward(const Particle *start);

    const Particle *FindConstConnected(const Particle *start,bool forward=true);

  public:

    // constructor
    Parton_Finder();

    // member funnctions
    void Clear();

    Particle *FindConnected(const Particle *start,
			    bool forward=true,unsigned int index=0);

    const Particle *FindConstConnected(); 

    // inline functions
    inline void SetStart(Particle *const start) { m_start=start; }
    inline void SetEnd(Particle *const end)     { m_end=end;     }

    inline void Exclude(const btp::code code) { m_excludeblobs.insert(code);    }
    inline void Exclude(const kf::code code)  { m_excludeflavours.insert(code); }

    inline Particle *Start() const { return (Particle*)m_start; }
    inline Particle *End() const   { return (Particle*)m_end;   }

    inline const Particle *ConstStart() const { return m_start; }
    inline const Particle *ConstEnd() const   { return m_end;   }

    inline void SetColour(const unsigned int i,const unsigned int code)
    { m_color=std::pair<unsigned int,unsigned int>(i,code); }

  };// end of class Parton_Finder

}

#endif
