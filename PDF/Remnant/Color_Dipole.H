#ifndef PDF_Remnant_Color_Dipole_H
#define PDF_Remnant_Color_Dipole_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob.H"
#include <map>

namespace PDF {

  struct qri {

    enum type {
      real = 0,
      anti = 1,
      size = 2
    };

  };// end of struct qri

  std::ostream &operator<<(std::ostream &str,const qri::type type);

  class Color_Dipole {
  public:

    typedef std::vector<Color_Dipole*>  Dipole_Vector;
    typedef std::set<ATOOLS::Particle*> Particle_Set;

    typedef std::map<ATOOLS::Particle *const,
		     ATOOLS::Flow *> Particle_Flow_Map;

  private:

    static Dipole_Vector s_dipoles;
    static Particle_Set  s_partons;

    static Particle_Flow_Map s_flows, s_oldflows;
   
    Particle_Flow_Map      m_flows[qri::size];
    ATOOLS::Particle_List *p_companions;

    ATOOLS::Particle *p_begin[qri::size], *p_end[qri::size];
    Color_Dipole     *p_next[qri::size];

    Color_Dipole();

    void SelectCompanion(ATOOLS::Particle *const begin);
    void CollectString(const qri::type type); 

    bool AssignColor(Particle_Flow_Map::iterator fit,
		     const unsigned int oldc,const unsigned int newc);
    bool AssignColors(const qri::type type,const int color);

    void DetectLoop(const qri::type type);

    friend std::ostream &operator<<(std::ostream &str,
				    const Color_Dipole &info);

  public:

    // constructor
    Color_Dipole(ATOOLS::Particle *const begin,
		 ATOOLS::Particle_List *const companions=NULL);
    
    // destructor
    ~Color_Dipole();

    // member functions
    void CollectStrings(); 
    void DetectLoops();
    
    bool Cat(Color_Dipole *const dipole,const qri::type type);

    bool Singlet(const qri::type type) const;
    bool Connected(const Color_Dipole *dipole,const qri::type type,
		   const size_t catcher=0) const;

    static void SetAllColors();
    static void ResetAllColors();
    static void ClearAll();

    // inline functions
    inline ATOOLS::Particle *const Begin(const qri::type type) const     
    { return p_begin[type]; }
    inline ATOOLS::Particle *const End(const qri::type type) const       
    { return p_end[type]; }

    inline Color_Dipole *const Next(const qri::type type) const       
    { return p_next[type]; }

  };// end of class Color_Dipole

  std::ostream &operator<<(std::ostream &str,const Color_Dipole &info);

}// end of namespace PDF

#define ANTI(TYPE) (qri::type)(1-(int)TYPE)

#define COLOR(TYPE) 1+(int)TYPE
#define TYPE(COLOR) (qri::type)((int)COLOR-1)

#endif
