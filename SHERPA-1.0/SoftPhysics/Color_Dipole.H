#ifndef Color_Dipole_H
#define Color_Dipole_H

#include "Particle_List.H"
#include "Blob.H"
#include <map>

namespace SHERPA {

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

    typedef std::map<ATOOLS::Particle *const,
		     ATOOLS::Flow *> Particle_Flow_Map;

    static std::set<ATOOLS::Particle*> s_partons;

  private:
   
    Particle_Flow_Map      m_flows[qri::size];
    ATOOLS::Particle_List *p_companions;

    ATOOLS::Particle *p_begin[qri::size], *p_end[qri::size];

    bool m_status[qri::size];

    Color_Dipole();

    void SelectCompanion(ATOOLS::Particle *const begin);
    void CollectString(const qri::type type); 
    
    void SetColors(const qri::type type);
    bool AssignColor(Particle_Flow_Map::iterator fit,
		     const unsigned int oldc,const unsigned int newc);

    bool Cat(Color_Dipole *const dipole,const qri::type type);
    void Split(const qri::type type);

    void UnDo(const qri::type type);
    void Clear(const qri::type type);

    friend std::ostream &operator<<(std::ostream &str,
				    const Color_Dipole &info);

  public:

    // constructor
    Color_Dipole(ATOOLS::Particle *const begin,
		 ATOOLS::Particle_List *const companions=NULL);
    
    // destructor
    ~Color_Dipole();

    // member functions
    bool Insert(Color_Dipole *const info,const qri::type type);
    bool Cross(Color_Dipole *const info,const qri::type type);

    bool AssignColors(const qri::type type,const int color);
    void SetColors();

    void UnDo();

    bool Includes(ATOOLS::Particle *const part);

    void Prepend(ATOOLS::Particle *const part,
		 const qri::type type,const bool same=true);
    void Append(ATOOLS::Particle *const part,
		const qri::type type,const bool same=true);

    bool Cat(Color_Dipole *const dipole);

    bool Singlet(const qri::type type) const;
    
    // inline functions
    inline ATOOLS::Particle *const Begin(const qri::type type) const     
    { return p_begin[type]; }
    inline ATOOLS::Particle *const End(const qri::type type) const       
    { return p_end[type]; }

    inline bool Status(const qri::type type) const 
    { return m_status[type]; }

  };// end of class Color_Dipole

  std::ostream &operator<<(std::ostream &str,const Color_Dipole &info);

}// end of namespace SHERPA

#define ANTI(TYPE) (qri::type)(1-(int)TYPE)

#define COLOR(TYPE) 1+(int)TYPE
#define TYPE(COLOR) (qri::type)((int)COLOR-1)

#endif
