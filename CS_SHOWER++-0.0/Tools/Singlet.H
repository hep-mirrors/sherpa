#ifndef Singlet_h
#define Singlet_h

#include "Parton.H"
#include "Blob.H"

namespace CS_SHOWER {
  class Singlet : public Parton_List {
    int  m_col;
    bool m_active;
  public :
    Singlet(int col=-1,bool act=true) : m_col(col), m_active(act) {}
    ~Singlet();

    int       SplitParton(PLiter &, Parton *, Parton *);
    Singlet * SplitList(PLiter);
    void      ReshuffleList(PLiter);
    void      Swap(PLiter,PLiter);
    void      ExtractPartons(ATOOLS::Blob *,ATOOLS::Blob *);

    inline bool const IsActive()  const;
    inline void SetActive();
    inline void SetPassive();

    friend std::ostream& operator<<(std::ostream &,Singlet &);
  };

  bool const Singlet::IsActive() const { return m_active; }
  void Singlet::SetActive()      { m_active = true;  }
  void Singlet::SetPassive()     { m_active = false; }

  typedef std::list<Singlet *>   All_Singlets;
  typedef All_Singlets::iterator ASiter;
  std::ostream& operator<<(std::ostream &,All_Singlets &); 
}

#endif
