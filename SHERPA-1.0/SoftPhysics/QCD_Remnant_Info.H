#ifndef QCD_Remnant_Info_H
#define QCD_Remnant_Info_H

#include "Particle.H"

namespace SHERPA {

  class QCD_Remnant_Info {
  private:
    
    ATOOLS::Particle *p_this[2], *p_finalstate;
    QCD_Remnant_Info *p_connected[2];

    bool m_anti;

    unsigned int m_oldcolor[2];
    
    static std::set<QCD_Remnant_Info*> m_tested;

    void SelectCompanion();

    bool Find(const QCD_Remnant_Info *info,const bool forward);
    bool Find(const ATOOLS::Particle *particle,const bool forward);

  public:

    // constructor
    QCD_Remnant_Info(ATOOLS::Particle *const particle);
   
    // member functions
    int Find(const QCD_Remnant_Info *info);
    int Find(const ATOOLS::Particle *particle);

    // inline functions
    inline ATOOLS::Particle *&operator()() { return p_this[m_anti]; }
    inline ATOOLS::Particle *&operator->() { return p_this[m_anti]; }

    inline ATOOLS::Particle *&operator++()  { return p_finalstate; }

    inline ATOOLS::Particle *&operator()(const size_t i) { return p_this[i];      }
    inline QCD_Remnant_Info *&operator[](const size_t i) { return p_connected[i]; }

    inline unsigned int &operator-(const size_t i) { return m_oldcolor[i];           }
    inline unsigned int  operator+(const size_t i) { return p_this[i]->GetFlow(1+i); }

    inline bool operator!() const { return m_anti; }


  };// end of class QCD_Remnant_Info

  std::ostream &operator<<(std::ostream &ostr,QCD_Remnant_Info &info);

}// end of namespace SHERPA

#endif
