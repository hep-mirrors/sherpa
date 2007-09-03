#ifndef Jet_Finder_Base_H
#define Jet_Finder_Base_H

#include "Selector.H"

namespace ATOOLS {
  
  std::vector<int> ID(size_t id);

  size_t IDCount(size_t id);

  template <typename Type>
  size_t ID(const std::vector<Type> &ids);

  class Jet_Finder_Base: public Selector_Base {
  protected:

    double m_ycut, m_gycut, m_cycut, m_gcycut, m_delta_r;
    int    m_type;

    std::string m_cuttag;

  public:

    // constructor
    Jet_Finder_Base();

    // destructor
    ~Jet_Finder_Base();

    // member functions
    void SetDeltaR(double dr);

    virtual void FillCombinations();

    // inline functions
    inline void SetType(const int mode) { m_type=mode; }

    inline int  Type() const { return m_type; }

    inline void SetYcut(const std::string &ycut) { m_cuttag=ycut; }

    inline double DeltaR() const         { return m_delta_r; }
    inline double Ycut() const           { return m_ycut;    }
    inline double GlobalYcut() const     { return m_gycut;   }
    inline double CoreYcut() const       { return m_cycut;   }
    inline double GlobalCoreYcut() const { return m_gcycut;  }

  };

}

#endif








