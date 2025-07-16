#ifndef AHADIC_Tools_Hadron_Multiplet_H
#define AHADIC_Tools_Hadron_Multiplet_H

#include "AHADIC++/Tools/Wave_Function.H"

namespace AHADIC {
  // This is just a container class - it does nothing.
  class Hadron_Multiplet {
  private:
    std::set<ATOOLS::Flavour> m_elements;
    double                    m_spinweight,m_extraweight,m_weight;
  public:
    Hadron_Multiplet() : 
      m_spinweight(1.),m_extraweight(1.),m_weight(1.) {}
    void AddToElements(ATOOLS::Flavour fl) { m_elements.insert(fl); }
    void SetSpinWeight(double wt)          { m_spinweight  = wt; }
    void SetExtraWeight(double wt)         { m_extraweight = wt; }
    void SetWeight(double wt=0.) { 
      if (wt==0.) m_weight = m_spinweight * m_extraweight;
             else m_weight = wt;
    }
    std::set<ATOOLS::Flavour> * GetElements() { return &m_elements;       }
    int                         Size()        { return m_elements.size(); }
    double                      SpinWeight()  { return m_spinweight;      }
    double                      ExtraWeight() { return m_extraweight;     }
    double                      Weight()      { return m_weight;          }
  };

  typedef std::map<std::string,Hadron_Multiplet *>  Hadron_Multiplet_Map;
}

#endif
