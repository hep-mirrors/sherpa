#ifndef Hadron_Multiplet_H
#define Hadron_Multiplet_H

#include "Hadron_Wave_Function.H"
#include "Flavour.H"
#include <map>
#include <string>


namespace AHADIC {

  class Hadron_Multiplet {
  private:
    std::string        m_name;
    ATOOLS::FlavourSet m_elements;
    double             m_spinweight,m_extraweight,m_weight;
  public:
    Hadron_Multiplet() : 
      m_name(std::string("")),
      m_spinweight(1.),m_extraweight(1.),m_weight(1.) {}
    void SetName(std::string name)                 { m_name        = name;  }
    void AddToElements(ATOOLS::Flavour & fl)       { m_elements.insert(fl); }
    void AddToElements(const ATOOLS::Flavour & fl) { m_elements.insert(fl); }
    void SetSpinWeight(double wt)                  { m_spinweight  = wt;    }
    void SetExtraWeight(double wt)                 { m_extraweight = wt;    }
    void SetWeight(double wt=0.)                   { 
      if (wt==0.) m_weight = m_spinweight * m_extraweight;
             else m_weight = wt;
    }
    std::string          Name()        { return m_name;            }
    ATOOLS::FlavourSet * GetElements() { return &m_elements;       }
    int                  Size()        { return m_elements.size(); }
    double               SpinWeight()  { return m_spinweight;      }
    double               ExtraWeight() { return m_extraweight;     }
    double               Weight()      { return m_weight;          }
  };

  typedef std::map<int,Hadron_Multiplet *>  Hadron_Multiplet_Map;
  typedef Hadron_Multiplet_Map::iterator    Hadron_Multiplet_Miter;

  class All_Hadron_Multiplets {
  private:
    Hadron_WF_Map         * p_wavefunctions;
    Hadron_Multiplet_Map  * p_multiplets;

    void ConstructWaveFunctions();
    Hadron_Wave_Function * ConstructMesonWaveFunction(int,int,int,int);
    Hadron_Wave_Function * ConstructBaryonWaveFunction(int,int,int,int,int);
    void ConstructAntiWaveFunctions();
    void CreateMultiplets();
    void AddMultipletWeights();

    void LookUpAngles(const int angular,const int spin,double & costh,double & sinth);
    bool LookUpWaveFunction(const int kfcode,const int fl1,const int fl2,
			    const Hadron_Wave_Function * wf);
  public:
    All_Hadron_Multiplets();
    ~All_Hadron_Multiplets();

    Hadron_WF_Map        * GetWaveFunctions() { return p_wavefunctions; }

    Hadron_Wave_Function * GetWaveFunction(ATOOLS::Flavour);
    void                   PrintWaveFunctions();
  };
}

#endif
