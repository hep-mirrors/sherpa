#ifndef Hadron_Multiplet_H
#define Hadron_Multiplet_H

#include "Hadron_Wave_Function.H"
#include "Flavour.H"
#include <map>
#include <string>


namespace AHADIC {

  struct Hadron_Multiplet {
    std::string        m_name;
    ATOOLS::FlavourSet m_elements;
    double             m_weight;
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
