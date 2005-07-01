#ifndef Hadronisation_Parameters_H
#define Hadronisation_Parameters_H

#include "Constituents.H"
#include "Hadron_Multiplet.H"
#include "Transitions.H"
#include "Vector.H"

#include <map>
#include <string>




namespace AHADIC {
  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class Momenta_Stretcher {
  public:
    Momenta_Stretcher() {}
    bool MassThem(const int,ATOOLS::Vec4D *,const double *);
    bool ZeroThem(const int,ATOOLS::Vec4D *);
  };

  class Hadronisation_Parameters {
  private:
    double                   m_offset;
    Constituents           * p_constituents;
    All_Hadron_Multiplets  * p_multiplets;
    All_Single_Transitions * p_transitions1;
    All_Double_Transitions * p_transitions2;
    Momenta_Stretcher      * p_stretcher;
    StringDoubleMap          m_parametermap;
    SDiter                   m_piter;

    void   ReadParameters(std::string,std::string);
  public:
    Hadronisation_Parameters();
    ~Hadronisation_Parameters();

    void   Init(std::string,std::string);
    double Get(std::string keyword);

    Constituents           * GetConstituents()      { return p_constituents; }
    All_Hadron_Multiplets  * GetMultiplets()        { return p_multiplets;   }
    All_Single_Transitions * GetSingleTransitions() { return p_transitions1; } 
    All_Double_Transitions * GetDoubleTransitions() { return p_transitions2; } 

    bool AdjustMomenta(const int,ATOOLS::Vec4D *,const double *);  
  };

  extern Hadronisation_Parameters hadpars;
}

#endif
