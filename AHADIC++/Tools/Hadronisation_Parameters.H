#ifndef AHADIC_Tools_Hadronisation_Parameters_H
#define AHADIC_Tools_Hadronisation_Parameters_H

#include "AHADIC++/Tools/Constituents.H"
#include "AHADIC++/Tools/Double_Transitions.H"
#include "AHADIC++/Tools/Single_Transitions.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"

#include <map>
#include <string>
#include <cassert>

namespace AHADIC {
  static ATOOLS::Vec4D s_AxisP(1.,0.,0.,1.);
  static ATOOLS::Vec4D s_AxisM(1.,0.,0.,-1.);

  struct control {
  public:
    static long int s_AHAparticles;
    static long int s_AHAclusters;
  };

  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class Hadronisation_Parameters {
  private:
    double               m_offset;
    Constituents       * p_constituents;
    Single_Transitions * p_stransitions;
    Double_Transitions * p_dtransitions;
    
    StringDoubleMap      m_parametermap;
    SDiter               m_piter;
    bool                 m_ana;
    
    void ReadParameters(std::string,std::string);
    void ReadGeneralSwitches(ATOOLS::Data_Reader & dataread);
    void ReadMassParameters(ATOOLS::Data_Reader & dataread);
    void ReadPoppingParameters(ATOOLS::Data_Reader & dataread);
    void ReadMesonWeights(ATOOLS::Data_Reader & dataread);
    void ReadGluonSplittingParameters(ATOOLS::Data_Reader & dataread);
    void ReadClusterDecayParameters(ATOOLS::Data_Reader & dataread);
    void ReadClusterToMesonParameters(ATOOLS::Data_Reader & dataread);
    void ReadDeprecatedParameters(ATOOLS::Data_Reader & dataread);
  public:
    Hadronisation_Parameters();
    ~Hadronisation_Parameters();

    void Init(std::string,std::string);
    double Get(std::string keyword);
    const bool & AnaOn() const { return m_ana; }

    Constituents          * GetConstituents()       { return p_constituents; }
    //All_Hadron_Multiplets * GetMultiplets()         { return p_multiplets; }
    Single_Transitions    * GetSingleTransitions()  { return p_stransitions; }
    Double_Transitions    * GetDoubleTransitions()  { return p_dtransitions; }
    
    bool AdjustMomenta(const int,ATOOLS::Vec4D *,const double *);  
  };

  extern Hadronisation_Parameters * hadpars;

  inline double Lambda(const double & ma2,const double & mb2,
		       const double & mc2,const double & kt2=0.)  {
    assert(ATOOLS::sqr(ma2-mb2-(mc2+kt2))>4.*mb2*(mc2*kt2));
    return sqrt(ATOOLS::sqr(ma2-mb2-(mc2+kt2))-4.*mb2*(mc2+kt2))/(2.*ma2);
  };
}

#endif
