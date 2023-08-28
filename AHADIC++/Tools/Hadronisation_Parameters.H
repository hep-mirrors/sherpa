#ifndef AHADIC_Tools_Hadronisation_Parameters_H
#define AHADIC_Tools_Hadronisation_Parameters_H

#include "AHADIC++/Tools/Constituents.H"
#include "AHADIC++/Tools/Double_Transitions.H"
#include "AHADIC++/Tools/Single_Transitions.H"
#include "ATOOLS/Math/Vector.H"

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

  class Hadronisation_Parameters {
  private:
    int                  m_shower, m_nvariations;
    double               m_offset;
    Constituents       * p_constituents;
    Single_Transitions * p_stransitions;
    Double_Transitions * p_dtransitions;

    std::map<std::string,int>    m_switchmap;
    std::map<std::string,double> m_parametermap;
    std::map<std::string,std::vector<double>> m_parametermap_vecs;

    void ReadParameters();
    void ReadGeneralSwitches();
    void ReadPoppingParameters();
    void ReadMesonWeights();
    void ReadSplittingParameters();
    void ReadClusterToMesonPSParameters();
    void CheckAndPad();
  public:
    Hadronisation_Parameters();
    ~Hadronisation_Parameters();

    void Init(std::string);
    const double Get(std::string keyword) const;
    const std::vector<double> GetVec(std::string) const;
    const int    Switch(std::string keyword) const;

    Constituents          * GetConstituents()       { return p_constituents; }
    //All_Hadron_Multiplets * GetMultiplets()         { return p_multiplets; }
    Single_Transitions    * GetSingleTransitions()  { return p_stransitions; }
    Double_Transitions    * GetDoubleTransitions()  { return p_dtransitions; }

    bool AdjustMomenta(const int,ATOOLS::Vec4D *,const double *);
    int NumberOfVariations() { return m_nvariations; }
  };

  extern Hadronisation_Parameters * hadpars;

  inline double Lambda(const double & ma2,const double & mb2,
		       const double & mc2,const double & kt2=0.)  {
    assert(ATOOLS::sqr(ma2-mb2-(mc2+kt2))>4.*mb2*(mc2*kt2));
    return sqrt(ATOOLS::sqr(ma2-mb2-(mc2+kt2))-4.*mb2*(mc2+kt2))/(2.*ma2);
  };
}

#endif
