#ifndef AHADIC_Tools_Hadronisation_Parameters_H
#define AHADIC_Tools_Hadronisation_Parameters_H

#include "AHADIC++/Tools/Constituents.H"
#include "AHADIC++/Tools/Hadron_Multiplet.H"
#include "AHADIC++/Tools/Transitions.H"
#include "ATOOLS/Math/Vector.H"

#include <map>
#include <string>


namespace ATOOLS {
  class Default_Reader;
}

namespace AHADIC {
  class Soft_Cluster_Handler;

  class Hadronisation_Parameters {
  private:
    double                   m_offset;
    Constituents           * p_constituents;
    All_Hadron_Multiplets  * p_multiplets;
    Single_Transitions     * p_singletransitions;
    Double_Transitions     * p_doubletransitions;
    Soft_Cluster_Handler   * p_softclusters;
    bool                     m_ana;

    StringDoubleMap          m_parametermap;
    SDiter                   m_piter;

    void ReadParameters(std::string,std::string);
    void ReadGeneralSwitches(ATOOLS::Default_Reader & reader);
    void ReadMassParameters(ATOOLS::Default_Reader & reader);
    void ReadPoppingParameters(ATOOLS::Default_Reader & reader);
    void ReadMesonWeights(ATOOLS::Default_Reader & reader);
    void ReadGluonSplittingParameters(ATOOLS::Default_Reader & reader);
    void ReadClusterDecayParameters(ATOOLS::Default_Reader & reader);
    void ReadClusterToMesonParameters(ATOOLS::Default_Reader & reader);
    void ReadDeprecatedParameters(ATOOLS::Default_Reader & reader);
  public:
    Hadronisation_Parameters();
    ~Hadronisation_Parameters();

    void   Init(std::string,std::string);
    double Get(std::string keyword);
    const bool & AnaOn() const { return m_ana; }

    Constituents           * GetConstituents()       { return p_constituents;      }
    All_Hadron_Multiplets  * GetMultiplets()         { return p_multiplets;        }
    Single_Transitions     * GetSingleTransitions()  { return p_singletransitions; }
    Double_Transitions     * GetDoubleTransitions()  { return p_doubletransitions; }
    Soft_Cluster_Handler   * GetSoftClusterHandler() { return p_softclusters;      }

    bool AdjustMomenta(const int,ATOOLS::Vec4D *,const double *);  
  };

  extern Hadronisation_Parameters *hadpars;
}

#endif
