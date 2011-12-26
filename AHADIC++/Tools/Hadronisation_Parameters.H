#ifndef AHADIC_Tools_Hadronisation_Parameters_H
#define AHADIC_Tools_Hadronisation_Parameters_H

#include "AHADIC++/Tools/Constituents.H"
#include "AHADIC++/Tools/Hadron_Multiplet.H"
#include "AHADIC++/Tools/Transitions.H"
#include "AHADIC++/Tools/Strong_Coupling.H"
#include "AHADIC++/Tools/Dipole_Splitter.H"
#include "ATOOLS/Math/Vector.H"

#include <map>
#include <string>




namespace AHADIC {
  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class Soft_Cluster_Handler;

  class Hadronisation_Parameters {
  private:
    double                   m_offset;
    Constituents           * p_constituents;
    All_Hadron_Multiplets  * p_multiplets;
    Single_Transitions     * p_singletransitions;
    Double_Transitions     * p_doubletransitions;
    Soft_Cluster_Handler   * p_softclusters;
    leading::code            m_leading;
    PTOrder::code            m_ptorder;
    ZForm::code              m_zform;
    asform::code             m_asform;
    Strong_Coupling        * p_coupling; 
    Dipole_Splitter        * p_splitter;

    StringDoubleMap          m_parametermap;
    SDiter                   m_piter;

    void   ReadParameters(std::string,std::string);
  public:
    Hadronisation_Parameters();
    ~Hadronisation_Parameters();

    void   Init(std::string,std::string);
    double Get(std::string keyword);

    Constituents           * GetConstituents()       { return p_constituents;      }
    All_Hadron_Multiplets  * GetMultiplets()         { return p_multiplets;        }
    Single_Transitions     * GetSingleTransitions()  { return p_singletransitions; }
    Double_Transitions     * GetDoubleTransitions()  { return p_doubletransitions; }
    Strong_Coupling        * GetCoupling()           { return p_coupling;          }
    Dipole_Splitter        * GetSplitter()           { return p_splitter;          }
    Soft_Cluster_Handler   * GetSoftClusterHandler() { return p_softclusters;      }

    bool AdjustMomenta(const int,ATOOLS::Vec4D *,const double *);  
  };

  extern Hadronisation_Parameters *hadpars;

  /*!
    \file Hadronisation_Parameters.H
    \brief Contains the class AHADIC::Hadronisation_Parameters.
  */
  /*!
    \class Hadronisation_Parameters
    \brief This class stores all important parameters of the model. 

    This class is made external, and thus is initialised only once in the AHADIC module. It includes a 
    parameter map (a map connecting word-strings with double numbers), the constituent characteristics 
    (accessible through AHADIC::Constituents), the hadron multiplets (AHADIC::All_Hadron_Multiplets), 
    the static single and double hadron transition probabilities of clusters with a given flavour 
    content (realised by AHADIC::Single_Transitions and AHADIC::Double_Transitions), the
    AHADIC::Dipole_Splitter relevant for the breakup of clusters and the decay of the incident
    gluons into constituents, and the AHADIC::Soft_Cluster_Handler which chacks for and manages
    the transition of clusters into hadrons.  These classes are all initialised only once,
    made accessible through the Hadronisation_Parameters, and, after the run, deleted by the 
    Hadronisation_Parameters.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Hadronisation_Parameters.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var double Hadronisation_Parameters::m_offset
    \brief Not used at the moment.
  */
  /*!
    \var Constituents * Hadronisation_Parameters::p_constituents
    \brief Pointer to the AHADIC::Constituents. They are initialised once in the Hadronisation_Parameters
    and then made available to all other classes by the corresponding Get-method. To create them with
    proper masses etc., before hadronisation parameters are read in and then used.
  */
  /*!
    \var All_Hadron_Multiplets * Hadronisation_Parameters::p_multiplets
    \brief Pointer to the AHADIC::All_Hadron_Multiplets object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var All_Single_Transitions * Hadronisation_Parameters::p_singletransitions
    \brief Pointer to the AHADIC::All_Single_Transitions object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var All_Double_Transitions * Hadronisation_Parameters::p_doubletransitions
    \brief Pointer to the AHADIC::All_Double_Transitions object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var Dipole_Splitter * Hadronisation_Parameters::p_splitter
    \brief Pointer to the AHADIC::Dipole_Splitter object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var Soft_Cluster_Handler * Hadronisation_Parameters::p_softclusters
    \brief Pointer to the AHADIC::Soft_Cluster_Handler object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var StringDoubleMap Hadronisation_Parameters::m_parametermap
    \brief A map connecting keyword (wordstrings) with parameters of the model, which in turn are
    accessible through the method Hadronisation_Parameters::Get(std::string keyword).  For a list
    of parameters, cf. AHADI::Ahadic. 
  */
  /*!
    \var SDiter Hadronisation_Parameters::m_piter
    \brief The iterator corresponding to Hadronisation_Parameters::m_parametermap.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Hadronisation_Parameters.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Hadronisation_Parameters::Hadronisation_Parameters()
    \brief A void constructor.
  */
  /*!
    \fn Hadronisation_Parameters::~Hadronisation_Parameters()
    \brief Destructor, deleting the Hadronisation_Parameters::p_constituents, 
    Hadronisation_Parameters::p_multiplets, Hadronisation_Parameters::p_singletransiions,
    Hadronisation_Parameters::p_doubletransitions, the Hadronisation_Parameters::p_splitter, 
    and the Hadronisation_Parameters::p_softclusters.
  */
  /*!
    \fn void Hadronisation_Parameters::Init(std::string,std::string)
    \brief Calls the Hadronisation_Parameters::ReadParameters(std::string,std::string)
    method with the arguments of path and filename, and initialises the 
    Hadronisation_Parameters::p_constituents, Hadronisation_Parameters::p_multiplets, 
    Hadronisation_Parameters::p_singletransiions, Hadronisation_Parameters::p_doubletransitions, 
    the Hadronisation_Parameters::p_splitter,
    and the Hadronisation_Parameters::p_softclusters.
  */
  /*!
    \fn void Hadronisation_Parameters::ReadParameters(std::string,std::string)
    \brief Reading the parameter file and filling the Hadronisation_Parameters::m_parametermap.
  */
  /*!
    \fn double Hadronisation_Parameters::Get(std::string keyword)
    \brief Returns the actual vlaue of a parameter, given by the keyword.
  */
  /*!
    \fn Constituents * Hadronisation_Parameters::GetConstituents()
    \brief Returns the pointer to the AHADIC::Constituents, 
    Hadronisation_Parameters::p_constituents.
  */
  /*!
    \fn All_Hadron_Multiplets * Hadronisation_Parameters::GetMultiplets()
    \brief  Returns the pointer to the AHADIC::All_Hadron_Multiplets,
    Hadronisation_Parameters::p_multiplets.
  */
  /*!
    \fn Single_Transitions * Hadronisation_Parameters::GetSingleTransitions()
    \brief  Returns the pointer to the AHADIC::Single_Transitions,
    Hadronisation_Parameters::p_singletransitions.
  */
  /*!
    \fn Double_Transitions * Hadronisation_Parameters::GetDoubleTransitions()
    \brief  Returns the pointer to the AHADIC::Double_Transitions, 
    Hadronisation_Parameters::p_singletransitions.
  */
  /*!
    \fn Dipole_Splitter * Hadronisation_Parameters::GetSplitter()
    \brief  Returns the pointer to the AHADIC::Dipole_Splitter, Hadronisation_Parameters::p_splitter.
  */
  /*!
    \fn Soft_Cluster_Handler * Hadronisation_Parameters::GetSoftclusterHandler()
    \brief  Returns the pointer to the AHADIC::Soft_Cluster_Handler, 
    Hadronisation_Parameters::p_softclusters.
  */
  /*!
    \fn bool Hadronisation_Parameters::AdjustMomenta(const int,ATOOLS::Vec4D *,const double *)
    \brief Puts a number of vectors on their new mass-shells.  In so doing, the 
    ATOOLS::Momenta_Stretcher is being used. 

    The algorithm for this putting of momenta onto a new mass-shell is to first put them
    all on a mass-shell of mass 0 (by inverting the Rambo->Massive trick) and then put them
    onto their new masses (by the direct Rambo->Massive trick).  Of course, this is an
    iterative process with an error of roughly 1.e-6, because square roots are involved.
  */
}

#endif
