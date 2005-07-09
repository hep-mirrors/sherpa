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

  /*!
    \file Hadronisation_Parameters.H
    \brief Contains the classes AHADIC::Momenta_Stretcher and AHADIC::Hadronisation_Parameters.
  */

  /*!
    \class Momenta_Stretcher
    \brief Utility class to bring sets of four vectors on different mass-shells.
  */
  /*!
    \class Hadronisation_Parameters
    \brief This class stores all important parameters of the model. 

    This class is made external, and thus is initialised only once in the AHADIC module. It includes a 
    parameter map (a map connecting word-strings with double numbers), the constituent characteristics 
    (accessible through AHADIC::Constituents), the hadron multiplets (AHADIC::All_Hadron_Multiplets), 
    and the static single and double hadron transition probabilities of clusters with a given flavour 
    content (realised by AHADIC::All_Single_Transitions and AHADIC::All_Double_Transitions). These 
    classes are all initialised and, after the run, deleted by the Hadronisation_Parameters.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Momenta_Stretcher.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn bool Momenta_Stretcher::MassThem(const int,ATOOLS::Vec4D *,const double *)
    \brief Puts the momenta on new mass-shells, respecting four-momentum conservation of the overall
    ensemble and the original directions as far as possible. If this succeeds, "true" is returned; if, 
    in contrast, the sum of all new masses is larger than the overall energy available, such a 
    transformation cannot work out and the method yields "false".

    In general the momenta set need to be in its overall c.m. frame for this method to work.
    Then the corresponding algorithm of Rambo is employed. It shifts the masses and uses
    scale transformations in an iterative procedure until convergence is reached. 
    
    For two momenta, the system does not need to be in its c.m. frame. There, the system is boosted 
    into its c.m. frame, the overall energy is calculated and the three-momenta are adjusted, 
    respecting the original directions. Finally, the systemv is boosted back.
  */
  /*!
    \fn bool Momenta_Stretcher::ZeroThem(const int,ATOOLS::Vec4D *)
    \brief Puts the momenta on zero mass, respecting four-momentum conservation of the overall
    ensemble and the original directions as far as possible.

    For two momenta this is straight-forwardly done after going to the two vectors c.m. system
    and calculating the total energy. For an arbitrary number of momenta, the inverse of the
    Rambo algorithm used in Momenta_Stretcher::MassThem(const int,ATOOLS::Vec4D *,const double *) is
    employed.
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
    \var All_Single_Transitions * Hadronisation_Parameters::p_transitions1
    \brief Pointer to the AHADIC::All_Single_Transitions object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var All_Double_Transitions * Hadronisation_Parameters::p_transitions2
    \brief Pointer to the AHADIC::All_Double_Transitions object. It is also initialised only once, in
    the Hadronisation_Parameters and then made available to the rest of the code through a corresponding
    Get-method.
  */
  /*!
    \var Momenta_Stretcher * Hadronisation_Parameters::p_stretcher
    \brief Pointer to the Momenta_Stretcher. It seems at this point this is a bit over-engineered.
  */
  /*!
    \var StringDoubleMap Hadronisation_Parameters::m_parametermap
    \brief 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Hadronisation_Parameters.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Hadronisation_Parameters::
    \brief
  */
}

#endif
