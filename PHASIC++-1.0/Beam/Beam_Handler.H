#ifndef Beam_Handler_h
#define Beam_Handler_h

#include "Vector.H"
#include "Poincare.H"
#include "Beam_Base.H"
#include "PDF_Base.H"

#include "Pol_Info.H"

namespace BEAM {
  class Beam_Handler {
  protected:
    Beam_Base          ** BeamBase;
    AMATOOLS::Poincare    CMSBoost;
    /*
      mode   
      0          no beam
      1          beam 1
      2          beam 2
      3          beam 1 & beam 2 
    */
    int                   mode; 
    double                mass12, mass22, x1, x2;
    double                exponent[2];
    double                smin,smax,splimits[3], ylimits[2];
    AMATOOLS::Vec4D       fiXVECs[2];
    std::string           name, type;
  public:
      Beam_Handler(int *,APHYTOOLS::Flavour *,double *,APHYTOOLS::Flavour *,double *,double *);
      ~Beam_Handler();

    bool   CheckConsistency(APHYTOOLS::Flavour *,APHYTOOLS::Flavour *);
    void   SetBeamMasses(APHYTOOLS::Flavour * beams); 

    bool   MakeBeams(AMATOOLS::Vec4D *,double,double);
    bool   CalculateWeight(double);
    double Weight(APHYTOOLS::Flavour *);
    double Weight();

    void   BoostInCMS(AMATOOLS::Vec4D *,int);
    void   BoostInLab(AMATOOLS::Vec4D *,int);

    void   SetName(std::string _name) { name         = _name; }
    void   SetSprimeMin(double _spl)  { splimits[0]  = _spl; }
    void   SetSprimeMax(double _spl)  { splimits[1]  = _spl; }
    void   SetYMin(double _yl)        { ylimits[0]   = _yl; }
    void   SetYMax(double _yl)        { ylimits[1]   = _yl; }


    std::string     Name()            { return name; }
    std::string     Type()            { return type; }
    double        * SprimeRange()     { return splimits; }
    double        * YRange()          { return ylimits; }
    double          SprimeMin()       { return splimits[0]; }
    double          SprimeMax()       { return splimits[1]; }
    double          YMin()            { return ylimits[0]; }
    double          YMax()            { return ylimits[1]; }
    double          Upper1()          { return BeamBase[0]->Xmax(); }
    double          Upper2()          { return BeamBase[1]->Xmax(); }
    double          Peak()            { return BeamBase[0]->Peak() *
					       BeamBase[1]->Peak(); }
    double          Exponent(int i)   { return exponent[i]; }
    int             On()              { return mode; }
  };
}





/*!
  \class Beam_Handler 
  \brief Manager of all Beam.
  This class manages all initial state radiation (Beam) according to the strategy 
  named in Beam.dat. It is initialized from Amegic and handed over to all Processes,
  both groups and individual ones. Before oming into full efect during integration
  or event generation, it initalises suitable Beam treatment (Beam_Bases) for both 
  beams and uses them to generate corresponding Beam weights. 
*/

/*!
  \var   Beam_Base * BeamBase[2]
  \brief The two Beam bases, one for each beam.
*/

/*!
  \fn    Beam_Handler(APHYTOOLS::Flavour *)
  \brief The explicit constructor doing a proper initialisation
  The constructor is a pretty simple object. Entering with the two beam flavours and
  a default s-range as arguments, the actual s-range will be deduced by comparing with
  the range given in the "Beam.dat" file. Additionally, the Beam strategy (none or SF)
  is read off there and corresponding Beam_Bases are initialized.
  We might have to think about "hybrid" stuff, e.g. DIS with no Beam off the electron
  and the like ...
*/

/*!
  \fn     bool CalculateWeight(double)
  This method calculates the parton densities according to the two pdf's depending
  on the external x-values and the scale, which is passed as an argument.
*/

/*!
  \fn     bool CalculateWeight2(double)
  This method calculates the parton densities according to the two pdf's depending
  on the external x-values and the scale, which is passed as an argument. The main point
  here is that the x-values are interchanged to allow for flin[0] coming from Beam[1]
  and vice versa.
*/

/*!
  \fn    double Weight(APHYTOOLS::Flavour * fl)
  \brief The weight corresponding to CalculateWeight(double).
  Here the a-priori weight (=1) is multiplied with the actual value of the pdf's for the 
  specific initial state partons and the flux.
*/

/*!
  \fn    double Weight2(APHYTOOLS::Flavour * fl)
  \brief The weight corresponding to CalculateWeight2(double).
  Here the a-priori weight (=1) is multiplied with the actual value of the pdf's for the 
  specific initial state partons and the flux.
*/
 
#endif



