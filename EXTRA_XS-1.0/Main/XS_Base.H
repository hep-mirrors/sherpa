#ifndef XS_Base_H
#define XS_Base_H


#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"
#include "Beam_Handler.H"
#include "Selector.H"
#include "Standard_Selector.H"

#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"

#include <string>

namespace EXTRAXS { class XS_Base; }

namespace EXTRAXS {
  class XS_Base {
  protected:
    int                   nin, nout;
    APHYTOOLS::Flavour  * fl;
    std::string           name;

    std::vector<int>      isr_types;
    std::vector<double>   isr_masses;
    std::vector<double>   isr_widths;

    long int              n;
    double                totalxs, totalerr, totalsum, totalsumsqr, max;
    double                last, lastdxs, lastlumi;
    double                scale, s,t,u;
    int                ** colours;

    double                taumin,taumax;

    APHYTOOLS::Selector_Base    * sel;
    BEAM::Beam_Handler          * beam;
    ISR::ISR_Handler            * isr;
    PHASIC::Phase_Space_Handler * ps; 
    XS_Base                     * selected;
    AMATOOLS::vec4d             * moms;

    void Init(int,int,APHYTOOLS::Flavour * =0);
    virtual void CreateSelector();
  public:
    XS_Base() {};
    virtual ~XS_Base();
    bool SetUpIntegrator(ISR::ISR_Handler *,BEAM::Beam_Handler *);
    bool CalculateTotalXSec();
    bool OneEvent();

    virtual void AddPoint(double);
    virtual void SetTotalXS();
    virtual void SelectOne();
    virtual void DeSelect();

    virtual double Differential(AMATOOLS::vec4d *);
    virtual double Differential2();
    virtual double Differential(double,double,double);
    virtual double DSigma(double,double,double);
    virtual double DSigma2();
    virtual double operator() (double,double,double);

    virtual bool   SetColours(AMATOOLS::vec4d *);
    virtual bool   SetColours(double,double,double);

    virtual void   ISRInfo(int,int &,double &,double &);
    virtual double Scale(AMATOOLS::vec4d *);
    virtual double KFactor(double);

    virtual void SetBeam(BEAM::Beam_Handler *);
    virtual void SetISR(ISR::ISR_Handler *);
    virtual XS_Base * Selected();

    int                        Nin();
    int                        Nout();
    std::string                Name();
    APHYTOOLS::Flavour       * Flavs();
    int                     ** Colours();
    AMATOOLS::vec4d          * Momenta();
    double                     Scale();

    double                     Total();
    double                     Max();
    double                     Last();
    double                     LastXS();
    double                     LastLumi();

    APHYTOOLS::Selector_Base * Selector();
    int                        ISRNumber();

    void SetMomenta(AMATOOLS::vec4d *);
    void SetMax(double);
  };


  inline int                        XS_Base::Nin()       { return nin; }
  inline int                        XS_Base::Nout()      { return nout; }  
  inline std::string                XS_Base::Name()      { return name; }
  inline APHYTOOLS::Flavour       * XS_Base::Flavs()     { return fl; }
  inline int                     ** XS_Base::Colours()   { return colours; }
  inline AMATOOLS::vec4d          * XS_Base::Momenta()   { return moms; }
  inline double                     XS_Base::Scale()     { return scale; }
  inline APHYTOOLS::Selector_Base * XS_Base::Selector()  { return sel; }
  inline int                        XS_Base::ISRNumber() { return isr_types.size(); }


  inline double                     XS_Base::Total()     { return totalxs; }
  inline double                     XS_Base::Max()       { return max; }
  inline double                     XS_Base::Last()      { return last; }
  inline double                     XS_Base::LastXS()    { return lastdxs; }
  inline double                     XS_Base::LastLumi()  { return lastlumi; }


  inline void XS_Base::SetMomenta(AMATOOLS::vec4d * _moms) { moms = _moms; }
  inline void XS_Base::SetMax(double _max)                 { max  = _max; }
  inline void XS_Base::SelectOne()                         { selected = this; }
  inline void XS_Base::DeSelect()                          { selected = 0; }
}





/*! 
  \class XS_Base 
  \brief The mother class of all processes : Single_XS, XS_Group, All_XSes.
  
  Every class, which is to be integrated by the Phase_Space_Integrator is derived 
  from this class or from Process_Base in AMEGIC. 
*/

/*! 
  \var   double   totalxs, totalerr, totalsum, totalsumsqr, max, long int n
  \brief Cross sections and their errors. 

  There's a double strategy : 
  First of all the relative contribution of individual channels is optimized within 
  the multichannel such that the overall variance is minimized. To do this, individual 
  results, variance for each /sa Single_Channel are to be stored within each 
  Single_Channel. 
  The total result and its error is accumulated in the /sa Phase_Space_Integrator and 
  in the full /sa Multi_Channel as well. It should be noted that this /sa Multi_Channel
  can also be a /sa Single_Channel only, both are derived from the /sa Channel_Base. 
  It should be stressed here that, if ISR is on, there are two Multi_Channels,
  one for th final state momenta and one for the initial state piece of the
  integration.
  
  On the other hand, results have to be distributed over the processes that are
  integrated. This is especially true for /sa XS_Groups. There the total result
  coincides with the one of the integration stored in /sa Phase_Space_Integrator
  and its specific /sa Channel_Base. However, this total result has to be distributed
  over the /sa Single_XS(es), therefore they have to accumulate their own
  result, variance, number of points etc..
*/

/*!
  \var double last, lastdxs, lastlumi;
  The results for the partial cross section at the last phase space point.
  This is important for instance for the integration, if one /sa Single_XS
  can be mapped onto another one. Then the last amplitude squared lastdxs of this 
  partner will be used and folded with the PDF (lastlumi) that depends on the 
  particular initial state. We have   last = lastlumi * lastdxs.
*/


#endif
