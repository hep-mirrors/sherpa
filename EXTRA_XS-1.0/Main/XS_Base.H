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

namespace AMEGIC { class Process_Group; }

namespace EXTRAXS {

  class Single_XS;

  class XS_Base {
  protected:

    int                    nin, nout;

    APHYTOOLS::Flavour   * fl;
    std::string            name;

    int                    massive;         

    double                 s, t, u;
    double                 scale;

    int                 ** colours;

    XS_Base              * selected;

    std::vector<int>    isr_types;
    std::vector<double> isr_masses;
    std::vector<double> isr_widths;

    void Init(int,int,APHYTOOLS::Flavour * =0);

  public:

    XS_Base() {};
    virtual ~XS_Base();

    virtual double operator()(double,double,double);

    virtual bool SetColours(AMATOOLS::Vec4D *);
    virtual bool SetColours(double, double, double);

    virtual XS_Base * Selected();

    virtual int Size() { return 1; }
    virtual XS_Base * operator[](int idx) { return this; }

    int                  Nin();
    int                  Nout();
    std::string          Name();
    APHYTOOLS::Flavour * Flavs();
    void                 SetFlav(int,APHYTOOLS::Flavour);
    int               ** Colours();
    double               Scale();
    virtual double       KFactor(double _scale) { return 1; }

    int  ISRNumber() { return isr_types.size(); }
    void ISRInfo(int,int &,double &,double &);
    void SetISRTypes(APHYTOOLS::Flavour *);

    virtual void MakeBroker(ISR::ISR_Handler *, BEAM::Beam_Handler *,
			    APHYTOOLS::Selector_Data *,AMEGIC::Process_Group *) { }
    virtual void Initialize(ISR::ISR_Handler *, BEAM::Beam_Handler *,
			    APHYTOOLS::Selector_Data *, AMEGIC::Process_Group *) { }

  };

  inline int                    XS_Base::Nin()     { return nin;     }
  inline int                    XS_Base::Nout()    { return nout;    }  
  inline std::string            XS_Base::Name()    { return name;    }
  inline APHYTOOLS::Flavour   * XS_Base::Flavs()   { return fl;      }
  inline int                 ** XS_Base::Colours() { return colours; }
  inline double                 XS_Base::Scale()   { return scale;   }

  inline void                   XS_Base::SetFlav(int i,APHYTOOLS::Flavour _fl) { fl[i] = _fl; }

}

#endif

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


