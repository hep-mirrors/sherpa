#ifndef Single_Broker_H
#define Single_Broker_H

#include "Message.H"
#include <string>

#include "XS_Base.H"
#include "Process_Base.H"

using namespace AMEGIC;
using namespace PHASIC;

namespace PHASIC {
  
  class Single_Broker: public AMEGIC::Process_Base {
  protected:
    
    EXTRAXS::XS_Base   * xsec;

    AMEGIC::Point       * fake;
    double                s,t,u;
    double                scale;

    double                taumin, taumax;

  public:
    
    Single_Broker(EXTRAXS::XS_Base *);
    Single_Broker(int, int, APHYTOOLS::Flavour *, ISR::ISR_Handler *, BEAM::Beam_Handler *,
		  APHYTOOLS::Selector_Data *, int, int); 
    ~Single_Broker();

    EXTRAXS::XS_Base * GetXS();

    // redefined methods of AMEGIC::Process_Base 

    void   Initialize(APHYTOOLS::Selector_Data * _seldata);
    int    InitAmplitude(AMEGIC::Topology *,AMATOOLS::Vec4D *&,std::vector<double> &,
			 std::vector<AMEGIC::Single_Process *> &)            { return 0; }
    void   InitAnalysis(std::vector<APHYTOOLS::Primitive_Observable_Base *>) { return;   }
    void   UpdateCuts(double,double) { return; }
    void   InitCuts()                { return; }

    virtual void SetName(std::string _name) { name = _name;        }
    virtual void SetResDir(std::string)     { return;              }
    virtual void SetAtoms(bool _atoms)      { atoms = _atoms;      }
    virtual void SetTables(bool _tables)    { tables = 0;  return; }

    void SetNStrong(int) { return; }
    void SetNEWeak(int)  { return; }
    void Empty()         { return; }

    virtual bool LookUpXSec(double,bool,std::string) { return 0; }
    virtual void RescaleXSec(double)                 { return;   }    
    virtual bool PrepareXSecTables()                 { return 0; }
    
    int             NStrong()                              { return 0;          }
    int             NEWeak()                               { return 0;          }
    int             NumberOfDiagrams()                     { return IS_XS_FLAG; }
    AMEGIC::Point * Diagram(int i)                         { return fake;       }
    bool            IsFreeOfFourVertex(AMEGIC::Point * _p) { return 0;          }
    void            PrintDifferential()                    { return;            }

    AMEGIC::Process_Base  * Partner()  { return this; }
    AMEGIC::Process_Base  * Selected() { return this; }

    void Add(AMEGIC::Process_Base *) { return;   }
    bool SetUpIntegrator()           { return 0; }
    int  Size()                      { return 0; }
 
    void SetScale(double _scale)      { scale = _scale; }

    AMEGIC::Process_Base * operator[](int idx) { return 0; }

    void SelectOne() { selected = this; } 
    void DeSelect()  { selected = 0;    } 

    // interfaces to EXTRAXS::XS_Base

//     bool SetColours(AMATOOLS::Vec4D *);
//     bool SetColours(double, double, double);

    void CreateSelector() { sel = new APHYTOOLS::No_Selector(); }
    bool SetUpIntegrator(ISR::ISR_Handler * ,BEAM::Beam_Handler * );

    void AddPoint(const double);   
    bool CalculateTotalXSec();
    void SetTotalXS(int);
    bool OneEvent();

    double Differential(AMATOOLS::Vec4D *);
    double Differential2();
    double DSigma(AMATOOLS::Vec4D *,bool);
    double DSigma2();
    bool   SameEvent();
    double WeightedEvent();

    double KFactor(double);
    
  };

  inline EXTRAXS::XS_Base * Single_Broker::GetXS() { return xsec; }

}

#endif

/*!
  This class is meant to be an interface between AMEGIC::Single_Process and 
  EXTRAXS::Single_XS. At the moment it is directly set up by PHASIC::Broker_Group.
  
  A reference to the correspnding Single_XS is stored in 'xsec' and is used 
  to evaluate differential cross sections, scales and K-factors.

  All methods are described within AMEGIC::Process_Base while not self - explanatory.
 */
