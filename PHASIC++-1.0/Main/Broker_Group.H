#ifndef Broker_Group_H
#define Broker_Group_H

#include "Process_Base.H"
#include "Single_Broker.H"
#include "XS_Base.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <vector>

namespace PHASIC {
  
  class Broker_Group : public AMEGIC::Process_Base {
  protected:   

    std::vector<AMEGIC::Process_Base *> brokers;
    
    EXTRAXS::XS_Base * xsecs;
    AMEGIC::Point    * fake;

    double s, t, u;
    double taumin, taumax;

  public:
    
    Broker_Group(EXTRAXS::XS_Base *);
    ~Broker_Group();

    EXTRAXS::XS_Base * GetXS();

    // redefined methods of AMEGIC::Process_Base()

    void   Initialize(APHYTOOLS::Selector_Data * _seldata);
    int    InitAmplitude(AMEGIC::Topology *,AMATOOLS::vec4d *&,std::vector<double> &,
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
    AMEGIC::Process_Base  * Selected();

    void CreateSelector() { sel = new APHYTOOLS::No_Selector(); }

    void Add(AMEGIC::Process_Base *) { return;   }
    bool SetUpIntegrator()           { return 0; } // this might be changed
 
    int                    Size()              { return brokers.size(); }
    AMEGIC::Process_Base * operator[](int idx) { return brokers[idx];   }

    void   SetScale(double _scale);
    double KFactor(double);

    // interfaces to EXTRAXS::XS_Base

    void Add(EXTRAXS::XS_Base *); 

    void AddPoint(const double); 
    bool CalculateTotalXSec();
    void SetTotalXS(int); 
    void SelectOne();
    void DeSelect();

    double Differential(AMATOOLS::vec4d *);
//     double Differential(double, double, double); 
    double Differential2();
    double DSigma(AMATOOLS::vec4d *, bool);
//     double DSigma(double, double, double);
    double DSigma2();

    void   SetISR(ISR::ISR_Handler *);
    void   SetBeam(BEAM::Beam_Handler *);
    bool   SetUpIntegrator(ISR::ISR_Handler *,BEAM::Beam_Handler *);
    bool   OneEvent();
    bool   SameEvent();
    double WeightedEvent();

  };

  inline EXTRAXS::XS_Base * Broker_Group::GetXS() { return xsecs; }

}

#endif

/*!
  This class is meant to be an interface between AMEGIC::Process_Base and EXTRAXS::XS_Base.
  It stores multiple PHASIC::Single_Brokers within 'brokers'.

  At the moment the construction of Broker_Group is explicitly managed by EXTRAXS::XS_Group.
  This could be wise since XS_Group decides by usage of XS_Selector which process is 
  to be considered.
  Thus the extention to a class storing multiple Broker_Groups should possibly be done 
  by defining 'PHASIC::All_Brokers' managed by XS_Group. 

  All routines are self - explanatory or quoted as taken over directly from Process_Base. 
*/
