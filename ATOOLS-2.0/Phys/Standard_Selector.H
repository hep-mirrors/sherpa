#ifndef Standard_Selector_h
#define Standard_Selector_h

#include "Selector.H"



namespace APHYTOOLS {
  class No_Selector : public Selector_Base {
  public:
    No_Selector()                       { name = std::string ("No_Selector"); sel_log = 0; }
    No_Selector(int,int,Flavour *)      { name = std::string ("No_Selector"); sel_log = 0; }
    bool Trigger(const AMATOOLS::vec4d*)                    { return 1; }
    void BuildCuts(Cut_Data * cuts)                         { }
    void UpdateCuts(double sprime,double y,Cut_Data * cuts) { }
    void Output()                                           { }
  };

  class Energy_Selector : public Selector_Base {
    double * emin, * emax, * value;
  public:
    Energy_Selector(int,int,Flavour *);
    ~Energy_Selector();
    void     SetRange(std::vector<APHYTOOLS::Flavour>,double,double);
    bool     Trigger(const AMATOOLS::vec4d * );
    void     BuildCuts(Cut_Data *);
    void     UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();
  };

  class PT_Selector : public Selector_Base {
    double * ptmin, * ptmax, * value;
  public:
    PT_Selector(int,int,Flavour *);
    ~PT_Selector();
    void     SetRange(std::vector<APHYTOOLS::Flavour>,double,double);
    bool     Trigger(const AMATOOLS::vec4d * );
    void     BuildCuts(Cut_Data *);
    void     UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();
  };

  class Rapidity_Selector : public Selector_Base {
    double * ymin, * ymax, * value;
  public:
    Rapidity_Selector(int,int,Flavour *);
    ~Rapidity_Selector();
    void     SetRange(std::vector<APHYTOOLS::Flavour>,double,double);
    bool     Trigger(const AMATOOLS::vec4d * );
    void     BuildCuts(Cut_Data *);
    void     UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();
  };

  class Angle_Selector : public Selector_Base {
    double ** cosmin, ** cosmax, * value;
  public:
    Angle_Selector(int,int,Flavour *);
    ~Angle_Selector();
    void     SetRange(std::vector<APHYTOOLS::Flavour>,double,double);
    void     SetRange(std::vector<APHYTOOLS::Flavour>,int,double,double);
    bool     Trigger(const AMATOOLS::vec4d * );
    void     BuildCuts(Cut_Data *);
    void     UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();
  };

  class Mass_Selector : public Selector_Base {
    double ** massmin, ** massmax, * value;
  public:
    Mass_Selector(int,int,Flavour *);
    ~Mass_Selector();
    void     SetRange(std::vector<APHYTOOLS::Flavour>,double,double);
    bool     Trigger(const AMATOOLS::vec4d * );
    void     BuildCuts(Cut_Data *);
    void     UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();
  };

  class Summed_PT_Selector : public Selector_Base {
    double             ptmin, ptmax;
    APHYTOOLS::Flavour crit;
  public:
    Summed_PT_Selector(int,int,Flavour *);
    ~Summed_PT_Selector();
    void     SetRange(std::vector<APHYTOOLS::Flavour>,double,double);
    bool     Trigger(const AMATOOLS::vec4d * );
    void     BuildCuts(Cut_Data *);
    void     UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();
  };
}




/*!
  \class No_Selector
  This is a dummy selector that does not apply any cuts on the
  phase space. Consequently, the methods are all empty.
*/

/*!
  \class Energy_Selector
  Selects according to energy ranges for the particles.
  These ranges are set with the method SetRange for each flavour, the 
  default choice is from 0 to the c.m. energy. However, all energy selection 
  takes place within one selector, this results in the following trigger 
  algorithm: The trigger checks via a loop for all particles whether they 
  are inside their energy range. If not, it returns a 0.
*/

/*!
  \class Angle_Selector
  Selects according to angular ranges (cosines) for pairs of particles.
  These ranges are set with the method SetRange for each flavour pair
  or for each flavour w.r.t. a beam specified by its number, the 
  default choice is from -1 to 1. Again, all angular selection 
  takes place within one selector, this results in the following trigger 
  algorithm: The trigger checks via a loop for all pairs of particles whether 
  they are inside their angular range. If not, it returns a 0.
*/


/*!
  \class Mass_Selector
  Selects according to invariant-mass ranges for pairs of particles.
  These ranges are set with the method SetRange for each flavour pair, the default 
  choice is from 0 to the square of the c.m. energy. Again, all mass selection 
  takes place within one selector, this results in the following trigger 
  algorithm: The trigger checks via a loop for all pairs of particles whether 
  they are inside their mass range. If not, it returns a 0.
*/
#endif
