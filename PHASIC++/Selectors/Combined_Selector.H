#ifndef PHASIC_Selectors_Combined_Selector_H
#define PHASIC_Selectors_Combined_Selector_H

#include "PHASIC++/Selectors/Selector.H"

namespace PHASIC {

  class Combined_Selector: public Selector_Base {
  protected: 

    long int m_count;
    int      m_update;
    std::vector<Selector_Base *> m_sels;
    std::vector<std::pair<std::string,double> > m_osc;

    Selector_Key FindArguments(const Selector_Key &strings,
			       size_t &starty,size_t &startx);

  public:

    Combined_Selector(Process_Integrator *const proc);

    ~Combined_Selector();

    bool Initialize(const Selector_Key &key);

    bool Trigger(const ATOOLS::Vec4D_Vector &);
    bool JetTrigger(const ATOOLS::Vec4D_Vector &,
                    const ATOOLS::Flavour_Vector &,int);
    bool JetTrigger(const ATOOLS::Vec4D_Vector &p);
    bool NoJetTrigger(const ATOOLS::Vec4D_Vector &);

    void BuildCuts(Cut_Data *);
    void UpdateCuts(double,double,Cut_Data *);
    void AddOnshellCondition(std::string,double);
    int  NeedUpdate();

    void Output();

    Selector_Base *GetSelector(const std::string &name) const;

  };

}

#endif
