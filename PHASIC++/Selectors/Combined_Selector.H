#ifndef PHASIC_Selectors_Combined_Selector_H
#define PHASIC_Selectors_Combined_Selector_H

#include "PHASIC++/Selectors/Selector.H"

namespace PHASIC {

  class Combined_Selector: public Selector_Base {
  protected: 

    long int m_count;
    int      m_on, m_res, m_rsres;
    std::vector<Selector_Base *> m_sels;
    std::vector<std::pair<size_t,double> > m_osc;

  public:

    Combined_Selector(Process_Base *const proc);

    ~Combined_Selector();

    bool Initialize(const Selector_Key &key);

    bool Trigger(const ATOOLS::Vec4D_Vector &p,
                 const ATOOLS::Flavour *fl=NULL, size_t n=0);
    bool RSTrigger(ATOOLS::NLO_subevtlist *const subs);
    bool Trigger(ATOOLS::Selector_List &sl);

    bool Pass() const;

    void BuildCuts(Cut_Data *);
    void AddOnshellCondition(size_t,double);

    void Output();

    Selector_Base *GetSelector(const std::string &name) const;
    void ListSelectors() const;

    inline void SetOn(const int on) { m_on=on; }

    inline int On() const { return m_on; }

    inline int Result() const { return m_on?m_res:1; }
    inline int RSResult() const { return m_on?m_rsres:1; }

    inline void SetResult(const int res) { m_res=res; }

    /// same as Result, but with results for all parameter variaitons (e.g.
    /// when varying the QCUT scale)
    std::vector<ATOOLS::Event_Weights> CombinedResults() const;

  };

}

#endif
