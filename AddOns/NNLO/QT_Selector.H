#ifndef NNLO__src__QT_Selector_H
#define NNLO__src__QT_Selector_H

#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/Process_Base.H"

namespace SHNNLO {

  class QT_Selector : public PHASIC::Selector_Base {
    PHASIC::Process_Base *p_proc;
    double m_qtmin, m_cqtmin;
    int m_type;
  public:
    QT_Selector(const PHASIC::Selector_Key &key);
    bool Trigger(ATOOLS::Selector_List &);
    void BuildCuts(PHASIC::Cut_Data *) {}
    inline double QTCut() const { return m_cqtmin; }
  };

}// end of namespace SHNNLO

#endif
