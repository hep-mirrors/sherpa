#ifndef QED_Process_H
#define QED_Process_H

#include "XS_Group.H"

namespace EXTRAXS {

  class QED_Processes: public XS_Group {
    void CreateSelector();
  public:
    QED_Processes();
  };
}

#endif
