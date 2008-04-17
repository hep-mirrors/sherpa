#ifndef XS_Model_Handler_H
#define XS_Model_Handler_H

#include "XS_Model_Base.H"

namespace EXTRAXS {

  class XS_Model_Handler {
  public:

    static XS_Model_Base *GetModel(const std::string &name);

  };// end of class XS_Model_Handler

}// end of namespace EXTRAXS

#endif

