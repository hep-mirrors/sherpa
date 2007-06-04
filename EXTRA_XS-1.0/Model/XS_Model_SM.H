#ifndef XS_Model_SM_H
#define XS_Model_SM_H

#include "XS_Model_Base.H"

namespace EXTRAXS {

  class XS_Model_SM: public XS_Model_Base {
  public:

    // constructor
    XS_Model_SM();
    
    // member functions
    void Initialize(MODEL::Model_Base *const model,
		    const std::string &file);

    bool IncludesModel(const std::string &name) const;

    std::vector<std::string> IncludedModels() const;

  };// end of class XS_Model_SM

}// end of namespace EXTRAXS

#endif
