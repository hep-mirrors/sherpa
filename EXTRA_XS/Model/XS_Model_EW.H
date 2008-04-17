#ifndef XS_Model_EW_H
#define XS_Model_EW_H

#include "XS_Model_Base.H"

namespace EXTRAXS {

  class XS_Model_EW: public XS_Model_Base {
  public:

    // constructor
    XS_Model_EW();
    
    // member functions
    void Initialize(MODEL::Model_Base *const model,
		    const std::string &file);

    bool IncludesModel(const std::string &name) const;

    std::vector<std::string> IncludedModels() const;

  };// end of class XS_Model_EW

}// end of namespace EXTRAXS

#endif
