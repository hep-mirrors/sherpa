#ifndef XS_Selector_H
#define XS_Selector_H

#include "Single_XS.H"
#include "XS_Group.H"

namespace EXTRAXS {
  
  class XS_Group;

  class XS_Selector {
  private:

    XS_Base *p_owner;

    bool m_offshell;

  public:

    // constructor
    XS_Selector(XS_Base *const owner);
    
    // member functions
    Single_XS *GetXS(const size_t nin,const size_t nout,
		     const ATOOLS::Flavour *flavours);
    
    static size_t FindInGroup(XS_Group *const group,XS_Base *&base,
			      const size_t nin,const size_t nout,
			      const ATOOLS::Flavour *flavours);
    
    // inline functions
    inline void SetOffShell(const bool offshell)   { m_offshell=offshell;   }

    inline bool OffShell() const  { return m_offshell;  }

  };// end of class XS_Selector

}//end of namespace EXTRAXS

#endif
