#ifndef Simple_Chain_inl_H
#define Simple_Chain_inl_H

#include "Simple_Chain.H"

namespace AMISIC {
	
  inline void Simple_Chain::SetEnvironmentFile(std::string _m_environmentfile)
  { m_environmentfile=_m_environmentfile; }

  inline std::string Simple_Chain::EnvironmentFile()
  { return m_environmentfile; }

  inline void Simple_Chain::SetNorm(GridResultType _m_norm)
  { m_norm=_m_norm; }

  inline Simple_Chain::GridResultType Simple_Chain::Norm()
  { return m_norm; }

  inline Profile_Function_Base *Simple_Chain::ProfileFunction()
  { return p_profile; }

} // end of namespace AMISIC

#endif
