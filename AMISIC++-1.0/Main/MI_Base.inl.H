#ifndef MI_Base_inl_H
#define MI_Base_inl_H

#include "MI_Base.H"

namespace AMISIC {

  inline void MI_Base::SetStart(const double _m_start,unsigned int i)
  { if (i<m_nparameter) m_start[i]=_m_start; }

  inline const double MI_Base::Start(unsigned int i) const
  { return m_start[i]; }

  inline void MI_Base::SetStop(const double _m_stop,unsigned int i)
  { if (i<m_nparameter) m_stop[i]=_m_stop; }

  inline const double MI_Base::Stop(unsigned int i) const
  { return m_stop[i]; }

  inline bool MI_Base::StopGeneration(TypeID type)
  { 
    switch (type) {
    case HardEvent:
      return s_stophard;
    case SoftEvent:
      return s_stopsoft;
    case Unknown:
      return s_stophard&&s_stopsoft;
    default:	
      return false;
    }
    return true;
  }

#ifdef USING__Sherpa
  inline SHERPA::Matrix_Element_Handler *MI_Base::MEHandler()
  { return p_mehandler; }
#endif

  inline EXTRAXS::XS_Base *MI_Base::XS()
  { return p_xs; }

  inline void MI_Base::SetRemnantHandler(ATOOLS::Remnant_Info *const remnanthandler,
					 const size_t i)
  { s_remnanthandlers[i]=remnanthandler; }

  inline ATOOLS::Remnant_Info *MI_Base::RemnantHandler(const size_t i)
  { return s_remnanthandlers[i]; }

} // end of namespace AMISIC

#endif
