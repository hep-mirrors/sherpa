#ifndef Amisic_inl_H
#define Amisic_inl_H

#include "Amisic.H"

namespace AMISIC {

  inline MI_Base *Amisic::HardBase()
  { return p_hardbase; }

  inline MI_Base *Amisic::SoftBase()
  { return p_softbase; }

#ifdef USING__SHERPA
  inline SHERPA::Matrix_Element_Handler *Amisic::HardMEHandler()
  { return p_hardbase->MEHandler(); }

  inline SHERPA::Matrix_Element_Handler *Amisic::SoftMEHandler()
  { return p_softbase->MEHandler(); }
#endif

  inline EXTRAXS::XS_Base *Amisic::HardXS()
  { return p_hardbase->XS(); }

  inline EXTRAXS::XS_Base *Amisic::SoftXS()
  { return p_softbase->XS(); }

  inline const Amisic::ModelID Amisic::HardModel() const
  { return m_hardmodel; }

  inline const Amisic::ModelID Amisic::SoftModel() const
  { return m_softmodel; }

} // end of namespace AMISIC

#endif
