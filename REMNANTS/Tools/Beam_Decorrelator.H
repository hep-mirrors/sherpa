#ifndef REMNANTS_Tools_Beam_Decorrelator_H
#define REMNANTS_Tools_Beam_Decorrelator_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Primordial_KPerp.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include <list>

namespace REMNANTS {
  static ATOOLS::Vec4D s_AxisP(1.,0.,0.,1.);
  static ATOOLS::Vec4D s_AxisM(1.,0.,0.,-1.);

  class Remnant_Handler;
  class Beam_Decorrelator {
  private:
    bool               m_on;
    double             m_expo, m_xiP, m_invxiP, m_maxeta, m_mass2, m_deltaM;
    Remnant_Handler  * p_rhandler;
    Primordial_KPerp*  p_kperpGenerator;

    ATOOLS::Blob     * p_softblob;
    ATOOLS::Particle * p_beam, * p_spect;
    std::list<ATOOLS::Particle * > m_softgluons;

    ATOOLS::Vec4D      m_pbeam, m_pspect, m_ktvec;
    ATOOLS::Poincare   m_boost, m_rotat;
    double             m_Q2, m_Q, m_mbeam, m_mbeam2, m_mspect, m_mspect2;
    double             m_minMbeam, m_minMbeam2, m_minMspect, m_minMspect2;
    double             m_x, m_kt2;

    bool MustEmit(ATOOLS::Particle * pi, ATOOLS::Particle * pj);
    bool SoftEmission();
    void InitSoftEmission();
    bool DefineKinematics();
    bool MakeKinematics();
  public:
    Beam_Decorrelator();
    ~Beam_Decorrelator() = default;

    void Initialize(Remnant_Handler * const rhandler);
    bool operator()(ATOOLS::Blob * softblob);
    void Reset();
  };
}

#endif
