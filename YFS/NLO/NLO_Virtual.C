#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "YFS/NLO/NLO_Base.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include <functional>
#include <map>
#include <cstdlib>
#include <iostream>
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/NLO/Photon_Counterterm.H"
#include "MODEL/Main/Model_Base.H"
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include "YFS/NLO/NLO_Base_Internal.H"

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

double NLO_Base::CalculateVirtual() {
  m_lhelok = false;
  m_vborn_ok = false;
  if (CeexSuppliesVirtual())
    return (m_ceexvirt - 1.) * m_born;
  if (m_eex_virt) {
    // subtract born to avoid double counting
    // already present in eex!!
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  if (!m_looptool)
    return 0;
  double virt;
  double sub;
  p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  CheckMassReg();
  /*
    YFS: VIRTUAL_CANONICAL_FRAME (default false): evaluate the Born's one-loop
    in the rest frame of the beams with the beams exactly on the z axis
    (CanonicalBeamFrame, which documents why). The event's Born point after
    ISR has two beams with the same tiny p_T, and there OpenLoops' 2 -> 2
    virtual is off by up to 8% of itself (Z-pole mu mu, full EW, alpha(0):
    0.085 against 0.091 in every rotated frame). Not for CEEX_Virtual:
    helicity, whose helicity labels are the lab's.
  */
  static const bool vcanon(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["VIRTUAL_CANONICAL_FRAME"].SetDefault(false).Get<bool>());
  if (vcanon && !(m_useceex && m_ceexvirtsrc == ceexvirt::helicity))
    virt = p_virt->CalcInFrame(m_plab, CanonicalBeamFrame(m_plab), m_born);
  else
    virt = p_virt->Calc(m_plab, m_born);
  if (m_check_virt_born) {
    // the provider's Born is pointlike, m_born is dressed with the pion form
    // factor, so compare against the dressed provider Born
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()
                         * ExternalFormFactor(m_plab, m_flavs), 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "\nSherpa" << ":  " << m_born << std::endl
                  << "PhaseSpace Point = ";
      for (auto _p : m_plab)
        msg_Error() << _p << std::endl;
    }
  }
  if (p_virt->FailCut())
    return 0;
  if (m_virt_sub && p_virt->p_loop_me->Mode() != 1)
    sub = p_dipoles->CalculateVirtualSub();
  else
    sub = 0;
  m_virt_raw = virt; m_virt_subval = sub * m_born / m_rescale_alpha;
  m_oneloop = (virt - sub * m_born / m_rescale_alpha);
  // YFS: CEEX_Virtual: helicity - the same loop call, resolved by helicity
  if (m_useceex && m_ceexvirtsrc == ceexvirt::helicity &&
      p_virt->p_loop_me->Mode() == 0 && !IsZero(virt) && m_born != 0.)
    BuildLoopHelicityFactors(virt, sub);
  if (IsZero(virt)){
    m_zeroV++;
    return 0;
  }
  if (p_virt->p_loop_me->Mode() == 1)
    m_oneloop /= m_rescale_alpha;
  if (IsBad(m_oneloop) || IsBad(sub)) {
    msg_Error() << "YFS Virtual is NaN" << std::endl
                << "Virtual:  " << virt << std::endl
                << "Subtraction: " << sub * m_born << std::endl
                << "PhaseSpace Point: " << std::endl
                << m_plab << std::endl;
  }
  if (m_check_poles == 1) {
    if (!m_virt_sub)
      sub = p_dipoles->CalculateVirtualSub();
    double p1 = p_virt->p_loop_me->ME_E1() * p_virt->m_factor;
    double yfspole = p_dipoles->Get_E1();
    int ncorrect = ::countMatchingDigits(p1, -yfspole);
    double reldiff = (p1 + yfspole) / p1;
    if (!IsEqual(p1, -yfspole, 1e-4)) {
      msg_Error() << "Poles do not cancel in YFS Virtuals" << std::endl
                  << "Correct digits =  " << ncorrect << std::endl
                  << "Relative diff =  " << reldiff
                  << std::endl
                  // <<"Process =  "<<p_virt->p_loop_me->Name()<<std::endl
                  << "One-Loop Provider V eps^{-1}  = " << p1 << std::endl
                  << "Sherpa V eps^{-1} = " << yfspole << std::endl
                  << "Sherpa/One-Loop = " << yfspole / p1 << std::endl;
      return 0;
    } else {
      int i = 0;
      msg_Debugging() << std::setprecision(32);
      msg_Debugging() << "Poles cancel in YFS Virtuals to " << ncorrect
                      << " digits" << std::endl
                      << "Relative diff =  " << reldiff << std::endl;
      m_histograms1d["SinglePoleCD"]->Insert(ncorrect);
      m_histograms1d["OneLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
      m_histograms1d["OneLoopEpsLP"]->Insert(log10(fabs(p1)));
      m_histograms1d["relativediff"]->Insert(log10(fabs(reldiff)));
      msg_Debugging() << std::setprecision(32)
                      << "One-Loop Provider V eps^{-1}  = " << p1 << std::endl
                      << "Sherpa V eps^{-1}  = " << yfspole << std::endl;
    }
  }
  // v of this event for the RV_MODE 1 remainder v_{n+1} - v
  if (m_born != 0. && !IsBad(m_oneloop)) {
    m_vborn = m_oneloop/m_born;
    m_vborn_ok = true;
  }
  return m_oneloop;
}

/*
  RV_MODE 1: the Born's v with the subtraction on the Born point's OWN legs,
  v_B = V_fin(m_plab)/B - B~(m_plab legs)/kappa, from the loop CalculateVirtual
  already evaluated. The event's v (m_vborn) subtracts B~
  on the legs the YFS form factor is built on: the II dipole and the initial
  legs of the IF dipoles at the FULL beams (YFS_Handler::MakeYFS), while the
  loop is at m_plab, whose beams are reduced by the ISR photons. For a
  one-photon event the (n+1)-body point's legs are those event legs, and
  v_{n+1} - v would be exact given that convention; but with a hard ISR
  companion the point of a soft photon sits at the reduced beams, and
  v_{n+1} - v -> B~(reduced legs) - B~(full-beam legs) != 0 as the photon
  goes soft: 1e-4 .. 6e-3 per soft photon, i.e. RV_SOFT_CUT dependence.
  v_{n+1} - v_B -> 0 for every photon (both IR subtracted on their own legs:
  the object the standalone referee olref.py computes), and differs from the
  exact one-photon identity only by rho (v - v_B), the Born virtual's own
  leg convention.
*/
bool NLO_Base::BornVirtualOnOwnLegs(double &vb) {
  if (!m_looptool || m_born == 0.) return false;
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  const double sub(p_nlodipoles->CalculateVirtualSub());
  vb = m_virt_raw/m_born - sub/m_rescale_alpha;
  return !IsBad(vb);
}

double NLO_Base::CalculateVV() {
  if (!m_vvtool)
    return 0;
  if (m_eex_virt) {
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  double virt;
  double sub;
  // CheckMassReg();
  if (!HasISR())
    virt = p_vv->Calc(m_bornMomenta, m_born);
  else
    virt = p_vv->Calc(m_plab, m_born);
  if (m_check_virt_born) {
    // the provider's Born is pointlike, m_born is dressed with the pion form
    // factor, so compare against the dressed provider Born
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()
                         * ExternalFormFactor(m_plab, m_flavs), 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "\nSherpa" << ":  " << m_born << std::endl
                  << "PhaseSpace Point = ";
      for (auto _p : m_plab)
        msg_Error() << _p << std::endl;
    }
  }
  if (p_vv->FailCut())
    return 0;
  if (m_virt_sub && p_virt->p_loop_me->Mode() != 1)
    sub = p_dipoles->CalculateVirtualSub();
  else
    sub = 0;
  double sub2 = p_dipoles->CalculateVVSubEps();
  // m_oneloop = (virt- sub * m_born/m_rescale_alpha );
  m_oneloop = (virt - sub * CalculateVirtual() / m_rescale_alpha -
               0.5 * sub * sub * m_born / m_rescale_alpha);
  if (p_virt->p_loop_me->Mode() == 1) {
    m_oneloop /= m_rescale_alpha;
  }
  if (IsBad(m_oneloop) || IsBad(sub)) {
    msg_Error() << "YFS Virtual is NaN" << std::endl
                << "Virtual:  " << virt << std::endl
                << "Subtraction: " << sub * m_born << std::endl
                << "PhaseSpace Point: " << std::endl
                << m_plab << std::endl;
  }
  double loope1 =
      p_vv->p_loop_me->ME_E1() *
      p_vv->m_factor; //*p_vv->m_factor;//+p_virt->p_loop_me->ME_E1()*p_virt->m_factor;;
  double loope2 =
      2. * p_vv->p_loop_me->ME_E2() * p_vv->m_factor * p_vv->m_factor;
  double yfse1 = p_dipoles->Get_E1();
  double yfse2 = p_dipoles->GetVV_E2();
  if (m_check_poles == 1) {
    if (!m_virt_sub)
      sub = p_dipoles->CalculateVirtualSub();
    const double p1 = p_vv->p_loop_me->ME_E1() * p_vv->m_factor;
    const double p2 =
        2. * p_vv->p_loop_me->ME_E2() * p_vv->m_factor * p_vv->m_factor;
    const double yfspole1 = (p_dipoles->Get_E1());
    const double yfspole2 = p_dipoles->GetVV_E2();
    PRINT_VAR(p1 / yfspole1);
    int ncorrect1 = ::countMatchingDigits(p1, yfspole1, 32);
    int ncorrect2 = ::countMatchingDigits(p2, -yfspole2, 32);
    if (!IsEqual(p2, -yfspole2, 1e-6) || ncorrect1 < 10) {
      msg_Error() << "Poles do not cancel in YFS Double Virtuals" << std::endl
                  << "Correct digits \epsion^{-1} =  " << ncorrect1 << std::endl
                  << "Correct digits \epsion^{-2} =  " << ncorrect2
                  << std::endl;
      return 0;
    } else {
      int i = 0;
      msg_Debugging() << std::setprecision(32);
      msg_Out() << "Poles cancel in YFS double Virtuals to " << ncorrect2
                << " digits" << std::endl;
      m_histograms1d["SinglePoleVV"]->Insert(ncorrect1);
      m_histograms1d["DoublePoleVV"]->Insert(ncorrect2);
    }
  }
  return 0;
}
