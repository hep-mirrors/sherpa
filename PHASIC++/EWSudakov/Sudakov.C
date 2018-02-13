#include "PHASIC++/EWSudakov/Sudakov.H"
#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <cassert>

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

Sudakov::Sudakov(Process_Base& proc):
  m_proc{ proc },
  p_ampl{ Sudakov::CreateAmplitude(m_proc) },
  m_ci{ m_proc, p_ampl },
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  m_sw{ sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) },
  // TODO: set default to false
  m_check{ Default_Reader().Get<bool>("CHECK_EWSUDAKOV", true) }
{
}

Cluster_Amplitude* Sudakov::CreateAmplitude(Process_Base& proc)
{
  Cluster_Amplitude* ampl{ Cluster_Amplitude::New() };
  ampl->SetNIn(proc.NIn());
  ampl->SetOrderQCD(proc.MaxOrder(0));
  for (size_t i(1);i<proc.MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+proc.MaxOrder(i));
  for(int i(0);i<proc.NIn()+proc.NOut();++i)
    if (i<proc.NIn()) ampl->CreateLeg(Vec4D(),proc.Flavours()[i].Bar());
    else ampl->CreateLeg(Vec4D(),proc.Flavours()[i]);
  ampl->SetProc(&proc);
  ampl->SetProcs(proc.AllProcs());
  return ampl;
}

Cluster_Amplitude* Sudakov::CreateSU2RotatedAmplitude(size_t legindex) const
{
  auto* ampl = p_ampl->Copy();
  auto* leg = ampl->Leg(legindex);
  auto flav = leg->Flav();
  Flavour newflav;
  // TODO: generalise to other flavours (preferredly, add an WeakIsoPartner
  // function to the Flavour class
  if (flav.IsPhoton())
    newflav = Flavour{kf_Z};
  else if (flav.Kfcode() == kf_Z)
    newflav = Flavour{kf_photon};
  leg->SetFlav(newflav);
  return ampl;
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  UpdateAmplitude(mom);
  m_SU2rotatedspinampls.clear();  // they will be calculated on demand
  m_ci.FillSpinAmplitudes(m_spinampls, p_ampl);
  CalculateSpinAmplitudeCoeffs();
  if (m_check) {
    if (CheckCoeffs()) {
      msg_Debugging() << METHOD << "(): Everything's fine ..." << std::endl;
    } else {
      THROW(fatal_error, "Something wrong with the coeff...");
    }
  }
  // TODO: dress with logs, calculate factor for squared amplitude
  return 1.0;
}

void Sudakov::UpdateAmplitude(const ATOOLS::Vec4D_Vector& mom)
{
  p_ampl->SetProcs(m_proc.AllProcs());
  // TODO: this weirdly makes a difference... check
  for(int i(0); i < m_proc.NIn()+m_proc.NOut();++i)
    p_ampl->Leg(i)->SetMom(mom[i]);//i<m_proc.NIn()?-mom[i]:mom[i]);
}

void Sudakov::CalculateSpinAmplitudeCoeffs()
{
  const auto& ampls = m_spinampls[0];
  const auto spinamplnum = ampls.size();
  m_coeffs = std::vector<Complex>(spinamplnum, 0.0);
  for (size_t i{ 0 }; i < spinamplnum; ++i) {
    const auto value = ampls.Get(i);
    if (value == 0.0)
      continue;
    m_coeffs[i] = DoubleLogCoeff(m_spinampls[0], i);
    // TODO: add other coefficients
  }
  for (const auto& coeff : m_coeffs)
    DEBUG_VAR(coeff);
}

Complex Sudakov::DoubleLogCoeff(const Spin_Amplitudes& ampls, size_t spinidx)
{
  const auto spincombination = m_spinampls[0].GetSpinCombination(spinidx);
  if (spincombination.size() != p_ampl->Legs().size())
    THROW(fatal_error, "Inconsistent state");
  Complex coeff{ 0.0 };
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ p_ampl->Leg(i)->Flav() };
    if (flav.IsBoson() && flav.Charge() == 0) {
      // only transverly polarized 
      // TODO: maybe we should only check if the current leg has longitudinal
      // polarisation, not if any has?
      const auto it = std::find(
          spincombination.cbegin(), spincombination.cend(), 2);
      if(it == spincombination.cend()){
	// mixing between neutral gauge bosons: non-diagonal terms appear
	const auto from = flav.Kfcode();
	coeff -= NondiagonalCew(from, from) / 2.0;
	const auto to = (from == kf_photon) ? kf_Z : kf_photon;
	const auto prefactor = -NondiagonalCew(from, to) / 2.0;
	const auto amplratio = 0.0;
	auto it = m_SU2rotatedspinampls.find(i);
	if (it == m_SU2rotatedspinampls.end()) {
	  auto* ampl = CreateSU2RotatedAmplitude(i);
	  m_ci.FillSpinAmplitudes(m_SU2rotatedspinampls[i], ampl);
	  it = m_SU2rotatedspinampls.find(i);
	  ampl->Delete();
	}
	const auto rotated = it->second[0].Get(spincombination);
        // TODO: check if we can just use the spinidx to get the unrotated
        // result
	const auto unrotated = ampls.Get(spincombination);
        // unrotated is guaranteed to be non-zero by virtue of caller logic in
        // CalculateSpinAmplitudeCoeffs
        assert(unrotated != 0.0);
	coeff += prefactor * rotated / unrotated;
      }
    } else {
      // only diagonal terms appear
      coeff -= DiagonalCew(flav, spincombination[i]) / 2.0;
    }
  }
  return coeff;
}

double Sudakov::DiagonalCew(const Flavour& flav, int pol) const
{
  // pol is either chirality or polarisation:
  // 0: + (right-handed or transverse polarisation)
  // 1: - (left-handed or transverse polarisation)
  // 2: 0 (longitudinal polarisation)
  static auto CewLefthandedLepton = (1 + 2*m_cw2) / (4*m_sw2*m_cw2);
  if (flav.IsLepton()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return 1/m_cw2;
    } else {
      return CewLefthandedLepton;
    }
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 1) {
      return (m_sw2 + 27*m_cw2) / (36*m_sw2*m_cw2);
    } else {
      if (flav.IsUptype())
        return 4 / (9*m_cw2);
      else
        return 1 / (9*m_cw2);
    }
  } else if (flav.IsScalar()) {  // cf. eq. (B.18) and (B.16)
    return CewLefthandedLepton;
  } else if (flav.Kfcode() == kf_Wplus) {
    // NOTE: for longitudinal bosons, use the Goldstone equivalence theorem
    if (pol == 2)
      return CewLefthandedLepton;
    else
      return 2/m_sw2;
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double Sudakov::NondiagonalCew(kf_code from, kf_code to) const
{
  if ((from != kf_Z && from != kf_photon) || (to != kf_Z && to != kf_photon))
    THROW(fatal_error, "Only neutral gauge bosons are supported");
  if (from != to)
    return -2.0 * m_cw/m_sw;
  if (from == kf_photon)
    return 2.0;
  if (from == kf_Z)
    return 2.0 * m_cw2/m_sw2;
  THROW(fatal_error, "Logic error");
}

bool Sudakov::CheckCoeffs()
{
  /*
    This simply grabs the process' name and check if it
    finds the right coefficients in.
   */
  bool res(false); size_t count(0);

  const auto denners_coeff = ReferenceCoeffs();
  for(const auto cc : m_coeffs){
    for(const auto dc : denners_coeff){
      if(std::abs(cc.real()-dc) < 1.e-1){
	msg_Debugging() << om::red 
			<< "  Calculated coeff: " << cc.real() << "\t vs \t  Paper's : " << dc
			<< "\n \t Test = " << (std::abs(cc.real()-dc) < 1.e-2)
			<< om::reset << std::endl;
	count += 1;
	break;
      }
    }
  }
  if(count >= denners_coeff.size()) res = true;
  return res;
}

std::vector<double> Sudakov::ReferenceCoeffs()
{
  const auto pname(m_proc.Name());
  std::map<std::string, std::vector<double> > _procs;

  _procs["2_2__e-__e+__mu-__mu+"] = {-2.58,-4.96,-7.35};
  _procs["2_2__e-__e+__u__ub"]    = {-1.86,-4.68,-4.25,-7.07};
  _procs["2_2__e-__e+__d__db"]    = {-1.43,-4.68,-3.82,-7.07};
  _procs["2_2__e-__e+__W+__W-"]   = {-7.35,-4.96,-12.6};
  _procs["2_2__e-__e+__P__P"]     = {-1.29,-8.15};
  _procs["2_2__e-__e+__Z__P"]     = {-1.29,-12.2};
  _procs["2_2__e-__e+__Z__Z"]     = {-1.29,-16.2};

  // check proc name is inside the few we have
  size_t check_name(0);
  for(std::map<std::string, std::vector<double> >::iterator it = _procs.begin();
      it != _procs.end(); ++it){
    if(it->first == pname) check_name = 1;
  }
  if(!check_name) THROW(not_implemented, "No test for proc: " + pname);
 
  return _procs[pname];
}
