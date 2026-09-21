#include "AddOns/Griffin/Griffin_Interface.H"

#include <iostream>

using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace griffin;

Griffin::Griffin_Interface *Griffin::Griffin_Interface::p_instance = nullptr;

Griffin::Griffin_Interface::Griffin_Interface() :
      ME_Generator_Base("Griffin")
{
  RegisterDefaults();
  p_instance = this;
}

Griffin::Griffin_Interface::~Griffin_Interface()
{
  if (p_instance == this) p_instance = nullptr;
}

Griffin::Griffin_Interface &Griffin::Griffin_Interface::Instance()
{
  if (!p_instance)
    THROW(fatal_error, "Griffin_Interface::Instance(): not yet constructed "
                        "-- is \"Griffin\" listed in ME_GENERATORS?");
  return *p_instance;
}

void Griffin::Griffin_Interface::RegisterDefaults()
{
   Scoped_Settings s{ Settings::GetMainSettings()["GRIFFIN"] };
   s["Order"].SetDefault(griffinorder::nnlo);
   s["Delta_Alpha"].SetDefault(0.06);
}

bool Griffin::Griffin_Interface::Initialize(MODEL::Model_Base *const model,
        BEAM::Beam_Spectra_Handler *const beam,
        PDF::ISR_Handler *const isr,
        YFS::YFS_Handler *const yfs)
{
    PrintLogo(msg->Info());
    rpa->gen.AddCitation(
        1, "The Griffin library is described in \\cite{Chen:2022dow}.");

    Settings& s = Settings::GetMainSettings();
    Scoped_Settings ss{ Settings::GetMainSettings()["GRIFFIN"] };

    const double GF    = s["GF"].Get<double>();
    m_order            = ss["Order"].Get<griffinorder::code>();
    const double delap = ss["Delta_Alpha"].Get<double>();

    m_ewscheme = s["EW_SCHEME"].Get<ew_scheme::code>();

    if (m_ewscheme == ew_scheme::alphamZ) {
        griffin::SMval inv;

        inv.set(MW,  Flavour(kf_Wplus).Mass());
        inv.set(MZ,  Flavour(kf_Z).Mass());
        inv.set(MH,  Flavour(kf_h0).Mass());
        inv.set(ME,  Flavour(kf_e).Mass());
        inv.set(MM,  Flavour(kf_mu).Mass());
        inv.set(ML,  Flavour(kf_tau).Mass());
        inv.set(MD,  Flavour(kf_d).Mass());
        inv.set(MS,  Flavour(kf_s).Mass());
        inv.set(MB,  Flavour(kf_b).Mass());
        inv.set(MU,  Flavour(kf_u).Mass());
        inv.set(MC,  Flavour(kf_c).Mass());
        inv.set(MT,  Flavour(kf_t).Mass());

        inv.set(al,  s_model->ScalarConstant("alpha_QED"));
        inv.set(als, s_model->ScalarConstant("alpha_S"));

        inv.set(Delal, delap);
        inv.set(Gmu,   GF);
        inv.set(GamW,  Flavour(kf_Wplus).Width());
        inv.set(GamZ,  Flavour(kf_Z).Width());
        m_griffin.emplace<griffin::SMval>(inv);
        std::visit([](auto &gx) {
            std::cout << "\n";
            std::cout << "complex-pole mass: mw = " << gx.get(MWc) << "\n";
            std::cout << "PDG mass:          mw = " << gx.get(MW)  << "\n";
            std::cout << "complex-pole mass: mZ = " << gx.get(MZc) << "\n";
            std::cout << "PDG mass:          mZ = " << gx.get(MZ)  << "\n";
            std::cout << "alpha(0):          al = " << gx.get(al)  << "\n";
            std::cout << "                   1/al = " << 1.0/gx.get(al) << "\n";
            std::cout << std::endl;
        }, m_griffin);
    }

    if (m_ewscheme == static_cast<ew_scheme::code>(3)) {
        griffin::inval tmp;
        tmp.set(MZ,   Flavour(kf_Z).Mass());
        tmp.set(MW,   Flavour(kf_Wplus).Mass());
        tmp.set(als,  s_model->ScalarConstant("alpha_S"));
        tmp.set(GamZ, Flavour(kf_Z).Width());
        tmp.set(GamW, Flavour(kf_Wplus).Width());
        tmp.set(MH,   Flavour(kf_h0).Mass());
        tmp.set(MT,   Flavour(kf_t).Mass());
        tmp.set(MB,   Flavour(kf_b).Mass());     // MSbar at scale mu=MZ
        tmp.set(Delal, delap);
        tmp.set(Gmu,   GF);

        m_griffin.emplace<griffin::SMvalGMwMz>(tmp);
        std::visit([](auto &gx) {
            std::cout << "\n";
            std::cout << "complex-pole mass: mw = " << gx.get(MWc) << "\n";
            std::cout << "PDG mass:          mw = " << gx.get(MW)  << "\n";
            std::cout << "complex-pole mass: mZ = " << gx.get(MZc) << "\n";
            std::cout << "PDG mass:          mZ = " << gx.get(MZ)  << "\n";
            std::cout << "alpha(0):          al = " << gx.get(al)  << "\n";
            std::cout << "                   1/al = " << 1.0/gx.get(al) << "\n";
            std::cout << std::endl;
        }, m_griffin);
    }

    return true;
}

Griffin::GriffinProcess
Griffin::Griffin_Interface::MakeProcess(const PHASIC::Process_Info &pi)
{
  // Griffin is 2->2 only, so the incoming/outgoing pair is
  // (flavour, antiflavour) and one PDG code per side is enough.
  GriffinProcess proc;
  proc.initial_pdg = static_cast<int>(pi.m_ii.GetExternal()[0]);
  proc.final_pdg   = static_cast<int>(pi.m_fi.GetExternal()[0]);
  if (proc.initial_pdg == proc.final_pdg) proc.norm = 2.;  // Bhabha
  return proc;
}

Griffin::GriffinProcess
Griffin::Griffin_Interface::MakeProcess(const External_ME_Args &args)
{
  GriffinProcess proc;
  proc.initial_pdg = static_cast<int>(args.m_inflavs[0]);
  proc.final_pdg   = static_cast<int>(args.m_outflavs[0]);
  if (proc.initial_pdg == proc.final_pdg) proc.norm = 2.;  // Bhabha
  return proc;
}

Griffin::Griffin_Interface::Kinematics
Griffin::Griffin_Interface::KinVars(const Vec4D_Vector &momenta)
{
  const double s = (momenta[0]+momenta[1]).Abs2();
  const double t = (momenta[0]-momenta[2]).Abs2();
  return {s, 1. + 2.*t/s};
}

template <class FA, class SW, class Mat, class SM>
double Griffin::Griffin_Interface::DiffXSecTerm(int initial_pdg, int final_pdg,
                                                 double s, double cost,
                                                 const SM &gx)
{
  FA FAi(initial_pdg, gx), FAf(final_pdg, gx);
  SW SWi(initial_pdg, gx), SWf(final_pdg, gx);
  Mat M(initial_pdg, final_pdg, VEC, VEC,
        FAi.result().real(), FAf.result().real(),
        SWi.result().real(), SWf.result().real(),
        s, cost, gx);
  M.setkinvar(s, cost);

  Cplx resvv, resav, resva, resaa;
  M.setform(VEC, VEC);   resvv = M.result();
  M.setform(AXV, VEC);   resav = M.result();
  M.setform(VEC, AXV);   resva = M.result();
  M.setform(AXV, AXV);   resaa = M.result();

  return real((1. + cost*cost) * (resvv*conj(resvv) + resav*conj(resav)
                                 + resva*conj(resva) + resaa*conj(resaa))
              + 4.*cost * (resvv*conj(resaa) + resva*conj(resav)).real());
}

void Griffin::Griffin_Interface::EvaluateLoop(const Vec4D_Vector &momenta,
                                               const GriffinProcess &proc,
                                               METOOLS::DivArrD &virt) const
{
  if (momenta.size() != 4)
    THROW(fatal_error, "Griffin library is for 2->2 scattering only");
  const double born = EvaluateLO(momenta, proc);
  double higher = 0.;
  switch (m_order) {
    case griffinorder::nlo:  higher = EvaluateNLO (momenta, proc); break;
    case griffinorder::nloe: higher = EvaluateNLOE(momenta, proc); break;
    case griffinorder::nnlo: higher = EvaluateNNLO(momenta, proc); break;
    default: THROW(not_implemented, "GRIFFIN: Order " + ToString(m_order));
  }
  virt.Finite() = higher - born;
}

double Griffin::Griffin_Interface::EvaluateLO(const Vec4D_Vector &momenta,
                                               const GriffinProcess &proc) const
{
  const Kinematics kv = KinVars(momenta);
  if (kv.cost > 1. || kv.cost < -1.)
    msg_Error() << "CosTheta out of range in " << METHOD << std::endl;
  return std::visit([&](auto &gx) {
    return DiffXSecTerm<FA_SMLO, SW_SMLO, matel>(
        proc.initial_pdg, proc.final_pdg, kv.s, kv.cost, gx);
  }, m_griffin);
}

double Griffin::Griffin_Interface::EvaluateNLO(const Vec4D_Vector &momenta,
                                                const GriffinProcess &proc) const
{
  const Kinematics kv = KinVars(momenta);
  return std::visit([&](auto &gx) {
    return DiffXSecTerm<FA_SMNLO, SW_SMNLO, mat_SMNLO>(
        proc.initial_pdg, proc.final_pdg, kv.s, kv.cost, gx);
  }, m_griffin);
}

double Griffin::Griffin_Interface::EvaluateNNLO(const Vec4D_Vector &momenta,
                                                 const GriffinProcess &proc) const
{
  const Kinematics kv = KinVars(momenta);
  return std::visit([&](auto &gx) {
    return DiffXSecTerm<FA_SMNNLO, SW_SMNNLO, mat_SMNNLO>(
        proc.initial_pdg, proc.final_pdg, kv.s, kv.cost, gx);
  }, m_griffin);
}

double Griffin::Griffin_Interface::EvaluateNLOE(const Vec4D_Vector &momenta,
                                                 const GriffinProcess &proc) const
{
  // The exact-NLO "improved Born" combination genuinely differs from the
  // shared LO/NLO/NNLO helicity sum above (it interferes the LO and NLOE
  // amplitudes rather than squaring one amplitude), so it keeps its own
  // formula rather than going through DiffXSecTerm.
  const Kinematics kv = KinVars(momenta);
  const double sw = s_model->ComplexConstant("csin2_thetaW").real();
  return std::visit([&](auto &gx) {
      FA_SMLO FAi0(proc.initial_pdg, gx), FAf0(proc.final_pdg, gx);
      SW_SMLO SWi0(proc.initial_pdg, gx), SWf0(proc.final_pdg, gx);
      matel M0(proc.initial_pdg, proc.final_pdg, VEC, VEC,
               FAi0.result().real(), FAf0.result().real(),
               sw, sw, kv.s, kv.cost, gx);

      FA_SMNLO FAi(proc.initial_pdg, gx), FAf(proc.final_pdg, gx);
      SW_SMNLO SWi(proc.initial_pdg, gx), SWf(proc.final_pdg, gx);
      mat_SMeNLO M(proc.initial_pdg, proc.final_pdg, VEC, VEC,
                    FAi.result().real(), FAf.result().real(),
                    sw, sw, kv.s, kv.cost, gx);

      M0.setkinvar(kv.s, kv.cost);
      M.setkinvar(kv.s, kv.cost);

      Cplx res0vv, res0va, res0av, res0aa;
      Cplx resvv, resva, resav, resaa;

      M0.setform(VEC, VEC);   res0vv = M0.result();
      M0.setform(AXV, VEC);   res0av = M0.result();
      M0.setform(VEC, AXV);   res0va = M0.result();
      M0.setform(AXV, AXV);   res0aa = M0.result();

      M.setform(VEC, VEC);    resvv = M.result();
      M.setform(AXV, VEC);    resav = M.result();
      M.setform(VEC, AXV);    resva = M.result();
      M.setform(AXV, AXV);    resaa = M.result();

      return real((1. + kv.cost*kv.cost) *
          ((res0vv + 2.*(resvv - res0vv)) * conj(res0vv)
           + (res0av + 2.*(resav - res0av)) * conj(res0av)
           + (res0va + 2.*(resva - res0va)) * conj(res0va)
           + (res0aa + 2.*(resaa - res0aa)) * conj(res0aa))
          + 4.*kv.cost*(res0vv*conj(res0aa) + (resvv - res0vv)*conj(res0aa)
                    + res0vv*conj(resaa - res0aa)
                    + res0va*conj(res0av) + (resva - res0va)*conj(res0av)
                    + res0va*conj(resav - res0av))
      );
    }, m_griffin);
}

void Griffin::Griffin_Interface::EvaluateBorn(const Vec4D_Vector &momenta,
                                               const GriffinProcess &proc,
                                               double &born) const
{
  const Kinematics kv = KinVars(momenta);
  born = std::visit([&](auto &gx) {
        FA_SMLO FAi(proc.initial_pdg, gx), FAf(proc.final_pdg, gx);
        SW_SMLO SWi(proc.initial_pdg, gx), SWf(proc.final_pdg, gx);

        matel M(proc.initial_pdg, proc.final_pdg, VEC, VEC, FAi, FAf, SWi, SWf, kv.s, kv.cost, gx);
        M.setkinvar(kv.s, kv.cost);

        Cplx resvv, resva, resav, resaa;

        M.setform(VEC, VEC);   resvv = M.result();
        M.setform(AXV, VEC);   resav = M.result();
        M.setform(VEC, AXV);   resva = M.result();
        M.setform(AXV, AXV);   resaa = M.result();

        return real((1. + kv.cost*kv.cost) * (resvv*conj(resvv) + resav*conj(resav)
                                       + resva*conj(resva) + resaa*conj(resaa))
                    + 4.*kv.cost * (resvv*conj(resaa) + resva*conj(resav)).real()
                    - 2.*(1. + kv.cost*kv.cost) * (resvv*conj(resav) + resva*conj(resaa)).real()
                    - 4.*kv.cost * (resvv*conj(resva) + resav*conj(resaa)));
    }, m_griffin);
}

std::istream &Griffin::operator>>(std::istream &str, griffinorder::code &mode)
{
  std::string tag;
  str>>tag;
  mode=griffinorder::nnlo;
  if      (tag == "NNLO")    mode = griffinorder::nnlo;
  else if (tag == "NLOE")    mode = griffinorder::nloe;
  else if (tag == "NLO")     mode = griffinorder::nlo;
  else if (tag == "LO")      mode = griffinorder::lo;
  else THROW(fatal_error, "Unknown GRIFFIN: Order = " + tag);
  return str;
}

std::ostream &Griffin::operator<<(std::ostream &str,const griffinorder::code &ym)
{
  if      (ym==griffinorder::lo)     return str<<"LO";
  else if (ym==griffinorder::nlo)    return str<<"NLO";
  else if (ym==griffinorder::nnlo)   return str<<"NNLO";
  else if (ym==griffinorder::nloe)   return str<<"Exact NLO";
  return str<<"unknown";
}

void Griffin::Griffin_Interface::PrintLogo(std::ostream &s){
  s<<"======================================================"<<std::endl;
  s<<"                                                      "<<std::endl;
  s<<"======================================================"<<std::endl;

  s<<"     ______ ____   ____ ______ ______ ____ _   __"<<std::endl;
  s<<"    / ____// __ \\ /  _// ____// ____//  _// | / /"<<std::endl;
  s<<"   / / __ / /_/ / / / / /_   / /_    / / /  |/ / "<<std::endl;
  s<<"  / /_/ // _, _/_/ / / __/  / __/  _/ / / /|  /  "<<std::endl;
  s<<"  \\____//_/ |_|/___//_/    /_/    /___//_/ |_/   "<<std::endl;

  s<<"======================================================"<<std::endl;
  s<<"                 version 1.1                          "<<std::endl;
  s<<"           Lisong Chen and Ayres Freitas             "<<std::endl;
  s<<"          https://arxiv.org/abs/2211.16272            "<<std::endl;
  s<<"======================================================"<<std::endl;
}


DECLARE_GETTER(Griffin::Griffin_Interface,"Griffin",PHASIC::ME_Generator_Base,PHASIC::ME_Generator_Key);

PHASIC::ME_Generator_Base *ATOOLS::Getter<PHASIC::ME_Generator_Base,PHASIC::ME_Generator_Key,
                                  Griffin::Griffin_Interface>::
operator()(const PHASIC::ME_Generator_Key &key) const
{
  return new Griffin::Griffin_Interface();
}

void ATOOLS::Getter<PHASIC::ME_Generator_Base,PHASIC::ME_Generator_Key,Griffin::Griffin_Interface>::
PrintInfo(std::ostream &str,const std::size_t width) const
{
  str<<"Interface to the Griffin loop ME generator";
}
