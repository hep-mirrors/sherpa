#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Scales/Pilot_Scale_Setter.H"

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Pilot_Scale_Setter, "PILOT", Scale_Setter_Base,
               Scale_Setter_Arguments);

Scale_Setter_Base *
ATOOLS::Getter<Scale_Setter_Base, Scale_Setter_Arguments, Pilot_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const {
  return new Pilot_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base, Scale_Setter_Arguments,
                    Pilot_Scale_Setter>::PrintInfo(std::ostream &str,
                                                   const size_t width) const {
  str << "pilot scale scheme";
}

Scale_Setter_Arguments METS_Args(Scale_Setter_Arguments args)
{
  args.m_scale = "STRICT_METS";
  return args;
}

Scale_Setter_Arguments VAR_Args(Scale_Setter_Arguments args)
{
  args.m_scale = StringReplace(args.m_scale, "PILOT", "VAR");
  return args;
}

Pilot_Scale_Setter::Pilot_Scale_Setter(const Scale_Setter_Arguments &args,
                                       const int mets_mode)
    : Scale_Setter_Base(args), mets{METS_Args(args), mets_mode}, var{VAR_Args(args)} {}

double Pilot_Scale_Setter::CalculateScale(const ATOOLS::Vec4D_Vector &p,
                                          const size_t mode) {
  if (rpa->gen.IsPilotRun())
    return var.Calculate(p, mode);
  else
    return mets.Calculate(p, mode);
}

PDF::CParam
Pilot_Scale_Setter::CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const {
  if (rpa->gen.IsPilotRun())
    return var.CoreScale(ampl);
  else
    return mets.CoreScale(ampl);
}

void Pilot_Scale_Setter::SetCaller(Process_Base *const proc) {
  var.SetCaller(proc);
  mets.SetCaller(proc);
}

double Pilot_Scale_Setter::Scale(const ATOOLS::stp::id type,
                                 const int mode) const {
  if (rpa->gen.IsPilotRun())
    return var.Scale(type, mode);
  else
    return mets.Scale(type, mode);
}

const std::vector<double> &Pilot_Scale_Setter::Scales() const {
  if (rpa->gen.IsPilotRun())
    return var.Scales();
  else
    return mets.Scales();
}

const std::vector<double> &Pilot_Scale_Setter::FixedScales() const {
  if (rpa->gen.IsPilotRun())
    return var.FixedScales();
  else
    return mets.FixedScales();
}

ATOOLS::ClusterAmplitude_Vector &Pilot_Scale_Setter::Amplitudes() {
  if (rpa->gen.IsPilotRun())
    return var.Amplitudes();
  else
    return mets.Amplitudes();
}

const ATOOLS::Vec4D_Vector &Pilot_Scale_Setter::Momenta() const {
  if (rpa->gen.IsPilotRun())
    return var.Momenta();
  else
    return mets.Momenta();
}

void Pilot_Scale_Setter::PreCalc(const ATOOLS::Vec4D_Vector &p,
                                 const size_t &mode) {
  if (rpa->gen.IsPilotRun())
    return var.PreCalc(p, mode);
  else
    return mets.PreCalc(p, mode);
}

double Pilot_Scale_Setter::Calculate(const ATOOLS::Vec4D_Vector &p,
                                     const size_t &mode) {
  if (rpa->gen.IsPilotRun())
    return var.Calculate(p, mode);
  else
    return mets.Calculate(p, mode);
}
