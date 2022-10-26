#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class Pilot_Scale_Setter: public Scale_Setter_Base {

  public:
    Pilot_Scale_Setter(const Scale_Setter_Arguments &args);
    ~Pilot_Scale_Setter();

    void PreCalc(const ATOOLS::Vec4D_Vector &p, const size_t &mode) override;

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
                     const size_t &mode) override;

    double CalculateScale(const ATOOLS::Vec4D_Vector &p,const size_t mode) override;

    PDF::CParam CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const override;

    void SetFixedScale(const std::vector<double> &s) override;

    void SetCaller(Process_Base *const proc) override;

    double Scale(const ATOOLS::stp::id type,const int mode) const override;

    const std::vector<double> &Scales() const override;

    const std::vector<double> &FixedScales() const override;

    ATOOLS::ClusterAmplitude_Vector &Amplitudes() override;

    const ATOOLS::Vec4D_Vector &Momenta() const override;

  private:
    Scale_Setter_Base* ActiveScaleSetter() const;
    Scale_Setter_Base* p_var_scale_setter;
    Scale_Setter_Base* p_mets_scale_setter;

  };

}

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Pilot_Scale_Setter,"PILOT",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,Pilot_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Pilot_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    Pilot_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"pilot scale scheme";
}

Pilot_Scale_Setter::Pilot_Scale_Setter
(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args)
{
  // Construct VAR scale setter from string argument passed to the PILOT scale
  // setter.
  Scale_Setter_Arguments varargs = args;
  varargs.m_scale = StringReplace(args.m_scale, "PILOT", "VAR");
  p_var_scale_setter = Scale_Setter_Base::Scale_Getter_Function::GetObject(
      varargs.m_scale, varargs);
  if (p_var_scale_setter == NULL)
    THROW(fatal_error,
          "Could not construct VAR scale setter within PILOT one.");

  // Construct STRICT_METS scale setter.
  Scale_Setter_Arguments metsargs = args;
  metsargs.m_scale = "STRICT_METS";
  p_mets_scale_setter = Scale_Setter_Base::Scale_Getter_Function::GetObject(
      metsargs.m_scale, metsargs);
  if (p_mets_scale_setter == NULL)
    THROW(fatal_error,
          "Could not construct STRICT_METS scale setter within PILOT one.");

  m_scale.reserve(std::max(p_var_scale_setter->m_scale.size(),
                           p_mets_scale_setter->m_scale.size()));
  SetCouplings();
}

void Pilot_Scale_Setter::PreCalc
(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  ActiveScaleSetter()->PreCalc(momenta, mode);
}

double Pilot_Scale_Setter::Calculate
(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  double ret = ActiveScaleSetter()->Calculate(momenta, mode);
  m_scale = ActiveScaleSetter()->m_scale;
  p_cpls->Calculate();
  return ret;
}

double Pilot_Scale_Setter::CalculateScale(const ATOOLS::Vec4D_Vector &p,const size_t mode)
{
  m_scale.resize(ActiveScaleSetter()->m_scale.size());
  if (p_subs) {
    for (size_t i(0);i<p_subs->size();++i) {
      NLO_subevt *sub((*p_subs)[i]);
      sub->m_mu2.resize(m_scale.size());
    }
  }
  double ret = ActiveScaleSetter()->CalculateScale(p, mode);
  m_scale = ActiveScaleSetter()->m_scale;
  p_cpls->Calculate();
  return ret;
}

PDF::CParam Pilot_Scale_Setter::CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const
{
  return ActiveScaleSetter()->CoreScale(ampl);
}

void Pilot_Scale_Setter::SetFixedScale(const std::vector<double> &s)
{
  ActiveScaleSetter()->SetFixedScale(s);
}

void Pilot_Scale_Setter::SetCaller(Process_Base *const proc)
{
  p_var_scale_setter->SetCaller(proc);
  p_mets_scale_setter->SetCaller(proc);
}

double Pilot_Scale_Setter::Scale(const ATOOLS::stp::id type,const int mode) const
{
  return ActiveScaleSetter()->Scale(type, mode);
}

const std::vector<double> &Pilot_Scale_Setter::Scales() const
{
  return ActiveScaleSetter()->Scales();
}

const std::vector<double> &Pilot_Scale_Setter::FixedScales() const
{
  return ActiveScaleSetter()->FixedScales();
}

ATOOLS::ClusterAmplitude_Vector &Pilot_Scale_Setter::Amplitudes()
{
  return ActiveScaleSetter()->Amplitudes();
}

const ATOOLS::Vec4D_Vector &Pilot_Scale_Setter::Momenta() const
{
  return ActiveScaleSetter()->Momenta();
}

Scale_Setter_Base* Pilot_Scale_Setter::ActiveScaleSetter() const
{
  // Use the cheaper scale setter within the pilot run, switch to the better and
  // more costly one for normal runs.
  if (rpa->gen.PilotRun())
    return p_var_scale_setter;
  else
    return p_mets_scale_setter;
}

Pilot_Scale_Setter::~Pilot_Scale_Setter()
{
  delete p_var_scale_setter;
  delete p_mets_scale_setter;
}
