#include "AddOns/GoSam/GoSam_Variation_Generator.H"

#include "AddOns/GoSam/GoSam_Interface.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Weights.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <ostream>

using namespace GoSam;
using namespace PHASIC;
using namespace ATOOLS;

using Base = Hard_Process_Variation_Generator_Base;
using Args = Hard_Process_Variation_Generator_Arguments;

Variation_Generator::Variation_Generator(const Args& args)
{
  p_proc = dynamic_cast<Single_Process*>(args.p_proc);
  if (!p_proc)
    THROW(fatal_error, "VARIATIONS: GoSam requires a Single_Process.");

  Scoped_Settings s = Settings::GetMainSettings()["GOSAM"];
  Scoped_Settings nominalparams = s["Parameters"];
  Scoped_Settings variations = s["ParameterVariations"];
  for (const std::string& name : variations.GetKeys()) {
    Variation var;
    var.m_name = name;
    Scoped_Settings varparams = variations[name];
    for (const std::string& pname : varparams.GetKeys()) {
      var.m_params[pname] = varparams[pname].SetDefault(0.).Get<double>();
      if (!m_nominal.count(pname)) {
        if (!nominalparams[pname].IsSetExplicitly())
          THROW(fatal_error, "GOSAM: ParameterVariations refers to parameter \""
                +pname+"\", which has no nominal value in GOSAM: Parameters.");
        m_nominal[pname] = nominalparams[pname].SetDefault(0.).Get<double>();
      }
    }
    m_variations.push_back(var);
  }
  if (m_variations.empty())
    THROW(fatal_error, "VARIATIONS: GoSam requires GOSAM: ParameterVariations to be set.");
}

void Variation_Generator::GenerateAndFillWeightsMap(Weights_Map& wgtmap)
{
  double nominal = p_proc->LastXS();
  bool lookup = p_proc->LookUp();
  p_proc->SetLookUp(false);
  ME_Weight_Info *mewgt = p_proc->GetMEwgtinfo();
  const double nominalB = mewgt->m_B, nominalK = mewgt->m_K;
  for (const Variation& var : m_variations) {
    if (nominal == 0.0) {
      wgtmap["GoSam"][var.m_name] = 0.0;
      continue;
    }
    for (const auto& param : var.m_params)
      GoSam_Interface::SetParameter(param.first, param.second);
    double varied = p_proc->Partonic(p_proc->Integrator()->Momenta());
    wgtmap["GoSam"][var.m_name] = varied/nominal;
    for (const auto& param : var.m_params)
      GoSam_Interface::SetParameter(param.first, m_nominal.at(param.first));
  }
  p_proc->SetLookUp(lookup);
  mewgt->m_B = nominalB;
  mewgt->m_K = nominalK;
}

void Variation_Generator::ResetWeightsMap(Weights_Map& wgtmap)
{
  for (const Variation& var : m_variations)
    wgtmap["GoSam"][var.m_name] = 1.0;
}

DECLARE_GETTER(Variation_Generator, "GoSam", Base, Args);

Base* ATOOLS::Getter<Base, Args, Variation_Generator>::
operator()(const Args& args) const
{
  return new Variation_Generator(args);
}

void ATOOLS::Getter<Base, Args, Variation_Generator>::
PrintInfo(std::ostream& str, const size_t width) const
{
  str << "GoSam parameter variations: reweight GoSam-evaluated matrix elements\n"
      << "for alternative values of GOSAM: Parameters, e.g. EFT couplings.\n";
}
