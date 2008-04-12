#include "XS_Model_EW.H"

#include "Run_Parameter.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;

XS_Model_EW::XS_Model_EW(): XS_Model_Base("EW") {}

void XS_Model_EW::Initialize(MODEL::Model_Base *const model,
			     const std::string &file)
{
  XS_Model_Base::Initialize(model,file);
  m_consts["g_1"]=sqrt(4.*M_PI*ScalarFunction("alpha_QED",rpa.gen.CplScale()));
  m_consts["\\sin\\theta_W"]=sqrt(ScalarConstant("sin2_thetaW"));
  m_consts["\\cos\\theta_W"]=sqrt(1.0-ScalarConstant("sin2_thetaW"));
  m_consts["P_L"]=1.0;
  m_consts["P_R"]=1.0;
  m_consts["v_{EW}"]=ScalarConstant("vev");
}

bool XS_Model_EW::IncludesModel(const std::string &name) const
{
  if (name=="QED") return true;
  if (name=="EW") return true;
  return false;
}

std::vector<std::string> XS_Model_EW::IncludedModels() const
{
  std::vector<std::string> models(1);
  models[0]="EW";
  return models;
}
