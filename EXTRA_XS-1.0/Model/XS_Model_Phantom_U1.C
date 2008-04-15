#include "XS_Model_Phantom_U1.H"

#include "Run_Parameter.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;

XS_Model_Phantom_U1::XS_Model_Phantom_U1(): XS_Model_Base("SM+Phantom_U1") {}

void XS_Model_Phantom_U1::Initialize(MODEL::Model_Base *const model,
			     const std::string &file)
{
  XS_Model_Base::Initialize(model,file);
  double scale(rpa.gen.CplScale());
  m_consts["g_1"]=sqrt(4.*M_PI*ScalarFunction("alpha_QED",scale));
  m_consts["\\sin\\theta_W"]=sqrt(ScalarConstant("sin2_thetaW"));
  m_consts["\\cos\\theta_W"]=sqrt(1.0-ScalarConstant("sin2_thetaW"));
  m_consts["P_L"]=1.0;
  m_consts["P_L"]=1.0;
  m_consts["v_{EW}"]=ScalarConstant("vev");
  
  m_consts["g_3"]=sqrt(4.*M_PI*ScalarFunction("alpha_S",scale));  
}

bool XS_Model_Phantom_U1::IncludesModel(const std::string &name) const
{
  if (name=="QED") return true;
  if (name=="EW") return true;
  if (name=="QCD") return true;
  if (name=="SM+PHANTOM_U1") return true;
  return false;
}

std::vector<std::string> XS_Model_Phantom_U1::IncludedModels() const
{
  std::vector<std::string> models(3);
  models[0]="QCD";
  models[1]="EW";
  models[2]="SM+PHANTOM_U1";
  return models;
}
