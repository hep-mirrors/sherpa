#include "XS_Model_MSSM.H"

#include "Run_Parameter.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;

XS_Model_MSSM::XS_Model_MSSM(): XS_Model_Base("MSSM") {}

void XS_Model_MSSM::Initialize(MODEL::Model_Base *const model,
			     const std::string &file)
{
  XS_Model_Base::Initialize(model,file);
  double scale(rpa.gen.CplScale());
  m_consts["g_1"]=sqrt(4.*M_PI*ScalarFunction("alpha_QED",scale));
  m_consts["\\sin\\theta_W"]=sqrt(ScalarConstant("sin2_thetaW"));
  m_consts["\\cos\\theta_W"]=sqrt(1.0-ScalarConstant("sin2_thetaW"));
  m_consts["P_L"]=1.0;
  m_consts["P_R"]=1.0;
  m_consts["v_{EW}"]=ScalarConstant("vev");
  m_consts["g_3"]=sqrt(4.*M_PI*ScalarFunction("alpha_S",scale));  
}

bool XS_Model_MSSM::IncludesModel(const std::string &name) const
{
  if (name=="QED") return true;
  if (name=="EW") return true;
  if (name=="QCD") return true;
  if (name=="SM") return true;
  if (name=="MSSM") return true;
  return false;
}

std::vector<std::string> XS_Model_MSSM::IncludedModels() const
{
  std::vector<std::string> models(5);
  models[0]="QCD";
  models[1]="EW";
  models[2]="QED";
  models[3]="SM";
  models[4]="MSSM";
  return models;
}
