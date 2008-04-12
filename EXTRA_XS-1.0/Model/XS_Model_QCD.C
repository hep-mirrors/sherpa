#include "XS_Model_QCD.H"

#include "Run_Parameter.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;

XS_Model_QCD::XS_Model_QCD(): XS_Model_Base("QCD") {}

void XS_Model_QCD::Initialize(MODEL::Model_Base *const model,
			     const std::string &file)
{
  XS_Model_Base::Initialize(model,file);
  m_consts["g_3"]=sqrt(4.*M_PI*ScalarFunction("alpha_S",rpa.gen.CplScale()));  
}

bool XS_Model_QCD::IncludesModel(const std::string &name) const
{
  if (name=="QCD") return true;
  return false;
}

std::vector<std::string> XS_Model_QCD::IncludedModels() const
{
  std::vector<std::string> models(1);
  models[0]="QCD";
  return models;
}
