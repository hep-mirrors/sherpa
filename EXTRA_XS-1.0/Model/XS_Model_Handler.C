#include "XS_Model_Handler.H"

#include "XS_Model_MSSM.H"
#include "XS_Model_SM.H"
#include "XS_Model_QCD.H"
#include "XS_Model_EW.H"

using namespace EXTRAXS;

XS_Model_Base *XS_Model_Handler::GetModel(const std::string &name)
{
  if (name=="MSSM") return new XS_Model_MSSM();
  if (name=="pure_QCD") return new XS_Model_QCD();
  if (name=="pure_EW") return new XS_Model_EW();
  return new XS_Model_SM();
}
