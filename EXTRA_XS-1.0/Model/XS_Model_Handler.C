#include "XS_Model_Handler.H"

#include "XS_Model_SM.H"

using namespace EXTRAXS;

XS_Model_Base *XS_Model_Handler::GetModel(const std::string &name)
{
  if (name=="SM") return new XS_Model_SM();
  return NULL;
}
