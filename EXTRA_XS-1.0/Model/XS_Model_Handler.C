#include "XS_Model_Handler.H"

#include "XS_Model_SM.H"
#include "XS_Model_Phantom_U1.H"

using namespace EXTRAXS;

XS_Model_Base *XS_Model_Handler::GetModel(const std::string &name)
{
  if (name=="SM") return new XS_Model_SM();
  if (name=="SM+Phantom_U1") return new XS_Model_Phantom_U1();
  return NULL;
}
