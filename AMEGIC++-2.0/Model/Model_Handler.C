#include "Model_Handler.H"
//SM
#include "Model_QCD.H"
#include "Model_EE_QCD.H"
#include "Model_EW.H"
#include "Model_SM.H"
//SUSY
#include "Model_MSSM.H"
#include "Model_THDM.H"
//Large Extra Dimensions
#include "Model_LED.H"

#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

Model_Handler::Model_Handler() {
  AORGTOOLS::msg.Tracking()<<"Initialise Model_Handler"<<endl;
}


Model* Model_Handler::GetModel()
{ 
  Model* m = 0;

  switch (rpa.me.Model()) {
  case Model_Type::pure_QCD : m = new Model_QCD;break;
  case Model_Type::QCD      : m = new Model_EE_QCD;break;
  case Model_Type::EW       : m = new Model_EW;break;
  case Model_Type::SM       : m = new Model_SM;break;
  case Model_Type::MSSM     : m = new Model_MSSM;break;
  case Model_Type::THDM     : m = new Model_THDM;break;
  case Model_Type::LED      : m = new Model_LED;break;
  default:
    msg.Error()<<"Model does not exist !!!"<<endl;
    abort();  
  }

  switch (rpa.me.Model()) {
  case Model_Type::MSSM :
  case Model_Type::THDM : m->Init_Isajet();
  }

  return m;
}


