#include "Run_Parameter.H"
#include "Model_Handler.H"
#include "Standard_Model.H"
#include "MSSM.H"
#include "ADD.H"
#include "Message.H"


using namespace MODEL;
using namespace ATOOLS;

Model_Base * Model_Handler::GetModel(Data_Read * _dataread,std::string _path,std::string _file) {
  std::string model     = _dataread->GetValue("MODEL",std::string("SM"));
  
  Model_Base * modelbase = 0;
  if (model==std::string("MSSM")) {
    modelbase = new MSSM(_path,_file);
  }
  if (model==std::string("ADD")) {
    modelbase = new ADD(_path,_file);
  }
  if (model!=std::string("SM") && modelbase==0) { 
    msg.Error()<<"Error in Model_Handler::GetModel :"<<std::endl
	       <<"   Tried to initialize model : "<<model<<std::endl
	       <<"   Option not available. Initialize Standard Model instead."<<std::endl;
  }
  modelbase = new Standard_Model(_path,_file);
  rpa.gen.SetModel(modelbase);
  return modelbase;
}
