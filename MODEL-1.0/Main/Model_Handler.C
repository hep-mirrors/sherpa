#include "Model_Handler.H"
#include "Standard_Model.H"
#include "MSSM.H"
#include "ADD.H"
#include "Message.H"


using namespace MODEL;
using namespace AORGTOOLS;
using namespace APHYTOOLS;

Model_Base * Model_Handler::GetModel(Data_Read * _dataread,std::string _path) {
  std::string model     = _dataread->GetValue("MODEL",std::string("SM"));
  std::string modelfile = _dataread->GetValue("CONST_FILE",std::string("Model.dat"));
  
  if (model==std::string("MSSM")) {
    msg.Debugging()<<"Initialize MSSM through "<<_path<<modelfile<<std::endl;
    return new MSSM(_path,modelfile);
  }
  if (model==std::string("ADD")) {
    msg.Debugging()<<"Initialize ADD through "<<_path<<modelfile<<std::endl;
    return new ADD(_path,modelfile);
  }
  if (model!=std::string("SM")) { 
    msg.Error()<<"Error in Model_Handler::GetModel :"<<std::endl
	       <<"   Tried to initialize model : "<<model<<std::endl
	       <<"   Option not available. Initialize Standard Model instead."<<std::endl;
  }
  msg.Debugging()<<"Initialize SM through "<<_path<<modelfile<<std::endl;
  return new Standard_Model(_path,modelfile);
}
