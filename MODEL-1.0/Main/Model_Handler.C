#include "Model_Handler.H"
#include "Standard_Model.H"
#include "MSSM.H"
//#include "ADD.H"
#include "Message.H"


using namespace MODEL;
using namespace AORGTOOLS;
using namespace APHYTOOLS;

Model_Base * Model_Handler::GetModel(Data_Read * _dataread,string _path) {
  string model     = _dataread->GetValue<string>("MODEL",string("SM"));
  string modelfile = _dataread->GetValue<string>("CONST_FILE",string("Model.dat"));
  
  if (model==string("MSSM")) {
    msg.Debugging()<<"Initialize MSSM through "<<_path<<modelfile<<endl;
    return new MSSM(_path,modelfile);
  }
  /*
  if (model==string("ADD")) {
    msg.Debugging()<<"Initialize ADD through "<<_path<<modelfile<<endl;
    return new ADD(_path,modelfile);
  }
  */
  if (model!=string("SM")) { 
    msg.Error()<<"Error in Model_Handler::GetModel :"<<endl
	       <<"   Tried to initialize model : "<<model<<endl
	       <<"   Option not available. Initialize Standard Model instead."<<endl;
  }
  msg.Debugging()<<"Initialize SM through "<<_path<<modelfile<<endl;
  return new Standard_Model(_path,modelfile);
}
