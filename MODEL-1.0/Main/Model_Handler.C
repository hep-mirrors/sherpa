#include "Run_Parameter.H"
#include "Model_Handler.H"
#include "Standard_Model.H"
#include "THDM.H"
#include "MSSM.H"
//#include "MUED.H"
#include "Fourth_Generation_Leptons.H"
#include "SM_Phantom_U1.H"
#include "ADD.H"
#include "Message.H"


using namespace MODEL;
using namespace ATOOLS;

Model_Base * Model_Handler::GetModel(Data_Read * _dataread,std::string _path,std::string _file) {
  std::string model     = _dataread->GetValue("MODEL",std::string("SM"));
  
  Model_Base * modelbase = 0;
  if (model==std::string("PHANTOM_U1")) {
    modelbase = new SM_Phantom_U1(_path,_file);
  }
//   else if (model==std::string("MUED")) {
//     modelbase = new MUED(_path,_file);
//   }
  else if (model==std::string("MSSM")) {
    modelbase = new MSSM(_path,_file);
  }
  else if (model==std::string("THDM")) {
    modelbase = new THDM(_path,_file);
  }
  else if (model==std::string("FOURTH_GEN_LEPTONS")) {
    modelbase = new Fourth_Generation_Leptons(_path,_file);
  }
  else if (model==std::string("ADD")) {
    modelbase = new ADD(_path,_file);
  }
  else if (model==std::string("SM")) { 
    modelbase = new Standard_Model(_path,_file);
  }
  else {
    msg_Error()<<METHOD<<"(): Tried to initialize model : "<<model<<std::endl
	       <<"   Option not available. Initialize Standard Model instead."<<std::endl;
    modelbase = new Standard_Model(_path,_file);
  }
  if (!modelbase->RunSpectrumGenerator()) {
    msg_Error()<<METHOD<<"(): RunSpectrumGenerator() delivered false. Abort."<<std::endl;
    abort();
  }

  modelbase->InitializeInteractionModel();
  modelbase->FillDecayTables();
  rpa.gen.SetModel(modelbase);
  return modelbase;
}
