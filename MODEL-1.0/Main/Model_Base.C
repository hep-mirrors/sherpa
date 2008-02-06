#include "Model_Base.H"
#include "Spectrum_Generator_Base.H"
#include "Interaction_Model_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Exception.H"


using namespace MODEL;
using namespace ATOOLS;
using namespace std;


Model_Base::Model_Base(std::string _dir,std::string _file) :
  m_dir(_dir), m_file(_file), p_dataread(NULL),
  p_numbers(NULL), p_constants(NULL), p_functions(NULL), p_matrices(NULL),
  p_spectrumgenerator(NULL), 
  p_vertex(NULL), p_vertextable(NULL), 
  p_decays(NULL)
{
}

Model_Base::~Model_Base() 
{
  if (p_numbers!=NULL) delete p_numbers;
  if (p_functions!=NULL) {
    while (!p_functions->empty()) {
      delete p_functions->begin()->second;
      p_functions->erase(p_functions->begin());
    }
    delete p_functions;
  }
  if (p_constants!=NULL)         delete p_constants;
  if (p_matrices!=NULL)          delete p_matrices;
  if (p_dataread!=NULL)          delete p_dataread;
  if (p_spectrumgenerator!=NULL) delete p_spectrumgenerator;
  if (p_vertex!=NULL)            delete p_vertex;
  if (p_vertextable!=NULL)       delete p_vertextable;
  if (p_decays!=NULL)            delete p_decays;
}

void Model_Base::InitializeInteractionModel()
{
  Data_Read read(m_dir+rpa.gen.Variable("ME_DATA_FILE"));
  string modeltype   = read.GetValue<string>("SIGNAL_MODEL",string("SM"));
  string cplscheme   = read.GetValue<string>("COUPLING_SCHEME",string("Running"));
  string massscheme  = read.GetValue<string>("YUKAWA_MASSES",string("Running"));
  string widthscheme = read.GetValue<string>("WIDTH_SCHEME",string("Fixed"));

  Interaction_Model_Handler mh(this);
  Interaction_Model_Base * model = mh.GetModel(modeltype,cplscheme,massscheme);

  p_vertex        = new Vertex(model);
  p_vertextable   = new Vertex_Table;
  for (int i=0;i<p_vertex->MaxNumber();++i) {
    if ((*p_vertex)[i]->on) {
      (*p_vertextable)[(*p_vertex)[i]->in[0]].push_back((*p_vertex)[i]);
    }
  }

  delete model;
}

void Model_Base::FillDecayTables() {
  p_decays = new All_Decays(this);
  Flavour flav;
  for (std::map<ATOOLS::Flavour, Vertex_List>::iterator vit=p_vertextable->begin();
       vit!=p_vertextable->end();vit++) {
    flav = vit->first;
    if (!flav.IsStable() && flav.Width()<0.) {
      std::cout<<METHOD<<" : "<<flav<<" : "<<flav.Width()<<std::endl;
      p_decays->AddToDecays(flav); 
    }
  }
  p_decays->InitializeDecayTables();
  p_decays->CalculateWidths();
}

int Model_Base::ScalarNumber(const std::string _name) {
  if (p_numbers->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No numbers stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0;
  }
  if (p_numbers->count(_name)>0) return (*p_numbers)[_name];

  msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0;
}


double Model_Base::ScalarConstant(const std::string _name) {
  if (p_constants->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<std::endl
	       <<"   No constants stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_constants->count(_name)>0) return (*p_constants)[_name];

  msg_Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


Function_Base * Model_Base::GetScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return NULL;
  }
  if (p_functions->count(_name)>0) return (*p_functions)[_name];

  msg_Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return NULL;
}


double Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_functions->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])(_t);

  msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


double Model_Base::ScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])();

  msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


CMatrix Model_Base::ComplexMatrix(const std::string _name) {
  if (p_matrices->empty()) {
    msg_Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<std::endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<std::endl;
    return CMatrix(1);
  }
  if (p_matrices->count(_name)>0) return (*p_matrices)[_name];

  msg_Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return CMatrix(1);
}


Complex Model_Base::ComplexMatrixElement(const std::string _name,const int _i,const int _j) {
  if (p_matrices->empty()) {
    msg_Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<")("<<_i<<","<<_j<<") : "<<std::endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0;
  }
  if (p_matrices->count(_name)>0) {
    int rank = (*p_matrices)[_name].Rank();
    if (_i<rank && _j<rank && 0<=_i && 0<=_j) return (*p_matrices)[_name][_i][_j];
  }

  msg_Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<")("<<_i<<","<<_j<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0;
}

