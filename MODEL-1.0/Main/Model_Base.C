#include "Model_Base.H"
#include "Spectrum_Generator_Base.H"
#include "Message.H"


using namespace MODEL;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;


Model_Base::Model_Base(std::string _dir,std::string _file) :
  m_dir(_dir), m_file(_file), p_dataread(NULL),
  p_numbers(NULL), p_functions(NULL), p_constants(NULL), p_matrices(NULL),
  p_spectrumgenerator(NULL)
{ 
  //  APHYTOOLS::ParticleInit(m_dir); 
}

Model_Base::~Model_Base() 
{
  if (!p_numbers->empty()) {
    p_numbers->erase(p_numbers->begin(),p_numbers->end());
    p_numbers = NULL;
  }
  if (!p_functions->empty()) {
    p_functions->erase(p_functions->begin(),p_functions->end());
    p_functions = NULL;
  }
  if (!p_constants->empty()) {
    p_constants->erase(p_constants->begin(),p_constants->end());
    p_constants = NULL;
  }
  if (!p_matrices->empty()) {
    p_matrices->erase(p_matrices->begin(),p_matrices->end());
    p_matrices = NULL;
  }
  if (p_dataread) {
    delete p_dataread;
    p_dataread = NULL;
  }
  if (p_spectrumgenerator) {
    delete p_spectrumgenerator;
    p_spectrumgenerator = NULL;
  }
}

int Model_Base::ScalarNumber(const std::string _name) {
  if (p_numbers->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<endl
	       <<"   No numbers stored in model "<<m_name<<". Return 0."<<endl;
    return 0;
  }
  if (p_numbers->count(_name)>0) return (*p_numbers)[_name];

  msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return 0;
}


double Model_Base::ScalarConstant(const std::string _name) {
  if (p_constants->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<endl
	       <<"   No constants stored in model "<<m_name<<". Return 0."<<endl;
    return 0.;
  }
  if (p_constants->count(_name)>0) return (*p_constants)[_name];

  msg.Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return 0.;
}


Function_Base * Model_Base::GetScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<endl;
    return NULL;
  }
  if (p_functions->count(_name)>0) return (*p_functions)[_name];

  msg.Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return NULL;
}


double Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_functions->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])(_t);

  msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return 0.;
}


double Model_Base::ScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])();

  msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return 0.;
}


CMatrix Model_Base::ComplexMatrix(const std::string _name) {
  if (p_matrices->empty()) {
    msg.Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<endl;
    return CMatrix(1);
  }
  if (p_matrices->count(_name)>0) return (*p_matrices)[_name];

  msg.Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return CMatrix(1);
}


Complex Model_Base::ComplexMatrixElement(const std::string _name,const int _i,const int _j) {
  if (p_matrices->empty()) {
    msg.Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<") : "<<endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<endl;
    return 0;
  }
  if (p_matrices->count(_name)>0) {
    int rank = (*p_matrices)[_name].Rank();
    if (_i<rank && _j<rank) return (*p_matrices)[_name][_i][_j];
  }

  msg.Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<") : "<<endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<endl;
  return 0;
}

