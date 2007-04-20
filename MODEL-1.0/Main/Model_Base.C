#include "Model_Base.H"
#include "Spectrum_Generator_Base.H"
#include "Message.H"


using namespace MODEL;
using namespace ATOOLS;


Model_Base::Model_Base(std::string _dir,std::string _file) :
  m_dir(_dir), m_file(_file), p_dataread(NULL),
  p_numbers(NULL), p_constants(NULL), p_functions(NULL), p_matrices(NULL),
  p_spectrumgenerator(NULL)
{ }

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
  if (p_constants!=NULL) delete p_constants;
  if (p_matrices!=NULL) delete p_matrices;
  if (p_dataread!=NULL) delete p_dataread;
  if (p_spectrumgenerator!=NULL) delete p_spectrumgenerator;
}

int Model_Base::ScalarNumber(const std::string _name) {
  if (p_numbers->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No numbers stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0;
  }
  if (p_numbers->count(_name)>0) return (*p_numbers)[_name];

  msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0;
}


double Model_Base::ScalarConstant(const std::string _name) {
  if (p_constants->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<std::endl
	       <<"   No constants stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_constants->count(_name)>0) return (*p_constants)[_name];

  msg.Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


Function_Base * Model_Base::GetScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return NULL;
  }
  if (p_functions->count(_name)>0) return (*p_functions)[_name];

  msg.Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return NULL;
}


double Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_functions->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])(_t);

  msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


double Model_Base::ScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])();

  msg.Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


CMatrix Model_Base::ComplexMatrix(const std::string _name) {
  if (p_matrices->empty()) {
    msg.Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<std::endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<std::endl;
    return CMatrix(1);
  }
  if (p_matrices->count(_name)>0) return (*p_matrices)[_name];

  msg.Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return CMatrix(1);
}


Complex Model_Base::ComplexMatrixElement(const std::string _name,const int _i,const int _j) {
  if (p_matrices->empty()) {
    msg.Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<")("<<_i<<","<<_j<<") : "<<std::endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0;
  }
  if (p_matrices->count(_name)>0) {
    int rank = (*p_matrices)[_name].Rank();
    if (_i<rank && _j<rank && 0<=_i && 0<=_j) return (*p_matrices)[_name][_i][_j];
  }

  msg.Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<")("<<_i<<","<<_j<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0;
}

