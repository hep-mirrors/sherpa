#include "MODEL/Main/Spectrum_Generator_Base.H"

using namespace MODEL;
using namespace ATOOLS;

Spectrum_Generator_Base::Spectrum_Generator_Base(ATOOLS::Data_Reader * _dataread,Model_Base * _model) :
  p_dataread(_dataread), p_model(_model)  {}

Spectrum_Generator_Base::~Spectrum_Generator_Base() 
{ 
}




