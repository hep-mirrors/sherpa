#include "BeamRemnant_Initialization.H"

BeamRemnant_Initialization::BeamRemnant_Initialization(std::string _file,std::string _path,
						       MODEL::Model_Base * _model,
						       ISR::ISR_Handler * _isr,
						       BEAM::Beam_Spectra_Handler * _beam) :
  m_dir(_dir), m_file(_file), p_isr(_isr), p_beam(_beam)
{
}

BeamRemnant_Initialization::~BeamRemnant_Initialization();

void BeamRemnant_Initialization::ReadInFile();

bool BeamRemnant_Initialization::InitializeTheBeamRemnantHandler();

void BeamRemnant_Initialization::Output();
