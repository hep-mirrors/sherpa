#include "Hard_Processes.H"

// global variables
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Run_Parameter.H"

#include "QED_Processes.H"
#include "QCD_Processes.H"

using namespace APACIC;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace ISR;

using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;


Hard_Processes::Hard_Processes(ISR::ISR_Handler * _isr,
			       BEAM::Beam_Handler * _beam,
			       bool & success) :
  isr(_isr), beam(_beam)
{
  success = 0;
  if (ProcessesInit()) success = 1;
}

Hard_Processes::~Hard_Processes() {}

bool Hard_Processes::ProcessesInit() {
  if ( ( (rpa.gen.Beam1() == Flavour(kf::e)) &&
	 (rpa.gen.Beam2() == (Flavour(kf::e).Bar())) ) ||
       ( (rpa.gen.Beam1() == (Flavour(kf::e)).Bar()) &&
	 (rpa.gen.Beam2() == Flavour(kf::e)) ) ) {
    two2two = (new QED_Processes())->CreateBroker();
    msg.Debugging()<<"In Hard_Processes::Process_Init : "<<std::endl;
    if (two2two) msg.Debugging()<<" Initialised new Broker " 
				<<two2two->Name()<<std::endl;
    else msg.Debugging()<<" Cannot initialise new Broker ! "<<std::endl;
    return 1;
  }
  if ( ( (rpa.gen.Beam1() == Flavour(kf::p_plus))         &&
	 (rpa.gen.Beam2() == (Flavour(kf::p_plus).Bar())) )   ||
       ( (rpa.gen.Beam1() == (Flavour(kf::p_plus)).Bar()) &&
	 (rpa.gen.Beam2() == Flavour(kf::p_plus)) )           ||
       ( (rpa.gen.Beam1() == (Flavour(kf::p_plus)))       &&
	 (rpa.gen.Beam2() == Flavour(kf::p_plus)) )           ||
       ( (rpa.gen.Beam1() == (Flavour(kf::p_plus).Bar())) &&
	 (rpa.gen.Beam2() == Flavour(kf::p_plus).Bar()) )    ) {
    two2two = (new QCD_Processes())->CreateBroker();
    msg.Debugging()<<"In Hard_Processes::Process_Init : "<<std::endl;
    if (two2two) msg.Debugging()<<" Initialised new Broker " 
				<<two2two->Name()<<std::endl;
    else msg.Debugging()<<" Cannot initialise new Broker ! "<<std::endl;
    return 1;
  }
  return 0;
};

bool Hard_Processes::PrepareCalculation() {
  msg.Debugging()<<"Hard_Processes::PrepareCalculation() : "<<std::endl
		 <<" SetISR  : "<<isr<<" for "<<two2two->Name()<<std::endl
		 <<" SetBeam : "<<beam<<" for "<<two2two->Name()<<std::endl;
  return two2two->SetUpIntegrator(isr,beam);
};


bool Hard_Processes::CalculateCrossSections() {
  return two2two->CalculateTotalXSec();
};


bool Hard_Processes::OneEvent() {
  return two2two->OneEvent();
}













