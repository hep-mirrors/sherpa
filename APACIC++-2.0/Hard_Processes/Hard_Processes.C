#include "Hard_Processes.H"

// global variables
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Run_Parameter.H"

#include "QED_Processes.H"
#include "QCD_Processes.H"

using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace APACIC;
using namespace AMEGIC;
using namespace EXTRAXS;

Hard_Processes::Hard_Processes(APHYTOOLS::Selector_Data * _seldata,
			       ISR::ISR_Handler * _isr,BEAM::Beam_Handler * _beam,
			       bool & success) : 
  seldata(_seldata), isr(_isr), beam(_beam)
{
  success = 0;
  broker  = new All_Processes(AMEGIC::XS_MODE);
  if (ProcessesInit()) success = 1;
}

Hard_Processes::~Hard_Processes() 
{
  if (broker) delete broker;
}

bool Hard_Processes::ProcessesInit() 
{
  if ( ((rpa.gen.Beam1() == Flavour(kf::e)) && 
	(rpa.gen.Beam2() == (Flavour(kf::e).Bar())) ) ||
       ((rpa.gen.Beam1() == (Flavour(kf::e)).Bar()) && 
	(rpa.gen.Beam2() == Flavour(kf::e)) ) ) 
    {
      two2two = new QED_Processes(0);
      two2two->Initialize(isr, beam, seldata, broker);
      return 1;
    }
  if ( ( (rpa.gen.Beam1() == Flavour(kf::p_plus)) && 
	 (rpa.gen.Beam2() == (Flavour(kf::p_plus).Bar())) )   ||
       ( (rpa.gen.Beam1() == (Flavour(kf::p_plus)).Bar()) &&
	 (rpa.gen.Beam2() == Flavour(kf::p_plus)) )           ||
       ( (rpa.gen.Beam1() == (Flavour(kf::p_plus)))       &&
	 (rpa.gen.Beam2() == Flavour(kf::p_plus)) )           ||
       ( (rpa.gen.Beam1() == (Flavour(kf::p_plus).Bar())) &&
	 (rpa.gen.Beam2() == Flavour(kf::p_plus).Bar()) )    ) 
    {
      two2two = new QCD_Processes();
      two2two->Initialize(isr, beam, seldata, broker);
      return 1;
    }
  return 0;
};

bool Hard_Processes::PrepareCalculation() 
{
  AMEGIC::Topology              * top = 0;
  Vec4D                         * moms = 0;
  std::vector<double>             results;
  std::vector<Single_Process *>   links;

  if (broker->InitAllProcesses(top,moms,results,links)) return 1;
  else {
      msg.Error()<<"Error in Hard_Processes::PrepareCalculation() : Broker Initialization failed !"
		 <<std::endl;
      return 0;
  }
};

bool Hard_Processes::CalculateCrossSections() {
  return broker->CalculateTotalXSec();
};

bool Hard_Processes::OneEvent() {
  return broker->OneEvent();
}













