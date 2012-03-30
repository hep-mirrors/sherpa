#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace SHRIMPS;

Event_Generator::Event_Generator(const run_mode::code & runmode,
				 const weight_mode::code & weightmode) :
  m_runmode(runmode), m_thisevent(m_runmode), m_weightmode(weightmode), 
  p_cross(NULL),
  p_elastic(NULL), p_sdiff(NULL), p_ddiff(NULL), 
  p_qelastic(NULL), p_inelastic(NULL), 
  p_active(NULL),
  m_minkt2(MBpars("min_kt2")),
  m_done(false)
{ }

Event_Generator::~Event_Generator() 
{   
  if (p_inelastic) delete p_inelastic; p_inelastic=NULL;
  if (p_qelastic)  delete p_qelastic; p_qelastic=NULL;
  if (p_sdiff)     delete p_sdiff; p_sdiff=NULL;
  if (p_ddiff)     delete p_ddiff; p_ddiff=NULL;
  if (p_elastic)   delete p_elastic; p_elastic=NULL;
}

void Event_Generator::
Initialise(Cross_Sections * cross,Beam_Remnant_Handler * beams,
	   const int & test) {
  p_cross = cross;

  switch (m_runmode) {
  case run_mode::xsecs_only:
    break;
  case run_mode::elastic_events:
    p_elastic   = new Elastic_Event_Generator(p_cross->GetSigmaElastic(),beams);
    m_xsec      = p_elastic->XSec();
    break;
  case run_mode::single_diffractive_events:
    p_sdiff     = new Single_Diffractive_Event_Generator(p_cross->GetSigmaSD(),beams);
    m_xsec      = p_sdiff->XSec();
    break;
  case run_mode::double_diffractive_events:
    p_ddiff     = new Double_Diffractive_Event_Generator(p_cross->GetSigmaDD(),beams);
    m_xsec      = p_ddiff->XSec();
    break;
  case run_mode::quasi_elastic_events:
    p_elastic   = new Elastic_Event_Generator(p_cross->GetSigmaElastic(),beams);
    p_sdiff     = new Single_Diffractive_Event_Generator(p_cross->GetSigmaSD(),beams);
    p_ddiff     = new Double_Diffractive_Event_Generator(p_cross->GetSigmaDD(),beams);
    p_qelastic  = new Quasi_Elastic_Event_Generator(p_elastic,p_sdiff,p_ddiff);
    m_xsec      = p_qelastic->XSec();
    break;
  case run_mode::inelastic_events:
  case run_mode::underlying_event:
    p_inelastic = new Inelastic_Event_Generator(p_cross->GetSigmaInelastic(),
						p_cross->GetEikonals(),beams,
						test);
    m_xsec      = p_inelastic->XSec();
    break;
  case run_mode::all_min_bias:
    p_elastic   = new Elastic_Event_Generator(p_cross->GetSigmaElastic(),beams);
    p_sdiff     = new Single_Diffractive_Event_Generator(p_cross->GetSigmaSD(),beams);
    p_ddiff     = new Double_Diffractive_Event_Generator(p_cross->GetSigmaDD(),beams);
    p_inelastic = new Inelastic_Event_Generator(p_cross->GetSigmaInelastic(),
						p_cross->GetEikonals(),beams,
						test);
    m_xsec = p_cross->SigmaTot();
    break;
  default:
    break;
  }
} 

bool Event_Generator::DressShowerBlob(ATOOLS::Blob * blob) {
  if (m_runmode!=run_mode::underlying_event) {
    msg_Error()<<"Error in "<<METHOD<<" for run mode = "<<m_runmode<<".\n";
    return false;
  }
  msg_Out()<<METHOD<<" for run mode = "<<m_runmode<<".\n";
  return p_inelastic->DressShowerBlob(blob);
}

int Event_Generator::MinimumBiasEvent(ATOOLS::Blob_List * blobs) {
  if (m_done) return 0;
  if (m_thisevent==run_mode::all_min_bias) 
    m_thisevent=p_cross->SelectCollisionMode();
  if (blobs->size()==1) {
    (*blobs)[0]->AddData("Weight",new ATOOLS::Blob_Data<double>(m_xsec));
  }
  switch (m_thisevent) {
  case run_mode::elastic_events:
    m_done   = true;
    p_active = p_elastic;
    return p_elastic->ElasticEvent(blobs,m_xsec);
  case run_mode::single_diffractive_events:
    m_done   = true;
    p_active = p_sdiff;
    return p_sdiff->SingleDiffractiveEvent(blobs,m_xsec);
  case run_mode::double_diffractive_events:
    m_done   = true;
    p_active = p_ddiff;
    return p_ddiff->DoubleDiffractiveEvent(blobs,m_xsec);
  case run_mode::quasi_elastic_events:
    m_done   = true;
    p_active = p_qelastic;
    return p_qelastic->QuasiElasticEvent(blobs,m_xsec);
  case run_mode::inelastic_events:
    p_active = p_inelastic;
    return p_inelastic->InelasticEvent(blobs,m_xsec,false,
				       m_weightmode==weight_mode::weighted);
  case run_mode::all_min_bias:
  default:
    msg_Error()<<"Error in "<<METHOD<<" (event omde = "<<m_thisevent<<"):\n"
	       <<"   No meaningful mode for this event selected.\n"
	       <<"   Will exit.\n";
    exit(1);
  }
  return -1;
}
