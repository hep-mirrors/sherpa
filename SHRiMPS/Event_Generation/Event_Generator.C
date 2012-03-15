#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace SHRIMPS;

Event_Generator::Event_Generator(const run_mode::code & runmode,
				 const weight_mode::code & weightmode) :
  m_runmode(runmode), m_weightmode(weightmode), 
  p_cross(NULL),
  p_elastic(NULL), p_sdiff(NULL), p_ddiff(NULL), 
  p_qelastic(NULL), p_inelastic(NULL), 
  p_active(NULL),
  m_minkt2(MBpars("min_kt2"))
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
  //msg_Out()<<METHOD<<": "<<blobs->size()<<".\n";
  if (blobs->size()==1) {
    (*blobs)[0]->AddData("Weight",new ATOOLS::Blob_Data<double>(m_xsec));
    //msg_Out()<<METHOD<<": put xsec = "<<m_xsec<<" in |"<<(*blobs)[0]<<"|\n";
  }
  switch (m_runmode) {
  case run_mode::elastic_events:
    p_active = p_elastic;
    return p_elastic->ElasticEvent(blobs,m_xsec);
  case run_mode::single_diffractive_events:
    p_active = p_sdiff;
    return p_sdiff->SingleDiffractiveEvent(blobs,m_xsec);
  case run_mode::double_diffractive_events:
    p_active = p_ddiff;
    return p_ddiff->DoubleDiffractiveEvent(blobs,m_xsec);
  case run_mode::quasi_elastic_events:
    p_active = p_qelastic;
    return p_qelastic->QuasiElasticEvent(blobs,m_xsec);
  case run_mode::inelastic_events:
    p_active = p_inelastic;
    return p_inelastic->InelasticEvent(blobs,m_xsec,false,
				       m_weightmode==weight_mode::weighted);
  case run_mode::all_min_bias:
    switch (p_cross->SelectCollisionMode()) {
    case 0:
      // elastic collision
      p_active = p_elastic;
      return p_elastic->ElasticEvent(blobs,m_xsec);
    case 1:
      // single diffractive collision
      p_active = p_sdiff;
      return p_sdiff->SingleDiffractiveEvent(blobs,m_xsec);
    case 2:
      // double diffractive collision
      p_active = p_ddiff;
      return p_ddiff->DoubleDiffractiveEvent(blobs,m_xsec);
    case 10:
      // inelastic collision
      p_active = p_inelastic;
      return p_inelastic->InelasticEvent(blobs,m_xsec,false,
					 m_weightmode==weight_mode::weighted);
      break;
    case -1:
    default:
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Did not generate a meaningful collision mode."<<std::endl
		 <<"   Return 'false' and hope for the best."<<std::endl;
    }
    return -1;
  default:
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Event mode "<<m_runmode<<" not initialised yet."
		 <<std::endl
		 <<"   Return 'false' and hope for the best."<<std::endl;
  }
  return -1;
}
