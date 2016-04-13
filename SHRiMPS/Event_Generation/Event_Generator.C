#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace SHRIMPS;

Event_Generator::Event_Generator() :
  m_runmode(MBpars.RunMode()),m_thisevent(m_runmode),
  p_inelastic(NULL), p_active(NULL), m_xsec(0.)
{ }

Event_Generator::~Event_Generator() 
{   
  if (p_inelastic) delete p_inelastic; p_inelastic=NULL;
}

void Event_Generator::Initialise() {
  p_inelastic = new Inelastic_Event_Generator();
  m_xsec += p_inelastic->XSec();
} 

void Event_Generator::Reset() {
  if (p_active) p_active->Reset();
  m_thisevent = m_runmode;
}


bool Event_Generator::DressShowerBlob(ATOOLS::Blob * blob) {
  if (m_runmode!=run_mode::underlying_event) {
    msg_Error()<<"Error in "<<METHOD<<" for run mode = "<<m_runmode<<".\n";
    return false;
  }
  msg_Out()<<METHOD<<" for run mode = "<<m_runmode<<".\n";
  return false; 
}

int Event_Generator::MinimumBiasEvent(ATOOLS::Blob_List * blobs) {
  msg_Tracking()<<"======================================================\n"
		<<METHOD<<": "<<blobs->size()<<" blobs.\n"
		<<"======================================================\n";
  if (blobs->size()==1) {
    (*blobs)[0]->AddData("Weight",new ATOOLS::Blob_Data<double>(m_xsec));
    (*blobs)[0]->AddData("Weight_Norm",new ATOOLS::Blob_Data<double>(m_xsec));
    (*blobs)[0]->AddData("Trials",new ATOOLS::Blob_Data<double>(1.));
  }
  p_active = p_inelastic;
  return p_inelastic->GenerateEvent(blobs,false);
}

void Event_Generator::Test(const std::string & dirname) {
  msg_Info()<<METHOD<<": Starting.\n";
  p_inelastic->Test(dirname);
}
