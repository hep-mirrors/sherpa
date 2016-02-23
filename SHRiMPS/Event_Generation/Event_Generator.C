#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace SHRIMPS;

Event_Generator::Event_Generator() :
  m_runmode(MBpars.RunMode()),m_thisevent(m_runmode),
  p_inelastic(NULL), 
  m_done(false)
{ }

Event_Generator::~Event_Generator() 
{   
  if (p_inelastic) delete p_inelastic; p_inelastic=NULL;
}

void Event_Generator::Initialise(Beam_Remnant_Handler * beams) {
  p_inelastic = new Inelastic_Event_Generator(beams);
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
  msg_Out()<<"======================================================\n"
	   <<"======================================================\n"
	   <<METHOD<<"(done = "<<m_done<<"):\n"<<(*blobs)<<"\n";
  if (m_done) return 0;
  if (blobs->size()==1) {
    (*blobs)[0]->AddData("Weight",new ATOOLS::Blob_Data<double>(m_xsec));
    (*blobs)[0]->AddData("Weight_Norm",new ATOOLS::Blob_Data<double>(m_xsec));
    (*blobs)[0]->AddData("Trials",new ATOOLS::Blob_Data<double>(1.));
  }
  return p_inelastic->GenerateEvent(blobs,false);
}

void Event_Generator::Test(const std::string & dirname) {
  msg_Info()<<METHOD<<": Starting.\n";
  p_inelastic->Test(dirname);
}
