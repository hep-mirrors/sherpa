#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Event_Generator::Event_Generator(Cross_Sections * xsecs,const bool & test) :
  m_runmode(MBpars.RunMode()), m_thisevent(m_runmode),
  p_inelastic(NULL), p_elastic(NULL), p_soft_diffractive(NULL),
  p_active(NULL), m_xsec(0.), m_xsec_inel(0.), m_xsec_elas(0.), m_xsec_diff(0.),
  m_eventsuccessful(1), m_mustinit(1)
{
  InitGenerator(xsecs,test);
}

Event_Generator::~Event_Generator() 
{   
  if (p_inelastic) delete p_inelastic; p_inelastic=NULL;
  if (p_elastic) delete p_elastic; p_elastic=NULL;
  if (p_soft_diffractive) delete p_soft_diffractive; p_soft_diffractive=NULL;
}

void Event_Generator::InitGenerator(Cross_Sections * xsecs,const bool & test) {
  switch (m_runmode) {
  case run_mode::inelastic_events:
    p_inelastic = new Inelastic_Event_Generator(xsecs->GetSigmaInelastic(),test);
    break; 
  case run_mode::elastic_events:
    p_elastic = new Elastic_Event_Generator(xsecs->GetSigmaElastic(),test);
    break;
  case run_mode::soft_diffractive_events:
    p_soft_diffractive = new Soft_Diffractive_Event_Generator(xsecs->GetSigmaSD(),test);
    break;
  case run_mode::all_min_bias:
    p_inelastic = new Inelastic_Event_Generator(xsecs->GetSigmaInelastic(),test);
    p_elastic = new Elastic_Event_Generator(xsecs->GetSigmaElastic(),test);
    p_soft_diffractive = new Soft_Diffractive_Event_Generator(xsecs->GetSigmaSD(),test);
    break;
  }
} 

void Event_Generator::
Initialise(Remnant_Handler * remnants,Cluster_Algorithm * cluster) {
  //m_xsec = 0.;
  if (p_inelastic) {
    p_inelastic->Initialise(remnants,cluster);
    m_xsec += p_inelastic->XSec();
    m_xsec_inel += p_inelastic->XSec();
  }
  if (p_elastic) {
    p_elastic->Initialise();
    m_xsec += p_elastic->XSec();
    m_xsec_elas += p_elastic->XSec();
  }
  if (p_soft_diffractive) {
    p_soft_diffractive->Initialise();
    m_xsec += p_soft_diffractive->XSec();
    m_xsec_diff += p_soft_diffractive->XSec();
  }
  msg_Info()<<METHOD<<" with sigma = "<<m_xsec<<" pb\n";
}

void Event_Generator::Reset() {
  if (p_active) p_active->Reset();
  m_thisevent = m_runmode;
  m_mustinit = 1;
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
  if (blobs->size()==1) {
    (*blobs)[0]->AddData("WeightsMap",new ATOOLS::Blob_Data<ATOOLS::Weights_Map>(m_xsec));
    (*blobs)[0]->AddData("Weight_Norm",new ATOOLS::Blob_Data<double>(1.));
    (*blobs)[0]->AddData("Trials",new ATOOLS::Blob_Data<double>(1));
  }

  switch (MBpars.RunMode()) {
  case run_mode::inelastic_events:
    p_active = p_inelastic;
    break; 
  case run_mode::elastic_events:
    p_active = p_elastic;
    break;
  case run_mode::soft_diffractive_events:
    p_active = p_soft_diffractive;
  case run_mode::all_min_bias:
    if (m_mustinit && m_eventsuccessful) {
        double R(ran->Get());
        if (R < m_xsec_inel/m_xsec) {
            p_active = p_inelastic;
            //msg_Out()<<METHOD<<" This event will be inelastic.\n";
        }
        else if (R < (m_xsec_inel+m_xsec_elas)/m_xsec) {
          p_active = p_elastic;
          //msg_Out()<<METHOD<<" This event will be elastic.\n";
        }
        else {
          p_active = p_soft_diffractive;
          //msg_Out()<<METHOD<<" This event will be diffractive.\n";
        }
    }
    m_mustinit = 0;
  }
  int outcome = p_active->GenerateEvent(blobs,false);
  if (outcome == -1) {
      m_eventsuccessful = false;
      return 0;
  }
  else {
      m_eventsuccessful = true;
      return outcome;
  }
}

void Event_Generator::Test(const std::string & dirname) {
  msg_Info()<<METHOD<<": Starting.\n";
  p_inelastic->Test(dirname);
}
