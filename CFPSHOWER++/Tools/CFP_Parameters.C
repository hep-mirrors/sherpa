#include "CFPSHOWER++/Tools/CFP_Parameters.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;


CFP_Parameters * CFPSHOWER::cfp_pars = NULL;

CFP_Parameters::CFP_Parameters() {}

bool CFP_Parameters::Init()
{
  msg_Out()<<"Entering "<<METHOD<<"==================================\n";
  m_switches[string("kinematics")]      = 2;  //CS Kinematics
  m_switches[string("kfactor")]         = 1; 
  m_switches[string("endpoint")]        = 0;
  m_switches[string("softcorrections")] = 0;
  m_switches[string("couplings")]       = 1;
  m_switches[string("ME_corrections")]  = 0;
  m_switches[string("SF_order")]        = 1;
  m_switches[string("Log_Type")]        = 2;
  m_switches[string("max_emissions")]   = 1;
  m_switches[string("max_particles")]   = 100000;
  m_switches[string("muR_scheme")]      = 1;
  m_parameters[string("PDF_min")]       = 1.e-4; 
  m_parameters[string("PDF_min_X")]     = 1.e-4;
  m_parameters[string("pt2min(FS)")]    = 1.;
  m_parameters[string("pt2min(IS)")]    = 1.;
  m_parameters[string("k_alpha(FS)")]   = 1.;
  m_parameters[string("k_alpha(IS)")]   = 1.;
  m_parameters[string("k_muR")]         = 1.;
  m_parameters[string("k_muF")]         = 1.;
  Output();
  msg_Out()<<"Leaving "<<METHOD<<"==================================\n";
  return true;
  //auto s = Settings::GetMainSettings();
  /*
  m_switches["kinematics"]      = s["CSS_KIN_SCHEME"].SetDefault(1).Get<int>();
  m_switches["kfactor"]         = s["CSS_KFACTOR_SCHEME"].SetDefault(0).Get<int>();
  m_switches["endpoint"]        = s["CSS_ENDPOINT_SCHEME"].SetDefault(0).Get<int>();
  m_switches["softcorrections"] = s["CSS_SOFTCORRECTION_SCHEME"].SetDefault(0).Get<int>();
  m_switches["couplings"]       = s["CSS_COUPLING_SCHEME"].SetDefault(1).Get<int>();
  m_switches["ME_corrections"]  = s["CSS_ME_CORRECTION"].SetDefault(0).Get<int>();
  m_switches["SF_order"]        = s["CSS_SF_ORDER"].SetDefault(1).Get<int>();
  m_switches["Log_Type"]        = s["CSS_LOG_TYPE"].SetDefault(2).Get<int>();
  m_switches["max_emissions"]   = s["CSS_MAXEM"].SetDefault(10000).Get<int>();
  m_switches["max_particles"]   = s["CSS_MAXPART"].SetDefault(10000).Get<int>();
  //m_switches["muR_scheme"]    = s["CSS_MUR_SCHEME"].SetDefault(1).Get<int>();
  m_parameters["recalc_fac"]    = s["CSS_RECALC_FACTOR"].SetDefault(4.0).Get<double>();
  m_parameters["PDF_min"]       = s["CSS_PDF_MIN"].SetDefault(1.0e-4).Get<double>();
  m_parameters["PDF_min_X"]     = s["CSS_PDF_MIN_X"].SetDefault(1.0e-2).Get<double>();
  m_parameters["NLO_enhance"]   = s["CSS_TC_ENHANCE"].SetDefault(1.0).Get<double>();
  */
}

int CFP_Parameters::operator[](string keyword) 
{
  map<string,int>::iterator siter = m_switches.find(keyword);
  if (siter!=m_switches.end()) return siter->second;
  msg_Error()<<"Error in CFP_Switches("<<keyword<<"), "
	     <<"in total "<<m_switches.size()<<" switches.\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  exit(1);
  return 0;
}

double CFP_Parameters::operator()(string keyword) 
{
  map<string,double>::iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  msg_Error()<<"Error in CFP_Parameters("<<keyword<<") "
	     <<"in "<<m_parameters.size()<<".\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  exit(1);
  return 0.;
}

void CFP_Parameters::Output() {
  msg_Out()<<"*********************************************\n"
	   <<METHOD<<" has "<<m_switches.size()<<" switches and "
	   <<m_parameters.size()<<" parameters:\n";
  for (map<string,int>::iterator siter=m_switches.begin();
       siter!=m_switches.end();siter++) {
    msg_Out()<<"  *** "<<siter->first<<" = "<<siter->second<<" (switch)\n";
  }
  for (map<string,double>::iterator piter=m_parameters.begin();
       piter!=m_parameters.end();piter++) {
    msg_Out()<<"  *** "<<piter->first<<" = "<<piter->second<<" (parameter)\n";
  }
  msg_Out()<<"*********************************************\n";
}
