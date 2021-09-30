#include "CFPSHOWER++/Tools/CFP_Parameters.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;


CFP_Parameters * CFPSHOWER::cfp_pars = NULL;

CFP_Parameters::CFP_Parameters() {}

bool CFP_Parameters::Init(Default_Reader * dataread) {
  msg_Out()<<"Entering "<<METHOD<<"==================================\n";

  //                                      1 = Alaric, 2 = CS
  m_switches[string("kinematics")]      = dataread->GetValue<int>("CSS_KIN_SCHEME",2);  
  m_switches[string("kfactor")]         = dataread->GetValue<int>("CSS_KFACTOR_SCHEME",1); 
  m_switches[string("endpoint")]        = dataread->GetValue<int>("CSS_ENDPOINT_SCHEME",0);
  m_switches[string("softcorrections")] = dataread->GetValue<int>("CSS_SOFTCORRECTION_SCHEME",0);
  m_switches[string("couplings")]       = dataread->GetValue<int>("CSS_COUPLING_SCHEME",1);
  m_switches[string("ME_corrections")]  = dataread->GetValue<int>("CSS_ME_CORRECTION",0);
  m_switches[string("SF_order")]        = dataread->GetValue<int>("CSS_SF_ORDER",1);
  m_switches[string("max_emissions")]   = dataread->GetValue<int>("CSS_MAXEM",100000);
  m_switches[string("max_particles")]   = dataread->GetValue<int>("CSS_MAXPART",100000);
  // kt2_pipj_for_gqq = 11, kt2_all          = 10
  // T_pipj_for_gqq   = 1,  T_all            = 0
  m_switches[string("muR_scheme")]      = dataread->GetValue<int>("CSS_MUR_SCHEME",1);
  //                                      1 = soft, 2 = coll, 3 = both
  m_switches[string("Log_Type")]        = dataread->GetValue<int>("CSS_LOG_TYPE",3);
  //                                      1 = q->qg, 2 = g->gg, 4 = g->qq, 7 = all
  m_switches[string("SF_Type")]         = dataread->GetValue<int>("CSS_SF_TYPE",8);
                                             
  m_parameters[string("PDF_min")]       = dataread->GetValue<double>("CSS_PDF_MIN",1.e-4); 
  m_parameters[string("PDF_min_X")]     = dataread->GetValue<double>("CSS_PDF_MIN_X",1.e-4);
  m_parameters[string("pt2min(FS)")]    = dataread->GetValue<double>("CSS_FS_PT2MIN",1.);
  m_parameters[string("pt2min(IS)")]    = dataread->GetValue<double>("CSS_IS_PT2MIN",1.);
  m_parameters[string("k_alpha(FS)")]   = dataread->GetValue<double>("CSS_FS_AS_FAC",1.);
  m_parameters[string("k_alpha(IS)")]   = dataread->GetValue<double>("CSS_IS_AS_FAC",1.);
  m_parameters[string("k_muF")]         = dataread->GetValue<double>("CSS_IS_PDF_FAC",1.);
  Output();
  msg_Out()<<"Leaving "<<METHOD<<"==================================\n";
  return true;
  //auto s = Settings::GetMainSettings();
  /*
  m_parameters["recalc_fac"]    = s["CSS_RECALC_FACTOR"].SetDefault(4.0).Get<double>();
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
