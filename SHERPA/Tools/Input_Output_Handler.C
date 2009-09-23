#include "Input_Output_Handler.H"

#include "Blob_List.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Output_Base.H"
#include "Output_Sherpa.H"
#include "Output_HepEvt.H"
#include "Output_D0_HepEvt.H"
#include "HepEvt_Interface.H"
#ifdef USING__HEPMC2
#include "HepMC2_Interface.H"
#include "HepMC/GenEvent.h"
#include "Output_HepMC2_Genevent.H"
#endif

#include <stdio.h>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

extern "C" {
  void outhepevt_();
}

Input_Output_Handler::Input_Output_Handler(Data_Reader* dr) :
  m_evtnumber(0), m_evtcount(0),
#ifdef USING__HEPMC2
  p_hepmc2(NULL),
#endif
  p_hepevt(NULL)
{
  InitialiseOutput(dr);
}

Input_Output_Handler::~Input_Output_Handler() 
{
#ifdef USING__HEPMC2
  if (p_hepmc2)       { delete p_hepmc2; p_hepmc2 = NULL; }
#endif
  if (p_hepevt!=NULL) { delete p_hepevt; p_hepevt=NULL; }

  for (map<string,Output_Base *>::iterator oit=m_outmap.begin();
       oit!=m_outmap.end();oit++) {
    delete oit->second;
  }
  m_outmap.clear();
}

bool Input_Output_Handler::InitialiseOutput(Data_Reader* dr) {
  string sherpaoutput=dr->GetValue<string>("SHERPA_OUTPUT",string(""));
  string hepevtoutput=dr->GetValue<string>("HEPEVT_OUTPUT",string(""));
  string d0hepevtoutput=dr->GetValue<string>("D0_HEPEVT_OUTPUT",string(""));
  string hepmc2genevent=dr->GetValue<string>("HEPMC2_GENEVENT_OUTPUT",string(""));
  string evtpath = dr->GetValue<string>
    ("EVT_FILE_PATH",rpa.gen.Variable("SHERPA_RUN_PATH"));
  int precision       = dr->GetValue<int>("OUTPUT_PRECISION",12);
  m_outmode = dr->GetValue<string>("EVENT_MODE",string("Sherpa"));
  m_filesize = dr->GetValue<int>("FILE_SIZE",1000);

  if (!sherpaoutput.empty()) {
    m_outmap["SHERPA"]=
      new Output_Sherpa(evtpath+"/"+sherpaoutput,".evts", precision);
  }
  if (!hepevtoutput.empty()) {
    m_outmap["HEPEVT"]=
      new Output_HepEvt(evtpath+"/"+hepevtoutput,".hepevt", precision);
  }
  if (!d0hepevtoutput.empty()) {
    m_outmap["D0_HEPEVT"]=
      new Output_D0_HepEvt(evtpath+"/"+d0hepevtoutput,".d0.hepevt", precision);
  }
  if (!hepmc2genevent.empty()) {
#ifdef USING__HEPMC2
    m_outmap["HEPMC2_GENEVENT"]= new
      Output_HepMC2_Genevent(evtpath+"/"+hepmc2genevent,".hepmc2g",precision);
#else
    THROW(fatal_error,"HepMC format can only be created when Sherpa was linked "
          +string("with HepMC2, please read our Howto for more information."));
#endif
  }

  if (m_outmode=="HepMC") {
#ifdef USING__HEPMC2
    if (p_hepmc2==NULL) p_hepmc2 = new HepMC2_Interface();
#else
    THROW(fatal_error,"HepMC format can only be created when Sherpa was linked "
          +string("with HepMC2, please read our Howto for more information."));
#endif
  }
  else if (m_outmode=="HepEvt") {
    if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface();
  }

  for (map<string,Output_Base *>::iterator oit=m_outmap.begin();
       oit!=m_outmap.end();oit++) {
    oit->second->Header();
  }

  return true;
}

void Input_Output_Handler::PrintEvent(ATOOLS::Blob_List *const blobs) {
  if (m_outmode=="HepMC") {
#ifdef USING__HEPMC2
    p_hepmc2->Sherpa2HepMC(blobs);
#else
    THROW(fatal_error,"HepMC format can only be created when Sherpa was linked "
          +string("with HepMC2, please read our Howto for more information."));
#endif
  }
  if (!msg_LevelIsEvents()) return;
  if (m_outmode=="Sherpa") {
    if (!blobs->empty()) {
      msg_Out()<<"  -------------------------------------------------  "<<std::endl;
      for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) 
        msg_Out()<<*(*blit)<<std::endl;
      msg_Out()<<"  -------------------------------------------------  "<<std::endl;
    }
    else msg_Out()<<"  ******** Empty event ********  "<<std::endl;
  }
  else if (m_outmode=="HepMC") {
#ifdef USING__HEPMC2
    p_hepmc2->GenEvent()->print(msg->Out());
#else
    THROW(fatal_error, "HepMC format can only be created when Sherpa was linked"
          +string(" with HepMC2, please read our Howto for more information."));
#endif
  }
  else if (m_outmode=="HepEvt") {
    p_hepevt->Sherpa2HepEvt(blobs);
    p_hepevt->WriteFormatedHepEvt(msg->Out(),p_hepevt->Nhep());
  }
  else {
    msg_Error()<<"Error in "<<METHOD<<std::endl
               <<"   Unknown Output format : "<<m_outmode<<std::endl
               <<"   No output, continue run ..."<<std::endl;
  }
}

bool Input_Output_Handler::OutputToFormat(ATOOLS::Blob_List *const blobs)
{
  ResetInterfaces();
  
  double weight=1.0;
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit){
    if ((*blit)->Type()==btp::Signal_Process) {
      Blob_Data_Base *info((**blit)["ME_Weight"]);
      if (info!=NULL) weight*=info->Get<double>();
      info=(**blit)["Process_Weight"];
      if (info!=NULL) weight/=info->Get<double>();
    }
  }

  for (map<string,Output_Base *>::iterator oit=m_outmap.begin();
       oit!=m_outmap.end();oit++) {
    oit->second->Output(blobs, weight);
  }
  
  m_evtnumber++;
  m_evtcount++;
  
  if (m_evtcount%m_filesize==0) {
    m_evtcount = 0;
    string number(ToString(int(m_evtnumber/m_filesize)));
    for (map<string,Output_Base *>::iterator oit=m_outmap.begin();
         oit!=m_outmap.end();oit++) {
      oit->second->Footer(number);
      oit->second->ChangeFile(number);
      oit->second->Header();
    }
  }
  return true;
}

void Input_Output_Handler::ResetInterfaces() {
}
