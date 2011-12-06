#include "SHERPA/Tools/Input_Output_Handler.H"

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/Output_Sherpa.H"
#include "SHERPA/Tools/Output_RootNtuple.H"
#include "SHERPA/Tools/Output_HepEvt.H"
#include "SHERPA/Tools/Output_D0_HepEvt.H"
#include "SHERPA/Tools/HepEvt_Interface.H"
#include "SHERPA/Tools/Output_LHEF.H"
#ifdef USING__HEPMC2
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "HepMC/GenEvent.h"
#include "SHERPA/Tools/Output_HepMC2_Genevent.H"
#include "SHERPA/Tools/Output_HepMC2_Short.H"
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
  p_hepevt(NULL),
  p_eventhandler(NULL), p_mehandler(NULL)
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
  std::string stag(rpa->gen.Variable("RNG_SEED"));
  while (stag.find(' ')!=std::string::npos) stag.replace(stag.find(' '),1,"-");
  dr->AddTag("RNG_SEED",stag);
  string sherpaoutput=dr->GetValue<string>("SHERPA_OUTPUT",string(""));
  string rootntupleoutput=dr->GetValue<string>("ROOTNTUPLE_OUTPUT",string(""));
  string hepevtoutput=dr->GetValue<string>("HEPEVT_OUTPUT",string(""));
  string d0hepevtoutput=dr->GetValue<string>("D0_HEPEVT_OUTPUT",string(""));
  string hepmc2genevent=dr->GetValue<string>("HEPMC2_GENEVENT_OUTPUT",string(""));
  string lheoutput=dr->GetValue<string>("LHEF_OUTPUT",string(""));
  string hepmc2short=dr->GetValue<string>("HEPMC2_SHORT_OUTPUT",string(""));
  string evtpath = dr->GetValue<string>
    ("EVT_FILE_PATH",rpa->gen.Variable("SHERPA_RUN_PATH"));
  int precision       = dr->GetValue<int>("OUTPUT_PRECISION",12);
  m_outmode = dr->GetValue<string>("EVENT_MODE",string("Sherpa"));
  dr->SetAllowUnits(true);
  m_filesize = dr->GetValue<long int>("FILE_SIZE",std::numeric_limits<int>::max());

  if (!sherpaoutput.empty()) {
    m_outmap["SHERPA"]=
      new Output_Sherpa(evtpath+"/"+sherpaoutput,".evts", precision);
  }
  if (!rootntupleoutput.empty()) {
#ifdef USING__ROOT
    if (rootntupleoutput[0] == '/'){
      m_outmap["ROOTNTUPLE"]= new Output_RootNtuple(rootntupleoutput,".root", precision,m_filesize);
    } else {
      m_outmap["ROOTNTUPLE"]= new Output_RootNtuple(evtpath+"/"+rootntupleoutput,".root", precision,m_filesize);
    }
#else
    THROW(fatal_error,"ROOTNTUPLE format can only be created when Sherpa was linked "
          +string("with root, please read our Manual for more information."));
#endif
  }
  if (!hepevtoutput.empty()) {
    m_outmap["HEPEVT"]=
      new Output_HepEvt(evtpath+"/"+hepevtoutput,".hepevt", precision);
  }
  if (!d0hepevtoutput.empty()) {
    m_outmap["D0_HEPEVT"]=
      new Output_D0_HepEvt(evtpath+"/"+d0hepevtoutput,".d0.hepevt", precision);
  }
  if (!lheoutput.empty()) {
    m_outmap["LHEF"]=
      new Output_LHEF(evtpath+"/"+lheoutput,".lhe", precision);
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
  if (!hepmc2short.empty()) {
#ifdef USING__HEPMC2
    m_outmap["HEPMC2_SHORT"]= new
      Output_HepMC2_Short(evtpath+"/"+hepmc2short,".hepmc",precision);
#else
    THROW(fatal_error,"HepMC format can only be created when Sherpa was linked "
          +string("with HepMC2, please read our Howto for more information."));
#endif
  }

  if (m_outmode=="HepMC" || m_outmode=="HepMC_Short") {
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
  else if (m_outmode=="HepMC_Short") {
#ifdef USING__HEPMC2
    p_hepmc2->Sherpa2ShortHepMC(blobs);
#else
    THROW(fatal_error,"HepMC format can only be created when Sherpa was linked "
          +string("with HepMC2, please read our Howto for more information."));
#endif
  }
  if (!msg_LevelIsEvents()) return;
  if (m_outmode=="Sherpa" || m_outmode=="LHEF") {
    if (!blobs->empty()) {
      msg_Out()<<"  -------------------------------------------------  "<<std::endl;
      for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) 
        msg_Out()<<*(*blit)<<std::endl;
      msg_Out()<<"  -------------------------------------------------  "<<std::endl;
    }
    else msg_Out()<<"  ******** Empty event ********  "<<std::endl;
  }
  else if (m_outmode=="HepMC" || m_outmode=="HepMC_Short") {
#ifdef USING__HEPMC2
    // if sub list present, print only those
    if (p_hepmc2->GenSubEventList().size()>0) {
      for (size_t i=0;i<p_hepmc2->GenSubEventList().size();++i) {
        p_hepmc2->GenSubEvent(i)->print(msg->Out());
      }
    }
    else {
      p_hepmc2->GenEvent()->print(msg->Out());
    }
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
  double weight=blobs->Weight();
  if (p_mehandler) weight/=p_mehandler->TotalXS()*rpa->Picobarn();

  double xs(1.0), xserr(0.0);
  if (p_eventhandler) {
    xs=p_eventhandler->TotalXS();
    xserr=p_eventhandler->TotalErr();
  }
  for (map<string,Output_Base *>::iterator oit=m_outmap.begin();
       oit!=m_outmap.end();oit++) {
    oit->second->SetXS(xs, xserr);
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
