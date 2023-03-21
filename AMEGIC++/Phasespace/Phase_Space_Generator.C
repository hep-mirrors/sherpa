#include <stdio.h>
#include <stdlib.h>
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "AMEGIC++/Phasespace/Channel_Generator.H"
#include "AMEGIC++/Phasespace/Channel_Generator_NPV.H"
#include "AMEGIC++/Phasespace/Channel_Generator3V.H"
#include "AMEGIC++/Phasespace/Channel_Generator_UniV.H"
#include "AMEGIC++/Phasespace/Channel_Generator3_NPV.H"
#include "AMEGIC++/Phasespace/Channel_Generator_KK.H"
#include "AMEGIC++/Phasespace/Channel_Generator_Decays.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "AMEGIC++/String/String_Library.H"

#include "MODEL/Main/Running_AlphaQED.H"


using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS; 
using namespace std;

Phase_Space_Generator::Phase_Space_Generator(int _nin, int _nout) :
  nin(_nin), nout(_nout), m_mode{ 1 }
{
  static int ishadronic{ -1 };
  if (ishadronic == -1) {
    ishadronic = (rpa->gen.Beam1().IsHadron() && rpa->gen.Beam2().IsHadron());
    auto s = Settings::GetMainSettings()["AMEGIC"];
    s["INTEGRATOR"].SetDefault((ishadronic) ? 6 : 7);
    s["RS_INTEGRATOR"].SetDefault(7);
  }
}

bool Phase_Space_Generator::Construct(std::list<std::string>* liblist,string _pathID,string _pID,
				      ATOOLS::Flavour* fl,Process_Base * proc)
{
  auto s = Settings::GetMainSettings()["AMEGIC"];
  path   = _pathID;
  pathID = _pathID + string("/") + _pID;
  pID    = string("P")+_pID;

  int ngraph = proc->NumberOfDiagrams();
  if (ngraph<=0) {
    msg_Error()<<"Error in Phase_Space_Generator::Construct for "<<proc->Name()<<endl;
    Abort();
  }

  string lmapname = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+pathID+string("/fsrchannels");
  string mapname  = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+path+string("/fsrchannels.map");

  int inttype{ -1 };
  if (proc->Info().Has(nlo_type::real)) {
    inttype = s["RS_INTEGRATOR"].Get<int>();
  } else {
    inttype = s["INTEGRATOR"].Get<int>();
  }
  if (nout==1) return 0;
  if (nin==1&&nout==2) return 0;
  if (inttype<4 && !(inttype>1 && nout==2)) return 0;
  if (inttype==2) inttype=6;
  if (inttype==3) inttype=7;
  if (inttype>20) return 0;

  if (My_In_File::FileInDB(lmapname)) return 1-GetLibList(liblist);

  int newchannels = 0;
  //int extrachannel = 0;
  My_Out_File lmf(lmapname);
  lmf.Open();
  int cnt=0;
  My_Out_File mf(mapname);
  My_In_File imf(mapname);
  if (!My_In_File::FileInDB(mapname)) {
    mf.Open();
    mf.Close();
  }
  else {
    char buffer[buffersize];
    imf.Open();
    for (;*imf;cnt++) imf->getline(buffer,buffersize);
    imf.Close();
    cnt--;
  }

  string fsrpath0= string("fsrchannels");
  string fsrpath = fsrpath0;
  char hlp[4];
  snprintf(hlp,4,"%i",nout);
  fsrpath += string(hlp);
  fsrpath0 = fsrpath;
  if (cnt>=maxchannels) {
    snprintf(hlp,4,"_%i",cnt/maxchannels);
    fsrpath = fsrpath+string(hlp);
  }
  string fsrp = path+string("/")+fsrpath;

  bool kk_fs=false;
  for (int i=0;i<nout;i++){
    if (fl[i+nin].IsKK()) kk_fs=true;
  }
  int ng = 2;
  if (inttype==4 || kk_fs || inttype==7) ng=1;
  if (proc->OSDecays()>0) ng=1;

  for (int i=0;i<ngraph;i++) {
    if (proc->IsFreeOfFourVertex(proc->Diagram(i))) {
      for(int j=0;j<ng;j++) {
        Channel_Generator_Base *cg(NULL);
	if (nin==1 && nout>2) cg = new Channel_Generator_Decays(nin,nout,proc->Diagram(i),0);
	else {
	  if (kk_fs) {
	    cg = new Channel_Generator_KK(nin,nout,proc->Diagram(i),0);
	  }
	  else {
	    if (inttype==6) {
	      if (j==0) cg = new Channel_Generator3V(nin,nout,proc->Diagram(i),0);
	      else cg = new Channel_Generator3_NPV(nin,nout,proc->Diagram(i),0);
	    }
	    else {
	      if (inttype==7) {
		cg = new Channel_Generator_UniV(nin,nout,proc->Diagram(i),0);
	      }
	      else {
		if (j==0) cg = new Channel_Generator(nin,nout,proc->Diagram(i),0);
		else cg = new Channel_Generator_NPV(nin,nout,proc->Diagram(i),0);
	      }
	    }
	  }
	}
	if (!cg) THROW(fatal_error, "No channel generated.");
	for (int k=0;k<cg->NumberOfChannels();k++) {
	  string chID = cg->CreateChannelID(k);
	  *lmf<<chID<<endl;
	  if (!RSearchInDB(mapname,chID)) {
	    bool hit;
	    do {
	      if (nin==2) snprintf(procname,12,"C%i_%i",nout,cnt);
	      else snprintf(procname,12,"CD%i_%i",nout,cnt);
	      string help = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+fsrp+string(procname);
	      hit = My_In_File::FileInDB(help);
	      if (hit) cnt++;
	    } while (hit);
	    
	  // making directory
	    if (cnt%maxchannels==0) {
	      if (cnt>0) {
		snprintf(hlp,4,"_%i",cnt/maxchannels);
		fsrpath = fsrpath0 + string(hlp);
		fsrp = path+string("/")+fsrpath;
	      }
	      ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+fsrp,true);
	      String_Library slib(1);
	      slib.InitMakefile(fsrp);
	    }
	    
	    int  rannumber;
	    //cg->SetName(string(procname));
	    rannumber    = cg->MakeChannel(k,cnt,fsrp,pID);
	    if (nout==1) rannumber=1;
	    if (!cg->Valid() && msg_LevelIsTracking()) PRINT_INFO("Channel "<<procname<<" kicked because of decoupled particle");
	    if (rannumber>0 && cg->Valid()) {
	      string makefilename = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+fsrp+string("/CMakeLists.txt");
	      AddToCMakefile(makefilename,fsrp,procname);
	      cnt++;
	      newchannels = 1;
	    }
            std::string mapstr, maptmp;
            if (imf.Open())
              for (getline(*imf,maptmp);imf->good(); getline(*imf,maptmp))
                mapstr+=maptmp+"\n";
            imf.Close();
            mf.Open();
            *mf<<mapstr;
            *mf<<chID<<": "<<fsrpath<<"/"<<string(procname)<<endl;
	    mf.Close();
	  }
	  else {
	    if (!newchannels) {
	      if (chID[0]!='%') {
		int pos = chID.find(string(": "));
		chID = chID.substr(pos+2);
		liblist->push_back(chID);
	      }	
	    }
	  }
	}
	delete cg;
      }
    }
  }
  lmf.Close();
  return newchannels;
}



void Phase_Space_Generator::AddToCMakefile(string makefilename,string pathID,string fileID)
{
  size_t hit=pathID.find("/");
  string base=pathID.substr(0,hit);
  string subdirname=pathID.substr(hit+1);

  if (!IsFile(makefilename)) {
    ofstream file(makefilename.c_str());

    file<<"set(libProc_"<<subdirname<<"_la_SOURCES"<<endl;
    file<<fileID<<".C"<<endl;
    file<<"CG.C"<<endl;
    file<<")"<<endl;
    file<<"add_library(Proc_"<<subdirname<<" SHARED ${libProc_"<<subdirname<<"_la_SOURCES})"<<endl;
    file<<"amegic_handle_shared_library(Proc_"<<subdirname<<")"<<endl;
    file<<"install(TARGETS Proc_"<<subdirname<<" DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})"<<endl;
    
    
    ofstream cgfile((rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+pathID+"/CG.C").c_str());
    cgfile<<"#include \"PHASIC++/Channels/Channel_Generator.H\"\n"
	  <<"#include \"PHASIC++/Channels/Multi_Channel.H\"\n"
	  <<"#include \"PHASIC++/Process/Process_Base.H\"\n"
	  <<"#include \"PHASIC++/Process/ME_Generator_Base.H\"\n"
	  <<"#include \"PHASIC++/Channels/Single_Channel.H\"\n"
	  <<"#include \"ATOOLS/Org/Library_Loader.H\"\n"
	  <<"#include \"ATOOLS/Org/Run_Parameter.H\"\n"
	  <<"#include \"PHASIC++/Main/Phase_Space_Handler.H\"\n"
	  <<"#include \"PHASIC++/Main/Process_Integrator.H\"\n";
    cgfile<<"\nusing namespace PHASIC;\nusing namespace ATOOLS;\n\n";
    cgfile<<"#define PTS long unsigned int\n#define PT(ARG) (PTS)(ARG)\n"
	  <<"typedef PHASIC::Single_Channel *(*Lib_Getter_Function)\n"
	  <<"  (int nin,int nout,ATOOLS::Flavour* fl,\n"
	  <<"   ATOOLS::Integration_Info * const info,PHASIC::Phase_Space_Handler *psh);\n";
    cgfile<<"\nnamespace PHASIC {\n";
    cgfile<<"  class "<<subdirname<<"_Channel_Generator: public Channel_Generator {\n";
    cgfile<<"  public:\n    "<<subdirname<<"_Channel_Generator(const Channel_Generator_Key &key):\n"
	  <<"      Channel_Generator(key) {}\n";
    cgfile<<"    Single_Channel *LoadChannel(int nin,int nout,Flavour* fl,"
	  <<"const std::string &pID,Phase_Space_Handler *psh)\n"
	  <<"    {\n      size_t pos(pID.find(\"/\"));\n"
	  <<"      s_loader->AddPath(rpa->gen.Variable(\"SHERPA_LIB_PATH\"));\n"
	  <<"      Lib_Getter_Function gf = (Lib_Getter_Function)\n"
	  <<"        PT(s_loader->GetLibraryFunction(\"Proc_"<<subdirname<<"\",\"Getter_\"+pID));\n"
	  <<"      if (gf==NULL) return NULL;\n"
	  <<"      return gf(nin,nout,fl,psh->GetInfo(),psh);\n    }\n";
    cgfile<<"    int GenerateChannels()\n"
	  <<"    {\n      int nin=p_proc->NIn(), nout=p_proc->NOut();\n"
	  <<"      Flavour *fl=(Flavour*)&p_proc->Flavours().front();\n"
	  <<"      Phase_Space_Handler *psh=&*p_proc->Integrator()->PSHandler();\n"
	  <<"      p_mc->Add(LoadChannel(nin,nout,fl,\""<<fileID<<"\",psh));\n"
	  <<"      return 0;\n    }\n";
    cgfile<<"  };\n}\n\n";
    cgfile<<"DECLARE_GETTER("<<subdirname<<"_Channel_Generator,\""
	  <<subdirname<<"\",Channel_Generator,Channel_Generator_Key);\n";
    cgfile<<"Channel_Generator *ATOOLS::Getter<Channel_Generator,"
	  <<"Channel_Generator_Key,"<<subdirname<<"_Channel_Generator>::\n"
	  <<"operator()(const Channel_Generator_Key &args) const "
	  <<"{ return new "<<subdirname<<"_Channel_Generator(args); }\n";
    cgfile<<"void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,"
	  <<subdirname<<"_Channel_Generator>::\n"
	  <<"PrintInfo(std::ostream &str,const size_t width) const "
	  <<"{ str<<\""<<subdirname<<"\"; }\n"<<std::flush;
  }
  else {
    ifstream from(makefilename.c_str());
    ofstream to((makefilename+string(".tmp")).c_str());  

    string buffer;
    for (;from;) {
      getline(from,buffer);
      to<<buffer<<endl;
      if (buffer.find(std::string("set(libProc_")+subdirname+"_la_SOURCES")!=string::npos) {to<<fileID<<".C"<<endl;}
    }
    from.close();
    to.close();

    Move(makefilename+".tmp",makefilename);
    {
      std::string fname=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+pathID+"/CG.C";
      ifstream from(fname.c_str());
      ofstream to((fname+".tmp").c_str());
      bool first=true;
      string buffer;
      for (;from;) {
        getline(from,buffer);
        if (first && buffer.find("p_mc->Add")!=string::npos) {
          to<<"      p_mc->Add(LoadChannel(nin,nout,fl,\""<<fileID<<"\",psh));\n";
          first=false;
        }
        to<<buffer<<endl;
      }
      from.close();
      to.close();

      Move(fname+".tmp",fname);
    }
  }
}

bool Phase_Space_Generator::GetLibList(std::list<std::string>* liblist)
{
  string chlname   = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+pathID + string("/fsrchannels");
  string chmapname = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+path   + string("/fsrchannels.map");

  My_In_File chlist(chlname);
  chlist.Open();
  if (!My_In_File::FileInDB(chmapname)) {
    msg_Error()<<"Error in Phase_Space_Generator:" <<chmapname<<" not found."<<endl;
    return 0;
  }

  char buffer[buffersize];
  string libname;
  for(;*chlist;) {
    chlist->getline(buffer,buffersize);
    libname = string(buffer);
    if (*chlist && libname[0]!='%') {
      if (!RSearchInDB(chmapname,libname) || libname.find(": ")==string::npos) {
	msg_Error()<<"Error in Phase_Space_Generator:"  <<"Mapping for "<<libname<<" not found."<<endl;
	return 0;
      }

      if (libname[0]!='%') {
	int pos = libname.find(string(": "));
	libname = libname.substr(pos+2);
      
	liblist->push_back(libname);
      }
    }
  }
  chlist.Close();
  return 1;
}

bool Phase_Space_Generator::IsFile(string &filename)
{
  ifstream from;
  from.open(filename.c_str());
  if (from) return 1;
  return 0;
}

bool Phase_Space_Generator::Search(ifstream &from,string search)
{
  char buffer[buffersize];
  for(;from;) {
    from.getline(buffer,buffersize);    
    if (string(buffer).find(string(search))!=string::npos) return 1;
  }
  return 0;
}

bool Phase_Space_Generator::RSearch(ifstream &from,string &search)
{
  char buffer[buffersize];
  for(;from;) {
    from.getline(buffer,buffersize);    
    if (string(buffer).find(string(search))!=string::npos) {
      search = string(buffer);
      return 1;
    }
  }
  return 0;
}

int  Phase_Space_Generator::Search(string file,string search)
{  

  ifstream from;
  //search name  
  from.open(file.c_str());

  char buffer[buffersize];

  for(;from;) {
    from.getline(buffer,buffersize);    
    if (string(buffer).find(string(search))!=string::npos) {
      from.close();
      return 1;
    }
  }
  from.close();
  return 0;
}

int  Phase_Space_Generator::RSearchInDB(string file,string &search)
{
  My_In_File from(file);

  if (from.Open()) {
    return RSearch(*from, search);
  }
  return 0;
}

void Phase_Space_Generator::Copy(string sfrom,string sto)
{
  ifstream from;
  ofstream to;
  
  from.open(sfrom.c_str());
  to.open(sto.c_str()); 

  char ch {0};
  while(from.get(ch)) to.put(ch);
  from.close();
  to.close();  

  remove(sfrom.c_str());
}



