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
#include "AMEGIC++/String/String_Library.H"

#include "MODEL/Main/Running_AlphaQED.H"


using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS; 
using namespace std;

Phase_Space_Generator::Phase_Space_Generator(int _nin,int _nout) : nin(_nin), nout(_nout) 
{
  m_mode=1;
}

bool Phase_Space_Generator::Construct(std::list<std::string>* liblist,string _pathID,string _pID,
				      ATOOLS::Flavour* fl,Process_Base * proc)
{ 
  path   = _pathID;
  pathID = _pathID + string("/") + _pID;
  pID    = string("P")+_pID;

  int ngraph = proc->NumberOfDiagrams();
  if (ngraph<=0) {
    msg_Error()<<"Error in Phase_Space_Generator::Construct for "<<proc->Name()<<endl;
    abort();
  }

  string lmapname = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pathID+string("/fsrchannels");
  string mapname  = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+path+string("/fsrchannels.map");

  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(rpa.GetPath());
  dr.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  int inttype  = dr.GetValue<int>("INTEGRATOR",6);
  if (nin==1&&nout==2) return 0;
  if (inttype<4 && !(inttype>1 && nout==2)) return 0;
  if (inttype==2) inttype=6;
  if (inttype==3) inttype=7;
  if (inttype>20) return 0;

  if (IsFile(lmapname)) return 1-GetLibList(liblist);

  ATOOLS::MakeDir(rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+path);
  int newchannels = 0;
  //int extrachannel = 0;
  ofstream lmf;
  lmf.open(lmapname.c_str());
  int cnt=0;
  ofstream mf;
  ifstream imf;
  if (!IsFile(mapname)) {
    mf.open(mapname.c_str(),ios::out);
    mf.close();
  }
  else {
    char buffer[buffersize];
    imf.open(mapname.c_str());
    for (;imf;cnt++) imf.getline(buffer,buffersize);
    imf.close();
    cnt--;
  }

  string fsrpath0= string("fsrchannels");
  string fsrpath = fsrpath0;
  char hlp[4];
  sprintf(hlp,"%i",nout);
  fsrpath += string(hlp);
  fsrpath0 = fsrpath;
  if (cnt>=maxchannels) {
    sprintf(hlp,"_%i",cnt/maxchannels);
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
      for(int j=0;j<ng;j++){
	Channel_Generator_Base *cg;
	if (nin==1 && nout>2) cg = new Channel_Generator_Decays(nin,nout,proc->Diagram(i),0);
	else {
	  if (kk_fs){
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
	for (int k=0;k<cg->NumberOfChannels();k++) {
	  string chID = cg->CreateChannelID(k);
	  lmf<<chID<<endl;
	  if (!RSearch(mapname,chID)) {
	    bool hit;
	    do {
	      if (nin==2) sprintf(procname,"C%i_%i",nout,cnt);
	      else sprintf(procname,"CD%i_%i",nout,cnt);
	      string help = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+fsrp+string(procname);
	      hit = IsFile(help);
	      if (hit) cnt++;
	    } while (hit);
	    
	  // making directory
	    if (cnt%maxchannels==0) {
	      if (cnt>0) {
		sprintf(hlp,"_%i",cnt/maxchannels);
		fsrpath = fsrpath0 + string(hlp);
		fsrp = path+string("/")+fsrpath;
	      }
	      ATOOLS::MakeDir(rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+fsrp);
	      String_Library slib(1);
	      slib.InitMakefile(fsrp);
	    }
	    
	    int  rannumber;
	    //cg->SetName(string(procname));
	    rannumber    = cg->MakeChannel(k,cnt,fsrp,pID);
	    if (rannumber>0) {
	      string makefilename = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+fsrp+string("/Makefile.am");
	      AddToMakefileAM(makefilename,fsrp,procname);
	      cnt++;
	      newchannels = 1;
	    }
	    mf.open(mapname.c_str(),ios::out|ios::app);
	    mf<<chID<<": "<<fsrpath<<"/"<<string(procname)<<endl;
	    mf.close();
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
  lmf.close();
  return newchannels;
}



void Phase_Space_Generator::AddToMakefileAM(string makefilename,string pathID,string fileID)
{
  size_t hit=pathID.find("/");
  string base=pathID.substr(0,hit);
  string subdirname=pathID.substr(hit+1);

  if (!IsFile(makefilename)) {
    ofstream file(makefilename.c_str());

    file<<"lib_LTLIBRARIES = libProc_"<<subdirname<<".la"<<endl;
    file<<"libProc_"<<subdirname<<"_la_SOURCES = "<<'\\'<<endl;
    file<<"\t"<<fileID<<".C"<<endl;
    file<<"CURRENT_SHERPASYS = "<<ATOOLS::rpa.gen.Variable("SHERPA_INC_PATH")<<endl;
    file<<"INCLUDES = -I$(CURRENT_SHERPASYS)"<<endl;
    file<<"DEFS     = "<<endl;
  }
  else {
    ifstream from(makefilename.c_str());
    ofstream to((makefilename+string(".tmp")).c_str());  

    string buffer;
    string key=string("libProc_"+subdirname+"_la_SOURCES");
    for (;from;) {
      getline(from,buffer);
      to<<buffer<<endl;
      if (buffer.find(key)!=string::npos) {
	to<<"\t"<<fileID<<".C"<<'\\'<<endl;
      }
    }
    from.close();
    to.close();

    MoveFile(makefilename+".tmp",makefilename);
  }
}


void  Phase_Space_Generator::AddToMakefile(string makefilename,string pathID,string fileID)
{
  if (IsFile(makefilename)==0) {
    cerr<<makefilename.c_str()<<" is not available !"<<endl;
    return;
  }

  if (Search(makefilename,string(fileID)+string(".C"))) return;

  AddToMakefileAM(makefilename+string(".am"),pathID,fileID);

  ofstream to;  
  ifstream from;


  from.open(makefilename.c_str()); 
  to.open((makefilename+string(".tmp")).c_str());

  char buffer[buffersize];

  string pID;
  pID=pathID;
  for (short int i=pathID.length()-1;i>=0;i--) {
    if (pathID[i]=='/') {
      pID    = pathID.substr(i+1);
      break;
    }
  }

  string lib = string("libProc_")+pID;

  for(;from;) {
    from.getline(buffer,buffersize);
    if (string(buffer).find(lib+string("_la_SOURCES"))==0) {
      if (string(buffer).find(string("\\"))==string::npos) {
	//no backslash
	to<<buffer<<"\\"<<endl
	  <<"\t"<<(fileID+string(".C")).c_str()<<endl; 
      }
      else {
	to<<buffer<<endl
	  <<"\t"<<(fileID+string(".C")).c_str()<<" \\"<<endl; 
      }
    }
    else {
      if (string(buffer).find(lib+string("_la_OBJECTS"))==0) {
	if (string(buffer).find(string("\\"))==string::npos) {
	  //no backslash
	  to<<buffer<<"\\"<<endl
	    <<"\t"<<(fileID+string(".lo")).c_str()<<endl; 
	}
	else {
	  to<<buffer<<endl
	    <<"\t"<<(fileID+string(".lo"))<<" \\"<<endl; 
	}
      }
      else to<<buffer<<endl;
    }
  }
  from.close();
  to.close();

  //copy back
  Copy(makefilename+string(".tmp"),makefilename);
}



bool Phase_Space_Generator::GetLibList(std::list<std::string>* liblist)
{
  string chlname   = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pathID + string("/fsrchannels");
  string chmapname = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+path   + string("/fsrchannels.map");
 
 ifstream chlist;
  chlist.open(chlname.c_str());
  if (!IsFile(chmapname)) {
    msg_Error()<<"Error in Phase_Space_Generator:"
		       <<chmapname<<" not found."<<endl;
    return 0;
  }

  char buffer[buffersize];
  string libname;
  for(;chlist;) {
    chlist.getline(buffer,buffersize);    
    libname = string(buffer);
    if (chlist && libname[0]!='%') {
      ifstream chmap(chmapname.c_str());
      if (!RSearch(chmap,libname) || libname.find(": ")==string::npos) {
	msg_Error()<<"Error in Phase_Space_Generator:"
			   <<"Mapping for "<<libname<<" not found."<<endl;	
	return 0;
      }
      chmap.close();

      if (libname[0]!='%') {
	int pos = libname.find(string(": "));
	libname = libname.substr(pos+2);
      
	liblist->push_back(libname);
      }
    }
  }
  chlist.close();
  return 1;
}
/*bool Phase_Space_Generator::LoadChannels(ATOOLS::Flavour * fl,Multi_Channel * Ch,Process_Base * proc)
{
  string chlname   = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pathID + string("/fsrchannels");
  string chmapname = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+path   + string("/fsrchannels.map");
 
 ifstream chlist;
  chlist.open(chlname.c_str());
  if (!IsFile(chmapname)) {
    msg_Error()<<"Error in Phase_Space_Generator:"
		       <<chmapname<<" not found."<<endl;
    return 0;
  }

  char buffer[buffersize];
  string libname;
  for(;chlist;) {
    chlist.getline(buffer,buffersize);    
    libname = string(buffer);
    if (chlist && libname[0]!='%') {
      ifstream chmap(chmapname.c_str());
      if (!RSearch(chmap,libname) || libname.find(string(": "))==string::npos) {
	msg_Error()<<"Error in Phase_Space_Generator:"
			   <<"Mapping for "<<libname<<" not found."<<endl;	
	return 0;
      }
      chmap.close();

      if (libname[0]!='%') {
	int pos = libname.find(string(": "));
	libname = libname.substr(pos+2);
      
	Single_Channel * sc = SetChannel(nin,nout,fl,libname);
	if (sc==0) {
	  msg_Error()<<"Phase_Space_Generator:"
			     <<"Channels are not compiled and linked yet."<<endl
			     <<"Type 'make install' and run again."<<endl;
	  return 0;
	}
	else {
	  sc->SetName(pID+string("--")+libname);
	  Ch->Add(sc);
	}
      }
    }
  }
  chlist.close();
  return 1;
  }*/

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

int  Phase_Space_Generator::RSearch(string file,string &search)
{  

  ifstream from;
  //search name  
  from.open(file.c_str());

  char buffer[buffersize];

  for(;from;) {
    from.getline(buffer,buffersize);    
    if (string(buffer).find(string(search))!=string::npos) {
      search = string(buffer);
      from.close();
      return 1;
    }
  }
  from.close();
  return 0;
}

void Phase_Space_Generator::Copy(string sfrom,string sto)
{
  ifstream from;
  ofstream to;
  
  from.open(sfrom.c_str());
  to.open(sto.c_str()); 

  char ch;
  while(from.get(ch)) to.put(ch);
  from.close();
  to.close();  

  remove(sfrom.c_str());
}



