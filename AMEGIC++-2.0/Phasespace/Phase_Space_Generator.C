#include <stdlib.h>
#include "Phase_Space_Generator.H"
#include "Channel_Generator.H"
#include "Process_Base.H"
#include "Run_Parameter.H"
#include "Message.H"

#include "Running_AlphaQED.H"

#include <stdio.h>


using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS; 
using namespace std;

Phase_Space_Generator::Phase_Space_Generator(int _nin,int _nout) : nin(_nin), nout(_nout) 
{
  m_mode=1;
}

bool Phase_Space_Generator::Construct(Multi_Channel * Ch,string _pathID,string _pID,
				      ATOOLS::Flavour* fl,Process_Base * proc)
{ 
  //  pathID = _pathID;
  pathID = _pathID + string("/") + _pID;
  pID    = string("P")+_pID;

  int ngraph = proc->NumberOfDiagrams();
  if (ngraph<=0) {
    msg.Error()<<"Error in Phase_Space_Generator::Construct for "<<proc->Name()<<endl;
    abort();
  }
  bool newchannels = 0;
  int extrachannel = 0;
  for (int i=0;i<ngraph;) {
    if (proc->IsFreeOfFourVertex(proc->Diagram(i))) {
      if (CheckForOldChannels(i,fl,Ch)) { i++; continue; } 
      Channel_Generator * cg = new Channel_Generator(nin,nout,fl,proc->Diagram(i));

      int  rannumber;
      // Channel construction starts here.
      if (extrachannel) sprintf(procname,"C%ia",i);
                   else sprintf(procname,"C%i",i);
      cg->SetName(pID+string("--")+string(procname));
      rannumber    = cg->MakeChannel(extrachannel,i,pathID,pID);
      if (newchannels == 0) {
	if ( (!MakeHeader(headername,procname,rannumber)) ||
	     (!MakeCfile(cfilename,procname,rannumber,Ch->Number()))  )
	{
	  msg.Error()<<"Error in Phase_Space_Generator::Construct for "<<proc->Name()<<endl
		     <<"  Header or C-File already existent."<<endl;
	  abort();
	} 
      }
      else {
	AddToHeader(headername,procname,1);
	AddToCfile(cfilename,procname,Ch->Number(),1);
      }
      if (rannumber>0) {
	sprintf(filename,"C%i",i);
	string makefilename = string("Process/")+pathID+string("/Makefile");
	AddToMakefile(makefilename,pathID,string("P"));	
	AddToMakefile(makefilename,pathID,filename);	
	AddToSetChannel();
	Ch->Add(cg);
      }
      newchannels = 1;
      if (!extrachannel) i++;
    }
    else i++;
  }
  return newchannels;
}

bool Phase_Space_Generator::MakeHeader(string &headername,char* name,int rannum)
{
  if (IsFile(headername)) return 0;

  //  system((string("cd Process;ln -s ")+pathID+string("/")+pID+string(".H ")).c_str());

  msg.Debugging()<<"Creating the headerfile "<<headername<<endl;
  ofstream header;

  header.open(headername.c_str());

  header<<"//Header file for process "<<pID.c_str()<<endl<<endl
	<<"#ifndef "<<pID.c_str()<<"_on"<<endl
	<<"#define "<<pID.c_str()<<"_on"<<endl
	<<"#include "<<'"'<<"Single_Channel.H"<<'"'<<endl<<endl
	<<"using namespace PHASIC;"<<endl<<endl;

  header<<"extern "<<'"'<<"C"<<'"'<<" Single_Channel * Getter_"<<pID
	<<"(int nin,int nout,ATOOLS::Flavour* fl, int chn);"<<endl<<endl;

  header<<"namespace PHASIC {"<<endl
	<<"  class "<<pID.c_str()<<" : public Single_Channel {"<<endl
	<<"    int     chnumber;"<<endl;

  //actual Channel
  header<<"    void   "<<name<<"Momenta(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *);"<<endl
	<<"    double "<<name<<"Weight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);"<<endl
	<<"    int    "<<name<<"Resonances(ATOOLS::Flavour*&);"<<endl
	<<"    void   "<<name<<"ISRtype(int &,double &,double &);"<<endl
	<<"  public:"<<endl
	<<"    "<<pID.c_str()<<"(int nin,int nout,ATOOLS::Flavour* fl,int _chn)"<<endl
	<<"       : Single_Channel(nin,nout,fl), chnumber(_chn)"<<endl
	<<"    { name = std::string(\""<<pID.c_str()<<"\"); };"<<endl
	<<"    void   GenerateWeight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);"<<endl
	<<"    void   GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *);"<<endl
	<<"    int    CountResonances(ATOOLS::Flavour*&);"<<endl
	<<"    void   ISRInfo(int &,double &,double &);"<<endl
	<<"    int    ChNumber()           { return chnumber; };"<<endl
	<<"    void   SetChNumber(int _ch) { chnumber = _ch; };"<<endl
	<<"  };"<<endl
	<<"}"<<endl<<endl
	<<"#endif"<<endl;

  header.close();
  msg.Debugging()<<headername.c_str()<<" saved."<<endl;
  return 1;
}

bool Phase_Space_Generator::MakeCfile(string &cfilename,char* name,
				      int rannum,int chnumber)
{
  if (IsFile(cfilename)) return 0;

  //  system((string("cd Process;ln -s ")+pathID+string("/")+pID+string(".C ")).c_str());

  msg.Debugging()<<"Creating the c file "<<cfilename<<endl;
  ofstream cfile;

  cfile.open(cfilename.c_str());

  cfile<<"//Process "<<pID<<endl<<endl
  //     <<"#include "<<'"'<<pID<<".H"<<'"'<<endl
       <<"#include "<<'"'<<"P.H"<<'"'<<endl
       <<"#include "<<'"'<<"Random.H"<<'"'<<endl<<endl
       <<"using namespace PHASIC;"<<endl<<endl;

  cfile<<"extern "<<'"'<<"C"<<'"'<<" Single_Channel * Getter_"<<pID
       <<"(int nin,int nout,ATOOLS::Flavour* fl, int chn) {"<<endl
       <<"  return new "<<pID<<"(nin,nout,fl,chn);"<<endl
       <<"}"<<endl<<endl;


  cfile<<"void "<<pID<<"::GeneratePoint("
       <<"ATOOLS::Vec4D * p,ATOOLS::Cut_Data * cuts,double* ran)"<<endl 
       <<"{"<<endl  
       <<"  switch (chnumber) {"<<endl
       <<"    case "<<chnumber<<": "<<name<<"Momenta(p,cuts,ran);break;"<<endl
       <<"    default:std::cerr<<"<<'"'
       <<"Channel Number"<<'"'<<"<<chnumber<<"<<'"'<<" not found !"<<'"'<<"<<std::endl;"<<endl
       <<"  }"<<endl
       <<"}"<<endl<<endl
       <<"void "<<pID<<"::GenerateWeight("
       <<"ATOOLS::Vec4D * p,ATOOLS::Cut_Data * cuts)"<<endl 
       <<"{"<<endl  
       <<"  switch (chnumber) {"<<endl
       <<"    case "<<chnumber<<": weight = "<<name<<"Weight(p,cuts);break;"<<endl
       <<"    default:std::cerr<<"<<'"'
       <<"Channel Number"<<'"'<<"<<chnumber<<"<<'"'<<" not found !"<<'"'<<"<<std::endl;"<<endl
       <<"  }"<<endl
       <<"}"<<endl<<endl
       <<"int "<<pID<<"::CountResonances(ATOOLS::Flavour*& fl_res)"<<endl 
       <<"{"<<endl  
       <<"  switch (chnumber) {"<<endl
       <<"    case "<<chnumber<<": return "<<name<<"Resonances(fl_res);"<<endl
       <<"    default:std::cerr<<"<<'"'
       <<"Channel Number"<<'"'<<"<<chnumber<<"<<'"'<<" not found !"<<'"'<<"<<std::endl;"<<endl
       <<"  }"<<endl
       <<"  return 0;"<<endl
       <<"}"<<endl<<endl
       <<"void "<<pID<<"::ISRInfo(int & type,double & mass,double & width)"<<endl 
       <<"{"<<endl
       <<"  switch (chnumber) {"<<endl
       <<"    case "<<chnumber<<": "<<name<<"ISRtype(type,mass,width); return;"<<endl
       <<"    default:std::cerr<<"<<'"'
       <<"Channel Number"<<'"'<<"<<chnumber<<"<<'"'<<" not found !"<<'"'<<"<<std::endl;"<<endl
       <<"  }"<<endl
       <<"  return;"<<endl
       <<"}"<<endl<<endl;
  
  cfile.close();
  msg.Debugging()<<cfilename.c_str()<<" saved."<<endl;
  return 1;
}

int Phase_Space_Generator::AddToHeader(string &headername,char* name,bool construct)
{
  if (IsFile(headername)==0) return 0;

  ifstream from;
  //search name  
  from.open(headername.c_str());
  if (Search(from,string(name))) return 1;
  from.close();
  
  if (!construct) return 0;
  //add Channel
  ofstream to;  


  ifstream from2;
  from2.open(headername.c_str()); 
  to.open((headername+string(".tmp")).c_str());
  
  char buffer[buffersize];
  
  for(;from2;) {
    from2.getline(buffer,buffersize);
    if (string(buffer).find(string("public:"))!=-1) {
      to<<"    void   "<<name<<"Momenta(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *);"<<endl
	<<"    double "<<name<<"Weight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);"<<endl
	<<"    int    "<<name<<"Resonances(ATOOLS::Flavour*&);"<<endl
	<<"    void   "<<name<<"ISRtype(int &,double &,double &);"<<endl;
    }
    to<<buffer<<endl;
  }
  from2.close();
  to.close();
  
  //copy back
  Copy((headername+string(".tmp")).c_str(),headername.c_str());
  
  return 0;
}

int Phase_Space_Generator::AddToCfile(string &cfilename,char* name,int chnumber,bool construct)
{
  if (IsFile(cfilename)==0) return 0;

  //search name
  ifstream from;
  //search name  
  from.open(cfilename.c_str());
  if (Search(from,string(name))) return 1;
  from.close();
  
  if (!construct) return 0;
  //add Channel
  ofstream to;  

  ifstream from2;
  from2.open(cfilename.c_str()); 
  to.open((cfilename+string(".tmp")).c_str());
  
  char buffer[buffersize];
  int mom = 1;

  for(;from2;) {
    from2.getline(buffer,buffersize);

    if (string(buffer).find(string("GenerateWeight"))!=-1) mom = 0;   
    if (string(buffer).find(string("GeneratePoint"))!=-1) mom = 1;
    if (string(buffer).find(string("CountResonances"))!=-1) mom = 2;   
    if (string(buffer).find(string("ISRtype"))!=-1) mom = 3;   

    if (string(buffer).find(string("default"))!=-1) {
      switch (mom) {
      case 0:
        to<<"    case "<<chnumber<<": weight = "<<name<<"Weight(p,cuts);break;"<<endl;
	break;
      case 1:
	to<<"    case "<<chnumber<<": "<<name<<"Momenta(p,cuts,ran);break;"<<endl;
	break;
      case 2:
	to<<"    case "<<chnumber<<": return "<<name<<"Resonances(fl_res);"<<endl;
	break;
      case 3:
	to<<"    case "<<chnumber<<": "<<name<<"ISRtype(type,mass,width); return;"<<endl;
	break;
      }

    }
    to<<buffer<<endl;
  }
  from2.close();
  to.close();

  //copy back
  Copy((cfilename+string(".tmp")).c_str(),cfilename.c_str());

  return 0;
}


void Phase_Space_Generator::AddToMakefileAM(string makefilename,string pathID,string fileID)
{
  cout<<"Phase_Space_Generator::AddToMakefileAM("<<makefilename<<","<<pathID<<","<<fileID<<")"<<endl;

  unsigned int hit=pathID.find("/");
  string base=pathID.substr(0,hit);
  string subdirname=pathID.substr(hit+1);

  ifstream from(makefilename.c_str());
  ofstream to((makefilename+string(".tmp")).c_str());  

  string buffer;
  string key=string("libProc_"+subdirname+"_la_SOURCES");
  string key_h=string("noinst_HEADERS");
  for (;from;) {
    getline(from,buffer);
    if (buffer.find(key_h)!=string::npos) {
      if (!(buffer.find("P.H")!=string::npos)) {
	to<<buffer<<" P.H"<<endl;
      }
      else {
	to<<buffer<<endl;
      }
    }
    else {
      to<<buffer<<endl;
    }
    if (buffer.find(key)!=string::npos) {
      to<<"\t"<<fileID<<".C"<<'\\'<<endl;
    }
  }
  from.close();
  to.close();
  
  string mv=string("mv ")+makefilename+".tmp "+makefilename;
  system(mv.c_str());
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
      if (string(buffer).find(string("\\"))==-1) {
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
	if (string(buffer).find(string("\\"))==-1) {
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

void Phase_Space_Generator::AddToSetChannel()
{
  // only include in Set_Channel.C if neccessary
  if (m_mode==1) return;

  ifstream from;
  ofstream to;

  from.open(string("Process/Set_Channel.C").c_str());
  to.open(string("Process/Set_Channel.C.tmp").c_str());

  int hit = 0;

  string buffer;

  //include into first line

  //  to<<"#include "<<'"'<<(pID+string(".H")).c_str()<<'"'<<endl;  
  to<<"#include "<<'"'<<(pathID+string("/P.H")).c_str()<<'"'<<endl;  

  for(;from;) {
    getline(from,buffer);

    if (buffer.find(pID)!=-1) break;

    if (buffer.find(string("return 0"))!=-1 || buffer.find(string("libname"))!=-1) {
      hit = 1;
      to<<"#ifdef "<<pID<<"_on"<<endl
	<<"  if (pID==string("<<'"'<<pID<<'"'<<")) "
	<<"return (new "<<pID<<"(nin,nout,fl,chn));"<<endl
	<<"#endif"<<endl;
    }
    to<<buffer<<endl;
  }
  from.close();
  to.close();

  if (hit) 
    Copy(string("Process/Set_Channel.C.tmp"),string("Process/Set_Channel.C"));
  else 
    remove("Process/Set_Channel.C.tmp");
}




bool Phase_Space_Generator::CheckForOldChannels(int n,ATOOLS::Flavour * fl,Multi_Channel * Ch)
{
  int  hit;
  bool oldch = 0;

  // Loop for Channel construction starts here.
  sprintf(procname,"C%i",n);
  
  headername = string("Process/")+pathID + string("/P.H");
  cfilename  = string("Process/")+pathID + string("/P.C");
  
  msg.Debugging()<<"   Checking for : "<<headername<<", "<<endl;
  hit  = 0;
  hit += AddToHeader(headername,procname,0); 
  msg.Debugging()<<"                  "<<cfilename<<endl;
  hit += AddToCfile (cfilename,procname,Ch->Number(),0); 
  
  if (hit==2) {
    Single_Channel * sc = SetChannel(nin,nout,fl,Ch->Number(),pID);
    if (sc==0) {
      ATOOLS::msg.Error()<<"Phase_Space_Generator:"
			    <<"Channels are not compiled and linked yet."<<endl
			    <<"Type 'make install' and run again."<<endl;
      return 0;
    }
    else {
      sc->SetName(pID+string("--")+string(procname));
      Ch->Add(sc);
    }
    oldch = 1;
  }
  if (hit==1) {
    ATOOLS::msg.Error()<<"Error in Phase_Space_Generator:"
			  <<"C- or H-file not found. Abort."<<endl;
    return 0;
  }

  if (oldch==1) return 1;
  return 0;
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
    if (string(buffer).find(string(search))!=-1) return 1;
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
    if (string(buffer).find(string(search))!=-1) {
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



