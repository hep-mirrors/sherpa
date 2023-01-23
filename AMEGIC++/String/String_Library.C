#include "AMEGIC++/String/String_Library.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <fstream>
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;


String_Library::String_Library(int mode):m_mode(mode)
{

}

void String_Library::UpdateConfigure(std::string pathID)
{
  msg_Debugging()<<"String_Library::UpdateConfigure("<<pathID<<") called :"<<std::endl;
  string cnf("/CMakeLists.txt");
  unsigned int hit=pathID.find("/");
  string base=pathID.substr(0,hit);
  string subdirname=pathID.substr(hit+1);
  string name=rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+base+cnf;
  if (!IsFile(name)) {
    msg_Tracking()<<"   file "<<name<<" does not exist, create it."<<endl;
    ofstream file(name.c_str());
    file<<"add_subdirectory("<<subdirname<<")"<<endl;
  } 
  else {
    ifstream from(name.c_str());
    ofstream to((name+string(".tmp")).c_str());
    string buffer;
    bool present=false;
    for (;from;) {
      getline(from,buffer);
      if (buffer.find(std::string("add_subdirectory(")+subdirname+")")!=string::npos) present=true;
      to<<buffer<<endl;
    }
    from.close();
    if (!present) to<<"add_subdirectory("<<subdirname<<")"<<endl;
    to.close();

    Move(name+".tmp",name);
  }


}

void String_Library::AddToCMakefile(string makefilename,string pathID,string fileID)
{
  msg_Debugging()<<"String_Library::AddToCMakefile("<<makefilename<<","<<pathID<<","<<fileID<<")"<<endl;

  unsigned int hit=pathID.find("/");
  string base=pathID.substr(0,hit);
  string subdirname=pathID.substr(hit+1);

  if (!IsFile(makefilename)) {
    ofstream file(makefilename.c_str());
    file<<"set(libProc_"<<subdirname<<"_la_SOURCES "<<endl;
    file<<fileID<<".C"<<endl;
    file<<")"<<endl;
    file<<"add_library(Proc_"<<subdirname<<"SHARED" << "${libProc_"<<subdirname<<"_la_SOURCES})"<<endl;
    file<<"amegic_handle_shared_library(Proc_"<<subdirname<<")"<<endl;
    file<<"install(TARGETS Proc_"<<subdirname<<" DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})"<<endl;
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
  }
}



void String_Library::InitMakefile(string pathID)
{
  UpdateConfigure(pathID);
  return;
}

void String_Library::Replace(string& buffer,const string& search,const string& replace)
{
  int minpos=0;
  while (SingleReplace(buffer,search,replace,minpos));
}

int String_Library::SingleReplace(string& buffer, const string& search,const string& replace,int& minpos)
{
  int pos= buffer.find(search,minpos);
  if (pos==-1) return 0;
  minpos=pos+replace.length();
  buffer = buffer.substr(0,pos)+replace+buffer.substr(pos+search.length());
  return 1;
}

void String_Library::Copy(string sfrom,string sto)
{
  ifstream from;
  ofstream to;
  
  from.open(sfrom.c_str());
  to.open(sto.c_str()); 

  char ch {0};
  while(from.get(ch)) to.put(ch);
  from.close();
  to.close();  

  //kill tmp
  remove(sfrom.c_str());
}

int String_Library::IsFile(string &filename)
{
  ifstream from;
  from.open(filename.c_str());
  if (from) return 1;
  return 0;
}

