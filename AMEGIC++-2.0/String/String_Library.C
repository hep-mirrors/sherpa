#include "String_Library.H"
#include <fstream>
#include "Message.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

void String_Library::InitMakefile(string pathID)
{
  string newMakefile = string("Process/")+pathID+string("/Makefile");

  if (IsFile(newMakefile)) return;

  msg.Out()<<newMakefile<<" is not available !"<<endl;
    

  system((string("cd Process;cp Dummy/Makefile ")+pathID).c_str());

  string pID;
  pID=pathID;
  for (short int i=pathID.length()-1;i>=0;i--) {
    if (pathID[i]=='/') {
      pID    = pathID.substr(i+1);
      break;
    }
  }
 

    
  ofstream to;  
  ifstream from;
  
  from.open(newMakefile.c_str()); 
  to.open((newMakefile+string(".tmp")).c_str());
  
  char buffer[buffersize];
  
  for(;from;) {
    from.getline(buffer,buffersize);
    
    string buf = string(buffer);
    
    //change subdir structure   (AS)
    Replace(buf,string("../../.."),string("../../../.."));
 
    //replace libDummy
    Replace(buf,string("libDummy"),string("libProc_")+pID);
    //kill Dummy.H & Dummy.C & Dummy.lo
    Replace(buf,string("Dummy.C"),string(""));
    Replace(buf,string("Dummy.H"),string(""));
    Replace(buf,string("Dummy.lo"),string(""));
    //DEP_FILES clear, include Line kill
    Replace(buf,string("include $(DEPDIR)/Dummy.Plo"),string(""));
    Replace(buf,string("$(DEPDIR)/Dummy.Plo"),string(""));
    //subDir & SOURCES '/Dummy'
    Replace(buf,string("/Dummy"),string("/")+pID);
    //Replace Makefile:
    Replace(buf,string("Makefile:"),string("Makefile.bla:"));
    if (buf.find(string("depcomp"))==-1) to<<buf<<endl;
  }
  from.close();
  to.close();
  //copy back
  Copy(newMakefile+string(".tmp"),newMakefile);
  
  //changing SUBDIRS from Process/Makefile
  //======================================================================
  ofstream to2;  
  ifstream from2;
  
  from2.open("Process/Makefile");
  to2.open("Process/Makefile.tmp");
  
  for(;from2;) {
    from2.getline(buffer,buffersize);
    string buf = string(buffer);
    //replace Dummy in SUBDIRS with DUMMY Proc
    int start=0;
    SingleReplace(buf,string("Dummy"),string("Dummy ")+pathID,start);
    to2<<buf<<endl;
  }
  from2.close();
  to2.close();
  //copy back
  Copy(string("Process/Makefile.tmp"),string("Process/Makefile"));
  
  //adding new library to linker
  //============================
  ofstream to3;  
  ifstream from3;

  from3.open("Makefile");
  to3.open("Makefile.tmp");
  
  for(;from3;) {
    from3.getline(buffer,buffersize);
    string buf = string(buffer);
    int pos = buf.find(string("libProcess"));
    if (pos!=-1) {
      if (buf.find(string("\\"))==-1) {
	//no backslash
	to3<<buffer<<"\\"<<endl;
	to3<<"\t"<<"${exec_prefix}/lib/libProc_"<<pID<<".la"<<endl; 
      }
      else {
	to3<<buffer<<endl;
	to3<<"\t"<<"${exec_prefix}/lib/libProc_"<<pID<<".la \\"<<endl; 
      }
    }
    else to3<<buffer<<endl;
  }
  from3.close();
  to3.close();
  //copy back
  Copy(string("Makefile.tmp"),string("Makefile"));
}

void String_Library::Replace(string& buffer,const string& search,const string& replace)
{
  int minpos=0;
  while (SingleReplace(buffer,search,replace,minpos));
}

int String_Library::SingleReplace(string& buffer,
				  const string& search,const string& replace,int& minpos)
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

  char ch;
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

int String_Library::Search(string file,string search)
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

void String_Library::AddToMakefile(string Makefile,string pathID,string fileID)
{
  if (IsFile(Makefile)==0) {
    cerr<<Makefile.c_str()<<" is not available !"<<endl;
    return;
  }

  if (Search(Makefile,string(fileID)+string(".C"))) return;

  ofstream to;  
  ifstream from;


  from.open(Makefile.c_str()); 
  to.open((Makefile+string(".tmp")).c_str());

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
	to<<buffer<<"\\"<<endl;
	to<<"\t"<<(fileID+string(".C")).c_str()<<endl; 
      }
      else {
	to<<buffer<<endl;
	to<<"\t"<<(fileID+string(".C")).c_str()<<" \\"<<endl; 
      }
    }
    else {
      if (string(buffer).find(lib+string("_la_OBJECTS"))==0) {
	if (string(buffer).find(string("\\"))==-1) {
	  //no backslash
	  to<<buffer<<"\\"<<endl;
	  to<<"\t"<<(fileID+string(".lo")).c_str()<<endl; 
	}
	else {
	  to<<buffer<<endl;
	  to<<"\t"<<(fileID+string(".lo"))<<" \\"<<endl; 
	}
      }
      else to<<buffer<<endl;
    }
  }
  from.close();
  to.close();

  //copy back
  Copy(Makefile+string(".tmp"),Makefile);
}
