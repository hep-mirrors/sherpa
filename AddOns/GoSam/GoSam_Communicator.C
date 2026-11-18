#include "GoSam_Communicator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include <fstream>
#include <sstream>

using namespace GoSam;
using namespace ATOOLS;
using namespace std;

GoSam_Communicator::GoSam_Communicator(string name) :
  m_filestatus(0), m_name(name)
{
  ifstream ifile;
  ifile.open(m_name.c_str());
  if (ifile) {
    ifile.close();
    m_filestatus=1;
    return;
  }
  ofstream ofile;
  ofile.open(m_name.c_str(),ios::trunc);
  ofile<<"# "<<m_name<<endl;
  ofile<<"# Created by Sherpa-"<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<endl<<endl;
  ofile.close();
}


GoSam_Communicator::~GoSam_Communicator()
{}

bool GoSam_Communicator::GetPLine(std::ifstream& file,std::string& s1,std::string&s2)
{
  s1="";s2="";
  string buffer;
  if (!file) return 0;
  getline(file,buffer);
  buffer=buffer.substr(0,buffer.find("#"));
  size_t pos=buffer.find("|");
  s1=buffer.substr(0,pos);
  if (pos!=string::npos) s2=buffer.substr(pos+1);
  while (s1.length()>0 && s1[0]==' ') s1.erase(0,1);
  while (s2.length()>0 && s2[0]==' ') s2.erase(0,1);
  return 1;
}

string GoSam_Communicator::StripProcessName(const string &s)
{
  istringstream iss(s);
  string tok, res;
  while (iss>>tok) {                // split at whitespace
    res+=(res.empty()?"":" ")+tok;  // join with exactly 1 whitespace
  }
  return res;
}


void GoSam_Communicator::AddParameter(string param)
{
  ofstream ofile;
  ofile.open(m_name.c_str(),ios::app);
  ofile<<param<<endl;
  ofile.close();
}

int  GoSam_Communicator::CheckParameterStatus()
{
  ifstream ifile;
  ifile.open(m_name.c_str());
  string s1,s2;
  int status(1);
  while (GetPLine(ifile,s1,s2)) {
    if (s1.length()>0&&s2.length()>0) 
      if (s1.find("->")==string::npos)
	if (s2.find("OK")==string::npos)
	  {
	    cout<<endl<<"Warning: OLE returned \""<<s2<<"\" for parameter "<<s1<<endl;
	    status=-1;
	  }
  }
  return status;
}

int GoSam_Communicator::CheckProcess(int nin,int nout,const Flavour_Vector& flavs)
{
  string pstr("");
  for (int i=0;i<nin;i++) pstr+=ToString((long int)flavs[i])+" ";
  pstr+="->";
  for (int i=nin;i<nin+nout;i++) pstr+=" "+ToString((long int)flavs[i]);
  ifstream ifile;
  ifile.open(m_name.c_str());
  string s1,s2;
  int status(-1);
  while (GetPLine(ifile,s1,s2)) {
    if (StripProcessName(s1)==pstr) {
      if (s2.length()>0) {
	s2=s2.substr(0,s2.find(" "));
	status=ToType<int>(s2);
	if (ToString(status)!=s2) {
	    cout<<endl<<"Warning: OLE returned \""<<s2<<"\" for process "<<s1<<endl;
	    status=-2;
	}
	break;
      }
      else {
	status=0;
	break;
      }
    }
  }
  ifile.close();
  return status;
}

void GoSam_Communicator::AddProcess(int nin,int nout,const Flavour_Vector& flavs)
{
  string pstr("");
  for (int i=0;i<nin;i++) pstr+=ToString((long int)flavs[i])+" ";
  pstr+="->";
  for (int i=nin;i<nin+nout;i++) pstr+=" "+ToString((long int)flavs[i]);
  AddParameter(pstr);
}

int GoSam_Communicator::GetID(int nin,int nout,const Flavour_Vector& flavs,int n)
{
  string pstr("");
  for (int i=0;i<nin;i++) pstr+=ToString((long int)flavs[i])+" ";
  pstr+="->";
  for (int i=nin;i<nin+nout;i++) pstr+=" "+ToString((long int)flavs[i]);
  ifstream ifile;
  ifile.open(m_name.c_str());
  string s1,s2;
  int id(-1);
  while (GetPLine(ifile,s1,s2)) {
    if (StripProcessName(s1)==pstr) {
      if (s2.length()>0) {
	string hstr=s2.substr(0,s2.find(" "));
	int nb=ToType<int>(hstr);
	if (n<nb) {
	  for(int i=0;i<n+1;i++) {
	    s2=s2.substr(s2.find(" ")+1);
	    while (s2.length()>0 && s2[0]==' ') s2.erase(0,1);
	  }
	  s2=s2.substr(0,s2.find(" "));
	  id=ToType<int>(s2);
	}
      }    
      break;
    }
  }
  ifile.close();
  return id;
}
