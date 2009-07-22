#include "EXTRA_XS/NLO/LH_OLE_Communicator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include <fstream>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace std;


LH_OLE_Communicator::LH_OLE_Communicator(string name) :
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
  ofile<<"# "<<m_name<<endl<<endl;
  ofile.close();
}


LH_OLE_Communicator::~LH_OLE_Communicator()
{}

void LH_OLE_Communicator::AddParameter(string param)
{
  ofstream ofile;
  ofile.open(m_name.c_str(),ios::app);
  ofile<<param<<endl;
  ofile.close();
}

int  LH_OLE_Communicator::CheckParameterStatus()
{
  ifstream ifile;
  ifile.open(m_name.c_str());
  string buffer;
  int status(1);
  while (ifile) {
    getline(ifile,buffer);
    if (buffer.find("| -")!=string::npos) if (buffer[0]!='#')
      {
	cout<<endl<<"Warning: OLE returned "<<buffer.substr(buffer.find("| -")+2)<<endl;
	status=-1;
      }
  }
  return status;
}

int LH_OLE_Communicator::CheckProcess(int nin,int nout,const Flavour_Vector& flavs)
{
  string pstr(ToString(nin)+" -> "+ToString(nout));
  for (int i=0;i<nin+nout;i++) pstr+=" "+ToString((long int)flavs[i]);
  ifstream ifile;
  ifile.open(m_name.c_str());
  string buffer;
  int status(-1);
  while (ifile) {
    getline(ifile,buffer);
    if (buffer.find(pstr)!=string::npos) {
      if (buffer.find("|")!=string::npos) {
	buffer=buffer.substr(buffer.find("|")+2);
	buffer=buffer.substr(0,buffer.find(" "));
	status=ToType<int>(buffer);
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

void LH_OLE_Communicator::AddProcess(int nin,int nout,const Flavour_Vector& flavs)
{
  string pstr(ToString(nin)+" -> "+ToString(nout));
  for (int i=0;i<nin+nout;i++) pstr+=" "+ToString((long int)flavs[i]);
  AddParameter(pstr);
}

int LH_OLE_Communicator::GetID(int nin,int nout,const Flavour_Vector& flavs,int n)
{
  string pstr(ToString(nin)+" -> "+ToString(nout));
  for (int i=0;i<nin+nout;i++) pstr+=" "+ToString((long int)flavs[i]);
  ifstream ifile;
  ifile.open(m_name.c_str());
  string buffer;
  int id(-1);
  while (ifile) {
    getline(ifile,buffer);
    if (buffer.find(pstr)!=string::npos) {
      if (buffer.find("|")!=string::npos) {
	buffer=buffer.substr(buffer.find("|")+2);
	for(int i=0;i<n+1;i++) buffer=buffer.substr(buffer.find(" ")+1);
	buffer=buffer.substr(0,buffer.find(" "));
	id=ToType<int>(buffer);
      }    
      break;
    }
  }
  ifile.close();
  return id;

}
