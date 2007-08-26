#include "Data_Collector.H"
#include "Message.H"

using namespace ATOOLS;

ATOOLS::String_BlobDataBase_Map Data_Collector::s_datacontainer;
int Data_Collector::s_counter=0;

Data_Collector global_data_collector;
 
Data_Collector::Data_Collector()
{
  ++s_counter;
}

void Data_Collector::AddData(const std::string name, Blob_Data_Base * data) 
{
  String_BlobDataBase_Map::iterator it=s_datacontainer.find(name);
  if (it==s_datacontainer.end()) {
    s_datacontainer[name]=data;
  }
  else {
    delete it->second;
    it->second=data;
  }
}

ATOOLS::Blob_Data_Base * Data_Collector::operator[](const std::string name) 
{
  ATOOLS::String_BlobDataBase_Map::const_iterator cit=s_datacontainer.find(name);
  if (cit==s_datacontainer.end()) return 0;
  return cit->second;
} 

void Data_Collector::Print() 
{
  msg_Out()<<"Data_Collector:"<<std::endl;
  for (String_BlobDataBase_Map::iterator it=s_datacontainer.begin();
       it!=s_datacontainer.end(); ++it) {
    msg_Out()<<"   * "<<it->first<<" ("<<*(it->second)<<")"<<std::endl;
  }
}

Data_Collector::~Data_Collector()
{
  // here the map should be cleared!!!
  --s_counter;
  if (s_counter>0) return;
  for (String_BlobDataBase_Map::iterator it=s_datacontainer.begin();
       it!=s_datacontainer.end(); ++it) delete it->second;
  s_datacontainer.clear();
}

// ======================================================================

template Process_Info &ATOOLS::Blob_Data_Base::Get<Process_Info>();


std::ostream & ATOOLS::operator<<(std::ostream & s, const Process_Info & wi)
{
  s<<" name="<<wi.texname<<"   xsec="<<wi.xsec<<" +- "<<wi.xsec_err<<" pb ";
  return s;
}

namespace ATOOLS {

template <> Blob_Data<Process_Info>::~Blob_Data() {}

template class Blob_Data<Process_Info>;

}
