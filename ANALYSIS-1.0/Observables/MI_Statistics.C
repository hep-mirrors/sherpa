#include "MI_Statistics.H"

using namespace ANALYSIS;

DECLARE_GETTER(MI_Statistics_Getter,"MIStats",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *const 
MI_Statistics_Getter::operator()(const String_Matrix &parameters) const
{
  std::string listname="Analysed";
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new MI_Statistics(listname);
}

void MI_Statistics_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

#include "Primitive_Analysis.H"

#include <fstream>

using namespace ATOOLS;

MI_Statistics::MI_Statistics(const std::string & listname, int type):
  Primitive_Observable_Base(type,0,50,50,NULL) 
{
  m_name  = "MI_Statistics.dat";
  m_type  = type;
  m_listname    = listname;
  m_splitt_flag = false;
}

void MI_Statistics::Evaluate(const Blob_List &  blobs,double weight,int ncount)
{
  unsigned int number=0;
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision) {
      ++number;
    }
  }
  p_histo->Insert(number,weight,ncount);
}

Primitive_Observable_Base * MI_Statistics::Copy() const 
{
  return new MI_Statistics(m_listname,m_type);
}

