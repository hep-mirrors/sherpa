#include "MI_Statistics.H"
#include "Primitive_Analysis.H"

#include <fstream>

using namespace ANALYSIS;
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

