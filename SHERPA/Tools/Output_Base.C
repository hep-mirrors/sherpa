#include "SHERPA/Tools/Output_Base.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_Base::Output_Base(string basename,string ext,int precision) :
  m_basename(basename), m_ext(ext)
{
#ifdef USING__GZIP
  m_ext += ".gz";
#endif
  m_outstream.open((m_basename+m_ext).c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
  m_outstream.precision(precision);
}

Output_Base::~Output_Base()
{
  m_outstream.close();
}

void Output_Base::ChangeFile(string number)
{
  string newfilename=m_basename+"."+number+m_ext;
  m_outstream.close();
  m_outstream.open(newfilename.c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+newfilename+".")
}
