#include "Output_HepMC2_Ascii.H"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_Ascii.h"
#include "Exception.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepMC2_Ascii::Output_HepMC2_Ascii(string basename,string ext,
                                         int precision) :
      Output_Base(basename, ext, precision)
{
  m_outstream.close();
  p_ioascii = new HepMC::IO_Ascii((m_basename+ext).c_str());
  if (p_ioascii->rdstate() != fstream::goodbit)
    THROW(fatal_error, "Could not open event file "+m_basename+ext+".")
}

Output_HepMC2_Ascii::~Output_HepMC2_Ascii()
{
  if(p_ioascii) delete p_ioascii; p_ioascii=NULL;
}

void Output_HepMC2_Ascii::Output(Blob_List* blobs, const double weight) 
{
  m_hepmc2.Sherpa2HepMC(blobs, weight);
  p_ioascii->write_event(m_hepmc2.GenEvent());
}

void Output_HepMC2_Ascii::ChangeFile(string number)
{
  string newfilename=m_basename+"."+number+m_ext;
  delete p_ioascii;
  p_ioascii = new HepMC::IO_Ascii(newfilename.c_str());
  if (p_ioascii->rdstate() != fstream::goodbit)
    THROW(fatal_error, "Could not open event file "+newfilename+".")
}
