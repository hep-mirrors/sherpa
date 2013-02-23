#include "SHERPA/Tools/Output_HepEvt.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_HepEvt::Output_HepEvt(const Output_Arguments &args):
  Output_Base("HEPEVT")
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".hepevt";
  int precision       = args.p_reader->GetValue<int>("OUTPUT_PRECISION",12);
#ifdef USING__GZIP
  m_ext += ".gz";
#endif
  m_outstream.open((m_basename+m_ext).c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
  m_outstream.precision(precision);
}

Output_HepEvt::~Output_HepEvt()
{
  m_outstream.close();
}

void Output_HepEvt::ChangeFile()
{
  m_outstream.close();
  std::string newname(m_basename+m_ext);
  for (size_t i(0);FileExists(newname);
       newname=m_basename+"."+ToString(++i)+m_ext);
  m_outstream.open(newname.c_str());
}

void Output_HepEvt::Output(Blob_List* blobs, const double weight) 
{
  m_hepevt.Sherpa2HepEvt(blobs);
  m_hepevt.SetWeight(weight);
  m_hepevt.WriteFullHepEvt(m_outstream,m_hepevt.Nhep());
}

DECLARE_GETTER(Output_HepEvt,"HEPEVT",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepEvt>::
operator()(const Output_Arguments &args) const
{
  return new Output_HepEvt(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HepEvt>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"HEPEVT output";
}

