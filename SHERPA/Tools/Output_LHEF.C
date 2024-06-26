#include "SHERPA/Tools/Output_LHEF.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/My_MPI.H"

#include <iomanip>
#include <stdio.h>
#include <cassert>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_LHEF::Output_LHEF(const Output_Arguments &args):
  Output_Base{ "LHEF" },
  m_xs(1.0), m_xserr(1.0),
  m_max(1.0)
{
  Settings& s = Settings::GetMainSettings();
  RegisterDefaults();
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".lhe";
  const int precision = s["EVENT_OUTPUT_PRECISION"].Get<int>();
  m_bntp = s["LHEF_BNTP"].Get<int>();
#ifdef USING__GZIP
  m_ext += ".gz";
#endif
#ifdef USING__MPI
  if (mpi->Size()>1)
    m_basename+="_"+rpa->gen.Variable("RNG_SEED");
#endif
  m_outstream.open((m_basename+m_ext).c_str());
  if (!m_outstream.good())
    THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
  m_outstream.precision(precision);
}

void Output_LHEF::RegisterDefaults() const
{
  Settings& s = Settings::GetMainSettings();

  s["LHEF_BNTP"].SetDefault(0);
  s["LHEF_IDWTUP"].SetDefault(0);

  int PDFSUP1_default{ rpa->gen.PDF(0)?rpa->gen.PDF(0)->LHEFNumber():-1 };
  int PDFSUP2_default{ rpa->gen.PDF(1)?rpa->gen.PDF(1)->LHEFNumber():-1 };
  s["LHEF_PDF_NUMBER_1"].SetDefault(PDFSUP1_default);
  s["LHEF_PDF_NUMBER_2"].SetDefault(PDFSUP2_default);
}

Output_LHEF::~Output_LHEF()
{
  m_outstream.close();
}

void Output_LHEF::ChangeFile()
{
  m_outstream.close();
  std::string newname(m_basename+m_ext);
  for (size_t i(0);FileExists(newname);
       newname=m_basename+"."+ToString(++i)+m_ext);
  m_outstream.open(newname.c_str());
}

void Output_LHEF::Header()
{
  Settings& s = Settings::GetMainSettings();

  const auto path = s.GetPath();
  const auto files = s.GetConfigFiles();

  m_outstream<<"<LesHouchesEvents version=\"1.0\">"<<std::endl;
  m_outstream<<"<header>"<<std::endl;
  m_outstream<<"<!-- "<<std::endl;
  m_outstream<<"# created by SHERPA "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION
             <<endl;
  for (const auto& file : files) {
    if (file.rfind("Decaydata.yaml") != std::string::npos) continue;
    My_In_File infile{ path, file };
    if (infile.Open()) {
      m_outstream<<"# Run data extracted from : "<<file<<std::endl;
      m_outstream<<"--> "<<std::endl;
      m_outstream<<"<SHRunCard> "<<std::endl;
      m_outstream<<infile.FileContent();
      m_outstream<<"</SHRunCard> "<<std::endl;
    }
  }
  m_outstream<<"</header>"<<std::endl;

  m_outstream<<"<init>"<<std::endl;
  //run info to be dumped here
  Flavour Beam1 = rpa->gen.Beam1();
  Flavour Beam2 = rpa->gen.Beam2();
  int IDBMUP1 = (long int)Beam1;
  int IDBMUP2 = (long int)Beam2;
  double EBMUP1 = rpa->gen.PBeam(0)[0];
  double EBMUP2 = rpa->gen.PBeam(1)[0];

  int IDWTUP{ s["LHEF_IDWTUP"].Get<int>() };
  if (IDWTUP==0) IDWTUP=ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"))==0?1:3;
  int NPRUP = 1;
  int PDFGUP1 = 0;
  int PDFGUP2 = 0;
  const int PDFSUP1{ s["LHEF_PDF_NUMBER_1"].Get<int>() };
  const int PDFSUP2{ s["LHEF_PDF_NUMBER_2"].Get<int>() };

  m_outstream<<std::setprecision(10);
  m_outstream<<std::setw(6)<<IDBMUP1<<" "
	     <<std::setw(6)<<IDBMUP2<<" "
	     <<std::setw(11)<<EBMUP1<<" "
	     <<std::setw(11)<<EBMUP2<<" "
	     <<std::setw(3)<<PDFGUP1<<" "
	     <<std::setw(3)<<PDFGUP2<<" "
	     <<std::setw(6)<<PDFSUP1<<" "
	     <<std::setw(6)<<PDFSUP2<<" "
	     <<std::setw(4)<<IDWTUP<<" "
	     <<std::setw(4)<<NPRUP<<std::endl;

  //process information
  m_outstream<<std::setw(18)<<m_xs<<" "
	     <<std::setw(18)<<m_xserr<<" "
	     <<std::setw(18)<<m_max<<" "
	     <<std::setw(4)<<NPRUP<<std::endl;
  m_outstream<<"</init>"<<std::endl;
}

void Output_LHEF::Output(Blob_List* blobs)
{
  const auto weight(blobs->Weight());
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  m_outstream<<"<event trials='"<<(int)(*sp)["Trials"]->Get<double>();
  if ((*sp)["MC@NLO_KT2_Start"])
    m_outstream<<"' kt2_start='"<<(*sp)["MC@NLO_KT2_Start"]->Get<double>()
	       <<"' kt2_stop='"<<(*sp)["MC@NLO_KT2_Stop"]->Get<double>();
  if ((*sp)["Factorisation_Scale"])
    m_outstream<<"' muf2='"<<(*sp)["Factorisation_Scale"]->Get<double>();
  double mur2=0.0;
  if ((*sp)["Renormalization_Scale"]) {
    mur2=(*sp)["Renormalization_Scale"]->Get<double>();
    m_outstream<<"' mur2='"<<mur2;
  }
  m_outstream<<"'>"<<std::endl;
  Poincare cms(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit){
    if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
      //LHE event information
      int NUP = (*blit)->NInP()+(*blit)->NOutP();
      int IDPRUP = 1;
      double XWGTUP = weight;
      double SCALUP = sqrt((*(*blit))["Factorisation_Scale"]->Get<double>());
      if (SCALUP==0.0) SCALUP=sqrt((*(*blit))["Renormalization_Scale"]->Get<double>());
      if ((*(*blit))["Resummation_Scale"])
	SCALUP=sqrt((*(*blit))["Resummation_Scale"]->Get<double>());
      if (SCALUP==0.0) SCALUP=sqrt((*(*blit))["Renormalization_Scale"]->Get<double>());
      double AQEDUP = -1.;
      double AQCDUP = mur2?(*MODEL::as)(mur2):-1.0;
      m_outstream<<std::setprecision(10);
      m_outstream<<std::setiosflags(std::ios::scientific);
      m_outstream<<std::setw(4)<<NUP<<" "
		 <<std::setw(4)<<IDPRUP<<" "
		 <<std::setw(18)<<XWGTUP<<" "
		 <<std::setw(18)<<SCALUP<<" "
		 <<std::setw(18)<<AQEDUP<<" "
		 <<std::setw(18)<<AQCDUP<<std::endl;
      for (int i=0;i<(*blit)->NInP();i++) {
	Vec4D p((*blit)->InParticle(i)->Momentum());
	if (m_bntp) cms.Boost(p);
	m_outstream<<std::setw(8)<<(long int)(*blit)->InParticle(i)->Flav()<<" -1  0  0 "
		   <<std::setw(4)<<(*blit)->InParticle(i)->GetFlow(1)<<" "
		   <<std::setw(4)<<(*blit)->InParticle(i)->GetFlow(2)<<" "
		   <<std::setw(18)<<p[1]<<" "
		   <<std::setw(18)<<p[2]<<" "
		   <<std::setw(18)<<p[3]<<" "
		   <<std::setw(18)<<p[0]<<" "
		   <<std::setw(18)<<(*blit)->InParticle(i)->FinalMass()<<" "
		   <<" 0  9"<<std::endl;
      }
      for (int i=0;i<(*blit)->NOutP();i++) {
	Vec4D p((*blit)->OutParticle(i)->Momentum());
	if (m_bntp) cms.Boost(p);
	m_outstream<<std::setw(8)<<(long int)(*blit)->OutParticle(i)->Flav()<<"  1  1  2 "
		   <<std::setw(4)<<(*blit)->OutParticle(i)->GetFlow(1)<<" "
		   <<std::setw(4)<<(*blit)->OutParticle(i)->GetFlow(2)<<" "
		   <<std::setw(18)<<p[1]<<" "
		   <<std::setw(18)<<p[2]<<" "
		   <<std::setw(18)<<p[3]<<" "
		   <<std::setw(18)<<p[0]<<" "
		   <<std::setw(18)<<(*blit)->OutParticle(i)->FinalMass()<<" "
		   <<" 0  9"<<std::endl;
      }
      m_outstream<<std::resetiosflags(std::ios::fixed);
    }
  }
  m_outstream<<"</event>"<<std::endl;
}

void Output_LHEF::Footer()
{
  string footer = std::string("</LesHouchesEvents>");
  m_outstream<<footer<<std::endl;
}

void Output_LHEF::SetXS(const ATOOLS::Weights_Map& xs,
			const ATOOLS::Weights_Map& xserr)
{
  m_xs = xs.Nominal();
  m_xserr = xserr.Nominal();
  m_max =1.;
}

DECLARE_GETTER(Output_LHEF,"LHEF",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_LHEF>::
operator()(const Output_Arguments &args) const
{
  return new Output_LHEF(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_LHEF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"LHEF output";
}

