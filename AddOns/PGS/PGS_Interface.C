#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"

inline void MakeFortranString
(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) 
    output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) 
    output[j]=(char)input[j];
}

const int nmxhep = HEPEVT_CB_SIZE;

extern "C" {

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevt_;
#define hepevt hepevt_

extern struct {
  int nevsha, istream;
  char shafile[80];
} pgspars_;
#define pgspars pgspars_

  void pgsxxx_(int *mode);

  inline void pgsxxx(int mode) { pgsxxx_(&mode); }

}

struct MD_Info {
  int m_n, m_c, m_mo[2], m_da[2];
  inline MD_Info(const int n=0): m_n(n), m_c(0)
  { m_mo[1]=m_mo[0]=m_da[1]=m_da[0]=0; }
};// end of struct MD_Info

typedef std::map<ATOOLS::Particle*,MD_Info> MotherDaughter_Map;

class PGS_Interface: public SHERPA::Analysis_Interface {
private:

  std::string m_inpath, m_infile, m_outfile;

  void ConvertParticle(ATOOLS::Particle *const cp,const int sc,
		       const int mode,MotherDaughter_Map &mdmap);
  void Convert(ATOOLS::Blob_List *const bl,
	       MotherDaughter_Map &mdmap);

  void CheckParticle(ATOOLS::Particle *const cp,const int sc,
		     const int mode,MotherDaughter_Map &mdmap);
  void Check(ATOOLS::Blob_List *const bl,
	     MotherDaughter_Map &mdmap);

public:

  inline PGS_Interface(const std::string &inpath,
			  const std::string &infile,
			  const std::string &outpath):
    Analysis_Interface("PGS"),
    m_inpath(inpath), m_infile(infile) {}

  bool Init();
  bool Run(ATOOLS::Blob_List *const bl);
  bool Finish();

  void ShowSyntax(const int i);

};// end of class PGS_Interface

using namespace SHERPA;
using namespace ATOOLS;

void PGS_Interface::ShowSyntax(const int i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   BEGIN_PGS {\n\n"
	   <<"     FILE_NAME <filename>\n";
  msg_Out()<<"\n   } END_PGS\n\n"
	   <<"}"<<std::endl;
}

void PGS_Interface::ConvertParticle
(ATOOLS::Particle *const cp,const int sc,
 const int mode,MotherDaughter_Map &mdmap)
{
  Blob *pb(cp->ProductionBlob());
  if (pb && pb->NInP()==1 && mode==0)
    ConvertParticle(pb->InParticle(0),3,0,mdmap);
  if (mdmap.find(cp)==mdmap.end()) {
  int n=hepevt.nhep;
  MD_Info &cmd(mdmap[cp]=MD_Info(n));
  hepevt.jmohep[n][0]=hepevt.jmohep[n][1]=0;
  hepevt.jdahep[n][0]=hepevt.jdahep[n][1]=0;
  hepevt.idhep[n]=cp->Flav().HepEvt();
  for (short int j=1;j<4;++j) 
    hepevt.phep[n][j-1]=cp->Momentum()[j];
  hepevt.phep[n][3]=cp->Momentum()[0];
  hepevt.phep[n][4]=Max(0.0,cp->Momentum().Abs2());
  if (pb) {
    for (short int j=1;j<4;++j)
      hepevt.vhep[n][j-1]=pb->Position()[j];
    hepevt.vhep[n][3]=pb->Position()[0];
    if (pb->NInP()==1) {
      MotherDaughter_Map::iterator mdit(mdmap.find(pb->InParticle(0)));
      if (mdit==mdmap.end()) {
	msg_Error()<<METHOD<<"(): Convert error JMOHEP."<<std::endl;
      }
      else {
	MD_Info &mmd(mdit->second);
	int m(mmd.m_n);
	cmd.m_mo[0]=hepevt.jmohep[n][0]=m+1;
	cmd.m_mo[1]=hepevt.jmohep[n][1]=m+1;
	if (hepevt.jdahep[m][0]==0) {
	  mmd.m_da[0]=hepevt.jdahep[m][0]=n+1;
	  mmd.m_da[1]=hepevt.jdahep[m][1]=n+1;
	}
	else {
	  if (n+1<hepevt.jdahep[m][0]-1 || n+1>hepevt.jdahep[m][1]+1)
	    msg_Error()<<METHOD<<"(): Convert error JDAHEP."<<std::endl;
	  mmd.m_da[0]=hepevt.jdahep[m][0]=Min(n+1,hepevt.jdahep[m][0]);
	  mmd.m_da[1]=hepevt.jdahep[m][1]=Max(n+1,hepevt.jdahep[m][1]);
	}
      }
    }
  }
  else {
    for (short int j=0;j<4;++j) hepevt.vhep[n][j]=0.0;
  }
  hepevt.isthep[n]=sc;
  ++hepevt.nhep;
  }
  Blob *db(cp->DecayBlob());
  if (db && db->NInP()==1 && mode==0)
    for (int i(0);i<db->NOutP();++i) {
      Particle *np(db->OutParticle(i));
      ConvertParticle(np,np->DecayBlob()?3:1,1,mdmap);
    }
}

void PGS_Interface::CheckParticle
(ATOOLS::Particle *const cp,const int sc,
 const int mode,MotherDaughter_Map &mdmap)
{
}

void PGS_Interface::Convert(ATOOLS::Blob_List *const bl,
			       MotherDaughter_Map &mdmap)
{
  hepevt.nhep=0;
  Blob_List bb(bl->Find(btp::Bunch));
  ConvertParticle(bb[0]->InParticle(0),3,-1,mdmap);
  ConvertParticle(bb[1]->InParticle(0),3,-1,mdmap);
  Particle_List pl(bl->ExtractParticles(1));
  for (size_t n(0);n<pl.size();++n)
    if (pl[n]->DecayBlob()==NULL)
      ConvertParticle(pl[n],1,0,mdmap);
}

void PGS_Interface::Check(ATOOLS::Blob_List *const bl,
			     MotherDaughter_Map &mdmap)
{
}

bool PGS_Interface::Init()
{
  Data_Reader reader(" ",";","//");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.SetInputPath(m_inpath);
  std::string infile(m_infile);
  if (infile.find('|')!=std::string::npos)
    infile=infile.substr(0,infile.find('|'));
  reader.SetInputFile(infile);
  reader.AddComment("#");
  reader.SetFileBegin("BEGIN_PGS");
  reader.SetFileEnd("END_PGS");
  m_outfile=reader.GetValue<std::string>("FILE_NAME","test.lhe");
  pgspars.nevsha=rpa->gen.NumberOfEvents();
  MakeFortranString(pgspars.shafile,m_outfile,80);
  msg_Info()<<"\n"<<METHOD<<"(): {"<<std::endl;
  pgsxxx(1);
  msg_Info()<<"}"<<std::endl;
  return true;
}

bool PGS_Interface::Run(ATOOLS::Blob_List *const bl)
{
  if (!bl->FourMomentumConservation())
    msg_Error()<<METHOD<<"(): Four momentum not conserved."<<std::endl;
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  Blob_Data_Base *xs((*sp)["Weight"]);
  if (xs==NULL) THROW(fatal_error,"No weight information");
  MotherDaughter_Map mdmap;
  Convert(bl,mdmap);
// #ifdef USING__HERACMN
//   heracmn.wtx=xs->Get<double>();
// #endif
  pgsxxx(2);
  Check(bl,mdmap);
  return true;
}

bool PGS_Interface::Finish()
{
  msg_Info()<<METHOD<<"(): {\n  Total xs is "
	    <<p_eventhandler->TotalXS()<<" pb +- "
	    <<p_eventhandler->TotalErr()<<" pb.\n";
  pgsxxx(3);
  msg_Info()<<"}"<<std::endl;
  return true;
}

DECLARE_GETTER(PGS_Interface_Getter,"PGS",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *PGS_Interface_Getter::operator()
(const Analysis_Arguments &args) const
{
  return new PGS_Interface
    (args.m_inpath,args.m_infile,args.m_outpath);
}

void PGS_Interface_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"PGS interface";
}
