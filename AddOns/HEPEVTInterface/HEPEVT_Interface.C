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

#ifdef USING__H1_Mode
#ifndef USING__HERACMN
#define USING__HERACMN
#endif
#endif
#ifdef USING__HZTOOL
#ifndef USING__HERACMN
#define USING__HERACMN
#endif
#endif

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

#ifdef USING__HERACMN
extern struct {
  double xsec;
  char gen[8];
  int ihchrg[nmxhep];
  float ntot, wtx;
  char beams[8][2];
} heracmn_;
#define heracmn heracmn_
#endif

#ifdef USING__HZTOOL
void hzfilhep_();
void hzlihep_(int *);

void hzdebug_()
{
  int iprint(1);
  if (msg_LevelIsDebugging()) hzlihep_(&iprint);
}

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevtp_;
#define hepevtp hepevtp_
#endif

}

struct MD_Info {
  int m_n, m_c, m_mo[2], m_da[2];
  inline MD_Info(const int n=0): m_n(n), m_c(0)
  { m_mo[1]=m_mo[0]=m_da[1]=m_da[0]=0; }
};// end of struct MD_Info

typedef std::map<ATOOLS::Particle*,MD_Info> MotherDaughter_Map;

class HEPEVT_Interface: public SHERPA::Analysis_Interface {
private:

  std::string m_inpath, m_infile;

  int m_rotate;

  void RotateEvent(ATOOLS::Blob_List *const bl);

  void ConvertParticle(ATOOLS::Particle *const cp,const int sc,
		       const int mode,MotherDaughter_Map &mdmap);
  void Convert(ATOOLS::Blob_List *const bl,
	       MotherDaughter_Map &mdmap);

  void CheckParticle(ATOOLS::Particle *const cp,const int sc,
		     const int mode,MotherDaughter_Map &mdmap);
  void Check(ATOOLS::Blob_List *const bl,
	     MotherDaughter_Map &mdmap);

public:

  inline HEPEVT_Interface(const std::string &inpath,
			  const std::string &infile,
			  const std::string &outpath):
    Analysis_Interface("HEPEVT"),
    m_inpath(inpath), m_infile(infile) {}

  bool Init();
  bool Run(ATOOLS::Blob_List *const bl);
  bool Finish();

  void ShowSyntax(const int i);

};// end of class HEPEVT_Interface

using namespace SHERPA;
using namespace ATOOLS;

void HEPEVT_Interface::ShowSyntax(const int i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   BEGIN_HEPEVT {\n\n"
	   <<"     ROTATE_EVENTS 0|1\n";
  msg_Out()<<"\n   } END_HEPEVT\n\n"
	   <<"}"<<std::endl;
}

void HEPEVT_Interface::ConvertParticle
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

void HEPEVT_Interface::CheckParticle
(ATOOLS::Particle *const cp,const int sc,
 const int mode,MotherDaughter_Map &mdmap)
{
#ifdef USING__HZTOOL
  if (mdmap.find(cp)==mdmap.end())
    THROW(fatal_error,"Internal error.");
  MD_Info &cmd(mdmap[cp]);
  Blob *pb(cp->ProductionBlob());
  if (pb && pb->NInP()==1 && mode==0)
    CheckParticle(pb->InParticle(0),3,0,mdmap);
  if (cmd.m_c==0) {
  cmd.m_c=1;
  int n=hepevt.nhep;
  if (cmd.m_n!=n) THROW(fatal_error,"Internal error");
  if (hepevtp.jmohep[n][0]!=cmd.m_mo[0] ||
      hepevtp.jmohep[n][1]!=cmd.m_mo[1])
    THROW(fatal_error,"HZFILHEP error JMOHEP");
  if (hepevtp.jdahep[n][0]!=cmd.m_da[0] ||
      hepevtp.jdahep[n][1]!=cmd.m_da[1])
    THROW(fatal_error,"HZFILHEP error JDAHEP");
  if (hepevtp.idhep[n]!=cp->Flav().HepEvt())
    THROW(fatal_error,"HZFILHEP error IDHEP");
  static double cacc(1.0e-6);
  for (short int j=1;j<4;++j)
    if (!IsEqual(hepevtp.phep[n][j-1],cp->Momentum()[j],cacc))
      msg_Error()<<"HZFILHEP error PHEP"<<std::endl;
  if (!IsEqual(hepevtp.phep[n][3],cp->Momentum()[0],cacc))
    msg_Error()<<"HZFILHEP error PHEP4"<<std::endl;
  if (!IsEqual(hepevtp.phep[n][4],
	       Max(0.0,cp->Momentum().Abs2()),cacc))
    msg_Error()<<"HZFILHEP error PHEP5"<<std::endl;
  if (cp->ProductionBlob()) {
    for (short int j=1;j<4;++j)
      if (!IsEqual(hepevtp.vhep[n][j-1],cp->XProd()[j],cacc))
	msg_Error()<<"HZFILHEP error VHEP"<<std::endl;
    if (!IsEqual(hepevtp.vhep[n][3],cp->XProd()[0],cacc))
      msg_Error()<<"HZFILHEP error VHEP4"<<std::endl;
  }
  else {
    for (short int j=0;j<4;++j) 
      if (hepevtp.vhep[n][j]!=0.0)
	THROW(fatal_error,"HZFILHEP error VHEP");
  }
  if (hepevtp.isthep[n]!=sc)
    THROW(fatal_error,"HZFILHEP error ISTHEP");
  ++hepevt.nhep;
  }
  Blob *db(cp->DecayBlob());
  if (db && db->NInP()==1 && mode==0)
    for (int i(0);i<db->NOutP();++i) {
      Particle *np(db->OutParticle(i));
      CheckParticle(np,np->DecayBlob()?3:1,1,mdmap);
    }
#endif
}

void HEPEVT_Interface::Convert(ATOOLS::Blob_List *const bl,
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

void HEPEVT_Interface::Check(ATOOLS::Blob_List *const bl,
			     MotherDaughter_Map &mdmap)
{
#ifdef USING__HZTOOL
  hepevt.nhep=0;
  Blob_List bb(bl->Find(btp::Bunch));
  CheckParticle(bb[0]->InParticle(0),3,-1,mdmap);
  CheckParticle(bb[1]->InParticle(0),3,-1,mdmap);
  Particle_List pl(bl->ExtractParticles(1));
  for (size_t n(0);n<pl.size();++n)
    if (pl[n]->DecayBlob()==NULL)
      CheckParticle(pl[n],1,0,mdmap);
  if (hepevt.nhep!=hepevtp.nhep)
    THROW(fatal_error,"HZFILHEP error NHEP");
#endif
}

void HEPEVT_Interface::RotateEvent(Blob_List *const bl)
{
  if (!m_rotate) return;
  for (size_t i(0);i<bl->size();++i) {
    for (int j(0);j<(*bl)[i]->NInP();++j) {
      Particle *cp((*bl)[i]->InParticle(j));
      if (cp->ProductionBlob()==NULL) {
	Vec4D p(cp->Momentum());
	p=Vec4D(p[0],-p[1],-p[2],-p[3]);
	cp->SetMomentum(p);
      }
    }
    for (int j(0);j<(*bl)[i]->NOutP();++j) {
      Particle *cp((*bl)[i]->OutParticle(j));
      Vec4D p(cp->Momentum());
      p=Vec4D(p[0],-p[1],-p[2],-p[3]);
      cp->SetMomentum(p);
    }
  }
}

bool HEPEVT_Interface::Init()
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
  reader.SetFileBegin("BEGIN_HEPEVT");
  reader.SetFileEnd("END_HEPEVT");
#ifdef USING__H1_Mode
  m_rotate=reader.GetValue<int>("ROTATE_EVENTS",1);
#else
  m_rotate=reader.GetValue<int>("ROTATE_EVENTS",0);
#endif
  msg_Info()<<METHOD<<"(): Set up HEPEVT interface."
	    <<" Rotate = "<<m_rotate<<"."<<std::endl;
#ifdef USING__HERACMN
  MakeFortranString(heracmn.gen,"SHA",8);
  heracmn.xsec=-1.0;
#endif
  return true;
}

bool HEPEVT_Interface::Run(ATOOLS::Blob_List *const bl)
{
  RotateEvent(bl);
  if (!bl->FourMomentumConservation())
    msg_Error()<<METHOD<<"(): Four momentum not conserved."<<std::endl;
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  Blob_Data_Base *xs((*sp)["Weight"]);
  if (xs==NULL) THROW(fatal_error,"No weight information");
  MotherDaughter_Map mdmap;
  Convert(bl,mdmap);
#ifdef USING__HERACMN
  heracmn.wtx=xs->Get<double>();
#endif
#ifdef USING__HZTOOL
  hzfilhep_();
  hzdebug_();
#endif
  Check(bl,mdmap);
  RotateEvent(bl);
  return true;
}

bool HEPEVT_Interface::Finish()
{
  msg_Info()<<METHOD<<"(): Total xs is "
	    <<p_eventhandler->TotalXS()<<" pb +- "
	    <<p_eventhandler->TotalErr()<<" pb."<<std::endl;
#ifdef USING__HERACMN
  heracmn.xsec=p_eventhandler->TotalXS();
#endif
  return true;
}

DECLARE_GETTER(HEPEVT_Interface_Getter,"HEPEVT",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *HEPEVT_Interface_Getter::operator()
(const Analysis_Arguments &args) const
{
  return new HEPEVT_Interface
    (args.m_inpath,args.m_infile,args.m_outpath);
}

void HEPEVT_Interface_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"HEPEVT interface";
}
