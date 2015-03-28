#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "AddOns/HZTool/HZTool_Wrapper.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"

namespace HZTOOL {

typedef void (*HZTool_Analysis)(int *flag);

class HZTool_Interface: public SHERPA::Analysis_Interface {
public:

  static HZTool_Interface *s_hztool;

  class Destructor {
  public:
    inline ~Destructor() { if (s_hztool) delete s_hztool; }
  };// end of class Destructor
  static Destructor s_destructor;

private:

  std::string m_inpath, m_infile, m_outpath, m_hzfile;

  size_t m_nevt, m_xsnevt, m_flags[4];
  double m_xssum, m_nchsum, m_nsum;
  bool   m_finished, m_check;

  std::vector<HZTool_Analysis> m_analyses;
  std::vector<std::string>     m_tags;

  void ConvertParticle(ATOOLS::Particle *const cp,const int sc);
  void Convert(ATOOLS::Blob_List *const bl);

  void CheckParticle(ATOOLS::Particle *const cp,const int sc);
  void Check(ATOOLS::Blob_List *const bl);

public:

  inline HZTool_Interface(const std::string &inpath,
			  const std::string &infile,
			  const std::string &outpath):
    Analysis_Interface("HZTool"),
    m_inpath(inpath), m_infile(infile), m_outpath(outpath),
    m_nevt(0), m_xsnevt(10000),
    m_xssum(0.0), m_nchsum(0.0), m_nsum(0.0),
    m_finished(false) {}
  ~HZTool_Interface() 
  { if (!m_finished) Finish(); s_hztool=NULL; }

  bool Init();
  bool Run(ATOOLS::Blob_List *const bl);
  bool Finish();

  void ShowSyntax(const int i);

  void HZxxxx(int *flag);

};// end of class HZTool_Interface

extern "C" {
  void structm_() {}
  void kpisflag_() {}
  void pdfset_() {}
  void vp2_() {}
  void pylist_() {}
  void lulist_() {}
  void luexec_() {}
}

using namespace SHERPA;
using namespace ATOOLS;

HZTool_Interface *HZTool_Interface::s_hztool(NULL);
HZTool_Interface::Destructor HZTool_Interface::s_destructor;

extern "C" void hzxxxx_(int *flag) 
{ 
  HZTool_Interface::s_hztool->HZxxxx(flag);
}

void HZTool_Interface::HZxxxx(int *flag)
{
  for (size_t i(0);i<m_analyses.size();++i) {
    if (*flag==3) msg_Info()<<"Writing out '"
			   <<m_tags[i]<<"'."<<std::endl;
    int cflag(m_flags[*flag]);
    m_analyses[i](&cflag);
  }
}

void HZTool_Interface::ShowSyntax(const int i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   BEGIN_HZTOOL {\n\n"
	   <<"     HISTO_NAME histogram name\n"
	   <<"     XS_EVENTS  events for estimating\n"
	   <<"                \\sigma_{tot} and <N_{chg}>\n"
	   <<"     EVT_CHECK  checkflag\n"
	   <<"     HZ_FLAGS   initflag runflag finishflag\n"
	   <<"     HZ_ENABLE  hzxxxx\n";
  msg_Out()<<"\n   } END_HZTOOL\n\n"
	   <<"}"<<std::endl;
}

void HZTool_Interface::ConvertParticle
(ATOOLS::Particle *const cp,const int sc)
{
  int n=hepevt.nhep;
  hepevt.jmohep[n][0]=hepevt.jmohep[n][1]=0;
  hepevt.jdahep[n][0]=hepevt.jdahep[n][1]=0;
  hepevt.idhep[n]=(long int) cp->Flav();
  for (short int j=1;j<4;++j) 
    hepevt.phep[n][j-1]=cp->Momentum()[j];
  hepevt.phep[n][3]=cp->Momentum()[0];
  hepevt.phep[n][4]=Max(0.0,cp->Momentum().Abs2());
  if (cp->ProductionBlob()) {
    for (short int j=1;j<4;++j)
      hepevt.vhep[n][j-1]=cp->XProd()[j];
    hepevt.vhep[n][3]=cp->XProd()[0];
  }
  else {
    for (short int j=0;j<4;++j) hepevt.vhep[n][j]=0.0;
  }
  hepevt.isthep[n]=sc;
  ++hepevt.nhep;
}

void HZTool_Interface::CheckParticle
(ATOOLS::Particle *const cp,const int sc)
{
  int n=hepevt.nhep;
  if (hepevtp.jmohep[n][0]!=0 ||
      hepevtp.jmohep[n][1]!=0 ||
      hepevtp.jdahep[n][0]!=0 ||
      hepevtp.jdahep[n][1]!=0)
    THROW(fatal_error,"HZFILHEP error JMO/DAHEP");
  if (hepevtp.idhep[n]!=(long int) cp->Flav())
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

void HZTool_Interface::Convert(ATOOLS::Blob_List *const bl)
{
  hepevt.nhep=0;
  Blob_List bb(bl->Find(btp::Bunch));
  ConvertParticle(bb[0]->InParticle(0),3);
  ConvertParticle(bb[1]->InParticle(0),3);
  if (bb[0]->InParticle(0)->Flav().IsLepton() &&
      bb[1]->InParticle(0)->Flav().IsHadron()) {
    // add virtual photon for crappy hztool routines
    // assume only one hard lepton is generated in signal
    Vec4D pp(bb[0]->InParticle(0)->Momentum());
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    for (int i(0);i<sp->NOutP();++i)
      if (sp->OutParticle(i)->Flav()==
	  bb[0]->InParticle(0)->Flav()) {
	pp-=sp->OutParticle(i)->Momentum();
	break;
      }
    Particle vp(1,Flavour(kf_photon),pp);
    ConvertParticle(&vp,3);
  }
  Particle_List pl(bl->ExtractParticles(1));
  for (size_t n(0);n<pl.size();++n)
    if (pl[n]->DecayBlob()==NULL) ConvertParticle(pl[n],1);
}

void HZTool_Interface::Check(ATOOLS::Blob_List *const bl)
{
  if (!m_check) return;
  hepevt.nhep=0;
  Blob_List bb(bl->Find(btp::Bunch));
  CheckParticle(bb[0]->InParticle(0),3);
  CheckParticle(bb[1]->InParticle(0),3);
  if (bb[0]->InParticle(0)->Flav().IsLepton() &&
      bb[1]->InParticle(0)->Flav().IsHadron()) {
    // add virtual photon for crappy hztool routines
    // assume only one hard lepton is generated in signal
    Vec4D pp(bb[0]->InParticle(0)->Momentum());
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    for (int i(0);i<sp->NOutP();++i)
      if (sp->OutParticle(i)->Flav()==
	  bb[0]->InParticle(0)->Flav()) {
	pp-=sp->OutParticle(i)->Momentum();
	break;
      }
    Particle vp(1,Flavour(kf_photon),pp);
    CheckParticle(&vp,3);
  }
  Particle_List pl(bl->ExtractParticles(1));
  for (size_t n(0);n<pl.size();++n)
    if (pl[n]->DecayBlob()==NULL) CheckParticle(pl[n],1);
  if (hepevt.nhep!=hepevtp.nhep)
    THROW(fatal_error,"HZFILHEP error");
}

bool HZTool_Interface::Init()
{
  if (m_nevt==0) {
    msg_Info()<<METHOD<<"(): {"<<std::endl;
    {
      msg_Indent();
      while (m_outpath.length() &&
	     m_outpath[m_outpath.length()-1]=='/')
	m_outpath.erase(m_outpath.length()-1,1);
      std::string outpath(m_outpath);
      for (size_t i(0);i<m_outpath.length();++i)
	m_outpath[i]=tolower(m_outpath[i]);
      if (m_outpath!=outpath)
	msg_Info()<<"Changed path name '"<<outpath<<"' -> '"<<m_outpath
		  <<"' to comply with hbook rules."<<std::endl;
      MakeDir(m_outpath);
      Data_Reader reader(" ",";","//");
      reader.AddWordSeparator("\t");
      reader.SetAddCommandLine(false);
      reader.SetInputPath(m_inpath);
      std::string infile(m_infile);
      if (infile.find('|')!=std::string::npos)
	infile=infile.substr(0,infile.find('|'));
      reader.SetInputFile(infile);
      reader.AddComment("#");
      reader.SetFileBegin("BEGIN_HZTOOL");
      reader.SetFileEnd("END_HZTOOL");
      reader.AddFileBegin("BEGIN_HZTOOL{");
      reader.AddFileEnd("}END_HZTOOL");
      m_hzfile=reader.GetValue<std::string>("HISTO_NAME","histo");
      MakeFortranString(hzhname.hname,m_outpath+"/"+m_hzfile,128);
      m_xsnevt=reader.GetValue<int>("XS_EVENTS",10000);
      m_check=reader.GetValue<int>("EVT_CHECK",1);
      msg_Info()<<"Using "<<m_xsnevt
		<<" events to estimate cross section."<<std::endl;
      std::vector<int> helpiv;
      if (reader.VectorFromFile(helpiv,"HZ_FLAGS") && 
	  helpiv.size()==3)
	for (size_t i(1);i<=3;++i) m_flags[i]=helpiv[i-1];
      else for (size_t i(1);i<=3;++i) m_flags[i]=i;
      std::vector<std::vector<std::string> > helpsvv;
      reader.MatrixFromFile(helpsvv,"HZ_ENABLE");
      for (size_t i(0);i<helpsvv.size();++i) {
	for (size_t j(0);j<helpsvv[i].size();++j) {
	  msg_Info()<<"Set up '"<<helpsvv[i][j]<<"' ... "<<std::flush;
	  void *func(s_loader->GetLibraryFunction
		     ("SherpaHZToolAnalysis",helpsvv[i][j]+"_"));
	  if (func==NULL) {
	    msg_Info()<<"not found."<<std::endl;
	  }
	  else {
	    m_analyses.push_back((HZTool_Analysis)func);
	    m_tags.push_back(helpsvv[i][j]);
	    msg_Info()<<"found."<<std::endl;
	  }
	}
      }
    }
    msg_Info()<<"}"<<std::endl;
  }
  if (m_nevt<=m_xsnevt) return false;
  msg_Info()<<"\n"<<METHOD<<"(): {"<<std::endl;
  hzinit(m_xssum/m_nsum,m_nchsum/m_xssum);
  if (!ChMod(m_outpath+"/"+m_hzfile,0666))
    msg_Info()<<"  Cannot chmod '"<<m_outpath<<"/"<<m_hzfile<<"'.\n";
  msg_Info()<<"}"<<std::endl;
  return true;
}

bool HZTool_Interface::Run(ATOOLS::Blob_List *const bl)
{
  s_hztool=this;
  if (m_nevt<=m_xsnevt) {
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    double cxs((*sp)["Weight"]->Get<double>());
    m_nsum+=(*sp)["Trials"]->Get<double>();
    m_xssum+=cxs;
    int nch=0;
    Particle_List pl(bl->ExtractParticles(1));
    for (size_t i(0);i<pl.size();++i)
      if (pl[i]->DecayBlob()==NULL && 
	  pl[i]->Flav().IntCharge()!=0) ++nch;
    m_nchsum+=cxs*nch;
    ++m_nevt;
    Convert(bl);
    heracmn.wtx=cxs;
    return true;
  }
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
  if (!bl->FourMomentumConservation())
    msg_Error()<<METHOD<<"(): Four momentum not conserved."<<std::endl;
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  Blob_Data_Base *xs((*sp)["Weight"]);
  if (xs==NULL) THROW(fatal_error,"No weight information");
  double wgt(xs->Get<double>());
  Convert(bl);
  hzevnt(wgt);
  Check(bl);
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
  return true;
}

bool HZTool_Interface::Finish()
{
  if (m_nevt<=m_xsnevt) return true;
  msg_Info()<<METHOD<<"(): {"<<std::endl;
  {
    msg_Indent();
    hzfinl();
  }
  msg_Info()<<"}"<<std::endl;
  return m_finished=true;
}

}// end of namespace HZTOOL

using namespace HZTOOL;
using namespace SHERPA;

DECLARE_GETTER(HZTool_Interface,"HZTool",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,HZTool_Interface>::
operator()(const Analysis_Arguments &args) const
{
  return new HZTool_Interface
    (args.m_inpath,args.m_infile,args.m_outpath);
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,
		    HZTool_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HZTool interface";
}
