#include "SHERPA/Tools/Output_Base.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"

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
  double eventweightlh, alphaqedlh, alphaqcdlh;
  double scalelh[10], spinlh[nmxhep][3];
  int icolorflowlh[nmxhep][2], idruplh;
} hepev4_;
#define hepev4 hepev4_

extern struct {
  int nevsha, istream;
  char shafile[80];
} pgspars_;
#define pgspars pgspars_

  void pgsxxx_(int *mode,int *imode);

  inline void pgsxxx(int mode,int imode=1)
  { pgsxxx_(&mode,&imode); }

}

namespace PGS {

struct MD_Info {
  int m_n, m_c, m_mo[2], m_da[2];
  inline MD_Info(const int n=0): m_n(n), m_c(0)
  { m_mo[1]=m_mo[0]=m_da[1]=m_da[0]=0; }
};// end of struct MD_Info

typedef std::map<ATOOLS::Particle*,MD_Info> MotherDaughter_Map;

class Output_PGS: public SHERPA::Output_Base {
private:

  std::string m_basename, m_ext;

  int m_imode;

  void ConvertParticle(ATOOLS::Particle *const cp,const int sc,
		       const int mode,MotherDaughter_Map &mdmap);
  void Convert(ATOOLS::Blob_List *const bl,
	       MotherDaughter_Map &mdmap);

public:

  Output_PGS(const SHERPA::Output_Arguments &args,const int mode);

  void Header();
  void Footer();

  void Output(ATOOLS::Blob_List *bl,const double weight);

  void ShowSyntax(const int i);

};// end of class Output_PGS

class Output_PGSW: public Output_PGS {};

using namespace SHERPA;
using namespace ATOOLS;
  
Output_PGS::Output_PGS(const SHERPA::Output_Arguments &args,const int mode):
  Output_Base("PGS"), m_imode(mode)
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=".lhe";
}

void Output_PGS::Header()
{
  pgspars.nevsha=rpa->gen.NumberOfEvents();
  MakeFortranString(pgspars.shafile,m_basename+m_ext,80);
  msg_Info()<<METHOD<<"("<<m_imode<<"): {"<<std::endl;
  pgsxxx(1);
  msg_Info()<<"}"<<std::endl;
  hepev4.alphaqedlh=hepev4.alphaqcdlh=0.0;
  for (size_t i(0);i<10;++i) hepev4.scalelh[i]=0.0;
  hepev4.idruplh=0;
}

void Output_PGS::ConvertParticle
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
  hepev4.icolorflowlh[n][0]=hepev4.icolorflowlh[n][1]=0;
  hepev4.spinlh[n][0]=hepev4.spinlh[n][1]=hepev4.spinlh[n][2]=0;
  hepevt.idhep[n]=(long int)cp->Flav();
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

void Output_PGS::Convert(ATOOLS::Blob_List *const bl,
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

void Output_PGS::Output(Blob_List *bl,const double weight)
{
  if (!bl->FourMomentumConservation())
    msg_Error()<<METHOD<<"(): Four momentum not conserved."<<std::endl;
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  Blob_Data_Base *xs((*sp)["Weight"]);
  if (xs==NULL) THROW(fatal_error,"No weight information");
  MotherDaughter_Map mdmap;
  Convert(bl,mdmap);
  hepev4.eventweightlh=xs->Get<double>();
  hepevt.nevhep=rpa->gen.NumberOfGeneratedEvents();
  pgsxxx(2,m_imode);
}

void Output_PGS::Footer()
{
  msg_Info()<<METHOD<<"(): {\n  Total xs is "
	    <<p_eventhandler->TotalXS()<<" pb +- "
	    <<p_eventhandler->TotalErr()<<" pb.\n";
  pgsxxx(3);
  msg_Info()<<"}"<<std::endl;
}

}// end of namespace PGS

using namespace PGS;
using namespace SHERPA;

DECLARE_GETTER(Output_PGS,"PGS",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_PGS>::
operator()(const Output_Arguments &args) const
{
  return new Output_PGS(args,1);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_PGS>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"PGS output unweighted";
}

DECLARE_GETTER(Output_PGSW,"PGS_Weighted",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_PGSW>::
operator()(const Output_Arguments &args) const
{
  return new Output_PGS(args,4);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_PGSW>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"PGS output weighted";
}
