#include "AddOns/Root/RootNtuple_Reader.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>

#ifdef USING__ROOT
#include "TChain.h"
#endif

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace SHERPA {

  class Dummy_Process: public PHASIC::Process_Base {
  public:
    void SetScale(const Scale_Setter_Arguments &args) {}
    void SetKFactor(const KFactor_Setter_Arguments &args) {}
    size_t Size() const { return 0; }
    Process_Base *operator[](const size_t &i) { return NULL; }

    Weight_Info *OneEvent(const int wmode,const int mode=0) { return NULL; }

    double Differential(const ATOOLS::Vec4D_Vector &p) { return 0.0; }
    double Differential2() { return 0.0; }
    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create=false) { return false; }
    void SetLookUp(const bool lookup) {}
  };

  struct RootNTupleReader_Variables {
#ifdef USING__ROOT
    static const Int_t s_kMaxParticle = 100;
    UInt_t m_id;
    Int_t m_nparticle;
    Float_t p_px[s_kMaxParticle];
    Float_t p_py[s_kMaxParticle];
    Float_t p_pz[s_kMaxParticle];
    Float_t p_E[s_kMaxParticle];
    Int_t p_kf[s_kMaxParticle];

    Double_t m_wgt,m_wgt2,m_mewgt,m_mewgt2;
    Double_t m_x1,m_x2,m_x1p,m_x2p,m_mur,m_muf,m_as;
    Int_t m_id1,m_id2,m_nuwgt;
    Double_t p_uwgt[18];
    Short_t m_oqcd;
    Char_t m_type[2];
    TChain* p_f;
#endif
  };
}

RootNtuple_Reader::RootNtuple_Reader(const Input_Arguments &args) :
  Event_Reader_Base(args), 
  m_evtid(0), m_subevtid(0), m_evtcnt(0), m_entries(0), m_evtpos(0),
  p_isr(args.p_isr), m_sargs(NULL,"","")
{
  std::string filename=m_path+m_file+".root";
  msg_Out()<<" Reading from "<<filename<<"\n";
  m_calc=args.p_reader->GetValue<int>("ROOTNTUPLE_CALC",1);
  if (m_calc) msg_Info()<<METHOD<<"(): Ntuple calc mode set to "<<m_calc<<"."<<std::endl;
  m_check=args.p_reader->GetValue<int>("ROOTNTUPLE_CHECK",m_calc&2?1:0);
  if (m_check) msg_Info()<<METHOD<<"(): Ntuple check mode set to "<<m_check<<"."<<std::endl;
  args.p_reader->SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  args.p_reader->RereadInFile();
  std::string scale=args.p_reader->GetValue<std::string>
    ("SCALES","VAR{sqr("+ToString(rpa->gen.Ecms())+")}");
  std::vector<std::string> helpsv;
  if (!args.p_reader->VectorFromFile(helpsv,"COUPLINGS")) helpsv.push_back("Alpha_QCD 1");
  std::string coupling(helpsv.size()?helpsv[0]:"");
  for (size_t i(1);i<helpsv.size();++i) coupling+=" "+helpsv[i];
  m_sargs=Scale_Setter_Arguments(args.p_model,scale,coupling);
  m_sargs.m_nin=2;
#ifdef USING__ROOT
  p_vars = new RootNTupleReader_Variables();
  p_vars->p_f=new TChain("t3");
  size_t bpos(filename.find("[")), epos(filename.find("]",bpos));
  if (bpos==std::string::npos ||
      epos==std::string::npos) {
    p_vars->p_f->Add(filename.c_str());
  }
  else {
    std::string basename(filename.substr(0,bpos));
    std::string suffix(filename.substr(epos+1));
    std::string range(filename.substr(bpos+1,epos-bpos-1));
    size_t spos(range.find("-"));
    if (spos==std::string::npos) THROW(fatal_error,"Ivalid syntax");
    size_t i(ToType<size_t>(range.substr(0,spos)));
    size_t e(ToType<size_t>(range.substr(spos+1)));
#ifdef USING__MPI
    int size=MPI::COMM_WORLD.Get_size();
    int rank=MPI::COMM_WORLD.Get_rank();
    int le=e, nact=size, values[2];
    if (size>1) {
      if (rank==0) {
	msg_Info()<<"MPI Analysis {\n";
	int inc=Max(1,(int)((e-i+1)/nact));
	e=i+inc-1;
	msg_Info()<<"  Rank 0 analyzes "<<basename
		  <<"["<<i<<"-"<<e<<"].\n";
	for (int tag=1;tag<size;++tag) {
	  values[0]=i+tag*inc;
	  values[1]=i+(tag+1)*inc-1;
	  if (tag==nact-1) values[1]=le;
	  MPI::COMM_WORLD.Send(&values,2,MPI::INT,tag,tag);
	  MPI::COMM_WORLD.Recv(&values,2,MPI::INT,MPI::ANY_SOURCE,size+tag);
	  msg_Info()<<"  Rank "<<tag<<" analyzes "
		    <<basename<<"["<<values[0]<<"-"<<values[1]<<"].\n";
	}
	msg_Info()<<"}\n";
      }
      else {
	MPI::COMM_WORLD.Recv(&values,2,MPI::INT,0,rank);
	i=values[0];
	e=values[1];
	MPI::COMM_WORLD.Send(&values,2,MPI::INT,0,size+rank);
      }
    }
#endif
    for (;i<=e;++i) {
      std::string lfile(basename+ToString(i)+suffix);
      if (FileExists(lfile)) p_vars->p_f->Add(lfile.c_str());
    }
  }
  m_entries=p_vars->p_f->GetEntries();
  if(m_entries==0) {
    msg_Error()<<"ERROR: Event file "<<filename<<" does not contain any event."<<std::endl;
    THROW(fatal_error,"Missing input");
  }
  p_vars->p_f->SetBranchAddress("id",&p_vars->m_id);
  p_vars->p_f->SetBranchAddress("nparticle",&p_vars->m_nparticle);
  p_vars->p_f->SetBranchAddress("px",p_vars->p_px);
  p_vars->p_f->SetBranchAddress("py",p_vars->p_py);
  p_vars->p_f->SetBranchAddress("pz",p_vars->p_pz);
  p_vars->p_f->SetBranchAddress("E",p_vars->p_E);

  p_vars->p_f->SetBranchAddress("kf",p_vars->p_kf);
  p_vars->p_f->SetBranchAddress("weight",&p_vars->m_wgt);
  p_vars->p_f->SetBranchAddress("weight2",&p_vars->m_wgt2);
  p_vars->p_f->SetBranchAddress("alphas",&p_vars->m_as);
  p_vars->p_f->SetBranchAddress("me_wgt",&p_vars->m_mewgt);
  p_vars->p_f->SetBranchAddress("me_wgt2",&p_vars->m_mewgt2);
  p_vars->p_f->SetBranchAddress("ren_scale",&p_vars->m_mur);
  p_vars->p_f->SetBranchAddress("fac_scale",&p_vars->m_muf);
  p_vars->p_f->SetBranchAddress("x1",&p_vars->m_x1);
  p_vars->p_f->SetBranchAddress("x2",&p_vars->m_x2);
  p_vars->p_f->SetBranchAddress("x1p",&p_vars->m_x1p);
  p_vars->p_f->SetBranchAddress("x2p",&p_vars->m_x2p);
  p_vars->p_f->SetBranchAddress("id1",&p_vars->m_id1);
  p_vars->p_f->SetBranchAddress("id2",&p_vars->m_id2);
  p_vars->p_f->SetBranchAddress("nuwgt",&p_vars->m_nuwgt);
  p_vars->p_f->SetBranchAddress("usr_wgts",p_vars->p_uwgt);
  p_vars->p_f->SetBranchAddress("alphasPower",&p_vars->m_oqcd);
  p_vars->p_f->SetBranchAddress("part",p_vars->m_type);
#else
  msg_Error()<<"Sherpa must be linked with root to read in root files!"<<endl;
#endif
  Read_Write_Base::AddCommandLine("K_PERP_MEAN_1 0;");
  Read_Write_Base::AddCommandLine("K_PERP_MEAN_2 0;");
  Read_Write_Base::AddCommandLine("K_PERP_SIGMA_1 0;");
  Read_Write_Base::AddCommandLine("K_PERP_SIGMA_2 0;");
}

RootNtuple_Reader::~RootNtuple_Reader()
{
  for (std::map<int,PHASIC::Process_Base*>::iterator
	 sit(m_procs.begin());sit!=m_procs.end();++sit) delete sit->second;
  for (std::map<int,PHASIC::Scale_Setter_Base*>::iterator
	 sit(m_scales.begin());sit!=m_scales.end();++sit) delete sit->second;
}


void RootNtuple_Reader::CloseFile() {
#ifdef USING__ROOT
  delete p_vars->p_f;
  delete p_vars;
#endif
}





bool RootNtuple_Reader::FillBlobs(Blob_List * blobs) 
{
  bool result=ReadInFullEvent(blobs);
  if (result==0) rpa->gen.SetNumberOfEvents(rpa->gen.NumberOfGeneratedEvents());
  
  long nev=rpa->gen.NumberOfEvents();
  if(nev==rpa->gen.NumberOfGeneratedEvents()) CloseFile();
  return result;
}

bool RootNtuple_Reader::ReadInEntry() 
{
  if (m_evtpos>=m_entries) return 0;
#ifdef USING__ROOT
  p_vars->p_f->GetEntry(m_evtpos);
  m_evtpos++;
  m_evtid=p_vars->m_id;
#endif
  return 1;
}

double RootNtuple_Reader::CalculateWeight
(const double &mur2,const double &muf2,const int mode) const
{
#ifdef USING__ROOT
  Flavour fl1((kf_code)abs(p_vars->m_id1),p_vars->m_id1<0);
  Flavour fl2((kf_code)abs(p_vars->m_id2),p_vars->m_id2<0);
  p_isr->PDF(0)->Calculate(p_vars->m_x1,muf2);
  p_isr->PDF(1)->Calculate(p_vars->m_x2,muf2);
  double fa=p_isr->PDF(0)->GetXPDF(fl1)/p_vars->m_x1;
  double fb=p_isr->PDF(1)->GetXPDF(fl2)/p_vars->m_x2;
  double asf=pow((*MODEL::as)(mur2)/p_vars->m_as,p_vars->m_oqcd);
  if (mode==0) {
    return p_vars->m_mewgt*asf*fa*fb;
  }
  else if (mode==2) {
    return p_vars->m_mewgt2*asf*fa*fb;
  }
  else if (mode==1) {
    double w[9];
    double lr=log(mur2/sqr(p_vars->m_mur)), lf=log(muf2/sqr(p_vars->m_muf));
    w[0]=p_vars->m_mewgt+p_vars->p_uwgt[0]*lr+p_vars->p_uwgt[1]*lr*lr/2.0;
    for (int i(1);i<9;++i) w[i]=p_vars->p_uwgt[i+1]+p_vars->p_uwgt[i+9]*lf;
    double wgt=w[0]*fa*fb;
    double faq=0.0, faqx=0.0, fag=0.0, fagx=0.0;
    double fbq=0.0, fbqx=0.0, fbg=0.0, fbgx=0.0;
    Flavour quark(kf_quark), gluon(kf_gluon);
    if (fl1.IsQuark()) {
      faq=fa;
      fag=p_isr->PDF(0)->GetXPDF(gluon)/p_vars->m_x1;
      p_isr->PDF(0)->Calculate(p_vars->m_x1/p_vars->m_x1p,muf2);
      faqx=p_isr->PDF(0)->GetXPDF(fl1)/p_vars->m_x1;
      fagx=p_isr->PDF(0)->GetXPDF(gluon)/p_vars->m_x1;
    }
    else {
      fag=fa;
      for (size_t i=0;i<quark.Size();++i)
	faq+=p_isr->PDF(0)->GetXPDF(quark[i])/p_vars->m_x1;
      p_isr->PDF(0)->Calculate(p_vars->m_x1/p_vars->m_x1p,muf2);
      fagx=p_isr->PDF(0)->GetXPDF(fl1)/p_vars->m_x1;
      for (size_t i=0;i<quark.Size();++i)
	faqx+=p_isr->PDF(0)->GetXPDF(quark[i])/p_vars->m_x1;
    }
    if (fl2.IsQuark()) {
      fbq=fb;
      fbg=p_isr->PDF(1)->GetXPDF(gluon)/p_vars->m_x2;
      p_isr->PDF(1)->Calculate(p_vars->m_x2/p_vars->m_x2p,muf2);
      fbqx=p_isr->PDF(1)->GetXPDF(fl2)/p_vars->m_x2;
      fbgx=p_isr->PDF(1)->GetXPDF(gluon)/p_vars->m_x2;
    }
    else {
      fbg=fb;
      for (size_t i=0;i<quark.Size();++i)
	fbq+=p_isr->PDF(1)->GetXPDF(quark[i])/p_vars->m_x2;
      p_isr->PDF(1)->Calculate(p_vars->m_x2/p_vars->m_x2p,muf2);
      fbgx=p_isr->PDF(1)->GetXPDF(fl2)/p_vars->m_x2;
      for (size_t i=0;i<quark.Size();++i)
	fbqx+=p_isr->PDF(1)->GetXPDF(quark[i])/p_vars->m_x2;
    }
    wgt+=(faq*w[1]+faqx*w[2]+fag*w[3]+fagx*w[4])*fb;
    wgt+=(fbq*w[5]+fbqx*w[6]+fbg*w[7]+fbgx*w[8])*fa;
    return wgt*asf;
  }
  else {
    THROW(not_implemented,"Invalid option");
  }
#else
  return -1.0;
#endif
}

bool RootNtuple_Reader::ReadInFullEvent(Blob_List * blobs) 
{
  if (!m_nlos.empty()) {
    for (size_t i=0;i<m_nlos.size();i++) {
      delete[] m_nlos[i]->p_fl;
      delete[] m_nlos[i]->p_mom;
      delete m_nlos[i];
    }
    m_nlos.clear();
  }
  if (m_evtid==0) if (!ReadInEntry()) return 0;
  Blob         *signalblob=blobs->FindFirst(btp::Signal_Process);
  m_weight = 0.;
  signalblob->SetTypeSpec("NLO");
  signalblob->SetId();
  signalblob->SetPosition(Vec4D(0.,0.,0.,0.));
  signalblob->SetStatus(blob_status::code(30));
  
#ifdef USING__ROOT
  size_t currentid=m_evtid;
  int id1, id2;
  double x1, x2;
  while (currentid==m_evtid) {
    Vec4D *moms = new Vec4D[2+p_vars->m_nparticle];
    Flavour *flav = new Flavour[2+p_vars->m_nparticle];
    for (int i=0;i<p_vars->m_nparticle;i++) {
      moms[i+2]=Vec4D(p_vars->p_E[i],p_vars->p_px[i],
		      p_vars->p_py[i],p_vars->p_pz[i]);
      flav[i+2]=Flavour(abs(p_vars->p_kf[i]),p_vars->p_kf[i]<0);
    }
    m_nlos.push_back(new NLO_subevt(p_vars->m_nparticle+2,NULL,flav,moms));
    m_nlos.back()->m_result=p_vars->m_wgt2;
    id1=p_vars->m_id1;
    id2=p_vars->m_id2;
    x1=p_vars->m_x1;
    x2=p_vars->m_x2;
    if (m_calc) {
      Vec4D_Vector p(2+p_vars->m_nparticle);
      for (int i=0;i<p_vars->m_nparticle;++i) {
	p[2+i]=Vec4D(p_vars->p_E[i],p_vars->p_px[i],
		     p_vars->p_py[i],p_vars->p_pz[i]);
      }
      if (m_scales.find(p_vars->m_nparticle)==m_scales.end()) {
	Process_Info pi;
	Flavour fl1((kf_code)abs(p_vars->m_id1),p_vars->m_id1<0);
	Flavour fl2((kf_code)abs(p_vars->m_id2),p_vars->m_id2<0);
	pi.m_ii.m_ps.push_back(Subprocess_Info(fl1));
	pi.m_ii.m_ps.push_back(Subprocess_Info(fl2));
	pi.m_oqcd=p_vars->m_oqcd;
	for (int i=0;i<p_vars->m_nparticle;++i)
	  pi.m_fi.m_ps.push_back(Subprocess_Info(flav[i+2]));
	m_sargs.p_proc=m_procs[p_vars->m_nparticle] = new Dummy_Process();
	m_sargs.p_proc->Init(pi,NULL,NULL);
	m_sargs.m_nout=p_vars->m_nparticle;
	m_scales[p_vars->m_nparticle] =
	  Scale_Setter_Base::Scale_Getter_Function::
	  GetObject(m_sargs.m_scale,m_sargs);
	if (m_scales[p_vars->m_nparticle]==NULL)
	  THROW(fatal_error,"Invalid scale scheme");
      }
      Scale_Setter_Base *scale(m_scales[p_vars->m_nparticle]);
      scale->CalculateScale(p);
      double weight=CalculateWeight
	(scale->Scale(stp::ren),scale->Scale(stp::fac),
	 p_vars->m_nuwgt?1:2);
      m_nlos.back()->m_result=weight;
      m_weight+=weight;
      if (m_check) {
	msg_Debugging()<<METHOD<<"(): "<<p_vars->m_type<<" computed "
		       <<weight<<", stored "<<p_vars->m_wgt2
		       <<", rel. diff. "<<weight/p_vars->m_wgt2-1.0<<".\n";
	if (!IsEqual(weight,p_vars->m_wgt2,rpa->gen.Accu()))
	  msg_Error()<<METHOD<<"(): "<<p_vars->m_type<<" weights differ by "
		     <<(weight/p_vars->m_wgt2-1.0)<<".\n  computed "
		     <<weight<<", stored "<<p_vars->m_wgt2<<"."<<std::endl;
      }
    }
    else if (m_check) {
      double weight=CalculateWeight
	(sqr(p_vars->m_mur),sqr(p_vars->m_muf),p_vars->m_nuwgt?1:2);
      m_weight+=weight;
      msg_Debugging()<<METHOD<<"(): "<<p_vars->m_type<<" computed "
		     <<weight<<", stored "<<p_vars->m_wgt2
		     <<", rel. diff. "<<weight/p_vars->m_wgt2-1.0<<".\n";
      if (!IsEqual(weight,p_vars->m_wgt2,rpa->gen.Accu()))
	msg_Error()<<METHOD<<"(): "<<p_vars->m_type<<" weights differ by "
		   <<(weight/p_vars->m_wgt2-1.0)<<".\n  computed "
		   <<weight<<", stored "<<p_vars->m_wgt2<<"."<<std::endl;
    }
    else {
      m_weight+=p_vars->m_wgt2;
    }
    if (!ReadInEntry()) m_evtid=0;
  }  
  Flavour fl1((kf_code)abs(id1),id1<0);
  Flavour fl2((kf_code)abs(id2),id2<0);
  Particle *part1=new Particle(0,fl1,x1*rpa->gen.PBeam(0));
  Particle *part2=new Particle(1,fl2,x2*rpa->gen.PBeam(1));
  signalblob->AddToInParticles(part1);
  signalblob->AddToInParticles(part2);
  for (size_t i=2;i<m_nlos.back()->m_n;++i) {
    Particle *part=new Particle
      (i,m_nlos.back()->p_fl[i],m_nlos.back()->p_mom[i]);
    signalblob->AddToOutParticles(part);
  }
#endif
  signalblob->SetStatus(blob_status::needs_beams);
  signalblob->SetWeight(m_weight);
  signalblob->AddData("Weight",new Blob_Data<double>(m_weight));
  signalblob->AddData("MEWeight",new Blob_Data<double>
		      (p_vars->m_nuwgt?p_vars->m_mewgt:p_vars->m_mewgt2));
  signalblob->AddData("Trials",new Blob_Data<double>(1.));
  signalblob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(&m_nlos));
  signalblob->AddData("Weight_Norm",new Blob_Data<double>(1.0));
  m_evtcnt++;  
  return 1;
}

DECLARE_GETTER(RootNtuple_Reader,"Root",
	       Event_Reader_Base,Input_Arguments);

Event_Reader_Base *ATOOLS::Getter
<Event_Reader_Base,Input_Arguments,RootNtuple_Reader>::
operator()(const Input_Arguments &args) const
{
  return new RootNtuple_Reader(args);
}

void ATOOLS::Getter
<Event_Reader_Base,Input_Arguments,RootNtuple_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Root NTuple input";
}

