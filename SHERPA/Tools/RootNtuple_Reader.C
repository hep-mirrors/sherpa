#include "ATOOLS/Org/CXXFLAGS.H"
#include "SHERPA/Tools/RootNtuple_Reader.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>

#ifdef USING__ROOT
#include "TChain.h"
#endif

#ifdef USING__MPI
#include "mpi.h"
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
    TChain* p_f;
#endif
  };
}

RootNtuple_Reader::RootNtuple_Reader(const std::string & path,const std::string & file,int mode,
				     MODEL::Model_Base *const model,PDF::ISR_Handler *const isr) :
  Event_Reader_Base(path,file), 
  m_eventmode(mode), m_evtid(0), m_subevtid(0), m_evtcnt(0), m_entries(0), m_evtpos(0),
  p_isr(isr), m_sargs(NULL,"","")
{
  std::string filename=m_path+m_file;
  msg_Out()<<" Reading from "<<filename<<"\n";
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  m_calc=dataread.GetValue<int>("ROOTNTUPLE_CALC",0);
  if (m_calc) msg_Info()<<METHOD<<"(): Ntuple calc mode set to "<<m_calc<<"."<<std::endl;
  m_check=dataread.GetValue<int>("ROOTNTUPLE_CHECK",m_calc&2?1:0);
  if (m_check) msg_Info()<<METHOD<<"(): Ntuple check mode set to "<<m_check<<"."<<std::endl;
  m_oqcd=dataread.GetValue<int>("ROOTNTUPLE_OQCD",0);
  if (m_oqcd) msg_Info()<<METHOD<<"(): Ntuple O(QCD) set to "<<m_oqcd<<"."<<std::endl;
  dataread.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  dataread.RereadInFile();
  std::string scale=dataread.GetValue<std::string>
    ("SCALES","VAR{sqr("+ToString(rpa->gen.Ecms())+")}");
  std::vector<std::string> helpsv;
  if (!dataread.VectorFromFile(helpsv,"COUPLINGS")) helpsv.push_back("Alpha_QCD 1");
  std::string coupling(helpsv.size()?helpsv[0]:"");
  for (size_t i(1);i<helpsv.size();++i) coupling+=" "+helpsv[i];
  m_sargs=Scale_Setter_Arguments(model,scale,coupling);
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
    exh->MPISync();
    int size=MPI::COMM_WORLD.Get_size();
    int rank=MPI::COMM_WORLD.Get_rank();
    int le=e, nact=1, values[2];
    if (size>1) {
      if (rank==0) {
	msg_Info()<<"MPI Analysis {\n";
	for (int tag=1;tag<size;++tag) if (exh->MPIStat(tag)) ++nact;
	int inc=Max(1,(int)((e-i+1)/nact));
	e=i+inc-1;
	msg_Info()<<"  Rank 0 analyzes "<<basename
		  <<"["<<i<<"-"<<e<<"].\n";
	for (int tag=1;tag<size;++tag) {
	  if (!exh->MPIStat(tag)) continue;
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

  msg_Out()<<"Event mode = "<<m_eventmode<<std::endl;
#else
  msg_Error()<<"Sherpa must be linked with root to read in root files!"<<endl;
#endif
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
  bool result;
  switch (m_eventmode) {
  case 1: result=ReadInFullEvent(blobs);
    break;
  default: result=ReadInSubEvent(blobs);
  }
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
  double asf=pow((*MODEL::as)(mur2)/p_vars->m_as,m_oqcd);
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

bool RootNtuple_Reader::ReadInSubEvent(Blob_List * blobs) 
{
  if (!ReadInEntry()) return 0;
  Blob         *signalblob=blobs->FindFirst(btp::Signal_Process);
  signalblob->SetTypeSpec("NLOsubevt");
  signalblob->SetId();
  signalblob->SetPosition(Vec4D(0.,0.,0.,0.));
  signalblob->SetStatus(blob_status::code(30));
#ifdef USING__ROOT
  if (!m_calc) m_weight=p_vars->m_wgt;
  else {
    Vec4D_Vector p(2+p_vars->m_nparticle);
    for (int i=0;i<p_vars->m_nparticle;++i)
      p[2+i]=Vec4D(p_vars->p_E[i],p_vars->p_px[i],
		   p_vars->p_py[i],p_vars->p_pz[i]);
    if (m_scales.find(p_vars->m_nparticle)==m_scales.end()) {
      Process_Info pi;
      Flavour fl1((kf_code)abs(p_vars->m_id1),p_vars->m_id1<0);
      Flavour fl2((kf_code)abs(p_vars->m_id2),p_vars->m_id2<0);
      pi.m_ii.m_ps.push_back(Subprocess_Info(fl1));
      pi.m_ii.m_ps.push_back(Subprocess_Info(fl2));
      for (int i=0;i<p_vars->m_nparticle;++i)
	pi.m_fi.m_ps.push_back
	  (Subprocess_Info(Flavour(abs(p_vars->p_kf[i]),
				   p_vars->p_kf[i]<0)));
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
    scale->CalculateScale(p,0);
    double weight=CalculateWeight
      (scale->Scale(stp::ren),scale->Scale(stp::fac),
       p_vars->m_nuwgt?1:0);
    if (m_check) {
      msg_Debugging()<<METHOD<<"(): computed "<<weight
		     <<", stored "<<p_vars->m_wgt
		     <<", rel. diff. "<<weight/p_vars->m_wgt-1.0<<".\n";
      if (!IsEqual(weight,p_vars->m_wgt,rpa->gen.Accu()))
	msg_Error()<<METHOD<<"(): Weights differ by "
		   <<(weight/p_vars->m_wgt-1.0)<<".\n  computed "
		   <<weight<<", stored "<<p_vars->m_wgt<<"."<<std::endl;
    }
    m_weight=weight;
  }
  if (m_check && !m_calc) {
    double weight=CalculateWeight
      (sqr(p_vars->m_mur),sqr(p_vars->m_muf),p_vars->m_nuwgt?1:0);
    msg_Debugging()<<METHOD<<"(): computed "<<weight
		   <<", stored "<<m_weight
		   <<", rel. diff. "<<weight/m_weight-1.0<<".\n";
    if (!IsEqual(weight,m_weight,rpa->gen.Accu()))
      msg_Error()<<METHOD<<"(): Weights differ by "
		 <<(weight/m_weight-1.0)<<".\n  computed "
		 <<weight<<", stored "<<m_weight<<"."<<std::endl;
  }
  signalblob->SetWeight(m_weight);
  for (int i=0;i<p_vars->m_nparticle;i++) {
    Particle *part=new Particle(i,Flavour(abs(p_vars->p_kf[i]),
					  p_vars->p_kf[i]<0),
				Vec4D(p_vars->p_E[i],
				      p_vars->p_px[i],
				      p_vars->p_py[i],
				      p_vars->p_pz[i]));
    signalblob->AddToOutParticles(part);
  }
  signalblob->AddData("Weight",new Blob_Data<double>(m_weight));
  signalblob->AddData("Trials",new Blob_Data<double>(1.));
#endif


  m_evtcnt++;
  return 1;
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
    m_nlos.back()->m_flip=0;
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
      scale->CalculateScale(p,0);
      double weight=CalculateWeight
	(scale->Scale(stp::ren),scale->Scale(stp::fac),
	 p_vars->m_nuwgt?1:2);
      m_nlos.back()->m_result=weight;
      m_weight+=weight;
      if (m_check) {
	msg_Debugging()<<METHOD<<"(): computed "<<weight
		       <<", stored "<<p_vars->m_wgt2
		       <<", rel. diff. "<<weight/p_vars->m_wgt2-1.0<<".\n";
	if (!IsEqual(weight,p_vars->m_wgt2,rpa->gen.Accu()))
	  msg_Error()<<METHOD<<"(): Weights differ by "
		     <<(weight/p_vars->m_wgt2-1.0)<<".\n  computed "
		     <<weight<<", stored "<<p_vars->m_wgt2<<"."<<std::endl;
      }
    }
    else if (m_check) {
      double weight=CalculateWeight
	(sqr(p_vars->m_mur),sqr(p_vars->m_muf),p_vars->m_nuwgt?1:2);
      m_weight+=weight;
      msg_Debugging()<<METHOD<<"(): computed "<<weight
		     <<", stored "<<p_vars->m_wgt2
		     <<", rel. diff. "<<weight/p_vars->m_wgt2-1.0<<".\n";
      if (!IsEqual(weight,p_vars->m_wgt2,rpa->gen.Accu()))
	msg_Error()<<METHOD<<"(): Weights differ by "
		   <<(weight/p_vars->m_wgt2-1.0)<<".\n  computed "
		   <<weight<<", stored "<<p_vars->m_wgt2<<"."<<std::endl;
    }
    else {
      m_weight+=p_vars->m_wgt2;
    }
    if (!ReadInEntry()) m_evtid=0;
  }  
  for (size_t i=2;i<m_nlos.back()->m_n;++i) {
    Particle *part=new Particle
      (i,m_nlos.back()->p_fl[i],m_nlos.back()->p_mom[i]);
    signalblob->AddToOutParticles(part);
  }
#endif
  signalblob->SetWeight(m_weight);
  signalblob->AddData("Weight",new Blob_Data<double>(m_weight));
  signalblob->AddData("Trials",new Blob_Data<double>(1.));
  signalblob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(&m_nlos));
  m_evtcnt++;  
  return 1;
}
