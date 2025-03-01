#include "AddOns/Root/RootNtuple_Reader.H"
#include "AddOns/Root/Output_RootNtuple.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/MINLO_Scale_Setter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Phys/Variations.H"
#include <iostream>

#ifdef USING__ROOT
#include "TChain.h"
#include "TFile.h"
#include "TLeaf.h"
#endif

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace SHERPA {

  bool RR_Process_Info::operator<(const RR_Process_Info pi) const
  {
    if (m_type[0]<pi.m_type[0]) return true;
    if (m_type[0]>pi.m_type[0]) return false;
    if (m_fl.size()<pi.m_fl.size()) return true;
    if (m_fl.size()>pi.m_fl.size()) return false;
    return m_fl<pi.m_fl;
  }

  std::ostream &operator<<(std::ostream &s,const RR_Process_Info &pi)
  {
    return s<<"{"<<pi.m_type<<","<<pi.m_fl.size()<<","<<pi.m_fl<<"}";
  }

  class Dummy_Process: public PHASIC::Process_Base {
  public:
    void SetScale(const Scale_Setter_Arguments& args) override
    {
      Scale_Setter_Arguments cargs(args);
      cargs.p_proc=this;
      p_scale = Scale_Setter_Base::Scale_Getter_Function::
	GetObject(cargs.m_scale,cargs);
      if (p_scale==NULL) THROW(fatal_error,"Invalid scale scheme");
    }
    void SetKFactor(const KFactor_Setter_Arguments& args) override
    {
      KFactor_Setter_Arguments cargs(args);
      cargs.p_proc=this;
      m_pinfo.m_kfactor=cargs.m_kfac;
      p_kfactor = KFactor_Setter_Base::KFactor_Getter_Function::
	GetObject(cargs.m_kfac,cargs);
      if (p_kfactor==NULL) THROW(fatal_error,"Invalid kfactor scheme");
    }
    size_t Size() const override { return 0; }
    Process_Base* operator[](const size_t& i) override { return NULL; }

    Weight_Info* OneEvent(const int wmode,
                          Variations_Mode varmode=Variations_Mode::all,
                          const int mode=0)
    {
      return NULL;
    }

    ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector&,
                                     ATOOLS::Variations_Mode) override
    {
      return 0.0;
    }
    double Differential2() { return 0.0; }
    bool CalculateTotalXSec(const std::string& resultpath,
                            const bool create = false) override
    {
      return false;
    }
    void SetLookUp(const bool lookup) override {}
    void InitializeTheReweighting(ATOOLS::Variations_Mode) override {}
  };

  struct RootNTupleReader_Variables {
#ifdef USING__ROOT
    static const Int_t s_kMaxParticle = 100;
    Int_t m_id;
    Int_t m_ncount, m_nparticle;
    Float_t p_px[s_kMaxParticle];
    Float_t p_py[s_kMaxParticle];
    Float_t p_pz[s_kMaxParticle];
    Float_t p_E[s_kMaxParticle];
    Double_t p_pxd[s_kMaxParticle];
    Double_t p_pyd[s_kMaxParticle];
    Double_t p_pzd[s_kMaxParticle];
    Double_t p_Ed[s_kMaxParticle];
    Int_t p_kf[s_kMaxParticle];

    Double_t m_wgt,m_wgt2,m_pswgt,m_mewgt,m_mewgt2;
    Double_t m_x1,m_x2,m_x1p,m_x2p,m_mur,m_muf,m_as,m_kfac;
    Int_t m_id1,m_id2,m_id1p,m_id2p,m_nuwgt;
    Double_t p_uwgt[18];
    Short_t m_oqcd;
    Char_t m_type[2], m_coqcd;
    TChain* p_f;
#endif
  };
}

RootNtuple_Reader::RootNtuple_Reader(const Input_Arguments &args,int exact,int ftype) :
  Event_Reader_Base(args), m_ftype(ftype), m_otype(0),
  m_evtid(0), m_subevtid(0), m_evtcnt(0), m_entries(0), m_evtpos(0),
  p_isr(args.p_isr), p_yfs(args.p_yfs), m_sargs(NULL,"",""), m_kargs(""), m_xf1(0.), m_xf2(0.)
{
  Settings& s = Settings::GetMainSettings();
  Common_Root_Settings().RegisterDefaults();
  std::string filename=m_path+m_file+".root";
  msg_Out()<<" Reading from "<<filename<<"\n";
  m_ecms = s["ROOTNTUPLE_ECMS"].Get<double>();
  m_calc = s["ROOTNTUPLE_CALC"].Get<int>();
  if (m_calc) msg_Info()<<METHOD<<"(): Ntuple calc mode set to "<<m_calc<<"."<<std::endl;
  m_treename = s["ROOTNTUPLE_TREENAME"].Get<std::string>();
  m_check = s["ROOTNTUPLE_CHECK"].Get<int>();
  if (m_check) msg_Info()<<METHOD<<"(): Ntuple check mode set to "<<m_check<<"."<<std::endl;
  std::string scale = s["SCALES"].Get<std::string>();
  m_lomode = s["ROOTNTUPLE_LO_MODE"].Get<int>();
  if (m_lomode) msg_Info()<<METHOD<<"(): Ntuple LO mode set to "<<m_lomode<<"."<<std::endl;
  std::string kfactor = s["KFACTOR"].Get<std::string>();
  const std::vector<std::string> helpsv{
    s["COUPLINGS"].GetVector<std::string>() };
  std::string coupling(helpsv.size()?helpsv[0]:"");
  for (size_t i(1);i<helpsv.size();++i) coupling+=" "+helpsv[i];
  m_sargs=Scale_Setter_Arguments(args.p_model,scale,coupling);
  m_sargs.m_nin=2;
  m_kargs=KFactor_Setter_Arguments(kfactor);
#ifdef USING__ROOT
  p_vars = new RootNTupleReader_Variables();
  p_vars->p_f=new TChain(m_treename.c_str());
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
    int size=mpi->Size();
    int rank=mpi->Rank();
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
	  mpi->Send(&values,2,MPI_INT,tag,tag);
	  mpi->Recv(&values,2,MPI_INT,MPI_ANY_SOURCE,size+tag);
	  msg_Info()<<"  Rank "<<tag<<" analyzes "
		    <<basename<<"["<<values[0]<<"-"<<values[1]<<"].\n";
	}
	msg_Info()<<"}\n";
      }
      else {
	mpi->Recv(&values,2,MPI_INT,0,rank);
	i=values[0];
	e=values[1];
	mpi->Send(&values,2,MPI_INT,0,size+rank);
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
  msg_Info()<<METHOD<<"(): Found "<<m_entries<<" entries."<<std::endl;
  if (s["ROOTNTUPLE_SET_NEVT"].Get<bool>())
    rpa->gen.SetNumberOfEvents(m_entries);
  p_vars->p_f->SetBranchAddress("id",&p_vars->m_id);
  p_vars->m_ncount=1.0;
  if (p_vars->p_f->GetBranch("ncount")) {
    p_vars->p_f->SetBranchAddress("ncount",&p_vars->m_ncount);
    msg_Info()<<METHOD<<"(): Using exact event recovery."<<std::endl;
  }
  p_vars->p_f->SetBranchAddress("nparticle",&p_vars->m_nparticle);
  m_ftype=p_vars->p_f->GetLeaf("E")->GetTypeName()[0]=='D';
  if (m_ftype) {
    msg_Info()<<METHOD<<"(): Found float type double."<<std::endl;
    p_vars->p_f->SetBranchAddress("px",p_vars->p_pxd);
    p_vars->p_f->SetBranchAddress("py",p_vars->p_pyd);
    p_vars->p_f->SetBranchAddress("pz",p_vars->p_pzd);
    p_vars->p_f->SetBranchAddress("E",p_vars->p_Ed);
  }
  else {
    p_vars->p_f->SetBranchAddress("px",p_vars->p_px);
    p_vars->p_f->SetBranchAddress("py",p_vars->p_py);
    p_vars->p_f->SetBranchAddress("pz",p_vars->p_pz);
    p_vars->p_f->SetBranchAddress("E",p_vars->p_E);
  }

  p_vars->p_f->SetBranchAddress("kf",p_vars->p_kf);
  p_vars->p_f->SetBranchAddress("weight",&p_vars->m_wgt);
  p_vars->p_f->SetBranchAddress("weight2",&p_vars->m_wgt2);
  p_vars->p_f->SetBranchAddress("alphas",&p_vars->m_as);
  p_vars->m_kfac=1.0;
  if (p_vars->p_f->GetBranch("kfactor"))
    p_vars->p_f->SetBranchAddress("kfactor",&p_vars->m_kfac);
  p_vars->m_pswgt=0.0;
  if (p_vars->p_f->GetBranch("ps_wgt"))
    p_vars->p_f->SetBranchAddress("ps_wgt",&p_vars->m_pswgt);
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
  p_vars->m_id1p=0;
  if (p_vars->p_f->GetBranch("id1p"))
    p_vars->p_f->SetBranchAddress("id1p",&p_vars->m_id1p);
  p_vars->m_id2p=0;
  if (p_vars->p_f->GetBranch("id2p"))
    p_vars->p_f->SetBranchAddress("id2p",&p_vars->m_id2p);
  p_vars->p_f->SetBranchAddress("nuwgt",&p_vars->m_nuwgt);
  p_vars->p_f->SetBranchAddress("usr_wgts",p_vars->p_uwgt);
  m_otype=p_vars->p_f->GetLeaf("alphasPower")->GetTypeName()[0]=='C';
  if (m_otype) {
    msg_Info()<<METHOD<<"(): Found order QCD type char."<<std::endl;
    p_vars->p_f->SetBranchAddress("alphasPower",&p_vars->m_coqcd);
  }
  else {
    p_vars->p_f->SetBranchAddress("alphasPower",&p_vars->m_oqcd);
  }
  p_vars->p_f->SetBranchAddress("part",p_vars->m_type);
#else
  msg_Error()<<"Sherpa must be linked with root to read in root files!"<<endl;
#endif
  s["INTRINSIC_KPERP"].OverrideScalar<bool>(false);
}

RootNtuple_Reader::~RootNtuple_Reader()
{
  for (std::map<RR_Process_Info,PHASIC::Process_Base*>::iterator
	 sit(m_procs.begin());sit!=m_procs.end();++sit) delete sit->second;
}

void RootNtuple_Reader::RegisterDefaults() const
{
  Settings& s = Settings::GetMainSettings();
  s["ROOTNTUPLE_ECMS"].SetDefault(rpa->gen.Ecms());
  s["ROOTNTUPLE_CALC"].SetDefault(1);
  const bool calc{ s["ROOTNTUPLE_CALC"].Get<bool>() };
  s["ROOTNTUPLE_CHECK"].SetDefault((calc & 2) ? 1 : 0);
  if (!s["SCALES"].IsSetExplicitly()) {
    s["SCALES"].OverrideScalar("VAR{sqr(" + ToString(rpa->gen.Ecms()) + ")}");
  }
  s["ROOTNTUPLE_LO_MODE"].SetDefault(0);
  s["ROOTNTUPLE_SET_NEVT"].SetDefault(false);
}

void RootNtuple_Reader::CloseFile() {
#ifdef USING__ROOT
  delete p_vars->p_f->GetCurrentFile();
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
  if (m_otype) p_vars->m_oqcd=p_vars->m_coqcd;
#endif
  return 1;
}

double RootNtuple_Reader::CalculateWeight
(const RootNtuple_Reader::Weight_Calculation_Args &args, MODEL::One_Running_AlphaS *as)
{
  const double mur2(args.m_mur2);
  const double muf2(args.m_muf2);
  const int mode(args.m_mode);
#ifdef USING__ROOT
  Flavour fl1((kf_code)abs(p_vars->m_id1),p_vars->m_id1<0);
  Flavour fl2((kf_code)abs(p_vars->m_id2),p_vars->m_id2<0);
  double sf(m_ecms/rpa->gen.Ecms());
  p_isr->PDF(0)->Calculate(p_vars->m_x1*sf,muf2);
  p_isr->PDF(1)->Calculate(p_vars->m_x2*sf,muf2);
  m_xf1=p_isr->PDF(0)->GetXPDF(fl1);
  m_xf2=p_isr->PDF(1)->GetXPDF(fl2);
  double fa=m_xf1/p_vars->m_x1;
  double fb=m_xf2/p_vars->m_x2;
  double asf=pow((*as)(mur2)/p_vars->m_as,p_vars->m_oqcd);
  MINLO_Scale_Setter *minlo(dynamic_cast<MINLO_Scale_Setter*>(args.p_scale));
  if (minlo && minlo->Amplitudes().size()) {
    asf=1.0;
    Cluster_Amplitude *ampl(minlo->Amplitudes().back());
    double moqcd(0);
    if (p_vars->m_type[0]=='R') {
      if (ampl->Next()) ampl=ampl->Next();
      else moqcd=1;
    }
    for (;ampl->Next();ampl=ampl->Next()) {
      double oqcd=ampl->OrderQCD()-ampl->Next()->OrderQCD();
      double casf=pow((*as)(ampl->KT2()*minlo->RSF()*args.m_mur2f)/p_vars->m_as,oqcd);
      asf*=casf;
#ifdef DEBUG__MINLO
      msg_Debugging()<<"DEBUG MINLO   local \\alpha_s weight "<<casf<<"  <-  ( "
		     <<casf*p_vars->m_as<<" / "<<p_vars->m_as<<" ) ^ "<<oqcd
		     <<"  <-  ( k_T = "<<sqrt(minlo->RSF()*args.m_mur2f)
		     <<" * "<<sqrt(ampl->KT2())<<" ) \n";
#endif
    }
    int oqcd(ampl->OrderQCD()-moqcd);
    if (p_vars->m_type[0]=='V' || p_vars->m_type[0]=='I') oqcd-=1;
    if (oqcd>0) {
      double casf=pow((*as)(ampl->KT2()*minlo->RSF()*args.m_mur2f)/p_vars->m_as,oqcd);
      asf*=casf;
#ifdef DEBUG__MINLO
      msg_Debugging()<<"DEBUG MINLO   local \\alpha_s weight "<<casf<<"  <-  ( "
		     <<casf*p_vars->m_as<<" / "<<p_vars->m_as<<" ) ^ "<<oqcd
		     <<"  <-  ( k_T = "<<sqrt(minlo->RSF()*args.m_mur2f)
		     <<" * "<<sqrt(ampl->KT2())<<" )\n";
#endif
    }
    if (p_vars->m_type[0]!='B' && !m_lomode) {
      double casf=(*as)(minlo->MuRAvg(1))/p_vars->m_as;
      asf*=casf;
#ifdef DEBUG__MINLO
      msg_Debugging()<<"DEBUG MINLO   nlo \\alpha_s weight "<<casf<<"  <-  "
		     <<casf*p_vars->m_as<<" / "<<p_vars->m_as<<"  <-  ( k_T = "
		     <<sqrt(minlo->MuRAvg(1))<<" )\n";
#endif
    }
#ifdef DEBUG__MINLO
    msg_Debugging()<<"DEBUG MINLO   full \\alpha_s weight "<<asf<<"\n";
  }
  else {
    msg_Debugging()<<"DEBUG MINLO   \\alpha_s weight "<<asf<<"  <-  ( "
		   <<(*as)(mur2)<<" / "<<p_vars->m_as<<" ) ^ "<<p_vars->m_oqcd<<"\n";
#endif
  }
#ifdef DEBUG__MINLO
  msg_Debugging()<<"DEBUG MINLO   \\mu_F = "<<sqrt(muf2)<<", \\mu_R = "<<sqrt(mur2)<<"\n";
#endif
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
    bool wnz=false;
    for (int i(1);i<9;++i) {
      w[i]=p_vars->m_nuwgt<=2?0:p_vars->p_uwgt[i+1]+p_vars->p_uwgt[i+9]*lf;
      if (w[i]) wnz=true;
    }
    double wgt=w[0]*fa*fb;
    if (wnz) {
    if (sf!=1.0) THROW(not_implemented,"I-term rescaling not supported.");
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
    }
    return wgt*asf;
  }
  else {
    THROW(not_implemented,"Invalid option");
  }
#else
  return -1.0;
#endif
}

double RootNtuple_Reader::CalculateWeight(const Weight_Calculation_Args& args,
                                          const QCD_Variation_Params& varparams)
{
  DEBUG_FUNC("R = " << sqrt(varparams.m_muR2fac)
                    << ", F = " << sqrt(varparams.m_muF2fac));
  // temporarily replace PDFs
  PDF::PDF_Base *nominalpdf1 = p_isr->PDF(0);
  PDF::PDF_Base *nominalpdf2 = p_isr->PDF(1);
  p_isr->SetPDF(varparams.p_pdf1, 0);
  p_isr->SetPDF(varparams.p_pdf2, 1);

  Weight_Calculation_Args varargs(args.m_mur2 * varparams.m_muR2fac,
                                  args.m_muf2 * varparams.m_muF2fac,
                                  args.m_mode,
                                  args.p_scale,
                                  args.p_kfac,
                                  args.m_K,
                                  varparams.m_muR2fac);
  if (args.p_scale && args.p_scale->UpdateScale(varparams)) {
    varargs.m_mur2=args.p_scale->Scale(stp::ren);
    varargs.m_muf2=args.p_scale->Scale(stp::fac);
  }
  double weight(CalculateWeight(varargs, varparams.p_alphas->GetAs()));
  if (args.p_kfac && args.p_kfac->UpdateKFactor(varparams))
    weight*=args.p_kfac->LastKFactor()/args.m_K;

  // reset PDFs
  p_isr->SetPDF(nominalpdf1, 0);
  p_isr->SetPDF(nominalpdf2, 1);
  p_isr->SetMuF2(args.m_muf2, 0);
  p_isr->SetMuF2(args.m_muf2, 1);

  return weight;
}

bool RootNtuple_Reader::ReadInFullEvent(Blob_List * blobs)
{
  m_mewgtinfo.Reset();
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
  signalblob->SetTypeSpec("NLO");
  signalblob->SetId();
  signalblob->SetPosition(Vec4D(0.,0.,0.,0.));
  signalblob->SetStatus(blob_status::code(30));
#ifdef USING__ROOT
  size_t currentid=m_evtid;
  int id1(0), id2(0);
  double x1(1.), x2(1.);
  double muR2(0.), muF2(0.);
  Vec4D bm[2]={rpa->gen.PBeam(0),rpa->gen.PBeam(1)};
  for (int i(0);i<2;++i) bm[i][3]=bm[i][3]<0.0?-bm[i][0]:bm[i][0];
  RootNTupleReader_Variables vars;
  while (currentid==m_evtid) {
    Vec4D sum;
    Vec4D *moms = new Vec4D[2+p_vars->m_nparticle];
    Flavour *flav = new Flavour[2+p_vars->m_nparticle];
    id1=p_vars->m_id1p?p_vars->m_id1p:p_vars->m_id1;
    id2=p_vars->m_id2p?p_vars->m_id2p:p_vars->m_id2;
    flav[0]=Flavour((long int)id1);
    flav[1]=Flavour((long int)id2);
    for (int i=0;i<p_vars->m_nparticle;i++) {
      if (m_ftype)
	moms[i+2]=Vec4D(p_vars->p_Ed[i],p_vars->p_pxd[i],
			p_vars->p_pyd[i],p_vars->p_pzd[i]);
      else
	moms[i+2]=Vec4D(p_vars->p_E[i],p_vars->p_px[i],
			p_vars->p_py[i],p_vars->p_pz[i]);
      flav[i+2]=Flavour(abs(p_vars->p_kf[i]),p_vars->p_kf[i]<0);
      sum+=moms[i+2];
    }
    m_nlos.push_back(new NLO_subevt(p_vars->m_nparticle+2,NULL,flav,moms));
    m_nlos.back()->m_result=p_vars->m_wgt2;
    m_nlos.back()->m_mu2[stp::fac]=sqr(p_vars->m_muf);
    m_nlos.back()->m_mu2[stp::ren]=sqr(p_vars->m_mur);
    m_nlos.back()->m_stype=sbt::qcd;
    // double sf(m_ecms/rpa->gen.Ecms());
    // x1=p_vars->m_x1*sf;
    // x2=p_vars->m_x2*sf;
    x1=sum.PPlus()/rpa->gen.PBeam(0).PPlus();
    x2=sum.PMinus()/rpa->gen.PBeam(1).PMinus();
    moms[0]=x1*bm[0];
    moms[1]=x2*bm[1];
    if (m_calc) {
      Vec4D_Vector p(moms,&moms[p_vars->m_nparticle+2]);
      RR_Process_Info info(p_vars->m_type,p_vars->m_nparticle+2,flav);
      if (m_procs.find(info)==m_procs.end()) {
	Process_Info pi;
	if (p_vars->m_type[0]!='B' || !m_lomode)
	  pi.m_fi.m_nlotype=ToType<nlo_type::code>(p_vars->m_type);
	pi.m_ii.m_ps.push_back(Subprocess_Info(flav[0]));
	pi.m_ii.m_ps.push_back(Subprocess_Info(flav[1]));
	pi.m_maxcpl[0]=pi.m_mincpl[0]=p_vars->m_oqcd;
	for (int i=0;i<p_vars->m_nparticle;++i)
	  pi.m_fi.m_ps.push_back(Subprocess_Info(flav[i+2]));
	m_sargs.p_proc=m_procs[info] = new Dummy_Process();
	m_sargs.p_proc->Init(pi,NULL,p_isr,p_yfs,1);
	m_sargs.p_proc->SetMaxOrder
	  (0,p_vars->m_oqcd-(p_vars->m_type[0]=='S'?1:0));
	m_sargs.m_nout=p_vars->m_nparticle;
	m_sargs.p_proc->SetScale(m_sargs);
	m_sargs.p_proc->SetKFactor(m_kargs);
      }
      signalblob->SetTypeSpec(m_procs[info]->Name());
      Scale_Setter_Base *scale(&*m_procs[info]->ScaleSetter());
      KFactor_Setter_Base *kfac(&*m_procs[info]->KFactorSetter());
      scale->CalculateScale(p);
      muR2=m_nlos.back()->m_mu2[stp::ren]=scale->Scale(stp::ren);
      muF2=m_nlos.back()->m_mu2[stp::fac]=scale->Scale(stp::fac);
      const double K(kfac->KFactor(p_vars->m_type[0]=='B'?1:0));
      const Weight_Calculation_Args args(muR2,muF2,p_vars->m_nuwgt?1:2,scale,kfac,K,1.);
      double weight=CalculateWeight(args, MODEL::as->GetAs());
      weight*=K/p_vars->m_kfac;
      m_nlos.back()->m_results = weight;
      ATOOLS::Reweight(
          m_nlos.back()->m_results["Main"],
          [this, &args, K, weight](double varweight,
                           const QCD_Variation_Params& varparams) -> double {
            varweight = CalculateWeight(args, varparams);
            return varweight / weight * K / p_vars->m_kfac;
          });
      m_nlos.back()->m_results["All"] = m_nlos.back()->m_results["Main"];

#ifdef DEBUG__MINLO
      msg_Debugging()<<"DEBUG MINLO   total weight "<<weight/p_vars->m_wgt2<<"\n";
#endif
      m_nlos.back()->m_result=weight;
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
	(Weight_Calculation_Args(sqr(p_vars->m_mur),sqr(p_vars->m_muf),
				 p_vars->m_nuwgt?1:2,NULL,NULL,1.0,1.0),
         MODEL::as->GetAs());
      RR_Process_Info info(p_vars->m_type,p_vars->m_nparticle+2,flav);
      KFactor_Setter_Base *kfac(&*m_procs[info]->KFactorSetter());
      weight*=kfac->KFactor(p_vars->m_type[0]=='B'?1:0)/p_vars->m_kfac;
      msg_Debugging()<<METHOD<<"(): "<<p_vars->m_type<<" computed "
		     <<weight<<", stored "<<p_vars->m_wgt2
		     <<", rel. diff. "<<weight/p_vars->m_wgt2-1.0<<".\n";
      if (!IsEqual(weight,p_vars->m_wgt2,rpa->gen.Accu()))
	msg_Error()<<METHOD<<"(): "<<p_vars->m_type<<" weights differ by "
		   <<(weight/p_vars->m_wgt2-1.0)<<".\n  computed "
		   <<weight<<", stored "<<p_vars->m_wgt2<<"."<<std::endl;
    }
    vars=*p_vars;
    if (!ReadInEntry()) m_evtid=0;
  }
  Particle *part1=new Particle(0,m_nlos.back()->p_fl[0],x1*bm[0]);
  Particle *part2=new Particle(1,m_nlos.back()->p_fl[1],x2*bm[1]);
  signalblob->AddToInParticles(part1);
  signalblob->AddToInParticles(part2);
  for (size_t i=2;i<m_nlos.back()->m_n;++i) {
    Particle *part=new Particle
      (i,m_nlos.back()->p_fl[i],m_nlos.back()->p_mom[i]);
    signalblob->AddToOutParticles(part);
  }
#endif
  m_pdfinfo=PDF_Info((long int)vars.m_id1,(long int)vars.m_id2,
		     x1,x2,muF2,muF2,m_xf1,m_xf2);
  // only reliable in SM,
  // HEFT breaks this, counting Yukawa's separately breaks this, etc.
  bool onemoreas(vars.m_type[0]=='V' || vars.m_type[0]=='I' ||
                 vars.m_type[0]=='S');
  int oew(m_nlos.back()->m_n-2+(onemoreas?1:0)-vars.m_oqcd);
  signalblob->SetStatus(blob_status::needs_beams);
  signalblob->AddStatus(blob_status::needs_harddecays);
  Weights_Map wgtmap {0.0};
  for (const auto* sub : m_nlos) {
    wgtmap += sub->m_results;
  }
  signalblob->AddData("WeightsMap",new Blob_Data<Weights_Map>(wgtmap));
  signalblob->AddData("MEWeight",new Blob_Data<double>
		      ((vars.m_nuwgt?vars.m_mewgt:vars.m_mewgt2)/
		       (vars.m_pswgt?vars.m_pswgt:1.0)));
  if (vars.m_pswgt)
    signalblob->AddData("PSWeight",new Blob_Data<double>(vars.m_pswgt));
  signalblob->AddData("Trials",new Blob_Data<double>(vars.m_ncount));
  if (vars.m_type[0]=='S' || vars.m_type[0]=='R')
    signalblob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(&m_nlos));
  signalblob->AddData("Weight_Norm",new Blob_Data<double>(1.0));
  signalblob->AddData("OQCD",new Blob_Data<int>(vars.m_oqcd));
  signalblob->AddData("OEW",new Blob_Data<int>(oew));
  signalblob->AddData("Renormalization_Scale", new Blob_Data<double>(muR2));
  signalblob->AddData("Factorization_Scale", new Blob_Data<double>(muF2));
  signalblob->AddData("PDFInfo", new Blob_Data<ATOOLS::PDF_Info>(m_pdfinfo));
  signalblob->AddData("MEWeightInfo",new Blob_Data<ATOOLS::ME_Weight_Info*>(&m_mewgtinfo));
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

DECLARE_GETTER(ERootNtuple_Reader,"ERoot",
	       Event_Reader_Base,Input_Arguments);

Event_Reader_Base *ATOOLS::Getter
<Event_Reader_Base,Input_Arguments,ERootNtuple_Reader>::
operator()(const Input_Arguments &args) const
{
  return new RootNtuple_Reader(args,1);
}

void ATOOLS::Getter
<Event_Reader_Base,Input_Arguments,ERootNtuple_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Root NTuple input";
}

DECLARE_GETTER(EDRootNtuple_Reader,"EDRoot",
	       Event_Reader_Base,Input_Arguments);

Event_Reader_Base *ATOOLS::Getter
<Event_Reader_Base,Input_Arguments,EDRootNtuple_Reader>::
operator()(const Input_Arguments &args) const
{
  return new RootNtuple_Reader(args,1,1);
}

void ATOOLS::Getter
<Event_Reader_Base,Input_Arguments,EDRootNtuple_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Root NTuple input";
}

