#include "PHASIC++/Process/Single_Process.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "PHASIC++/Channels/CS_Dipole.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "METOOLS/Explicit/Collinear_Splitting_Tools.H"
#include "MODEL/Main/Coupling_Data.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Single_Process::Single_Process(): m_lastbxs(0.0), m_zero(false)
{
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  m_nloct=read.GetValue<int>("SP_NLOCT",m_pinfo.m_ckkw&1);
  if (!m_nloct && m_pinfo.m_ckkw&1) {
    static bool print(false);
    if (!print) msg_Info()
      <<METHOD<<"(): Switch off NLO counterterms.\n";
    print=true;
  }
}

Single_Process::~Single_Process()
{
  for (Coupling_Map::const_iterator
	 cit(m_cpls.begin());cit!=m_cpls.end();++cit)
    delete cit->second;
}

size_t Single_Process::Size() const
{
  return 1;
}

Process_Base *Single_Process::operator[](const size_t &i)
{
  if (i==0) return this;
  return NULL;
}

Weight_Info *Single_Process::OneEvent(const int wmode,const int mode)
{
  p_selected=this;
  return p_int->PSHandler()->OneEvent(this,mode);
}

double Single_Process::KFactor() const
{
  if (p_kfactor) return p_kfactor->KFactor();
  return 1.0;
}

namespace PHASIC {

  double Beta0()
  {
    return 11.0/6.0*3.0-2.0/3.0*0.5*(Flavour(kf_jet).Size()/2);
  }

}

double Single_Process::NLOCounterTerms() const
{
  static double th(1.0e-12);
  if (!m_use_biweight) return 0.0;
  DEBUG_FUNC(m_name);
  if (p_scale->Scales().size()<stp::size+2)
    THROW(fatal_error,"Invalid number of scales: "+
	  ToString(p_scale->Scales().size())+
	  " < "+ToString(stp::size+2));
  double lmuf2(p_scale->Scale(stp::id(stp::fac)));
  double lmur2(p_scale->Scale(stp::id(stp::ren)));
  double muf2(p_scale->Scale(stp::id(stp::size+stp::fac)));
  double mur2(p_scale->Scale(stp::id(stp::size+stp::ren)));
  msg_Debugging()<<"\\mu_F = "<<sqrt(muf2)<<" -> "<<sqrt(muf2/lmuf2)<<"\n";
  msg_Debugging()<<"\\mu_R = "<<sqrt(mur2)<<" -> "<<sqrt(mur2/lmur2)<<"\n";
  MODEL::Coupling_Data *cpl(m_cpls.Get("Alpha_QCD"));
  double as(cpl->Default()*cpl->Factor());
  double ct(0.0);
  if (!IsEqual(mur2,lmur2)) {
    if (m_oqcd>1) ct-=double(m_oqcd-1)*as/M_PI*MODEL::as->Beta0(lmur2)
                                      *log(mur2/lmur2);
    msg_Debugging()<<"\\alpha_s term (nf="<<MODEL::as->Nf(lmur2)<<"): "
                   <<(m_oqcd-1)<<" * "<<as
                   <<"/2\\pi * "<<MODEL::as->Beta0(lmur2)
                   <<" * log("<<sqrt(mur2)<<"/"
                   <<sqrt(lmur2)<<") -> "<<ct<<"\n";
  }
  double z[2]={m_mewgtinfo.m_y1,m_mewgtinfo.m_y2};
  for (size_t i(0);i<2;++i)
    ct+=CollinearCounterTerms
      (i,m_flavs[i],p_int->Momenta()[i],z[i],muf2,lmuf2);
  msg_Debugging()<<"C = "<<ct<<"\n";
  return ct;
}

double Single_Process::CollinearCounterTerms
(const int i,const Flavour &fl,const Vec4D &p,
 const double &z,const double &t1,const double &t2) const
{
  if (!(p_int->ISR() && p_int->ISR()->On()&(1<<i))) return 0.0;
  static double th(1.0e-12);
  DEBUG_FUNC("Q = "<<sqrt(t1)<<" / "<<sqrt(t2));
  if (IsEqual(t1,t2)) return 0.0;
  double lmuf2(p_scale->Scale(stp::fac));
  msg_Debugging()<<"\\mu_F = "<<sqrt(lmuf2)<<"\n";
  msg_Debugging()<<"\\mu_R = "<<sqrt(p_scale->Scale(stp::ren))<<"\n";
  MODEL::Coupling_Data *cpl(m_cpls.Get("Alpha_QCD"));
  double as(cpl->Default()*cpl->Factor());
  double ct(0.0), lt(log(t1/t2)), x(p_int->ISR()->CalcX(p));
  msg_Debugging()<<as<<"/(2\\pi) * log("<<sqrt(t1)<<"/"
		 <<sqrt(t2)<<") = "<<as/(2.0*M_PI)*lt<<"\n";
  Flavour jet(kf_jet);
  double fb=p_int->ISR()->PDFWeight(1<<(i+1),p,p,lmuf2,lmuf2,fl,fl,0);
  if (IsZero(fb,th)) {
    msg_Tracking()<<METHOD<<"(): Zero xPDF ( f_{"<<fl<<"}("
		  <<x<<","<<sqrt(lmuf2)<<") = "<<fb<<" ). Skip.\n";
    return 0.0;
  }
  msg_Debugging()<<"Beam "<<i<<": z = "<<z<<", f_{"<<fl
		 <<"}("<<x<<","<<sqrt(lmuf2)<<") = "<<fb<<" {\n";
  for (size_t j(0);j<jet.Size();++j) {
    double Pf(METOOLS::FPab(jet[j],fl,z));
    double Ps(METOOLS::SPab(jet[j],fl,z));
    if (Pf+Ps==0.0) continue;
    double Pi(METOOLS::IPab(jet[j],fl,x));
    double H(METOOLS::Hab(jet[j],fl));
    double fa=p_int->ISR()->PDFWeight
      (1<<(i+1),p/z,p/z,lmuf2,lmuf2,jet[j],jet[j],0);
    double fc=p_int->ISR()->PDFWeight
      (1<<(i+1),p,p,lmuf2,lmuf2,jet[j],jet[j],0);
    msg_Debugging()<<"  P_{"<<jet[j]<<","<<fl
		   <<"}("<<z<<") = {F="<<Pf<<",S="<<Ps
		   <<",I="<<Pi<<"}, f_{"<<jet[j]<<"}("
		   <<x/z<<","<<sqrt(lmuf2)<<") = "<<fa
		   <<", f_{"<<jet[j]<<"}("<<x<<","
		   <<sqrt(lmuf2)<<") = "<<fc<<"\n";
    if (IsZero(fa,th)||IsZero(fc,th)) {
      msg_Tracking()<<METHOD<<"(): Zero xPDF. Skip.\n";
      return 0.0;
    }
    ct+=as/(2.0*M_PI)*lt*
      ((fa/z*Pf+(fa/z-fc)*Ps)*(1.0-x)+fc*(H-Pi))/fb;
  }
  msg_Debugging()<<"} -> "<<ct<<"\n";
  return ct;
}

ATOOLS::Cluster_Sequence_Info Single_Process::BeamISRWeight
(const double& Q2,const int imode,
 const ClusterAmplitude_Vector &ampls) const
{
  int mode(imode&1);
  if (mode) msg_Out()<<"Flipped initial states.\n";
  if (!m_use_biweight) return 1.;
  if (m_nin==1) return 0.5/p_int->Momenta()[0].Mass();
  else if (m_nin>2) THROW(not_implemented,"More than two incoming particles.");
  Cluster_Sequence_Info csi;
  if (p_int->ISR()) {
    // external PDFs contain flux
    double pdfext(p_int->ISR()->PDFWeight(mode,p_int->Momenta()[0],
                                               p_int->Momenta()[1],
                                               Q2,Q2,m_flavs[0],m_flavs[1]));
    msg_Debugging()<<"PDF(fla="<<m_flavs[0]
                   <<", xa="<<p_int->ISR()->CalcX(p_int->Momenta()[0])
                   <<", ta="<<Q2<<") * "
                   <<"PDF(flb="<<m_flavs[1]
                   <<", xb="<<p_int->ISR()->CalcX(p_int->Momenta()[1])
                   <<", tb="<<Q2<<") -> "<<pdfext<<std::endl;
    csi.AddWeight(pdfext);
    if (ampls.size() && (m_pinfo.m_ckkw&1)) {
      DEBUG_FUNC(m_name<<", \\mu_F = "<<sqrt(Q2)<<", mode = "<<mode
                 <<", #ampls="<<ampls.size());
      // add outer splitting
      csi.AddSplitting(Q2,p_int->ISR()->CalcX(p_int->Momenta()[mode]),
                          p_int->ISR()->CalcX(p_int->Momenta()[1-mode]),
                          m_flavs[mode],m_flavs[1-mode]);
      csi.AddPDFRatio(pdfext,1.);
      Cluster_Amplitude *ampl(ampls.front());
      msg_IODebugging()<<*ampl<<"\n";
      if (imode&2) {
	ampl=ampl->Next();
	msg_IODebugging()<<*ampl<<"\n";
      }
      int set(false);
      double LQ2(Q2);
      double pdfnum(pdfext), pdfden(pdfext);
      for (;ampl;ampl=ampl->Next()) {
        // skip: decays, equal scales, unordered configs, quarks below threshold
        // add into ClusterSteps everything but decays
        msg_IODebugging()<<*ampl<<"\n";
	if (ampl->Next()) {
	  if (ampl->Next()->Splitter()->Stat()==3) {
	    msg_Debugging()<<"Skip. Decay "<<
	      ID(ampl->Next()->Splitter()->Id())<<"\n";
	    continue;
	  }
	}
	Flavour f1(ampl->Leg(0)->Flav().Bar());
	Flavour f2(ampl->Leg(1)->Flav().Bar());
	if (MapProc() && LookUp() && !(imode&2)) {
	  f1=ReMap(f1,ampl->Leg(0)->Id());
	  f2=ReMap(f2,ampl->Leg(1)->Id());
	}
	csi.AddSplitting(ampl->KT2(),
			 p_int->ISR()->CalcX(-ampl->Leg(0)->Mom()),
			 p_int->ISR()->CalcX(-ampl->Leg(1)->Mom()),
			 f1,f2);
	if (IsEqual(LQ2,ampl->Next()?ampl->KT2():ampl->MuF2())) {
	  msg_Debugging()<<"Skip. Scales equal: t_i="<<LQ2
			 <<", t_{i+1}="<<(ampl->Next()?ampl->KT2():ampl->MuF2())
			 <<std::endl;
	  if (ampl->Next()!=NULL) csi.AddPDFRatio(pdfnum,pdfden);
	  else                    csi.AddPDFRatio(1.,pdfden);
	  continue;
	}
	if (set && LQ2>ampl->KT2()) {
	  msg_Debugging()<<"Skip. Unordered history "<<
	    sqrt(LQ2)<<" > "<<sqrt(ampl->KT2())<<"\n";
	  LQ2=sqrt(std::numeric_limits<double>::max());
	  continue;
	}
	if (LQ2<sqr(2.0*f1.Mass(true)) || LQ2<sqr(2.0*f2.Mass(true))) {
	  msg_Debugging()<<"Skip. Quarks below threshold: t="<<LQ2
			 <<" vs. "<<sqr(2.0*f1.Mass(true))
			 <<" / "<<sqr(2.0*f2.Mass(true))<<std::endl;
	  continue;
	}
	// denominators
	double wd1=p_int->ISR()->PDFWeight
	  (mode|2,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	double wd2=p_int->ISR()->PDFWeight
	  (mode|4,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	double LLQ2=LQ2;
	LQ2=ampl->KT2();
	if (ampl->Next()==NULL) LQ2=ampl->MuF2();
	// numerators
	double wn1=p_int->ISR()->PDFWeight
	  (mode|2,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	double wn2=p_int->ISR()->PDFWeight
	  (mode|4,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	if (!IsZero(wn1) && !IsZero(wd1)) csi.AddWeight(wn1/wd1);
	if (!IsZero(wn2) && !IsZero(wd2)) csi.AddWeight(wn2/wd2);
	// book-keep PDF ratios excl.
	//   a) first one correcting outer PDF from muF to t
	//   b) last numerator taken at muF (this one is to be varied)
	// use the following identity with i=0 -> core and i=N -> ext
	// wn-ext * [\prod_{i=0}^{N-1} wn_i/wd_i]
	// = [wn-ext * \prod_{i=1}^{N-1} wn_i/wd_i * 1/wd_0] * wn-core
	// = [\prod_{i=1}^N wn_i/wd_{i-1}] * wn-core
	pdfnum=wn1*wn2;
	pdfden=wd1*wd2;
	if (ampl->Next()!=NULL) csi.AddPDFRatio(pdfnum,pdfden);
	else                    csi.AddPDFRatio(1.,pdfden);
	msg_Debugging()<<"* [  "
		       <<"PDF(fla="<<f1
		       <<", xa="<<p_int->ISR()->CalcX(-ampl->Leg(0)->Mom())
		       <<", ta="<<LQ2<<") * "
		       <<"PDF(flb="<<f2
		       <<", xb="<<p_int->ISR()->CalcX(-ampl->Leg(1)->Mom())
		       <<", tb="<<LQ2<<") -> "<<wn1*wn2<<"\n"
		       <<"   / "
		       <<"PDF(fla="<<f1
		       <<", xa="<<p_int->ISR()->CalcX(-ampl->Leg(0)->Mom())
		       <<", ta="<<LLQ2<<") * "
		       <<"PDF(flb="<<f2
		       <<", xb="<<p_int->ISR()->CalcX(-ampl->Leg(1)->Mom())
		       <<", tb="<<LLQ2<<") -> "<<wd1*wd2
		       <<" ] = "<<wn1*wn2/wd1/wd2<<std::endl;
	if (m_pinfo.Has(nlo_type::born)) {
	  double rn[2]={ran->Get(),ran->Get()};
	  for (int i(0);i<2;++i) {
	    if (i==0 && (IsZero(wn1) || IsZero(wd1))) continue;
	    if (i==1 && (IsZero(wn2) || IsZero(wd2))) continue;
	    Vec4D p(-ampl->Leg(i)->Mom());
	    double x(p_int->ISR()->CalcX(p)), z(x+(1.0-x)*rn[i]);
	    csi.AddCounterTerm(CollinearCounterTerms(i,i?f2:f1,p,z,LQ2,LLQ2));
	  }
	}
	set=true;
      }
    }
  }
  if (p_int->Beam() && p_int->Beam()->On()) {
    p_int->Beam()->CalculateWeight(Q2);
    csi.AddWeight(p_int->Beam()->Weight());
  }
  return csi;
}

void Single_Process::BeamISRWeight
(NLO_subevtlist *const subs,const int mode) const
{
  double muf2(subs->back()->m_mu2[stp::fac]);
  double flux(p_int->ISR()->Flux(p_int->Momenta()[0],p_int->Momenta()[1]));
  if (m_nin==2 && p_int->ISR()) {
    size_t nscales(0);
    for (size_t i(0);i<subs->size();++i) {
      NLO_subevt *sub((*subs)[i]);
      if (sub->m_me==0.0) sub->m_result=0.0;
      if ((!IsEqual(sub->m_mu2[stp::fac],muf2) ||
	   m_pinfo.m_nlomode!=1) && sub->m_me!=0.0) {
	ClusterAmplitude_Vector ampls(sub->p_ampl?1:0,sub->p_ampl);
	if (ampls.size()) ampls.front()->SetProc(sub->p_proc);
	sub->m_result=sub->m_me
		      *BeamISRWeight(sub->m_mu2[stp::fac],mode|2,ampls).m_pdfwgt
		      *flux;
	sub->m_xf1=p_int->ISR()->XF1(0);
	sub->m_xf2=p_int->ISR()->XF2(0);
	++nscales;
      }
    }
    if (nscales<subs->size() && m_pinfo.m_nlomode==1) {
      double lumi(BeamISRWeight(muf2,mode,ClusterAmplitude_Vector()).m_pdfwgt);
      lumi*=flux;
      for (size_t i(0);i<subs->size();++i) {
	if ((*subs)[i]->m_me==0.0) (*subs)[i]->m_result=0.0;
	if (IsEqual((*subs)[i]->m_mu2[stp::fac],muf2) &&
	    (*subs)[i]->m_me!=0.0) {
          (*subs)[i]->m_result=(*subs)[i]->m_me*lumi;
	  (*subs)[i]->m_xf1=p_int->ISR()->XF1(0);
	  (*subs)[i]->m_xf2=p_int->ISR()->XF2(0);
        }
      }
    }
  }
  else {
    for (size_t i(0);i<subs->size();++i) {
      ClusterAmplitude_Vector ampls(1,(*subs)[i]->p_ampl);
      if (ampls.size()) ampls.front()->SetProc((*subs)[i]->p_proc);
      (*subs)[i]->m_result=(*subs)[i]->m_me*
        BeamISRWeight((*subs)[i]->m_mu2[stp::fac],mode,ampls).m_pdfwgt*flux;
    }
  }
}

double Single_Process::Differential(const Vec4D_Vector &p)
{
  DEBUG_FUNC(Name()<<", RS:"<<GetSubevtList());
  m_lastb=m_last=0.0;
  m_mewgtinfo.Reset();
  m_mewgtinfo.m_oqcd=OrderQCD();
  m_mewgtinfo.m_oew=OrderEW();
  m_mewgtinfo.m_fl1=(int)(Flavours()[0]);
  m_mewgtinfo.m_fl2=(int)(Flavours()[1]);
  p_int->SetMomenta(p);
  if (IsMapped()) p_mapproc->Integrator()->SetMomenta(p);
  double flux(p_int->ISR()->Flux(p[0],p[1]));
  if (GetSubevtList()==NULL) {
    if (m_zero) return 0.0;
    Scale_Setter_Base *scs(ScaleSetter(1));
    scs->SetCaller(Proc());
    if (Partonic(p,0)==0.0) return 0.0;
    if (m_mewgtinfo.m_type==mewgttype::none) m_mewgtinfo.m_B = m_lastxs;
    m_mewgtinfo*=flux;
    m_mewgtinfo.m_muf2=scs->Scale(stp::fac);
    m_mewgtinfo.m_mur2=scs->Scale(stp::ren);
    if (m_lastxs==0.0) return m_last=0.0;
    m_last=m_lastxs;
    if (m_nloct && m_pinfo.Has(nlo_type::born))
      m_last+=m_lastbxs*NLOCounterTerms();
    ATOOLS::Cluster_Sequence_Info csi=BeamISRWeight
      (scs->Scale(stp::fac),0,scs->Amplitudes().size()?
       scs->Amplitudes():ClusterAmplitude_Vector());
    csi.AddFlux(flux);
    m_mewgtinfo.m_clusseqinfo=csi;
    msg_Debugging()<<m_mewgtinfo;
    m_last=(m_last-m_lastbxs*csi.m_ct)*csi.m_pdfwgt*csi.m_flux;
    m_lastb=m_lastbxs*csi.m_pdfwgt*csi.m_flux;
    if (p_mc==NULL) return m_last;
    // calculate DADS for MC@NLO, one PS point, many dipoles
    Dipole_Params dps(p_mc->Active(this));
    std::vector<double> x(2,-1.0);
    for (size_t j(0);j<2;++j) x[j]=Min(dps.m_p[j][3]/rpa->gen.PBeam(j)[3],1.);
    for (size_t i(0);i<dps.m_procs.size();++i) {
      Process_Base *cp(dps.m_procs[i]);
      size_t mcmode(cp->SetMCMode(m_mcmode));
      bool lookup(cp->LookUp());
      cp->SetLookUp(false);
      double dadswgt(cp->Differential(dps.m_p)*dps.m_weight);
      double dadsmewgt(cp->GetMEwgtinfo()->m_B*dps.m_weight);
      DADS_Info dads(-dadsmewgt,
                     cp->Flavours()[0],cp->Flavours()[1],
                     x[0],x[1],
                     0.,0.,
                     cp->ScaleSetter()->Scale(stp::ren),
                     cp->ScaleSetter()->Scale(stp::fac));
      msg_Debugging()<<dads<<std::endl;
      m_mewgtinfo.m_dadsinfos.push_back(dads);
      m_last-=dadswgt;
      cp->SetLookUp(lookup);
      cp->SetMCMode(mcmode);
    }
    return m_last;
  }
  Partonic(p,0);
  NLO_subevtlist *subs(GetSubevtList());
  for (size_t i=0;i<subs->size();++i) {
    m_mewgtinfo.m_RS+=(*subs)[i]->m_mewgt;
  }
  m_mewgtinfo*=flux;
  Scale_Setter_Base *scs(ScaleSetter(1));
  if (scs!=NULL) {
    m_mewgtinfo.m_muf2=scs->Scale(stp::fac);
    m_mewgtinfo.m_mur2=scs->Scale(stp::ren);
  }
  BeamISRWeight(subs,0);
  for (size_t i=0;i<subs->size();++i) {
    m_last+=(*subs)[i]->m_result;
    (*subs)[i]->m_mewgt*=flux;
  }
  return m_last;
}

bool Single_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  p_int->Reset();
  SP(Phase_Space_Handler) psh(p_int->PSHandler());
  if (p_int->ISR()) {
    if (m_nin==2) {
      if (m_flavs[0].Mass()!=p_int->ISR()->Flav(0).Mass() ||
          m_flavs[1].Mass()!=p_int->ISR()->Flav(1).Mass()) {
        p_int->ISR()->SetPartonMasses(m_flavs);
      }
    }
  }
  psh->InitCuts();
  if (p_int->ISR())
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  exh->AddTerminatorObject(p_int);
  psh->InitIncoming();
  double var(p_int->TotalVar());
  msg_Info()<<METHOD<<"(): Calculate xs for '"
            <<m_name<<"' ("<<(p_gen?p_gen->Name():"")<<")"<<std::endl;
  double totalxs(psh->Integrate()/rpa->Picobarn());
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->Points()) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    exh->RemoveTerminatorObject(p_int);
    return 1;
  }
  exh->RemoveTerminatorObject(p_int);
  return 0;
}

void Single_Process::SetScale(const Scale_Setter_Arguments &args)
{
  if (IsMapped()) return;
  Scale_Setter_Arguments cargs(args);
  cargs.p_proc=this;
  cargs.p_cpls=&m_cpls;
  p_scale = Scale_Setter_Base::Scale_Getter_Function::
    GetObject(m_pinfo.m_scale=cargs.m_scale,cargs);
  if (p_scale==NULL) THROW(fatal_error,"Invalid scale scheme");
}

void Single_Process::SetKFactor(const KFactor_Setter_Arguments &args)
{
  if (IsMapped()) return;
  KFactor_Setter_Arguments cargs(args);
  cargs.p_proc=this;
  m_pinfo.m_kfactor=cargs.m_kfac;
  p_kfactor = KFactor_Setter_Base::KFactor_Getter_Function::
    GetObject(m_pinfo.m_kfactor=cargs.m_kfac,cargs);
  if (p_kfactor==NULL) THROW(fatal_error,"Invalid kfactor scheme");
}

void Single_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
}

bool Single_Process::Combinable
(const size_t &idi,const size_t &idj)
{
  return true;
}

const Flavour_Vector &Single_Process::
CombinedFlavour(const size_t &idij)
{
  static Flavour_Vector fls(1,kf_none);
  return fls;
}

ATOOLS::Flavour Single_Process::ReMap
(const ATOOLS::Flavour &fl,const size_t &id) const
{
  return fl;
}

Cluster_Amplitude *Single_Process::Cluster
(const Vec4D_Vector &p,const size_t &mode)
{
  MCatNLO_Process *mp(dynamic_cast<MCatNLO_Process*>(Parent()));
  if (mp) {
    Cluster_Amplitude *ampl(mp->GetAmplitude());
    if (ampl) return ampl;
  }
  if (!(mode&256)) {
    ClusterAmplitude_Vector &ampls(ScaleSetter(1)->Amplitudes());
    if (ampls.size()) {
      msg_Debugging()<<METHOD<<"(): Found "
		     <<ampls.size()<<" amplitude(s) ... ";
      msg_Debugging()<<"select 1st.\n";
      return ampls.front()->CopyAll();
    }
    if (mode&2048) return NULL;
  }
  PDF::Cluster_Definitions_Base* cd=p_shower->GetClusterDefinitions();
  int amode=0, cmode=mode;
  if (cd) {
    amode=cd->AMode();
    if (amode) cmode|=512;
    if (mode&512) cd->SetAMode(1);
  }
  p_gen->SetClusterDefinitions(cd);
  p_gen->PreCluster(this,p);
  Cluster_Amplitude* ampl(p_gen->ClusterConfiguration(this,p,cmode));
  if (ampl) ampl->Decays()=m_pinfo.m_fi.GetDecayInfos();
  if (cd) cd->SetAMode(amode);
  return ampl;
}
