#include "Tools.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

#include <algorithm>

using namespace ATOOLS;
using namespace PHASIC;

extern "C" {
  double li2(double x);
  double li3(double x);
  double li4(double x);
  double s2p2(double x);
  double h0p0m(double x); 
  double h00mp(double x);
  double h0mpp(double x);
  double h0mpm(double x);
  double h0mmp(double x);
  double hmpmp(double x);
}

PDF::PDF_Base *SHNNLO::s_pdf(NULL);
MODEL::Running_AlphaS *SHNNLO::s_as(NULL);
AMEGIC::Basic_Sfuncs *SHNNLO::s_bs;

long int SHNNLO::s_ntrials=-1;
double SHNNLO::s_disc, SHNNLO::s_p1, SHNNLO::s_p2;
double SHNNLO::s_z[4]={0.,0.,0.,0.};

double SHNNLO::pmap[6]={0,1,2,3,4,5};
double SHNNLO::s_pdfmin[2]={1.e-4,1.e-2};

Complex SHNNLO::spa(int i,int j) { return s_bs->S0(pmap[i],pmap[j]); }
Complex SHNNLO::spb(int i,int j) { return s_bs->S1(pmap[i],pmap[j]); }
double SHNNLO::sij(int i,int j) { return (spa(i,j)*spb(j,i)).real(); }

double SHNNLO::Li2(double z) {
  if(z>1.0)
    THROW(fatal_error, "Called Li2 with argument larger than one");
  return li2(z);
}
double SHNNLO::Li3(double z) {
  if(z>1.0)
    THROW(fatal_error, "Called Li3 with argument larger than one");
  return li3(z);
}
double SHNNLO::Li4(double z) {
  if(z>1.0)
    THROW(fatal_error, "Called Li4 with argument larger than one");
  return li4(z);
}
double SHNNLO::S2p2(double z) { 
  if(z>1.0)
    THROW(fatal_error, "Called S2p2 with argument larger than one");
  return s2p2(z);
}
double SHNNLO::H0p0m(double z) {
  if(z>1.0||z<-1.0)
    THROW(fatal_error, "Called H0p0m with argument larger than one or smaller than minus one");
  return h0p0m(z);
}
double SHNNLO::H00mp(double z) {
  if(z>1.0||z<-1.0)
    THROW(fatal_error, "Called H00mp with argument larger than one or smaller than minus one");
  return h00mp(z);
}
double SHNNLO::H0mpp(double z) {
  if(z>1.0||z<-1.0)
    THROW(fatal_error, "Called H0mpp with argument larger than one or smaller than minus one");
  return h0mpp(z);
}
double SHNNLO::H0mpm(double z) {
  if(z>1.0||z<-1.0)
    THROW(fatal_error, "Called H0mpm with argument larger than one or smaller than minus one");
  return h0mpm(z);
}
double SHNNLO::H0mmp(double z) {
  if(z>1.0||z<-1.0)
    THROW(fatal_error, "Called H0mmp with argument larger than one or smaller than minus one");
  return h0mmp(z);
}
double SHNNLO::Hmpmp(double z) {
  if(z>1.0||z<-1.0)
    THROW(fatal_error, "Called Hmpmp with argument larger than one or smaller than minus one");
  return hmpmp(z);
}

double SHNNLO::GetX(const ATOOLS::Vec4D &p,const int i)
{
  if (p[3]>0.0) return p.PPlus()/rpa->gen.PBeam(0).PPlus();
  return p.PMinus()/rpa->gen.PBeam(1).PMinus();
}

double SHNNLO::PDF(int i,double x,double muf)
{
  if (x>1.0) return 0.0;
  Flavour fl(i==0?Flavour(kf_gluon):Flavour(abs(i),i<0));
  s_pdf->Calculate(x,muf*muf);
  return s_pdf->GetXPDF(fl)/x;
}

double SHNNLO::GetPDF(const Flavour &fl,double x,double muf2)
{
  s_pdf->Calculate(x,muf2);
  return s_pdf->GetXPDF(fl)/x;
}

double SHNNLO::GetXPDF(Cluster_Leg *l,double muf2)
{
  double x(GetX(-l->Mom(),l->Id()));
  s_pdf->Calculate(x,muf2);
  double f(s_pdf->GetXPDF(l->Flav().Bar()));
  double q2min(sqr(2.0*l->Flav().Mass(true)));
  double min=s_pdfmin[0]*log(1.0-x)/log(1.0-s_pdfmin[1]);
  msg_Debugging()<<"f_{"<<l->Flav().Bar()<<"}("<<x<<","
		 <<sqrt(muf2)<<") = "<<f/x<<" <-> min = "
		 <<min<<", Q_{min} = "<<sqrt(q2min)<<"\n";
  if (dabs(f)<min || muf2<q2min) return 0.0;
  return f;
}

Cluster_Leg *SHNNLO::GetSplitter(const Cluster_Amplitude &ampl)
{
  for (size_t i(0);i<ampl.Legs().size();++i)
    if (ampl.Leg(i)->K()) return ampl.Leg(i);
  return NULL;
}

double SHNNLO::Beta0(const double &nf)
{
  return 11.0/6.0*3.0-2.0/3.0*0.5*nf;
}

double SHNNLO::Hab(const Flavour &a,const Flavour &b)
{
  if (a.IsQuark()) {
    if (b.IsQuark()) return a==b?4.0/3.0*3.0/2.0:0.0;
    return 0.0;
  }
  else {
    if (b.IsQuark()) return 0.0;
    return 11.0/6.0*3.0-2.0/3.0*0.5*(Flavour(kf_jet).Size()/2);
  }
}

double SHNNLO::FPab(const Flavour &a,const Flavour &b,const double &z)
{
  if (a.IsQuark()) {
    if (b.IsQuark()) return a==b?-4.0/3.0*(1.0+z):0.0;
    return 4.0/3.0*(1.0+sqr(1.0-z))/z;
  }
  else {
    if (b.IsQuark()) return 1.0/2.0*(z*z+sqr(1.0-z));
    return 3.0*2.0*((1.0-z)/z-1.0+z*(1.0-z));
  }
}

double SHNNLO::SPab(const Flavour &a,const Flavour &b,const double &z)
{
  if (a.IsQuark()) {
    if (b.IsQuark()) return a==b?4.0/3.0*2.0/(1.0-z):0.0;
    return 0.0;
  }
  else {
    if (b.IsQuark()) return 0.0;
    return 3.0*2.0/(1.0-z);
  }
}

double SHNNLO::IPab(const Flavour &a,const Flavour &b,const double &x)
{
  if (a.IsQuark()) {
    if (b.IsQuark() && a==b)
      return 4.0/3.0*2.0*log(1.0/(1.0-x));
    return 0.0;
  }
  else {
    if (b.IsQuark()) return 0.0;
    return 3.0*2.0*log(1.0/(1.0-x));
  }
}

double SHNNLO::CollinearCounterTerms
(Cluster_Leg *const l,
 const double &t1,const double &t2,
 const double &lmur2,const double &lmuf2,
 const double &ran)
{
  DEBUG_FUNC("Q = "<<sqrt(t1)<<" / "<<sqrt(t2));
  if (IsEqual(t1,t2)) return 0.0;
  msg_Debugging()<<"\\mu_F = "<<sqrt(lmuf2)<<"\n";
  msg_Debugging()<<"\\mu_R = "<<sqrt(lmur2)<<"\n";
  double as((*s_as)(lmur2)), x(GetX(-l->Mom(),l->Id()));
  double z(x+(1.0-x)*ran);
  double ct(0.0), lt(log(t1/t2));
  msg_Debugging()<<as<<"/(2\\pi) * log("<<sqrt(t1)<<"/"
		 <<sqrt(t2)<<") = "<<as/(2.0*M_PI)*lt<<"\n";
  Flavour jet(kf_jet), fl(l->Flav().Bar());
  double fb=GetPDF(fl,x,lmuf2);
  if (IsZero(fb)) {
    msg_Tracking()<<METHOD<<"(): Zero xPDF ( f_{"<<fl<<"}("
		  <<x<<","<<sqrt(lmuf2)<<") = "<<fb<<" ). Skip.\n";
    return 0.0;
  }
  msg_Debugging()<<"Beam "<<ID(l->Id())<<": z = "<<z<<", f_{"<<fl
		 <<"}("<<x<<","<<sqrt(lmuf2)<<") = "<<fb<<" {\n";
  for (size_t j(0);j<jet.Size();++j) {
    double Pf(FPab(jet[j],fl,z));
    double Ps(SPab(jet[j],fl,z));
    if (Pf+Ps==0.0) continue;
    double Pi(IPab(jet[j],fl,x));
    double H(Hab(jet[j],fl));
    double fa=GetPDF(jet[j],x/z,lmuf2);
    double fc=GetPDF(jet[j],x,lmuf2);
    msg_Debugging()<<"  P_{"<<jet[j]<<","<<fl
		   <<"}("<<z<<") = {F="<<Pf<<",S="<<Ps
		   <<",I="<<Pi<<"}, f_{"<<jet[j]<<"}("
		   <<x/z<<","<<sqrt(lmuf2)<<") = "<<fa
		   <<", f_{"<<jet[j]<<"}("<<x<<","
		   <<sqrt(lmuf2)<<") = "<<fc<<"\n";
    if (IsZero(fa)||IsZero(fc)) {
      msg_Tracking()<<METHOD<<"(): Zero xPDF. Skip.\n";
      return 0.0;
    }
    ct+=as/(2.0*M_PI)*lt*
      ((fa/z*Pf+(fa/z-fc)*Ps)*(1.0-x)+fc*(H-Pi))/fb;
  }
  msg_Debugging()<<"} -> "<<ct<<"\n";
  return ct;
}

double SHNNLO::Weight
(double &ct,Cluster_Amplitude *const ampl,
 const double &mur2,const double &muf2,
 const double *k0sq,const int mode)
{
  ct=0.0;
  if (ampl->Next()==NULL) return 1.0;
  if (ampl->OrderQCD()-1!=ampl->Next()->OrderQCD())
    THROW(fatal_error,"Invalid shower history");
  DEBUG_FUNC(ampl);
  msg_Debugging()<<*ampl->Next()<<"\n";
  double kt2c(k0sq[0]), kt2(ampl->KT2());
  Cluster_Leg *split(GetSplitter(*ampl->Next()));
  if (split->Id()&3) kt2c=k0sq[1];
  if (kt2<kt2c) {
    if ((mode&2) && ampl->Prev()->KT2()>kt2) return 1.0;
    return sqrt(-1.0);
  }
  if (kt2>ampl->Next()->KT2()) return 1.0;
  if (mode!=1) return 1.0;
  double Q2[2]={muf2,muf2};
  double asmu((*s_as)(mur2));
  double askt((*s_as)(ampl->Mu2()));
  if (mode&1) {
    double cpl(asmu);
    std::vector<double> ths(MODEL::as->Thresholds(mur2,ampl->Mu2()));
    if (ampl->Mu2()<mur2) std::reverse(ths.begin(),ths.end());
    if (ths.empty() || !IsEqual(ampl->Mu2(),ths.back())) ths.push_back(ampl->Mu2());
    if (!IsEqual(mur2,ths.front())) ths.insert(ths.begin(),mur2);
    for (size_t i(1);i<ths.size();++i) {
      double nf=MODEL::as->Nf((ths[i]+ths[i-1])/2.0);
      double L=log(ths[i]/ths[i-1]), lct=cpl/(2.0*M_PI)*Beta0(nf)*L;
      cpl*=1.0-lct;
      ct+=lct;
    }
  }
  double K=askt/asmu;
  msg_Debugging()<<"\\alpha_s weight: "<<askt<<" / "
		 <<asmu<<" = "<<askt/asmu<<", ct: "<<ct
		 <<" -> rel. diff "<<asmu*(1.0-ct)/askt-1.0<<"\n";
  for (int i(0);i<2;++i) {
    Cluster_Leg *li(ampl->Leg(i));
    if ((li->Id()!=split->K() || split->Id()&3) &&
	(li->Id()&split->Id())==0) continue;
    double fkt=GetXPDF(li,kt2);
    double fmu=GetXPDF(li,Q2[i]);
    if (IsZero(fmu)) { ct=0.0; return 1.0; }
    double cct=(mode&1)?CollinearCounterTerms(li,kt2,Q2[i],mur2,muf2,s_z[i]):0.0;
    ct-=cct;
    K*=fkt/fmu;
    Q2[i]=kt2;
    msg_Debugging()<<"PDF term("<<i<<"): "<<fkt<<" / "
		   <<fmu<<" = "<<fkt/fmu<<", ct: "<<-cct<<"\n";
  }
  msg_Debugging()<<"K = "<<K<<", CT = "<<ct
		 <<" -> "<<(K+ct)<<"\n";
  if (ampl->Next()->Next()==NULL) {
    double muc2(ampl->Next()->KT2());
    for (int i(0);i<2;++i) {
      Cluster_Leg *li(ampl->Next()->Leg(i));
      if (IsEqual(muc2,Q2[i])) continue;
      double fkt=GetXPDF(li,muc2);
      double fmu=GetXPDF(li,Q2[i]);
      if (IsZero(fmu)) { ct=0.0; return 1.0; }
      double cct=(mode&1)?CollinearCounterTerms(li,muc2,Q2[i],mur2,muf2,s_z[i+2]):0.0;
      ct-=cct;
      K*=fkt/fmu;
      msg_Debugging()<<"PDF term("<<i<<"): "<<fkt<<" / "
		     <<fmu<<" = "<<fkt/fmu<<", ct: "<<-cct<<"\n";
    }
  }
  msg_Debugging()<<"K = "<<K<<", CT = "<<ct
		 <<" -> "<<(K+ct)<<"\n";
  return K;
}

double SHNNLO::NLODiffWeight
(Process_Base *const proc,double &wgt,
 const double &mur2,const double &muf2,
 const double *k0sq,const int fomode,
 const int umode,const std::string &varid)
{
  if (fomode) return wgt;
  DEBUG_FUNC(proc->Name());
  Cluster_Amplitude *ampl(NULL);
  Scale_Setter_Base *sc(proc->ScaleSetter());
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  if (ampl==NULL || ampl->Next()==NULL) return wgt;
  msg_Debugging()<<*ampl<<"\n";
  double w1, K=Weight(w1,ampl,mur2,muf2,k0sq,0);
  msg_Debugging()<<"K = "<<K<<"\n";
  if (IsBad(K)) {
    ampl->Next()->SetNLO(128);
    return umode?0.0:wgt;
  }
  bool gen=rpa->gen.NumberOfTrials()>s_ntrials;
  s_ntrials=rpa->gen.NumberOfTrials();
  if (gen) s_disc=ran->Get();
  else msg_Debugging()<<"keep random point\n";
  double p1=s_p1;
  if (gen) p1=s_p1=1.0/(1.0+2.0*dabs(K-1.0));
  if (umode) p1=0.0;
  if (s_disc<=p1) {
    wgt*=1.0/p1;
    ampl->Next()->SetNLO(16|64);
  }
  else {
    double pw(0.5*(1.0-p1));
    wgt*=(umode?K:K-1.0)/pw;
    ampl->Next()->SetNLO(16);
    if (s_disc>p1+pw) {
      ampl->Next()->SetNLO(16|32);
      wgt*=-1.0;
    }
  }
  msg_Debugging()<<"K = "<<wgt<<"\n";
  return wgt;
}

double SHNNLO::NNLODiffWeight
(Process_Base *const proc,double &wgt,
 const double &mur2,const double &muf2,
 const double *k0sq,const int mode,
 const int fomode,const int umode,
 const std::string &varid)
{
  nlo_type::code nlot=proc->Info().m_fi.m_nlotype;
  DEBUG_FUNC(proc->Name()<<", wgt = "<<wgt
	     <<", type "<<nlot<<", mode = "<<mode);
  if (fomode || rpa->gen.NumberOfTrials()<1) return wgt;
  Cluster_Amplitude *ampl(NULL);
  Scale_Setter_Base *sc(proc->ScaleSetter(1));
  if (sc->Amplitudes().size()) ampl=sc->Amplitudes().front();
  NLO_subevtlist *subs(proc->GetRSSubevtList());
  if (subs) ampl=subs->back()->Proc<Single_Process>()->
	      ScaleSetter(1)->Amplitudes().front()->Next();
  if (ampl==NULL || ampl->Next()==NULL)
    return wgt=SetWeight(NULL,mode,wgt,1.0,0.0,umode,varid);
  msg_Debugging()<<*ampl<<"\n";
  double w1, w=Weight(w1,ampl,mur2,muf2,k0sq,mode|(subs?2:0));
  msg_Debugging()<<"w = "<<w<<", w1 = "<<w1<<"\n";
  if (IsBad(w) || IsBad(w1)) {
    ampl->Next()->SetNLO(128);
    if (umode) return 0.0;
    return wgt=SetWeight(NULL,mode,wgt,1.0,0.0,umode,varid);
  }
  return wgt=SetWeight(ampl->Next(),mode,wgt,w,w1,umode,varid);
}

double SHNNLO::NNLODeltaWeight
(Process_Base *const proc,const double &wgt,const int fomode)
{
  nlo_type::code nlot=proc->Info().m_fi.m_nlotype;
  DEBUG_FUNC(proc->Name()<<", 0j type "<<nlot);
  if (fomode) return wgt;
  NLO_subevtlist *subs(proc->GetSubevtList());
  if (subs==NULL) subs=proc->GetRSSubevtList();
  if (subs==NULL) return wgt;
  if (subs->back()->p_ampl && subs->back()->p_ampl->Next())
    subs->back()->p_ampl->Next()->SetNLO(128);
  return wgt;
}

double SHNNLO::SetWeight
(ATOOLS::Cluster_Amplitude *const ampl,const int mode,
 double wgt,const double &w,const double &w1,
 const int umode,const std::string &varid)
{
  bool gen=rpa->gen.NumberOfTrials()>s_ntrials;
  s_ntrials=rpa->gen.NumberOfTrials();
  if (gen) s_disc=ran->Get();
  else msg_Debugging()<<"keep random point\n";
  double p1=s_p1;
  if (gen) {
    s_p2=0.0;
    p1=s_p1=1.0/(1.0+2.0*dabs(w-1.0)+2.0*dabs(w*w1));
  }
  if (umode) p1=0.0;
  if (s_disc<=p1) {
    wgt*=1.0/p1;
    if (ampl) {
      if (mode==1) ampl->SetFlag(2);
      ampl->SetNLO(16|64);
    }
  }
  else {
    double p2=s_p2;
    if (gen) {
      p2=s_p2=1.0-p1;
      if (w!=1.0) p2=s_p2/=(1.0+dabs(w*w1/(w-1.0)));
    }
    if (umode) p2=1.0;
    if (s_disc<=p1+p2) {
      double pw(0.5*p2);
      wgt*=(umode?w:w-1.0)/pw;
      if (ampl) {
	if (mode==1) ampl->SetFlag(2);
	ampl->SetNLO(16);
      }
      if (s_disc>p1+pw) {
	if (ampl) {
	  if (mode==1) ampl->SetFlag(2);
	  ampl->SetNLO(16|32);
	}
	wgt*=-1.0;
      }
    }
    else {
      double pw(0.5*(1.0-p1-p2));
      wgt*=w*w1/pw;
      if (ampl)	ampl->SetNLO(16);
      if (s_disc>p1+p2+pw) {
	if (ampl) ampl->SetNLO(16|32);
	wgt*=-1.0;
      }
    }
  }
  msg_Debugging()<<"wgt = "<<wgt<<" (w = "<<w<<", w1 = "<<w1
		 <<", p1 = "<<s_p1<<", p2 = "<<s_p2
		 <<", # = "<<s_disc<<") "<<varid<<"\n";
  return wgt;
}
