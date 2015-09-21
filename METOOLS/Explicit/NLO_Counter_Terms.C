#include "METOOLS/Explicit/NLO_Counter_Terms.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace METOOLS;
using namespace ATOOLS;

double METOOLS::Hab(const Flavour &a,const Flavour &b)
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

double METOOLS::FPab(const Flavour &a,const Flavour &b,const double &z)
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

double METOOLS::SPab(const Flavour &a,const Flavour &b,const double &z)
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

double METOOLS::IPab(const Flavour &a,const Flavour &b,const double &x)
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

double METOOLS::AlphaSCounterTerm(const double& muR2,const double& muR2ref,
                                  const double& as,
                                  MODEL::One_Running_AlphaS * oras,
                                  const size_t oqcd)
{
  DEBUG_FUNC(muR2<<" -> "<<muR2ref);
  if (oqcd==0) return 0.;
  if (IsEqual(muR2,muR2ref)) return 0.;
  // if flavour threshold between muR2 and muR2ref split beta0*log(muR2/muR2ref)
  // into regions with constant nf
  std::vector<double> thrs(oras->Thresholds(muR2,muR2ref));
  msg_Debugging()<<"Flavour thresholds in range ["<<muR2<<","<<muR2ref<<"]: "
                 <<thrs<<std::endl;
  thrs.push_back((muR2>muR2ref)?muR2:muR2ref);
  thrs.insert(thrs.begin(),(muR2>muR2ref)?muR2ref:muR2);
  double betalog(0.);
  msg_Debugging()<<"\\sum_{\\mu_{th}} \\beta_0(nf(\\mu_i)) "
                 <<"log(\\mu_{i+1}/\\mu_i) = "<<std::endl;
  for (size_t i(0);i<thrs.size()-1;++i) {
    msg_Debugging()<<(i==0?"    ":"  + ")<<oras->Beta0(thrs[i+1])
                   <<" * "<<log(thrs[i+1]/thrs[i])<<"  (nf="
                   <<oras->Nf(thrs[i+1])<<", "<<thrs[i]<<".."<<thrs[i+1]
                   <<") \n";
    betalog+=oras->Beta0(thrs[i+1])*log(thrs[i+1]/thrs[i]);
  }
  if (muR2>muR2ref) betalog*=-1.;
  msg_Debugging()<<"  = "<<betalog<<std::endl;
  msg_Debugging()<<"\\alpha_s term: "<<oqcd<<" * "<<as
                 <<"/2\\pi * \\sum_{\\mu_{th}} \\beta_0(n_f(\\mu_i)) "
                 <<"log(\\mu_{i+1}/\\mu_i) = "
                 <<oqcd*as/M_PI*betalog<<std::endl;
  return double(oqcd)*as/M_PI*betalog;
}

double METOOLS::CollinearCounterTerms
(const Flavour &fl,const double &x,const double &z,const double &as,
 const double &t1,const double &t2,const double &tref,PDF::PDF_Base * pdf)
{
  if (!pdf) return 0.;
  static double th(1.0e-12);
  DEBUG_FUNC("Q = "<<sqrt(t1)<<" / "<<sqrt(t2));
  if (IsEqual(t1,t2)) return 0.0;
  msg_Debugging()<<"\\mu_F = "<<sqrt(tref)<<"\n";
  double ct(0.0), lt(log(t1/t2));
  msg_Debugging()<<as<<"/(2\\pi) * log("<<sqrt(t1)<<"/"
                 <<sqrt(t2)<<") = "<<as/(2.0*M_PI)*lt<<"\n";
  Flavour jet(kf_jet);
  if (!pdf->Contains(fl)) return 0.;
  pdf->Calculate(x,tref);
  double fb(pdf->GetXPDF(fl)/x);
  if (IsZero(fb,th)) {
    msg_Tracking()<<METHOD<<"(): Zero xPDF ( f_{"<<fl<<"}("
                  <<x<<","<<sqrt(tref)<<") = "<<fb<<" ). Skip.\n";
    return 0.0;
  }
  msg_Debugging()<<"z = "<<z<<", f_{"<<fl
                 <<"}("<<x<<","<<sqrt(tref)<<") = "<<fb<<" {\n";
  for (size_t j(0);j<jet.Size();++j) {
    if (!pdf->Contains(jet[j])) continue;
    double Pf(METOOLS::FPab(jet[j],fl,z));
    double Ps(METOOLS::SPab(jet[j],fl,z));
    if (Pf+Ps==0.0) continue;
    double Pi(METOOLS::IPab(jet[j],fl,x));
    double H(METOOLS::Hab(jet[j],fl));
    pdf->Calculate(x/z,tref);
    double fa(pdf->GetXPDF(jet[j])/x);
    pdf->Calculate(x,tref);
    double fc(pdf->GetXPDF(jet[j])/x);
    msg_Debugging()<<"  P_{"<<jet[j]<<","<<fl
                   <<"}("<<z<<") = {F="<<Pf<<",S="<<Ps
                   <<",I="<<Pi<<"}, f_{"<<jet[j]<<"}("
                   <<x/z<<","<<sqrt(tref)<<") = "<<fa
                   <<", f_{"<<jet[j]<<"}("<<x<<","
                   <<sqrt(tref)<<") = "<<fc<<"\n";
    ct+=as/(2.*M_PI)*lt*(fa/z*(Pf+(1.-x)*Ps)+fc*(H-Pi-(1.-x)*Ps))/fb;
  }
  msg_Debugging()<<"} -> "<<ct<<"\n";
  return ct;
}
