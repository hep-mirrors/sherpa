#include "PHASIC++/Process/Single_Process.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
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
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "METOOLS/Explicit/NLO_Counter_Terms.H"
#include "MODEL/Main/Coupling_Data.H"
#include "MODEL/Main/Running_AlphaS.H"

#include <algorithm>

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Single_Process::Single_Process():
  m_lastxs(0.0), m_lastbxs(0.0), m_dsweight(1.0), m_lastflux(0.0),
  m_zero(false), m_dads(true), m_pdfcts(true), m_nfconvscheme(0)
{
  Settings& s = Settings::GetMainSettings();
  m_pdfcts = s["MEPSNLO_PDFCT"].SetDefault(true).Get<bool>();
  m_dads = s["MCNLO_DADS"].SetDefault(true).Get<bool>();

  std::string ncs{ s["NLO_NF_CONVERSION_TERMS"]
    .SetDefault("None")
    .UseNoneReplacements()
    .Get<std::string>() };
  if      (ncs=="None" || ncs=="0") m_nfconvscheme=0;
  else if (ncs=="On"   || ncs=="1") m_nfconvscheme=1;
  else if (ncs=="Only" || ncs=="2") m_nfconvscheme=2;
  else THROW(fatal_error,"Unknown NLO_NF_CONVERSION_TERMS scheme.");
  static bool printed(false);
  if (!printed && m_nfconvscheme>0) {
    if (m_nfconvscheme==1)
      msg_Info()<<"NLO n_f scheme conversion terms added to I-operator."<<std::endl;
    if (m_nfconvscheme==2)
      msg_Info()<<"NLO n_f scheme conversion terms computed only."<<std::endl;
    printed=true;
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

double Single_Process::KFactor(const int mode) const
{
  if (p_kfactor && p_kfactor->On()) return p_kfactor->KFactor(mode);
  return 1.0;
}

double Single_Process::NfSchemeConversionTerms() const
{
  if (m_nfconvscheme==0) return 0.;
  DEBUG_FUNC("scheme="<<m_nfconvscheme<<", nlo-type="<<m_pinfo.m_nlomode);
  if (m_nfconvscheme==1 && !(m_pinfo.Has(nlo_type::vsub))) return 0.;
  // determine nf of PDF, if no PDF take from alphaS
  // nfl from number of light quarks
  int nfa(p_int->ISR()->PDF(0)?p_int->ISR()->PDF(0)->ASInfo().m_nf:-1);
  int nfb(p_int->ISR()->PDF(1)?p_int->ISR()->PDF(1)->ASInfo().m_nf:-1);
  if (nfa!=nfb)       THROW(fatal_error,"Nf(a) and Nf(b) differ. Abort.");
  if (nfb==-1 && nfb!=MODEL::as->Nf(1.e20))
    THROW(fatal_error,"Inconsistent Nf between PDFs and \\alpha_s.");
  if (nfb==-1) nfa=MODEL::as->Nf(1.e20);
  const size_t nf(nfb);
  const size_t nfl(Flavour(kf_quark).Size()/2);
  if (nfl==nf) return 0.;
  // debugging output
  msg_Debugging()<<"a="<<m_flavs[0]<<", b="<<m_flavs[1]
                 <<", nf(a)="<<nfa<<", nf(b)="<<nfb<<", nfl="<<nfl<<std::endl;
  // compute the scheme conversion terms
  MODEL::Coupling_Data *qcdcpl(m_cpls.Get("Alpha_QCD"));
  double as(qcdcpl->Default()*qcdcpl->Factor());
  MODEL::Coupling_Data *qedcpl(m_cpls.Get("Alpha_QED"));
  double aqed(qedcpl->Default()*qedcpl->Factor());
  double res(0.);
  double p(m_maxcpl[0]-m_pinfo.m_fi.m_nlocpl[0]);
  double facqcd(as/3./M_PI*MODEL::as->TR()), facqed(aqed/3/M_PI);
  msg_Debugging()<<"p="<<p<<", facqcd="<<facqcd<<", facqed="<<facqed<<std::endl;
  double muR2(ScaleSetter(1)->Scale(stp::ren));
  double muF2(ScaleSetter(1)->Scale(stp::fac));
  for (size_t i(nfl+1);i<=nf;++i) {
    Flavour fl((kf_code)(i));
    double m2(sqr(fl.Mass()));
    msg_Debugging()<<fl<<": m="<<sqrt(m2)<<", \\mu_R="<<sqrt(muR2)
                   <<", \\mu_F="<<sqrt(muF2)<<std::endl;
    // logs only if scale larger than mass
    double logm2muR2(m2<muR2?log(m2/muR2):0.);
    double logm2muF2(m2<muF2?((muR2==muF2)?logm2muR2:log(m2/muF2)):0.);
    msg_Debugging()<<"log(m2/muR2)="<<logm2muR2
                   <<", log(m2/muF2)="<<logm2muF2<<std::endl;
    // alphaS term
    if (p && logm2muR2 && m_flavs[0].Strong()) res+=facqcd*p*logm2muR2;
    // PDF terms
    for (size_t j(0);j<2;++j) {
      if (logm2muF2 &&
          p_int->ISR()->PDF(j) && p_int->ISR()->PDF(j)->Contains(m_flavs[j])) {
        if      (m_flavs[j].IsGluon())  res-=facqcd*logm2muF2;
        else if (m_flavs[j].IsPhoton()) res-=facqed*sqr(fl.Charge())*logm2muF2;
      }
    }
  }
  // scheme dependent return value
  msg_Debugging()<<"B = "<<m_lastbxs<<" ,  C_nf = "<<res<<std::endl
                 <<" => "<<(m_nfconvscheme==2?0.:m_last)+m_lastbxs*res<<std::endl;
  return (m_nfconvscheme==2?-m_last:0.)+m_lastbxs*res;
}

double Single_Process::CollinearCounterTerms
(const int i,const Flavour &fl,const Vec4D &p,
 const double &z,const double &t1,const double &t2,
 const double &muf2fac, const double &mur2fac,
 MODEL::Running_AlphaS * as) const
{
  if (!(p_int->ISR() && p_int->ISR()->On()&(1<<i))) return 0.0;
  static double th(1.0e-12);
  DEBUG_FUNC("Q = "<<sqrt(t1)<<" / "<<sqrt(t2));
  if (IsEqual(t1,t2)) return 0.0;

  // determine scales
  double lmuf2(p_scale->Scale(stp::fac)*muf2fac);
  double lmur2(p_scale->Scale(stp::ren)*mur2fac);
  msg_Debugging()<<"\\mu_F = "<<sqrt(lmuf2)<<"\n";
  msg_Debugging()<<"\\mu_R = "<<sqrt(lmur2)<<"\n";

  // determine running AlphaS and evaluate at lmur2
  // if as is not given, the nominal results will be used
  double asvalue(0.0);
  if (as) {
    asvalue = (*as)(lmur2);
 } else {
    MODEL::Coupling_Data *cpl(m_cpls.Get("Alpha_QCD"));
    asvalue = cpl->Default() * cpl->Factor();
  }

  // determine counter-term
  double ct(0.0), lt(log(t1/t2)), x(p_int->ISR()->CalcX(p));
  msg_Debugging()<<asvalue<<"/(2\\pi) * log("<<sqrt(t1)<<"/"
		 <<sqrt(t2)<<") = "<<asvalue/(2.0*M_PI)*lt<<"\n";
  Flavour jet(kf_jet);
  double fb=p_int->ISR()->PDFWeight(1<<(i+1),p,p,lmuf2,lmuf2,fl,fl,0);
  if (IsZero(fb,th)) {
    msg_Tracking()<<METHOD<<"(): Zero xPDF ( f_{"<<fl<<"}("
		  <<x<<","<<sqrt(lmuf2)<<") = "<<fb<<" ). Skip.\n";
    return 0.0;
  }

  // skip PDF ratio if high-x sanity condition not fullfilled
  if (dabs(fb)<1.0e-4*log(1.0 - x)/log(1.0 - 1.0e-2)){
    msg_Debugging() << "Invalid pdf ratio, ct set to zero." << std::endl;
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
      msg_Tracking()<<METHOD<<"(): Zero xPDF. No contrib from "<<j
                    <<". Skip .\n";
    }
    ct+=asvalue/(2.0*M_PI)*lt*
      ((fa/z*Pf+(fa/z-fc)*Ps)*(1.0-x)+fc*(H-Pi))/fb;
  }
  msg_Debugging()<<"} -> "<<ct<<"\n";
  return ct;
}

ATOOLS::Cluster_Sequence_Info Single_Process::ClusterSequenceInfo(
    const ATOOLS::ClusterAmplitude_Vector &ampls,
    const double &Q2,
    const double &muf2fac,
    const double &mur2fac,
    const double &showermuf2fac,
    MODEL::Running_AlphaS * as,
    const ATOOLS::Cluster_Sequence_Info * const nominalcsi)
{
  if (!m_use_biweight) {
    return 1.;
  }
  if (m_nin == 1) {
    return 1.0;
  } else if (m_nin > 2) {
    THROW(not_implemented, "More than two incoming particles.");
  }
  Cluster_Sequence_Info csi;
  AddISR(csi, ampls, Q2, muf2fac, mur2fac, showermuf2fac, as, nominalcsi);
  AddBeam(csi, Q2);
  return csi;
}

void Single_Process::AddISR(ATOOLS::Cluster_Sequence_Info &csi,
            const ATOOLS::ClusterAmplitude_Vector &ampls,
            const double &Q2,
            const double &muf2fac, const double &mur2fac,
            const double &showermuf2fac,
            MODEL::Running_AlphaS * as,
            const ATOOLS::Cluster_Sequence_Info * const nominalcsi)
{
  DEBUG_FUNC(Name());
  if (p_int->ISR()) {
    // add external PDF weight (before clustering)
    double pdfext(p_int->ISR()->PDFWeight(0,
                                          p_int->Momenta()[0],
                                          p_int->Momenta()[1],
                                          Q2, Q2,
                                          m_flavs[0], m_flavs[1]));
    msg_Debugging()<<"PDF(fla="<<m_flavs[0]
                   <<", xa="<<p_int->ISR()->CalcX(p_int->Momenta()[0])
                   <<", Qa="<<sqrt(Q2)<<") * "
                   <<"PDF(flb="<<m_flavs[1]
                   <<", xb="<<p_int->ISR()->CalcX(p_int->Momenta()[1])
                   <<", Qb="<<sqrt(Q2)<<") -> "<<pdfext<<std::endl;
    csi.AddWeight(pdfext);

    // add splittings and their PDF weight ratios from clustering
    if (ampls.size() && (m_pinfo.m_ckkw&1)) {
      DEBUG_FUNC(m_name<<", \\mu_F = "<<sqrt(Q2));
      m_mewgtinfo.m_type|=mewgttype::METS;

      // add external splitting
      csi.AddSplitting(Q2,
                       p_int->ISR()->CalcX(p_int->Momenta()[0]),
                       p_int->ISR()->CalcX(p_int->Momenta()[1]),
                       m_flavs[0], m_flavs[1]);
      csi.AddPDFRatio(pdfext, 1.);

      Cluster_Amplitude *ampl(ampls.front());
      msg_IODebugging()<<*ampl<<"\n";

      // add subsequent splittings
      bool addedfirstsplitting(false);
      double currentQ2(Q2);
      double currentscalefactor(1.0);
      double pdfnum(pdfext), pdfden(pdfext);
      for (; ampl; ampl = ampl->Next()) {
        // skip decays, equal scales, unordered configs,
        // and quarks below threshold

        // skip decays (they are not even added to the splittings)
        msg_IODebugging()<<*ampl<<"\n";
	if (ampl->NLO()&1) {
          msg_Debugging()<<"Skip. NLO clustering "<<
            ID(ampl->Next()->Splitter()->Id())<<"\n";
          continue;
	}
	if (ampl->Next() && ampl->Next()->Splitter()->Stat() == 3) {
          msg_Debugging()<<"Skip. Decay "<<
            ID(ampl->Next()->Splitter()->Id())<<"\n";
          continue;
	}

        // add splitting
	Flavour f1(ampl->Leg(0)->Flav().Bar());
	Flavour f2(ampl->Leg(1)->Flav().Bar());
	if (MapProc() && LookUp() && !m_pinfo.Has(nlo_type::real)) {
	  f1=ReMap(f1, ampl->Leg(0)->Id());
	  f2=ReMap(f2, ampl->Leg(1)->Id());
	}
	// ampl->KT2() has actual splitting scale
	// ampl->MuF() is fac-scale as given in scale setter
	csi.AddSplitting(ampl->KT2(),
			 p_int->ISR()->CalcX(-ampl->Leg(0)->Mom()),
			 p_int->ISR()->CalcX(-ampl->Leg(1)->Mom()),
			 f1, f2);

        // skip equal scales
        if (IsEqual(currentQ2 / currentscalefactor, ampl->Next() ? showermuf2fac * ampl->KT2() : Q2)) {
            msg_Debugging()<<"Skip. Scales equal: t_i="<<currentQ2 / currentscalefactor
                           <<", t_{i+1}="<<(ampl->Next()?ampl->KT2():Q2)
                           <<std::endl;
          if (ampl->Next() == NULL) {
            csi.AddPDFRatio(1., pdfden);
          } else {
            csi.AddPDFRatio(pdfnum, pdfden);
          }
          continue;
        }

        // skip unordered configuration
	if (addedfirstsplitting && (currentQ2 / currentscalefactor > ampl->KT2())) {
	  msg_Debugging()<<"Skip. Unordered history "<<
	    sqrt(currentQ2 / currentscalefactor)<<" > "<<sqrt(ampl->KT2())<<"\n";
	  currentQ2 = sqrt(std::numeric_limits<double>::max());
	  csi.AddPDFRatio(1.,1.);
	  continue;
	}

        // skip when a scale is below a (quark) mass threshold
        if (currentQ2 / currentscalefactor < sqr(2.0 * f1.Mass(true))
            || currentQ2 / currentscalefactor < sqr(2.0 * f2.Mass(true))) {
          msg_Debugging()<<"Skip. Quarks below threshold: t="
                         <<currentQ2 / currentscalefactor
                         <<" vs. "<<sqr(2.0*f1.Mass(true))
                         <<" / "<<sqr(2.0*f2.Mass(true))<<std::endl;
          continue;
        }

	// denominators
	double wd1(p_int->ISR()->PDFWeight(2,
                                           -ampl->Leg(0)->Mom(),
                                           -ampl->Leg(1)->Mom(),
                                           currentQ2, currentQ2, f1, f2,
                                           0));
	double wd2(p_int->ISR()->PDFWeight(4,
                                           -ampl->Leg(0)->Mom(),
                                           -ampl->Leg(1)->Mom(),
                                           currentQ2, currentQ2, f1, f2,
                                           0));

        // new scale (note: for the core scale we use Q2 instead of ampl->MuF2
        // because we might be reweighting and Q2 could have been multiplied
        // by a scaling factor, whereas ampl->MuF2 would not reflect this)
        double lastQ2 = currentQ2;
        if (ampl->Next() == NULL) {
          currentQ2 = Q2;
          currentscalefactor = 1.0;
        } else {
          currentQ2 = showermuf2fac * ampl->KT2();
          currentscalefactor = showermuf2fac;
        }

        // skip when a scale is below a (quark) mass threshold, new scale
        if (currentQ2 < sqr(2.0 * f1.Mass(true)) || currentQ2 < sqr(2.0 * f2.Mass(true))) {
          msg_Debugging()<<"Skip. Quarks below threshold: t="<<currentQ2
                         <<" vs. "<<sqr(2.0*f1.Mass(true))
                         <<" / "<<sqr(2.0*f2.Mass(true))<<std::endl;
          continue;
        }

	// numerators
	double wn1(p_int->ISR()->PDFWeight(2,
                                           -ampl->Leg(0)->Mom(),
                                           -ampl->Leg(1)->Mom(),
                                           currentQ2, currentQ2, f1, f2,
                                           0));
	double wn2(p_int->ISR()->PDFWeight(4,
                                           -ampl->Leg(0)->Mom(),
                                           -ampl->Leg(1)->Mom(),
                                           currentQ2, currentQ2, f1, f2,
                                           0));

        double x1=p_int->ISR()->CalcX(-ampl->Leg(0)->Mom());
        double x2=p_int->ISR()->CalcX(-ampl->Leg(1)->Mom());

        // skip PDF ratio if high-x sanity condition not fullfilled
        auto validratio1 = (!IsZero(wn1) && !IsZero(wd1));
        if (validratio1 && x1 < 1.0)
          validratio1 = !(dabs(wd1)<1.0e-4*log(1.0 - x1)/log(1.0 - 1.0e-2));
        if (validratio1) {
          csi.AddWeight(wn1 / wd1);
        } else {
          msg_Debugging() << "invalid pdf ratio in beam 0," << std::endl;
          msg_Debugging() << "skip weight." << std::endl;
        }
        auto validratio2 = (!IsZero(wn2) && !IsZero(wd2));
        if (validratio2 && x2 < 1.0)
          validratio2 = !(dabs(wd2)<1.0e-4*log(1.0 - x2)/log(1.0 - 1.0e-2));
        if (validratio2) {
          csi.AddWeight(wn2 / wd2);
        } else {
          msg_Debugging() << "invalid pdf ratio in beam 1," << std::endl;
          msg_Debugging() << "skip weight." << std::endl;
        }

	// book-keep PDF ratios excl.
	//   a) first one correcting outer PDF from muF to t
	//   b) last numerator taken at muF (this one is to be varied)
	// use the following identity with i=0 -> core and i=N -> ext
	// wn-ext * [\prod_{i=0}^{N-1} wn_i/wd_i]
	// = [wn-ext * \prod_{i=1}^{N-1} wn_i/wd_i * 1/wd_0] * wn-core
	// = [\prod_{i=1}^N wn_i/wd_{i-1}] * wn-core
	pdfnum = wn1 * wn2;
	pdfden = wd1 * wd2;
	if (ampl->Next() == NULL) {
          csi.AddPDFRatio(1., pdfden);
        } else {
          csi.AddPDFRatio(pdfnum, pdfden);
        }
	msg_Debugging()<<"* [  "
		       <<"PDF(fla="<<f1
		       <<", xa="<<p_int->ISR()->CalcX(-ampl->Leg(0)->Mom())
		       <<", ta="<<currentQ2<<") * "
		       <<"PDF(flb="<<f2
		       <<", xb="<<p_int->ISR()->CalcX(-ampl->Leg(1)->Mom())
		       <<", tb="<<currentQ2<<") -> "<<wn1*wn2<<"\n"
		       <<"   / "
		       <<"PDF(fla="<<f1
		       <<", xa="<<p_int->ISR()->CalcX(-ampl->Leg(0)->Mom())
		       <<", ta="<<lastQ2<<") * "
		       <<"PDF(flb="<<f2
		       <<", xb="<<p_int->ISR()->CalcX(-ampl->Leg(1)->Mom())
		       <<", tb="<<lastQ2<<") -> "<<wd1*wd2
		       <<" ] = "<<wn1*wn2/wd1/wd2<<std::endl;

        // add collinear counterterm
        if (m_pdfcts && m_pinfo.Has(nlo_type::born)) {
          for (int i(0); i < 2; ++i) {
            // skip PDF ratio if high-x sanity condition not fullfilled
            if (i == 0 && (IsZero(wn1) || IsZero(wd1) || (dabs(wd1)<1.0e-4*log(1.0 - x1)/log(1.0 - 1.0e-2)) )) continue;
            if (i == 1 && (IsZero(wn2) || IsZero(wd2) || (dabs(wd2)<1.0e-4*log(1.0 - x2)/log(1.0 - 1.0e-2)) )) continue;
            Vec4D p(-ampl->Leg(i)->Mom());
            const double x(p_int->ISR()->CalcX(p));
            double z(-1.0);
            if (nominalcsi) {
              const size_t currentsplittingindex = csi.m_txfl.size() - 1;
              // this might return -1.0 if there were no counterterms
              // calculated for this splitting with nominal parameters
              z = (i == 0) ?
                nominalcsi->m_txfl[currentsplittingindex].m_xap :
                nominalcsi->m_txfl[currentsplittingindex].m_xbp;
            }
            if (z == -1.0) {
              z = x + (1.0 - x) * ran->Get();
            }
            const double ct(CollinearCounterTerms(i, i ? f2 : f1, p, z,
                                                  currentQ2, lastQ2,
                                                  muf2fac, mur2fac, as));
            csi.AddCounterTerm(ct, z, i);
          }
        }
        addedfirstsplitting = true;
      }
    }
  }
}

void Single_Process::AddBeam(ATOOLS::Cluster_Sequence_Info& csi,
                             const double& Q2)
{
  if (p_int->Beam() && p_int->Beam()->On()) {
    p_int->Beam()->CalculateWeight(Q2);
    csi.AddWeight(p_int->Beam()->Weight());
  }
}

Event_Weights Single_Process::Differential(const Vec4D_Vector& p,
                                           Weight_Type type)
{
  DEBUG_FUNC(Name()<<", RS:"<<GetSubevtList());

  ResetResultsForDifferential(type);
  InitMEWeightInfo();
  UpdateIntegratorMomenta(p);
  CalculateFlux(p);

  if (m_zero) {
    m_eventweights = 0.0;
    return m_eventweights;
  }

  if (IsMapped()) {
    // NOTE: this needs to be reset before returning
    p_mapproc->SetCaller(this);
  }

  Scale_Setter_Base* scales {ScaleSetter(1)};

  Partonic(p);

  if (GetSubevtList() == nullptr) {

    if (type == Weight_Type::all) {
      assert(Selector()->CombinedResults().size() == 1);
      m_eventweights *= Selector()->CombinedResults()[0];
    }

    // calculate ISR weight
    m_csi = ClusterSequenceInfo(scales->Amplitudes(), scales->Scale(stp::fac));
    m_csi.AddFlux(m_lastflux);

    // update results
    m_last = m_lastxs + NfSchemeConversionTerms() - m_lastbxs * m_csi.m_ct;
    m_lastb = m_lastbxs;
    if (m_use_biweight) {
      double prefac {m_csi.m_pdfwgt * m_csi.m_flux};
      m_last *= prefac;
      m_lastb *= prefac;
    }

    if (p_mc != nullptr && m_dsweight && m_pinfo.Has(nlo_type::vsub)) {
      // calculate DADS term for MC@NLO: one PS point, many dipoles
      m_mewgtinfo.m_type |= mewgttype::DADS;

      if (m_dads) {

        // ask BBar emission generator for dipole parameters
        const Dipole_Params dps {p_mc->Active(this)};

        if (dps.p_dip != nullptr) {

          // calculate incoming parton longitudinal momentum fractions
          std::array<double, 2> x;
          for (size_t i {0}; i < 2; ++i) {
            x[i] = Min(p_int->ISR()->CalcX(dps.m_p[i]), 1.0);
          }

          for (Process_Base* proc : dps.m_procs) {

            // NOTE: the following adjustments need to be reset before returning
            const size_t mcmode {proc->SetMCMode(m_mcmode)};
            const bool lookup {proc->LookUp()};
            proc->SetLookUp(false);

            m_dadseventweights +=
                proc->Differential(dps.m_p, type) * dps.m_weight * m_dsweight;

            double dadsmewgt {proc->GetMEwgtinfo()->m_B * dps.m_weight *
                              m_dsweight};
            DADS_Info dads {-dadsmewgt,
                            x[0],
                            x[1],
                            proc->Flavours()[0],
                            proc->Flavours()[1]};
            m_mewgtinfo.m_dadsinfos.push_back(dads);

            // NOTE: here we reset the adjustments we have done above
            proc->SetLookUp(lookup);
            proc->SetMCMode(mcmode);
          }
        }
      }
    }

  } else {

    const auto triggers = Selector()->CombinedResults();

    for (int i {0}; i < GetSubevtList()->size(); ++i) {
      auto sub = (*GetSubevtList())[i];

      // append RDA info to m_mewgtinfo
      if (m_mewgtinfo.m_type & mewgttype::H) {
        RDA_Info rda {sub->m_mewgt,
                      sub->m_mu2[stp::ren],
                      sub->m_mu2[stp::fac],
                      sub->m_mu2[stp::fac],
                      sub->m_i,
                      sub->m_j,
                      sub->m_k};
        m_mewgtinfo.m_rdainfos.push_back(rda);
        msg_Debugging()<<i<<": wgt="<<m_mewgtinfo.m_rdainfos.back().m_wgt
                       <<std::endl;
      }

      // calculate weight for each subevent
      if (sub->m_me == 0.0) {
        sub->m_result = 0.0;
        sub->m_results = 0.0;
      } else {
        // calculate ISR weight
        ClusterAmplitude_Vector ampls(sub->p_real->p_ampl ? 1 : 0,
                                      sub->p_real->p_ampl);
        if (!ampls.empty()) {
          ampls.front()->SetProc(sub->p_proc);
        }
        Cluster_Sequence_Info csi {
            ClusterSequenceInfo(ampls, sub->m_mu2[stp::fac])};
        csi.AddFlux(m_lastflux);
        if (m_mewgtinfo.m_type & mewgttype::H) {
          m_mewgtinfo.m_rdainfos.back().m_csi = csi;
        }

        // update subevent information
        sub->m_result = sub->m_me * csi.m_pdfwgt * csi.m_flux;
        assert(!triggers.empty());
        const auto& jet_trigger_weights =
            (triggers.size() == 1) ? triggers[0] : triggers[i];
        sub->m_results = Event_Weights {};
        if (type == Weight_Type::all) {
          sub->m_results *= jet_trigger_weights;
        } else {
          sub->m_results *= jet_trigger_weights.Nominal();
        }
        sub->m_results *= sub->m_result;
        sub->m_mewgt *= m_lastflux;
        sub->m_xf1 = p_int->ISR()->XF1(0);
        sub->m_xf2 = p_int->ISR()->XF2(0);

        // update result
        m_last += (sub->m_trig & 1) ? sub->m_result : 0.0;
      }
    }

  }

  if (IsMapped()) {
    // NOTE: here we reset the mapped process caller we have modified above
    p_mapproc->SetCaller(p_mapproc);
  }

  UpdateMEWeightInfo(scales);

  // perform on-the-fly reweighting
  m_eventweights *= m_last;
  if (type != Weight_Type::nominal && m_eventweights.ContainsVariations()) {
    if (GetSubevtList() == nullptr) {
      RefactoredReweightBVI(type, scales->Amplitudes());
    } else {
      RefactoredReweightRS(type, scales->Amplitudes());
    }
  }
  m_last -= m_dadseventweights.Nominal();
  m_eventweights -= m_dadseventweights;

  // propagate (potentially) re-clustered momenta
  if (GetSubevtList() == nullptr) {
    UpdateIntegratorMomenta(scales->Amplitudes());
  } else {
    for (NLO_subevt* sub : *GetSubevtList()) {
      UpdateSubeventMomenta(*sub);
    }
  }

  return m_eventweights;
}

void Single_Process::ResetResultsForDifferential(Weight_Type type)
{
  m_lastb=m_last=m_lastflux=0.0;
  m_mewgtinfo.Reset();
  if (type == Weight_Type::all) {
    m_dadseventweights = m_eventweights = Event_Weights {};
  } else {
    m_dadseventweights = m_eventweights = Event_Weights {0, 1.0};
  }
  m_dadseventweights = 0.0;
}

void Single_Process::UpdateIntegratorMomenta(const Vec4D_Vector& p)
{
  p_int->SetMomenta(p);
  if (IsMapped())
    p_mapproc->Integrator()->SetMomenta(p);
}

void Single_Process::UpdateIntegratorMomenta(ClusterAmplitude_Vector& ampls)
{
  if (ampls.size()) {
    Cluster_Amplitude* ampl {ampls.front()->Last()};
    if (ampl->NLO() & 256) {
      Vec4D_Vector p(m_nin + m_nout);
      for (size_t i(0); i < ampl->Legs().size(); ++i)
        p[i] = i < m_nin ? -ampl->Leg(i)->Mom() : ampl->Leg(i)->Mom();
      p_int->SetMomenta(p);
    }
  }
}

void Single_Process::UpdateSubeventMomenta(NLO_subevt& sub)
{
  if (sub.p_real->p_ampl == nullptr) {
    return;
  }
  Cluster_Amplitude* ampl {sub.p_real->p_ampl->Last()};
  if (ampl->NLO() & 256) {
    for (size_t i(0); i < ampl->Legs().size(); ++i) {
      *((Vec4D*)&sub.p_mom[i]) = ampl->Leg(i)->Mom();
    }
    for (size_t i(ampl->Legs().size()); i < sub.m_n; ++i) {
      *((Vec4D*)&sub.p_mom[i]) = Vec4D();
    }
  }
}

void Single_Process::CalculateFlux(const Vec4D_Vector& p)
{
  if (m_nin == 1)
    m_lastflux = p_int->ISR()->Flux(p[0]);
  else
    m_lastflux = p_int->ISR()->Flux(p[0], p[1]);
  m_lastflux /= m_issymfac;
}

void Single_Process::RefactoredReweightBVI(Weight_Type type,
                                           ClusterAmplitude_Vector& ampls)
{
  BornLikeReweightingInfo info {m_mewgtinfo, ampls, m_last};
  m_eventweights.Apply(
      [this, &ampls, &info](double varweight,
                            size_t varindex,
                            Variation_Parameters& varparams) -> double {
        if (varweight == 0.0) {
          return 0.0;
        }
        double K {1.0};
        if (!m_mewgtinfo.m_bkw.empty()) {
          K = m_mewgtinfo.m_bkw[varindex];
        }
        if (m_mewgtinfo.m_type == mewgttype::none ||
            m_mewgtinfo.m_type == mewgttype::METS) {

          return K * ReweightBornLike(varparams, info);

        } else {

          const double muR2new {MuR2(varparams, info)};
          const Cluster_Sequence_Info csi {
              ClusterSequenceInfo(varparams, info, muR2new / info.m_muR2,
                                  &m_mewgtinfo.m_clusseqinfo)};

          double BVIKPnew {0.0};

          if (csi.m_pdfwgt != 0.0) {

            // calculate AlphaS factors (for Born and non-Born contributions)
            const double alphasratio {
              AlphaSRatio(info.m_muR2, muR2new, varparams.p_alphas)};
            const double alphasfac {pow(alphasratio, info.m_orderqcd)};
            double bornalphasfac {1.0};
            if (alphasfac != 1.0) {
              // for the Born contribution within BVIKP, we need to evaluate at
              // the lower order
              const bool needslowerorderqcd {
                  m_mewgtinfo.m_type & mewgttype::VI ||
                  m_mewgtinfo.m_type & mewgttype::KP};
              bornalphasfac =
                  needslowerorderqcd ? alphasfac / alphasratio : alphasfac;
            }

            // B term
            const double Bnew {m_mewgtinfo.m_B * bornalphasfac};

            // VI term
            const double logR {log(muR2new / info.m_muR2)};
            const double VInew {
                (m_mewgtinfo.m_VI + m_mewgtinfo.m_wren[0] * logR +
                 m_mewgtinfo.m_wren[1] * 0.5 * ATOOLS::sqr(logR)) *
                alphasfac};

            // KP terms
            const double KPnew {KPTerms(&varparams) * alphasfac};

            // Calculate K1
            double K1 {1.0};
            if (m_mewgtinfo.m_bkw.size() > s_variations->Size()) {
              K1 = m_mewgtinfo
                       .m_bkw[s_variations->Size() + varindex];
            }

            // Calculate final reweighted BVIKP result
            BVIKPnew = (Bnew * K * (1.0 - csi.m_ct) + (VInew + KPnew) * K1) *
                       csi.m_pdfwgt;
          }

          return BVIKPnew;
        }
      });
}

void Single_Process::RefactoredReweightRS(Weight_Type type,
                                          ClusterAmplitude_Vector& ampls)
{
  BornLikeReweightingInfo info {m_mewgtinfo, ampls, m_last};
  auto last_subevt_idx = GetSubevtList()->size() - 1;
  m_eventweights.Apply([this, &info, &last_subevt_idx](
                           double varweight,
                           size_t varindex,
                           Variation_Parameters& varparams) -> double {
    double K {1.0};
    if (!m_mewgtinfo.m_bkw.empty()) {
      K = m_mewgtinfo.m_bkw[varindex];
    }
    double weight {0.0};
    for (int i {0}; i <= last_subevt_idx; ++i) {
      auto sub = (*GetSubevtList())[i];
      if (sub->m_results.Variation(varindex) != 0.0) {
        info.m_wgt = sub->m_mewgt;
        info.m_muR2 = sub->m_mu2[stp::ren];
        info.m_muF2 = sub->m_mu2[stp::fac];
        info.m_ampls = ClusterAmplitude_Vector(sub->p_real->p_ampl ? 1 : 0,
                                               sub->p_real->p_ampl);
        info.m_fallbackresult = sub->m_result;
        auto contrib = K * ReweightBornLike(varparams, info);
        sub->m_results.Variation(varindex) = contrib;
        weight += contrib;
      }
    }
    return weight;
  });
}

void Single_Process::InitMEWeightInfo()
{
  m_mewgtinfo.m_oqcd = MaxOrder(0);
  m_mewgtinfo.m_oew  = MaxOrder(1);
  m_mewgtinfo.m_fl1  = (int)(Flavours()[0]);
  m_mewgtinfo.m_fl2  = (int)(Flavours()[1]);
  m_mewgtinfo.m_x1   = p_int->ISR()->X1();
  m_mewgtinfo.m_x2   = p_int->ISR()->X2();
}

void Single_Process::UpdateMEWeightInfo(Scale_Setter_Base* scales)
{
  m_mewgtinfo *= m_lastflux;
  if (scales != nullptr) {
    m_mewgtinfo.m_muf2 = scales->Scale(stp::fac);
    m_mewgtinfo.m_mur2 = scales->Scale(stp::ren);
  }
  m_mewgtinfo.m_clusseqinfo = m_csi;
  msg_Debugging()<<m_mewgtinfo;
}

double
Single_Process::ReweightBornLike(ATOOLS::Variation_Parameters& varparams,
                                 Single_Process::BornLikeReweightingInfo& info)
{
  if (info.m_wgt == 0.0) {
    return 0.0;
  }
  const double muR2new(MuR2(varparams, info));
  ATOOLS::Cluster_Sequence_Info csi(
      ClusterSequenceInfo(varparams, info, muR2new / info.m_muR2,
                          &m_mewgtinfo.m_clusseqinfo));
  if (csi.m_pdfwgt == 0.0) {
    return 0.0;
  } 
  const double alphasratio(AlphaSRatio(info.m_muR2, muR2new, varparams.p_alphas));
  const double alphasfac(pow(alphasratio, info.m_orderqcd));
  const double newweight(info.m_wgt * alphasfac * csi.m_pdfwgt);
  return newweight;
}

ATOOLS::Cluster_Sequence_Info Single_Process::ClusterSequenceInfo(
    ATOOLS::Variation_Parameters& varparams,
    Single_Process::BornLikeReweightingInfo & info,
    const double &mur2fac,
    const ATOOLS::Cluster_Sequence_Info * const nominalcsi)
{
  const double Q2(info.m_muF2 * varparams.m_muF2fac);

  // insert target PDF into ISR_Handler, such that ClusterSequenceInfo uses
  // them through the ISR_Handler instead of the nominal PDF
  PDF::PDF_Base *nominalpdf1 = p_int->ISR()->PDF(0);
  PDF::PDF_Base *nominalpdf2 = p_int->ISR()->PDF(1);
  p_int->ISR()->SetPDF(varparams.p_pdf1, 0);
  p_int->ISR()->SetPDF(varparams.p_pdf2, 1);

  ATOOLS::Cluster_Sequence_Info csi(
      ClusterSequenceInfo(info.m_ampls,
                          Q2, varparams.m_muF2fac, mur2fac,
                          varparams.m_showermuF2fac,
                          varparams.p_alphas,
                          nominalcsi));

  if (csi.m_pdfwgt == 0.0 && m_mewgtinfo.m_clusseqinfo.m_pdfwgt != 0.0) {
    varparams.IncrementOrInitialiseWarningCounter("Target PDF weight is zero, nominal is non-zero");
  } else if (csi.m_pdfwgt != 0.0 && m_mewgtinfo.m_clusseqinfo.m_pdfwgt == 0.0) {
    varparams.IncrementOrInitialiseWarningCounter("Target PDF weight is non-zero, nominal is zero");
  }

  // reset
  p_int->ISR()->SetPDF(nominalpdf1, 0);
  p_int->ISR()->SetPDF(nominalpdf2, 1);
  p_int->ISR()->SetMuF2(info.m_muF2, 0);
  p_int->ISR()->SetMuF2(info.m_muF2, 1);

  return csi;
}

double Single_Process::KPTerms(const ATOOLS::Variation_Parameters * varparams)
{
  double KP(KPTerms(0, varparams->p_pdf1,
                       varparams->p_pdf2, varparams->m_muF2fac) * m_lastflux);

  return KP;
}

double Single_Process::KPTerms
(int mode, PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb, double scalefac2)
{
  THROW(fatal_error,"Virtual function not reimplemented.");
  return 0.;
}

double Single_Process::MuR2(
  const ATOOLS::Variation_Parameters& varparams,
  Single_Process::BornLikeReweightingInfo & info) const
{
  double mu2new(info.m_muR2 * varparams.m_muR2fac);
  if ((varparams.m_showermuR2fac != 1.0) && (m_pinfo.m_ckkw & 1)) {
    ATOOLS::ClusterAmplitude_Vector &ampls = info.m_ampls;
    if (ampls.size()) {
      // go through cluster sequence
      double alphasnewproduct(1.0);
      double oqcdsum(0.0);
      double minmu2(1.0);
      Cluster_Amplitude *ampl(ampls.front());
      for (; ampl->Next(); ampl = ampl->Next()) {
        if (m_pinfo.Has(nlo_type::real) && ampl->Prev() == NULL) {
          continue;
        }
        double oqcd(ampl->OrderQCD() - ampl->Next()->OrderQCD());
        if (oqcd > 0.0) {
          double mu2(Max(ampl->Mu2(), MODEL::as->CutQ2()));
          mu2 = Min(mu2, sqr(rpa->gen.Ecms()));
          const double mu2new(mu2 * varparams.m_showermuR2fac);
          minmu2 = Min(minmu2, mu2new);
          const double alphasnew(varparams.p_alphas->BoundedAlphaS(mu2new));
          alphasnewproduct *= pow(alphasnew, oqcd);
          oqcdsum += oqcd;
        }
      }
      const double oqcdremainder(ampl->OrderQCD() - (m_pinfo.Has(nlo_type::vsub) ? 1 : 0));
      if (oqcdremainder) {
        const double mu2(Max(ampl->Mu2(), MODEL::as->CutQ2()));
        const double mu2new(mu2 * varparams.m_showermuR2fac);
        minmu2 = Min(minmu2, mu2new);
        const double alphasnew(varparams.p_alphas->BoundedAlphaS(mu2new));
        alphasnewproduct *= pow(alphasnew, oqcdremainder);
        oqcdsum += oqcdremainder;
      }
      if (oqcdsum) {
        // solve for new mu2
        const double alphasnewaverage(pow(alphasnewproduct, 1.0 / oqcdsum));
        const double maxmu2(varparams.m_showermuR2fac * 1.01 * sqr(rpa->gen.Ecms()));
        mu2new = MODEL::as->WDBSolve(alphasnewaverage, minmu2, maxmu2);
      }
    }
  }
  return mu2new;
}


double Single_Process::AlphaSRatio(
    double mur2old, double mur2new,
    MODEL::Running_AlphaS * asnew)
{
  const double alphasnew((*asnew)(mur2new));
  const double alphasold((*MODEL::as)(mur2old));
  return alphasnew / alphasold;
}


bool Single_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  p_int->Reset();
  auto psh = p_int->PSHandler();
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
  THROW(not_implemented, "To be implemented by child classes");
  return false;
}

const Flavour_Vector &Single_Process::
CombinedFlavour(const size_t &idij)
{
  THROW(not_implemented, "To be implemented by child classes");
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
  if (!(mode&1)) {
  MCatNLO_Process *mp(dynamic_cast<MCatNLO_Process*>(Parent()));
  if (mp) {
    Cluster_Amplitude *ampl(mp->GetAmplitude());
    return ampl;
  }
  }
  NLO_subevtlist *subs(GetRSSubevtList());
  if (subs==NULL) subs=GetSubevtList();
  if (subs && subs->back()->p_ampl) return subs->back()->p_ampl->CopyAll();
  ClusterAmplitude_Vector &ampls(ScaleSetter(1)->Amplitudes());
  if (ampls.size()) {
    msg_Debugging()<<METHOD<<"(): Found "
		   <<ampls.size()<<" amplitude(s) ... ";
    msg_Debugging()<<"select 1st.\n";
    for (Cluster_Amplitude *ampl(ampls.front());ampl;ampl=ampl->Next())
      msg_Debugging()<<*ampl<<"\n";
    return ampls.front()->CopyAll();
  }
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->SetNIn(m_nin);
  ampl->SetOrderQCD(m_maxcpl[0]);
  ampl->SetOrderEW(m_maxcpl[1]);
  ampl->SetProc(this);
  ampl->SetMS(p_gen);
  ampl->SetProcs(p_apmap);
  if (p_int->ColorIntegrator()==NULL) {
    for (size_t i(0);i<m_nin+m_nout;++i)
      ampl->CreateLeg(i<m_nin?-p_int->Momenta()[i]:p_int->Momenta()[i],
		      i<m_nin?m_flavs[i].Bar():m_flavs[i]);
    std::vector<int> tids, atids;
    for (size_t i(0);i<ampl->Legs().size();++i)
      if (ampl->Leg(i)->Flav().StrongCharge()>0) {
	tids.push_back(i);
	if (ampl->Leg(i)->Flav().StrongCharge()==8)
	  atids.push_back(i);
	else ampl->Leg(i)->SetCol(ColorID(ampl->Leg(i)->Col().m_i,0));
      }
      else if (ampl->Leg(i)->Flav().StrongCharge()<0) {
	ampl->Leg(i)->SetCol(ColorID(0,ampl->Leg(i)->Col().m_j));
	atids.push_back(i);
      }
      else {
	ampl->Leg(i)->SetCol(ColorID(0,0));
      }
    while (true) {
      std::random_shuffle(atids.begin(),atids.end(),*ran);
      size_t i(0);
      for (;i<atids.size();++i) if (atids[i]==tids[i]) break;
      if (i==atids.size()) break;
    }
    for (size_t i(0);i<tids.size();++i) {
      int cl(Flow::Counter());
      ampl->Leg(tids[i])->SetCol(ColorID(cl,ampl->Leg(tids[i])->Col().m_j));
      ampl->Leg(atids[i])->SetCol(ColorID(ampl->Leg(atids[i])->Col().m_i,cl));
    }
  }
  else {
    Int_Vector ci(p_int->ColorIntegrator()->I());
    Int_Vector cj(p_int->ColorIntegrator()->J());
    for (size_t i(0);i<m_nin+m_nout;++i)
      ampl->CreateLeg(i<m_nin?-p_int->Momenta()[i]:p_int->Momenta()[i],
		      i<m_nin?m_flavs[i].Bar():m_flavs[i],
		      ColorID(ci[i],cj[i]));
  }  
  ampl->SetMuF2(ScaleSetter(1)->Scale(stp::fac));
  ampl->SetMuR2(ScaleSetter(1)->Scale(stp::ren));
  ampl->SetMuQ2(ScaleSetter(1)->Scale(stp::res));
  return ampl;
}
