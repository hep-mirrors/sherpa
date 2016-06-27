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
#include "METOOLS/Explicit/NLO_Counter_Terms.H"
#include "MODEL/Main/Coupling_Data.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Single_Process::Single_Process() :
  m_lastbxs(0.0), m_lastflux(0.0), m_zero(false)
{
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

double Single_Process::CollinearCounterTerms(
    const int i,
    const ATOOLS::Flavour &fl,
    const ATOOLS::Vec4D &p,
    const double &z,
    const double &t1, const double &t2,
    const double &muf2fac,
    const double &mur2fac,
    MODEL::One_Running_AlphaS * as) const
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

  // determine ct
  double ct(0.0), lt(log(t1/t2)), x(p_int->ISR()->CalcX(p));
  msg_Debugging()<<asvalue<<"/(2\\pi) * log("<<sqrt(t1)<<"/"
		 <<sqrt(t2)<<") = "<<asvalue/(2.0*M_PI)*lt<<"\n";
  Flavour jet(kf_jet);
  double fb=p_int->ISR()->PDFWeight((1<<(i+1))|8,p,p,lmuf2,lmuf2,fl,fl,0);
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
    bool skipsfirstampl, const double &Q2,
    const double &muf2fac,
    const double &mur2fac,
    MODEL::One_Running_AlphaS * as)
{
  if (!m_use_biweight) {
    return 1.;
  }
  if (m_nin == 1) {
    return 0.5 / p_int->Momenta()[0].Mass();
  } else if (m_nin > 2) {
    THROW(not_implemented, "More than two incoming particles.");
  }
  Cluster_Sequence_Info csi;
  AddISR(csi, ampls, skipsfirstampl, Q2, muf2fac, mur2fac, as);
  AddBeam(csi, Q2);
  return csi;
}

void Single_Process::AddISR(ATOOLS::Cluster_Sequence_Info &csi,
            const ATOOLS::ClusterAmplitude_Vector &ampls,
            bool skipsfirstampl, const double &Q2,
            const double &muf2fac, const double &mur2fac,
            MODEL::One_Running_AlphaS * as)
{
  if (p_int->ISR()) {
    // add external PDF weight (before clustering)
    double pdfext(p_int->ISR()->PDFWeight(0,
                                          p_int->Momenta()[0],
                                          p_int->Momenta()[1],
                                          Q2, Q2,
                                          m_flavs[0], m_flavs[1]));
    msg_Debugging()<<"PDF(fla="<<m_flavs[0]
                   <<", xa="<<p_int->ISR()->CalcX(p_int->Momenta()[0])
                   <<", ta="<<Q2<<") * "
                   <<"PDF(flb="<<m_flavs[1]
                   <<", xb="<<p_int->ISR()->CalcX(p_int->Momenta()[1])
                   <<", tb="<<Q2<<") -> "<<pdfext<<std::endl;
    csi.AddWeight(pdfext);

    // add splittings and their PDF weight ratios from clustering
    if (ampls.size() && (m_pinfo.m_ckkw&1)) {
      DEBUG_FUNC(m_name<<", \\mu_F = "<<sqrt(Q2)
                 <<", #ampls="<<ampls.size());
      m_mewgtinfo.m_type|=mewgttype::METS;

      // add external splitting
      csi.AddSplitting(Q2,
                       p_int->ISR()->CalcX(p_int->Momenta()[0]),
                       p_int->ISR()->CalcX(p_int->Momenta()[1]),
                       m_flavs[0], m_flavs[1]);
      csi.AddPDFRatio(pdfext, 1.);

      Cluster_Amplitude *ampl(ampls.front());
      msg_IODebugging()<<*ampl<<"\n";
      if (skipsfirstampl) {
	ampl = ampl->Next();
	msg_IODebugging()<<*ampl<<"\n";
      }

      // add subsequent splittings
      bool addedfirstsplitting(false);
      double currentQ2(Q2);
      double pdfnum(pdfext), pdfden(pdfext);
      for (; ampl; ampl = ampl->Next()) {
        // skip decays, equal scales, unordered configs,
        // and quarks below threshold

        // skip decays (they are not even added to the splittings)
        msg_IODebugging()<<*ampl<<"\n";
	if (ampl->Next() && ampl->Next()->Splitter()->Stat() == 3) {
          msg_Debugging()<<"Skip. Decay "<<
            ID(ampl->Next()->Splitter()->Id())<<"\n";
          continue;
	}

        // add splitting
	Flavour f1(ampl->Leg(0)->Flav().Bar());
	Flavour f2(ampl->Leg(1)->Flav().Bar());
	if (MapProc() && LookUp() && !skipsfirstampl) {
	  f1=ReMap(f1, ampl->Leg(0)->Id());
	  f2=ReMap(f2, ampl->Leg(1)->Id());
	}
	csi.AddSplitting(ampl->KT2(),
			 p_int->ISR()->CalcX(-ampl->Leg(0)->Mom()),
			 p_int->ISR()->CalcX(-ampl->Leg(1)->Mom()),
			 f1, f2);

        // skip equal scales
	if (IsEqual(currentQ2, ampl->Next() ? ampl->KT2() : ampl->MuF2())) {
	  msg_Debugging()<<"Skip. Scales equal: t_i="<<currentQ2
			 <<", t_{i+1}="<<(ampl->Next()?ampl->KT2():ampl->MuF2())
			 <<std::endl;
	  if (ampl->Next() == NULL) {
            csi.AddPDFRatio(1., pdfden);
          } else {
            csi.AddPDFRatio(pdfnum, pdfden);
          }
	  continue;
	}

        // skip unordered configuration
	if (addedfirstsplitting && currentQ2 > ampl->KT2()) {
	  msg_Debugging()<<"Skip. Unordered history "<<
	    sqrt(currentQ2)<<" > "<<sqrt(ampl->KT2())<<"\n";
	  currentQ2 = sqrt(std::numeric_limits<double>::max());
	  continue;
	}

        // skip when a scale is below a (quark) mass threshold
	if (currentQ2 < sqr(2.0 * f1.Mass(true)) || currentQ2 < sqr(2.0 * f2.Mass(true))) {
	  msg_Debugging()<<"Skip. Quarks below threshold: t="<<currentQ2
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
	double lastQ2=currentQ2;
        currentQ2 = (ampl->Next() == NULL) ? Q2 : ampl->KT2();

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

        // add weights
	if (!IsZero(wn1) && !IsZero(wd1)) csi.AddWeight(wn1 / wd1);
	if (!IsZero(wn2) && !IsZero(wd2)) csi.AddWeight(wn2 / wd2);

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
	if (m_pinfo.Has(nlo_type::born)) {
	  double rn[2] = {ran->Get(), ran->Get()};
	  for (int i(0); i < 2; ++i) {
	    if (i == 0 && (IsZero(wn1) || IsZero(wd1))) continue;
	    if (i == 1 && (IsZero(wn2) || IsZero(wd2))) continue;
	    Vec4D p(-ampl->Leg(i)->Mom());
	    double x(p_int->ISR()->CalcX(p)), z(x + (1.0 - x) * rn[i]);
	    csi.AddCounterTerm(CollinearCounterTerms(i, i ? f2 : f1, p, z, currentQ2, lastQ2, muf2fac, mur2fac, as),
			       z, i);
	  }
	}
	addedfirstsplitting = true;
      }
    }
  }
}

void Single_Process::AddBeam(ATOOLS::Cluster_Sequence_Info &csi,
                             const double &Q2)
{
  if (p_int->Beam() && p_int->Beam()->On()) {
    p_int->Beam()->CalculateWeight(Q2);
    csi.AddWeight(p_int->Beam()->Weight());
  }
}

double Single_Process::Differential(const Vec4D_Vector &p)
{
  DEBUG_FUNC(Name()<<", RS:"<<GetSubevtList());
  m_lastb=m_last=m_lastflux=0.0;
  m_mewgtinfo.Reset();
  m_mewgtinfo.m_oqcd=MaxOrder(0);
  m_mewgtinfo.m_oew=MaxOrder(1);
  m_mewgtinfo.m_fl1=(int)(Flavours()[0]);
  m_mewgtinfo.m_fl2=(int)(Flavours()[1]);
  p_int->SetMomenta(p);
  if (IsMapped()) p_mapproc->Integrator()->SetMomenta(p);
  m_lastflux = p_int->ISR()->Flux(p[0],p[1]);
  if (GetSubevtList()==NULL) {
    if (m_zero) return 0.0;
    Scale_Setter_Base *scs(ScaleSetter(1));
    scs->SetCaller(Proc());
    if (Partonic(p,0)==0.0) return 0.0;
    m_mewgtinfo*=m_lastflux;
    m_mewgtinfo.m_muf2=scs->Scale(stp::fac);
    m_mewgtinfo.m_mur2=scs->Scale(stp::ren);
    if (m_lastxs==0.0) return m_last=0.0;
    m_last=m_lastxs;
    ClusterAmplitude_Vector ampls = scs->Amplitudes().size() ? 
        scs->Amplitudes() : ClusterAmplitude_Vector();
    const double facscale(scs->Scale(stp::fac));
    ATOOLS::Cluster_Sequence_Info csi(
        ClusterSequenceInfo(ampls, false, facscale));
    csi.AddFlux(m_lastflux);
    m_mewgtinfo.m_clusseqinfo=csi;
    msg_Debugging()<<m_mewgtinfo;
    m_last=(m_last-m_lastbxs*csi.m_ct)*
      (m_use_biweight?csi.m_pdfwgt*csi.m_flux:1.0);
    m_lastb=m_lastbxs*
      (m_use_biweight?csi.m_pdfwgt*csi.m_flux:1.0);
    if (p_mc!=NULL) {
      // calculate DADS for MC@NLO, one PS point, many dipoles
      msg_Debugging()<<"Calculating DADS terms"<<std::endl;
      m_mewgtinfo.m_type|=mewgttype::DADS;
      Dipole_Params dps(p_mc->Active(this));
      if (dps.p_dip!=NULL) {
        std::vector<double> x(2,-1.0);
        for (size_t j(0);j<2;++j) x[j]=Min(p_int->ISR()->CalcX(dps.m_p[j]),1.);
        for (size_t i(0);i<dps.m_procs.size();++i) {
          Process_Base *cp(dps.m_procs[i]);
          size_t mcmode(cp->SetMCMode(m_mcmode));
          bool lookup(cp->LookUp());
          cp->SetLookUp(false);
          // Set up reweighting within this dipole process, the result can
          // be retrieved later in the Reweighting function of this (main)
          // process
          if (p_variationweights) {
            SHERPA::Variations *v = p_variationweights->GetVariations();
            SHERPA::Variation_Weights * vw = new SHERPA::Variation_Weights(v);
            cp->SetOwnedVariationWeights(vw);
          }
          double dadswgt(cp->Differential(dps.m_p)*dps.m_weight);
          msg_Debugging()<<"DADS_"<<i<<" = "<<-dadswgt<<std::endl;
          double dadsmewgt(cp->GetMEwgtinfo()->m_B*dps.m_weight);
          DADS_Info dads(-dadsmewgt,x[0],x[1],
                         cp->Flavours()[0],cp->Flavours()[1]);
          msg_Debugging()<<dads<<std::endl;
          m_mewgtinfo.m_dadsinfos.push_back(dads);
          if (p_variationweights && dadsmewgt != 0.0) {
            *cp->VariationWeights() *= dps.m_weight;
          }
          m_last-=dadswgt;
          cp->SetLookUp(lookup);
          cp->SetMCMode(mcmode);
        }
      }
    }
    if (p_variationweights) {
      p_variationweights->UpdateOrInitialiseWeights(&Single_Process::ReweightWithoutSubevents,
                                                    *this, ampls);
    }
    return m_last;
  }
  else {
    Partonic(p,0);
    NLO_subevtlist *subs(GetSubevtList());
    for (size_t i(0);i<subs->size();++i) {
      NLO_subevt *sub((*subs)[i]);
      if (sub->m_me==0.0) sub->m_result=0.0;
      if (m_mewgtinfo.m_type&mewgttype::H) {
        RDA_Info rda(sub->m_mewgt,sub->m_mu2[stp::ren],
                     sub->m_mu2[stp::fac],sub->m_mu2[stp::fac],
                     sub->m_i,sub->m_j,sub->m_k);
        m_mewgtinfo.m_rdainfos.push_back(rda);
        msg_Debugging()<<i<<": wgt="<<m_mewgtinfo.m_rdainfos.back().m_wgt
                       <<std::endl;
      }
      if (sub->m_me!=0.0) {
        ClusterAmplitude_Vector ampls(sub->p_ampl?1:0,sub->p_ampl);
        if (ampls.size()) ampls.front()->SetProc(sub->p_proc);
        const double facscale(sub->m_mu2[stp::fac]);
        ATOOLS::Cluster_Sequence_Info csi(
            ClusterSequenceInfo(ampls, true, facscale));
        msg_Debugging()<<"fa*fb="<<csi.m_pdfwgt<<std::endl;
        csi.AddFlux(m_lastflux);
        if (m_mewgtinfo.m_type&mewgttype::H)
          m_mewgtinfo.m_rdainfos.back().m_csi=csi;
        sub->m_result=sub->m_me*csi.m_pdfwgt*csi.m_flux;
        sub->m_mewgt*=m_lastflux;
        // store core PDF values
        sub->m_xf1=p_int->ISR()->XF1(0);
        sub->m_xf2=p_int->ISR()->XF2(0);
      }
    }
    Scale_Setter_Base *scs(ScaleSetter(1));
    if (scs!=NULL) {
      m_mewgtinfo.m_muf2=scs->Scale(stp::fac);
      m_mewgtinfo.m_mur2=scs->Scale(stp::ren);
    }
    m_mewgtinfo*=m_lastflux;
    for (size_t i=0;i<subs->size();++i) {
      m_last+=(*subs)[i]->m_result;
    }
    if (p_variationweights) {
      const long emptyadditionaldata(NULL);
      p_variationweights->InitialiseWeights(&Single_Process::ReweightSubevents,
                                            *this, emptyadditionaldata);
    }
    return m_last;
  }
  THROW(fatal_error,"Internal error.");
  return 0.;
}

SHERPA::Subevent_Weights_Vector Single_Process::ReweightSubevents(
  SHERPA::Variation_Parameters * varparams,
  SHERPA::Variation_Weights * varweights,
  const long &additionaldata)
{
  NLO_subevtlist *sevtlist = GetSubevtList();
  if (!sevtlist) THROW(fatal_error, "Missing subevents.");

  // obtain ME wgt info and subevents fetch common info
  BornLikeReweightingInfo info;
  info.m_orderqcd = m_mewgtinfo.m_oqcd;
  info.m_fl1 = m_mewgtinfo.m_fl1;
  info.m_fl2 = m_mewgtinfo.m_fl2;
  info.m_x1 = p_int->ISR()->X1();
  info.m_x2 = p_int->ISR()->X2();
  SHERPA::Subevent_Weights_Vector weights;
  for (size_t i(0); i < sevtlist->size(); ++i) {
    NLO_subevt *sub((*sevtlist)[i]);

    // fetch subevent-specific info, then reweight
    info.m_wgt = sub->m_mewgt;
    info.m_muR2 = sub->m_mu2[stp::ren];
    info.m_muF12 = sub->m_mu2[stp::fac];
    info.m_muF22 = info.m_muF12;
    info.m_ampls = ClusterAmplitude_Vector(sub->p_ampl ? 1 : 0,
                                           sub->p_ampl);
    info.m_skipsfirstampl = true;
    info.m_fallbackresult = sub->m_result;
    weights.push_back(ReweightBornLike(varparams, info));
  }
  return weights;
}

double Single_Process::ReweightWithoutSubevents(
  SHERPA::Variation_Parameters * varparams,
  SHERPA::Variation_Weights * varweights,
  ATOOLS::ClusterAmplitude_Vector & ampls)
{
  // build type minus METS
  mewgttype::code nometstype((m_mewgtinfo.m_type & mewgttype::METS) ?
                             m_mewgtinfo.m_type ^ mewgttype::METS : m_mewgtinfo.m_type);

  if (GetSubevtList()) {
    THROW(fatal_error, "Unexpected subevents.");
  }

  // fetch common info
  BornLikeReweightingInfo info;
  info.m_orderqcd = m_mewgtinfo.m_oqcd;
  info.m_fl1 = m_mewgtinfo.m_fl1;
  info.m_fl2 = m_mewgtinfo.m_fl2;
  info.m_x1 = p_int->ISR()->X1();
  info.m_x2 = p_int->ISR()->X2();
  info.m_fallbackresult = m_last; 

  if (nometstype==mewgttype::none) { // non-NLO Born
    info.m_wgt = m_mewgtinfo.m_B;
    info.m_muR2 = m_mewgtinfo.m_mur2;
    info.m_muF12 = m_mewgtinfo.m_muf2;
    info.m_muF22 = info.m_muF12;
    info.m_ampls = ampls;
    info.m_skipsfirstampl = false;
    return ReweightBornLike(varparams, info);

  } else { // NLO Born, Virtual Integrated, KP, DADS
    // calculate factors common to all these contributions
    info.m_muR2 = m_mewgtinfo.m_mur2;
    info.m_muF12 = m_mewgtinfo.m_muf2;
    info.m_muF22 = info.m_muF12;
    info.m_ampls = ampls;
    info.m_skipsfirstampl = false;

    const ATOOLS::Cluster_Sequence_Info csi(ClusterSequenceInfo(varparams, info));
    if (csi.m_pdfwgt == 0.0) {
      varweights->IncrementOrInitialiseWarningCounter("Single process different PDF cut-off");
      return info.m_fallbackresult;
    } 

    // calculate AlphaS factors (for Born and non-Born contributions)
    const std::vector<double> alphasratios(AlphaSRatios(varparams, info));
    double alphasfac(1.0);
    for (std::vector<double>::const_iterator it(alphasratios.begin());
        it != alphasratios.end(); it++) {
      alphasfac *= *it;
    }
    double bornalphasfac(1.0);
    if (alphasfac != 1.0) {
      // for the Born contribution within BVIKP, we need to evaluate at the lower order
      // divide out one of the core AlphaS ratios
      // NOTE: The assumption that the last amplitude in a cluster
      // sequence has to be divided out, is most certainly wrong and this will
      // lead to an error in MEPS@NLO/MENLOPS runs, where we reweight the
      // AlphaS argument for individual splittings
      const bool needslowerorderqcd(nometstype & mewgttype::VI || nometstype & mewgttype::KP);
      bornalphasfac = needslowerorderqcd ? alphasfac / alphasratios.back() : alphasfac;
    }

    // Born
    const double Bnew(m_mewgtinfo.m_B * bornalphasfac);

    // Virtual Integrated
    double VInew(0.0);
    {
      const double logR(log(varparams->m_muR2fac));
      VInew = (m_mewgtinfo.m_VI
               + m_mewgtinfo.m_wren[0]*logR
               + m_mewgtinfo.m_wren[1]*0.5*ATOOLS::sqr(logR)) * alphasfac;
    }

    // KP
    const double KPnew(KPTerms(varparams) * alphasfac);

    // apply NLO counterterms
    double BVIKPnew(((Bnew + VInew + KPnew) - Bnew * csi.m_ct) * csi.m_pdfwgt);

    // DADS
    double DADSnew(0.0);
    {
      if (p_mc!=NULL) {
        Dipole_Params dps(p_mc->Active(this));
        if (dps.p_dip!=NULL) {
          for (size_t i(0);i<dps.m_procs.size();++i) {
            Process_Base *cp(dps.m_procs[i]);
            // when the dipole's Partonic returns 0, the variation weights will be
            // empty, therefore check before retrieval
            SHERPA::Variation_Weights *dipvarweights(cp->VariationWeights());
            if (dipvarweights->GetNumberOfVariations() > 0) {
              size_t paramindex = varweights->CurrentParametersIndex();
              DADSnew -= dipvarweights->GetVariationWeightAt(paramindex);
            }
          }
        }
      }
    }

    return BVIKPnew + DADSnew; 
  }
}

double Single_Process::ReweightBornLike(
  SHERPA::Variation_Parameters * varparams,
  Single_Process::BornLikeReweightingInfo & info)
{
  if (info.m_wgt == 0.0) {
    return 0.0;
  }
  ATOOLS::Cluster_Sequence_Info csi(ClusterSequenceInfo(varparams, info));
  if (csi.m_pdfwgt == 0.0) {
    p_variationweights->IncrementOrInitialiseWarningCounter("Single process different PDF cut-off");
    return info.m_fallbackresult;
  } 
  const std::vector<double> alphasratios(AlphaSRatios(varparams, info));
  double alphasfac(1.0);
  for (std::vector<double>::const_iterator it(alphasratios.begin());
      it != alphasratios.end(); it++) {
    alphasfac *= *it;
  }
  const double newweight(info.m_wgt * alphasfac * csi.m_pdfwgt);
  return newweight;
}

std::pair<double, double> Single_Process::GetPairOfPDFValuesOrOne(
    SHERPA::Variation_Parameters * varparams,
    Single_Process::BornLikeReweightingInfo & info) const
{
  const double muF12new(info.m_muF12 * varparams->m_muF2fac);
  const double muF22new(info.m_muF22 * varparams->m_muF2fac);
  double fa(1.0);
  if (varparams->p_pdf1) {
    varparams->p_pdf1->Calculate(info.m_x1, muF12new);
    fa = varparams->p_pdf1->GetXPDF(info.m_fl1) / info.m_x1;
  }
  double fb(1.0);
  if (varparams->p_pdf2) {
    varparams->p_pdf2->Calculate(info.m_x2, muF22new);
    fb = varparams->p_pdf2->GetXPDF(info.m_fl2) / info.m_x2;
  }
  return std::make_pair(fa, fb);
}

ATOOLS::Cluster_Sequence_Info Single_Process::ClusterSequenceInfo(
    SHERPA::Variation_Parameters * varparams,
    Single_Process::BornLikeReweightingInfo & info)
{
  // calculate Q2
  if (info.m_muF12 != info.m_muF22) {
    THROW(not_implemented, "Unequal factorisation scale factors.");
  }
  const double Q2(info.m_muF12 * varparams->m_muF2fac);

  // insert target PDF into ISR_Handler, such that ClusterSequenceInfo uses
  // them through the ISR_Handler instead of the nominal PDF
  PDF::PDF_Base *nominalpdf1 = p_int->ISR()->PDF(0);
  PDF::PDF_Base *nominalpdf2 = p_int->ISR()->PDF(1);
  p_int->ISR()->SetPDF(varparams->p_pdf1, 0);
  p_int->ISR()->SetPDF(varparams->p_pdf2, 1);

  ATOOLS::Cluster_Sequence_Info csi(
      ClusterSequenceInfo(info.m_ampls, info.m_skipsfirstampl,
                          Q2, varparams->m_muF2fac, varparams->m_muR2fac,
                          varparams->p_alphas));

  // reset
  p_int->ISR()->SetPDF(nominalpdf1, 0);
  p_int->ISR()->SetPDF(nominalpdf2, 1);
  p_int->ISR()->SetMuF2(info.m_muF12, 0);
  p_int->ISR()->SetMuF2(info.m_muF12, 1);

  return csi;
}

double Single_Process::KPTerms(SHERPA::Variation_Parameters * varparams)
{
  // insert target PDF into ISR_Handler, such that KP_Terms uses them through
  // the ISR_Handler instead of the nominal PDF
  PDF::PDF_Base *nominalpdf1 = p_int->ISR()->PDF(0);
  PDF::PDF_Base *nominalpdf2 = p_int->ISR()->PDF(1);
  p_int->ISR()->SetPDF(varparams->p_pdf1, 0);
  p_int->ISR()->SetPDF(varparams->p_pdf2, 1);

  double KP(KPTerms(0, varparams->m_muF2fac) * m_lastflux);

  // reset
  p_int->ISR()->SetPDF(nominalpdf1, 0);
  p_int->ISR()->SetPDF(nominalpdf2, 1);

  return KP;
}

std::vector<double> Single_Process::AlphaSRatios(
  SHERPA::Variation_Parameters * varparams,
  Single_Process::BornLikeReweightingInfo & info) const
{
  std::vector<double> ratios;
  if (varparams->m_showermuR2fac != 1.0) {
    ATOOLS::ClusterAmplitude_Vector &ampls = info.m_ampls;
    if (ampls.size() && (m_pinfo.m_ckkw & 1)) {
      // go through cluster sequence
      for (Cluster_Amplitude *ampl(ampls.front()); ampl; ampl = ampl->Next()) {
        const size_t power = ampl->Next() ?
          ampl->OrderQCD() - ampl->Next()->OrderQCD() : ampl->OrderQCD();
        if (power > 0) {
          const double mu2(Max(ampl->Mu2(), MODEL::as->CutQ2()));
          const double mu2new(mu2 * varparams->m_showermuR2fac);
          const double alphasold(MODEL::as->BoundedAlphaS(mu2));
          const double alphasnew(varparams->p_alphas->BoundedAlphaS(mu2new));
          const double alphasratio(alphasnew / alphasold);
          for (size_t i(0); i < power; i++) {
            ratios.push_back(alphasratio);
          }
        }
      }
    }
  } else {
    const double muR2new(info.m_muR2 * varparams->m_muR2fac);
    const double alphasnew((*varparams->p_alphas)(muR2new));
    const double alphasold((*MODEL::as)(info.m_muR2));
    const double alphasratio(alphasnew / alphasold);
    for (size_t i(0); i < info.m_orderqcd; i++) {
      ratios.push_back(alphasratio);
    }
  }
  return ratios;
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
