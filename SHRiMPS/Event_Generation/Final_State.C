#include "SHRiMPS/Event_Generation/Final_State.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

// asform:         constant=-1, frozen=0, smooth=1, IR0=2
// ktform:         frozen=0, smooth=1, cut=-1
// ordering:       rap_only=0, keep=1, ao=2, ao_keep=3, rap_phys=4, ao_phys=6

Final_State::Final_State(const int & test) :
  m_ktform(MBpars.KTForm()),
  m_ordering(MBpars.Ordering()), 
  m_resc_ktmin(MBpars.RescKTMin()),
  m_resc_nosing(MBpars.RescNoSing()),
  m_ktmin_mode(MBpars("KTMin_Mode")),
  //p_alphaS(static_cast<MODEL::Strong_Coupling *>
  //	   (MODEL::s_model->GetScalarFunction(std::string("strong_cpl")))), 
  p_alphaS(new Strong_Coupling(static_cast<Running_AlphaS *>
			       (s_model->GetScalarFunction(string("alpha_S"))),
				MBpars.AsForm(),MBpars("Q_as2"))),
  m_Q02(MBpars("Q02")), m_Q12(MBpars("Q12")), m_QN2(MBpars("QN2")), 
  m_d2(MBpars("Ddiff2")), m_kdiff(MBpars("kdiff")), 
  m_test(test), m_output(true)
{
  InitHistograms();
}


bool Final_State::FirstSinglet(const double & y1,const double & y2,
			       const double & sup,const int & nbeam) {
//  return false;
  if (p_ladder->IsRescatter() && m_resc_nosing==resc_nosing::on) return false;
  double wt1 = p_eikonal->SingletWeight(m_b1,m_b2,y1,y2,sup,nbeam); 
  double wt8 = p_eikonal->OctetWeight(m_b1,m_b2,y1,y2,sup,nbeam); 
  if (wt1>(wt1+wt8)*ran->Get()) {
    p_ladder->GetPropsBegin()->m_col = colour_type::singlet;
    p_ladder->SetDiffractive(true);
    return true;
  }
  return false;
}

double Final_State::
operator()(Ladder * ladder,const double & Deltay,
	   const bool & firstladder,const bool & firstattempt) {
  msg_Tracking()<<"--------------------------------------------\n"
		<<"--------------------------------------------\n"
		<<"--------------------------------------------\n"
		<<METHOD<<" for \n"<<(*ladder);
  p_ladder    = ladder;
  p_emissions = p_ladder->GetEmissions();
  m_plusiter  = p_emissions->begin();
  m_minusiter = p_emissions->end(); m_minusiter--;

  double y0(m_plusiter->first), y1(m_minusiter->first);
  m_lastwt = 1.;
  int nbeam(!p_ladder->IsRescatter()?2:
	    int(dabs(p_ladder->GetIn1()->m_mom.Y())>6.)+
	    int(dabs(p_ladder->GetIn2()->m_mom.Y())>6.));
  if (firstattempt && FirstSinglet(y0,y1,1.,nbeam)) {
    m_firstsing++;
    m_histomap[std::string("Delta_final")]->Insert(1./dabs(y0-y1)); 
    p_ladder->GetProps()->begin()->m_col = colour_type::singlet;
    p_ladder->SetDiffractive(true);
  }
  else {
    m_pprop    = -p_ladder->GetIn1()->m_mom;
    m_mprop    =  p_ladder->GetIn2()->m_mom;
    
    p_props    = p_ladder->GetProps();
    m_propiter = p_props->begin();
    
    if (m_resc_ktmin==resc_ktmin::off || !p_ladder->IsRescatter()) 
      m_kt1cut2 = 0.;
    else 
      m_kt1cut2 = sqrt(m_pprop.PPerp2()*m_mprop.PPerp2()); 
    m_kt1max2   = 0.;
    m_Deltay    = Deltay;
    m_lastwt    = GenerateEmissions();
  }
  return m_lastwt;
}

double Final_State::GenerateEmissions() {
  msg_Tracking()
    <<"### "<<METHOD<<"(kt1cut2="<<m_kt1cut2<<", kt1min2="<<m_kt1min2<<").\n";
  LadderMap::iterator split, spect;
  bool   run;
  double kt2(-1.);
  do {
    run = false;
    if (OneEmission(split,spect,kt2)) {
      split->second.m_mom    = m_k0;
      spect->second.m_mom    = m_k2;
      Ladder_Particle part(Flavour(kf_gluon),m_k1);
      p_ladder->AddParticle(m_k1.Y(),part);
      p_ladder->SetMaxKT2(kt2);
      if (FixPropColours(split,spect)) run = true;
      else {
	m_singexit++;
        m_histomap[std::string("Deltay_singexit")]->Insert(dabs(m_k2.Y()-m_k1.Y()));
	msg_Tracking()
	  <<"   active interval in singlet for "
	  <<p_ladder->Size()<<", lastwt="<<m_lastwt<<".\n";
	break;
      }
    }
  } while (run);

  double Delta = 
    double(p_ladder->Size()-1)/
    dabs(p_ladder->GetEmissionsBegin()->second.m_mom.Y()-
	 p_ladder->GetEmissionsRBegin()->second.m_mom.Y());

  m_histomap[std::string("Delta_final")]->Insert(Delta); 
  m_histomap[std::string("FSWt")]->Insert(m_lastwt);

  bool check(true);
  TPropList::iterator tpit=p_ladder->GetPropsBegin(),tbef(tpit);
  tpit++;
  while (tpit!=p_ladder->GetPropsEnd()) {
    if (tbef->m_col==colour_type::singlet && 
	tpit->m_col==colour_type::singlet) {
      msg_Error()<<"Error in "<<METHOD<<": "
		 <<"found two adjacent singlet propagators in\n"
		 <<(*p_ladder)
		 <<"   will replace one with an octet.\n";
      (tbef->m_q2<tpit->m_q2?tbef:tpit)->m_col=colour_type::octet;
      m_cols++;
      check = false;
    }
    tpit++;
    tbef++;
  }
  if (check) msg_Tracking()<<METHOD<<" yields ladder with tested colours:\n"
			   <<(*p_ladder);
  
  return m_lastwt;
}


bool Final_State::
OneEmission(LadderMap::iterator & split,LadderMap::iterator & spect,
	    double & kt2) {
  bool dir(SelectSplitterAndSpectator(split,spect)), test;
  if (dir) m_dir1++; else m_dir0++;
  m_k0 = split->second.m_mom;
  m_k2 = spect->second.m_mom;
  if (!dir) { 
    m_k0[3]*=-1.; m_k2[3]*=-1.; 
    m_q0[3]*=-1.; m_q2[3]*=-1; 
  }
  test = TryEmission(kt2,dir);
  if (!dir) { 
    m_k0[3]*=-1.; m_k1[3]*=-1.; m_k2[3]*=-1.; 
    m_q0[3]*=-1.; m_q2[3]*=-1; 
  }
  return test;
}

bool Final_State::
TryEmission(double & kt12,const bool & dir) {
  m_trials++;
  Vec4D Q(m_k0+m_k2);
  double Q2(Q.Abs2()),kt1max2(Q2/4.);
  switch (m_ktmin_mode) {
  case 1:
    m_kt1min2 = m_propiter->m_qt2;
    break;
  case 0:
  default:
    m_kt1min2 = 0.;
    break;
  }
  double kt1min2=Max(m_kt1min2,m_kt1cut2);
  if (kt1max2<kt1min2) {
    m_histomap[std::string("Delta_final")]->Insert(0.); 
    return false;
  }
  double ktexpo(1.);
  double QT(Q.PPerp()),MT(sqrt(Q.Abs2()+Q.PPerp2()));
  double Y(Q.Y()), Phi(Q.Phi()), cphi(cos(Phi)), sphi(sin(Phi));
  double expy(std::exp(Y));
  double kt2old(m_k2.PPerp()), y2old(m_k2.Y()), phi2old(m_k2.Phi());
  double kt0old(m_k0.PPerp()), y0old(m_k0.Y());
  double ystop(y2old);
  double kt2(kt2old), y2(y2old), phi2(phi2old), kt02, kt0, y0(y0old);

  double y1(m_k0.Y()),mu01_2(0.),mu12_2(0.),kt1,phi1;
  double term1(MT*cosh(Y-y2old)-QT*cos(Phi-phi2old)), term2, term3;
  double cphi0,sphi0,cphi1,sphi1,cphi2(cos(phi2)),sphi2(sin(phi2));
  double expy0,expy1,expy2(std::exp(y2));
  double wt,expo,rarg,reggewt(0.),sup(0.),recombwt(0.),diffwt(0.);
  double kmrwt(0.);
  size_t ntrials(0);
  
  double colfac(3.);
  double Delta(colfac*p_alphaS->MaxValue()/M_PI *
	       KT2integral(kt1max2,kt1min2,m_Q02eff,ktexpo));
  m_histomap[std::string("Delta_naive")]->Insert(Delta); 
  if (dir) m_histomap[std::string("Delta_naive1")]->Insert(Delta); 
  else m_histomap[std::string("Delta_naive0")]->Insert(Delta);

  Vec4D k_0, k_1, k_2, q01, q12;
  msg_Tracking()
    <<"------------- ["<<y0old<<", "<<y2old<<" - ["<<kt0old<<", "<<kt2old<<"]\n"
    <<"------------- cms: Y = "<<Y<<", MT = "<<MT<<", "
    <<"Q = "<<Q<<", Phi = "<<Phi<<".\n"
    <<"------------- Delta = "<<Delta<<"("<<kt1min2<<" -> "<<kt1max2<<"), "
    <<"alphaS = "<<p_alphaS->MaxValue()<<", "
    <<"colfac = "<<colfac<<", dir = "<<dir<<".\n"<<(*p_ladder)<<"\n";
  do {
    wt = 0.;
    if (ntrials++>10000) {
      m_histomap[std::string("Delta_final")]->Insert(0.); 
      return false;
    }
    double deltay1(-1.*log(ran->Get())/Delta);
    y1 += deltay1;
    m_histomap[std::string("Deltay_test")]->Insert(deltay1); 
    if (y1>ystop || y1<y0old) {
      m_ys++;
      m_histomap[std::string("Deltay_regexit")]->Insert(dabs(y2old-y0old));
      m_histomap[std::string("Delta_final")]->Insert(0.); 
      return false;
    }
    y0      = y0old;
    y2      = y2old;
    kt12    = SelectKT2(kt1max2,kt1min2,m_Q02eff,ktexpo);
    m_histomap[std::string("KT2_test1")]->Insert(kt12); 
    m_histomap[std::string("KT2_test2")]->Insert(kt12); 
    if (p_alphaS->Weight(kt12)<ran->Get()) {
      m_alphaS++;
      continue;
    }
    kt1     = sqrt(kt12);
    m_histomap[std::string("KT2_accept1")]->Insert(kt12); 
    m_histomap[std::string("KT2_accept2")]->Insert(kt12); 
    if (IsNan(kt1)) {
      m_rej_negkt1++;
      continue;
    }
    phi1    = ran->Get()*2.*M_PI;
    cphi1   = cos(phi1);
    sphi1   = sin(phi1);
    term2   = kt1*(cosh(y2old-y1)-cos(phi2old-phi1));
    term3   = MT*kt1*cosh(Y-y1)-QT*kt1*cos(Phi-phi1);
    if (dabs(term1-term2)<1.e-8) {
      m_rej_nokt2++;
      continue;
    }
    kt2     = (Q2/2.-term3)/(term1-term2);
    if (kt2<0) {
      m_rej_nokt2++;
      continue;
    }
    kt02    =
      QT*QT+kt2*kt2+kt1*kt1-
      2.*QT*kt2*cos(Phi-phi2)-2.*QT*kt1*cos(Phi-phi1)+
      2.*kt2*kt1*cos(phi2-phi1);
    if (kt02<1.e-12) {
      m_rej_nokt0++;
      continue;
    }
    kt0     = sqrt(kt02);
    cphi0   = (QT*cphi-kt2*cphi2-kt1*cphi1)/kt0;
    sphi0   = (QT*sphi-kt2*sphi2-kt1*sphi1)/kt0;
    
    if (dabs(cphi0)>1. || dabs(sphi0)>1.) {
      m_rej_nophi0++;
      continue;
    }
    if (dabs(cphi0*cphi0+sphi0*sphi0-1.)>1.e-6) {
      m_rej_nophi0++;
      continue;
    }
    expy1   = std::exp(y1);
    expy0   = (MT*expy-kt2*expy2-kt1*expy1)/kt0;
    if (expy0<0.) {
      m_rej_noy0++;
      continue;
    }
    y0      = log(expy0);
    double coshy0(cosh(y0)),coshy1(cosh(y1)),coshy2(cosh(y2)); 
    k_0     = kt0*Vec4D(coshy0,cphi0,sphi0,(y0>0.?1.:-1.)*sqrt((coshy0+1.)*(coshy0-1.)));
    k_1     = kt1*Vec4D(coshy1,cphi1,sphi1,(y1>0.?1.:-1.)*sqrt((coshy1+1.)*(coshy1-1.)));
    k_2     = kt2*Vec4D(coshy2,cphi2,sphi2,(y2>0.?1.:-1.)*sqrt((coshy2+1.)*(coshy2-1.)));
/*    k_0     = kt0*Vec4D(cosh(y0),cphi0,sphi0,sinh(y0));
    k_1     = kt1*Vec4D(cosh(y1),cphi1,sphi1,sinh(y1));
    k_2     = kt2*Vec4D(cosh(y2),cphi2,sphi2,sinh(y2));*/
    k_0[0]  = k_0.P();
    k_1[0]  = k_1.P();
    k_2[0]  = k_2.P();
    if ( k_0[3]>=k_0[0] || k_1[3]>=k_1[0] || k_2[3]>=k_2[0]) {
      msg_Error()<<METHOD<<": Reject emission due to inaccuracy in four-momentum\n"
                 <<"k_0 = "<<k_0<<" with y = "<<k_0.Y()<<"\n"
                 <<"k_1 = "<<k_1<<" with y = "<<k_1.Y()<<"\n"
                 <<"k_2 = "<<k_2<<" with y = "<<k_2.Y()<<std::endl;
      continue;
    }
    if (MomViolation(k_0,k_1,k_2,dir)) continue;
    m_histomap[std::string("Delta_kin")]->Insert(1./Max(1.e-2,dabs(y1-y0))); 
    q01   = dir?m_q0+k_0:m_q0-k_0;
    q12   = dir?-m_q2+k_2:-m_q2-k_2;

    m_q01_2 = dabs(q01.Abs2());
    m_q12_2 = dabs(q12.Abs2());
    if (dir) {
      m_histomap[std::string("q01_2_1")]->Insert(m_q01_2);
      m_histomap[std::string("q12_2_1")]->Insert(m_q12_2);
    }
    else {
      m_histomap[std::string("q01_2_0")]->Insert(m_q01_2);
      m_histomap[std::string("q12_2_0")]->Insert(m_q12_2);
    }
    m_histomap[std::string("Delta_order")]->Insert(1./Max(1.e-2,dabs(y1-y0))); 
    if (!IsOrdered(dir,k_0,k_1,k_2,m_q01_2)) continue;
    m_histomap[std::string("Delta_order")]->Insert(1./Max(1.e-2,dabs(y1-y0))); 
    wt   = 1.;
    double deltay(dabs(k_1.Y()-k_0.Y()));
    mu01_2 = Q02((k_0.Y()+k_1.Y())/2.);
    mu12_2 = Q02((k_1.Y()+k_2.Y())/2.);
    if (dir) m_histomap[std::string("ytest1")]->Insert(y1);
        else m_histomap[std::string("ytest0")]->Insert(y1);
    if ((MBpars.LadderWeight()==ladder_weight::Regge || 
	 MBpars.LadderWeight()==ladder_weight::ReggeDiffusion) && 
	deltay>m_Deltay) { 
      rarg = Min(mu01_2/q01.PPerp2(),q01.PPerp2()/mu01_2);
      expo = colfac*(*p_alphaS)(q01.PPerp2())*deltay/M_PI; 
//       rarg = Min(mu01_2/q01.Abs2(),q01.Abs2()/mu01_2);
//       expo = colfac*(*p_alphaS)(q01.Abs2())*deltay/M_PI; 
      //expo = colfac*p_alphaS->MaxValue()*deltay/M_PI; 
      wt  *= reggewt = pow(rarg,expo);
      m_histomap[std::string("ReggeWt")]->Insert(reggewt);
    }
    if (MBpars.LadderWeight()==ladder_weight::Regge && 
	dabs(m_kdiff)>1.e-6) {
      wt    *= diffwt =
	exp(-m_kdiff*sqrt(m_d2+sqr(log(Max(m_q01_2,m_Q02eff)/
				       Max(m_q12_2,m_Q02eff)))));
    }
    sup    = SuppressionTerm(m_q01_2,m_q12_2)*m_Q02eff/(m_Q02eff+kt12);
    wt    *= recombwt= 
      Min(1.,p_eikonal->EmissionWeight(m_b1,m_b2,dir?y1:-y1,sup));
    m_histomap[std::string("RecombWt")]->Insert(recombwt);
    if (dir) {
      m_histomap[std::string("ReggeWt1")]->Insert(reggewt);
      m_histomap[std::string("RecombWt1")]->Insert(recombwt);
      m_histomap[std::string("RecombSup1")]->Insert(sup);
    }
    else {
      m_histomap[std::string("ReggeWt0")]->Insert(reggewt);
      m_histomap[std::string("RecombWt0")]->Insert(recombwt);
      m_histomap[std::string("RecombSup0")]->Insert(sup);
    }
  } while (wt<ran->Get());

  m_k0 = k_0;
  m_k1 = k_1;
  m_k2 = k_2;

  if (MBpars.LadderWeight()==ladder_weight::Regge) {
    rarg = Min(mu12_2/q12.PPerp2(),q12.PPerp2()/mu12_2);
    expo = colfac*(*p_alphaS)(q12.PPerp2())*dabs(k_2.Y()-k_1.Y())/M_PI; 
//     rarg = Min(mu12_2/q12.Abs2(),q12.Abs2()/mu12_2);
//     expo = colfac*(*p_alphaS)(q12.Abs2())*dabs(k_2.Y()-k_1.Y())/M_PI; 
    m_lastwt  = pow(rarg,expo);
    //m_lastwt *= kmrwt = q12.PPerp2()/m_q12_2;
    //m_histomap[std::string("KMRWt")]->Insert(kmrwt);
  }

  m_histomap[std::string("Delta_final")]->Insert(1./(dabs(y1-y0))); 
  return true;
}

bool Final_State::MomViolation(ATOOLS::Vec4D & k_0, ATOOLS::Vec4D & k_1,
			       ATOOLS::Vec4D & k_2, bool dir) {
  //if (k_0[0]<0. || k_1[0]<0. || k_2[0]<0.) {
  //  m_rej_offshell++;
  //  return true;
  //}
  if (k_0[0]<0. || k_0.Nan() || dabs(k_0.Abs2())>1.e-3) {
    m_rej_offshell++;
    return true;
  }
  if (k_1[0]<0. || k_1.Nan() || dabs(k_1.Abs2())>1.e-3) {
    m_rej_offshell++;
    return true;
  }
  if (k_2[0]<0. || k_2.Nan() || dabs(k_2.Abs2())>1.e-3) {
    m_rej_offshell++;
    return true;
  }
  if (dabs(k_0.Abs2()>1.e-12)) k_0[0] = k_0.P();
  if (dabs(k_1.Abs2()>1.e-12)) k_1[0] = k_1.P();
  if (dabs(k_2.Abs2()>1.e-12)) k_2[0] = k_2.P();

  bool momviolation(false);
  Vec4D bef(m_k0+m_k2),aft(k_0+k_1+k_2);
  if (!IsEqual(aft.Abs2()/aft[0],bef.Abs2()/bef[0],1.e-6)) {
    momviolation = true;
  }
  if (dabs((bef-aft)[3]/bef[0])>1.e-6) {
    momviolation = true;
  }
//   k_2 = m_k0+m_k2-k_0-k_1;
  if (momviolation) m_rej_nofit++;
  return momviolation;
}

bool Final_State::
SelectSplitterAndSpectator(LadderMap::iterator & split,
			   LadderMap::iterator & spect) {
  //double pp2 = dabs(m_pprop.Abs2());
  //double pm2 = dabs(m_mprop.Abs2());
  //double yp  = dabs(m_plusiter->first);
  //double ym  = dabs(m_minusiter->first);
  bool   dir = ran->Get()>0.5; //(yp>ym)?true:false;
  if (dir) {
    split = m_plusiter;
    spect = m_minusiter;
    m_q0  = m_pprop;
    m_q2  = m_mprop;
  }
  else {
    split = m_minusiter;
    spect = m_plusiter;
    m_q0  = m_mprop;
    m_q2  = m_pprop;
  }
  return dir;
}  

void Final_State::
SwapSplitterAndSpectator(LadderMap::iterator & split,
			 LadderMap::iterator & spect) {
  if (split==m_minusiter) {
    split = m_plusiter;
    spect = m_minusiter;
    m_q0  = m_pprop;
    m_q2  = m_mprop;
  }
  else {
    split = m_minusiter;
    spect = m_plusiter;
    m_q0  = m_mprop;
    m_q2  = m_pprop;
  }
}

bool Final_State::FixPropColours(const LadderMap::iterator & split,
				 const LadderMap::iterator & spect) {
  bool prev(false),next(false),dir(true);
  TPropList::iterator pend(p_props->end());pend--;
  LadderMap::iterator emend(p_emissions->end());emend--;

  if (split==m_plusiter && spect==m_minusiter) dir = true;
  else if (split==m_minusiter && spect==m_plusiter) dir = false;
  else {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Do not understand orientation, will abort."<<std::endl;
    exit(1);
  }
    
  if (m_propiter!=p_props->begin()) {
    m_propiter--;
    if (m_propiter->m_col==colour_type::singlet){
      if(dir) prev = true;
      else next = true;
    }
    m_propiter++;
  }
  if (m_propiter!=pend) {
    m_propiter++;
    if (m_propiter->m_col==colour_type::singlet){
      if (dir) next = true;
      else prev = true;
    }  
    m_propiter--;
  }

  double y0(m_k0.Y()), y1(m_k1.Y()), y2(m_k2.Y());
  double wt81(0.), wt18(0.), wt88(0.);
  int    beam1 = int((split==p_emissions->begin()) || split==emend);
  int	 beam2 = int((spect==p_emissions->begin()) || spect==emend);

  double sup01(SuppressionTerm(m_q01_2,m_q01_2));
  double sup12(SuppressionTerm(m_q12_2,m_q12_2));

  double tot = wt18 = prev?0.:
    p_eikonal->SingletWeight(m_b1,m_b2,y0,y1,sup01,beam1)*
    p_eikonal->OctetWeight(m_b1,m_b2,y1,y2,sup12,beam2);
  tot +=wt81 = next?0.:
    p_eikonal->OctetWeight(m_b1,m_b2,y0,y1,sup01,beam1)*
    p_eikonal->SingletWeight(m_b1,m_b2,y1,y2,sup12,beam2);
  tot += wt88 =       
    p_eikonal->OctetWeight(m_b1,m_b2,y0,y1,sup01,beam1)*
    p_eikonal->OctetWeight(m_b1,m_b2,y1,y2,sup12,beam2);
        
  std::pair<colour_type::code,colour_type::code> cols;

  tot *= ran->Get();
  tot -= wt18;
  if (tot<0.) {
    cols.first  = colour_type::singlet;
    cols.second = colour_type::octet;
    p_ladder->SetDiffractive(true);
  }
  else {
    tot -= wt81;
    if (tot<0.) { 
      cols.first  = colour_type::octet;
      cols.second = colour_type::singlet;
      p_ladder->SetDiffractive(true);
    }
    else {
      cols.first  = colour_type::octet;
      cols.second = colour_type::octet;
    }
  }
  if (dir) {
    m_plusiter++;
    m_pprop          += m_k0;
    m_propiter->m_q   = m_q0+m_k0;
    m_propiter->m_q2  = m_propiter->m_q.Abs2();
    m_propiter->m_qt2 = m_propiter->m_q.PPerp2();
    m_propiter->m_q02 = Q02((y0+y1)/2.);
    m_propiter->m_col = cols.first;
    m_propiter++; 
    T_Prop prop(cols.second,-m_k2+m_q2,Q02((y1+y2)/2.));
    m_propiter        = p_props->insert(m_propiter,prop);
    if (m_propiter->m_col!=colour_type::singlet) return true;
  }
  else {
    m_minusiter--;
    m_mprop          -= m_k0;
    m_propiter->m_q   = -m_k0+m_q0;
    m_propiter->m_q2  = m_propiter->m_q.Abs2();
    m_propiter->m_qt2 = m_propiter->m_q.PPerp2();
    m_propiter->m_q02 = Q02((y0+y1)/2.);
    m_propiter->m_col = cols.first;
    T_Prop prop(cols.second,m_q2+m_k2,Q02((y1+y2)/2.));
    m_propiter        = p_props->insert(m_propiter,prop);
    if (m_propiter->m_col!=colour_type::singlet) return true;
  }
  return false;
}



bool Final_State::
IsOrdered(const bool & dir,Vec4D & k_0,
	  Vec4D & k_1,Vec4D & k_2,const double & t) {
  double xi0=k_0.PPerp2()/sqr(k_0[3]);
  double xi1=k_1.PPerp2()/sqr(k_1[3]);
  double xi2=k_2.PPerp2()/sqr(k_2[3]);
  double s  =(k_0+k_1).Abs2();
  double u  = s-t;

  double y0(k_0.Y()), y1(k_1.Y()), y2(k_2.Y());
  
  bool ordered(true);
  if (s<t || u<t) 
    msg_Tracking()<<"   --> s = "<<s<<", t = "<<t<<", u = "<<u<<std::endl
		  <<" from "<<k_0<<" "<<k_1<<std::endl
		  <<" with spectator "<<k_2<<", "<<std::endl
		  <<" and "<<m_q0<<" "<<m_q2<<"."<<std::endl; 
  switch (m_ordering) {
  case ordering::ao_phys:
    // physics s>t
    if (s<t) {
      ordered = false;
      break;
    }
    // plus
  case ordering::ao:
    // ao + rap ordering
    if (xi0>xi1 || xi2>xi1 || y0>y1 || y1>y2) ordered = false;
    break;
  case ordering::ao_keep:
    if (y0>y1) {
      Vec4D swap(k_0);
      k_0 = k_1;
      k_1 = swap;
    }
    if (xi0>xi1 || xi2>xi1) {
      ordered = false;
      break;
    }
    break;
  case ordering::keep:
    // allow swap for split <--> emit : reorder raps to enforce y0<y1
    // but insist on y1<y2
    if (y0>y1) {
      Vec4D swap(k_0);
      k_0 = k_1;
      k_1 = swap;
    }
    if (y1<y2 || y1>y2) {
      ordered = false;
      break;
    }
    break;
  case ordering::rap_phys:
    // physics s>t
    if (s<t) {
      ordered = false;
      break;
    }
    // plus 
  case ordering::rap_only:
  default:
    // rapidity ordering
    if (y0>y1 || y1>y2) ordered = false;
    break;
  }
  if (!ordered) {
    if (dir) m_rej_order1++; else m_rej_order0++;
    m_rej_order++;
  }
  return ordered;
}

double Final_State::SuppressionTerm(const double & q02,const double & q12) {
  return sqrt(p_alphaS->Weight(q02,true)*p_alphaS->Weight(q12,true));
}

double Final_State::Q02(const double & y,const bool & first) {
  return m_Q02eff; //+m_Q12*term;
  //double term = pow(p_eikonal->Sum(m_b1,m_b2,y),1.5);
}

double Final_State::SelectKT2(const double & kt2max,const double & kt2min,
			      const double & mu2,const double & ktexpo) {
  double fixterm(0.), logterm(0.), expo(1.-ktexpo), kt2(-1.);
  double kt2cut(Max(mu2,kt2min)),rand(ran->Get());
  if (kt2max>kt2min) {
    switch (m_ktform) {
    case ktform::cut:
      if (kt2max>mu2) {
	if (expo==0.) 
	  kt2 = kt2cut * pow(kt2max/kt2cut,rand);
	else 
	  kt2 = pow(rand*pow(kt2max,expo)+(1.-rand)*pow(kt2cut,expo),1./expo);
      }
      break;      
    case ktform::smooth:
      if (expo==0.) 
	kt2 = (kt2min+mu2)*pow((kt2max+mu2)/(kt2min+mu2),rand)-mu2; 
      else 
	kt2 = pow(rand*pow(kt2max+mu2,expo)+
		  (1.-rand)*pow(kt2cut+mu2,expo),1./expo)-mu2;
      break;      
    case ktform::IR0:
      if (kt2min<mu2) fixterm = 0.5*(1.-sqr(kt2min/mu2));
      if (expo==0.) logterm = log(Max(mu2,kt2max)/kt2cut);
      else logterm = (pow(kt2cut,expo)-pow(Max(mu2,kt2max),expo))/expo;
      if (fixterm>ran->Get()*(fixterm+logterm)) {
	kt2 = sqrt(kt2min+rand*(mu2-kt2min));
      }
      else {
	if (expo==0.) 
	  kt2 = kt2cut * pow(Max(kt2max,mu2)/kt2cut,rand);
	else 
	  kt2 = pow(rand*pow(kt2max,expo)+(1.-rand)*pow(kt2cut,expo),1./expo);
      }
      break;
    case ktform::frozen:
    default:
      if (kt2min<mu2) fixterm = 1.-kt2min/mu2;
      if (expo==0.) logterm = log(Max(mu2,kt2max)/kt2cut);
      else logterm = (pow(kt2cut,expo)-pow(Max(mu2,kt2max),expo))/expo;
      if (fixterm>ran->Get()*(fixterm+logterm)) {
	kt2 = kt2min+rand*(mu2-kt2min);
      }
      else {
	if (expo==0.) 
	  kt2 = kt2cut * pow(Max(kt2max,mu2)/kt2cut,rand);
	else 
	  kt2 = pow(rand*pow(kt2max,expo)+(1.-rand)*pow(kt2cut,expo),1./expo);
      }
      break;
    }
  }
  if (kt2<0.) {
  }
  return kt2;
}

double Final_State::KT2integral(const double & kt2max,const double & kt2min,
				const double & q02,const double & ktexpo) {
  double fixterm(0.), logterm(0.), expo(1.-ktexpo);
  if (kt2max>kt2min) {
    switch (m_ktform) {
    case ktform::cut:
      if (kt2max>q02) {
	if (expo==0.) logterm = log(kt2max/Max(q02,kt2min));
	else logterm = (pow(Max(q02,kt2min),expo)-pow(kt2max,expo))/expo;
      }
      break;      
    case ktform::smooth:
      if (expo==0.) logterm = log((kt2max+q02)/(kt2min+q02));
      else logterm = (pow(kt2max/q02+1.,expo)-pow(kt2min/q02+1.,expo))/expo;
      break;      
    case ktform::IR0:
      if (kt2min<q02) fixterm = 0.5*(1.-sqr(kt2min/q02));
      if (expo==0.) logterm = log(Max(q02,kt2max)/Max(q02,kt2min));
      else logterm = (pow(Max(q02,kt2min),expo)-pow(Max(q02,kt2max),expo))/expo;
      break;
    case ktform::frozen:
    default:
      if (kt2min<q02) fixterm = 1.-kt2min/q02;
      if (expo==0.) logterm = log(Max(q02,kt2max)/Max(q02,kt2min));
      else logterm = (pow(Max(q02,kt2min),expo)-pow(Max(q02,kt2max),expo))/expo;
      break;
    }
  }
  return fixterm + logterm;
}

void Final_State::InitHistograms() {
  m_histomap[std::string("Deltay_test")] = 
    new Histogram(0,0.0,10.0,100);
  m_histomap[std::string("KT2_test1")] = 
    new Histogram(0,0.0,100.0,100);
  m_histomap[std::string("KT2_accept1")] = 
    new Histogram(0,0.0,100.0,100);
  m_histomap[std::string("KT2_test2")] = 
    new Histogram(0,0.0,4.0,100);
  m_histomap[std::string("KT2_accept2")] = 
    new Histogram(0,0.0,4.0,100);
  m_histomap[std::string("Delta_naive")] = 
    new Histogram(0,0.0,20.0,100);
  m_histomap[std::string("Delta_naive0")] = 
    new Histogram(0,0.0,20.0,100);
  m_histomap[std::string("Delta_naive1")] = 
    new Histogram(0,0.0,20.0,100);
  m_histomap[std::string("Delta_kin")] = 
    new Histogram(0,0.0,5.0,100);
  m_histomap[std::string("Delta_order")] = 
    new Histogram(0,0.0,5.0,100);
  m_histomap[std::string("Delta_final")] = 
    new Histogram(0,0.0,5.0,100);
  //m_histomap[std::string("KMRWt")] = 
  //  new Histogram(0,0.0,20.0,200);
  m_histomap[std::string("ReggeWt")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("ReggeWt0")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("ReggeWt1")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("RecombWt")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("RecombWt0")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("RecombWt1")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("RecombSup0")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("RecombSup1")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("FSWt")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("FSWt0")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("FSWt1")] = 
    new Histogram(0,0.0,2.0,200);
  m_histomap[std::string("ytest1")] = 
    new Histogram(0,-10.0,10.0,100);
  m_histomap[std::string("ytest0")] = 
    new Histogram(0,-10.0,10.0,100);
  m_histomap[std::string("q01_2_1")] = 
    new Histogram(0,0.0,400.0,100);
  m_histomap[std::string("q01_2_0")] = 
    new Histogram(0,0.0,400.0,100);
  m_histomap[std::string("q12_2_1")] = 
    new Histogram(0,0.0,400.0,100);
  m_histomap[std::string("q12_2_0")] = 
    new Histogram(0,0.0,400.0,100);
  m_histomap[std::string("Deltay_regexit")] = 
    new Histogram(0,0.0,20.0,100);
  m_histomap[std::string("Deltay_singexit")] = 
    new Histogram(0,0.0,20.0,100);

  m_trials=m_rej_negkt1=m_rej_nokt2=m_rej_nokt0=m_rej_nophi0=m_rej_noy0=0;
  m_rej_offshell=m_rej_nofit=m_rej_order=m_alphaS=m_ys=0;
  m_dir0=m_dir1=m_rej_order1=m_rej_order0=m_cols=0;
  m_singexit=m_firstsing=0;
}

void Final_State::OutputHistograms() {
  msg_Info()
    <<METHOD<<" yields rejections from "<<m_trials<<" trial emissions:\n"
    <<"   regular exits       "<<m_ys<<"\n"
    <<"   alphaS weight       "<<m_alphaS<<"\n"
    <<"   negative kt1        "<<m_rej_negkt1<<"\n"
    <<"   no kt2 constructed  "<<m_rej_nokt2<<"\n"
    <<"   no kt0 constructed  "<<m_rej_nokt0<<"\n"
    <<"   no phi0 found       "<<m_rej_nophi0<<"\n"
    <<"   no y0 found         "<<m_rej_noy0<<"\n"
    <<"   offshell momenta    "<<m_rej_offshell<<"\n"
    <<"   mismatch kinematics "<<m_rej_nofit<<"\n"
    <<"   unordered emissions "<<m_rej_order<<".\n";
  msg_Info()<<"   A posteriori colour reorderings: "<<m_cols<<".\n";
  msg_Info()<<"   Ratio of forward/backward trials = "
	    <<double(m_dir1)/double(m_dir0)<<", "
	    <<" in ordering rejections = "
	    <<double(m_rej_order1)/double(m_rej_order)<<".\n";
  msg_Info()
    <<"   Naive Delta from Regge approximation    = "
    <<m_histomap[std::string("Delta_naive")]->Average()<<", "
    <<"asymmetry = "
    <<(m_histomap[std::string("Delta_naive1")]->Average()/
       m_histomap[std::string("Delta_naive0")]->Average())<<";\n"
    <<"   after taking only kinematically allowed = "
    <<m_histomap[std::string("Delta_kin")]->Average()<<";\n"
    <<"   after taking only ordered emissions     = "
    <<m_histomap[std::string("Delta_order")]->Average()<<";\n"
    <<"   after all weights, final Delta          = "
    <<m_histomap[std::string("Delta_final")]->Average()<<".\n";
  msg_Info()
    <<METHOD<<" ways of exiting:\n"
    <<"   regular exits       "<<m_ys<<"\n"
    <<"   singlet exits       "<<m_singexit<<"\n"
    <<"   first singlet       "<<m_firstsing<<"\n";
  
  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator 
	 hit=m_histomap.begin();hit!=m_histomap.end();hit++) {
    histo = hit->second;
    name  = std::string("Ladder_Analysis/")+hit->first+std::string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histomap.clear();
}
