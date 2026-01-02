#include "SHRiMPS/Ladders/Ladder_Generator_TMD.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Ladder_Generator_TMD::Ladder_Generator_TMD() :
  Ladder_Generator_Base(),
  m_S(sqr(rpa->gen.Ecms())), m_Ecms(rpa->gen.Ecms()/2.), m_pt02(0.25),
  m_Pplus(Vec4D(1.,0.,0.,1.)), m_Pminus(Vec4D(1.,0.,0.,-1.)),
  m_sigmaTMD(Sigma_TMD(m_pt02)),
  m_test(true)
{
  m_sigmaTMD.SetAlphaS(p_alphaS);
  if (m_test) {
    m_histos[std::string("Y_first")]        = new Histogram(0, -8.0,  8.0, 32);
    m_histos[std::string("Y_first_pt")]     = new Histogram(0, -8.0,  8.0, 32);
    m_histos[std::string("Y_first_highpt")] = new Histogram(0, -8.0,  8.0, 32);
  }
}

Ladder_Generator_TMD::~Ladder_Generator_TMD() {
  for (size_t beam=0;beam<2;beam++) { delete p_tmd[beam]; p_tmd[beam]=NULL; }
  if (m_test) {
    std::string name  = std::string("Ladder_Analysis/");
    for (std::map<std::string, Histogram * >::iterator hit=m_histos.begin();
	 hit!=m_histos.end();hit++) {
      hit->second->Finalize();
      hit->second->Output(name+hit->first);
      delete hit->second;
    }
  }
}

void Ladder_Generator_TMD::Initialise(PDF::ISR_Handler *const isr) {
  for (size_t beam=0;beam<2;beam++) {
    double Y = (beam==0?-1.:1.)*m_Ymax;
    p_tmd[beam] = new REMNANTS::Pseudo_TMD(isr->PDF(beam),Y,m_pt02);
  }
  m_sigmaTMD.Initialise(p_tmd);
}

Ladder * Ladder_Generator_TMD::operator()(const Vec4D & pos) {
  msg_Out()<<"     --- "<<METHOD<<" @ position = "<<pos<<"\n";
  InitialiseLadder(pos);
  if (m_ngluons>0) AddEmissions();
  AddInitialState();
  AddPropagators();
  SelectPropagatorColours();
  msg_Out()<<"\n"<<(*p_ladder)<<"\n";
  return p_ladder;
}

void Ladder_Generator_TMD::InitialiseLadder(const Vec4D & pos) {
  m_ylimits[1] = -m_Ymax;
  m_ylimits[0] =  m_Ymax;
  p_ladder     = new Ladder(pos);
  p_emissions  = p_ladder->GetEmissions();
  p_props      = p_ladder->GetProps();
  m_ngluons    = m_density.NGluons(-m_Ymax,m_Ymax);
  msg_Out()<<"     --- N_gluons =  "<<m_ngluons<<" in "
  	   <<"rapidity range = ["<<m_ylimits[1]<<", "<<m_ylimits[0]<<"], "
  	   <<"energy range ["<<m_E[0]<<", "<<m_E[1]<<"].\n";
  if (m_ngluons==0.) {
    FillZeroEmissionLadder();
    return;
  }
  double xsec  = m_sigmaTMD.MakeEvent();
  Vec4D  kvec  = m_sigmaTMD.GetEmit();
  double kt    = kvec.PPerp();
  m_Yhat       = kvec.Y();  
  m_ktsum      = kvec;
  p_ladder->AddRapidity(m_Yhat,Flavour(kf_gluon),kvec);
  for (size_t beam=0;beam<2;beam++) {
    m_qtprev[beam]    = m_sigmaTMD.GetProp(beam);
    m_qt2prev[beam]   = m_sigmaTMD.GetQT2(beam);
    m_xprev[beam]     = m_sigmaTMD.GetX(beam);
    m_xtmdprev[beam]  = m_sigmaTMD.xTMD(beam);
    m_thetaprev[beam] = kvec.Theta();
    for (size_t j=0;j<2;j++) m_y[beam][j] = m_Yhat;
  }
  if (m_test) FillAnalysis();
  msg_Out()<<"     --- start with emission: kt = "<<kt<<" @ y = "<<m_Yhat<<", "
  	   <<"theta = "<<(180.*kvec.Theta()/M_PI)<<"\n"
  	   <<"         from "<<m_sigmaTMD.GetProp(0)<<" + "<<m_sigmaTMD.GetProp(1)<<" "
  	   <<"--> "<<m_sigmaTMD.GetEmit()<<"\n";
}

void Ladder_Generator_TMD::FillZeroEmissionLadder() {
  double y[2], x[2], ypow(2.), qt, qt0(sqrt(m_pt02)), weight;
  double phi   = 2.*M_PI*ran->Get();
  do { qt  = qt0*dabs(ran->GetGaussian()); }
  while (2.*qt*cosh(Max(dabs(y[0]),dabs(y[1])))>Min(m_E[0],m_E[1]));
  m_qt2prev[0] = m_qt2prev[1] = sqr(qt);
  m_qtprev[0]  = qt*Vec4D(0.,cos(phi),sin(phi),0.);
  m_qtprev[1]  = -m_qtprev[0];
  m_ktsum      = m_qtprev[0];
}

void Ladder_Generator_TMD::AddEmissions() {
  msg_Out()<<"     ----- "<<METHOD<<" adds emissions starting from "<<m_y[0][0]<<" & "<<m_y[1][0]
  	   <<" in ["<<m_ylimits[0]<<", "<<m_ylimits[1]<<"]\n";
  double deltay[2];
  do {
    for (size_t j=0;j<2;j++) deltay[j] = dabs(m_ylimits[j]-m_y[j][0]);
    if (deltay[0]<m_deltaY && deltay[1]<m_deltaY) break;
    size_t beam = (deltay[0]>deltay[1] ? 0 : 1);
    msg_Out()<<"\n     ----- compare y = "<<m_y[0][0]<<", "<<m_y[1][0]<<" "
    	     <<"("<<deltay[0]<<" vs. "<<deltay[1]<<") --> beam = "<<beam<<"\n"; 
    if (TrialEmission(beam)) {
      m_ktsum         += m_ktvectest;
      msg_Out()<<"     ----- accepted: y = "<<m_y[beam][1]<<", "
	       <<"qt = "<<m_qtvectest<<" - qt_prev = "<<m_qtprev[beam]<<"\n"
	       <<"     -----           kt = "<<m_ktvectest<<" "<<"(kt = "<<m_kttest<<", "
	       <<"all kt = "<<m_ktsum<<")\n";
      m_y[beam][0]     = m_y[beam][1];
      m_xprev[beam]    = m_xtest;
      m_qtprev[beam]   = m_qtvectest;
      m_qt2prev[beam]  = m_qt2test;
      m_xtmdprev[beam] = m_xtmdtest;
      p_ladder->AddRapidity(m_y[beam][0],Flavour(kf_gluon),m_ktvectest);
    }
    else {
      m_y[beam][0] = m_ylimits[beam];
    }
  } while (m_ylimits[0]>m_y[0][1] || m_y[1][1]>m_ylimits[1]); 
}

bool Ladder_Generator_TMD::TrialEmission(const size_t beam) {
  msg_Out()<<"       --- "<<METHOD<<"(beam = "<<beam<<"): "
  	   <<m_y[beam][0]<<" --> "<<m_y[beam][1]<<"  ("<<m_ylimits[beam]<<")\n";
  double qt2ratio = log(m_S/(4.*m_pt02)+1.);
  double arg      = M_PI/(3.*AlphaSMax()*qt2ratio);
  double deltay, phitest;
  do {
    deltay         = arg*log(ran->Get());
    m_y[beam][1]  += (beam==0 ? -deltay : +deltay);
    m_qt2test      = m_pt02 *  (pow(qt2ratio,ran->Get())+1.);
    msg_Out()<<"\n       --- next trial at deltay = "<<deltay<<" --> y_test = "<<m_y[beam][1]<<", "
    	     <<"qt2_test = "<<m_qt2test<<".\n";
    // Check if trial emission still within rapidity range.
    if (beam==1 ? m_y[beam][1]<m_ylimits[beam] : m_y[beam][1]>m_ylimits[beam]) return false;
    // Reconstruct kinematics of emitted parton
    phitest        = 2.*M_PI*ran->Get();
    m_qtvectest    = sqrt(m_qt2test)*Vec4D(0.,cos(phitest),sin(phitest),0.);
    m_ktvectest    = m_qtvectest-m_qtprev[beam];
    m_kttest       = m_ktvectest.PPerp();
    m_ktvectest[0] = m_kttest*cosh(m_y[beam][1]);
    m_ktvectest[3] = m_kttest*sinh(m_y[beam][1]);
    m_qtvectest    = m_qtprev[beam]+m_ktvectest;
  } while (EmissionWeight(beam)<ran->Get());
  return true;
}

double Ladder_Generator_TMD::EmissionWeight(const size_t beam) {
  m_xtest    = (beam==0 ?
		(m_qtvectest[0]+m_qtvectest[3]) :
		(m_qtvectest[0]-m_qtvectest[3]))/(2.*m_Ecms);
  if (m_xtest<1.e-5 || m_xtest>0.95) {
    //msg_Out()<<"         -   "<<METHOD<<"(x = "<<m_xtest<<")\n";
    return 0.;
  }
  //if (sqr(m_kttest)<m_qt2test) return 0.;
  p_tmd[beam]->Calculate(m_xtest,sqr(m_kttest),m_qt2test,m_y[beam][1]);
  m_xtmdtest = p_tmd[beam]->XTMD(Flavour(kf_gluon));
  double tmdweight = ( 1. * /*m_xprev[beam]/m_xtest */
		       m_xtmdtest/m_xtmdprev[beam]);
  double absweight = AbsorptionWeight(beam);
  double weight    = AlphaSWeight(sqr(m_kttest)) * tmdweight; 
  msg_Out()<<"         -   "<<METHOD
  	   <<"(x = "<<m_xtest<<" vs. "<<m_xprev[beam]<<"): weight = "<<weight<<"\n"
  	   <<"         -   (as = "<<AlphaSWeight(sqr(m_kttest))<<", "
	   <<"tmd = "<<tmdweight<<", abs = "<<absweight<<")\n";
  return weight;
}

double Ladder_Generator_TMD::AbsorptionWeight(const size_t beam) {
  return m_density.AbsorptionWeight(m_y[beam][1]);
}

void Ladder_Generator_TMD::AddInitialState() {
  //msg_Out()<<"     - "<<METHOD<<"(ktsum = "<<m_ktsum<<",\n"
  //	   <<"     -             qt_0 = "<<m_qtprev[0]<<", "<<m_qtprev[1]<<")\n";
  double kt2, kt, arg, dy, y;
  Vec4D  out[2];
  double alpha = 0., beta = 0., PpPm=m_Pplus*m_Pminus;
  for (size_t beam=0;beam<2;beam++) {
    alpha += m_qtprev[beam]*m_Pminus/PpPm;
    beta  += m_qtprev[beam]*m_Pplus/PpPm;
  }
  for (size_t beam=0;beam<2;beam++) {
    kt2 = dabs(m_qtprev[beam].PPerp2());
    kt  = sqrt(kt2);
    arg = M_PI/(3.*AlphaSMax()*log(1.+kt*m_E[beam]/m_pt02));
    do {
      dy = -log(ran->Get())*arg;
      y  = m_ylimits[beam] + (beam==0?-1.:1.) * dy;
    } while (2.*kt*cosh(y)>m_E[beam]);
    out[beam] = ( Vec4D(kt*cosh(y), 0., 0., kt*sinh(y)) -
		  m_qtprev[beam].Perp() );
    alpha    += out[beam]*m_Pminus/PpPm;
    beta     += out[beam]*m_Pplus/PpPm;
    m_ktsum  += out[beam];
    //msg_Out()<<"     - out["<<beam<<"] = "<<out[beam]<<" @ "<<out[beam].Y()
    //	     <<" ("<<out[beam].Abs2()<<")\n";
    p_ladder->AddRapidity(y,Flavour(kf_gluon),out[beam]);
  }
  p_ladder->InPart(0)->SetMomentum(alpha*m_Pplus);
  p_ladder->InPart(1)->SetMomentum(beta*m_Pminus);
  //msg_Out()<<"kt sum = "<<m_ktsum<<" vs "
  //	   <<(p_ladder->InPart(0)->Momentum()+p_ladder->InPart(1)->Momentum())<<", "
  //	   <<"in^ = "<<p_ladder->InPart(0)->Momentum().Abs2()<<" & "
  //	   <<p_ladder->InPart(1)->Momentum().Abs2()<<"\n"
  //	   <<(*p_ladder)<<"\n";
  //exit(1);
}

void Ladder_Generator_TMD::AddPropagators() {
  Vec4D qtvec = p_ladder->InPart(0)->Momentum();
  LadderMap::iterator lit=p_emissions->begin(),lend=p_emissions->end();lend--;
  while (lit!=lend) {
    qtvec -= lit->second.Momentum(); 
    p_props->push_back(T_Prop(colour_type::octet,
			      (qtvec[0]<0?-1.:1.)*qtvec,m_qt2min));
    lit++;
  }
  //msg_Out()<<"     - final check: "<<qtvec<<" vs. "<<m_qtprev[1]<<".\n";
}

void Ladder_Generator_TMD::SelectPropagatorColours() {
  LadderMap::iterator lit1=p_emissions->begin(),  lit2=p_emissions->end();  lit2--;
  TPropList::iterator pit1=p_props->begin(), pit2=p_props->end(); pit2--;
  double y1, y2, wt1, wt8;
  size_t dir;
  while (lit1->first>lit2->first) {
    dir = (dabs(lit1->first) > dabs(lit2->first))?1:0;
    if (dir) { y1 = lit1->first; lit1++; y2 = lit1->first; }
    else     { y2 = lit2->first; lit2--; y1 = lit2->first; }
    wt1 = m_density.SingletWeight(y1,y2);
    wt8 = m_density.OctetWeight(y1,y2);
    if (wt1/(wt1+wt8)>ran->Get()) {
      if (dir) { pit1->SetCol(colour_type::singlet); pit1++; }
      else     { pit2->SetCol(colour_type::singlet); pit2--; }
    }
  }
  pit1 = p_props->begin(); pit2 = pit1; pit2++;
  while (pit2!=p_props->end()) {
    if (pit1->Col()==colour_type::singlet && pit2->Col()==colour_type::singlet) {
      if (ran->Get()>0.5) pit1->SetCol(colour_type::octet); 
      else pit2->SetCol(colour_type::octet);
    }
    pit1++; pit2++;
  }
}

void Ladder_Generator_TMD::FillAnalysis() {
  m_histos[std::string("Y_first")]->Insert(m_Yhat,1.);
  m_histos[std::string("Y_first_pt")]->Insert(m_Yhat,m_ktsum.PPerp(),1.);
  if (m_ktsum.PPerp()>2.5)
    m_histos[std::string("Y_first_highpt")]->Insert(m_Yhat,1.);
}
