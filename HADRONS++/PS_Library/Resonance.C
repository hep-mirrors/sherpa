#include "HADRONS++/PS_Library/Resonance.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


using namespace HADRONS;
using namespace ATOOLS;
using namespace std;


Resonance_Base::Resonance_Base(const Res_Params & params) :
  m_inflav(params.m_inflav), m_type(params.m_type),
  m_OSmass(params.m_OSmass<0.   ? m_inflav.Mass(true) : params.m_OSmass),
  m_OSmass2(m_OSmass*m_OSmass),
  m_OSwidth(params.m_OSwidth<0. ? m_inflav.Width() : params.m_OSwidth),
  m_OSwidth2(m_OSwidth*m_OSwidth),
  m_weight(params.m_weight), m_phase(params.m_phase),
  m_threshold(0.), m_threshold2(0.), m_lambda(0.),
  m_exponent(m_inflav.IsVector() ? 3 : (m_inflav.IsScalar() ? 1 : 0)) 
{
  m_name = std::string("R_")+m_inflav.IDName();
  if (m_type==resonance_type::bespoke)      m_name += std::string("bespoke");
  else if (m_type==resonance_type::running) m_name += std::string("running");
  else                                      m_name += std::string("fixed");
  m_decmasses2.resize(params.m_outflavs.size(),0.);
  for (size_t i=0;i<params.m_outflavs.size();i++) {
    double mass     = params.m_outflavs[i].Mass(true);
    m_threshold    += mass;
    m_decmasses2[i] = ATOOLS::sqr(mass);
  }
  m_threshold2  = ATOOLS::sqr(m_threshold);
  m_lambda      = Lambda(m_OSmass2,m_decmasses2[0],m_decmasses2[1]);
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

Complex Resonance_Tree::operator()(const double & s,const double & s1) {
  msg_Out()<<"--------------------------------------------------\n"
	   <<METHOD<<"("<<sqrt(s)<<", "<<sqrt(s1)<<")\n";
  return ( m_norm * ( p_start ? p_start->BreitWigner(s) : (1.,0.) ) *
	   p_nodes->BreitWigner(s1) );
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

RunningWidth_Resonance::RunningWidth_Resonance(const Res_Params & params) :
  Resonance_Base(params)
{}

double RunningWidth_Resonance::CalculateWidth(const double & s) {
  //msg_Out()<<METHOD<<": s = "<<sqrt(s)<<" vs. "<<sqrt(m_threshold2)
  //	   <<"[pass = "<<(s>m_threshold2)<<"]"
  //	   <<" -> "<<m_OSwidth<<"\n";
  if (s>m_threshold2) {
    if (dabs(1.-s/m_OSmass2)<1.e-12) return m_OSwidth;
    if (m_decmasses2.size()==2) {
      return ( m_OSwidth *
	       pow(Lambda(s,m_decmasses2[0],m_decmasses2[1])/m_lambda, m_exponent) *
	       sqrt(s)/m_OSmass
	       );
    }
  }
  return 0.;
}

Complex RunningWidth_Resonance::BreitWigner(const double & s) {
  double MW = sqrt(s)*CalculateWidth(s);
  //msg_Out()<<METHOD<<"("<<sqrt(s)<<"): "<<MW<<"\n";
  return m_OSmass2/Complex( m_OSmass2-s, -MW );
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


GS_Resonance::GS_Resonance(const Res_Params & params) : 
  Resonance_Base(params),
  m_m_pi(Flavour(kf_pi_plus).Mass(true))
{
  double mpi2 = sqr(m_m_pi);
  m_d   = ( 3/M_PI * mpi2/sqr(m_lambda) * log((m_OSmass+2.*m_lambda)/(2.*m_m_pi)) +
	    m_OSmass/(2.*M_PI*m_lambda) * (1. - 2.*mpi2/sqr(m_lambda)) );
  m_h0  = h(m_OSmass2);
  m_dh0 = m_h0 * (1./(8.*sqr(m_lambda))-1./(2.*m_OSmass2));
}

double GS_Resonance::CalculateWidth(const double & s) {
  if (s>m_threshold2) {
    if (dabs(1.-s/m_OSmass2)<1.e-12) return m_OSwidth;
    if (m_decmasses2.size()==2) {
      return ( m_OSwidth *
	       pow(Lambda(s,m_decmasses2[0],m_decmasses2[1])/sqrt(s) /
		   (m_lambda / m_OSmass), m_exponent)
	       );
    }
  }
  return m_OSwidth;
}

Complex GS_Resonance::BreitWigner(const double & s) {
  double Gamma = CalculateWidth(s);
  return ( (m_OSmass2+m_d*m_OSmass*m_OSwidth)/
	   Complex(m_OSmass2-s+f(s), -sqrt(s)*Gamma) );
}

double GS_Resonance::h(const double & s) {
  double ks = Lambda(s,m_decmasses2[0],m_decmasses2[1]), E = sqrt(s);
  return ( 2.*ks/(M_PI*E)*log((E+2.*ks)/(2.*m_m_pi)) );
}

double GS_Resonance::f(const double & s) { 
  return ( m_OSwidth*m_OSmass2/pow(m_lambda,3) *
	   ( sqr(Lambda(s,m_decmasses2[0],m_decmasses2[1])) * (h(s)-m_h0) +
	     (m_OSmass2-s) * sqr(m_lambda) * m_dh0 )  );
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

// Comment: we should generalise this to *** ALL*** members of the a1 multiplet,
// in particular K_1^+(1400), which decays predominantly into:
// K_1(1400) -> K^*(892) + pi (94%)
// K_1(1400) -> K rho(770)    (3%)

RunningWidth3_Resonance::RunningWidth3_Resonance(const Res_Params & params) :
  Resonance_Base(params), m_allm2(0.), m_nbins(100), m_nsbins(50)
{
  size_t rindex = 0;
  vector<Flavour> flavs; flavs.resize(2);
  Res_Params V;
  for (size_t i=0;i<3;i++) {
    m_allm2 += params.m_outflavs[i].Mass(true);
    for (size_t j=0;j<3;j++) {
      if (i==j) continue;
      if (params.m_outflavs[i].Charge()<=params.m_outflavs[j].Charge()) continue;
      flavs[0] = params.m_outflavs[i];
      flavs[1] = params.m_outflavs[j];
      if (GetParams(flavs,V)) {
	if (rindex==0) {
	  m_mi2       = m_decmasses2[i];
	  m_mj2       = m_decmasses2[j];
	  p_BW_ij     = new RunningWidth_Resonance(V);
	  m_delta2_ij = m_mi2-m_mj2;
	  m_sum2_ij   = m_mi2+m_mj2;
	  m_smin_ij   = sqr(flavs[0].Mass(true)+flavs[1].Mass(true));
	  //msg_Out()<<"Made new resonance ["<<i<<j<<"]: "<<p_BW_ij->Flav()<<", "
	  //	   <<"Delta = "<<m_delta2_ij<<", Sum = "<<m_sum2_ij<<", "
	  //	   <<"s > "<<m_smin_ij<<"\n";
	}
	else {
	  m_mk2       = m_decmasses2[i];
	  p_BW_kj     = new RunningWidth_Resonance(V);
	  m_delta2_kj = m_mk2-m_mj2;
	  m_sum2_kj   = m_mk2+m_mj2;
	  m_smin_kj   = sqr(flavs[0].Mass(true)+flavs[1].Mass(true));
	  m_smin_ik   = sqr(sqrt(m_mi2)+sqrt(m_mk2));
	  //msg_Out()<<"Made new resonance ["<<i<<j<<"]: "<<p_BW_kj->Flav()<<", "
	  //	   <<"Delta = "<<m_delta2_kj<<", Sum = "<<m_sum2_kj<<", "
	  //	   <<"s > "<<m_smin_kj<<"\n";
	}
	rindex++;
      }
    }
  }
  FillTable();
  m_gV = (*p_g)(sqr(m_inflav.Mass(true)));
}

RunningWidth3_Resonance::~RunningWidth3_Resonance() {
  if (p_BW_ij) { delete p_BW_ij; p_BW_ij = NULL; }
  if (p_BW_kj) { delete p_BW_kj; p_BW_kj = NULL; }
  if (p_g)     { delete p_g;     p_g     = NULL; }
}

bool RunningWidth3_Resonance::GetParams(vector<Flavour> & flavs,Res_Params & V) {
  if (flavs.size()!=2) return false;
  if ( (flavs[0].Kfcode()==kf_pi_plus && flavs[1].Kfcode()==kf_pi) ||
       (flavs[0].Kfcode()==kf_pi && flavs[1].Kfcode()==kf_pi_plus) ) {
    V.m_inflav   = Flavour(kf_rho_770_plus);
    V.m_outflavs = flavs;
    V.m_type     = resonance_type::running;
    V.m_OSmass   = 0.774;
    V.m_OSwidth  = 0.143;
    return true;
  }
  if ( (flavs[0]==Flavour(kf_pi_plus) && flavs[1]==Flavour(kf_pi_plus).Bar()) ||
       (flavs[0]==Flavour(kf_pi_plus).Bar() && flavs[1]==Flavour(kf_pi_plus)) ) {
    V.m_inflav   = Flavour(kf_rho_770);
    V.m_outflavs = flavs;
    V.m_type     = resonance_type::running;
    V.m_OSmass   = 0.774;
    V.m_OSwidth  = 0.143;
    return true;
  }
  if ( (flavs[0]==Flavour(kf_K_S) || flavs[0]==Flavour(kf_K_L) ||
	flavs[0]==Flavour(kf_K)) &&
       (flavs[1]==Flavour(kf_pi_plus) ||
	flavs[1]==Flavour(kf_pi_plus).Bar()) ) {
    V.m_inflav   = Flavour(kf_K_star_892_plus);
    V.m_outflavs = flavs;
    V.m_type     = resonance_type::running;
    V.m_OSmass   = 0.8955;
    V.m_OSwidth  = 0.051;
    return true;
  }
  return false;
}


double RunningWidth3_Resonance::
dg(const double & Q2,const double & sij,const double & skj) {
  // this is the implementation of the integrand of equations (3.8) and (3.9)
  // in the original Kuhn-Santamaria paper.  Note the two relative minus signs that
  // have been absorbed here.
  double  sik    = Q2-sij-skj+m_allm2;
  if (sik<=m_smin_ik || sik>=sqr(sqrt(Q2)-sqrt(m_mj2))) return 0.;
  double  Vij2  = ( ( sij-2.*m_sum2_ij + sqr(m_delta2_ij)/sij ) + 
		    sqr(sik-skj+m_delta2_ij-(Q2+sij-m_mk2)*sqr(m_delta2_ij)/sij)/(4.*Q2) );
  double  Vkj2  = ( ( skj-2.*m_sum2_kj + sqr(m_delta2_kj)/skj ) + 
		    sqr(sik-sij+m_delta2_kj-(Q2+skj-m_mi2)*sqr(m_delta2_kj)/skj)/(4.*Q2) );
  double  Vijkj = ( ( (Q2-2.*sik+m_allm2-3.*m_mj2) -
		      (Q2-m_allm2+m_mj2)*(m_delta2_ij*m_delta2_kj)/(sij*skj) +
		      (Q2-2.*skj-m_sum2_kj+m_mi2)*m_delta2_ij/sij +
		      (Q2-2.*sij-m_sum2_ij+m_mk2)*m_delta2_kj/skj )/2. +
		    ( (sik-skj+m_delta2_ij-(Q2-2.*skj-m_sum2_ij+m_mk2)*m_delta2_ij/sij) *
		      (sik-sij+m_delta2_kj-(Q2-2.*sij-m_sum2_kj+m_mi2)*m_delta2_kj/skj) )/
		    (4.*Q2) );
  Complex BWij = p_BW_ij->BreitWigner(sij), BWkj = p_BW_kj->BreitWigner(skj);
  double  val = (Vij2  * sqr(abs(BWij)) + Vkj2  * sqr(abs(BWkj)) +
		 2.*Vijkj * abs(BWij*conj(BWkj)));
  return val;
}

void RunningWidth3_Resonance::FillTable(const double & Qmax) {
  p_g = new OneDim_Table(axis(m_nbins,m_threshold,
			      sqr(Qmax<0.?Flavour(kf_tau).Mass(true):Qmax)));
  for (size_t q2=0;q2<m_nbins;q2++) {
    double Q2      = p_g->GetAxis().x(q2);
    double value   = 0.;
    double sijmax  = sqr(sqrt(Q2)-sqrt(m_mk2));
    double sijstep = (sijmax-m_smin_ij)/double(m_nsbins);
    double skjmax  = sqr(sqrt(Q2)-sqrt(m_mi2));
    double skjstep = (skjmax-m_smin_kj)/double(m_nsbins);
    double volume  = sijstep*skjstep, dgR;
    double sij     = m_smin_ij+sijstep/2., skj;
    do {
      skj = m_smin_kj+skjstep/2.;
      do {
	value += dgR = dg(Q2,sij,skj)/Q2 * volume;
	skj   += skjstep;
      } while(skj<skjmax);
      sij += sijstep;
    } while(sij<sijmax);
    //msg_Out()<<"Q2 = "<<q2<<": "<<value<<"\n";
    p_g->Fill(q2,value);
  }
}

double RunningWidth3_Resonance::CalculateWidth(const double & s) {
  return m_OSwidth * (*p_g)(s)/m_gV; }

Complex RunningWidth3_Resonance::BreitWigner(const double & s) {
  double MW = sqrt(s)*CalculateWidth(s);
  msg_Out()<<METHOD<<"("<<sqrt(s)<<"): for "<<sqrt(m_OSmass2)<<" --> "<<MW<<"\n";
  return m_OSmass2/Complex( m_OSmass2-s, -MW );
}
 
