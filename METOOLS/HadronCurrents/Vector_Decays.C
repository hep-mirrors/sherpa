#include "METOOLS/HadronCurrents/Vector_Decays.H"
#include "METOOLS/HadronCurrents/Line_Shapes.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <set>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//
// Vector -> Pseudoscalar + Photon
//
// Important: assume REAL photon is always 2nd particle (with mass 0)
//
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

V_PGamma::V_PGamma(const Flavour & inflav,const vector<Flavour> & outflavs,
		   const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR)
{
  FixPrefactor();
}

const double V_PGamma::Calculate(const double & s) {
  m_lambda2 = sqr(s-m_decmasses2[0]);
  return Flux(s) * ME2(s) * PS_2(s);
}

const double V_PGamma::ME2(const double & s)  const {
  return m_prefactor * m_lambda2/s;
}

const double V_PGamma::PS_2(const double & s) const {
  return sqrt(m_lambda2)/(8.*M_PI*s);
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//
// Vector -> Pseudoscalar + Pseudoscalar
//
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

V_PP::V_PP(const Flavour & inflav,const vector<Flavour> & outflavs,
	   const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR)
{
  FixPrefactor();
}

const double V_PP::Calculate(const double & s) {
  m_lambda2 = ( sqr(s-m_decmasses2[0]-m_decmasses2[1])-
		4.*m_decmasses2[0]*m_decmasses2[1] )/s;
  return Flux(s) * ME2(s) * PS_2(s);
}

const double V_PP::ME2(const double & s)  const {
  return m_prefactor/3. * m_lambda2;
}

const double V_PP::PS_2(const double & s) const {
  return sqrt(m_lambda2)/(8.*M_PI);
}

////////////////////////////////////////////////////////////////////
//
// Vector -> Off-shell (Pseudo-)Scalar + Pseudoscalar 
//
// Important:
// - assume on-shell pseudoscalar is always last flavour
// - assume off-shell particle decays into the first n-1 particle
//   flavours
//
////////////////////////////////////////////////////////////////////

V_PoffP::V_PoffP(const ATOOLS::Flavour & inflav,
		 const std::vector<ATOOLS::Flavour> & outflavs,
		 const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR),
  p_table(NULL)
{
  m_mPS    = m_decmasses[m_decmasses.size()-2];
  m_mPS2   = sqr(m_mPS);
  m_mPmin2 = sqr(m_mmin-m_mPS);
}
  
V_PoffP::~V_PoffP() { if (p_table) delete p_table; }

void V_PoffP::Init(const vector<Propagator_Base * > & props,
		   const axis & Qrange) {
  Propagator_Base * prop = props[0];
  p_table = new OneDim_Table(Qrange);
  for (size_t qbin=0;qbin<p_table->GetAxis().m_nbins;qbin++) {
    double Q        = p_table->GetAxis().x(qbin)+p_table->GetAxis().m_xstep/2.;
    double value    = 0.;
    if (Q>m_mmin) {
      double Q2     = sqr(Q);
      double mPmax2 = sqr(Q-m_mPS), mPstep = (mPmax2-m_mPmin2)/double(200);
      double mP2    = m_mPmin2 + mPstep/2.; 
      do {
	/////////////////////////////////////////////////////////////////////
	// Last factor below is absolute value squared of (summed) BW propagator
	// numerator, each term individually corrected to (M*Gamma) before squaring
	// the sum.
	// The factor mPstep is the "integration" volume over the propagator
	// invariant mass, which scales linearly in size with Q.
	/////////////////////////////////////////////////////////////////////
	value      += ( mPstep *
			(sqr(Q2-mP2-m_mPS2)-4.*m_mPS2*mP2)/Q2 *
			prop->Normalised2(mP2) );
	mP2        += mPstep;
      } while (mP2<mPmax2);
    }
    p_table->Fill(qbin,value);
  }
  delete prop;
  FixPrefactor();
}

const double V_PoffP::Calculate(const double & Q2) {
  if (!p_table)
    THROW(fatal_error,"Cannot call non-existing look-up table.");
  return m_prefactor*Flux(Q2)*(*p_table)(sqrt(Q2));  
}


////////////////////////////////////////////////////////////////////
//
// Vector -> Off-shell (Pseudo-)Vector + Pseudoscalar 
//
// Important:
// - assume on-shell pseudoscalar is always last flavour
// - assume off-shell particle decays into the first n-1 particle
//   flavours
//
////////////////////////////////////////////////////////////////////

V_VoffP::V_VoffP(const ATOOLS::Flavour & inflav,
		 const std::vector<ATOOLS::Flavour> & outflavs,
		 const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR),
  p_table(NULL)
{
  m_mPS    = m_decmasses[m_decmasses.size()-2];
  m_mPS2   = sqr(m_mPS);
  m_mVmin2 = sqr(m_mmin-m_mPS);
}
  
V_VoffP::~V_VoffP() { if (p_table) delete p_table; }

void V_VoffP::Init(const vector<Propagator_Base * > & props,
		   const axis & Qrange) {
  Propagator_Base * prop = props[0];
  p_table = new OneDim_Table(Qrange);
  for (size_t qbin=0;qbin<p_table->GetAxis().m_nbins;qbin++) {
    double Q        = p_table->GetAxis().x(qbin)+p_table->GetAxis().m_xstep/2.;
    double value    = 0.;
    if (Q>m_mmin) {
      double Q2     = sqr(Q);
      double mVmax2 = sqr(Q-m_mPS), mVstep = (mVmax2-m_mVmin2)/double(200);
      double mV2    = m_mVmin2 + mVstep/2.; 
      do {
	/////////////////////////////////////////////////////////////////////
	// Last factor below is absolute value squared of (summed) BW propagator
	// numerator, each term individually corrected to (M*Gamma) before squaring
	// the sum.
	// The factor mPstep is the "integration" volume over the propagator
	// invariant mass, which scales linearly in size with Q.
	/////////////////////////////////////////////////////////////////////
	value      += ( mVstep *
			(sqr(Q2-mV2-m_mPS2)-4.*m_mPS2*mV2)/Q2 *
			prop->Normalised2(mV2) );
	mV2        += mVstep;
      } while (mV2<mVmax2);
    }
    p_table->Fill(qbin,value);
  }
  delete prop;
  FixPrefactor();
}

const double V_VoffP::Calculate(const double & Q2) {
  if (!p_table)
    THROW(fatal_error,"Cannot call non-existing look-up table.");
  return m_prefactor*Flux(Q2)*(*p_table)(sqrt(Q2));  
}


////////////////////////////////////////////////////////////////////
//
// Vector -> Off-shell Axial-Vector + Pseudoscalar
//
// Important:
// - assume on-shell pseudoscalar is always last flavour
// - assume off-shell particle decays into the first n-1 particle
//   flavours
//
////////////////////////////////////////////////////////////////////

V_AoffP::V_AoffP(const ATOOLS::Flavour & inflav,
		 const std::vector<ATOOLS::Flavour> & outflavs,
		 const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR),
  p_table(NULL)
{
  m_mPS    = m_decmasses[m_decmasses.size()-2];
  m_mPS2   = sqr(m_mPS);
  m_mAmin2 = sqr(m_mmin-m_mPS);
}
  
V_AoffP::~V_AoffP() { if (p_table) delete p_table; }

void V_AoffP::Init(const vector<Propagator_Base * > & props,
		   const axis & Qrange) {
  Propagator_Base * prop = props[0];
  p_table = new OneDim_Table(Qrange);
  for (size_t qbin=0;qbin<p_table->GetAxis().m_nbins;qbin++) {
    double Q        = p_table->GetAxis().x(qbin)+p_table->GetAxis().m_xstep/2.;
    double value    = 0.;
    if (Q>m_mmin) {
      double Q2     = sqr(Q);
      double mAmax2 = sqr(Q-m_mPS), mAstep = (mAmax2-m_mAmin2)/double(200);
      double mA2    = m_mAmin2 + mAstep/2.; 
      do {
	/////////////////////////////////////////////////////////////////////
	// Last factor below is absolute value squared of (summed) BW propagator
	// numerator, each term individually corrected to (M*Gamma) before squaring
	// the sum.
	// The factor mPstep is the "integration" volume over the propagator
	// invariant mass, which scales linearly in size with Q.
	/////////////////////////////////////////////////////////////////////
	value      += ( mAstep *
			(sqr(Q2-mA2-m_mPS2)-4.*m_mPS2*mA2)/Q2 *
			prop->Normalised2(mA2) );
	mA2        += mAstep;
      } while (mA2<mAmax2);
    }
    p_table->Fill(qbin,value);
  }
  delete prop;
  FixPrefactor();
}

const double V_AoffP::Calculate(const double & Q2) {
  if (!p_table)
    THROW(fatal_error,"Cannot call non-existing look-up table.");
  return m_prefactor*Flux(Q2)*(*p_table)(sqrt(Q2));  
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//
// Vector -> Pseudoscalar + Pseudoscalar + Pseudoscalar
//
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

V_PPP::V_PPP(const Flavour & inflav,const vector<Flavour> & outflavs,
	     const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR),
  p_table(NULL), p_argument(NULL)
{ }

V_PPP::~V_PPP() {
  if (p_table) delete p_table;
}

void V_PPP::Init(const vector<Propagator_Base *> & props,
		 const axis & Qrange) {
  p_table = new OneDim_Table(Qrange);
  set<Flavour> vectors = { Flavour(kf_omega_782) };
  set<Flavour> axials  = { Flavour(kf_a_1_1260),
			   Flavour(kf_a_1_1260_plus) };
  if (vectors.find(m_inflav)!=vectors.end()) {
    p_argument = new V1minus_PPP_Arg(m_outflavs);
    p_argument->SetProps(props);
  }
  else if (axials.find(m_inflav)!=axials.end()) {
    p_argument = new V1plus_PPP_Arg(m_outflavs);
    p_argument->SetProps(props);
  }
  if (p_argument==NULL) {
    string error = "Fatal attempt at partial width for "+m_inflav.IDName();
    THROW(fatal_error, error);
  }
  FillTable();
  FixPrefactor();
}

const double V_PPP::Calculate(const double & s) {
  if (!p_table)
    THROW(fatal_error,"Cannot call non-existing look-up table.");
  return m_prefactor*Flux(s)*(*p_table)(sqrt(s));  
}

/////////////////////////////////////////////////////////////////////////
//
// Helper base class for the Dalitz decay V -> PPP mediated through
// vector propagators
//
/////////////////////////////////////////////////////////////////////////
V_PPP_Arg_Base::V_PPP_Arg_Base(const std::vector<Flavour> & outflavs) {
  if (outflavs.size()!=3) {
    string error = ( string("Wrong number of flavours for V_PPP: ")+
		     to_string(outflavs.size()));
    THROW(fatal_error,error); 
  }
  m_mi = outflavs[0].HadMass(); m_mi2 = sqr(m_mi);
  m_mj = outflavs[1].HadMass(); m_mj2 = sqr(m_mj);
  m_mk = outflavs[2].HadMass(); m_mk2 = sqr(m_mk);
  m_sij_min   = sqr(m_mi+m_mj);
  m_skj_min   = sqr(m_mk+m_mj);
  m_sik_min   = sqr(m_mi+m_mk);
  m_sum2      = m_mi2+m_mj2+m_mk2;
  m_sum2_ij   = m_mi2+m_mj2;    m_delta2_ij = m_mi2-m_mj2; 
  m_sum2_kj   = m_mk2+m_mj2;    m_delta2_kj = m_mk2-m_mj2; 
  m_sij_steps = m_skj_steps = 200;
}

void V_PPP_Arg_Base::FillTable(OneDim_Table * table) {
  for (size_t qbin=0;qbin<table->GetAxis().m_nbins;qbin++) {
    double Q       = table->GetAxis().x(qbin)+table->GetAxis().m_xstep/2.;
    double Q2      = sqr(Q);
    double value   = 0.;
    double sijmax  = sqr(Q-m_mk);
    double sijstep = (sijmax-m_sij_min)/double(m_sij_steps);
    double skjmax  = sqr(Q-m_mi);
    double skjstep = (skjmax-m_skj_min)/double(m_skj_steps);
    double volume  = sijstep*skjstep, dgR;
    double sij     = m_sij_min+sijstep/2., skj;
    do {
      skj = m_skj_min+skjstep/2.;
      do {
	value += dgR = dg(Q2,sij,skj)/Q2 * volume;
	skj   += skjstep;
      } while(skj<skjmax);
      sij += sijstep;
    } while(sij<sijmax);
    table->Fill(qbin,value);    
  }
}

/////////////////////////////////////////////////////////////////////////
//
// This helper encodes the structures relevant for a_1(1260)
//
/////////////////////////////////////////////////////////////////////////
void V1plus_PPP_Arg::SetProps(const vector<Propagator_Base *> & props) {
  if (props.size()!=2) {
    string error = ( string("Wrong number of flavours for V_PPP(1+): ")+
		     to_string(props.size()));
    THROW(fatal_error,error); 
  }
  p_propij  = props[0];
  p_propkj  = props[1];
}

double V1plus_PPP_Arg::dg(const double & Q2,
			  const double & sij,const double & skj) const {
  double  sik    = Q2-sij-skj+m_sum2;
  if (sik<=m_sik_min || sik>=sqr(sqrt(Q2)-sqrt(m_mj2))) return 0.;
  double  Vij2  = ( ( sij-2.*m_sum2_ij + sqr(m_delta2_ij)/sij ) + 
		    sqr(sik-skj+m_delta2_ij-
			(Q2+sij-m_mk2)*sqr(m_delta2_ij)/sij)/(4.*Q2) );
  double  Vkj2  = ( ( skj-2.*m_sum2_kj + sqr(m_delta2_kj)/skj ) + 
		    sqr(sik-sij+m_delta2_kj-
			(Q2+skj-m_mi2)*sqr(m_delta2_kj)/skj)/(4.*Q2) );
  double  Vijkj = ( ( (Q2-2.*sik+m_sum2-3.*m_mj2) -
		      (Q2-m_sum2+m_mj2)*(m_delta2_ij*m_delta2_kj)/(sij*skj) +
		      (Q2-2.*skj-m_sum2_kj+m_mi2)*m_delta2_ij/sij +
		      (Q2-2.*sij-m_sum2_ij+m_mk2)*m_delta2_kj/skj  )/2. +
		    ( (sik-skj+m_delta2_ij-
		       (Q2-2.*skj-m_sum2_ij+m_mk2)*m_delta2_ij/sij) *
		      (sik-sij+m_delta2_kj-
		       (Q2-2.*sij-m_sum2_kj+m_mi2)*m_delta2_kj/skj) )/
		    (4.*Q2) );
  Complex BWij = (*p_propij)(sij), BWkj = (*p_propkj)(skj);
  double  val = (Vij2  * norm(BWij) + Vkj2  * norm(BWkj) +
		 2.*Vijkj * abs(BWij*conj(BWkj)));
  return val;
}

/////////////////////////////////////////////////////////////////////////
//
// This helper encodes the structures relevant for omega(782)
//
/////////////////////////////////////////////////////////////////////////
void V1minus_PPP_Arg::SetProps(const vector<Propagator_Base *> & props) {
  if (props.size()!=3) {
    string error = ( string("Wrong number of flavours for V_PPP(1+): ")+
		     to_string(props.size()));
    THROW(fatal_error,error); 
  }
  p_propij  = props[0];
  p_propkj  = props[1];
  p_propik  = props[2];
}

double V1minus_PPP_Arg::dg(const double & Q2,
			   const double & sij,const double & skj) const {
  double  sik    = Q2-sij-skj+m_sum2;
  if (sik<=m_sik_min || sik>=sqr(sqrt(Q2)-sqrt(m_mj2))) return 0.;
  double pipj = (sij-m_mi2-m_mj2)/2.; 
  double pkpj = (skj-m_mk2-m_mj2)/2.; 
  double pipk = (sik-m_mi2-m_mk2)/2.;
  double num  = ( m_mi2*m_mj2*m_mk2 + 2.*pipj*pkpj*pipk -
		  ( m_mi2*sqr(pkpj) + m_mj2*sqr(pipk) + m_mk2*sqr(pipj) ) );
  double denom = norm((*p_propij)(sij)+(*p_propkj)(skj)+(*p_propik)(sik));
  return num*denom;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//
// (Axial-)Vector -> (Pseudo-)Scalar + Off-shell (Pseudo-)Scalar 
//
// We assume that the off-shell particle decays into two scalars,
// and we may capture interference effects between different
// resulting amplitudes.
// 
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

V_PPP_scalar::V_PPP_scalar(const Flavour & inflav,
			   const vector<Flavour> & outflavs,
			   const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR), p_table(NULL)
{
  for (size_t i=0;i<3;i++) p_props[i]=NULL;
}

  
V_PPP_scalar::~V_PPP_scalar() {
  if (p_table) delete p_table;
}
  
void V_PPP_scalar::Init(const vector<Propagator_Base *> & props,const axis & Qrange) {
  p_table = new OneDim_Table(Qrange);
  for (size_t i=0;i<props.size();i++) p_props[i] = props[i];
  FillTable();
  FixPrefactor();
}

void V_PPP_scalar::FillTable() {
  double sijmin    = sqr(m_decmasses[0]+m_decmasses[1]);
  double sjkmin    = sqr(m_decmasses[1]+m_decmasses[2]);
  m_sikmin         = sqr(m_decmasses[0]+m_decmasses[2]);
  m_sum2           = m_decmasses2[0]+m_decmasses2[1]+m_decmasses2[2];
  m_sum2_ij        = m_decmasses2[0]+m_decmasses2[1];
  m_sum2_jk        = m_decmasses2[1]+m_decmasses2[2];
  m_sum2_ik        = m_decmasses2[0]+m_decmasses2[2];
  size_t sij_steps = 200, sjk_steps = 200;
  for (size_t qbin=0;qbin<p_table->GetAxis().m_nbins;qbin++) {
    double Q        = p_table->GetAxis().x(qbin)+p_table->GetAxis().m_xstep/2.;
    double value    = 0.;
    if (Q>=m_decmasses[0]+m_decmasses[1]+m_decmasses[2]) {
      double Q2      = sqr(Q);
      double sijmax  = sqr(Q-m_decmasses[2]);
      double sijstep = (sijmax-sijmin)/double(sij_steps);
      double sjkmax  = sqr(Q-m_decmasses[0]);
      double sjkstep = (sjkmax-sjkmin)/double(sjk_steps);
      m_sikmax       = sqr(Q-m_decmasses[1]);
      double volume  = sijstep*sjkstep/Q2;
      double sij     = sijmin+sijstep/2., sjk, dgr;
      do {
	sjk = sjkmin+sjkstep/2.;
	do {
	  value += dgr = dg(Q2,sij,sjk) * volume/3.;
	  sjk   += sjkstep;
	} while (sjk<sjkmax);
	sij += sijstep;
      } while (sij<sijmax);
    }
    p_table->Fill(qbin,value);
  }
}

const double V_PPP_scalar::dg(const double & Q2,const double & sij,const double & sjk) const {
  double sik  = Q2-sij-sjk+m_sum2;
  if (sik<=m_sikmin || sik>=m_sikmax) return 0.;
  double dg = ( (sqr(Q2-m_decmasses2[0]-sjk)-4.*m_decmasses2[0]*sjk)/Q2 *
		norm((*p_props[0])(sjk)) );
  if (p_props[1]) {
    dg +=     ( (sqr(Q2-m_decmasses2[1]-sik)-4.*m_decmasses2[1]*sik)/Q2 *
		norm((*p_props[1])(sik)) );
  }
  if (p_props[2]) {
    dg +=     ( (sqr(Q2-m_decmasses2[2]-sij)-4.*m_decmasses2[2]*sij)/Q2 *
		norm((*p_props[2])(sij)) );
  }
  /*
  msg_Out()<<METHOD<<": "<<Q2<<" "<<sjk<<" ("<<m_decmasses2[0]<<") --> "
  	   <<(sqr(Q2-m_decmasses2[0]-sjk)-4.*m_decmasses2[0]*sjk)/Q2
  	   <<" * "<<norm((*p_props[0])(sjk))<<" from "<<(*p_props[0])(sjk)<<"\n";
  */
  return dg;
}

const double V_PPP_scalar::Calculate(const double & Q2) {
  if (!p_table)
    THROW(fatal_error,"Cannot call non-existing look-up table.");
  return m_prefactor*Flux(Q2)*(*p_table)(sqrt(Q2));  
}


