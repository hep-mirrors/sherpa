#include "HADRONS++/Current_Library/VA_0_PP.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//
// If you want to test the formfactors, you have to work with the .dat files
// for the individual decays - for example, if you want to test K+pi final states
// you need to replace the "VA_0_PP_strange" name tag in the .dat file with
// "VA_0_PP".
// Also, at the moment we do not yet transfer parameters from the dat files into
// the code - so you need to explicitly change the defaults in the model("name",default)
// structures - see the comments marked with *** in the code.
//
/////////////////////////////////////////////////////////////////////////////////

VA_0_PP::VA_0_PP(const ME_Parameters &parameters,const std::string& name) :
  Current_Base(parameters.flavs, parameters.indices, name)
{
  double CG = 1.;
  if (name==std::string("VA_0_PP_pipi")) {
    // *** make this model("Form_Factor",1) a model("Form_Factor",2)
    if (int(parameters.model("FORM_FACTOR",1))==1)
      p_ff = new VA_0_PP::FF_KS_pipi(parameters.model);
    else
      p_ff = new VA_0_PP::FF_RChT_pipi(parameters.model);
    if  (m_flavs[p_i[1]].Kfcode()==kf_pi_plus) CG = 1./SQRT_05;
    m_global  = parameters.model("Vud", Tools::Vud) * CG;
    m_DeltaM2 = 0.;
  }
  else if (name==std::string("VA_0_PP_Kpi") ||
	   name==std::string("VA_0_PP_piK") ) {
    if (int(parameters.model("FORM_FACTOR",1))==1)
      p_ff = new VA_0_PP::FF_KS_Kpi(parameters.model);
    /*else
      p_ff = new VA_0_PP::FF_RChT_Kpi(parameters.model);*/ //fix 
    // the Clebsch-gordon factos between pi pi and pi K look a bit odd.
    // we need to cross check this with literature.
    if  (m_flavs[p_i[1]].Kfcode()==kf_pi ||
	 m_flavs[p_i[0]].Kfcode()==kf_pi ) CG = SQRT_05;
    m_global  = parameters.model("Vus", Tools::Vus) * CG;
    m_DeltaM2 = dabs(sqr(m_flavs[p_i[1]].HadMass())-sqr(m_flavs[p_i[0]].HadMass())); 
  }
}
  
void VA_0_PP::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  double q2 = (moms[p_i[1]] + moms[p_i[0]] ).Abs2();
  // this is a two-component vector for left-handed current and right-handed current:
  // J^\mu = \langle P(1)P(0) | \bar u Gamma^{\mu L} d | 0\rangle +
  //         \langle P(1)P(0) | \bar u Gamma^{\mu R} d | 0\rangle
  //       ~ (p_1-p_0)^\mu * F_V(q^2) +
  //         |m_1^2-m_0^2|/q^2 * (p_1+p_0)^\mu * [F_S(q^2)-F_V(q^2)]
  // If we had more form factors we would have to expand this equation,
  // and produce the form factors in the FormFactor base class and their
  // KS and RChT realisations.
  // Of course, the right-handed component is always zero for EW interactions in
  // the Standard Model.
  Insert( m_global * ( (moms[p_i[1]]-moms[p_i[0]]) * p_ff->F_V(q2) +
		       m_DeltaM2/q2 * (moms[p_i[1]]+moms[p_i[0]]) *
		       (p_ff->F_S(q2)-p_ff->F_V(q2)  
		      ) ),
	  0. );
}


////////////////////////////////////////////////////////////////////
//
// Vector form factor in Kuehn-Santamaria model for pi+pi/K+K
// https://link.springer.com/article/10.1007/BF01572024
//
////////////////////////////////////////////////////////////////////

VA_0_PP::FF_KS_pipi::FF_KS_pipi(const GeneralModel & model) :
  VA_0_PP::FormFactor(model) {
  FillResonances(model);
}

void VA_0_PP::FF_KS_pipi::FillResonances(const GeneralModel & model) {
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  // constructor of the ResonanceFavour needs: PDG code, mass, width, and whether the
  // width is running.  we also need the "weight" or prefactor in the combination
  // of the resonance contributions in the form factor.
  // *** change the masses/widths/prefactors in the code if you want to.
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_770_plus,
						      model("Mass_rho(770)+",  0.7769 ),
						      model("Width_rho(770)+", 0.149 ),
						      running_width ),
				    1.));
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_1450_plus,
						      model("Mass_rho(1450)+",  1.363 ),
						      model("Width_rho(1450)+", 0.310 ),				     
						      running_width ),
				    model("beta", -0.167 ) ));
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_1700_plus,
						      model("Mass_rho(1700)+",  1.700 ),
						      model("Width_rho(1700)+", 0.235 ),				     
						      running_width ),
				    model("gamma", -0.050 ) ));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    m_Vresonance_norm += rit->second;
  }
}

Complex VA_0_PP::FF_KS_pipi::F_V(const double & q2) {
  Complex fv(0.,0.);
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    // each term is: weight/prefactor * Breit-Wigner "propagator" term
    // the calculaion of the second term depends on whether the width is
    // momentum-dependent and is performed in the "ResonanceFlavour" class
    fv += rit->second * rit->first.BreitWigner(q2);
  }
  return fv/m_Vresonance_norm;
}

////////////////////////////////////////////////////////////////////
//
// Vector and scalar form factors in Kuehn-Santamaria model for K+pi
// https://link.springer.com/article/10.1007/BF01572024
//
////////////////////////////////////////////////////////////////////

VA_0_PP::FF_KS_Kpi::FF_KS_Kpi(const GeneralModel & model) :
  VA_0_PP::FormFactor(model) {
  FillResonances(model);
}

void VA_0_PP::FF_KS_Kpi::FillResonances(const GeneralModel & model) {
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_K_star_892_plus,
						      model("Mass_K*(892)+",  0.8921 ),
						      model("Width_K*(892)+", 0.0513 ),
						      running_width ),
				    1.));
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_K_star_1680_plus,
						      model("Mass_K*(1680)+",  1.7000 ),
						      model("Width_K*(1680)+", 0.2350 ),
						      running_width ),
				    model("beta", -0.038)));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    m_Vresonance_norm += rit->second;
  }
  m_Sresonances.push_back(make_pair(ResonanceFlavour( kf_K_0_star_1430_plus,
						      model("Mass_K(0)*(1430)+",  1.396 ),
						      model("Width_K(0)*(1430)+", 0.294 ),				     
						      running_width ),
				    1.));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Sresonances.begin();
       rit!=m_Sresonances.end();rit++) {
    m_Sresonance_norm += rit->second;
  }
}

Complex VA_0_PP::FF_KS_Kpi::F_V(const double & q2) {
  Complex fv(0.,0.);
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    fv += rit->second * rit->first.BreitWigner(q2);
  }
  return fv/m_Vresonance_norm;
}

Complex VA_0_PP::FF_KS_Kpi::F_S(const double & q2) {
  Complex fs(0.,0.);
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Sresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    fs += rit->second * rit->first.BreitWigner(q2);
  }
  return fs/m_Sresonance_norm;
}

////////////////////////////////////////////////////////////////////
//
// Vector form factor in Resonance Chiral Theory for pi+pi/K+K
//
// we need to implement https://arxiv.org/pdf/1902.02273.pdf:
// this means more resonances in the for-loop - prepared in FillResonances
//
////////////////////////////////////////////////////////////////////

VA_0_PP::FF_RChT_pipi::FF_RChT_pipi(const GeneralModel & model) :
  VA_0_PP::FormFactor(model)
{
  // Funny definition of f_pi here.
  // I assume it is 130.7 MeV/\sqrt{2} (this is a question of different conventions) --
  // but we need to check this!
  m_fpi   = SQRT_05 * model("fpi", 0.1307 ); //### eq 2.1 and page 7
  m_m2_pi = sqr( Flavour(kf_pi_plus).HadMass() );
  m_m2_K  = sqr( Flavour(kf_K_plus).HadMass() );
  FillResonances(model);  
}

void VA_0_PP::FF_RChT_pipi::FillResonances(const GeneralModel & model) {
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_770_plus,
						      model("Mass_rho(770)+",  0.7769 ),
						      model("Width_rho(770)+", 0.149 ),
						      running_width ),
				    1.));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    msg_Out()<<"* "<<rit->first.KfCode()<<" --> "<<rit->second<<"\n";
    m_Vresonance_norm += rit->second;
  }
}

Complex VA_0_PP::FF_RChT_pipi::F_V(const double & q2) {
  Complex sum(0.,0.);
  double  prefactor = -q2/(96.*sqr(M_PI*m_fpi)); 
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    Complex pi_loop   = A_Loop(m_m2_pi/q2, m_m2_pi/rit->first.Mass2() ); 
    Complex K_loop    = A_Loop(m_m2_K/q2,  m_m2_K/rit->first.Mass2() );
    Complex sum_loop  = pi_loop + 0.5*K_loop;
    Complex term;
    if (rit->first.Running()) {
      double running_width = -rit->second * rit->first.Mass2()  * prefactor * sum_loop.imag();
      term = Tools::BreitWigner(q2, rit->first.Mass2(), running_width);
    }
    else {
      term = Tools::BreitWignerFix(q2, rit->first.Mass2(), rit->first.MassWidth() );
    }
    sum +=  term * exp(prefactor*sum_loop.real());
  }
  return sum;
}

Complex VA_0_PP::FF_RChT_pipi::A_Loop( double x, double y )
{
  Complex sigma = csqrt(1.-4.*x);
  return ( log(y) + 8.*x - 5./3. + pow(sigma,3.)*log((sigma+1.)/(sigma-1.)) ); 
}



////////////////////////////////////////////////////////////////////
//
// Vector and scalar form factors in Resonance Chiral Theory for K+pi
// https://arxiv.org/pdf/1902.02273.pdf
//
////////////////////////////////////////////////////////////////////

VA_0_PP::FF_RChT_Kpi::FF_RChT_Kpi(const GeneralModel & model) :
  VA_0_PP::FormFactor(model)
{
  m_m2_pi   = sqr( Flavour(kf_pi_plus).HadMass() );
  m_m2_K    = sqr( Flavour(kf_K_plus).HadMass() );
  m_m2_eta  = sqr( Flavour(kf_eta).HadMass() );
  m_DeltaM2 = m_m2_K - m_m2_pi; 
  m_SigmaM2 = m_m2_K + m_m2_pi; 
  m_piK_threshold  = Flavour(kf_pi).HadMass()  + Flavour(kf_K_plus).HadMass();
  m_etaK_threshold = Flavour(kf_eta).HadMass() + Flavour(kf_K_plus).HadMass();
  // Different definition of f_pi here compared to the pi pi formfactor above.
  // We need to check this!
  //m_fpi    = 1./SQRT_05 * model("fpi", 0.1307 );  
  m_fpi    = SQRT_05 * model("fpi", 0.1307);  //## eq. 2.1 https://arxiv.org/pdf/1902.02273.pdf
  m_fpi2   = sqr(m_fpi);
  m_mu2    = sqr( model("renorm", model("Mass_rho(770)+",
					Flavour(kf_rho_770_plus).HadMass())));
  m_cd     = model("const_cd", 0.014);
  m_cm     = m_fpi2/(4.*m_cd);
  FillResonances(model);  
}

void VA_0_PP::FF_RChT_Kpi::FillResonances(const GeneralModel & model) {
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_K_star_892_plus,
						      model("Mass_K*(892)+",  0.8921 ),
						      model("Width_K*(892)+", 0.0513 ),
						      running_width ),
				    1.));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    m_Vresonance_norm += rit->second;
  }
  m_Sresonances.push_back(make_pair(ResonanceFlavour( kf_K_0_star_1430_plus,
						      model("Mass_K(0)*(1430)+",  1.396 ),
						      model("Width_K(0)*(1430)+", 0.294 ),				     
						      running_width ),
				    1.));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Sresonances.begin();
       rit!=m_Sresonances.end();rit++) {
    m_Sresonance_norm += rit->second;
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////
//
// the form factors need to be filled with the terms in VA_0_PP_strange and checked against
// the underlying literature.
// A lot of the literature is authored by someone called P. Roig.
// https://inspirehep.net/literature?sort=mostrecent&size=25&page=1&q=f%20a%20roig%2C%20p&ui-citation-summary=true
//
//////////////////////////////////////////////////////////////////////////////////////////////
/*
//##
Complex HADRONS::VA_0_PP::FF_RChT_Kpi::F_V(const double & q2, double s) {
  Complex sum(1.,0.);
  double    MG_K = MassWidthVector(s);
  Complex     BW = Tools::BreitWigner( s, m_mK2, MG_K );
  Complex M_part = Mr(s,m_mK2,m_mPi2) + Mr(s,m_mK2,m_mEta2);
  Complex L_part = L(s,m_mK2,m_mPi2) + L(s,m_mK2,m_mEta2);
  double   expon = 3./2./m_fpi2 *( s*M_part.real() - L_part.real() );
  sum = BW * exp(expon);
  return sum;
}

//##
Complex HADRONS::VA_0_PP::FF_RChT_Kpi::F_S(const double & q2, double s) {
  Complex sum(1.,0.);
  double  MG_K0 = MassWidthScalar(s);
  double Mass2 = 3; //delete
  double MassWidth = 2; //delete
  Complex    BW = HADRONS::Tools::BreitWigner( s, Mass2, MassWidth);
  Complex    F4 = (
  	1./(8.*m_fpi2) *
      		( 5.*s - 2.*m_Sigma_KP - 3.*sqr(m_Delta_KP)/s ) *
     		JBar(s,m_mK2,m_mPi2,m_mK2+m_mPi2,m_mK2-m_mPi2)
	+ 1./(24.*m_fpi2) *
      		( 3.*s - 2.*m_Sigma_KP - sqr(m_Delta_KP)/s ) *
      		JBar(s,m_mK2,m_mEta2,m_mK2+m_mEta2,m_mK2-m_mEta2)
	+ s/(4.*m_Delta_KP)*(5.*MuOf(m_mPi2)-2.*MuOf(m_mK2)-3.*MuOf(m_mEta2)));
  double  inter = 1. - (1.-m_fpi2/4./sqr(m_cd))*m_Sigma_KP/m_MK02;
  Complex  expon = Complex( F4.real(), F4.imag()/(1.+sqr(F4.imag())) );
  sum = BW * inter * exp(expon);
  return sum;
}
*/

////////////////////////////////////////////////////////////////////
//
// Pion Vector Form Factor in Resonance Chiral Theory for pi+pi
// (A dispersive analysis of the pion vector form factor 
// and τ − → KK S ν τ deca(y)
//  https://arxiv.org/abs/1902.02273
////////////////////////////////////////////////////////////////////

//##
VA_0_PP::FF_RChT2_pipi::FF_RChT2_pipi(const GeneralModel& model) :
  VA_0_PP::FormFactor(model)
{
  m_fpi  = SQRT_05 * model("fpi", 0.1307 ); 
  m_m_pi = Flavour(kf_pi_plus).HadMass(); m_m2_pi = sqr(m_m_pi);
  m_m_K  = Flavour(kf_K_plus).HadMass();  m_m2_K  = sqr(m_m_K);
  FillResonances(model);  
}

 void VA_0_PP::FF_RChT2_pipi::FillResonances(const GeneralModel & model) { //2.35
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_770_plus,
						      model("Mass_rho(770)+",  0.7736 ),
						      model("Width_rho(770)+", 0.149 ),
						      running_width ),
				    1.));
  m_Vresonances.back().first.SetMassFactor(1.);
  m_Vresonances.back().first.SetSAmplitudesAndPhases(model("gamma",0.15),model("delta",-0.12),
						     model("phi1",-0.36),model("phi2",-0.02));
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_1450_plus,
						      model("Mass_rho(1450)+",  1.376 ),
						      model("Width_rho(1450)+", 0.603 ),
						      running_width ),
				    1.));
  m_Vresonances.back().first.SetMassFactor(0);
  m_Vresonances.back().first.SetSAmplitudesAndPhases(-model("gamma",0.15),0.,
						     model("phi1",-0.36),0.);
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_1700_plus,
						      model("Mass_rho(1700)+",  1.718 ),
						      model("Width_rho(1700)+", 0.465 ),
						      running_width ),
				    1.));
  m_Vresonances.back().first.SetMassFactor(0.);
  m_Vresonances.back().first.SetSAmplitudesAndPhases(0.,-model("delta",-0.12),
						     0.,model("phi2",-0.02));

  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    msg_Out()<<"* "<<rit->first.KfCode()<<" --> "<<rit->second<<"\n";
  }
  m_Vresonance_norm = 1.;
}

Complex VA_0_PP::FF_RChT2_pipi::Sigma(const double & m2,const double & s) {
  return sqrt(1.-4.*m2/s);
}

Complex VA_0_PP::FF_RChT2_pipi::CalcA(const double & m2,const double & s,const double & mu2) {
  Complex sigma = Sigma(m2,s);
  return (mu2<0. ? 0. : log(m2/mu2)) + 8.*m2/s - 5./3. + pow(sigma,3.)*log((sigma+1.)/(sigma-1.));
}

Complex VA_0_PP::FF_RChT2_pipi::F_V(const double & q2) {
  Complex sum(0.,0.);
  Complex A_K  = CalcA(m_m2_K,q2,q2);
  Complex A_pi = CalcA(m_m2_pi,q2,q2);

  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    Complex BW = rit->second * rit->first.BreitWignerPhases(q2), arg;
    if (rit->first.KfCode()==kf_rho_770_plus) {
      arg = real(q2/(96.*sqr(M_PI*m_fpi))*(A_pi+A_K/2.));
    }
    else {
      arg  = (q2*rit->first.TwoBodyResonanceMassWidth(q2,m_m_pi)/
	      (M_PI * sqr(rit->first.Mass2()) * pow(Sigma(rit->first.Mass2(),q2),3))*real(A_pi));
    }
    sum += BW * exp(-arg);
  }
  return sum;
}


////////////////////////////////////////////////////////////////////
//
// Gounaris-Sakurai parametrization Vector Form Factor in
// Resonance Chiral Theory for pi+pi https://arxiv.org/abs/1509.09140
//
////////////////////////////////////////////////////////////////////

//##


VA_0_PP::FF_GS_pipi::FF_GS_pipi(const GeneralModel & model) :
  VA_0_PP::FormFactor(model)
{
  m_fpi   = SQRT_05 * model("fpi", 0.1307 ); 
  m_m2_pi = sqr( Flavour(kf_pi_plus).HadMass() );
  m_m2_K  = sqr( Flavour(kf_K_plus).HadMass() );
  FillResonances(model);  
}

void VA_0_PP::FF_GS_pipi::FillResonances(const GeneralModel & model) {
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_770_plus,
						      model("Mass_rho(770)+",  0.7746 ), //table VII https://journals.aps.org/prd/pdf/10.1103/PhysRevD.78.072006
						      model("Width_rho(770)+", 0.148 ),
						      running_width ),
				    1.));
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_1450_plus,
						      model("Mass_rho(1450)+",  1.446 ),
						      model("Width_rho(1450)+", 0.434 ),				     
						      running_width ),
				    model("beta", 0.15 ) )); 
  m_Vresonances.push_back(make_pair(ResonanceFlavour( kf_rho_1700_plus,
						      model("Mass_rho(1700)+",  1.728 ),
						      model("Width_rho(1700)+", 0.164 ),				     
						      running_width ),
				    model("gamma", 0.037 ) ));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_Vresonances.begin();
       rit!=m_Vresonances.end();rit++) {
    msg_Out()<<"* "<<rit->first.KfCode()<<" --> "<<rit->second<<"\n";
    m_Vresonance_norm += rit->second;
  }
}
Complex VA_0_PP::FF_GS_pipi::F_V(const double & q2) {
  Complex sum(0.,0.);
  term1 = Tools::BreitWignerrhoprime(q2, 1.446, 0.434);
  term2 = Tools::BreitWignerrhobiprime(q2,1.728,0.164);
  complex sum += (1/(1+beta+gamma))*(BWGS+beta*term1+gamma*term2);
  return sum;
}




//DEFINE_CURRENT_GETTER(VA_0_PP,"VA_0_PP")

DECLARE_GETTER(VA_0_PP,"VA_0_PP",Current_Base,ME_Parameters);			

Current_Base *ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PP>::	
operator()(const ME_Parameters &parameters) const  {
  //msg_Out()<<"GETTER for VA_0_PP: "
  //	   <<parameters.flavs[parameters.indices[1]]<<" "
  //	   <<parameters.flavs[parameters.indices[0]]<<"\n"
  //	   <<"Formfactor = "<<int(parameters.model("FORM_FACTOR",1) )<<"\n";
  if ( (parameters.flavs[parameters.indices[0]].Kfcode()==kf_pi &&
	parameters.flavs[parameters.indices[1]].Kfcode()==kf_pi_plus) ||
       ( (parameters.flavs[parameters.indices[0]].Kfcode()==kf_K_S ||
	  parameters.flavs[parameters.indices[0]].Kfcode()==kf_K_L ||
	  parameters.flavs[parameters.indices[0]].Kfcode()==kf_K) &&
	 parameters.flavs[parameters.indices[1]].Kfcode()==kf_K_plus) ) {
    return new VA_0_PP(parameters, "VA_0_PP_pipi");
  }
  else if ( (parameters.flavs[parameters.indices[0]].Kfcode()==kf_pi_plus &&
	     (parameters.flavs[parameters.indices[1]].Kfcode()==kf_K_S ||
	      parameters.flavs[parameters.indices[1]].Kfcode()==kf_K_L ||
	      parameters.flavs[parameters.indices[1]].Kfcode()==kf_K)) ||
	    (parameters.flavs[parameters.indices[0]].Kfcode()==kf_pi &&
	     parameters.flavs[parameters.indices[1]].Kfcode()==kf_K_plus) ) { 
    return new VA_0_PP(parameters, "VA_0_PP_piK");
  }
  else if ( ((parameters.flavs[parameters.indices[0]].Kfcode()==kf_K_S ||
	      parameters.flavs[parameters.indices[0]].Kfcode()==kf_K_L ||
	      parameters.flavs[parameters.indices[0]].Kfcode()==kf_K) &&
	     parameters.flavs[parameters.indices[1]].Kfcode()==kf_pi_plus) ||
	    (parameters.flavs[parameters.indices[0]].Kfcode()==kf_K_plus &&
	     parameters.flavs[parameters.indices[1]].Kfcode()==kf_pi) ) { 
    return new VA_0_PP(parameters, "VA_0_PP_Kpi");
  }
  return NULL; 
}

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PP>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^\\pm$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}

