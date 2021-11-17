#include "HADRONS++/Current_Library/VA_0_PPP.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

VA_0_PPP::VA_0_PPP(const ME_Parameters &parameters,const std::string& name) :
  Current_Base(parameters.flavs, parameters.indices, name)
{
  msg_Out()<<METHOD<<": "<<m_flavs[1]<<" "<<m_flavs[2]<<" "<<m_flavs[3]<<"\n"
	   <<" ** "<<m_flavs[p_i[0]]<<" ("<<p_i[0]<<") "
	   <<m_flavs[p_i[1]]<<" ("<<p_i[1]<<") "
	   <<m_flavs[p_i[2]]<<" ("<<p_i[2]<<")\n";
  
  p_ff = SelectFormFactor(parameters,name);
}

void VA_0_PPP::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti) {
  // this is a two-component vector for left-handed current and right-handed current:
  // J^\mu = \langle P(3)P(2)P(1) | \bar u Gamma^{\mu L} d | 0\rangle +
  //         \langle P(3)P(2)P(1) | \bar u Gamma^{\mu R} d | 0\rangle
  //       ~ { Q^\mu - [Q * (p1-p2)]/Q2 * (p1-p2)^\mu } * F_12 +
  //         { Q^\mu - [Q * (p1-p3)]/Q2 * (p1-p3)^\mu } * F_13 +
  //         { Q^\mu - [Q * (p2-p3)]/Q2 * (p2-p3)^\mu } * F_23 +
  //         i B_123 * epsilon^{\mu\nu\rho\sigma} p_{1\nu} p_{2\rho} p_{3\sigma} * F_V +
  //         Q^\mu * F_S
  // Comments:
  // We choose to make F_12 = 0 for all processes for the time being, so we can ignore it
  // at the moment.
  // Also, we push the term i B_123 into the form factor calculation to save one global parameter
  // which would be zero for some of the channels.
  // Of course, the right-handed component is always zero for EW interactions in
  // the Standard Model.
  Vec4D p1(moms[p_i[0]]), p2(moms[p_i[1]]), p3(moms[p_i[2]]), Q(p1+p2+p3);
  double Q2   = Q.Abs2();
  double s12  = (p1+p2).Abs2(), s13 = (p1+p3).Abs2(), s23 = (p2+p3).Abs2();
  Insert( m_global *
	  ( (Q - (Q*(p1-p3))/Q2 * (p1-p3)) * p_ff->F13(s12,s13,s23,Q2) +
	    (Q - (Q*(p2-p3))/Q2 * (p2-p3)) * p_ff->F23(s12,s13,s23,Q2) +
	    cross(p1,p2,p3) * p_ff->FV(s12,s13,s23,Q2) + Q * p_ff->FS(s12,s13,s23,Q2) ),
	  0.);
}


VA_0_PPP::FormFactor * VA_0_PPP::SelectFormFactor(const ME_Parameters &parameters,
					     const std::string& name) {
  if (name==std::string("VA_0_PPP_pipipi")) {
    return ThreePiFF(parameters);
  }
  msg_Error()<<"Error in "<<METHOD<<":\n"
	     <<"   No form factor model found for "<<name<<"\n"
	     <<"   Return NULL and hope for the Best.\n";
  return NULL;
}


VA_0_PPP::FormFactor * VA_0_PPP::ThreePiFF(const ME_Parameters &parameters) {
  // global prefactor as product of CKM element, Clebsch-Gordon factor, and
  // prefactor involving pion decay constant
  m_global = ( parameters.model("Vud", Tools::Vud) * SQRT_05 *
	       4./(3. * parameters.model("fpi", 0.1307 )) );
  switch (int(parameters.model("FORM_FACTOR",1))) {
  case 1:
    // Simple Kuehn-Santamaria Model,
    // https://link.springer.com/content/pdf/10.1007/BF01572024.pdf
    // as implemented in tauola: DOI: 10.1016/0010-4655(93)90061-G
    // we did not have this implemented before, but instead the
    // KS model from Decker, Finkemeier, Mirkes, hep-ph/9303474
    return new VA_0_PPP::FF_KS_pipipi(parameters.model,1);
  case 2:
    // KS model from Decker, Finkemeier, Mirkes, hep-ph/9303474
    return new VA_0_PPP::FF_KS_pipipi(parameters.model,2);
  case 3:
    // KS model from Finkemeier, Mirkes, hep-ph/9503474
  case 4:
    // CLEO model: hep-ex/9902022
  case 5:
    // RChT model: Dumm, Pich, Portoles, hep-ph/0312183 and update eprint 0911.4436 [hep-ph]
  default:
    break;
  }
  msg_Error()<<"Error in "<<METHOD<<":\n"
	     <<"   No form factor found for ff model = "
	     <<int(parameters.model("FORM_FACTOR",1))<<" in 3 pion mode.\n"
	     <<"   Return NULL and hope for the Best.\n";
  return NULL;
}

////////////////////////////////////////////////////////////////////
//
// Form factors in Kuehn-Santamaria model for pi pi pi
//
////////////////////////////////////////////////////////////////////

VA_0_PPP::FF_KS_pipipi::FF_KS_pipipi(const GeneralModel & model,const int & version) :
  VA_0_PPP::FormFactor(model),
  m_version(version),
  m_a1s_norm(0.), m_rhos_norm(0.),
  m_mpi2(sqr(Flavour(kf_pi).Mass())),
  m_ma12(sqr(model("Mass_a1(1260)+",  1.251 )))
{
  FillResonances(model);
}

void VA_0_PPP::FF_KS_pipipi::FillResonances(const GeneralModel & model) {
  int running_width  = int( model("RUNNING_WIDTH", 1 ) );
  m_a1s.push_back(make_pair(ResonanceFlavour( kf_a_1_1260_plus,
					      model("Mass_a1(1260)+",  1.251 ),
					      model("Width_a1(1260)+", 0.599 ),
					      running_width ),
			      1.));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_a1s.begin();
       rit!=m_a1s.end();rit++) {
    m_a1s_norm += rit->second;
  }
  // we do not have to distinguish between rho^0 and \rho^+- --
  // their masses and widths etc. are practically identical ...
  // when summing over the lists, each term will become
  // weight/prefactor * Breit-Wigner "propagator" term
  // the calculaion of the Breit-Wigner term depends on whether the width is
  // momentum-dependent and is performed in the "ResonanceFlavour" class
  m_rhos.push_back(make_pair(ResonanceFlavour( kf_rho_770_plus,
					       model("Mass_rho(770)+",  0.7769 ),
					       model("Width_rho(770)+", 0.149 ),
					       running_width ),
			     1.));
  m_rhos.push_back(make_pair(ResonanceFlavour( kf_rho_1450_plus,
					       model("Mass_rho(1450)+",  1.363 ),
					       model("Width_rho(1450)+", 0.310 ),
					       running_width ),
			     model("beta", -0.167 ) ));
  m_rhos.push_back(make_pair(ResonanceFlavour( kf_rho_1700_plus,
					       model("Mass_rho(1700)+",  1.700 ),
					       model("Width_rho(1700)+", 0.235 ),
					       running_width ),
			     model("gamma", 0. ) ));
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_rhos.begin();
       rit!=m_rhos.end();rit++) {
    m_rhos_norm += rit->second;
  }
}

Complex VA_0_PPP::FF_KS_pipipi::F(const double & spipi,const double & Q2) {
  Complex fa1s(0.,0.);
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_a1s.begin();
       rit!=m_a1s.end();rit++) {
    fa1s += rit->second * rit->first.BreitWigner(Q2);
  }
  Complex frhos(0.,0.);
  for (list<pair<ResonanceFlavour, double> >::iterator rit=m_rhos.begin();
       rit!=m_rhos.end();rit++) {
    frhos += rit->second * rit->first.BreitWigner(spipi);
  }
  return fa1s*frhos/(m_a1s_norm*m_rhos_norm);
}

Complex VA_0_PPP::FF_KS_pipipi::
FS(const double & s12, const double & s13, const double & s23,const double & Q2) {
  Complex fs(0.,0.);
  if (m_version==2) {
    // this implements the scalar form factor in Eq (35) of hep-ph/9303474
    Complex a1term(0.,0.);
    for (list<pair<ResonanceFlavour, double> >::iterator rit=m_a1s.begin();
	 rit!=m_a1s.end();rit++) {
      a1term += rit->second * rit->first.BreitWigner(Q2);
    }
    Complex rho13term(0.,0.), rho23term(0.,0.);
    for (list<pair<ResonanceFlavour, double> >::iterator rit=m_rhos.begin();
	 rit!=m_rhos.end();rit++) {
      rho13term += rit->second * rit->first.BreitWigner(s13);
      rho23term += rit->second * rit->first.BreitWigner(s23);
    }    
    fs = m_mpi2/(Q2-m_mpi2) *
      ( Complex(1.,0.) +
	(Q2-m_ma12)/(m_ma12*Q2) * a1term * ((s12-s23)/2. * rho13term +
					    (s12-s13)/2. * rho23term) );
  }
  return fs;
}


//DEFINE_CURRENT_GETTER(VA_0_PPP,"VA_0_PPP")
DECLARE_GETTER(VA_0_PPP,"VA_0_PPP",Current_Base,ME_Parameters);			

Current_Base *ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PPP>::	
operator()(const ME_Parameters &parameters) const  {
  msg_Out()<<"GETTER for VA_0_PPP: "
  	   <<parameters.flavs[parameters.indices[0]]<<" "
  	   <<parameters.flavs[parameters.indices[1]]<<" "
  	   <<parameters.flavs[parameters.indices[2]]<<"\n"
  	   <<"Formfactor = "<<int(parameters.model("FORM_FACTOR",1) )<<"\n";
  if ( (parameters.flavs[parameters.indices[0]].Kfcode()==kf_pi &&
	parameters.flavs[parameters.indices[1]].Kfcode()==kf_pi &&
	parameters.flavs[parameters.indices[2]].Kfcode()==kf_pi_plus) ||
       (parameters.flavs[parameters.indices[0]].Kfcode()==kf_pi_plus &&
	parameters.flavs[parameters.indices[1]].Kfcode()==kf_pi_plus &&
	parameters.flavs[parameters.indices[2]].Kfcode()==kf_pi_plus)) {
    return new VA_0_PPP(parameters, "VA_0_PPP_pipipi");
  }
  if ( (parameters.flavs[parameters.indices[0]].Kfcode()==kf_pi_plus &&
	( parameters.flavs[parameters.indices[1]].Kfcode()==kf_K ||
	  parameters.flavs[parameters.indices[1]].Kfcode()==kf_K_S ||
	  parameters.flavs[parameters.indices[1]].Kfcode()==kf_K_L) &&
	parameters.flavs[parameters.indices[2]].Kfcode()==kf_pi)) {
    //return new VA_0_PPP(parameters, "VA_0_PPP_piKpi");
  }
  return new VA_0_PPP(parameters, "VA_0_PPP_pipipi");
  return NULL; 
}

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PPP>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi K $ \n\n"
    <<"Order: special, see tau decays. \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
