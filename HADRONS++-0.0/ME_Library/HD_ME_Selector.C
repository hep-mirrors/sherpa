#include "HD_ME_Selector.H"
#include "Message.H"

#include "Two_Body_MEs.H"
#include "Three_Body_MEs.H"


using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

HD_ME_Base * HD_ME_Selector::GetME(int nin,int nout,Flavour * flavs,
				   std::string met) {
  HD_ME_Base * hdme = NULL;
  double mass = flavs[0].Mass();
  for (int i=1;i<1+nout;i++) {
    mass-=flavs[i].Mass();
    if (mass<=0.) {
      msg.Error()<<"Error in HD_ME_Selector::GetME("<<nin<<"->"<<nout<<") : "
		 <<"   Masses don not match : "<<flavs[0].Mass()<<" -> ";
      for (int j=1;j<nout;j++) msg.Error()<<flavs[j].Mass()<<"+";
      msg.Error()<<flavs[nout].Mass()<<endl
		 <<"   Will return NULL and hope for the best."<<endl;
      return hdme;
    }
  }

  switch (flavs[0]) {
  case (kf::pi):
  case (kf::eta):
  case (kf::eta_prime_958):
    SelectLightPseudoScalarDecay(nout,flavs,hdme);
    break;
  case (kf::rho_770):
  case (kf::rho_770_plus):
  case (kf::omega_782):
  case (kf::Kstar_892):
  case (kf::Kstar_892_plus):
  case (kf::phi_1020):
    //SelectLightVectorDecay(nout,flavs,hdme);
    break;
  }

  if (hdme==NULL) hdme = new Isotropic(nout,flavs);
  return hdme;
}



void HD_ME_Selector::SelectLightPseudoScalarDecay(int nout,Flavour * flavs,
						  HD_ME_Base *& hdme)
{
  switch (nout) {
  case 2:
    if (flavs[1]==Flavour(kf::photon) &&
	flavs[2]==Flavour(kf::photon)) 
      hdme = new P_2Gamma(nout,flavs);
    break;
  case 3: 
    if (flavs[1].IsLepton() &&
	flavs[2].IsLepton() && flavs[2].IsAnti() &&
	flavs[3]==Flavour(kf::photon)) {
      hdme = new P_GammaFF(nout,flavs);
    }
    if (IsPseudoScalar(flavs[1]) &&
	flavs[2]==Flavour(kf::photon) &&
	flavs[3]==Flavour(kf::photon)) { 
      hdme = new P_P2Gamma(nout,flavs);
    }
    if (flavs[1]==Flavour(kf::pi_plus) &&
	flavs[2]==Flavour(kf::pi_plus).Bar() &&
	flavs[3]==Flavour(kf::photon)) {
      hdme = new P_2PGamma(nout,flavs);
      SetVector_For_2PS(nout,1,2,flavs,hdme);
    }
    if (IsPseudoScalar(flavs[1]) &&
	IsPseudoScalar(flavs[2]) &&
	IsPseudoScalar(flavs[3])) {
      hdme = new P_3P_Dalitz(nout,flavs);
      SetDalitz_For_3P(nout,flavs,hdme);
    }
    break;
  }
}

void HD_ME_Selector::SetVector_For_2PS(int nout,int PS1, int PS2,
				       Flavour * flavs,HD_ME_Base * hdme)
{
  double mass   = flavs[0].Mass();
  double charge = flavs[PS1].Charge()+flavs[PS2].Charge();
  if (charge==0.) {
    if (flavs[0].Mass()>Flavour(kf::rho_770).Mass())    
      hdme->AddVector(Flavour(kf::rho_770).Mass(),Flavour(kf::rho_770).Width());
    if (flavs[0].Mass()>Flavour(kf::omega_782).Mass())
      hdme->AddVector(Flavour(kf::omega_782).Mass(),Flavour(kf::omega_782).Width());
    if (flavs[0].Mass()>Flavour(kf::K_star_892).Mass())
      hdme->AddVector(Flavour(kf::K_star_892).Mass(),Flavour(kf::K_star_892).Width());
    if (flavs[0].Mass()>Flavour(kf::phi_1020).Mass())
      hdme->AddVector(Flavour(kf::phi_1020).Mass(),Flavour(kf::phi_1020).Width());
  }
  else if (charge!=0.) {
    if (flavs[0].Mass()>Flavour(kf::rho_770).Mass())
      hdme->AddVector(Flavour(kf::rho_770).Mass(),Flavour(kf::rho_770).Width());
    if (flavs[0].Mass()>Flavour(kf::K_star_892).Mass())
      hdme->AddVector(Flavour(kf::K_star_892).Mass(),Flavour(kf::K_star_892).Width());
  }
}

void HD_ME_Selector::SetDalitz_For_3P(int nout,Flavour * flavs,
				      HD_ME_Base * hdme)
{
}
