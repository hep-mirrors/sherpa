#include "Flavour.H"
#include "HD_ME_Selector.H"
#include "Message.H"

#include "Tau_Decay_MEs.H"
#include "Two_Body_MEs.H"
#include "Three_Body_MEs.H"
#include "Four_Body_MEs.H"


using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

HD_ME_Base * HD_ME_Selector::GetME(int nin,int nout,Flavour * flavs)
{
  HD_ME_Base * hdme = NULL;							// pointer on ME_Base
  double mass = flavs[0].Mass();					// mass of decaying particle
  if (flavs[0]==Flavour(kf_tau) || 
      flavs[0]==Flavour(kf_tau)) mass = flavs[0].PSMass();

  // sanity check if sum of outgoing masses > incoming mass
  for (int i=1;i<1+nout;i++) {
    mass-=flavs[i].Mass();
    if (mass<0.) {
      msg_Error()<<"Error in HD_ME_Selector::GetME("<<nin<<"->"<<nout<<") : "
		 <<"   Masses do not match : "<<flavs[0].Mass()<<" -> ";
      for (int j=1;j<nout;j++) msg_Error()<<flavs[j].Mass()<<"+";
      msg_Error()<<flavs[nout].Mass()<<endl
		 <<"   Will return NULL and hope for the best."<<endl;
      return hdme;
    }
  }

  // Select ME depending on decaying particle
  switch (flavs[0]) {
    case (kf_tau):
      SelectTauDecay(nout,flavs,hdme);
      break;
    case (kf_eta):
    case (kf_eta_prime_958):
      SelectLightPseudoScalarDecay(nout,flavs,hdme);
      break;
  }

  if (hdme==NULL) hdme = new Isotropic(nout,flavs);
  return hdme;
}


// Select the corresponding tau decay ME
void HD_ME_Selector::SelectTauDecay(int nout,Flavour * flavs,
				    HD_ME_Base *& hdme )
{
  msg_Tracking()<<"HD_ME_Selector::SelectTauDecay 1->"<<nout<<endl;
  switch( nout ) {
    case 2:
      {
        int nPion(0), nKaon(0);
        for( int i=1; i<3; i++ ) {
          if( flavs[i].Kfcode() == kf_pi_plus ) nPion++;
          if( flavs[i].Kfcode() == kf_K_plus ) nKaon++;
        }
        if( nPion == 1 ) {
          hdme = new Tau_Pseudo( nout, flavs );
        }
        if( nKaon == 1 ) {
          hdme = new Tau_Pseudo( nout, flavs );
        }
        break;
      }
    case 3:
      {
        int nLep(0), nPion(0), nKaon(0);
        for( int i=1; i<4; i++ ) {
          if( flavs[i].IsLepton() ) nLep++;
          if( flavs[i].Kfcode() == kf_pi_plus ||
              flavs[i].Kfcode() == kf_pi ) nPion++;
          if( flavs[i].Kfcode() == kf_K_plus ||
              flavs[i].Kfcode() == kf_K_S ||
              flavs[i].Kfcode() == kf_K_L ) nKaon++;
        }
        if( nLep == 3 ) {
          hdme = new Tau_Lepton( nout, flavs );
        }
        if( nPion == 2 ) {
          hdme = new Tau_Two_Pion( nout, flavs );
        }
        if( nKaon == 2 ) {
          hdme = new Tau_Two_Pion( nout, flavs );
        }
        if( nKaon == 1 && nPion == 1 ) {
          hdme = new Tau_Pion_Kaon( nout, flavs );
        }
        break;
      }
    case 4:
      {
        int nPseudo(0), nEta(0);
        // count number of pseudoscalars
        for( int i=1; i<5; i++ ) {
          if( flavs[i].Kfcode() == kf_pi_plus ||
              flavs[i].Kfcode() == kf_pi ||
              flavs[i].Kfcode() == kf_K_plus ||
              flavs[i].Kfcode() == kf_K_L ||
              flavs[i].Kfcode() == kf_K_S) nPseudo++;
          if( flavs[i].Kfcode() == kf_eta ) nEta++;
        }
        if( nPseudo==3 )
          hdme = new Tau_Three_Pseudo( nout, flavs );
        if( nEta==1 && nPseudo==2 )
          hdme = new Tau_Eta_Two_Pion( nout, flavs );
	break;
      }
    case 5:
      {
        int nPion_ch (0), nPion_0 (0);
        // count number of pions
        for( int i=1; i<6; i++ ) {
          if( flavs[i].Kfcode() == kf_pi_plus ) nPion_ch++;
          if( flavs[i].Kfcode() == kf_pi )      nPion_0++;
        }
        if( nPion_ch==3 && nPion_0==1 ) {
          hdme = new Tau_Four_Pion_3( nout, flavs );
        }
        if( nPion_ch==1 && nPion_0==3 ) {
          hdme = new Tau_Four_Pion_1( nout, flavs );
        }
	break;
      }
  }
}

// Select the corresponding ME of light pseudoscalar decay
void HD_ME_Selector::SelectLightPseudoScalarDecay(
						  int nout,Flavour * flavs,
						  HD_ME_Base *& hdme)
{
//   msg_Tracking()<<"HD_ME_Selector::SelectLightPseudoScalarDecay 1->"<<nout<<endl;
//   switch (nout) {
//   case 2:
//     if (flavs[1]==Flavour(kf_photon) &&
// 	flavs[2]==Flavour(kf_photon))
//       hdme = new P_2Gamma(nout,flavs);
//     break;
//   case 3:
//     if (flavs[1].IsLepton() &&
// 	flavs[2].IsLepton() && flavs[2].IsAnti() &&
// 	flavs[3]==Flavour(kf_photon)) {
//       hdme = new P_GammaFF(nout,flavs);
//     }
//     if (IsPseudoScalar(flavs[1]) &&
// 	flavs[2]==Flavour(kf_photon) &&
// 	flavs[3]==Flavour(kf_photon)) {
//       hdme = new P_P2Gamma(nout,flavs);
//     }
//     if (flavs[1]==Flavour(kf_pi_plus) &&
// 	flavs[2]==Flavour(kf_pi_plus).Bar() &&
// 	flavs[3]==Flavour(kf_photon)) {
//       hdme = new P_2PGamma(nout,flavs);
//       SetVector_For_2PS(nout,1,2,flavs,hdme);
//     }
//     if (IsPseudoScalar(flavs[1]) &&
// 	IsPseudoScalar(flavs[2]) &&
// 	IsPseudoScalar(flavs[3])) {
//       //hdme = new P_3P_Dalitz(nout,flavs);
//     }
//     break;
//   }
}

void HD_ME_Selector::SetVector_For_2PS(int nout,int PS1, int PS2,
				       Flavour * flavs,HD_ME_Base * hdme)
{
  double mass   = flavs[0].Mass();
  double charge = flavs[PS1].Charge()+flavs[PS2].Charge();
  if (charge==0.) {
    if (flavs[0].Mass()>Flavour(kf_rho_770).Mass())    
      hdme->AddVector(Flavour(kf_rho_770).Mass(),Flavour(kf_rho_770).Width());
    if (flavs[0].Mass()>Flavour(kf_omega_782).Mass())
      hdme->AddVector(Flavour(kf_omega_782).Mass(),Flavour(kf_omega_782).Width());
    if (flavs[0].Mass()>Flavour(kf_K_star_892).Mass())
      hdme->AddVector(Flavour(kf_K_star_892).Mass(),Flavour(kf_K_star_892).Width());
    if (flavs[0].Mass()>Flavour(kf_phi_1020).Mass())
      hdme->AddVector(Flavour(kf_phi_1020).Mass(),Flavour(kf_phi_1020).Width());
  }
  else if (charge!=0.) {
    if (flavs[0].Mass()>Flavour(kf_rho_770).Mass())
      hdme->AddVector(Flavour(kf_rho_770).Mass(),Flavour(kf_rho_770).Width());
    if (flavs[0].Mass()>Flavour(kf_K_star_892).Mass())
      hdme->AddVector(Flavour(kf_K_star_892).Mass(),Flavour(kf_K_star_892).Width());
  }
}

