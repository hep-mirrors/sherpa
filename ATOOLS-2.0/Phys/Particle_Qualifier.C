#include "Particle_Qualifier.H"

using namespace ATOOLS;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Particle_Qualifier_Base
#define PARAMETER_TYPE std::string
#include "Getter_Function.C"

template <class Class>
Particle_Qualifier_Base *const GetQualifier(const std::string &parameter)
{									
  return new Class();
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Particle_Qualifier_Base *const					\
  NAME::operator()(const std::string &parameter) const			\
  { return GetQualifier<CLASS>(parameter); }

#define DEFINE_PRINT_METHOD(NAME,PRINT)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<PRINT; }

#define DEFINE_QUALIFIER_GETTER(CLASS,NAME,TAG,PRINT)			\
  DECLARE_GETTER(NAME,TAG,Particle_Qualifier_Base,std::string);		\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME,PRINT)

#include "Message.H"

void Particle_Qualifier_Base::ShowQualifiers(const int mode)
{
  if (!msg.LevelIsInfo() || mode==0) return;
  msg.Out()<<"Particle_Qualifier_Base::ShowQualifiers(): {\n\n";
  Particle_Qualifier_Getter::PrintGetterInfo(msg.Out(),20);
  msg.Out()<<"\n}"<<std::endl;
}

DEFINE_QUALIFIER_GETTER(Is_Charged_Hadron,Is_Charged_Hadron_Getter,
			"1","charged hadron");
DEFINE_QUALIFIER_GETTER(Is_Charged_Hadron,Is_Charged_Hadron_Getter_,
			"ChargedHadron","charged hadron");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Hadron,Is_Neutral_Hadron_Getter,
			"2","neutral hadron");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Hadron,Is_Neutral_Hadron_Getter_,
			"NeutralHadron","neutral hadron");
DEFINE_QUALIFIER_GETTER(Is_Hadron,Is_Hadron_Getter,
			"3","hadron");
DEFINE_QUALIFIER_GETTER(Is_Hadron,Is_Hadron_Getter_,
			"Hadron","hadron");
DEFINE_QUALIFIER_GETTER(Is_Charged,Is_Charged_Getter,
			"4","charged");
DEFINE_QUALIFIER_GETTER(Is_Charged,Is_Charged_Getter_,
			"Charged","charged");
DEFINE_QUALIFIER_GETTER(Is_Neutral,Is_Neutral_Getter,
			"5","neutral");
DEFINE_QUALIFIER_GETTER(Is_Neutral,Is_Neutral_Getter_,
			"Neutral","neutral");
DEFINE_QUALIFIER_GETTER(Is_Charged_Pion,Is_Charged_Pion_Getter,
			"11","charged pion");
DEFINE_QUALIFIER_GETTER(Is_Charged_Pion,Is_Charged_Pion_Getter_,
			"ChargedPion","charged pion");
DEFINE_QUALIFIER_GETTER(Is_Charged_Kaon,Is_Charged_Kaon_Getter,
			"12","charged kaon");
DEFINE_QUALIFIER_GETTER(Is_Charged_Kaon,Is_Charged_Kaon_Getter_,
			"ChargedKaon","charged kaon");
DEFINE_QUALIFIER_GETTER(Is_Proton_Antiproton,Is_Proton_Antiproton_Getter,
			"13","proton antiproton");
DEFINE_QUALIFIER_GETTER(Is_Proton_Antiproton,Is_Proton_Antiproton_Getter_,
			"ProtonAntiproton","proton antiproton");
DEFINE_QUALIFIER_GETTER(Is_Parton,Is_Parton_Getter,
			"21","parton");
DEFINE_QUALIFIER_GETTER(Is_Parton,Is_Parton_Getter_,
			"Parton","parton");
DEFINE_QUALIFIER_GETTER(Is_There,Is_There_Getter,
			"22","there");
DEFINE_QUALIFIER_GETTER(Is_There,Is_There_Getter_,
			"There","there");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Pion,Is_Neutral_Pion_Getter,
			"101","neutral pion");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Pion,Is_Neutral_Pion_Getter_,
			"NeutralPion","neutral pion");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Kaon,Is_Neutral_Kaon_Getter,
			"102","neutral kaon");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Kaon,Is_Neutral_Kaon_Getter_,
			"NeutralKaon","neutral kaon");
DEFINE_QUALIFIER_GETTER(Is_Charged_KStar,Is_Charged_KStar_Getter,
			"103","charged kstar");
DEFINE_QUALIFIER_GETTER(Is_Charged_KStar,Is_Charged_KStar_Getter_,
			"ChargedKStar","charged kstar");
DEFINE_QUALIFIER_GETTER(Is_Neutral_KStar,Is_Neutral_KStar_Getter,
			"104","charged kstar");
DEFINE_QUALIFIER_GETTER(Is_Neutral_KStar,Is_Neutral_KStar_Getter_,
			"NeutralKStar","charged kstar");
DEFINE_QUALIFIER_GETTER(Is_Eta,Is_Eta_Getter,
			"105","eta");
DEFINE_QUALIFIER_GETTER(Is_Eta,Is_Eta_Getter_,
			"Eta","eta");
DEFINE_QUALIFIER_GETTER(Is_Rho0,Is_Rho0_Getter,
			"106","rho0");
DEFINE_QUALIFIER_GETTER(Is_Rho0,Is_Rho0_Getter_,
			"Rho0","rho0");
DEFINE_QUALIFIER_GETTER(Is_Omega,Is_Omega_Getter,
			"107","omega");
DEFINE_QUALIFIER_GETTER(Is_Omega,Is_Omega_Getter_,
			"Omega","omega");
DEFINE_QUALIFIER_GETTER(Is_EtaPrime,Is_EtaPrime_Getter,
			"108","eta prime");
DEFINE_QUALIFIER_GETTER(Is_EtaPrime,Is_EtaPrime_Getter_,
			"EtaPrime","eta prime");
DEFINE_QUALIFIER_GETTER(Is_Phi,Is_Phi_Getter,
			"109","phi");
DEFINE_QUALIFIER_GETTER(Is_Phi,Is_Phi_Getter_,
			"Phi","phi");
DEFINE_QUALIFIER_GETTER(Is_Lambda,Is_Lambda_Getter,
			"110","lambda");
DEFINE_QUALIFIER_GETTER(Is_Lambda,Is_Lambda_Getter_,
			"Lambda","lambda");
DEFINE_QUALIFIER_GETTER(Is_Charged_Sigma,Is_Charged_Sigma_Getter,
			"111","charged sigma");
DEFINE_QUALIFIER_GETTER(Is_Charged_Sigma,Is_Charged_Sigma_Getter_,
			"ChargedSigma","charged sigma");
DEFINE_QUALIFIER_GETTER(Is_Charged_Xi,Is_Charged_Xi_Getter,
			"112","charged xi");
DEFINE_QUALIFIER_GETTER(Is_Charged_Xi,Is_Charged_Xi_Getter_,
			"ChargedXi","charged xi");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Xi,Is_Neutral_Xi_Getter,
			"113","neutral xi");
DEFINE_QUALIFIER_GETTER(Is_Neutral_Xi,Is_Neutral_Xi_Getter_,
			"NeutralXi","neutral xi");
