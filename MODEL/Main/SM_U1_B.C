#include "MODEL/Main/SM_U1_B.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Effective_Higgs_Coupling.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(SM_U1_B,"SM+U1_B",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,SM_U1_B>::
operator()(const Model_Arguments &args) const
{
  return new SM_U1_B(args.m_path,args.m_file,args.m_elementary);
}

void Getter<Model_Base,Model_Arguments,SM_U1_B>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + U(1)_B \n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- M_Z', g'_1\n"
     <<std::setw(width+7)<<" "<<"- MIX = mixing matrix (1=diag, 0=demo, -1=user)\n"
     <<std::setw(width+7)<<" "<<"- with -1=user not implemented yet\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


SM_U1_B::SM_U1_B(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    CustomContainerInit();
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool SM_U1_B::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model plus U(1)_B from "
	      <<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+U1_B");

  p_sm->ModelInit(isr);
  p_numbers          = p_sm->ExtractScalarNumbers();
  p_constants        = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions        = p_sm->ExtractScalarFunctions();
  p_matrices         = p_sm->ExtractComplexMatrices();

  delete p_sm;

  FillSpectrum(isr);

  if (!SanityChecks()) {
    msg_Error()<<"Potential Error in "<<METHOD<<":"<<endl
	       <<"   Sanity checks not passed."<<endl
	       <<"   Continue and hope for the best."<<endl;
  }
  return true;
}

SM_U1_B::~SM_U1_B() 
{ }

void SM_U1_B::ParticleInit() {
  //add new gauge boson
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,
  //                                                                idname,tex_name
  s_kftable[kf_Z0_2] = 
    new Particle_Info(kf_Z0_2,1000.,10.0,0,0,0,2,-1,1,1,1,"Z'","Z'");
  ReadParticleData();
}

void SM_U1_B::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  p_dataread->RereadInFile();
  p_constants->insert(make_pair(string("M_Z'"),    
				p_dataread->GetValue<double>("M_Z'",-1.)));
  p_constants->insert(make_pair(string("g'_1"),    
				p_dataread->GetValue<double>("g'_1",0.)));

  double MZprime(ScalarConstant(string("M_Z'"))),GZprime(0.),inc;
  double g1prime(sqr(ScalarConstant(std::string("g'_1"))));
  FixMix();
  msg_Out()<<"Calculate width for Z' with mass = "<<MZprime<<" GeV.\n";

  for (short int i=1;i<=6;i++) {
    Flavour quark = Flavour((kf_code)(i));
    double  massq = quark.HadMass(); 
    if (!quark.IsOn() || !quark.Strong()) continue;
    if (i%2) {
      if (2.*massq>MZprime) continue;
      GZprime += inc = 
	3.*MZprime*g1prime/(24.*M_PI)*pow(1.-sqr(2.*massq/MZprime),3./2.);
      msg_Out()<<"   Add Z' --> "<<quark<<"+"<<quark.Bar()<<" ("<<i<<"), "
	       <<" add "<<inc<<" GeV.\n";
    }
    else {
      double M2 = MZprime*MZprime, mq2 = massq*massq;
      for (short int j=2;j<=7;j+=2) {
	Flavour anti = Flavour((kf_code)(j)).Bar();
	if (!anti.IsOn() || !anti.Strong()) continue;
	double massa  = anti.HadMass(), ma2 = massa*massa; 
	if (massq+massa>MZprime) continue;
	Complex mix   = ComplexMatrixElement(std::string("UpMix"),i/2-1,j/2-1);
	double absmix = ATOOLS::sqr(mix.real())+ATOOLS::sqr(mix.imag());
	if (ATOOLS::IsZero(absmix)) continue;
	GZprime += inc = 
	  3.*g1prime*absmix/(24.*M_PI*M2)*
	  ((2.*M2-mq2-ma2)/2.-3.*massq*massa-ATOOLS::sqr(ma2-mq2)/(2.*M2)) *
	  sqrt(ATOOLS::sqr(M2-mq2-ma2)-4.*ma2*mq2)/(2.*MZprime);
	msg_Out()<<"   Add Z' --> "<<quark<<"+"<<anti<<" ("<<i<<"), "
		 <<" add "<<inc<<" GeV.\n";
      }
    }
  }

  Flavour flav;
  flav = Flavour(kf_Z0_2);
  flav.SetMass(ScalarConstant(string("M_Z'")));
  flav.SetHadMass(MZprime);
  flav.SetMassOn(true);
  flav.SetWidth(GZprime);

  msg_Out()<<METHOD<<" initializes Z' boson with \n"
	   <<"    mass = "<<MZprime<<" GeV and width = "<<GZprime<<" GeV\n"
	   <<"    for g'_1 = "<<ScalarConstant(std::string("g'_1"))<<".\n";
}

void SM_U1_B::FixMix() {
  CMatrix Mix(3);
  switch (p_dataread->GetValue<int>("MIX",1)) {
  case 1:
    for (int i=0;i<3;i++) {
      for (int j=i;j<3;j++) Mix[i][j] = Mix[j][i] = Complex(0.,0.);
      Mix[i][i] = Complex(1.,0.);
    }
    break;
  case 0: {
    double  Theta(2.*M_PI/3.);
    Complex diag(Complex(1./sqrt(3.),0.));
    Complex offdiag(Complex(cos(Theta)/sqrt(3.),sin(Theta)/sqrt(3.)));
    Mix[0][0] = Mix[1][1] = Mix[2][2] = diag;
    Mix[0][1] = Mix[0][2] = Mix[1][2] = Mix[1][0] = Mix[2][0] = Mix[2][1] = offdiag;
    break;
  }
  case -1:
  default:
    msg_Error()<<"Error in "<<METHOD<<"(mix = "<<p_dataread->GetValue<int>("MIX",1)<<".\n"
	       <<"   Option not implemented yet.  Will exit the run.\n";
    exit(1);
  }
  p_matrices->insert(std::make_pair(std::string("UpMix"),Mix));

  int output=p_dataread->GetValue<int>("OUTPUT_MIXING",0);
  msg_Out()<<"Output = "<<output<<".\n";
  if (output>0) {
    unsigned int os(25);
    msg_Out()<<"quark mixing matrix:\n";
    for (int i=0;i<3;++i) {
      for (int j=0;j<3;++j) {
	msg_Out()<<std::setw(os)<<Mix[i][j];
      }
      msg_Out()<<"\n";
    }
    msg_Out()<<"Unitarity check:\n";
    CMatrix Mixconj = Mix.Conjugate();
    for (int i=0;i<3;++i) {
      for (int j=0;j<3;++j) {
	Complex entry=Complex(0.,0.);
	for (int k=0;k<3;++k) entry += Mix[i][k]*Mixconj[k][j];
	if (ATOOLS::dabs(entry.real())<1.e-12) entry = Complex(0., entry.imag());
	if (ATOOLS::dabs(entry.imag())<1.e-12) entry = Complex(entry.real(), 0.);
	if (ATOOLS::dabs(1.-entry.real())<1.e-12) entry = Complex(1., entry.imag());
	if (ATOOLS::dabs(1.-entry.imag())<1.e-12) entry = Complex(entry.real(), 1.);
	msg_Out()<<std::setw(os)<<entry;
      }
      msg_Out()<<"\n";
    }	
  }
}

bool SM_U1_B::SanityChecks() {
  if (ScalarConstant(string("g'_1"))<1.e-6 ||
      ScalarConstant(string("M_Z'"))<0.) return false;
  return true;
}



