#include "Hdecay_Fortran_Interface.H"
#include "Message.H"

using namespace HDECAY;
using namespace MODEL;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;



extern "C" {
  void hdecayinter_(double &,int *,double *,double *,double *,double *);
  void hdecaysm_(double *,double *,double &);
}

Hdecay_Fortran_Interface::Hdecay_Fortran_Interface(AORGTOOLS::Data_Read * _dataread,
						   Model_Base * _model) :
  Spectrum_Generator_Base(_dataread,_model) {
}


void Hdecay_Fortran_Interface::Run(std::string _mode) {
  double * couplings = new double[3];
  couplings[0]       = 1./p_model->ScalarConstant(std::string("alpha_QED(0)"));
  couplings[1]       = p_model->ScalarConstant(std::string("GF"));
  couplings[2]       = p_model->ScalarConstant(std::string("alpha_S(MZ)"));

  double * bosons    = new double[4];
  bosons[0]          = p_model->ScalarConstant(std::string("MW"));
  bosons[1]          = p_model->ScalarConstant(std::string("MZ"));
  bosons[2]          = Flavour(kf::W).Width();
  bosons[3]          = Flavour(kf::Z).Width();

  double * yukawas   = new double[9];
  yukawas[0]         = 0.;
  yukawas[1]         = 0.;
  yukawas[2]         = p_dataread->GetValue<double>("MSTRANGE_EFF",0.48);
  yukawas[3]         = p_model->ScalarConstant(std::string("Yukawa_c"));
  yukawas[4]         = p_model->ScalarConstant(std::string("Yukawa_b"));
  yukawas[5]         = p_model->ScalarConstant(std::string("Yukawa_t"));
  yukawas[6]         = 0.;
  yukawas[7]         = p_model->ScalarConstant(std::string("Yukawa_mu"));
  yukawas[8]         = p_model->ScalarConstant(std::string("Yukawa_tau"));

  if (yukawas[2]==0. || yukawas[3]==0. || yukawas[4]==0. || 
      yukawas[5]==0. || yukawas[7]==0. || yukawas[8]==0.)
    msg.Error()<<"Potential error in Hdecay : "<<endl
	       <<"   The yukawas for muon, tau, s, c, b, t MUST be larger than 0. !!!!"<<endl
	       <<"   Expect nonsense results from Hdecay."<<endl;


  double * ckms      = new double[3];
  ckms[0]            = p_model->ComplexMatrixElement(std::string("CKM"),0,1).real();
  ckms[1]            = p_model->ComplexMatrixElement(std::string("CKM"),0,2).real();
  ckms[2]            = p_model->ComplexMatrixElement(std::string("CKM"),1,2).real();

  int * flags        = new int[4];
  flags[0]           = p_dataread->GetValue<int>("HIGGS_NNLO(M)",1); 
  flags[1]           = p_dataread->GetValue<int>("HIGGS_ON_SHELL_WZ",1);
  flags[2]           = p_dataread->GetValue<int>("HIGGS_MASS_SCHEME",0);
  flags[3]           = p_dataread->GetValue<int>("HIGGS_NF_LIGHT_Q_IN_G",5);

  if (_mode==std::string("SM")) {
    hmass = Flavour(kf::h).Mass(); 
    hdecayinter_(hmass,flags,couplings,bosons,yukawas,ckms);
  }
  CalculateEffectiveCouplings(_mode);

  delete [] flags;
  delete [] ckms;
  delete [] yukawas;
  delete [] bosons;
  delete [] couplings;
}

void Hdecay_Fortran_Interface::CalculateEffectiveCouplings(std::string _mode) {
  double * brff = new double[9];
  double * brVV = new double[5];
  double hwidth;
  hdecaysm_(brff,brVV,hwidth);

  Flavour  higgs = Flavour(kf::h);
  widths.push_back(hwidth);
  masses.push_back(higgs.Mass());
  particles.push_back(higgs);

  double   yukawa;
  double   vev       = p_model->ScalarConstant(std::string("vev"));
  double   ferm_pref = 8.*M_PI*sqr(vev)*hwidth/hmass;


  DecayChannel * decay;
  for (int i=1;i<7;++i) {
    decay = new DecayChannel(higgs);
    decay->AddDecayProduct(Flavour(kf::code(i)));
    decay->AddDecayProduct(Flavour(kf::code(i)).Bar());
    decayproducts.push_back(decay);
  }
  for (int i=11;i<16;i+=2) {
    decay = new DecayChannel(higgs);
    decay->AddDecayProduct(Flavour(kf::code(i)));
    decay->AddDecayProduct(Flavour(kf::code(i)).Bar());
    decayproducts.push_back(decay);
  }

  if (_mode==std::string("SM")) {
    for (int i=0;i<9;i++) {
      decayproducts[i]->SetDecayWidth(brff[i]*hwidth);
      msg.Debugging()<<"Decay width for "<<decayproducts[i]->Flin()<<" -> ";
      for (int j=0;j<decayproducts[i]->Nout();j++) 
	msg.Debugging()<<decayproducts[i]->Flout(j)<<" ";
      msg.Debugging()<<": "<<brff[i]<<" * "<<hwidth<<" GeV."<<endl;

      /*
	if (flav1.PSMass()<hmass/2. && p_model->Yukawa(i)>0.) {
	yukawa = sqrt(ferm_pref*brff[i]/(1.*sqrt(pow(1.-4.*sqr(flav1.PSMass()/hmass),3.))));
	if (flav1.IsQuark()) yukawa /= 3.;
	msg.Debugging()<<"Maybe set effective coupling for h -> "<<flav1<<" : "
		       <<yukawa<<" for mass = "<<flav1.PSMass()<<endl;
		       }
      */
    }
  }
  
  delete [] brff;
  delete [] brVV;
}
