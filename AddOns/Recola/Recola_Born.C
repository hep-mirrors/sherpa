#include "AddOns/Recola/Recola_Born.H"

#include "AddOns/Recola/Recola_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace Recola {

Recola_Born::Recola_Born(const Process_Info& pi,
                         const Flavour_Vector& flavs,
                         unsigned int recola_id, int amptype) :
  Tree_ME2_Base(pi, flavs), m_recola_id(recola_id), m_amptype(amptype)
  {
    m_symfac=pi.m_fi.FSSymmetryFactor();
    m_symfac*=pi.m_ii.ISSymmetryFactor();
    m_eventcount=0;
  }


  double Recola_Born::Calc(const Vec4D_Vector& momenta) {
    if(Recola_Interface::checkProcGeneration() == false) {
      Data_Reader reader(" ",";","#","=");
      int ewscheme = reader.GetValue<int>("EW_SCHEME",1);
      if(ewscheme == 3){
	use_gfermi_scheme_and_set_alpha_rcl(AlphaQED());
      }
      else if (ewscheme == 2){
//	use_alphaz_scheme_and_set_alpha_rcl(AlphaQED());
	use_alphaz_scheme_rcl(AlphaQED());
      }
      else if(ewscheme == 1){
//	use_alpha0_scheme_and_set_alpha_rcl(AlphaQED());
	use_alpha0_scheme_rcl(AlphaQED());
      }
      else{
	cout << "The EW scheme "<<ewscheme<<" is not available with the SHERPA+Recola interface. Valid options are:" << endl;
	cout <<"   1. alpha_QED(0)"<<endl;
	cout <<"   2. alpha_QED(M_Z)"<<endl;
	cout <<"   3. GFermi"<<endl;
	exit(1);
      }
      
      int nlight=0;
      set_mu_ir_rcl(sqrt(m_mur2));
      set_mu_uv_rcl(sqrt(m_mur2));
      double default_alphas(Recola_Interface::GetDefaultAlphaS());
      double default_scale(Recola_Interface::GetDefaultScale());
      double default_flavscheme(Recola_Interface::GetDefaultFlav());
      double fixed(Recola_Interface::GetFixed());
      double alpha_mat;
      if (default_flavscheme==0)
	default_flavscheme=fixed;
      if (default_flavscheme==16) default_flavscheme=-1;
      if (fixed>0 && fixed<10){
	nlight=fixed;
      }
      else{
	if (default_flavscheme>10){
	  nlight=Recola_Interface::PDFnf(m_mur2,default_flavscheme-10);
	}
	if (default_flavscheme==-1)
	  nlight=-1;
	if (default_flavscheme==-2 || default_flavscheme==0){
	  if (Flavour(kf_c).Mass()!=0)
	    nlight=3;
	  else if (Flavour(kf_b).Mass()!=0)
	    nlight=4;
	  else if (Flavour(kf_t).Mass()!=0)
	    nlight=5;
	  else {
	    msg_Out()<<"WARNING: 6 light flavours detected.\n";
	    nlight=6;
	  }
	}
      }
      if (nlight==0){
	msg_Error()<<METHOD<<"(): Cannot determine number of flavours\n";
      }
      if (nlight>6){
	msg_Error()<<METHOD<<"(): Too many light flavours: "<<nlight<<"\n   Max is 6\n";
      }
      
      if (default_scale==-1. || default_alphas==-1. || Recola_Interface::GetPDFDefault()==0){
	set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),nlight); 
	msg_Debugging() << "use AlphaQCD\n";
      }
      else{	
	set_alphas_rcl(default_alphas,sqrt(default_scale),nlight); 
	msg_Debugging() << "use default_alphas\n";
      }
      
      cout << "processes in Recola are being generated..." << endl;
      Recola_Interface::setProcGenerationTrue();
      generate_processes_rcl();
      cout << "process generation in Recola completed..." << endl;  
      get_alpha_rcl(alpha_mat);
      Recola_Interface::SetDefaultFlav(nlight);
    }
    double alpha(0);
    get_alpha_rcl(alpha);
    
    msg_Debugging() << "default_alphas = " << Recola_Interface::GetDefaultAlphaS() << ", sqrt(default_scale) = " << sqrt(Recola_Interface::GetDefaultScale()) << endl;
    msg_Debugging() << "AlphaQCD() = " << AlphaQCD() << ", sqrt(m_mur2) = " << sqrt(m_mur2) << endl;
    
    MyTiming* timing;
    if (msg_LevelIsDebugging()) {
      timing = new MyTiming();
      timing->Start();
    }
    
    msg_Debugging() << "m_born = " << m_born << ", m_res = " << m_res << ", m_amptype" << m_amptype << "\n";
    set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),Recola_Interface::GetDefaultFlav());
    if (m_amptype==12) Recola_Interface::EvaluateLoop(m_recola_id, momenta, m_born, m_res);
    
    if (msg_LevelIsDebugging()) {
      timing->Stop();
      PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
		  <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    }
    
    msg_Debugging() << "m_born = " << m_born << ", m_res = " << m_res << ", m_symfac = " << m_symfac << ", m_res.Finite()*m_symfac" << m_res.Finite()*m_symfac <<"\n";
    //need probably some kind of symmetry factor for the result
    Data_Reader reader(" ",";","#","=");
    m_eventcount+=1;
    if (reader.GetValue<int>("JUST_ONE_POINT",0)==1){
      std::cout<<std::setprecision(15)<<"\nMomenta are:  "<<momenta<<std::endl;
      std::cout<<"flavour:  "<<m_flavs<<std::endl;
      std::cout<< std::setprecision(15)<<"    and born is:  "<<m_res.Finite()<<std::endl;
      std::cout<< std::setprecision(15)<<"    and born*sym_factor is:  "<<m_res.Finite()*m_symfac<<std::endl;
      std::cout<<"mu:  "<<sqrt(m_mur2)<<std::endl;
      std::cout<<"AlphaQED():  "<<AlphaQED()<<std::endl;
      std::cout<<"AlphaQCD():  "<<AlphaQCD()<<std::endl;
    
    if (m_eventcount==201)
      exit(1);
    }
    
    // Recola returns ME2 including 1/symfac, but Calc is supposed to return it
    // without 1/symfac, thus multiplying with symfac here
    return m_res.Finite()*m_symfac;
  }
  
}

using namespace Recola;

DECLARE_TREEME2_GETTER(Recola_Born,"Recola_Born")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,Recola_Born>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Recola") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo && pi.m_fi.m_nloqcdtype!=nlo_type::lo && pi.m_fi.m_nloewtype!=nlo_type::real && pi.m_fi.m_nloqcdtype!=nlo_type::real) return NULL;
  msg_Debugging() << "Getter Function called, Born check passed\n";

  int id(0);
  id = Recola_Interface::RegisterProcess(pi, 12);
  if (id>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new Recola_Born(pi, flavs, id, 12);
  }

  return NULL;
}
