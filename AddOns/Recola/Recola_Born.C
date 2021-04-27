#include "AddOns/Recola/Recola_Born.H"

#include "AddOns/Recola/Recola_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace Recola;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Recola_Born::Recola_Born(const Process_Info& pi,
                         const Flavour_Vector& flavs,
                         unsigned int recola_id) :
  Tree_ME2_Base(pi, flavs), m_recola_id(recola_id),
  m_ewscheme(ToType<int>(rpa->gen.Variable("EW_SCHEME"))),
  m_oqcd(pi.m_maxcpl[0])
{
  m_symfac=pi.m_fi.FSSymmetryFactor();
  m_symfac*=pi.m_ii.ISSymmetryFactor();
}


double Recola_Born::Calc(const Vec4D_Vector& momenta) {
  if (Recola_Interface::checkProcGeneration() == false) {
    // use different ew scheme in different interface functions
    switch (m_ewscheme) {
    case 3:
      use_gfermi_scheme_and_set_alpha_rcl(AlphaQED());
      break;
    case 2:
      //use_alphaz_scheme_and_set_alpha_rcl(AlphaQED());
      use_alphaz_scheme_rcl(AlphaQED());
      break;
    case 1:
      //use_alpha0_scheme_and_set_alpha_rcl(AlphaQED());
      use_alpha0_scheme_rcl(AlphaQED());
      break;
    default:
      msg_Out()<<"The EW scheme "<<m_ewscheme<<" is not available with the "
               <<"Sherpa+Recola interface. Valid options are:"<<std::endl
               <<"  1. alpha_QED(0)"<<std::endl
               <<"  2. alpha_QED(M_Z)"<<std::endl
               <<"  3. Gmu"<<std::endl;
      THROW(not_implemented,"Electroweak scheme not implemented in Sherpa+Recola.");
    }

    int nlight=0;
    set_mu_ir_rcl(Recola_Interface::IRScale());
    set_mu_uv_rcl(Recola_Interface::UVScale());
    int fixed(Recola_Interface::GetFixedFlav());
    if (Recola_Interface::GetDefaultFlav()==0) fixed=5;

    double alpha_mat;
    int default_flavscheme(fixed);
    if (default_flavscheme==16) default_flavscheme=-1;
    if (fixed>0 && fixed<10) {
      nlight=fixed;
    }
    else {
      if (default_flavscheme>10) {
        nlight=Recola_Interface::PDFnf(m_mur2,default_flavscheme-10);
      }
      if (default_flavscheme==-1)
        nlight=-1;
      if (default_flavscheme==-2 || default_flavscheme==0) {
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

    Recola_Interface::SetDefaultFlav(nlight);
    double default_alphaQCD=Recola_Interface::GetDefaultAlphaQCD();
    double default_scale=Recola_Interface::GetDefaultScale();
    set_alphas_rcl(default_alphaQCD,sqrt(default_scale),nlight);
    msg_Debugging() << "use AlphaQCD\n";


    msg_Out() << "processes in Recola are being generated..." << endl;
    Recola_Interface::setProcGenerationTrue();
    generate_processes_rcl();
    msg_Out() << "process generation in Recola completed..." << endl;
    get_alpha_rcl(alpha_mat);
    Recola_Interface::SetDefaultFlav(nlight);
  }
  double alpha(0);
  get_alpha_rcl(alpha);

  msg_Debugging() << "default_alphas = " << Recola_Interface::GetDefaultAlphaQCD() << ", sqrt(default_scale) = " << sqrt(Recola_Interface::GetDefaultScale()) << endl;
  msg_Debugging() << "AlphaQCD() = " << AlphaQCD() << ", sqrt(m_mur2) = " << sqrt(m_mur2) << endl;

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
  }

  METOOLS::DivArrD vres; std::vector<double> ass;
  set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),Recola_Interface::GetDefaultFlav());
  Recola_Interface::EvaluateProcess(m_recola_id, momenta, 0, m_oqcd,
                                    vres, m_born, ass, false);

  if (msg_LevelIsDebugging()) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
                <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
  }

  // Recola returns ME2 including 1/symfac, but Calc is supposed to return it
  // without 1/symfac, thus multiplying with symfac here
  return m_born*m_symfac;
}




DECLARE_TREEME2_GETTER(Recola_Born,"Recola_Born")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,Recola_Born>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Recola") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo && pi.m_fi.m_nloqcdtype!=nlo_type::lo && pi.m_fi.m_nloewtype!=nlo_type::real && pi.m_fi.m_nloqcdtype!=nlo_type::real) return NULL;
  msg_Debugging() << "Getter Function called, Born check passed\n";

  int id(0);
  id = Recola_Interface::RegisterProcess(pi, amptype::treetree);
  if (id>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new Recola_Born(pi, flavs, id);
  }

  return NULL;
}
