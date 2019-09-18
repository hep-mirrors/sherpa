#include "Recola_Virtual.H"

#include "AddOns/Recola/Recola_Interface.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace Recola;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Recola_Virtual::Recola_Virtual(const Process_Info& pi,
			      const Flavour_Vector& flavs,
			      unsigned int recola_id) :
  Virtual_ME2_Base(pi, flavs), m_recola_id(recola_id),
  m_ewscheme(ToType<int>(rpa->gen.Variable("EW_SCHEME"))),
  m_voqcd(pi.m_maxcpl[0]),
  m_boqcd(pi.m_maxcpl[0]-(pi.m_fi.m_nloqcdtype==nlo_type::loop))
{
  m_fixedIRscale=true;

  m_IRscale=Recola_Interface::IRScale();
  m_UVscale=Recola_Interface::UVScale();

  // init associated contribs
  size_t n(0);
  if (pi.m_fi.m_asscontribs&asscontrib::EW) {
    ++n;
    if (pi.m_fi.m_asscontribs&asscontrib::LO1) {
      ++n;
      if (pi.m_fi.m_asscontribs&asscontrib::LO2) {
        ++n;
        if (pi.m_fi.m_asscontribs&asscontrib::LO3) {
          ++n;
        }
      }
    }
  }
  m_asscontribs.resize(n);

  if (m_asscontribs.size()>0 && m_voqcd!=m_boqcd+1)
    THROW(fatal_error,"Associated contribs only implemented for NLO QCD.");
}

void Recola_Virtual::Calc(const Vec4D_Vector& momenta)
{
  DEBUG_FUNC("");
  if (!Recola_Interface::checkProcGeneration()) {
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
    set_mu_ir_rcl(m_IRscale);
    set_mu_uv_rcl(m_UVscale);
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

    // Recola_Interface::SetDefaultFlav(nlight);

    double default_alphaQCD=Recola_Interface::GetDefaultAlphaQCD();
    double default_scale=Recola_Interface::GetDefaultScale();
    set_alphas_rcl(default_alphaQCD,sqrt(default_scale),nlight);
    msg_Debugging()<<"use AlphaQCD="<<AlphaQCD()
                   <<",  sqrt(m_mur2)="<<sqrt(m_mur2)<<std::endl;

    msg_Out()<<"processes in Recola are being generated..."<<std::endl;
    generate_processes_rcl();
    Recola_Interface::setProcGenerationTrue();
    msg_Out()<<"process generation in Recola completed..."<<std::endl;
  }

  // calculate
  m_res*=0.; m_born=0.;
  for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]=0.;

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
    msg_Out()<<"virtual correction, ";
  }
  set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),Recola_Interface::GetDefaultFlav());
  Recola_Interface::EvaluateProcess(m_recola_id, momenta, m_voqcd, m_boqcd,
                                    m_res,m_born,m_asscontribs);

  if (msg_LevelIsDebugging()) {
    timing->Stop();
    msg_Out()<<"me="<<m_res.Finite()<<" "<<m_born<<std::endl;
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
               <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
  }

  // factor which by Sherpa convention has to be divided out at this stage
  if (m_born==0) m_born=1.;
  double factor=m_born*AlphaQCD()/2.0/M_PI;
  m_res.Finite()/=factor;
  m_res.IR()/=factor;
  m_res.IR2()/=factor;
  for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]/=factor;
  msg_Debugging()<<"V/B="<<m_res.Finite()<<std::endl;
  for (size_t i(0);i<m_asscontribs.size();++i)
    msg_Debugging()<<"ASS/B="<<m_asscontribs[i]<<std::endl;
}



DECLARE_VIRTUALME2_GETTER(Recola_Virtual,"Recola_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,Recola_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Recola") return NULL;

  if (pi.m_fi.m_nloqcdtype!=nlo_type::loop) return NULL;

  int procIndex=Recola_Interface::RegisterProcess(pi, amptype::treeloop);
  if (procIndex>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new Recola_Virtual(pi, flavs, procIndex);
  }
  else {
    return NULL;
  }
}
