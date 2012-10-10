#include "LH_OLE_Communicator.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Poincare.H"

using namespace OLE;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace OLE {
extern "C" void OLP_Start(const char * filename, int* success);
extern "C" void OLP_EvalSubProcess(int,double*,double,double*,double*);
extern "C" void OLP_Option(const char * assignment, int* success);

  class LH_OLE_Interface : public Virtual_ME2_Base {
    size_t m_pn;
    bool m_active, m_needcmsboost;
    int m_OLE_id;
    double* p_momenta;
    double p_result[4];
    static int s_oleinit;
    int m_nf;
  public:
    LH_OLE_Interface(const Process_Info& pi,const Flavour_Vector& flavs,bool active);
    ~LH_OLE_Interface() {
      if (p_momenta) delete[] p_momenta;
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    bool SetColours(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);

  };
}

int LH_OLE_Interface::s_oleinit=0;

LH_OLE_Interface::LH_OLE_Interface(const Process_Info& pi, 
                                   const Flavour_Vector& flavs,bool active) :
  Virtual_ME2_Base(pi, flavs), m_pn(flavs.size()), m_active(active), 
  m_needcmsboost(false), m_OLE_id(-1), p_momenta(NULL), m_nf(0)
{
  if (!m_active) return;
  m_mode = 0;
  m_drmode = 1;
  for (size_t i(0);i<Flavour(kf_quark).Size();++i) {
    if (Flavour(kf_quark)[i].Strong()) ++m_nf;
  }
  if (m_nf%2==1) THROW(fatal_error,"Uneven number of quark and anti-quark flavours.");
  m_nf/=2;
  msg_Debugging()<<METHOD<<"(): nf = "<<m_nf<<std::endl;

  bool contract(0);
  string orderfn("OLE_order.lh");
  string contractfn("OLE_contract.lh");
  string fname("");
  Data_Reader reader(" ",";","!","=");
  reader.SetInputPath(rpa->GetPath());
  if (reader.ReadFromFile(fname,"LHOLE_ORDERFILE")) {
    orderfn=fname;
  }
  if (reader.ReadFromFile(fname,"LHOLE_CONTRACTFILE")) {
    contractfn=fname;
  }
  m_needcmsboost=reader.GetValue<int>("LHOLE_CMSBOOST",0);
  ifstream ifile;
  ifile.open(contractfn.c_str());
  if (ifile) {
    contract=1;
    fname=contractfn;
    ifile.close();
  }
  else fname=orderfn;

  LH_OLE_Communicator lhfile(fname);
  if (!contract) {
    if (lhfile.FileStatus()==0) {
      lhfile.AddParameter("MatrixElementSquareType CHsummed");
      lhfile.AddParameter("CorrectionType          QCD");
      lhfile.AddParameter("IRregularisation        DRED");
      lhfile.AddParameter("AlphasPower             "+ToString(pi.m_oqcd-1));
      lhfile.AddParameter("AlphaPower              "+ToString(pi.m_oew));
      lhfile.AddParameter("OperationMode           CouplingsStrippedOff");
      std::string widthscheme("FixedWidthScheme");
      if (MODEL::s_model->ScalarNumber(std::string("WidthScheme")))
        widthscheme=std::string("ComplexMassScheme");
      lhfile.AddParameter("ResonanceTreatment      "+widthscheme);
      lhfile.AddParameter("EWRenormalisationScheme alphaMZ");
      lhfile.AddParameter("");
      lhfile.AddParameter("Z_mass                  "+ToString(Flavour(kf_Z).Mass()));
      lhfile.AddParameter("Z_width                 "+ToString(Flavour(kf_Z).Width()));
      lhfile.AddParameter("W_mass                  "+ToString(Flavour(kf_Wplus).Mass()));
      lhfile.AddParameter("W_width                 "+ToString(Flavour(kf_Wplus).Width()));
      std::string sin_th_2(ToString(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))));
      if (MODEL::s_model->ScalarNumber(std::string("WidthScheme")))
        sin_th_2=ToString(ToString(MODEL::s_model->ScalarConstant(std::string("csin2_thetaW"))));
      lhfile.AddParameter("sin_th_2                "+sin_th_2);
      lhfile.AddParameter("H_mass                  "+ToString(Flavour(kf_h0).Mass()));
      lhfile.AddParameter("H_width                 "+ToString(Flavour(kf_h0).Width()));
      lhfile.AddParameter("top_mass                "+ToString(Flavour(kf_t).Mass()));
      lhfile.AddParameter("top_width               "+ToString(Flavour(kf_t).Width()));
      lhfile.AddParameter("bottom_mass             "+ToString(Flavour(kf_b).Mass()));
      lhfile.AddParameter("bottom_width            "+ToString(Flavour(kf_b).Width()));
      lhfile.AddParameter("");
      lhfile.AddParameter("# process list");
    }
    if(lhfile.CheckProcess(2,m_pn-2,flavs)==-1) {
      lhfile.AddProcess(2,m_pn-2,flavs);
    }
    return;
  }

  if (lhfile.CheckParameterStatus()!=1) {
    THROW(fatal_error,"Bad OLE parameter");
  }

  int pstatus=lhfile.CheckProcess(2,m_pn-2,flavs);
  switch (pstatus) {
  case -2: 
  case 0:
    THROW(fatal_error,"OLE did not return valid process ID.");
  case -1: 
    THROW(fatal_error,"Process not found in contract file.");
  default:
    if (pstatus!=1) msg_Info()<<"Found "<<pstatus<<" subprocesses. Cannot "
			<<"handle this yet, only first ID will be used!"<<endl;
    m_OLE_id=lhfile.GetID(2,m_pn-2,flavs,0);
  }

  if (s_oleinit==0) {
    int check(0);
    // -- GoSam specific: --
    // Weak Gauge Bosons + Higgs
    string mZ_string("mZ="+ToString(Flavour(kf_Z).Mass()));
    string wZ_string("wZ="+ToString(Flavour(kf_Z).Width()));
    string mW_string("mW="+ToString(Flavour(kf_Wplus).Mass()));
    string wW_string("wW="+ToString(Flavour(kf_Wplus).Width()));
    string mH_string("mH="+ToString(Flavour(kf_h0).Mass()));
    OLE::OLP_Option(mZ_string.c_str(),&check);
    OLE::OLP_Option(wZ_string.c_str(),&check);
    OLE::OLP_Option(mW_string.c_str(),&check);
    OLE::OLP_Option(wW_string.c_str(),&check);
    OLE::OLP_Option(mH_string.c_str(),&check);
    // Quarks
    string mB_string("mB="+ToString(Flavour(kf_b).Mass()));
    string wB_string("wB="+ToString(Flavour(kf_b).Width()));
    string mT_string("mT="+ToString(Flavour(kf_t).Mass()));
    string wT_string("wT="+ToString(Flavour(kf_t).Width()));
    OLE::OLP_Option(mB_string.c_str(),&check);
    OLE::OLP_Option(wB_string.c_str(),&check);
    OLE::OLP_Option(mT_string.c_str(),&check);
    OLE::OLP_Option(wT_string.c_str(),&check);
    // -- GoSam specific end --
    OLE::OLP_Start(fname.c_str(),&check);
    if (check != 1) THROW(fatal_error,"OLP initialisation failed");
    s_oleinit=1;
  }
  p_momenta = new double[m_pn*5];
  for (size_t i=0;i<m_pn;i++) p_momenta[4+i*5]=flavs[i].Mass();

  for (size_t i=0;i<3;i++) p_result[i]=0.;
  p_result[3]=1.;
}

void LH_OLE_Interface::Calc(const Vec4D_Vector& pp) {
  if (!m_active) return;
  if (m_OLE_id<0) return;

  Vec4D_Vector momenta(pp);
  if (m_needcmsboost) {
    Poincare cms(momenta[0]+momenta[1]);
    for (size_t i(0);i<momenta.size();++i) cms.Boost(momenta[i]);
  }
  for (size_t i=0;i<m_pn;i++) {
    p_momenta[0+i*5]=momenta[i][0];
    p_momenta[1+i*5]=momenta[i][1];
    p_momenta[2+i*5]=momenta[i][2];
    p_momenta[3+i*5]=momenta[i][3];
  }
  double param(1.);

  OLE::OLP_EvalSubProcess(m_OLE_id,p_momenta,sqrt(m_mur2),&param,p_result);
  // correct normalization:
  double one_over_2pi = 0.15915494309189533577;
  for (size_t i=0;i<3;i++) p_result[i]/=one_over_2pi;
  // finite
  m_res.Finite() = p_result[2]/p_result[3];
  // 1/epsIR
  m_res.IR()     = p_result[1]/p_result[3];
  // 1/epsIR2
  m_res.IR2()    = p_result[0]/p_result[3];
}

bool LH_OLE_Interface::SetColours(const ATOOLS::Vec4D_Vector& momenta) {
  return true;
}

double LH_OLE_Interface::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{   
  return 4.*M_PI;
}

DECLARE_VIRTUALME2_GETTER(LH_OLE_Interface_Getter,"LH_OLE_Interface")
Virtual_ME2_Base *LH_OLE_Interface_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="LHOLE") return NULL;
  msg_Info()<<"Les Houches One-Loop Generator called.\n";
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    return new LH_OLE_Interface(pi, fl, true);
  }
  else if (pi.m_fi.m_nloqcdtype&nlo_type::vsub) {
    return new LH_OLE_Interface(pi, fl, false);
  }
  msg_Info()<<"Les Houches One-Loop Generator could not provide one-loop \n"
           <<"matrix element for "
           <<PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_ii)<<".\n";
  return NULL;
}
