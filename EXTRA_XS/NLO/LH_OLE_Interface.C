#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Exception.H"
#include "ATOOLS/Org/LH_OLE_Communicator.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace OLE {
  void Init(const char * filename) {}
  void EvalSubprocess(int,double*,double,double,double,double*) {}
}

namespace EXTRAXS {
  class LH_OLE_Interface : public Virtual_ME2_Base {
    double m_bf;
    size_t m_pn;
    bool m_active;
    int m_OLE_id;
    double* p_momenta;
    double p_result[4];
    double m_as,m_aqed;
    static int s_bhinit;
    double m_cpl;
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

int EXTRAXS::LH_OLE_Interface::s_bhinit=0;

LH_OLE_Interface::LH_OLE_Interface(const Process_Info& pi, const Flavour_Vector& flavs,bool active) :
  Virtual_ME2_Base(pi, flavs), m_OLE_id(-1), p_momenta(0)
{
  m_active=active;
  if (!m_active) return;
  m_needsborn = true;
  m_cpl = MODEL::s_model->ScalarFunction(std::string("alpha_S"),
					 sqr(rpa.gen.Ecms()));
  m_cpl /= 2.*M_PI;
  Flavour hfl(kf_quark);
  m_nf = hfl.Size()/2;

  m_pn=flavs.size();

  bool contract(0);
  string fname("OLE_order.lh");
  ifstream ifile;
  ifile.open(string("OLE_contract.lh").c_str());
  if (ifile) {
    contract=1;
    fname=string("OLE_contract.lh");
    ifile.close();
  }

  LH_OLE_Communicator lhfile(fname);
  if (!contract) {
    if (lhfile.FileStatus()==0) {
      lhfile.AddParameter("AmplitudeType CH_SUMMED");
      lhfile.AddParameter("Z_mass "+ToString(Flavour(kf_Z).Mass()));
      lhfile.AddParameter("Z_width "+ToString(Flavour(kf_Z).Width()));
      lhfile.AddParameter("W_mass "+ToString(Flavour(kf_Wplus).Mass()));
      lhfile.AddParameter("W_width "+ToString(Flavour(kf_Wplus).Width()));
      double sin_th_2=MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"));
      lhfile.AddParameter("sin_th_2 "+ToString(sin_th_2));
      lhfile.AddParameter("sin_2th "+ToString(sin(2.*asin(sqrt(sin_th_2)))));
      lhfile.AddParameter("");
    }
    if(lhfile.CheckProcess(2,m_pn-2,flavs)==-1) {
      lhfile.AddProcess(2,m_pn-2,flavs);
      m_newlibs=1;
    }
    return;
  }

  if (lhfile.CheckParameterStatus()!=1) {
    THROW(fatal_error,"Bad OLE parameter");
  }

  int pstatus=lhfile.CheckProcess(2,m_pn-2,flavs);
  switch (pstatus) {
  case -1: cout<<"Error: Process not found in contract file!"<<endl;
    abort();
  case 0: cout<<"Error: No OLE info"<<endl;
    abort();
  default:
    if (pstatus!=1) cout<<endl<<"Found "<<pstatus<<" subprocesses. Cannot handle this yet,"
			<<" only first ID will be used!"<<endl;
    m_OLE_id=lhfile.GetID(2,m_pn-2,flavs,0);
  }

  if (s_bhinit==0) {
    OLE::Init(fname.c_str());
    s_bhinit=1;
  }
  p_momenta = new double[m_pn*5];
  for (size_t i=0;i<m_pn;i++) p_momenta[4+i*5]=flavs[i].Mass();
  m_as   = m_cpl*2.*M_PI;
  m_aqed = MODEL::s_model->ScalarFunction(std::string("alpha_QED"),rpa.gen.CplScale());

  for (size_t i=0;i<3;i++) p_result[i]=0.;
  p_result[3]=1.;
}

void LH_OLE_Interface::Calc(const Vec4D_Vector& momenta) {
  if (!m_active) return;
  if (m_OLE_id<0) return;
  m_bf  = m_born*m_cpl;

  for (size_t i=0;i<m_pn;i++) {
    p_momenta[0+i*5]=momenta[i][0];
    p_momenta[1+i*5]=momenta[i][1];
    p_momenta[2+i*5]=momenta[i][2];
    p_momenta[3+i*5]=momenta[i][3];
  }

  OLE::EvalSubprocess(m_OLE_id,p_momenta,m_mur2,m_as,m_aqed,p_result);
  // finite
  m_res.Finite()= p_result[2]/p_result[3]*m_bf;
  // 1/epsIR
  m_res.IR()=  p_result[1]/p_result[3]*m_bf;
  // 1/epsIR2
  m_res.IR2()= p_result[0]/p_result[3]*m_bf;
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
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    return new LH_OLE_Interface(pi, fl, true);
  }
  else if (pi.m_fi.m_nloqcdtype&nlo_type::vsub) {
    return new LH_OLE_Interface(pi, fl, false);
  }
  return NULL;
}
