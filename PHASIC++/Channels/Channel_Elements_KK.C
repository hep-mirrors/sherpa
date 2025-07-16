#include "PHASIC++/Channels/Channel_Elements_KK.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_MPI.H"
#include "MODEL/Main/Model_Base.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Channel_Elements_KK PHASIC::CEKK;

Channel_Elements_KK::Channel_Elements_KK()
{
  m_ed    = 0;
  m_rs    = 0;
}

Channel_Elements_KK::~Channel_Elements_KK()
{
}

void Channel_Elements_KK::Init(int nin,int nout,ATOOLS::Flavour* fl)
{
  if (m_rs>0) return;
  m_nin = nin;
  m_nout = nout;
  m_kkp=-1;m_mpss=1.;
  int mode = MODEL::s_model->ScalarNumber(std::string("KK_mode"));
  for (int i=m_nin;i<m_nin+m_nout;i++) {
    if(fl[i].IsKK() && (mode==1 || mode==2 || mode==5)){
      if(ATOOLS::IsZero(ATOOLS::sqr(fl[i].Mass()))){
	msg_Error()<<"Error in Channel_Elements_KK: "<<endl
		   <<"   Please initialize with nonzero particle mass ("<<fl[i]<<") !"<<std::endl;
	Abort();
      }
      m_kkp=i;
      
      m_ed  = MODEL::s_model->ScalarNumber(std::string("ED"));
      m_r2  = sqr(MODEL::s_model->ScalarConstant(std::string("Radius")));
      m_gn  = MODEL::s_model->ScalarConstant(std::string("G_Newton"));

      //Calculation of Gamma(ed/2)
      if(m_ed%2==0) m_gam=1.;
      else m_gam=sqrt(M_PI);
      for(int k=2-m_ed%2;k<m_ed;k+=2)m_gam*=0.5*k;

      m_prevET = 0.;
      break;
    }
  }
  m_rs = 1;
}

void Channel_Elements_KK::SetKKmass(double *ms, double ET, Cut_Data* cuts,double ran)
{
  if (!IsEqual(ET,m_prevET) && m_kkp>-1) {
    double mm = m_prevET = ET;
    for(int j=m_nin;j<m_nin+m_nout;j++)
      if(j!=m_kkp) mm -= Max(sqrt(ms[j]),cuts->etmin[j]);
    if (m_nout==2) mm = Min(mm,sqrt(sqr(ET)-2*cuts->etmin[5-m_kkp]*ET));
    m_maxm2 = sqr(mm);
    m_maxn  = sqrt(m_maxm2*m_r2);
    m_mpss  = 2.*pow(sqrt(M_PI)*m_maxn,double(m_ed))/m_gam;
  }
  m_sran = ran;
  double ms2=sqr(ran)*m_maxm2;
  m_weight=pow(ran,double(m_ed-1));

  ms[m_kkp]=ms2;
}

double Channel_Elements_KK::GetWeightKK(double& ran){  
  ran = m_sran;
  return m_mpss*m_weight;
}
