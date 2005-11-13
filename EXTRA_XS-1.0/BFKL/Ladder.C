#include "Ladder.H"

#include "Phase_Space_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Run_Parameter.H"
#include "Random.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

Ladder::Ladder(const size_t nin,const size_t nout,
	       const ATOOLS::Flavour *flavours,
	       const int scalescheme,const int kfactorscheme,
	       BEAM::Beam_Spectra_Handler *const beamhandler,
	       PDF::ISR_Handler *const isrhandler,
	       ATOOLS::Selector_Data *const selectordata):
  XS_Group(nin,nout,flavours,scalescheme,kfactorscheme,
	   beamhandler,isrhandler,selectordata),
  p_sudakov(new BFKL_Sudakov()) 
{ 
  m_name="BFKL_ME";
}

Ladder::~Ladder()
{
  delete p_sudakov;
}

bool Ladder::Initialize()
{
  m_moms.resize(m_nin);
  m_kt2min=0.0;
  for (short unsigned int i(0);i<m_nin;++i) {
    m_flavs.push_back(p_flavours[i]);
    p_pdfs[i]=p_isrhandler->PDF(i);
    m_kt2min=ATOOLS::Max(m_kt2min,p_pdfs[i]->Cut("kp"));// is kt2 in fact :(
    if (p_pdfs[i]==NULL || 
	p_pdfs[i]->Type().find("DUPDF")==std::string::npos)
      THROW(fatal_error,"BFKL ME needs UPDF.");
  }
  p_sudakov->SetKT2Min(m_kt2min);
  delete [] p_momenta;
  delete [] p_flavours;
  p_momenta=&m_moms.front();
  p_flavours=&m_flavs.front();
  p_addflavours = new ATOOLS::Flavour[2];
  p_addflavours[0]=ATOOLS::kf::gluon;
  p_addflavours[1]=ATOOLS::kf::gluon;
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  return true;
}

void Ladder::CreateFSRChannels() 
{
  p_pshandler->FSRIntegrator()->DropAllChannels();
}

void Ladder::CreateISRChannels() 
{
  p_pshandler->ISRIntegrator()->DropAllChannels();
  p_pshandler->KMRKPIntegrator()->DropAllChannels();
  p_pshandler->KMRZIntegrator()->DropAllChannels();
}

double Ladder::Jacobian()
{
  double jacobian(1.0);
  jacobian*=sqrt(m_z1/(1.0-m_z1)*m_kt21/m_q2);
  jacobian*=m_k2[0]*m_k2[3]/(m_a2*m_z2*(1.0-m_z2));
  return jacobian;
}

bool Ladder::GeneratePSJet()
{
  double rn[4], Q(sqrt(m_q2));
  while (true) {
    for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
    double ymin(Min(m_yb,0.5*log(4.0*m_kt2min/m_q2)));
    double ymax(ATOOLS::Max(m_ya,-0.5*log(4.0*m_kt2min/m_q2)));
    // dice y
    m_y1=ymin+rn[0]*(ymax-ymin);
    double kt2max(sqrt(m_kt2min*m_q2)*exp(-dabs(m_y1)));
    kt2max=ATOOLS::Min(m_q2*exp(-2.0*dabs(m_y1)),kt2max);
    if (m_kt2min>=kt2max) continue;
    // dice kt2
    m_kt21=m_kt2min*pow(kt2max/m_kt2min,rn[1]);
    // dice phi
    m_phi1=2.0*M_PI*rn[2];
    double kt1(sqrt(m_kt21));
    double a(kt1/Q*exp(m_y1));
//     PRINT_INFO(kt1<<" "<<m_y1<<" -> "<<a);
    // construct emission
    double pp(a*m_pa.PPlus()), pm(m_kt21/pp);
    m_k1=Vec4D(0.5*(pp+pm),kt1*cos(m_phi1),kt1*sin(m_phi1),0.5*(pp-pm));
    // dice z
    m_z1=rn[3]*(1.0-a);
    m_a1=m_z1/(1.0-m_z1)*a;
    pp=m_a1*m_pa.PPlus();
    pm=m_pa.PMinus()-pm;
    m_moms[0]=m_q=Vec4D(0.5*(pp+pm),-m_k1[1],-m_k1[2],0.5*(pp-pm));
//     PRINT_INFO(m_kt21/m_k1.PPerp2()<<" & "<<m_y1/m_k1.Y());
//     PRINT_INFO(a<<" "<<m_a1<<" "<<m_z1<<" "<<m_k1<<" "<<m_k1.Abs2());
//     PRINT_INFO(a<<" "<<m_a1<<" "<<m_z1<<" "<<m_q<<" "<<m_q.Abs2());
    if (a>1.0 || m_a1>1.0 || m_z1>1.0) {
      msg.Error()<<METHOD<<"(): LCM out of range."<<std::endl;
      return false;
    }
    break;
  }
  if (m_k1.Nan() || m_q.Nan()) return false;
  return true;
}

bool Ladder::ConstructRung()
{
  return true;
}

bool Ladder::GenerateLadder()
{
  while (true) {
    m_moms.resize(m_nin);
    m_flavs.resize(m_nin);
    p_sudakov->Initialize();
    if (!GeneratePSJet()) return false;
    p_sudakov->SetKT2Max(sqr(rpa.gen.Ecms()));// subtract pdf jet kt
    while (p_sudakov->Dice()) {
      ConstructRung();
      //   if (!p_sudakov->CalculateWeight(m_moms)) {
      //     m_moms.pop_back();
      //     continue;
      //   }
      // m_weight*=p_sudakov->Weight();
      // p_sudakov->SetKT2Max(sqr(rpa.gen.Ecms()));// subtract pdf jet kt
    }
    double pp(m_q.PPlus()), pm(m_q.PPerp2()/pp);
    m_k2=Vec4D(0.5*(pp+pm),m_q[1],m_q[2],0.5*(pp-pm));
    m_a2=-m_q.PMinus()/m_pb.PMinus();
    m_kt22=m_k2.PPerp2();
    m_moms[1]=-1.0*m_q;
    m_q=m_q-m_k2;
    if (!m_q.Nan() && m_q[0]<0.0 && m_q[0]>-m_pb[0]) {
//       PRINT_INFO(pp<<" "<<pm<<" "<<m_k2<<" "<<m_k2.Abs2());
//       PRINT_INFO(pp<<" "<<pm<<" "<<m_moms[1]<<" "<<m_moms[1].Abs2());
      break; 
    }
  }
  m_z2=1.0+m_k2.PMinus()/m_q.PMinus();
  m_weight*=Jacobian();
  return true;
}

bool Ladder::SetScales()
{
  double fac(ATOOLS::rpa.gen.FactorizationScaleFactor());
  if (m_moms.size()>2) {
    m_mu21=fac*sqr(m_a1)*m_q2*m_moms[2].PMinus()/m_moms[2].PPlus();
    m_mu22=fac*sqr(m_a2)*m_q2*m_moms.back().PPlus()/m_moms.back().PMinus();
  }
  else {
    m_mu21=fac*sqr(m_a1)*m_q2*m_k2.PMinus()/m_k2.PPlus();
    m_mu22=fac*sqr(m_a2)*m_q2*m_k1.PPlus()/m_k1.PMinus();
  }
  return true;
}

double Ladder::Differential(const ATOOLS::Vec4D *momenta)
{
  m_weight=1.0;
  m_pa=p_beamhandler->GetBeam(0)->InMomentum();
  m_pb=p_beamhandler->GetBeam(1)->InMomentum();
  m_q2=2.0*m_pa*m_pb;
  p_sudakov->SetYMax(m_ya=m_pa.Y());
  p_sudakov->SetYMin(m_yb=m_pb.Y());
  if (!GenerateLadder()) {
    msg.Error()<<METHOD<<"(..): No ME generated. Set weight 0."<<std::endl;
    return 0.0;
  }
  if (!SetScales()) {
    msg.Error()<<METHOD<<"(..): Invalid scales. Set weight 0."<<std::endl;
    return 0.0;
  }
  p_pdfs[0]->Calculate(m_a1,m_z1,m_kt21,m_mu21);
  p_pdfs[1]->Calculate(m_a2,m_z2,m_kt22,m_mu22);
  m_weight*=p_pdfs[0]->GetXPDF(m_flavs[0]);
  m_weight*=p_pdfs[1]->GetXPDF(m_flavs[1]);
//   if (rpa.gen.NumberOfDicedEvents()>1) PRINT_INFO(m_k1<<" "<<m_k2);
  m_nvector=m_moms.size();
  m_nout=m_nvector-m_nin;
  p_momenta=&m_moms.front();
//   PRINT_INFO(p_momenta[0]<<p_momenta[1]);
  p_flavours=&m_flavs.front();
  m_nstrong=m_nin+m_nout;
  p_addmomenta[0]=m_k1;
  p_addmomenta[1]=m_k2;
  return m_weight;
}

double Ladder::Differential2()
{
  return 0.0;
}

bool Ladder::CalculateTotalXSec(const std::string &rpath,
				const bool create)
{
  return true;
}

bool Ladder::SetColours(const ATOOLS::Vec4D *momenta)
{
  return true;
}

bool Ladder::SetColours(const double s,const double t,const double u)
{
  return true;
}

bool Ladder::SelectOne()
{
  p_selected=this;
  return true;
}

void Ladder::DeSelect()
{
  p_selected=NULL;
}

