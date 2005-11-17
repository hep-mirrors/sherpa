#include "Ladder.H"

#include "Phase_Space_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Flow.H"

#define NC 3.0

using namespace EXTRAXS;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Ladder::Ladder(const size_t nin,const size_t nout,
	       const ATOOLS::Flavour *flavours,
	       const int scalescheme,const int kfactorscheme,
	       BEAM::Beam_Spectra_Handler *const beamhandler,
	       PDF::ISR_Handler *const isrhandler,
	       ATOOLS::Selector_Data *const selectordata):
  XS_Group(nin,nout,flavours,scalescheme,kfactorscheme,
	   beamhandler,isrhandler,selectordata),
  p_sudakov(new BFKL_Sudakov()),
  m_ncols(0)
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
    p_pdfs[i]=dynamic_cast<Doubly_Unintegrated_PDF*>(p_isrhandler->PDF(i));
    if (p_pdfs[i]==NULL || 
	p_pdfs[i]->Type().find("DUPDF")==std::string::npos)
      THROW(fatal_error,"BFKL ME needs UPDF.");
    m_kt2min=ATOOLS::Max(m_kt2min,p_pdfs[i]->Cut("kp"));// is kt2 in fact :(
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

double Ladder::Jacobian() const
{
  double jacobian(m_a1);
  jacobian*=(m_a2*m_z2*(1.0-m_z2))/(m_k2[0]*m_k2[3])/(2.0*M_PI);
  return jacobian;
}

double Ladder::Flux() const
{
  double flux(2.0*m_q2);
  double b1(m_kt21/(m_a1*m_q2)), b2(m_kt22/(m_a2*m_q2));
  flux*=m_a1*m_a2-b1*b2;
  return flux;
}

bool Ladder::GeneratePDFJet()
{
  double rn[4], Q(sqrt(m_q2)), mu(sqrt(m_kt2min));
  while (true) {
    for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
    double ymin(m_yb), ymax(m_ya);
    if (ymin>=ymax) THROW(fatal_error,"No allowed y range.");
    // dice y
    m_y1=ymin+rn[0]*(ymax-ymin);
    m_weight=ymax-ymin;
    double kt2max(ATOOLS::Min(m_kt2min*exp(ymax-ymin),m_q2/4.0));
    if (m_kt2min>=kt2max) continue;
    // dice kt2
    m_kt21=m_kt2min*pow(kt2max/m_kt2min,rn[1]);
    m_weight*=log(kt2max/m_kt2min)*m_kt21;
    // dice phi
    m_phi1=2.0*M_PI*rn[2];
    // no weight, averaging
    double kt1(sqrt(m_kt21)), a(kt1/Q*exp(m_y1));
    if (a>1.0) continue;
    // construct emission
    double pp(a*m_pa.PPlus()), pm(m_kt21/pp+m_pa.PMinus());
    m_k1=Vec4D(0.5*(pp+pm),kt1*cos(m_phi1),kt1*sin(m_phi1),0.5*(pp-pm));
    double zmin(mu/(mu+kt1)), zmax(Min(kt1/(mu+kt1),1.0-a));
    // dice z
    m_z1=zmin+rn[3]*(zmax-zmin);
    m_weight*=zmax-zmin;
    m_a1=m_z1/(1.0-m_z1)*a;
    pp=m_a1*m_pa.PPlus();
    pm=m_pa.PMinus()-pm;
    m_moms[0]=m_q=Vec4D(0.5*(pp+pm),-m_k1[1],-m_k1[2],0.5*(pp-pm));
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
  double y(p_sudakov->GetY());
  double kt(sqrt(p_sudakov->GetKT2())), phi(p_sudakov->GetPhi());
  Vec4D k(kt*cosh(y),kt*cos(phi),kt*sin(phi),kt*sinh(y));
  m_oldq=m_q;
  m_q=m_q-k;
  m_moms.push_back(k);
  m_flavs.push_back(p_sudakov->Selected()->GetB());
  return true;
}

bool Ladder::GenerateLadder()
{
  while (true) {
    m_moms.resize(m_nin);
    m_flavs.resize(m_nin);
    p_sudakov->Initialize();
    if (!GeneratePDFJet()) return false;
    p_sudakov->SetKT2Max(m_q2/4.0);
    p_sudakov->SetYMax(m_y1);
    while (false) {
//     while (p_sudakov->Dice()) {
      if (ConstructRung()) {
	if (p_sudakov->CalculateWeight(m_moms,m_q)) {
	  m_weight*=p_sudakov->Weight();
	}
	else {
	  m_moms.pop_back();
	  m_flavs.pop_back();
	  m_q=m_oldq;
	}
      }
    }
    double pp(m_q.PPlus()+m_pb.PPlus()), pm(m_q.PPerp2()/pp);
    m_k2=Vec4D(0.5*(pp+pm),m_q[1],m_q[2],0.5*(pp-pm));
    m_a2=-m_q.PMinus()/m_pb.PMinus();
    m_kt22=m_k2.PPerp2();
    m_moms[1]=-1.0*m_q;
    m_q=m_q-m_k2;
    m_z2=-m_moms[1].PMinus()/m_q.PMinus();
    double mu(sqrt(m_kt2min)), kt2(sqrt(m_kt22));
    if (m_z2<mu/(mu+kt2) || m_z2>kt2/(mu+kt2)) continue;
    if (m_moms.size()>m_nin) {
      if (m_moms.back().Y()<=m_k2.Y()) continue;
    }
    else {
      if (m_y1<=m_k2.Y()) continue;
    }
    if (!m_q.Nan() && m_q[0]<0.0 && m_q[0]>-m_pb[0]) break; 
  }
  m_weight*=Jacobian();
  return true;
}

bool Ladder::SetScales()
{
  double fac(ATOOLS::rpa.gen.FactorizationScaleFactor());
  if (m_moms.size()>m_nin) {
    m_mu21=fac*m_moms[2].PPerp2()*
      sqr(m_moms.front().PPlus()/m_moms[2].PPlus());
    m_mu22=fac*m_moms.back().PPerp2()*
      sqr(m_moms[1].PMinus()/m_moms.back().PMinus());
  }
  else {
    m_mu21=fac*m_kt22/sqr(1.0-m_z2);
    m_mu22=fac*m_kt21/sqr(1.0-m_z1);
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
  p_pdfs[1]->SetSudakovMode(0);
  p_pdfs[0]->Calculate(m_a1,m_z1,m_kt21,m_mu21);
  p_pdfs[1]->Calculate(m_a2,m_z2,m_kt22,m_mu22);
  m_weight*=p_pdfs[0]->GetXPDF(m_flavs[0])/m_a1;
  m_weight*=p_pdfs[1]->GetXPDF(m_flavs[1])/m_a2;
  m_weight*=4.0*sqr((m_moms[0]+m_moms[1]+m_k1+m_k2).Abs2())/(NC*NC-1.0);
  m_weight/=Flux();
  if (m_ncols<m_moms.size()) {
    for (size_t i(0);i<m_ncols;++i) delete [] p_colours[i]; 
    if (m_ncols>0) delete [] p_colours;
    m_ncols=m_moms.size();
    p_colours = new int*[m_ncols];
    for (size_t i(0);i<m_ncols;++i) p_colours[i] = new int[2];
  }
  p_colours[0][0]=p_colours[1][1]=Flow::Counter();
  for (size_t i(1);i<m_moms.size();++i)
    p_colours[i][0]=p_colours[i-1][1]=Flow::Counter();
  p_colours[m_moms.size()-1][1]=p_colours[0][0];
  std::swap<int>(p_colours[0][0],p_colours[0][1]);
  std::swap<int>(p_colours[1][0],p_colours[1][1]);
  m_nvector=m_moms.size();
  m_maxjetnumber=m_nout=m_nvector-m_nin;
  p_momenta=&m_moms.front();
  p_flavours=&m_flavs.front();
  m_nstrong=m_nin+m_nout;
  p_addmomenta[0]=m_k1;
  p_addmomenta[1]=m_k2;
  return dabs(m_weight);
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

