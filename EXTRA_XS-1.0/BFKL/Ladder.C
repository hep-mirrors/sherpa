#include "Ladder.H"

#include "Phase_Space_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Flow.H"

// #define USING__Two_To_Zero_ME
#define USING__Gluon_Kernel_Only

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
    m_kt2min=ATOOLS::Max(m_kt2min,p_pdfs[i]->Cut("kp"));
    p_pdfs[i]->SetSudakovMode(0);
#ifdef USING__Gluon_Kernel_Only
    p_pdfs[i]->SetSplitMode(1);
#endif
  }
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

void Ladder::AddEvent(const double xs,const double validxs,const int ncounts)
{ 
  Integrable_Base::AddEvent(xs,validxs,ncounts);
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
  Vec4D P(m_moms[0]+m_k1+m_moms[1]+m_k2);
  return 4.0*(m_k2.PPlus()*P.PMinus()+m_k2.PMinus()*P.PPlus())/P.Abs2();
}

double Ladder::Flux() const
{
  return 2.0*m_q2*m_a1*m_a2/(m_z1*m_z2);
}

bool Ladder::GeneratePDFJet()
{
  double rn[4];
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  // dice y
  m_y=m_y1=m_yb+rn[0]*(m_ya-m_yb);
  m_weight*=m_ya-m_yb;
  // dice kt2
  m_kt21=m_kt2min*pow(0.25*m_q2/m_kt2min,rn[1]);
  m_weight*=log(0.25*m_q2/m_kt2min)*m_kt21/m_q2;
  // dice phi
  m_phi1=2.0*M_PI*rn[2];
  m_weight*=2.0*M_PI;
  // construct jet
  double kt1(sqrt(m_kt21));
  m_k1=Vec4D(kt1*cosh(m_y1),kt1*cos(m_phi1),
	     kt1*sin(m_phi1),kt1*sinh(m_y1));
  m_q=m_k1;
  if (m_k1[0]*m_k1[0]>=0.25*m_q2) return false;
  msg_Debugging()<<"m_k1 = "<<m_k1<<"\n";
  // phase space weight
  m_weight/=32.0*sqr(M_PI);
  // dice y
  m_yn=m_yb+rn[3]*(m_y1-m_yb);
  m_weight*=m_y1-m_yb;
  if (!ConstructIncoming()) {
    msg.Error()<<METHOD<<"(): No kinematics found."<<std::endl;
    return false;
  }
  return true;
}

bool Ladder::TestEmission()
{
  if (m_p1[0]<0.0 || m_p2[0]<0.0 ||
      m_p1[0]>m_pa[0] || m_p2[0]>m_pb[0]) return false;
  return true;
}

bool Ladder::ConstructIncoming()
{
  msg_Debugging()<<"m_q     = "<<m_q<<"\n";
  msg_Debugging()<<"m_oldq  = "<<m_oldq<<"\n";
  double knt(sqrt(m_q.PPerp2())), yn(m_yn);
  m_kn=Vec4D(knt*cosh(yn),-m_q[1],-m_q[2],knt*sinh(yn));
  msg_Debugging()<<"k2 = "<<m_kn<<"\n";
  Vec4D cms(m_k1+m_kn);
  for (size_t i(m_nin);i<m_moms.size();++i)
    cms+=m_moms[i];
  double Q2(0.5*sqrt(cms.Abs2())), ey(exp(cms.Y()));
  msg_Debugging()<<"Q = "<<2.0*Q2<<" <- "<<cms.Abs2()<<", y = "<<cms.Y()<<"\n";
  m_p1=Vec4D(Q2*ey,0.0,0.0,Q2*ey);
  m_p2=Vec4D(Q2/ey,0.0,0.0,-Q2/ey);
  if (m_p1.Nan() || m_p2.Nan() || m_kn.Nan()) return false;
  return true;
}

bool Ladder::ConstructRung()
{
  double y(p_sudakov->GetY());
  double kt(sqrt(p_sudakov->GetKT2())), phi(p_sudakov->GetPhi());
  Vec4D k(kt*cosh(y),kt*cos(phi),kt*sin(phi),kt*sinh(y));
  m_moms.back()=k;
  m_q+=k;
  if (k.Nan()) return false;
  return ConstructIncoming();
}

bool Ladder::GenerateLadder()
{
  m_q=Vec4D();
  m_weight=1.0;
  m_moms.resize(m_nin);
  m_flavs.resize(m_nin);
  if (!GeneratePDFJet()) return false;
  m_k2=m_kn;
  m_moms[0]=m_p1-m_k1;
  m_moms[1]=m_p2-m_k2;
#ifndef USING__Two_To_Zero_ME
  p_sudakov->SetYMax(m_y1);
  p_sudakov->SetYMin(m_yn);
  p_sudakov->SetKT2Min(m_kt2min);
  p_sudakov->SetKT2Max(0.25*m_q2);
  msg_Debugging()<<"set ymin = "<<m_yn<<", ymax = "<<m_y1<<"\n";
  p_sudakov->Initialize();
  int cnt(0);
  while (p_sudakov->Dice()) {
    msg_Debugging()<<"test emission at y = "<<p_sudakov->GetY()
		   <<", qt = "<<sqrt(p_sudakov->GetKT2())<<"\n";
    m_y=p_sudakov->GetY();
    m_oldq=m_q;
    m_moms.push_back(Vec4D());
    m_flavs.push_back(p_sudakov->Selected()->GetB());
    if (ConstructRung() && TestEmission() && 
	p_sudakov->Approve(m_moms,m_oldq,m_q)) {
      msg_Debugging()<<"accept emission at y = "<<p_sudakov->GetY()
		     <<", qt = "<<sqrt(p_sudakov->GetKT2())<<"\n";
      m_k2=m_kn;
      m_moms[0]=m_p1-m_k1;
      m_moms[1]=m_p2-m_k2;
      ++cnt;
    }
    else {
      m_moms.pop_back();
      m_flavs.pop_back();
      m_q=m_oldq;
    }
  }
#endif
  msg_Debugging()<<"k1 = "<<m_k1<<"\n";
  msg_Debugging()<<"k2 = "<<m_k2<<"\n";
  msg_Debugging()<<"p1 = "<<m_moms[0]<<"\n";
  msg_Debugging()<<"p2 = "<<m_moms[1]<<"\n";
  return true;
}

bool Ladder::SetScales()
{
  m_a1=m_moms[0].PPlus()/m_pa.PPlus();
  m_a2=m_moms[1].PMinus()/m_pb.PMinus();
  m_z1=m_k1.PPlus()/m_moms[0].PPlus();
  m_z1=1.0/(1.0+m_z1);
  m_z2=m_k2.PMinus()/m_moms[1].PMinus();
  m_z2=1.0/(1.0+m_z2);
  m_kt22=m_k2.PPerp2();
  msg_Debugging()<<"a1 = "<<m_a1<<"\n";
  msg_Debugging()<<"a2 = "<<m_a2<<"\n";
  msg_Debugging()<<"z1 = "<<m_z1<<"\n";
  msg_Debugging()<<"z2 = "<<m_z2<<"\n";
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
  msg_Debugging()<<"====================\n";
  // set beam parameters
  m_pa=p_beamhandler->GetBeam(0)->InMomentum();
  m_pb=p_beamhandler->GetBeam(1)->InMomentum();
  m_q2=2.0*m_pa*m_pb;
  m_ya=m_pa.Y();
  m_yb=m_pb.Y();
  msg_Debugging()<<"ya/b set to "<<m_ya<<" "<<m_yb<<"\n";
  // construct me
  if (!GenerateLadder()) {
    msg_Debugging()<<METHOD<<"(..): No ME generated. Set weight 0.\n";
    return 0.0;
  }
  // set scales
  if (!SetScales()) {
    msg.Error()<<METHOD<<"(..): Invalid scales. Set weight 0."<<std::endl;
    return 0.0;
  }
  // calculate pdfs
  p_pdfs[0]->Calculate(m_a1,m_z1,m_kt21,m_mu21);
  p_pdfs[1]->Calculate(m_a2,m_z2,m_kt22,m_mu22);
  m_weight*=p_pdfs[0]->GetXPDF(m_flavs[0])/m_a1*m_z1*m_kt21;
  m_weight*=p_pdfs[1]->GetXPDF(m_flavs[1])/m_a2*m_z2*m_kt22;
  // add ll approximation of gg->gg me
  m_weight*=4.0*pow(M_PI,4.0);
  // add jacobian and flux
  m_weight*=Jacobian();
  m_weight/=Flux();
  // set colours
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
  // set momenta
  m_nvector=m_moms.size();
  m_maxjetnumber=m_nout=m_nvector-m_nin;
  p_momenta=&m_moms.front();
  p_flavours=&m_flavs.front();
  m_nstrong=m_nin+m_nout;
  p_addmomenta[0]=m_k1;
  p_addmomenta[1]=m_k2;
  Vec4D sum;
  for (size_t i(0);i<m_nvector;++i) {
    if (i>1) sum+=p_momenta[i];
    msg_Debugging()<<i<<" -> "<<p_momenta[i]<<" "<<p_momenta[i].Y()<<"\n";
  }
  msg_Debugging()<<" k1 -> "<<p_addmomenta[0]<<" "<<p_addmomenta[0].Y()<<"\n";
  msg_Debugging()<<" k2 -> "<<p_addmomenta[1]<<" "<<p_addmomenta[1].Y()<<"\n";
  msg_Debugging()<<" sum "<<sum<<"\n";
  msg_Debugging()<<" sum "<<m_moms[0]+m_moms[1]<<"\n";
  msg_Debugging()<<"weight = "<<m_weight<<"\n";
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

