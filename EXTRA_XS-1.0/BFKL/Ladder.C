#include "Ladder.H"

#include "Phase_Space_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Flow.H"

#define NC 3.0

// #define USING__Collinear_Factorization
#define USING__Gluon_Kernel_Only
#define USING__Two_To_Zero_ME

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
#ifdef USING__Gluon_Kernel_Only
    p_pdfs[i]->SetSplitMode(1);
#endif
  }
  p_pdfs[0]->SetSudakovMode(1);
  p_pdfs[1]->SetSudakovMode(0);
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
#ifdef USING__Collinear_Factorization
  // jacobian for ds' dct -> dkt2 dz 
  return 2.0*sqr(m_kt21)/(m_z1*(1.0-m_z1));
#else
  // jacobian for ds' -> dz 
  return 2.0*sqr(m_kt21)/(m_z1*(1.0-m_z1));
#endif
}

double Ladder::Flux() const
{
  return 2.0*m_q2*m_a1*m_a2/(m_z1*m_z2);
}

bool Ladder::CheckEnergy() const
{
  return -m_q.PMinus()/(1.0+m_q.PPerp2()/m_q.Abs2())<m_pb.PMinus();
}

bool Ladder::GeneratePDFJet()
{
  m_weight=1.0;
  double rn[4], Q(sqrt(m_q2)), kt2max(0.25*m_q2);
  double kt2min(0.25*m_q2*p_pdfs[0]->XMin()*p_pdfs[1]->XMin());
  kt2min=ATOOLS::Max(m_kt2min,kt2min);
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  // dice kt2
  m_kt21=kt2min*pow(kt2max/kt2min,rn[0]);
  m_weight*=log(kt2max/kt2min)*m_kt21/m_q2;
  // dice z
  double zmax(sqrt(0.25-m_kt21/m_q2)), zmin(0.5-zmax);
  zmax+=0.5;
  m_z1=zmin+rn[1]*(zmax-zmin);
  m_weight*=zmax-zmin;
  // dice phi
  m_phi1=2.0*M_PI*rn[2];
  m_weight*=2.0*M_PI;
  // dice y
  double shat(m_kt21/(m_z1*(1.0-m_z1))), Qhat(sqrt(shat));
  double ymin(log(Qhat/Q)), ymax(-ymin);
  if (ymin>=ymax) THROW(fatal_error,"No allowed y range.");
  double y(ymin+rn[3]*(ymax-ymin));
  m_weight*=ymax-ymin;
  // construct emission
  double kt1(sqrt(m_kt21)), E1(Qhat/2.0);
  m_k1=Vec4D(E1,kt1*cos(m_phi1),kt1*sin(m_phi1),sqrt(E1*E1-m_kt21));
  Vec4D cm(cosh(y),0.0,0.0,sinh(y));
  Poincare cms(cm);
  cms.BoostBack(m_k1);
  m_y1=m_k1.Y();
  // construct propagator
  Vec4D pi1(Qhat/2.,0.,0.,Qhat/2.);
  cms.BoostBack(pi1);
  m_moms[0]=m_q=pi1-m_k1;
  m_a1=m_q.PPlus()/m_pa.PPlus();
  m_z1=m_q.PPlus()/pi1.PPlus();
  // phase space weight
  m_weight/=32.0*sqr(M_PI);
  if (m_a1>1.0 || m_z1>1.0 || m_a1/m_z1>1.0) {
    msg_Error()<<METHOD<<"(): LCM out of range."<<std::endl;
    return false;
  }
  if (!CheckEnergy() || m_k1.Nan() || m_q.Nan()) return false;
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
  m_moms.resize(m_nin);
  m_flavs.resize(m_nin);
  p_sudakov->Initialize();
  if (!GeneratePDFJet()) return false;
#ifndef USING__Two_To_Zero_ME
  p_sudakov->SetKT2Max(m_kt21);
  p_sudakov->SetYMax(m_y1);
  while (p_sudakov->Dice()) {
    if (ConstructRung()) {
      if (CheckEnergy() && p_sudakov->CalculateWeight(m_moms,m_q)) {
	m_weight*=p_sudakov->Weight();
      }
      else {
	m_moms.pop_back();
	m_flavs.pop_back();
	m_q=m_oldq;
      }
    }
  }
#endif
  double pp(m_q.PPlus()), pm(m_q.PPerp2()/pp);
  m_k2=Vec4D(0.5*(pp+pm),m_q[1],m_q[2],0.5*(pp-pm));
  m_a2=-m_q.PMinus()/m_pb.PMinus();
  m_kt22=m_k2.PPerp2();
  m_moms[1]=-1.0*m_q;
  m_q=m_q-m_k2;
  m_z2=-m_moms[1].PMinus()/m_q.PMinus();
  if (m_moms.size()>m_nin && m_moms.back().Y()<=m_k2.Y()) return false;
  if (m_q.Nan() || m_q[0]>=0.0 || m_q[0]<=-m_pb[0]) return false; 
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
    msg_Debugging()<<METHOD<<"(..): No ME generated. Set weight 0.\n";
    return 0.0;
  }
  if (!SetScales()) {
    msg_Error()<<METHOD<<"(..): Invalid scales. Set weight 0."<<std::endl;
    return 0.0;
  }
#ifdef USING__Collinear_Factorization
  /*
    double s((m_k1+m_k2).Abs2()), t(m_moms[0].Abs2());
    double u((m_moms[0]+m_k1-m_k2).Abs2()), Q2(2.0*s*t*u/(s*s+u*u+t*t));
    p_pdfs[0]->GetBasicPDF()->Calculate(m_a1/m_z1,0.0,0.0,Q2);
    p_pdfs[1]->GetBasicPDF()->Calculate(m_a2/m_z1,0.0,0.0,Q2);
    m_weight*=p_pdfs[0]->GetBasicPDF()->GetXPDF(m_flavs[0])/m_a1*m_z1;
    m_weight*=p_pdfs[1]->GetBasicPDF()->GetXPDF(m_flavs[1])/m_a2*m_z2;
    double Ms(1.0-t*u/(s*s)), Mt(1.0-s*u/(t*t)), Mu(1.0-s*t/(u*u));
    m_weight*=sqr(4.0*M_PI*(*MODEL::as)(Q2))*9.0/4.0*(Mt+Mu+Ms)/sqr(m_kt21);
  */
  // ll approximation of gg->gg me
  /*
    m_weight*=4.0*sqr(M_PI)*sqr(M_PI)/sqr(m_kt21)*
      sqr((*MODEL::as)(Q2)/(2.0*M_PI*m_kt21)*2.0*3.0/(m_z1*(1.0-m_z1)))*
        4.0*sqr((m_moms[0]+m_moms[1]+m_k1+m_k2).Abs2())/(NC*NC-1.0)*
        sqr(m_z1*(1.0-m_z1));
  */
  // ll approximation of gg->gg me
  p_pdfs[0]->GetBasicPDF()->Calculate(m_a1/m_z1,0.0,0.0,m_kt21);
  p_pdfs[1]->GetBasicPDF()->Calculate(m_a2/m_z1,0.0,0.0,m_kt22);
  m_weight*=p_pdfs[0]->GetBasicPDF()->GetXPDF(m_flavs[0])/m_a1*m_z1;
  m_weight*=p_pdfs[1]->GetBasicPDF()->GetXPDF(m_flavs[1])/m_a2*m_z2;
  m_weight*=4.0*pow(M_PI,4.0)*
    sqr((*MODEL::as)(m_kt21)/(2.0*M_PI*m_kt21)*
	2.0*3.0*sqr(1.0-m_z1*(1.0-m_z1))/(m_z1*(1.0-m_z1)));
#else
  p_pdfs[0]->Calculate(m_a1,m_z1,m_kt21,m_mu21);
  p_pdfs[1]->Calculate(m_a2,m_z2,m_kt22,m_mu22);
  m_weight*=p_pdfs[0]->GetXPDF(m_flavs[0])/m_a1*m_z1;
  m_weight*=p_pdfs[1]->GetXPDF(m_flavs[1])/m_a2*m_z2;
  // ll approximation of gg->gg me
  m_weight*=4.0*pow(M_PI,4.0);
#endif  
  m_weight*=Jacobian();
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

