#include "Ladder.H"

#include "Phase_Space_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Flow.H"
#include "Data_Reader.H"
#include "MyStrStream.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Ladder::Ladder(const size_t nin,const size_t nout,
	       const ATOOLS::Flavour *flavours,
	       const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
	       BEAM::Beam_Spectra_Handler *const beamhandler,
	       PDF::ISR_Handler *const isrhandler,
	       ATOOLS::Selector_Data *const selectordata):
  XS_Group(nin,nout,flavours,scalescheme,kfactorscheme,
	   beamhandler,isrhandler,selectordata),
  p_sudakov(new BFKL_Sudakov()),
  m_ncols(0), m_sudmode(1), m_nfixed(std::numeric_limits<size_t>::max()),
  m_multimode(1), m_splitmode(1)
{ 
  m_name="BFKL_ME";
  Data_Reader read(" ",";","!","=");
  int ktscheme(1);
  if (!read.ReadFromFile(ktscheme,"BFKL_KT_SCHEME")) ktscheme=1;
  else msg_Info()<<METHOD<<"(): Set k_T-scheme "<<ktscheme<<".\n";
  p_sudakov->SetKTScheme(ktscheme);
  if (!read.ReadFromFile(m_nfixed,"BFKL_FIXED_MULTI"))
    m_nfixed=std::numeric_limits<size_t>::max();
  else msg_Info()<<METHOD<<"(): Fixed multiplicity "<<m_nfixed<<".\n";
  if (!read.ReadFromFile(m_multimode,"BFKL_MULTI_MODE")) m_multimode=1;
  else msg_Info()<<METHOD<<"(): Multiplicity selection mode "
		 <<m_multimode<<".\n";
  if (!read.ReadFromFile(m_sudmode,"BFKL_SUDAKOV_MODE")) m_sudmode=1;
  else msg_Info()<<METHOD<<"(): Set Sudakov mode "<<m_sudmode<<".\n";
  if (!read.ReadFromFile(m_splitmode,"BFKL_SPLIT_MODE")) m_splitmode=0;
  else msg_Info()<<METHOD<<"(): Set splitting mode "<<m_splitmode<<".\n";
  p_sudakov->SetSplitMode(m_splitmode);
  if (m_multimode==0) m_splitmode=0;
}

Ladder::~Ladder()
{
  delete p_sudakov;
  p_momenta=NULL;
  p_flavours=NULL;
}

bool Ladder::Initialize()
{
  m_moms.resize(m_nin);
  m_kt2min=0.0;
  for (short unsigned int i(0);i<m_nin;++i) {
    m_flavs.push_back(kf::gluon);
    p_pdfs[i]=dynamic_cast<Doubly_Unintegrated_PDF*>(p_isrhandler->PDF(i));
    if (p_pdfs[i]==NULL || 
	p_pdfs[i]->Type().find("DUPDF")==std::string::npos)
      THROW(fatal_error,"BFKL ME needs UPDF.");
    m_kt2min=ATOOLS::Max(m_kt2min,p_pdfs[i]->Cut("kp"));
    p_pdfs[i]->SetSudakovMode(0);
    if (m_splitmode<1) p_pdfs[i]->SetSplitMode(0);
  }
  delete [] p_momenta;
  delete [] p_flavours;
  p_momenta=&m_moms.front();
  p_flavours=&m_flavs.front();
  p_addflavours = new ATOOLS::Flavour[4];
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  return p_sudakov->Initialize();
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

double Ladder::Flux() const
{
  return 2.0*m_q2*m_a1*m_an/(m_z1*m_zn);
}

bool Ladder::GeneratePDFJet()
{
  double rn[4];
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  // dice y
  m_y=m_y1=m_yb+rn[0]*(m_ya-m_yb);
  m_weight*=m_ya-m_yb;
  // dice kt2
  m_kt2=m_kt21=m_kt2min*pow(0.25*m_q2/m_kt2min,rn[1]);
  m_weight*=log(0.25*m_q2/m_kt2min)*m_kt21;
  // dice phi
  double phi1(2.0*M_PI*rn[2]);
  m_weight*=2.0*M_PI;
  // construct jet
  double kt1(sqrt(m_kt21));
  m_k1=Vec4D(kt1*cosh(m_y1),kt1*cos(phi1),
	     kt1*sin(phi1),kt1*sinh(m_y1));
  m_q=m_k1;
  if (m_k1[0]*m_k1[0]>=0.25*m_q2) return false;
  // dice y
  m_yn=m_yb+rn[3]*(m_ya-m_yb);
  m_weight*=m_ya-m_yb;
  if (m_yn<m_y1) return false;
  // dice first flavour
  m_props.back()=kf::gluon;
  if (m_splitmode>0) {
    // select t-channel quark w/ probability T_R/C_A
    if (ran.Get()<0.5/3.0)
      m_props.back()=(kf::code)Min(p_sudakov->Nf(),(size_t)
				   (p_sudakov->Nf()*ran.Get()+1));
  }
  m_flavs[0]=m_props.back();
  // construct incoming
  if (!ConstructIncoming() || 
      !TestEmission()) return false;
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
  double knt(sqrt(m_q.PPerp2())), yn(m_yn);
  m_kn=Vec4D(knt*cosh(yn),-m_q[1],-m_q[2],knt*sinh(yn));
  Vec4D cms(m_k1+m_kn);
  for (size_t i(m_nin);i<m_moms.size();++i)
    cms+=m_moms[i];
  double Q2(0.5*sqrt(cms.Abs2())), ey(exp(cms.Y()));
  m_p1=Vec4D(Q2*ey,0.0,0.0,Q2*ey);
  m_p2=Vec4D(Q2/ey,0.0,0.0,-Q2/ey);
  if (m_p1.Nan() || m_p2.Nan() || m_kn.Nan() ||
      !TestEmission()) return false;
  m_moms[0]=m_p1-m_k1;
  m_moms[1]=m_p2-m_kn;
  m_flavs[1]=m_props.back().Bar();
  return true;
}

bool Ladder::ConstructRung()
{
  double kt(sqrt(m_kt2));
  Vec4D k(kt*cosh(m_y),kt*cos(m_phi),kt*sin(m_phi),kt*sinh(m_y));
  m_q+=m_moms.back()=k;
  if (k.Nan()) return false;
  return ConstructIncoming();
}

bool Ladder::DiceOneEmission()
{
  /* store old values */
  double lasty(m_y);
  m_oldq=m_q;
  /* dice y, ps weight */
  m_y=ran.Get()*(lasty-m_yn)+m_yn;
  m_weight*=dabs(lasty-m_yn);
  /* dice kt2, ps weight */
  double kt2max((m_pa+m_pb-m_q).Abs2());
  m_kt2=m_kt2min*pow(kt2max/m_kt2min,ran.Get());
  m_weight*=log(kt2max/m_kt2min);
  /* dice phi, ps weight is 1 */
  m_phi=ran.Get()*2.0*M_PI;
  /* add new rung */
  m_moms.push_back(Vec4D());
  m_flavs.push_back(kf::gluon);
  m_props.push_back(kf::gluon);
  if (!ConstructRung() ||
      m_q.PPerp2()<m_kt2min) return false;
  /* alpha_s & splitting weight */
  m_weight*=(*MODEL::as)(m_kt2)*3.0/M_PI;
  /* sudakov weight */
  if (m_sudmode>0)
    m_weight*=exp(-(*MODEL::as)(m_oldq.PPerp2())*3.0/M_PI
		  *log(m_oldq.PPerp2()/m_kt2min)*dabs(lasty-m_y));
  m_oldq=m_q;
  return true;
}

bool Ladder::GenerateLadder()
{
  m_q=Vec4D();
  m_weight=1.0;
  m_moms.resize(m_nin);
  m_flavs.resize(m_nin);
  m_props.resize(1);
  if (!GeneratePDFJet()) return false;
  if (m_multimode>0) {
    p_sudakov->SetYA(m_y1);
    p_sudakov->SetYB(m_yn);
    p_sudakov->SetKT2Min(m_kt2min);
    p_sudakov->SetKT2Max((m_pa+m_pb-m_q).Abs2());
    p_sudakov->SetInMomentum(m_q);
    p_sudakov->SetInFlavour(m_props.back());
    p_sudakov->Init();
    msg_Debugging()<<"init sud at y_a = "<<m_y1<<", y_b = "<<m_yn<<"\n";
    size_t cnt(0);
    Vec4D k1(m_k1);
    while (p_sudakov->Dice()) {
      msg_Debugging()<<"test emission at y = "<<p_sudakov->GetY()
		     <<", qt = "<<sqrt(p_sudakov->GetKT2())<<"\n";
      m_oldq=m_q;
      m_y=p_sudakov->GetY();
      m_kt2=p_sudakov->GetKT2();
      m_phi=p_sudakov->GetPhi();
      m_moms.push_back(Vec4D());
      m_flavs.push_back(p_sudakov->Selected()->GetC());
      m_props.push_back(p_sudakov->Selected()->GetB());
      if (!ConstructRung() ||
	  !p_sudakov->Approve(k1,m_oldq,m_moms.back(),m_q)) {
	return false;
      }
      else {
	msg_Debugging()<<"accept emission at y = "<<p_sudakov->GetY()
		       <<", qt = "<<sqrt(p_sudakov->GetKT2())<<"\n";
	m_weight*=p_sudakov->GetWeight();
	k1=m_moms.back();
	p_sudakov->SetKT2Max((m_pa+m_pb-m_q).Abs2());
	p_sudakov->SetInMomentum(m_q);
	p_sudakov->SetInFlavour(m_props.back());
	++cnt;
      }
    }
    if (m_nfixed<std::numeric_limits<size_t>::max() &&
	cnt!=m_nfixed) return false;
    ConstructIncoming();
    m_moms[0]=m_p1-m_k1;
    m_moms[1]=m_p2-m_kn;
  }
  else {
    for (size_t i(0);i<m_nfixed;++i) 
      if (!DiceOneEmission()) return false;
    /* last rung sudakov weight */
    if (m_sudmode>0)
      m_weight*=exp(-(*MODEL::as)(m_q.PPerp2())*3.0/M_PI
		    *log(m_q.PPerp2()/m_kt2min)*dabs(m_y-m_yn));
    m_moms[0]=m_p1-m_k1;
    m_moms[1]=m_p2-m_kn;
  }
  return true;
}

bool Ladder::SetScales()
{
  m_a1=m_moms[0].PPlus()/m_pa.PPlus();
  m_an=m_moms[1].PMinus()/m_pb.PMinus();
  m_z1=m_k1.PPlus()/m_moms[0].PPlus();
  m_z1=1.0/(1.0+m_z1);
  m_zn=m_kn.PMinus()/m_moms[1].PMinus();
  m_zn=1.0/(1.0+m_zn);
  m_kt21=m_k1.PPerp2();
  m_kt2n=m_kn.PPerp2();
  m_mu21=m_kt21/sqr(1.0-m_z1);
  m_mu2n=m_kt2n/sqr(1.0-m_zn);
  return true;
}

double Ladder::Differential(const ATOOLS::Vec4D *momenta)
{
  msg_Debugging()<<"====================\n";
  // set beam parameters
  m_pa=p_beamhandler->GetBeam(0)->InMomentum();
  m_pb=p_beamhandler->GetBeam(1)->InMomentum();
  m_q2=(m_pa+m_pb).Abs2();
  m_ya=m_pa.Y();
  m_yb=m_pb.Y();
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
  double fac(ATOOLS::rpa.gen.FactorizationScaleFactor());
  p_pdfs[0]->Calculate(m_a1,m_z1,m_kt21,fac*m_mu21);
  p_pdfs[1]->Calculate(m_an,m_zn,m_kt2n,fac*m_mu2n);
  m_weight*=p_pdfs[0]->GetXPDF(m_flavs[0])/m_a1*m_z1*m_kt21;
  m_weight*=p_pdfs[1]->GetXPDF(m_flavs[1])/m_an*m_zn*m_kt2n;
  // dice pdf jet flavours
  if (!p_pdfs[0]->SelectJetFlavour
      (p_addflavours[2],m_fl1,ran.Get())) return false;
  if (!p_pdfs[1]->SelectJetFlavour
      (p_addflavours[3],m_fln,ran.Get())) return false;
  // add ll approximation of gg->gg me and symmetry factor
  m_weight*=sqr(M_PI)/(2.0*m_q2);
  // add flux factor
  m_weight/=Flux();
  // set colours
  if (m_ncols<m_moms.size()) {
    for (size_t i(0);i<m_ncols;++i) delete [] p_colours[i]; 
    if (m_ncols>0) delete [] p_colours;
    m_ncols=m_moms.size();
    p_colours = new int*[m_ncols];
    for (size_t i(0);i<m_ncols;++i) p_colours[i] = new int[2];
  }
  SetColours();
  msg_Debugging()<<"c_0 = ("<<p_colours[0][0]
		 <<","<<p_colours[0][1]<<")\n";
  for (size_t i(2);i<m_ncols;++i) {
    msg_Debugging()<<"c_"<<i<<" = ("<<p_colours[i][0]
		   <<","<<p_colours[i][1]<<")\n";
  }
  msg_Debugging()<<"c_1 = ("<<p_colours[1][0]
		 <<","<<p_colours[1][1]<<")\n";
  // set momenta
  m_nvector=m_moms.size();
  m_maxjetnumber=m_nout=m_nvector-m_nin;
  p_momenta=&m_moms.front();
  p_flavours=&m_flavs.front();
  m_nstrong=m_nin+m_nout;
  p_addmomenta[0]=m_k1;
  p_addmomenta[1]=m_kn;
  p_addflavours[0]=m_fl1;
  p_addflavours[1]=m_fln;
  std::string name("BFKL__2_"+ToString(m_nout)+"__"+
		   ToString(m_flavs[0])+"_"+ToString(m_flavs[1]));
  if (m_nout>0) name+="_";
  for (size_t i(2);i<2+m_nout;++i) name+="_"+ToString(m_flavs[i]);
  msg_Debugging()<<METHOD<<"(): Generated '"<<name<<"'\n";
  m_name=name;
  return dabs(m_weight);
}

void Ladder::SetColour(const size_t &i,const bool fw,int *const lc)
{
  Flavour q(m_props[i-2]), f(m_flavs[i]);
  msg_Debugging()<<"step "<<i<<" "<<q<<" -> "<<f<<" "<<m_props[i-1]
		 <<": ("<<lc[0]<<","<<lc[1]<<") -> ";
  if (q.IsGluon()) {
    if (f.IsGluon()) {
      if (ran.Get()>0.5) {
	p_colours[i][0]=lc[0];
	lc[0]=p_colours[i][1]=Flow::Counter();
      }
      else {
	p_colours[i][1]=lc[1];
	lc[1]=p_colours[i][0]=Flow::Counter();
      }
    }
    else if (!f.IsAnti()) {
      p_colours[i][0]=lc[0];
      lc[0]=p_colours[i][1]=0;
    }
    else {
      p_colours[i][1]=lc[1];
      lc[1]=p_colours[i][0]=0;
    }  
  }
  else if (!q.IsAnti()) {
    if (f.IsGluon()) {
      p_colours[i][0]=lc[0];
      lc[0]=p_colours[i][1]=Flow::Counter();
    }
    else {
      lc[1]=p_colours[i][0]=Flow::Counter();
      p_colours[i][1]=0;
    }
  }
  else {
    if (f.IsGluon()) {
      p_colours[i][1]=lc[1];
      lc[1]=p_colours[i][0]=Flow::Counter();
    }
    else {
      lc[0]=p_colours[i][1]=Flow::Counter();
      p_colours[i][0]=0;
    }  
  }
  msg_Debugging()<<"("<<p_colours[i][0]<<","<<p_colours[i][1]
		 <<") ("<<lc[0]<<","<<lc[1]<<")\n";
}

void Ladder::SetColours()
{
  int lc[2]={0,0};
  if (m_flavs[0].IsGluon()||!m_flavs[0].IsAnti())
    lc[0]=Flow::Counter();
  if (m_flavs[0].IsGluon()||m_flavs[0].IsAnti())
    lc[1]=Flow::Counter();
  p_colours[0][0]=lc[0];
  p_colours[0][1]=lc[1];
  for (size_t i(2);i<m_flavs.size();++i)
    SetColour(i,true,lc);
  p_colours[1][0]=lc[1];
  p_colours[1][1]=lc[0];
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

