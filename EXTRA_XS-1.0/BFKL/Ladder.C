#include "Ladder.H"

#include "Phase_Space_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Flow.H"
#include "Data_Reader.H"
#include "MyStrStream.H"
#include "XS_Selector.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

double Factorial(const double &n) 
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

Ladder::Ladder(const size_t nin,const size_t nout,
	       const ATOOLS::Flavour *flavours,
	       const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
	       BEAM::Beam_Spectra_Handler *const beamhandler,
	       PDF::ISR_Handler *const isrhandler,
	       ATOOLS::Selector_Data *const selectordata,
	       XS_Model_Base *const model):
  XS_Group(nin,nout,flavours,scalescheme,kfactorscheme,
	   beamhandler,isrhandler,selectordata,model),
  p_sudakov(new BFKL_Sudakov()),
  m_ncols(0), m_sudmode(1), m_nfixed(std::numeric_limits<size_t>::max()),
  m_multimode(1), p_jf(NULL)
{ 
  m_name="BFKL_ME";
  Data_Reader read(" ",";","!","=");
  int ktscheme(1);
  if (!read.ReadFromFile(ktscheme,"BFKL_KT_SCHEME")) ktscheme=1;
  else msg_Info()<<METHOD<<"(): Set k_T-scheme "<<ktscheme<<".\n";
  p_sudakov->SetKTScheme(ktscheme);
  if (!read.ReadFromFile(m_nfixed,"BFKL_FIXED_MULTI"))
    m_nfixed=std::numeric_limits<size_t>::max();
  else msg_Info()<<METHOD<<"(): Set fixed multiplicity "<<m_nfixed<<".\n";
  if (!read.ReadFromFile(m_multimode,"BFKL_MULTI_MODE")) m_multimode=1;
  else msg_Info()<<METHOD<<"(): Set multiplicity selection mode "
		 <<m_multimode<<".\n";
  if (!read.ReadFromFile(m_sudmode,"BFKL_SUDAKOV_MODE")) m_sudmode=1;
  else msg_Info()<<METHOD<<"(): Set Sudakov mode "<<m_sudmode<<".\n";
  if (!read.ReadFromFile(m_splitmode,"BFKL_SPLIT_MODE")) m_splitmode=(1<<6)-1;
  else msg_Info()<<METHOD<<"(): Set splitting mode "<<m_splitmode<<".\n";
  p_sudakov->SetSplitMode(m_splitmode);
  if (m_multimode==0) m_splitmode=1;
  if (!read.ReadFromFile(m_memode,"BFKL_ME_MODE")) m_memode=0;
  else msg_Info()<<METHOD<<"(): Set me correction mode "<<m_memode<<".\n";
  if (!read.ReadFromFile(m_pnmode,"BFKL_PN_MODE")) m_pnmode=1;
  else msg_Info()<<METHOD<<"(): Set process name mode "<<m_pnmode<<".\n";
  if (!read.ReadFromFile(m_nmaxme,"BFKL_NMAX_ME")) m_nmaxme=4;
  else msg_Info()<<METHOD<<"(): Set nmax me "<<m_nmaxme<<".\n";
  if (!read.ReadFromFile(m_fname,"BFKL_PROCESS_NAME")) m_fname="";
  else msg_Info()<<METHOD<<"(): Set process name '"<<m_fname<<"'.\n";
  if (m_memode&1) m_nmaxme=2;
  m_asecms=(*MODEL::as)(sqr(rpa.gen.Ecms()));
}

Ladder::~Ladder()
{
  if (p_jf!=NULL) delete p_jf;
  for (std::map<std::string,XS_Base*>::const_iterator
	 mit(m_xsmes.begin());mit!=m_xsmes.end();++mit) delete mit->second;
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
    if (m_memode&4) p_pdfs[i]->SetOrderingMode(0);
    if (m_splitmode<2) p_pdfs[i]->SetSplitMode(0);
    p_pdfs[i]->SetWeightMode(0);
  }
  delete [] p_momenta;
  delete [] p_flavours;
  p_momenta=&m_moms.front();
  p_flavours=&m_flavs.front();
  p_addflavours = new ATOOLS::Flavour[4];
  p_addmomenta = new ATOOLS::Vec4D[2];
  m_naddout=2;
  if (m_memode>0) 
    p_jf = new Jet_Finder(ToString(m_kt2min/sqr(rpa.gen.Ecms())),4);
  p_sudakov->SetKT2Min(m_kt2min);
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

double Ladder::KT2Max(const double &y,const double &sab,
		      const double &s1,const double &sn) const
{
  return 0.25*sqr(sab+s1-sn)/(sab*sqr(cosh(y)))-s1;
}

double Ladder::YMax(const double &kt2,const double &sab,
		    const double &s1,const double &sn) const
{
  double pzr(sqrt(1.0-4.0*(s1+kt2)*sab/sqr(sab+s1-sn)));
  return 0.5*log((1.0+pzr)/(1.0-pzr));
}

bool Ladder::GeneratePDFJet()
{
  double rn[4];
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  double ymax(YMax(m_kt2min,m_q2));
  // dice y
  m_y=m_y1=rn[0]*2.0*ymax-ymax;
  m_weight*=2.0*ymax;
  double kt2max(KT2Max(m_y,m_q2));
  // dice kt2
  double qt2(m_q2), b0(MODEL::as->Beta0(qt2)/M_PI);
  double l2(qt2*exp(-1.0/(b0*(*MODEL::as)(qt2))));
  m_kt2=m_kt21=l2*exp(p_sudakov->DicePolynomial
		      (log(m_kt2min/l2),log(kt2max/l2),-1.0,rn[1]));
  m_weight*=m_kt21*p_sudakov->WeightPolynomial
    (log(m_kt2min/l2),log(kt2max/l2),-1.0,log(m_kt21/l2));
  // dice phi
  double phi1(2.0*M_PI*rn[2]);
  m_weight*=2.0*M_PI;
  // construct jet
  double kt1(sqrt(m_kt21));
  m_k1=Vec4D(kt1*cosh(m_y1),kt1*cos(phi1),
	     kt1*sin(phi1),kt1*sinh(m_y1));
  m_q=m_k1;
  // dice y
  m_yn=rn[3]*2.0*ymax-ymax;
  m_weight*=2.0*ymax;
  // dice first flavour
  m_props.back()=kf::gluon;
  if (m_splitmode>1) {
    // select t-channel quark w/ probability T_R/2C_A
    double orn(ran.Get()), ehf(p_sudakov->Nf()*0.5/3.0);
    double x[3]={0.0,1.0-ehf/(ehf+1.0),1.0};
    int iv((int)Min(2.0*orn,1.0));
    double rn(x[iv]+(orn-iv)*(x[iv+1]-x[iv]));
    double nf(p_sudakov->Nf()), pw(1.0+2.0*nf), enf(rn*pw);
    if (enf>1.0) {
      m_props.back()=Flavour(kf::code(0.5+0.5*enf));
      if (int(1.0+enf)%2==1) m_props.back()=m_props.back().Bar();
    }
    m_weight*=pw*(x[iv+1]-x[iv]);
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
  double yn(m_yn);
  double mnt(sqrt(m_q.PPerp2()+sqr(m_props.back().Mass())));
  m_kn=Vec4D(mnt*cosh(yn),-m_q[1],-m_q[2],mnt*sinh(yn));
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

bool Ladder::BoostToCMS()
{
  Vec4D pa(m_pa-m_q);
  m_cms=Poincare(pa+m_pb);
  m_cms.Boost(pa);
  m_zaxis=Poincare(pa,Vec4D::ZVEC);
  m_sab=(pa+m_pb).Abs2();
  return m_sab>0.0;
}

bool Ladder::ConstructRung(const bool in)
{
  double kt(sqrt(m_kt2)), mt(sqrt(m_kt2+sqr(m_flavs.back().Mass())));
  Vec4D k(mt*cosh(m_y),kt*cos(m_phi),kt*sin(m_phi),mt*sinh(m_y));
  m_q+=m_moms.back()=k;
  if (k.Nan()) return false;
  return in?ConstructIncoming():true;
}

bool Ladder::BoostToLab()
{
  m_zaxis.RotateBack(m_moms.back());
  m_cms.BoostBack(m_moms.back());
  m_q=m_oldq+m_moms.back();
  m_y=m_moms.back().Y();
  m_kt2=m_moms.back().PPerp2();
  return true;
}

bool Ladder::DiceFixedMulti()
{
  /* store old values */
  m_oldq=m_q;
  double lasty(m_y);
  /* boost and rotate into new cm frame */
  BoostToCMS();
  double kt2max(KT2Max(0.0,m_sab));
  /* dice kt2, ps weight */
  m_kt2=p_sudakov->DicePolynomial(0.0,kt2max,-0.99,ran.Get());
  m_weight*=p_sudakov->WeightPolynomial(0.0,kt2max,-0.99,m_kt2);
  double ymax(YMax(m_kt2,m_sab));
  /* dice y, ps weight */
  m_y=ran.Get()*2.0*ymax-ymax;
  m_weight*=2.0*ymax;
  /* dice phi, ps weight is 1 */
  m_phi=ran.Get()*2.0*M_PI;
  /* add new rung */
  m_moms.push_back(Vec4D());
  m_flavs.push_back(kf::gluon);
  m_props.push_back(kf::gluon);
  if (!ConstructRung() || !BoostToLab() || 
      !ConstructIncoming()) return false;
  if (m_kt2<m_kt2min || m_q.PPerp2()<m_kt2min) return false;
  if (m_y1>m_yn) {
    if (m_yn>m_y || m_y>lasty) return false;
  }
  else {
    if (m_yn<m_y || m_y<lasty) return false;
  }
  /* alpha_s & splitting weight */
  m_weight*=(*MODEL::as)(m_kt2)/m_kt2*3.0/M_PI;
  /* sudakov weight */
  if (m_sudmode>0)
    m_weight*=exp(-(*MODEL::as)(m_oldq.PPerp2())*3.0/M_PI
  		  *log(m_oldq.PPerp2()/m_kt2min)*dabs(lasty-m_y));
  m_oldq=m_q;
  return true;
}

int Ladder::DiceVariableMulti()
{
  /* store old data */
  m_oldq=m_q;
  /* boost and rotate into new cm frame */
  BoostToCMS();
  double kt2max(KT2Max(0.0,m_sab));
  /* set up sudakov for new branching */
  p_sudakov->SetKT2Max(kt2max);
  p_sudakov->SetInMomentum(m_q);
  p_sudakov->SetInFlavour(m_props.back());
  /* dice y */
  if (!p_sudakov->Dice()) return 0;
  m_y=p_sudakov->GetY();
  m_kt2=p_sudakov->GetKT2();
  m_phi=p_sudakov->GetPhi();
  Vec4D k1(m_moms.size()>2?m_moms.back():m_k1);
  /* add new rung */
  m_moms.push_back(Vec4D());
  m_flavs.push_back(p_sudakov->Selected()->GetC());
  m_props.push_back(p_sudakov->Selected()->GetB());
  if (!ConstructRung() || !ConstructIncoming()) return -1;
  if (!p_sudakov->Approve(k1,m_oldq,m_moms.back(),m_q)) return -1;
  m_weight*=p_sudakov->GetWeight();
  /* add splitting weight for me correction */
  m_splitweight*=p_sudakov->GetSplitWeight();
  return 1;
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
    /* initialize sudakov */
    p_sudakov->SetYA(m_y1);
    p_sudakov->SetYB(m_yn);
    p_sudakov->SetKT2Max(KT2Max(0.0,m_sab));
    p_sudakov->SetInMomentum(m_q);
    p_sudakov->SetInFlavour(m_props.back());
    p_sudakov->Init();
    size_t cnt(0);
    /* iterate */
    int stat(1);
    while ((stat=DiceVariableMulti())>0) {
      if (cnt>m_nfixed) return false;
      ++cnt;
    }
    if (stat<0) return false;
    if (m_nfixed<std::numeric_limits<size_t>::max() &&
	cnt!=m_nfixed) return false;
  }
  else {
    for (size_t i(0);i<m_nfixed;++i)
      if (!DiceFixedMulti()) return false;
    /* last rung sudakov weight */
    if (m_sudmode>0)
      m_weight*=exp(-(*MODEL::as)(m_q.PPerp2())*3.0/M_PI
 		    *log(m_q.PPerp2()/m_kt2min)*dabs(m_y-m_yn));
  }
  if (m_kn.PPerp2()<m_kt2min) return false;
  /* fix inital state */
  m_moms[0]=m_p1-m_k1;
  m_moms[1]=m_p2-m_kn;
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
  m_splitweight=1.0;
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
    msg_Error()<<METHOD<<"(..): Invalid scales. Set weight 0."<<std::endl;
    return 0.0;
  }
  // calculate pdfs
  double fac(ATOOLS::rpa.gen.FactorizationScaleFactor());
  int ids[2]={m_y1>m_yn?0:1,m_y1>m_yn?1:0};
  p_pdfs[ids[0]]->Calculate(m_a1,m_z1,m_kt21,fac*m_mu21);
  p_pdfs[ids[1]]->Calculate(m_an,m_zn,m_kt2n,fac*m_mu2n);
  m_weight*=p_pdfs[ids[0]]->GetXPDF(m_flavs[0])/m_a1*m_z1*m_kt21;
  m_weight*=p_pdfs[ids[1]]->GetXPDF(m_flavs[1])/m_an*m_zn*m_kt2n;
  if (m_weight==0.0) return 0.0;
  // dice pdf jet flavours
  if (!p_pdfs[ids[0]]->SelectJetFlavour
      (p_addflavours[2],m_fl1,ran.Get())) return false;
  if (!p_pdfs[ids[1]]->SelectJetFlavour
      (p_addflavours[3],m_fln,ran.Get())) return false;
  m_weight*=p_pdfs[0]->CorrectionWeight()*p_pdfs[1]->CorrectionWeight();
  // add ll approximation of gg->gg me
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
  double mass(m_fl1.PSMass()+m_fln.PSMass());
  if (m_k1[0]<m_fl1.PSMass() || m_kn[0]<m_fln.PSMass()) return 0.0;
  for (size_t i(2);i<2+m_nout;++i) {
    mass+=m_flavs[i].PSMass();
    if (m_moms[i][0]<m_flavs[i].PSMass()) return 0.0;
  }
  if (sqr(mass)>(m_p1+m_p2).Abs2()) return 0.0;
  mass=p_addflavours[2].PSMass()+p_addflavours[3].PSMass();
  if (sqr(mass)>(m_p1+m_p2).Abs2()) return 0.0;
  std::string name(p_addflavours[2].IDName()+"_"+
		   p_addflavours[3].IDName()+"_");
  name+="_"+m_fl1.IDName();
  for (size_t i(2);i<2+m_nout;++i) name+="_"+m_flavs[i].IDName();
  name+="_"+m_fln.IDName();
  msg_Debugging()<<METHOD<<"(): Generated '"<<name<<"'\n";
  if (m_fname!="" && name!=m_fname) return 0.0;
  m_name=name;
  if (m_memode>0 && m_nvector<=Min(m_nmaxme,(size_t)6)) 
    m_weight*=MECorrection();
  name="BFKL_2_0_j_j";
  if (m_pnmode&1) name+="{2_"+ToString(m_nout+2)+"_"+m_name+"}";
  msg_Debugging()<<METHOD<<"(): Hand out '"<<name<<"'\n";
  m_name=name;
  return dabs(m_weight);
}

double Ladder::MECorrection()
{
  // check jet measures
  for (size_t i(2);i<m_moms.size();++i) {
    if (p_jf->MTij2(m_k1,m_moms[i],m_fl1.Mass(),
		    m_flavs[i].Mass())<m_kt2min ||
	p_jf->MTij2(m_kn,m_moms[i],m_fln.Mass(),
		    m_flavs[i].Mass())<m_kt2min) return 0.0;
    for (size_t j(i+1);j<m_moms.size();++j) 
      if (p_jf->MTij2(m_moms[i],m_moms[j],m_flavs[i].Mass(),
		      m_flavs[j].Mass())<m_kt2min) return 0.0;
  }
  if (p_jf->MTij2(m_k1,m_kn,m_fl1.Mass(),m_fln.Mass())<m_kt2min) return 0.0;
  int ids[2]={m_y1>m_yn?0:1,m_y1>m_yn?1:0};
  double sf1(p_pdfs[ids[0]]->JetKernel()->Value(m_z1));
  double sfn(p_pdfs[ids[1]]->JetKernel()->Value(m_zn));
  m_bfkldxs=sqr(4.0*M_PI*m_asecms)/8.0*sf1*sfn;
  if (m_nvector>2 && (m_memode&1)) 
    THROW(fatal_error,"Analytic reweighting not implemented");
  if (m_nvector<=2 && (m_memode&1)) {
    XS_Base *me(NULL);
    double s((m_k1+m_kn).Abs2());
    double t((m_k1-m_p1).Abs2());
    double u((m_k1-m_p2).Abs2());
    std::string name(m_name);
    std::map<std::string,XS_Base*>::iterator mit(m_xsmes.find(name));
    if (mit!=m_xsmes.end()) {
      me=mit->second;
    }
    else {
      std::vector<Flavour> fl(4);
      fl.front()=p_addflavours[2];
      fl[1]=p_addflavours[3];
      fl[2]=m_fl1;
      fl.back()=m_fln;
      msg_Debugging()<<METHOD<<"(): Initialize '"<<name<<"'\n";
      XS_Selector select(this);
      m_xsmes[name] = me = select.GetXS
	(m_nin,m_naddout,&fl.front(),false,0,2,false);
      if (me==NULL) {
	msg_Error()<<METHOD<<"(): Invalid Process '"<<name<<"'."<<std::endl;
	abort();
      }
    }
    m_dxs=(*me)(s,t,u);
    msg_Debugging()<<METHOD<<"(): \\sigma/\\sigma_{BFKL} = "<<m_dxs
		   <<" / ( "<<m_bfkldxs<<" * "<<m_splitweight
		   <<" ) = "<<m_dxs/(m_bfkldxs*m_splitweight)<<"\n";
    if (m_memode&8) m_dxs=1.0;
    return m_dxs/(m_bfkldxs*m_splitweight);
  }
  std::vector<ATOOLS::Vec4D> moms(m_moms.size()+2);
  for (size_t i(2);i<m_moms.size();++i) {
    moms[i+1]=m_moms[i];
    m_bfkldxs*=16.0*M_PI*m_asecms*3.0/m_moms[i].PPerp2();
  }
  if (m_memode&8) return 1.0/m_bfkldxs;
  return 0.0;
}

std::string Ladder::MapProcess(const std::vector<ATOOLS::Flavour> &fl)
{
  std::vector<Flavour> mapfl(fl.size());
  std::map<Flavour,Flavour> flmap;
  int ckf(1);
  for (size_t i(0);i<fl.size();++i) {
    if (flmap.find(fl[i])==flmap.end()) {
      if (fl[i].IsQuark()) {
	Flavour qr(fl[i].Kfcode());
	flmap[qr]=Flavour((kf::code)ckf++);
	flmap[qr.Bar()]=flmap[qr].Bar();
      }
      else {
	flmap[fl[i]]=fl[i];
      }
    }
    mapfl[i]=flmap[fl[i]];
  }
  std::string name(mapfl[0].IDName()+"_"+mapfl[1].IDName()+"__");
  for (size_t i(2);i<mapfl.size();++i) name+="_"+mapfl[i].IDName();
  msg_Debugging()<<METHOD<<"(): Map '"<<m_name<<"' -> '"<<name<<"'\n";
  return name;
}

void Ladder::SetColour(const size_t &i,const bool fw,int *const lc)
{
  Flavour q(m_props[i-2]), f(m_flavs[i]);
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

