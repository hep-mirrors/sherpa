#include "BG_Amplitude.H"

#include "Message.H"
#include "Exception.H"
#include "STL_Tools.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

using namespace EXTRAXS;
using namespace ATOOLS;

static const double sqrttwo(sqrt(2.0));
static const double invsqrttwo(1.0/sqrttwo);

BG_Amplitude::BG_Amplitude():
  p_p(NULL), p_ps(NULL), p_ssp(NULL),
  p_ep(NULL), p_em(NULL), p_ssep(NULL), p_ssem(NULL), p_j(NULL),
  p_ch(NULL), p_id(NULL),
  m_n(0), m_nh(0)
{
  Vec4D k(1.0,1.0,0.0,0.0);
  m_kp=Spinor(1,k);
  m_km=Spinor(-1,k);
}

BG_Amplitude::~BG_Amplitude()
{
  CleanUp();
}

void BG_Amplitude::CleanUp()
{
  if (p_j!=NULL) {
    for (size_t i(0);i<m_n;++i) {
      for (size_t j(0);j<m_n;++j) delete p_j[i][j];
      delete p_j[i];
    }
    delete p_j;
  }
  if (p_ps!=NULL) {
    for (size_t i(0);i<m_n;++i) delete p_ps[i];
    delete p_ps;
  }
  if (p_p!=NULL) delete p_p;
  if (p_ssp!=NULL) delete p_ssp;
  if (p_ep!=NULL) delete p_ep;
  if (p_em!=NULL) delete p_em;
  if (p_ssep!=NULL) delete p_ssep;
  if (p_ssem!=NULL) delete p_ssem;
  if (p_ch!=NULL) delete p_ch;
  if (p_id!=NULL) delete p_id;
  m_nh=m_n=0;
  p_ssp=p_p=NULL;
  p_ssep=p_ep=NULL;
  p_ssem=p_em=NULL;
  p_ch=NULL;
  p_id=NULL;
  p_j=NULL;
  p_ps=NULL;
  m_hmap.clear();
  m_hamps.clear();
  m_maxid=0;
}

bool BG_Amplitude::Construct(const size_t &n)
{
  m_n=n;
  m_nh=1<<(m_n+1);
  p_p = new Vec4D[n];
  p_ssp = new Vec4D[n];
  p_em = new CVec4D[n];
  p_ep = new CVec4D[n];
  p_ssem = new CVec4D[n];
  p_ssep = new CVec4D[n];
  p_ch = new int[n];
  p_id = new int[n];
  p_j = new CVec4D**[n];
  p_ps = new Vec4D*[n];
  for (size_t i(0);i<n;++i) {
    p_j[i] = new CVec4D*[n];
    for (size_t j(0);j<n;++j) p_j[i][j] = new CVec4D[m_nh];
    p_ps[i] = new Vec4D[n];
  }
  return true;
}

CVec4D BG_Amplitude::EM(const Vec4D &p)
{
  Spinor pp(1,p);
  CVec4D e;
  e[0]=pp.U1()*m_km.U1()+pp.U2()*m_km.U2();
  e[3]=pp.U1()*m_km.U1()-pp.U2()*m_km.U2();
  e[1]=pp.U1()*m_km.U2()+pp.U2()*m_km.U1();
  e[2]=Complex(0.0,1.0)*(pp.U1()*m_km.U2()-pp.U2()*m_km.U1());
  return e/(sqrttwo*std::conj(m_kp*pp));
}

CVec4D BG_Amplitude::EP(const Vec4D &p)
{
  Spinor pm(-1,p);
  CVec4D e;
  e[0]=pm.U1()*m_kp.U1()+pm.U2()*m_kp.U2();
  e[3]=pm.U1()*m_kp.U1()-pm.U2()*m_kp.U2();
  e[1]=pm.U1()*m_kp.U2()+pm.U2()*m_kp.U1();
  e[2]=Complex(0.0,-1.0)*(pm.U1()*m_kp.U2()-pm.U2()*m_kp.U1());
  return e/(sqrttwo*std::conj(m_km*pm));
}

CVec4D BG_Amplitude::V3L(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2,
			 const ATOOLS::CVec4D &j1,const ATOOLS::CVec4D &j2)
{
  return invsqrttwo*(j1*j2)*(p1-p2)+sqrttwo*((j1*p2)*j2-(j2*p1)*j1);
}

CVec4D BG_Amplitude::V4L(const ATOOLS::CVec4D &j1,const ATOOLS::CVec4D &j2,
			 const ATOOLS::CVec4D &j3)
{
  return (j1*j3)*j2-0.5*((j2*j3)*j1+(j2*j1)*j3);
}

CVec4D BG_Amplitude::VL(const size_t &ei,const size_t &ej,const size_t &ih)
{
  p_ps[ei][ej]=p_ps[ei][ei]+p_ps[ei+1][ej];
  CVec4D sum;
  for (size_t i(ei);i<ej;++i) {
    size_t ihj(ih>>(i-ei+1)), ihi(ih-(ihj<<(i-ei+1)));
    if (i>ei) p_ps[ei][i]=p_ps[ei][ei]+p_ps[ei+1][i];
    if (ej>i+1) p_ps[i+1][ej]=p_ps[i+1][i+1]+p_ps[i+2][ej];
    sum+=V3L(p_ps[ei][i],p_ps[i+1][ej],p_j[ei][i][ihi],p_j[i+1][ej][ihj]);
    if (i<ej-1) {
      for (size_t j(i+1);j<ej;++j) {
	size_t ihjj(ihj>>(j-i)), ihji(ihj-(ihjj<<(j-i)));
	sum+=V4L(p_j[ei][i][ihi],p_j[i+1][j][ihji],p_j[j+1][ej][ihjj]);
      }
    }
  }
  return sum;
}

CVec4D BG_Amplitude::JL(const size_t &i,const size_t &j,const size_t &ih)
{
  if (j==i) return ih>0?p_ep[i]:p_em[i];
  if (i==1 && j==m_n-1) return VL(i,j,ih);
  return VL(i,j,ih)/p_ps[i][j].Abs2();
}

void BG_Amplitude::CalcJL()
{
  PROFILE_HERE;
  for (size_t d(0);d<m_n;++d) 
    for (size_t i(1);i<m_n;++i) {
      if (i+d<m_n) {
	size_t nh(1<<(d+1));
	for (size_t ih(0);ih<nh;++ih)
	  p_j[i][i+d][ih]=JL(i,i+d,ih);
      }
    }
  p_j[0][0][0]=JL(0,0,0);
  p_j[0][0][1]=JL(0,0,1);
}

void BG_Amplitude::SetMomenta(const Vec4D_Vector &moms)
{
  PROFILE_HERE;
  Vec4D sum;
  for (size_t i(0);i<m_n;++i) {
    p_ssp[i]=p_ps[i][i]=p_p[i]=moms[i];
    p_ssep[i]=p_ep[i]=EP(p_p[i]);
    p_ssem[i]=p_em[i]=p_ep[i].Conj();
    sum+=p_p[i];
  }
  static double accu(sqrt(Accu()));
  CVec4D::SetAccu(accu);
  if (!((CVec4D)sum).IsZero()) 
    msg.Error()<<METHOD<<"(): Four momentum not conserved. sum = "
	       <<sum<<"."<<std::endl;
  CVec4D::SetAccu(Accu());
}

size_t BG_Amplitude::MakeId(const Int_Vector &ids) const
{
  if (ids.size()!=m_n) 
    THROW(fatal_error,"Invalid particle number");
  size_t id(0);
  for (size_t i(0);i<ids.size();++i) 
    if (ids[i]>0) id+=1<<i;
  return id;
}

Int_Vector BG_Amplitude::MakeId(const size_t &id) const
{
  size_t ic(id);
  Int_Vector ids(m_n,-1);
  for (size_t i(0);i<ids.size();++i) {
    size_t c(1<<i);
    if (ic&c) {
      ids[i]=1;
      ic-=c;
    }
  }
  if (ic!=0) THROW(fatal_error,"Invalid particle number");
  return ids;
}

bool BG_Amplitude::EvaluateAll(const Int_Vector &perm)
{
  PROFILE_HERE;
  for (size_t j(0);j<m_n;++j) {
    p_ps[j][j]=p_p[j]=p_ssp[perm[j]];
    p_ep[j]=p_ssep[perm[j]];
    p_em[j]=p_ssem[perm[j]];
    for (size_t i(0);i<m_ress.size();++i)
      m_chirs[i][j]=m_schirs[i][perm[j]];
  }
  CalcJL();
  for (size_t i(0);i<m_ress.size();++i) {
    if (m_hmap[i]>=0) {
      m_ress[i]=std::conj(m_ress[m_hmap[i]]);
      continue;
    }
    size_t ih(MakeId(m_chirs[i])>>1);
    m_ress[i]=p_j[0][0][m_chirs[i][0]>0?1:0]*p_j[1][m_n-1][ih];
  }
  return true;
}

bool BG_Amplitude::CheckChirs(const Int_Vector &chirs)
{
  size_t p(0), m(0);
  for (size_t i(0);i<chirs.size();++i) {
    if (chirs[i]>0) ++p;
    else if (chirs[i]<0) ++m;
    else THROW(fatal_error,"Invalid helicities");
  }
  return p>1 && m>1;
}

bool BG_Amplitude::MapChirs(const Int_Vector &chirs)
{
  Int_Vector rchirs(chirs.size());
  for (size_t i(0);i<chirs.size();++i) rchirs[i]=-chirs[i];
  m_hmap.push_back(-1);
  for (size_t i(0);i<m_schirs.size();++i) {
    bool hit(true);
    for (size_t j(0);j<rchirs.size();++j)
      if (m_schirs[i][j]!=rchirs[j]) hit=false;
    if (!hit) continue;
    m_hmap.back()=i;
    return true;
  }
  return false;
}

bool BG_Amplitude::Construct(Int_Vector chirs,const size_t &i)
{
  if (i==chirs.size()) {
    if (CheckChirs(chirs)) {
      m_schirs.push_back(chirs);
      m_maxid=Max(m_maxid,MakeId(chirs));
      MapChirs(chirs);
      m_ress.push_back(0.0);
    }
    return true;
  }
  if (chirs[i]!=0) {
    Construct(chirs,i+1);
  }
  else {
    for (int ch(1);ch>=-1;ch-=2) {
      chirs[i]=ch;
      Construct(chirs,i+1);
    }
  }
  return true;
}

bool BG_Amplitude::Construct(Flavour_Vector flavs,Int_Vector chirs)
{
  CleanUp();
  if (!Construct(chirs.size())) return false;
  if (!Construct(chirs,0)) return false;
  m_hamps.resize(m_maxid+1,0);
  for (size_t i(0);i<m_schirs.size();++i)
    m_hamps[MakeId(m_schirs[i])]=1;
  m_chirs=m_schirs;
  return true;
}

bool BG_Amplitude::GaugeTest(const Vec4D_Vector &moms)
{
  msg_Info()<<METHOD<<"(): Performing gauge test ... "<<std::flush;
  Vec4D k(1.0,1.0,0.0,0.0);
  m_kp=Spinor(1,k);
  m_km=Spinor(-1,k);
  SetMomenta(moms);
  Int_Vector perm(m_n);
  for (size_t i(0);i<m_n;++i) perm[i]=i;
  if (!EvaluateAll(perm)) return false;
  std::vector<Complex> ress(m_ress);
  k=Vec4D(1.0,0.0,1.0,0.0);
  m_kp=Spinor(1,k);
  m_km=Spinor(-1,k);
  SetMomenta(moms);
  if (!EvaluateAll(perm)) return false;
  double mean(0.0);
  for (size_t i(0);i<m_ress.size();++i) mean+=std::abs(ress[i]);
  mean/=m_ress.size();
  for (size_t i(0);i<m_ress.size();++i) {
    msg_Debugging()<<"A("<<m_schirs[i]
		   <<") = "<<m_ress[i]<<" vs. "<<ress[i]<<" -> dev. "
		   <<m_ress[i].real()/ress[i].real()-1.0<<" "
		   <<m_ress[i].imag()/ress[i].imag()-1.0<<"\n";
    double accu(sqrt(Accu()));
    if (!IsEqual(m_ress[i].real(),ress[i].real(),accu) ||
	!IsEqual(m_ress[i].imag(),ress[i].imag(),accu)) {
      double rat(mean/std::abs(m_ress[i])*Accu());
      if (IsEqual(m_ress[i].real(),ress[i].real(),rat) &&
	  IsEqual(m_ress[i].imag(),ress[i].imag(),rat)) {
	msg.Error().precision(12);
	msg.Error()<<METHOD<<"(): Large deviation: "
		   <<m_ress[i]<<" vs. "<<ress[i]<<"\n  => ("
		   <<(m_ress[i].real()/ress[i].real()-1.0)<<","
		   <<(m_ress[i].imag()/ress[i].imag()-1.0)
		   <<") {"<<rat<<"}."<<std::endl;
	msg.Error().precision(6);
      }
      else {
	msg.Error().precision(12);
	msg.Error()<<METHOD<<"(): Gauge test failed. "
		   <<m_ress[i]<<" vs. "<<ress[i]<<"\n  => ("
		   <<(m_ress[i].real()/ress[i].real()-1.0)<<","
		   <<(m_ress[i].imag()/ress[i].imag()-1.0)
		   <<") {"<<rat<<"}."<<std::endl;
	msg.Error().precision(6);
	return false;
      }
    }
  }
  msg_Info()<<"satisfied."<<std::endl;
  return true;
}

