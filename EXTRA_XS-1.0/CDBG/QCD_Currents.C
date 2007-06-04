#include "QCD_Currents.H"

#include "Message.H"
#include "Exception.H"
#include "STL_Tools.H"

using namespace EXTRAXS;
using namespace ATOOLS;

const double sqrttwo(sqrt(2.0));

CVec4D QCD_G::EM(const Vec4D &p,const int cr,const int ca)
{
  Spinor pp(1,p);
  CVec4D e(cr,ca,1<<m_key,0);
  e[0]=pp.U1()*m_km.U1()+pp.U2()*m_km.U2();
  e[3]=pp.U1()*m_km.U1()-pp.U2()*m_km.U2();
  e[1]=pp.U1()*m_km.U2()+pp.U2()*m_km.U1();
  e[2]=Complex(0.0,1.0)*(pp.U1()*m_km.U2()-pp.U2()*m_km.U1());
  return e/(sqrttwo*std::conj(m_kp*pp));
}

CVec4D QCD_G::EP(const Vec4D &p,const int cr,const int ca)
{
  Spinor pm(-1,p);
  CVec4D e(cr,ca,0,1<<m_key);
  e[0]=pm.U1()*m_kp.U1()+pm.U2()*m_kp.U2();
  e[3]=pm.U1()*m_kp.U1()-pm.U2()*m_kp.U2();
  e[1]=pm.U1()*m_kp.U2()+pm.U2()*m_kp.U1();
  e[2]=Complex(0.0,-1.0)*(pm.U1()*m_kp.U2()-pm.U2()*m_kp.U1());
  return e/(sqrttwo*std::conj(m_km*pm));
}

void QCD_G::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		       const int cr,const int ca)
{
  m_p=p;
  ResetJ();
  if (ch>=0) {
    CVec4D j(EP(p,cr,ca));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): '+' "<<m_id<<" "<<j<<"\n";
#endif
    AddJ(j);
  }
  if (ch<=0) {
    CVec4D j(EM(p,cr,ca));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): '-' "<<m_id<<" "<<j<<"\n";
#endif
    AddJ(j);
  }
}

void QCD_G::SetGauge(const ATOOLS::Vec4D &k)
{
  m_k=k;
  m_kp=Spinor(1,m_k);
  m_km=Spinor(-1,m_k);
}

void QCD_G::AddPropagator()
{
  // add propagator for off-shell leg
  double prop(1.0/m_p.Abs2());
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop<<"\n";
#endif
  for (CVec_Vector::iterator jit(m_j.begin());
       jit!=m_j.end();++jit) *jit*=prop;
}

Complex QCD_G::Contract(const Current_Base &c,
			const long unsigned int &hm,
			const long unsigned int &hp) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"("<<hm<<","<<hp<<"): {\n";
  msg_Indent();
#endif
  if (c.Type()!='V') THROW(fatal_error,"Invalid current type.");
  Complex res(0.0);
  for (CVec_Vector::const_iterator jit1(m_j.begin());
       jit1!=m_j.end();++jit1)
    for (CVec_Vector::const_iterator jit2(c.J<CVec4D>().begin());
	 jit2!=c.J<CVec4D>().end();++jit2) 
      if ((*jit1)(0)==(*jit2)(1) && (*jit1)(1)==(*jit2)(0) &&
	  jit1->H(0)+jit2->H(0)==hm && jit1->H(1)+jit2->H(1)==hp) {
#ifdef DEBUG__BG
	msg_Debugging()<<"add "<<*jit1**jit2<<"\n";
	msg_Debugging()<<"    "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	res+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return res;
}

char QCD_G::Type() const
{
  return 'V';
}

std::string QCD_G::CLabel() const
{
  return "gluon,label.side=right,label.dist=1.5curly_len,label=$g$";
}

DECLARE_GETTER(QCD_G_Getter,"QCD_V",Current_Base,Current_Key);

Current_Base *QCD_G_Getter::operator()(const Current_Key &key) const
{ 
  if (!key.m_fl.IsGluon()) return NULL;
  return new QCD_G(); 
}

void QCD_G_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd gluon current"; 
}

void QCD_QGC::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			 const int cr,const int ca)
{
}

void QCD_QGC::SetGauge(const ATOOLS::Vec4D &k)
{
}

void QCD_QGC::AddPropagator()
{
}

Complex QCD_QGC::Contract(const Current_Base &c,
			  const long unsigned int &hm,
			  const long unsigned int &hp) const
{
  THROW(fatal_error,"Multiplication of tensor particles not allowed");
  return 0.0;
}

char QCD_QGC::Type() const
{
  return 'T';
}

std::string QCD_QGC::CLabel() const
{
  return "dashes,label.side=right,label=$G$";
}

DECLARE_GETTER(QCD_QGC_Getter,"QCD_T",Current_Base,Current_Key);

Current_Base *QCD_QGC_Getter::operator()(const Current_Key &key) const
{ 
  if (key.m_fl.Kfcode()!=kf::gluonqgc) return NULL;
  return new QCD_QGC(); 
}

void QCD_QGC_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd tensor current"; 
}

void QCD_Q::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		       const int cr,const int ca)
{
  m_p=p;
  ResetJ();
  if (ch>=0) {
    CSpinor j(m_fl.IsAnti()?CSpinor(-1,1,1,p,ca,0,1<<m_key):
	      CSpinor(1,-1,1,p,cr,0,1<<m_key));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): '+' "<<m_id<<" "<<j
		   <<", m = "<<p.Mass()<<"\n";
#endif
    AddJ(j);
  }
  if (ch<=0) {
    CSpinor j(m_fl.IsAnti()?CSpinor(-1,1,-1,p,ca,1<<m_key,0):
	      CSpinor(1,-1,-1,p,cr,1<<m_key,0));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): '-' "<<m_id<<" "<<j
		   <<", m = "<<p.Mass()<<"\n";
#endif
    AddJ(j);
  }
}

void QCD_Q::SetGauge(const ATOOLS::Vec4D &k)
{
}

void QCD_Q::AddPropagator()
{
  // add propagator for off-shell leg
  double m(m_fl.Mass()), prop(1.0/(m_p.Abs2()-m*m));
  CSpinor hs;
  Complex pp(hs.PPlus(m_p)), pm(hs.PMinus(m_p));
  Complex pt(hs.PT(m_p)), ptc(hs.PTC(m_p));
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop
		 <<" <- p^2 = "<<m_p.Abs2()<<", m = "<<m<<"\n";
  msg_Debugging()<<"pp = "<<pp<<", pm = "<<pm<<", pt = "<<pt<<"\n";
#endif
  if (m_fl.IsAnti()) {
    for (CSpinor_Vector::iterator jit(m_j.begin());
	 jit!=m_j.end();++jit) {
      CSpinor j(jit->R(),jit->B(),(*jit)(),jit->H(0),jit->H(1));
      j[0]=pm*(*jit)[2]-ptc*(*jit)[3];
      j[1]=-pt*(*jit)[2]+pp*(*jit)[3];
      j[2]=pp*(*jit)[0]+ptc*(*jit)[1];
      j[3]=pt*(*jit)[0]+pm*(*jit)[1];
      *jit=(m==0.0?j:j-*jit*m)*prop;
    }
  }
  else {
    for (CSpinor_Vector::iterator jit(m_j.begin());
	 jit!=m_j.end();++jit) {
      CSpinor j(jit->R(),jit->B(),(*jit)(),jit->H(0),jit->H(1));
      j[0]=(*jit)[2]*pp+(*jit)[3]*pt;
      j[1]=(*jit)[2]*ptc+(*jit)[3]*pm;
      j[2]=(*jit)[0]*pm-(*jit)[1]*pt;
      j[3]=-(*jit)[0]*ptc+(*jit)[1]*pp;
      *jit=-(m==0.0?j:j+*jit*m)*prop;
    }
  }
}

Complex QCD_Q::Contract(const Current_Base &c,
			const long unsigned int &hm,
			const long unsigned int &hp) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"("<<hm<<","<<hp<<"): {\n";
  msg_Indent();
#endif
  if (c.Type()!='S') THROW(fatal_error,"Invalid current type.");
  Complex res(0.0);
  for (CSpinor_Vector::const_iterator jit1(m_j.begin());
       jit1!=m_j.end();++jit1)
    for (CSpinor_Vector::const_iterator jit2(c.J<CSpinor>().begin());
	 jit2!=c.J<CSpinor>().end();++jit2) 
      if ((*jit1)()==(*jit2)() &&
	  jit1->H(0)+jit2->H(0)==hm && jit1->H(1)+jit2->H(1)==hp) {
#ifdef DEBUG__BG
	msg_Debugging()<<"add "<<*jit1**jit2<<"\n";
	msg_Debugging()<<"    "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	res+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return res;
}

char QCD_Q::Type() const
{
  return 'S';
}

std::string QCD_Q::CLabel() const
{
  return "fermion,label.side=right,label=$"+m_fl.TexName()+"$";
}

DECLARE_GETTER(QCD_Q_Getter,"QCD_S",Current_Base,Current_Key);

Current_Base *QCD_Q_Getter::operator()(const Current_Key &key) const
{ 
  if (!key.m_fl.IsQuark()) return NULL;
  return new QCD_Q(key.m_fl); 
}

void QCD_Q_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd quark current"; 
}

