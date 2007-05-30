#include "QCD_Vertices.H"

#include "Running_AlphaS.H"
#include "Run_Parameter.H"
#include "STL_Tools.H"

using namespace EXTRAXS;
using namespace ATOOLS;

const double sqrttwo(sqrt(2.0));
const double invsqrttwo(1.0/sqrttwo);

QCD_GGG::QCD_GGG(Current_Base *const c):
  Vertex(c),
  m_cpl(sqrt(4.0*M_PI*(*MODEL::as)(sqr(rpa.gen.Ecms())))) {}

void QCD_GGG::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<"\n";
  msg_Indent();
#endif
  // triple gluon vertex
  for (CVec_Vector::const_iterator ait(p_a->J<CVec4D>().begin());
       ait!=p_a->J<CVec4D>().end();++ait) {
    for (CVec_Vector::const_iterator bit(p_b->J<CVec4D>().begin());
	 bit!=p_b->J<CVec4D>().end();++bit) {
      if ((*ait)(0)==(*bit)(1) || (*ait)(1)==(*bit)(0)) {
	// sum V_3(\pi_1,\pi_2) & V_3(\pi_2,\pi_1)
	CVec4D j(invsqrttwo*(*ait**bit)*(p_a->P()-p_b->P())
		 +sqrttwo*((*ait*p_b->P())**bit-(*bit*p_a->P())**ait));
	j.SetH(ait->H(0)+bit->H(0),0);
	j.SetH(ait->H(1)+bit->H(1),1);
	if ((*ait)(0)==(*bit)(1)) {
	  if ((*ait)(1)==(*bit)(0) && (*ait)(0)==(*ait)(1))
	    // color singlet, yields zero
	    continue;
	  // case 1: original ordering, V_3(\pi_1,\pi_2)
	  j(0)=(*bit)(0);
	  j(1)=(*ait)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"'+' "<<*ait<<"\n";
	  msg_Debugging()<<"    "<<*bit<<"\n";
#endif
	  AddJ(m_cpl*j);
	}
	if ((*ait)(1)==(*bit)(0)) {
	  // case 2: reverse ordering, V_3(\pi_2,\pi_1)
	  j(0)=(*ait)(0);
	  j(1)=(*bit)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"'-' "<<*ait<<"\n";
	  msg_Debugging()<<"    "<<*bit<<"\n";
#endif
	  AddJ(-m_cpl*j);
	}
      }
    }
  }
}

DECLARE_GETTER(QCD_GGG_Getter,"QCD_VVV",Vertex,Vertex_Key);

Vertex *QCD_GGG_Getter::operator()(const Vertex_Key &key) const
{ 
  if (!key.p_a->Flav().IsGluon()||!key.p_b->Flav().IsGluon()||
      !key.p_c->Flav().IsGluon()) return NULL;
  return new QCD_GGG(key.p_c); 
}

void QCD_GGG_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd triple gluon vertex"; 
}

QCD_GGT::QCD_GGT(Current_Base *const c):
  Vertex(c),
  m_cpl(sqrt(4.0*M_PI*(*MODEL::as)(sqr(rpa.gen.Ecms())))) {}

void QCD_GGT::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<"\n";
  msg_Indent();
#endif
  for (CVec_Vector::const_iterator jit1(p_a->J<CVec4D>().begin());
       jit1!=p_a->J<CVec4D>().end();++jit1) 
    for (CVec_Vector::const_iterator jit2(p_b->J<CVec4D>().begin());
	 jit2!=p_b->J<CVec4D>().end();++jit2) {
      // sum T_3(\pi_1,\pi_2) & T_3(\pi_2,\pi_1)
      if ((*jit1)(0)==(*jit2)(1)) {
	if ((*jit1)(1)==(*jit2)(0) && (*jit1)(0)==(*jit1)(1))
	  // color singlet, yields zero
	  continue;
#ifdef DEBUG__BG
	msg_Debugging()<<"'+' "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	AddJ(m_cpl*CAsT4D(*jit1,*jit2));
      }
      if ((*jit1)(1)==(*jit2)(0)) {
#ifdef DEBUG__BG
	msg_Debugging()<<"'-' "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	AddJ(m_cpl*CAsT4D(*jit2,*jit1));
      }
    }
}

DECLARE_GETTER(QCD_GGT_Getter,"QCD_VVT",Vertex,Vertex_Key);

Vertex *QCD_GGT_Getter::operator()(const Vertex_Key &key) const
{ 
  if (!key.p_a->Flav().IsGluon()||!key.p_b->Flav().IsGluon()||
      key.p_c->Flav().Kfcode()!=kf::gluonqgc) return NULL;
  return new QCD_GGT(key.p_c); 
}

void QCD_GGT_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd tensor production vertex"; 
}

QCD_GTG::QCD_GTG(Current_Base *const c):
  Vertex(c),
  m_cpl(sqrt(4.0*M_PI*(*MODEL::as)(sqr(rpa.gen.Ecms())))) {}

void QCD_GTG::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<"\n";
  msg_Indent();
#endif
  // gluon/tensor vertex
  for (CVec_Vector::const_iterator ait(p_a->J<CVec4D>().begin());
       ait!=p_a->J<CVec4D>().end();++ait) {
    for (CAsT_Vector::const_iterator bit(p_b->J<CAsT4D>().begin());
	 bit!=p_b->J<CAsT4D>().end();++bit) {
      if ((*ait)(0)==(*bit)(1) || (*ait)(1)==(*bit)(0)) {
	// sum V_4(\pi_1,\{\pi_2,\pi_3\}) & V_4(\{\pi_2,\pi_3\},\pi_1)
	CVec4D j(*ait**bit/2.0);
	j.SetH(ait->H(0)+bit->H(0),0);
	j.SetH(ait->H(1)+bit->H(1),1);
	if ((*ait)(0)==(*bit)(1)) {
	  if ((*ait)(1)==(*bit)(0) && (*ait)(0)==(*ait)(1))
	    // color singlet, yields zero
	    continue;
	  // case 1: original ordering, V_4(\pi_1,\{\pi_2,\pi_3\})
	  j(0)=(*bit)(0);
	  j(1)=(*ait)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"'+' "<<*ait<<"\n";
	  msg_Debugging()<<"    "<<*bit<<"\n";
#endif
	  AddJ(m_cpl*j);
	}
	if ((*ait)(1)==(*bit)(0)) {
	  // case 2: reverse ordering, V_4(\{\pi_2,\pi_3\},\pi_1)
	  j(0)=(*ait)(0);
	  j(1)=(*bit)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"'-' "<<*ait<<"\n";
	  msg_Debugging()<<"    "<<*bit<<"\n";
#endif
	  AddJ(-m_cpl*j);
	}
      }
    }
  }
}

DECLARE_GETTER(QCD_GTG_Getter,"QCD_VTV",Vertex,Vertex_Key);

Vertex *QCD_GTG_Getter::operator()(const Vertex_Key &key) const
{ 
  if (!key.p_a->Flav().IsGluon()||!key.p_c->Flav().IsGluon()||
      key.p_b->Flav().Kfcode()!=kf::gluonqgc) return NULL;
  return new QCD_GTG(key.p_c); 
}

void QCD_GTG_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd tensor decay vertex"; 
}

QCD_QQG::QCD_QQG(Current_Base *const c):
  Vertex(c),
  m_cpl(sqrt(4.0*M_PI*(*MODEL::as)(sqr(rpa.gen.Ecms())))) {}

void QCD_QQG::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<"\n";
  msg_Indent();
#endif
  for (CSpinor_Vector::const_iterator ait(p_a->J<CSpinor>().begin());
       ait!=p_a->J<CSpinor>().end();++ait) {
    for (CSpinor_Vector::const_iterator bit(p_b->J<CSpinor>().begin());
	 bit!=p_b->J<CSpinor>().end();++bit) {
      bool singlet((*ait)()==(*bit)());
      Complex j01((*ait)[0]*(*bit)[2]+(*ait)[3]*(*bit)[1]);
      Complex j02((*ait)[1]*(*bit)[3]+(*ait)[2]*(*bit)[0]);
      Complex j11((*ait)[0]*(*bit)[3]-(*ait)[2]*(*bit)[1]);
      Complex j12((*ait)[1]*(*bit)[2]-(*ait)[3]*(*bit)[0]);
      CVec4D j(j01+j02,j11+j12,Complex(0.0,-1.0)*(j11-j12),j01-j02);
      j.SetH(ait->H(0)+bit->H(0),0);
      j.SetH(ait->H(1)+bit->H(1),1);
      j(0)=(*ait)();
      j(1)=(*bit)();
#ifdef DEBUG__BG
      msg_Debugging()<<"'+' "<<*ait<<"\n";
      msg_Debugging()<<"    "<<*bit<<"\n";
#endif
      AddJ(j*=invsqrttwo*m_cpl);
      if (singlet) {
	j*=-1.0/3.0;
  	for (short unsigned int i(1);i<=3;++i) {
  	  j(0)=j(1)=i;
	  AddJ(j);
  	}
      }
    }
  }
}

DECLARE_GETTER(QCD_QQG_Getter,"QCD_SSV",Vertex,Vertex_Key);

Vertex *QCD_QQG_Getter::operator()(const Vertex_Key &key) const
{ 
  if (!key.p_a->Flav().IsQuark()||!key.p_b->Flav().IsQuark()||
      !key.p_b->Flav().IsAnti()||!key.p_c->Flav().IsGluon()||
      !(key.p_a->Flav()==key.p_b->Flav().Bar())) return NULL;
  return new QCD_QQG(key.p_c); 
}

void QCD_QQG_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd gluon production vertex"; 
}

QCD_QGQ::QCD_QGQ(Current_Base *const c):
  Vertex(c),
  m_cpl(sqrt(4.0*M_PI*(*MODEL::as)(sqr(rpa.gen.Ecms())))) {}

void QCD_QGQ::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<"\n";
  msg_Indent();
#endif
  bool anti(p_a->Flav().IsAnti());
  for (CSpinor_Vector::const_iterator ait(p_a->J<CSpinor>().begin());
       ait!=p_a->J<CSpinor>().end();++ait) {
    for (CVec_Vector::const_iterator bit(p_b->J<CVec4D>().begin());
	 bit!=p_b->J<CVec4D>().end();++bit) {
      bool singlet((*bit)(0)==(*bit)(1));
      if (!anti && ((*ait)()==(*bit)(1) || singlet)) {
	CSpinor j(ait->R(),ait->B(),(*bit)(0),
		  ait->H(0)+bit->H(0),ait->H(1)+bit->H(1));
	Complex jp(ait->PPlus(*bit)), jm(ait->PMinus(*bit));
	Complex jt(ait->PT(*bit)), jtc(ait->PTC(*bit));
	j[0]=(*ait)[2]*jp+(*ait)[3]*jt;
	j[1]=(*ait)[2]*jtc+(*ait)[3]*jm;
	j[2]=(*ait)[0]*jm-(*ait)[1]*jt;
	j[3]=-(*ait)[0]*jtc+(*ait)[1]*jp;
#ifdef DEBUG__BG
	msg_Debugging()<<"'+' "<<*ait<<"\n";
	msg_Debugging()<<"    "<<*bit<<"\n";
#endif
	j*=invsqrttwo*m_cpl;
	if ((*ait)()==(*bit)(1)) AddJ(j);
	if (singlet) {
	  j()=(*ait)();
	  AddJ(-j/3.0);
	}
      }
      else if (anti && ((*ait)()==(*bit)(0) || singlet)) {
	CSpinor j(ait->R(),ait->B(),(*bit)(1),
		  ait->H(0)+bit->H(0),ait->H(1)+bit->H(1));
	Complex jp(ait->PPlus(*bit)), jm(ait->PMinus(*bit));
	Complex jt(ait->PT(*bit)), jtc(ait->PTC(*bit));
	j[0]=(*ait)[2]*jm-(*ait)[3]*jtc;
	j[1]=-(*ait)[2]*jt+(*ait)[3]*jp;
	j[2]=(*ait)[0]*jp+(*ait)[1]*jtc;
	j[3]=(*ait)[0]*jt+(*ait)[1]*jm;
#ifdef DEBUG__BG
	msg_Debugging()<<"'-' "<<*ait<<"\n";
	msg_Debugging()<<"    "<<*bit<<"\n";
#endif
	j*=invsqrttwo*m_cpl;
	if ((*ait)()==(*bit)(0)) AddJ(j);
	if (singlet) {
	  j()=(*ait)();
	  AddJ(-j/3.0);
	}
      }
    }
  }
}

DECLARE_GETTER(QCD_QGQ_Getter,"QCD_SVS",Vertex,Vertex_Key);

Vertex *QCD_QGQ_Getter::operator()(const Vertex_Key &key) const
{ 
  if (!key.p_a->Flav().IsQuark()||!key.p_b->Flav().IsGluon()||
      !(key.p_c->Flav()==key.p_a->Flav())) return NULL;
  return new QCD_QGQ(key.p_c); 
}

void QCD_QGQ_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd gluon decay vertex"; 
}


