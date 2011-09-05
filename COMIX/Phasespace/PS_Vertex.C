#include "COMIX/Phasespace/PS_Vertex.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#ifdef USING__MPI
#include "mpi.h"
#endif

using namespace COMIX;
using namespace ATOOLS;

PS_Vertex::PS_Vertex(const Vertex_Key &key):
  Vertex_Base(key,0,0), m_alpha(1.0), m_oldalpha(1.0), m_weight(1.0),
  m_np(0.0), m_sum(0.0), m_sum2(0.0), m_max(0.0),
  m_mnp(0.0), m_msum(0.0), m_msum2(0.0) {}

void PS_Vertex::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<"\n";
  msg_Indent();
#endif
  int sca(p_a->Flav().StrongCharge());
  int scb(p_b->Flav().StrongCharge());
  if (sca==0) {
    const PS_Info_Vector &ca(p_a->J<PS_Info>());
    const PS_Info_Vector &cb(p_b->J<PS_Info>());
    for (PS_Info_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (PS_Info_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit) {
	AddJ(PS_Info((*bit)(0),(*bit)(1),ait->H()+bit->H()));
	m_zero=false;
      }
  }
  else if (scb==0) {
    const PS_Info_Vector &ca(p_a->J<PS_Info>());
    const PS_Info_Vector &cb(p_b->J<PS_Info>());
    for (PS_Info_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (PS_Info_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit) {
	AddJ(PS_Info((*ait)(0),(*ait)(1),ait->H()+bit->H()));
	m_zero=false;
      }
  }
  else if (abs(sca)==3 && abs(scb)==3) {
    const PS_Info_Vector &ca((sca<0?p_b:p_a)->J<PS_Info>());
    const PS_Info_Vector &cb((sca<0?p_a:p_b)->J<PS_Info>());
    for (PS_Info_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (PS_Info_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit) {
	if ((*ait)(0)!=(*bit)(1))
	  AddJ(PS_Info((*ait)(0),(*bit)(1),ait->H()+bit->H()));
	else
	  for (size_t i(this->s_cimin);i<=this->s_cimax;++i) 
	    AddJ(PS_Info(i,i,ait->H()+bit->H()));
	m_zero=false;
      }
  }
  else if (abs(sca)==3 || abs(scb)==3) {
    const PS_Info_Vector &ca((abs(scb)==3?p_b:p_a)->J<PS_Info>());
    const PS_Info_Vector &cb((abs(scb)==3?p_a:p_b)->J<PS_Info>());
    if ((abs(scb)==3?scb:sca)<0) {
      for (PS_Info_Vector::const_iterator 
	     ait(ca.begin());ait!=ca.end();++ait)
	for (PS_Info_Vector::const_iterator 
	       bit(cb.begin());bit!=cb.end();++bit) {
	  if ((*ait)(1)==(*bit)(0)) {
	    AddJ(PS_Info(0,(*bit)(1),ait->H()+bit->H()));
	    m_zero=false;
	  }
	  else if ((*bit)(0)==(*bit)(1)) {
	    AddJ(PS_Info(0,(*ait)(1),ait->H()+bit->H()));
	    m_zero=false;
	  }
	}
    }
    else {
      for (PS_Info_Vector::const_iterator 
	     ait(ca.begin());ait!=ca.end();++ait)
	for (PS_Info_Vector::const_iterator 
	       bit(cb.begin());bit!=cb.end();++bit) {
	  if ((*ait)(0)==(*bit)(1)) {
	    AddJ(PS_Info((*bit)(0),0,ait->H()+bit->H()));
	    m_zero=false;
	  }
	  else if ((*bit)(0)==(*bit)(1)) {
	    AddJ(PS_Info((*ait)(0),0,ait->H()+bit->H()));
	    m_zero=false;
	  }
	}
    }
  }
  else {
    const PS_Info_Vector &ca(p_a->J<PS_Info>());
    const PS_Info_Vector &cb(p_b->J<PS_Info>());
    for (PS_Info_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
	for (PS_Info_Vector::const_iterator 
	       bit(cb.begin());bit!=cb.end();++bit) {
	  if ((*ait)(1)==(*bit)(0)) {
	    if ((*ait)(0)==(*bit)(1))
	      for (size_t i(this->s_cimin);i<=this->s_cimax;++i) 
		AddJ(PS_Info(i,i,ait->H()+bit->H()));
	    else
	      AddJ(PS_Info((*ait)(0),(*bit)(1),ait->H()+bit->H()));
	    m_zero=false;
	  }
	  else if ((*ait)(0)==(*bit)(1)) {
	    AddJ(PS_Info((*bit)(0),(*ait)(1),ait->H()+bit->H()));
	    m_zero=false;
	  }
	}
  }
}

void PS_Vertex::MPISync()
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    int rank=MPI::COMM_WORLD.Get_rank();
    double val[4];
    if (rank==0) {
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Recv(&val,4,MPI::DOUBLE,MPI::ANY_SOURCE,tag);
	m_mnp+=val[0];
	m_msum+=val[1];
	m_msum2+=val[2];
	m_max=ATOOLS::Max(m_max,val[3]);
      }
      val[0]=m_mnp;
      val[1]=m_msum;
      val[2]=m_msum2;
      val[3]=m_max;
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Send(&val,4,MPI::DOUBLE,tag,size+tag);
      }
    }
    else {
      val[0]=m_mnp;
      val[1]=m_msum;
      val[2]=m_msum2;
      val[3]=m_max;
      MPI::COMM_WORLD.Send(&val,4,MPI::DOUBLE,0,rank);
      MPI::COMM_WORLD.Recv(&val,4,MPI::DOUBLE,0,size+rank);
      m_mnp=val[0];
      m_msum=val[1];
      m_msum2=val[2];
      m_max=val[3];
    }
  }
  m_np+=m_mnp;
  m_sum+=m_msum;
  m_sum2+=m_msum2;
  m_mnp=m_msum=m_msum2=0.0;
#endif
}

void PS_Vertex::AddPoint(const double &weight)
{
  double wgt(sqr(weight)*m_weight);
  if (IsBad(wgt)) return;
#ifdef USING__MPI
  ++m_mnp;
  m_msum+=wgt;
  m_msum2+=sqr(wgt);
#else
  ++m_np;
  m_sum+=wgt;
  m_sum2+=sqr(wgt);
#endif
  m_max=ATOOLS::Max(m_max,dabs(wgt));
}

void PS_Vertex::Reset()
{
  m_mnp=m_np=0.0;
  m_sum=m_sum2=m_max=0.0;
  m_msum=m_msum2=0.0;
}

DECLARE_GETTER(PS_Vertex_Getter,"PSV",Vertex_Base,Vertex_Key);

Vertex_Base *PS_Vertex_Getter::operator()(const Vertex_Key &key) const
{
  return new PS_Vertex(key);
}

void PS_Vertex_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"phase space vertex";
}

