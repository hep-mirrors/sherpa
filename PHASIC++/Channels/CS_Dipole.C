#include "PHASIC++/Channels/CS_Dipole.H"

#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Channels/Vegas.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace ATOOLS;
using namespace PHASIC;

CS_Dipole::CS_Dipole(NLO_subevt *const sub,
		     Phase_Space_Handler *const psh,const bool bmcw):
  m_sub(*sub), p_vegas(NULL),
  m_alpha(1.0), m_oldalpha(1.0), m_weight(1.0),
  m_amin(0.0), m_type(0), m_on(false), m_bmcw(bmcw)
{
  m_sub.p_ampl=NULL;
  p_fsmc=psh->FSRIntegrator();
  p_ismc=psh->ISRIntegrator();
  m_isrspkey.Assign("s' isr",5,0,psh->GetInfo());
  m_isrykey.Assign("y isr",3,0,psh->GetInfo());
  m_isrmode=psh->Process()->ISR()->On();
  for (size_t i(0);i<m_sub.m_n;++i) {
    size_t k(ID(sub->p_id[i]).front());
    if (k==m_sub.m_i) m_ijt=i;
    else if (k==m_sub.m_k) m_kt=i;
    m_rbmap[m_brmap[i]=k]=i;
  }
  m_rbmap[m_sub.m_i]=m_rbmap[m_sub.m_j]=m_ijt;
  m_flij=m_sub.p_fl[m_ijt];
  m_fli=m_sub.p_real->p_fl[m_sub.m_i];
  m_flj=m_sub.p_real->p_fl[m_sub.m_j];
  m_flk=m_sub.p_real->p_fl[m_sub.m_k];
  m_id=m_ijt<psh->Process()->NIn()?"I":"F";
  m_id+=m_kt<psh->Process()->NIn()?"I":"F";
  m_id+="_"+ToString(m_ijt)
    +"["+ToString(m_sub.m_i)+","+ToString(m_sub.m_j)+"]_"
    +ToString(m_kt)+"["+ToString(m_sub.m_k)+"]";
  if (m_ijt<psh->Process()->NIn()) m_type|=1;
  if (m_kt<psh->Process()->NIn()) m_type|=2;
  Reset();
}

CS_Dipole::~CS_Dipole() 
{
  if (p_vegas) delete p_vegas;
}

bool CS_Dipole::IsMapped(CS_Dipole *const dip) const
{
  if (m_sub.m_i!=dip->m_sub.m_i ||
      m_sub.m_j!=dip->m_sub.m_j ||
      m_sub.m_k!=dip->m_sub.m_k) return false;
  if (m_ijt!=dip->m_ijt || m_kt!=dip->m_kt) return false;
  if (m_brmap!=dip->m_brmap) return false;
  return true;
}

void CS_Dipole::InitVegas(const std::string &pid)
{
  p_vegas = new Vegas(3,100,m_id,0);
}

double CS_Dipole::Alpha(const int mode) const
{
  if (m_on) return mode?m_alpha:1.0;
  return 0.0;
}

void CS_Dipole::AddPoint(const double &value,
			 const double &ewgt,const int mode)
{
  if (value==0.0 || m_weight==0.0) return;
  double wgt(sqr(value)*ewgt/m_weight);
#ifdef USING__MPI
  ++m_mnp;
  m_msum+=wgt;
  m_msum2+=sqr(wgt);
#else
  ++m_np;
  m_sum+=wgt;
  m_sum2+=sqr(wgt);
#endif
  if (mode==1) p_vegas->AddPoint(dabs(wgt),m_rn);
}

void CS_Dipole::Optimize()
{
  p_vegas->Optimize();
}

void CS_Dipole::EndOptimize()
{
  p_vegas->EndOptimize();
}

void CS_Dipole::MPISync()
{
  p_vegas->MPISync();
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    double val[3], rval[3];
    val[0]=m_mnp;
    val[1]=m_msum;
    val[2]=m_msum2;
    mpi->MPIComm()->Allreduce(&val,&rval,3,MPI::DOUBLE,MPI::SUM);
    m_mnp=rval[0];
    m_msum=rval[1];
    m_msum2=rval[2];
  }
  m_np+=m_mnp;
  m_sum+=m_msum;
  m_sum2+=m_msum2;
  m_mnp=m_msum=m_msum2=0.0;
#endif
}

void CS_Dipole::Reset()
{
  m_np=m_mnp=0.0;
  m_sum=m_sum2=m_msum=m_msum2=0.0;
}

void CS_Dipole::WriteOut(const std::string &pid,
			 std::vector<std::string> &info) const
{
  p_vegas->WriteOut(pid);
  info.resize(6);
  info[0]=m_id;
  info[1]=ToString(m_alpha,12);
  info[2]=ToString(m_oldalpha,12);
  info[3]=ToString(m_np,12);
  info[4]=ToString(m_sum,12);
  info[5]=ToString(m_sum2,12);
}

void CS_Dipole::ReadIn(const std::string &pid,
		       const std::vector<std::string> &info)
{
  p_vegas->ReadIn(pid);
  if (info.size()!=6 || info[0]!=m_id)
    THROW(fatal_error,"Corrupted input file");
  m_alpha=ToType<double>(info[1],12);
  m_oldalpha=ToType<double>(info[2],12);
  m_np=ToType<double>(info[3],12);
  m_sum=ToType<double>(info[4],12);
  m_sum2=ToType<double>(info[5],12);
}

double CS_Dipole::Lambda
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

double CS_Dipole::Phi
(Vec4D pijt,Vec4D pkt,Vec4D pi,const bool ii) const
{
  Vec4D ktt(0.0,cross(Vec3D(pijt),Vec3D(pkt)));
  Poincare cms(pijt+pkt);
  cms.Boost(pijt);
  cms.Boost(pi);
  Poincare zax(pijt,Vec4D::ZVEC);
  if (!ii && ktt.PSpat2()>1.0e-6) zax.Rotate(ktt);
  else ktt=Vec4D(0.0,1.0,1.0,0.0);
  zax.Rotate(pi);
  Poincare xax(ktt,Vec4D::XVEC);
  xax.Rotate(pi);
  return pi.Phi();
}

double CS_Dipole::GetS
(const double &Q2,const double &y,
 const double &mi2,const double &mj2,const double &mk2) const
{
  return y*(Q2-mk2)+(1.0-y)*(mi2+mj2);
}

double CS_Dipole::GetZ
(const double &Q2,const double &sij,const double &y,const double &zt,
 const double &mi2,const double &mk2) const
{
  double ecm=0.5*(Q2-sij-mk2), rtlam=sqr(ecm);
  if (rtlam<sij*mk2) return sqrt(-1.0);
  rtlam=sqrt(rtlam-sij*mk2);
  double gam=ecm+Sign(Q2-sij-mk2)*rtlam;
  return ecm/rtlam*(zt-mk2/dabs(gam)*(y/(1.0-y)+mi2/ecm));
}

double CS_Dipole::GetKT2
(const double &Q2,const double &y,const double &z,
 const double &mi2,const double &mj2,const double &mk2) const
{
  return (Q2-mi2-mj2-mk2)*y*z*(1.0-z)-sqr(1.0-z)*mi2-z*z*mj2;
}

double CS_Dipole::ConstructLN
(const double &Q2,const double &sij,
 const double &mij2,const double &mk2,
 const Vec4D &Q,Vec4D &pk,Vec4D &l,Vec4D &n) const
{
  double po=sqr(Q2-mij2-mk2)-4.0*mij2*mk2, pn=sqr(Q2-sij-mk2)-4.0*sij*mk2;
  if ((po<0.0)^(pn<0.0)) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return 0.0;
  }
  pk=(Q2+mk2-sij)/(2.0*Q2)*Q+(pk-(Q2+mk2-mij2)/(2.0*Q2)*Q)*sqrt(pn/po);
  Vec4D pij=Q-pk;
  double gam=pij*pk+Sign(Q2-sij-mk2)*sqrt(sqr(pij*pk)-sij*mk2);
  double a13=sij/gam, a2=mk2/gam, bet=1.0/(1.0-a13*a2);
  l=bet*(pij-a13*pk);
  n=bet*(pk-a2*pij);
  return gam;
}

std::ostream &PHASIC::operator<<(std::ostream &ostr,const CS_Dipole &dip) 
{
  return ostr<<"("<<&dip<<")'"<<dip.Id()<<"': m_a = "<<dip.Alpha()
	     <<" <- "<<*dip.GetSubEvt();
}

