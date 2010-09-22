#include "PHASIC++/Channels/CS_Dipoles.H"

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Vegas.H"
#include "PHASIC++/Channels/Channel_Basics.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;
using namespace PHASIC;

FF_Dipole::FF_Dipole(NLO_subevt *const sub,
		     Phase_Space_Handler *const psh):
  CS_Dipole(sub,psh), m_yexp(0.5), m_zexp(0.01)
{
  // read in y,z mode
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa.GetPath());
  read.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  double helpd;
  if (read.ReadFromFile(helpd,"EEG_FF_Y_EXPONENT")) m_yexp=helpd;
  if (read.ReadFromFile(helpd,"EEG_FF_Z_EXPONENT")) m_zexp=helpd;
}

FF_Dipole::~FF_Dipole() {}

Vec4D_Vector FF_Dipole::GeneratePoint
(const Vec4D_Vector &p,Cut_Data *const cuts,const double *rns)
{
  DEBUG_FUNC("");
  // massless y- and z-bounds so far
  double *rn(p_vegas->GeneratePoint(rns));
  msg_Debugging()<<"vegased :     ";
  msg_Debugging()<<"y = "<<rn[0]<<", z = "<<rn[1]
                 <<", phi = "<<rn[2]<<"\n";
  m_rn[0]=Channel_Basics::PeakedDist(0.0,m_yexp,m_amin,1.0,1,rn[0]);
  m_rn[1]=Channel_Basics::PeakedDist(0.0,m_zexp,0.0,1.0,1,rn[1]);
  m_rn[2]=rn[2]*2.0*M_PI;
  msg_Debugging()<<"transformed : ";
  msg_Debugging()<<"y = "<<m_rn[0]<<", z = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  Vec4D_Vector pp(p.size()+1);
  for (size_t i(0);i<p.size();++i) pp[m_brmap[i]]=p[i];
  Construct(pp[m_sub.m_i],pp[m_sub.m_j],pp[m_sub.m_k],
	    m_rn[0],m_rn[1],m_rn[2],p[m_ijt],p[m_kt]);
  return pp;
}

double FF_Dipole::GenerateWeight(const Vec4D_Vector &p,Cut_Data *const cuts)
{
  // massless y-/z-bounds and weight so far
  Vec4D_Vector pp(p.size()-1);
  for (size_t i(0);i<p.size();++i) pp[m_rbmap[i]]=p[i];
  Calculate(p[m_sub.m_i],p[m_sub.m_j],p[m_sub.m_k],
	    m_rn[0],m_rn[1],m_rn[2],pp[m_ijt],pp[m_kt]);
  p_fsmc->GenerateWeight(&pp.front(),cuts);
  if (p_ismc) {
    m_isrspkey[3]=(pp[0]+pp[1]).Abs2();
    m_isrykey[2]=(pp[0]+pp[1]).Y();
    p_ismc->GenerateWeight(m_isrmode);
  }
  if (m_rn[2]<0.0) m_rn[2]+=2.0*M_PI;
  msg_Debugging()<<"again :       ";
  msg_Debugging()<<"y = "<<m_rn[0]<<", z = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  if (m_rn[0]<m_amin) {
    m_weight=0.0;
    return 0.0;
  }
  m_weight=(pp[m_ijt]+pp[m_kt]).Abs2()/(16.0*sqr(M_PI))*(1.0-m_rn[0]);
  m_weight*=pow(m_rn[0],m_yexp)*pow(m_rn[1],m_zexp);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_yexp,m_amin,1.0,m_rn[0],1,m_rn[0]);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_zexp,0.0,1.0,m_rn[1],1,m_rn[1]);
  m_rn[2]/=2.0*M_PI;
  msg_Debugging()<<"recovered :   ";
  msg_Debugging()<<"y = "<<m_rn[0]<<", z = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  m_weight*=p_vegas->GenerateWeight(m_rn);
  if (p_ismc) m_weight*=p_ismc->Weight();
  return m_weight*=p_fsmc->Weight();
}

void FF_Dipole::Calculate
(const Vec4D &pi,const Vec4D &pj,const Vec4D &pk,
 double &y,double &z,double &phi,Vec4D &pijt,Vec4D &pkt)
{
  Vec4D  Q(pi+pj+pk);
  double pipj=pi*pj, sij=(pi+pj).Abs2(), pipk=pi*pk, pjpk=pj*pk;
  double Q2=Q*Q, mk2=sqr(m_flk.Mass()), mij2=sqr(m_flij.Mass());
  y=pipj/(pipj+pipk+pjpk);
  z=pipk/(pipk+pjpk);
  pkt=Vec4D(sqrt(Lambda(Q2,mij2,mk2)/Lambda(Q2,sij,mk2))*
	    (pk-(Q*pk)/Q2*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q);
  phi=Phi(Q-pkt,pkt,pi,0);
  pijt=Q-pkt;
}

void FF_Dipole::Construct
(ATOOLS::Vec4D &pi,ATOOLS::Vec4D &pj,ATOOLS::Vec4D &pk,
 const double &y,const double &z,const double &phi,
 const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &pkt)
{
  Vec4D rpijt(pijt), Q(pijt+pkt), l, n;
  Vec4D n_perp(0.0,cross(Vec3D(pijt),Vec3D(pkt)));
  Poincare cms(Q);
  cms.Boost(rpijt);
  Poincare zrot(rpijt,Vec4D::ZVEC);
  if (n_perp.PSpat2()<=1.0e-6) {
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
    zrot.RotateBack(n_perp);
  }
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(0.0,cross(Vec3D(rpijt),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();
  double mi2(sqr(m_fli.Mass())), mj2(sqr(m_flj.Mass()));
  double mk2(sqr(m_flk.Mass())), mij2(sqr(m_flij.Mass()));
  double Q2(Q.Abs2()), sij(GetS(Q2,y,mi2,mj2,mk2));
  double gam(ConstructLN(Q2,sij,mij2,mk2,Q,pk=pkt,l,n));
  double zt(GetZ(Q2,sij,y,z,mi2,mk2));
  double ktt(GetKT2(Q2,y,zt,mi2,mj2,mk2));
  if (ktt<0.0) THROW(fatal_error,"Invalid kinematics");
  ktt=sqrt(ktt);
  pi=ktt*sin(phi)*l_perp;
  cms.BoostBack(pi);
  pi+=zt*l+(mi2+ktt*ktt)/(gam*zt)*n+ktt*cos(phi)*n_perp;
  pj=Q-pi-pk;
}

FI_Dipole::FI_Dipole(ATOOLS::NLO_subevt *const sub,
		     Phase_Space_Handler *const psh):
  CS_Dipole(sub,psh), m_xexp(0.5), m_zexp(0.01)
{
  // read in x,z mode
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa.GetPath());
  read.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  double helpd;
  if (read.ReadFromFile(helpd,"EEG_FI_X_EXPONENT")) m_xexp=helpd;
  if (read.ReadFromFile(helpd,"EEG_FI_Z_EXPONENT")) m_zexp=helpd;
}

FI_Dipole::~FI_Dipole() {}

ATOOLS::Vec4D_Vector FI_Dipole::GeneratePoint
(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts,const double *rns)
{
  DEBUG_FUNC("");
  // massless x- and z-bounds so far
  double *rn(p_vegas->GeneratePoint(rns));
  msg_Debugging()<<"vegased :     ";
  if (m_kt==0) m_xmin=p[m_kt].PPlus()/rpa.gen.PBeam(0).PPlus();
  else m_xmin=p[m_kt].PMinus()/rpa.gen.PBeam(1).PMinus();
  msg_Debugging()<<"x = "<<rn[0]<<", z = "<<rn[1]
                 <<", phi = "<<rn[2]<<", xmin = "<<m_xmin<<"\n";
  m_rn[0]=Channel_Basics::PeakedDist(0.0,m_xexp,m_xmin,1.0,1,rn[0]);
  m_rn[1]=Channel_Basics::PeakedDist(0.0,m_zexp,0.0,1.0,1,rn[1]);
  m_rn[2]=rn[2]*2.0*M_PI;
  msg_Debugging()<<"transformed : ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", z = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  Vec4D_Vector pp(p.size()+1);
  for (size_t i(0);i<p.size();++i) pp[m_brmap[i]]=p[i];
  Construct(pp[m_sub.m_i],pp[m_sub.m_j],pp[m_sub.m_k],
            m_rn[0],m_rn[1],m_rn[2],p[m_ijt],p[m_kt]);
  return pp;
}

double FI_Dipole::GenerateWeight
(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts)
{
  // massless x-/z-bounds and weight so far
  Vec4D_Vector pp(p.size()-1);
  for (size_t i(0);i<p.size();++i) pp[m_rbmap[i]]=p[i];
  Calculate(p[m_sub.m_i],p[m_sub.m_j],p[m_sub.m_k],
	    m_rn[0],m_rn[1],m_rn[2],pp[m_ijt],pp[m_kt]);
  if (m_rn[2]<0.0) m_rn[2]+=2.0*M_PI;
  if (m_kt==0) m_xmin=pp[m_kt].PPlus()/rpa.gen.PBeam(0).PPlus();
  else m_xmin=pp[m_kt].PMinus()/rpa.gen.PBeam(1).PMinus();
  msg_Debugging()<<"again :       ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", z = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  if (m_rn[0]<m_xmin) {
    m_weight=0.0;
    return 0.0;
  }
  p_fsmc->GenerateWeight(&pp.front(),cuts);
  m_isrspkey[3]=(pp[0]+pp[1]).Abs2();
  m_isrykey[2]=(pp[0]+pp[1]).Y();
  p_ismc->GenerateWeight(m_isrmode);
  // 2(pijt*pk)/16pi^2
  m_weight=(pp[m_ijt]+p[m_sub.m_k]).Abs2()/
    (16.0*sqr(M_PI))/m_rn[0];
  m_weight*=pow(m_rn[0],m_xexp)*pow(m_rn[1],m_zexp);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_xexp,m_xmin,1.0,m_rn[0],1,m_rn[0]);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_zexp,0.0,1.0,m_rn[1],1,m_rn[1]);
  m_rn[2]/=2.0*M_PI;
  msg_Debugging()<<"recovered :   ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", z = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<", xmin = "<<m_xmin<<"\n";
  m_weight*=p_vegas->GenerateWeight(m_rn);
  return m_weight*=p_fsmc->Weight()*p_ismc->Weight();
}

void FI_Dipole::Calculate
(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj, const ATOOLS::Vec4D &pk,
 double &x, double &z, double &phi,
 ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt)
{
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  x=(pjpk+pipk-pipj)/(pipk+pjpk);
  z=pjpk/(pipk+pjpk);

  pijt=pi+pj-(1.-x)*pk;
  pkt=x*pk;

  double kp(sqrt(2.*pijt*pk*(1.-x)*z*(1.-z)));
  Vec4D kperp(pj-z*pijt-(1.-x)*(1.-z)*pk);
  Poincare bf(pijt+pkt);
  Poincare rot(bf*pijt,Vec4D(0.,0.,0.,1.));
  bf.Boost(kperp);
  rot.Rotate(kperp);
  if      ((kperp[1]>=0. && kperp[2]>=0.) || (kperp[1]>=0. && kperp[2]<0.))
    phi=acos(kperp[2]/kp);
  else if ((kperp[1]<0. && kperp[2]<0.) || (kperp[1]<0. && kperp[2]>=0.))
    phi=-acos(kperp[2]/kp)+2.*M_PI;
  else THROW(fatal_error,"Could not determine phi.");
}

void FI_Dipole::Construct
(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
 const double &x, const double &z, const double &phi,
 const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt)
{
  DEBUG_FUNC("");
  pk = 1./x*pkt;

  double kp(sqrt(2.*pijt*pk*(1.-x)*z*(1.-z)));
  Poincare bf(pkt+pijt);
  Poincare rot(bf*pijt,Vec4D(0.,0.,0.,1.));
  Vec4D kperp(0.,sin(phi)*kp,cos(phi)*kp,0.);
  rot.RotateBack(kperp);
  bf.BoostBack(kperp);

  pi = (1.-z)*pijt + (1.-x)*z*pk - kperp;
  pj = z*pijt + (1.-x)*(1.-z)*pk + kperp;

  msg_Debugging()<<"("<<m_ijt<<"):"<<pijt
                 <<" -> ("<<m_sub.m_i<<"):"<<pi
                 <<" ("<<m_sub.m_j<<"):"<<pj<<std::endl;
  msg_Debugging()<<"("<<m_kt<<"):"<<pkt
                 <<" -> ("<<m_sub.m_k<<"):"<<pk<<std::endl;
}

IF_Dipole::IF_Dipole(ATOOLS::NLO_subevt *const sub,
		     Phase_Space_Handler *const psh):
  CS_Dipole(sub,psh), m_xexp(0.5), m_uexp(0.5)
{
  // read in x,u mode
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa.GetPath());
  read.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  double helpd;
  if (read.ReadFromFile(helpd,"EEG_IF_X_EXPONENT")) m_xexp=helpd;
  if (read.ReadFromFile(helpd,"EEG_IF_U_EXPONENT")) m_uexp=helpd;
}

IF_Dipole::~IF_Dipole() {}

ATOOLS::Vec4D_Vector IF_Dipole::GeneratePoint
(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts,const double *rns)
{
  DEBUG_FUNC("");
  // massless x- and u-bounds so far
  double *rn(p_vegas->GeneratePoint(rns));
  if (m_ijt==0) m_xmin=p[m_ijt].PPlus()/rpa.gen.PBeam(0).PPlus();
  else m_xmin=p[m_ijt].PMinus()/rpa.gen.PBeam(1).PMinus();
  msg_Debugging()<<"vegased :     ";
  msg_Debugging()<<"x = "<<rn[0]<<", u = "<<rn[1]
                 <<", phi = "<<rn[2]<<", xmin = "<<m_xmin<<"\n";
  m_rn[0]=Channel_Basics::PeakedDist(0.0,m_xexp,m_xmin,1.0,1,rn[0]);
  m_rn[1]=Channel_Basics::PeakedDist(0.0,m_uexp,m_amin,1.0,1,rn[1]);
  m_rn[2]=rn[2]*2.0*M_PI;
  msg_Debugging()<<"transformed : ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", u = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  Vec4D_Vector pp(p.size()+1);
  for (size_t i(0);i<p.size();++i) pp[m_brmap[i]]=p[i];
  Construct(pp[m_sub.m_i],pp[m_sub.m_j],pp[m_sub.m_k],
            m_rn[0],m_rn[1],m_rn[2],p[m_ijt],p[m_kt]);
  if (m_ijt!=m_sub.m_i) {
    for (size_t i(0);i<pp.size();++i) pp[i]=Rotate(pp[i]);
  }
  return pp;
}

double IF_Dipole::GenerateWeight
(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts)
{
  // massless x-/v-bounds and weight so far
  Vec4D_Vector pp(p.size()-1);
  for (size_t i(0);i<p.size();++i) pp[m_rbmap[i]]=p[i];
  if (m_ijt==m_sub.m_i) {
  Calculate(p[m_sub.m_i],p[m_sub.m_j],p[m_sub.m_k],
	    m_rn[0],m_rn[1],m_rn[2],pp[m_ijt],pp[m_kt]);
  }
  else {
    for (size_t i(0);i<pp.size();++i) pp[i]=Rotate(pp[i]);
    Calculate(Rotate(p[m_sub.m_i]),Rotate(p[m_sub.m_j]),
	      Rotate(p[m_sub.m_k]),m_rn[0],m_rn[1],m_rn[2],
	      pp[m_ijt],pp[m_kt]);
  }
  if (m_rn[2]<0.0) m_rn[2]+=2.0*M_PI;
  if (m_ijt==0) m_xmin=pp[m_ijt].PPlus()/rpa.gen.PBeam(0).PPlus();
  else m_xmin=pp[m_ijt].PMinus()/rpa.gen.PBeam(1).PMinus();
  msg_Debugging()<<"again :       ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", u = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  if (m_rn[0]<m_xmin) {
    m_weight=0.0;
    return 0.0;
  }
  p_fsmc->GenerateWeight(&pp.front(),cuts);
  m_isrspkey[3]=(pp[0]+pp[1]).Abs2();
  m_isrykey[2]=(pp[0]+pp[1]).Y();
  p_ismc->GenerateWeight(m_isrmode);
  // 2(pa*pkt)/16pi^2
  m_weight=(pp[m_ijt]+pp[m_kt]).Abs2()/
    (16.0*sqr(M_PI))/sqr(m_rn[0]);
  m_weight*=pow(m_rn[1],m_uexp)*pow(m_rn[0],m_xexp);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_uexp,m_amin,1.0,m_rn[1],1,m_rn[1]);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_xexp,m_xmin,1.0,m_rn[0],1,m_rn[0]);
  m_rn[2]/=2.0*M_PI;
  msg_Debugging()<<"recovered :   ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", u = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<", xmin = "<<m_xmin<<"\n";
  m_weight*=p_vegas->GenerateWeight(m_rn);
  return m_weight*=p_fsmc->Weight()*p_ismc->Weight();
}

void IF_Dipole::Calculate
(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj, const ATOOLS::Vec4D &pk,
 double &x, double &u, double &phi,
 ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt)
{
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  x=(pipj+pipk-pjpk)/(pipj+pipk);
  u=pipj/(pipj+pipk);

  pijt=x*pi;
  pkt=pk+pj-(1.-x)*pi;

  double kp(sqrt(2.*pi*pkt*(1.-x)*u*(1.-u)));
  Vec4D kperp(pj-(1.-u)*(1.-x)/x*pijt-u*pkt);
  Poincare bf(pijt+pkt);
  Poincare rot(bf*pijt,Vec4D(0.,0.,0.,1.));
  bf.Boost(kperp);
  rot.Rotate(kperp);
  if      ((kperp[1]>=0. && kperp[2]>=0.) || (kperp[1]>=0. && kperp[2]<0.))
    phi=acos(kperp[2]/kp);
  else if ((kperp[1]<0. && kperp[2]<0.) || (kperp[1]<0. && kperp[2]>=0.))
    phi=-acos(kperp[2]/kp)+2.*M_PI;
  else THROW(fatal_error,"Could not determine phi.");
}

void IF_Dipole::Construct
(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
 const double &x, const double &u, const double &phi,
 const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt)
{
  DEBUG_FUNC("");
  pi = 1./x*pijt;

  double kp(sqrt(2.*pi*pkt*(1.-x)*u*(1.-u)));
  Poincare bf(pkt+pijt);
  Poincare rot(bf*pijt,Vec4D(0.,0.,0.,1.));
  Vec4D kperp(0.,sin(phi)*kp,cos(phi)*kp,0.);
  rot.RotateBack(kperp);
  bf.BoostBack(kperp);

  pj = (1.-u)*(1.-x)/x*pijt + u*pkt + kperp;
  pk = u*(1.-x)/x*pijt + (1.-u)*pkt - kperp;

  msg_Debugging()<<"("<<m_ijt<<"):"<<pijt
                 <<" -> ("<<m_sub.m_i<<"):"<<pi
                 <<" ("<<m_sub.m_j<<"):"<<pj<<std::endl;
  msg_Debugging()<<"("<<m_kt<<"):"<<pkt
                 <<" -> ("<<m_sub.m_k<<"):"<<pk<<std::endl;
}

II_Dipole::II_Dipole(ATOOLS::NLO_subevt *const sub,
		     Phase_Space_Handler *const psh):
  CS_Dipole(sub,psh), m_xexp(0.5), m_vexp(0.5)
{
  // read in x,v mode
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa.GetPath());
  read.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  double helpd;
  if (read.ReadFromFile(helpd,"EEG_II_X_EXPONENT")) m_xexp=helpd;
  if (read.ReadFromFile(helpd,"EEG_II_V_EXPONENT")) m_vexp=helpd;
}

II_Dipole::~II_Dipole() {}

ATOOLS::Vec4D_Vector II_Dipole::GeneratePoint
(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts,const double *rns)
{
  DEBUG_FUNC("");
  // massless x- and v-bounds so far
  double *rn(p_vegas->GeneratePoint(rns));
  if (m_ijt==0) m_xmin=p[m_ijt].PPlus()/rpa.gen.PBeam(0).PPlus();
  else m_xmin=p[m_ijt].PMinus()/rpa.gen.PBeam(1).PMinus();
  msg_Debugging()<<"vegased :     ";
  msg_Debugging()<<"x = "<<rn[0]<<", v = "<<rn[1]
                 <<", phi = "<<rn[2]<<", xmin = "<<m_xmin<<"\n";
  m_rn[0]=Channel_Basics::PeakedDist(0.0,m_xexp,m_xmin,1.0,1,rn[0]);
  double amin(Min(m_amin,0.5*(1.0-m_rn[0])));
  m_rn[1]=Channel_Basics::PeakedDist(0.0,m_vexp,amin,1.0-m_rn[0],1,rn[1]);
  m_rn[2]=rn[2]*2.0*M_PI;
  msg_Debugging()<<"transformed : ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", v = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  Vec4D_Vector pp(p.size()+1);
  for (size_t i(0);i<p.size();++i) pp[m_brmap[i]]=p[i];
  if (m_rn[1]>1.-m_rn[0]) {
    msg_Error()<<METHOD<<"(): v > 1-x, "<<m_rn[1]
	       <<" vs. "<<1.0-m_rn[0]<<"\n";
    m_rn[1]=(1.0-1.0e-6)*(1.0-m_rn[0]);
  }
  Construct(pp[m_sub.m_i],pp[m_sub.m_j],pp[m_sub.m_k],pp,
            m_rn[0],m_rn[1],m_rn[2],p[m_ijt],p[m_kt],p);
  return pp;
}

double II_Dipole::GenerateWeight
(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts)
{
  // massless x-/v-bounds and weight so far
  Vec4D_Vector pp(p.size()-1);
  Calculate(p[m_sub.m_i],p[m_sub.m_j],p[m_sub.m_k],p,
            m_rn[0],m_rn[1],m_rn[2],pp[m_ijt],pp[m_kt],pp);
  if (m_rn[2]<0.0) m_rn[2]+=2.0*M_PI;
  if (m_ijt==0) m_xmin=pp[m_ijt].PPlus()/rpa.gen.PBeam(0).PPlus();
  else m_xmin=pp[m_ijt].PMinus()/rpa.gen.PBeam(1).PMinus();
  msg_Debugging()<<"again :       ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", v = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<"\n";
  if (m_rn[0]<m_xmin) {
    m_weight=0.0;
    return 0.0;
  }
  p_fsmc->GenerateWeight(&pp.front(),cuts);
  m_isrspkey[3]=(pp[0]+pp[1]).Abs2();
  m_isrykey[2]=(pp[0]+pp[1]).Y();
  p_ismc->GenerateWeight(m_isrmode);
  // 2(pa*pb)/16pi^2
  m_weight=(p[m_sub.m_i]+p[m_sub.m_k]).Abs2()/
    (16.0*sqr(M_PI))/m_rn[0];
  double amin(Min(m_amin,0.5*(1.0-m_rn[0])));
  m_weight*=pow(m_rn[1],m_vexp)*pow(m_rn[0],m_xexp);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_vexp,amin,1.0-m_rn[0],m_rn[1],1,m_rn[1]);
  m_weight*=Channel_Basics::PeakedWeight
    (0.0,m_xexp,m_xmin,1.0,m_rn[0],1,m_rn[0]);
  m_rn[2]/=2.0*M_PI;
  msg_Debugging()<<"recovered :   ";
  msg_Debugging()<<"x = "<<m_rn[0]<<", v = "<<m_rn[1]
                 <<", phi = "<<m_rn[2]<<", xmin = "<<m_xmin<<"\n";
  m_weight*=p_vegas->GenerateWeight(m_rn);
  return m_weight*=p_fsmc->Weight()*p_ismc->Weight();
}

void II_Dipole::Calculate
(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj, const ATOOLS::Vec4D &pk,
 const ATOOLS::Vec4D_Vector& kj,
 double &x, double &v, double &phi,
 ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt,
 ATOOLS::Vec4D_Vector& kjt)
{
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  x=(pipk-pipj-pjpk)/pipk;
  v=pipj/pipk;

  pijt=x*Rotate(pi);
  pkt=Rotate(pk);

  double kp(sqrt(2.*pipk*v*(1.-x-v)));
  Vec4D kperp(Rotate(pj)-(1.-x-v)/x*pijt-v*pkt);
  if      ((kperp[1]>=0. && kperp[2]>=0.) || (kperp[1]>=0. && kperp[2]<0.))
    phi=acos(kperp[2]/kp);
  else if ((kperp[1]<0. && kperp[2]<0.) || (kperp[1]<0. && kperp[2]>=0.))
    phi=-acos(kperp[2]/kp)+2.*M_PI;
  else THROW(fatal_error,"Could not determine phi.");

  Vec4D K(pi-pj+pk), Kt(pijt+pkt);
  ATOOLS::Lorentz_Ten2D Lambda = MetricTensor()
                                 - 2./(K+Kt).Abs2()*BuildTensor(Kt+K,Kt+K)
                                 + 2./Kt.Abs2()*BuildTensor(Kt,K);

  for (size_t j(2), i(j);j<kjt.size();++i,++j) {
    if (i==m_sub.m_j) ++i;
    kjt[m_rbmap[i]] = Contraction(Lambda,2,kj[i]);
    msg_Debugging()<<"("<<i<<"):"<<kj[i]
		   <<" -> ("<<m_rbmap[i]<<"):"<<kjt[m_rbmap[i]]<<std::endl;
  }
}

void II_Dipole::Construct
(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
 ATOOLS::Vec4D_Vector& kj,
 const double &x, const double &v, const double &phi,
 const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt,
 const ATOOLS::Vec4D_Vector& kjt)
{
  DEBUG_FUNC("");
  pi=1./x*pijt;
  pk=pkt;

  double kp=sqrt(2.*pi*pk*v*(1.-x-v));
  Vec4D kperp(0.,sin(phi)*kp,cos(phi)*kp,0.);

  pj=(1.-x-v)/x*pijt + v*pkt + kperp;
  pi=Rotate(pi);
  pj=Rotate(pj);
  pk=Rotate(pk);
  msg_Debugging()<<"("<<m_ijt<<"):"<<pijt
                 <<" -> ("<<m_sub.m_i<<"):"<<pi
                 <<" ("<<m_sub.m_j<<"):"<<pj<<std::endl;
  msg_Debugging()<<"("<<m_kt<<"):"<<pkt
                 <<" -> ("<<m_sub.m_k<<"):"<<pk<<std::endl;

  Vec4D K(pi-pj+pk), Kt(pijt+pkt);
  ATOOLS::Lorentz_Ten2D Lambda = MetricTensor()
                                 - 2./(K+Kt).Abs2()*BuildTensor(Kt+K,Kt+K)
                                 + 2./Kt.Abs2()*BuildTensor(K,Kt);

  for (size_t i(0);i<kjt.size();++i) {
    if (i!=m_ijt && i!=m_kt) {
      kj[m_brmap[i]] = Contraction(Lambda,2,kjt[i]);
      msg_Debugging()<<"("<<i<<"):"<<kjt[i]
                     <<" -> ("<<m_brmap[i]<<"):"<<kj[m_brmap[i]]<<std::endl;
    }
  }
}

