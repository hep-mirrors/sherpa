#include "PHASIC++/Channels/CS_Dipoles.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace ATOOLS;
using namespace PHASIC;

II_Dipole::II_Dipole(const ATOOLS::Cluster_Leg &lij,
                     const ATOOLS::Cluster_Leg &lk,
                     const ATOOLS::Cluster_Leg &li,
                     const ATOOLS::Cluster_Leg &lj) :
  CS_Dipole(lij,lk,li,lj), m_xmode(1), m_vmode(1), m_exp(2)
{
  // read in y,z mode
  Data_Reader read(" ",";","!","=");
  read.SetInputFile("Integration.dat");
  std::string xmode,vmode;
  if (!read.ReadFromFile(xmode,"EEG_II_X_DICING")) xmode="logarithmic";
//   else msg_Info()<<METHOD<<"(): Set x mode in II dipole '"<<m_id
//                  <<"' to "<<xmode<<".\n";
  if (!read.ReadFromFile(vmode,"EEG_II_V_DICING")) vmode="linear";
//   else msg_Info()<<METHOD<<"(): Set v mode in II dipole '"<<m_id
//                  <<"' to "<<vmode<<".\n";
  if      (xmode=="linear")      m_xmode=0;
  else if (xmode=="logarithmic") m_xmode=1;
  else if (xmode=="quadratic")   m_xmode=2;
  else THROW(fatal_error,"unknown EEG_II_X_DICING mode");
  if      (vmode=="linear")      m_vmode=0;
  else if (vmode=="power")       m_vmode=1;
  else if (vmode=="exponential") m_vmode=2;
  else THROW(fatal_error,"unknown EEG_II_V_DICING mode");

  if (m_vmode==2) THROW(not_implemented,"not implemented yet");
  THROW(not_implemented,"II dipole splitting for the EEG not implemented yet");
}

II_Dipole::~II_Dipole() {}

CS_Dipole *II_Dipole::Copy() const
{
  return new II_Dipole(*this);
}

bool II_Dipole::GeneratePoint
(Cluster_Amplitude *const ampl,const double *rns)
{
//   DEBUG_FUNC("");
//   // massless y- and z-bounds so far
//   double *rn(p_vegas->GeneratePoint(rns));
//   msg_Debugging()<<"vegased :     ";
//   msg_Debugging()<<"x = "<<rn[0]<<", z = "<<rn[1]
//                  <<", phi = "<<rn[2]<<"\n";
// 
//   // x
//   // linear f(x)~1
//   if      (m_xmode==0) m_rn[0]=m_amin+rn[0]*(m_amax-m_amin);
//   // logarithmic f(x)~1/x
//   else if (m_xmode==1) m_rn[0]=m_amin*pow(m_amax/m_amin,rn[0]);
//   // quadratic f(x)~(x-a/2)
//   else if (m_xmode==2) m_rn[0]=m_amin/2.0+sqrt(sqr(m_amin/2.0)+rn[0]*m_amax*(m_amax-m_amin));
//   else THROW(fatal_error,"m_xmode changed to undefined value");
// 
//   // v
//   // linear f(x)~1
//   if      (m_vmode==0) m_rn[1]=rn[1];
//   // power f(x)~(x-1/2)^(2n)
//   else if (m_vmode==1)
//     if (rn[1]>0.5)
//       m_rn[1]=(1.+pow(2.*rn[1]-1.,1./(2.*m_exp+1.)))/2.;
//     else
//       m_rn[1]=(1.-pow(-2.*rn[1]+1.,1./(2.*m_exp+1.)))/2.;
//   // exponential
//   else if (m_vmode==2) m_rn[1]=rn[1];
//   else THROW(fatal_error,"m_vmode changed to undefined value");
// 
// 
//   // phi
//   // linear f(x)~1
//   m_rn[2]=rn[2]*2.0*M_PI;
// 
//   msg_Debugging()<<"transformed : ";
//   msg_Debugging()<<"x = "<<m_rn[0]<<", v = "<<m_rn[1]
//                  <<", phi = "<<m_rn[2]<<"\n";
//   Cluster_Leg *lij(ampl->IdLeg(m_lij.Id()));
//   Cluster_Leg *lk(ampl->IdLeg(m_lk.Id()));
//   Map(ampl);
//   Construct(ampl->IdLeg(m_li.Id()),ampl->IdLeg(m_lj.Id()),
//             lk,lij,m_rn[0],m_rn[1],m_rn[2]);
//   return true;
  return true;
}

double II_Dipole::GenerateWeight(const Cluster_Amplitude *ampl)
{
//   // massless y-/z-bounds and weight so far
//   Cluster_Amplitude *wampl(ampl->Copy());
//   ReMap(wampl);
//   Vec4D pijt, pkt;
//   Calculate(ampl->IdLeg(m_li.Id()),
//             ampl->IdLeg(m_lj.Id()),wampl->IdLeg(m_lk.Id()),
//             m_rn[0],m_rn[1],m_rn[2],pijt,pkt);
//   delete wampl;
//   if (m_rn[2]<0.0) m_rn[2]+=2.0*M_PI;
//   msg_Debugging()<<"again :       ";
//   msg_Debugging()<<"x = "<<m_rn[0]<<", v = "<<m_rn[1]
//                  <<", phi = "<<m_rn[2]<<"\n";
//   if (m_rn[0]<m_amin) {
//     m_weight=0.0;
//     return 0.0;
//   }
//   // does this read any different with alpha<1 ?
//   m_weight=(pijt+pkt).Abs2()/(16.0*sqr(M_PI))*(1.0-m_rn[0]);
// 
//   // x
//   // linear
//   if (m_xmode==0) {
//     m_weight*=(m_amax-m_amin);
//     m_rn[0]=(m_rn[0]-m_amin)/(m_amax-m_amin);
//   }
//   // logarithmic
//   else if (m_xmode==1) {
//     m_weight*=m_rn[0]*log(m_amax/m_amin);
//     m_rn[0]=log(m_rn[0]/m_amin)/log(m_amax/m_amin);
//   }
//   // quadratic
//   else if (m_xmode==2) {
//     m_weight*=m_amax*(m_amax-m_amin)/(2.*(m_rn[0]-m_amin/2.0));
//     m_rn[0]=m_rn[0]*(m_rn[0]-m_amin)/(m_amax*(m_amax-m_amin));
//   }
//   else THROW(fatal_error,"m_xmode changed to undefined value");
// 
//   // v
//   // power
//   if (m_vmode==1) {
//     double rn1old=m_rn[1];
//     m_rn[1]=(pow(2.*m_rn[1]-1.,2.*m_exp+1.)+1.)/2.;
//     m_weight*=1./(2.*(2.*m_exp+1.))*(rn1old-.5)/(m_rn[1]-.5);
//   }
//   // exponential
//   else if (m_vmode==2) {
//     THROW(fatal_error, "");
//   }
// 
//   m_rn[2]/=2.0*M_PI;
//   msg_Debugging()<<"recovered :   ";
//   msg_Debugging()<<"x = "<<m_rn[0]<<", v = "<<m_rn[1]
//                  <<", phi = "<<m_rn[2]<<"\n";
//   return m_weight*=p_vegas->GenerateWeight(m_rn);
  return 0.;
}

void II_Dipole::Calculate
(const Cluster_Leg *li,const Cluster_Leg *lj,
 const Cluster_Leg *lk,double &x,double &v,double &phi,
 Vec4D &pijt,Vec4D &pkt)
{
//   Vec4D  Q(li->Mom()+lj->Mom()+lk->Mom());
//   double pipj=li->Mom()*lj->Mom(), pipj2=(li->Mom()+lj->Mom()).Abs2();
//   double pipk=li->Mom()*lk->Mom(), pjpk=lj->Mom()*lk->Mom(), Q2=Q*Q;
//   double mk2=sqr(lk->Flav().Mass()), mij2=sqr(m_lij.Flav().Mass());
//   y=pipj/(pipj+pipk+pjpk);
//   z=pipk/(pipk+pjpk);
//   pkt=Vec4D(sqrt(Lambda(Q2,mij2,mk2)/Lambda(Q2,pipj2,mk2))*
//             (lk->Mom()-(Q*lk->Mom())/Q2*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q);
//   phi=Phi(Q-pkt,pkt,li->Mom(),0);
//   pijt=Q-pkt;
}

void II_Dipole::Construct
(Cluster_Leg *const li,Cluster_Leg *const lj,Cluster_Leg *const lk,
 Cluster_Leg *const lij,const double &x,const double &v,const double &phi)
{
//   Vec4D p1(lij->Mom()), p2(lk->Mom());
//   Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));
//   Poincare cms(p1+p2);
//   cms.Boost(p1);
//   cms.Boost(p2);
//   Poincare zrot(p1,Vec4D::ZVEC);
//   zrot.Rotate(p1);
//   zrot.Rotate(p2);
//   if (n_perp.PSpat2()>1.0e-6) zrot.Rotate(n_perp);
//   else n_perp=Vec4D(0.0,1.0,1.0,0.0);
//   Poincare xrot(n_perp,Vec4D::XVEC);
//   Vec4D m_p[3], Q(p1+p2);
//   double mi2=sqr(m_li.Flav().Mass()), mj2=sqr(m_lj.Flav().Mass());
//   double mk2=sqr(m_lk.Flav().Mass()), mij2=sqr(m_lij.Flav().Mass()), Q2=Q.Abs2();
//   double sij=y*(Q2-mk2)+(1.0-y)*(mi2+mj2);
//   double po=sqr(Q2-mij2-mk2)-4.0*mij2*mk2, pn=sqr(Q2-sij-mk2)-4.0*sij*mk2;
//   if (po<0.0 || pn<0.0) THROW(fatal_error,"Invalid kinematics");
//   m_p[2]=(Q2+mk2-sij)/(2.0*Q2)*Q+(p2-(Q2+mk2-mij2)/(2.0*Q2)*Q)*sqrt(pn/po);
//   Vec4D q13(Q-m_p[2]);
//   double gam=q13*m_p[2]+sqrt(sqr(q13*m_p[2])-sij*mk2);
//   double a13=sij/gam, a2=mk2/gam, bet=1.0/(1.0-a13*a2);
//   Vec4D l=bet*(q13-a13*m_p[2]), n=bet*(m_p[2]-a2*q13);
//   double zt=(z*(1.0+a13*a2)-a2*(2.0*mi2/gam+y/(1.0-y)*(1.0+a13*a2)))/
//     (1.0-a2*((mi2+mj2)/gam+y/(1.0-y)*(1.0+a13*a2)));
//   double ktt=gam*y/(1.0-y)*(1.0+a13*a2)*zt*(1.0-zt)-sqr(1.0-zt)*mi2-zt*zt*mj2;
//   if (ktt<0.0) THROW(fatal_error,"Invalid kinematics");
//   ktt=sqrt(ktt);
//   m_p[0]=zt*l+(mi2+ktt*ktt)/(gam*zt)*n
//     +ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0)
//     +ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0);
//   m_p[1]=Q-m_p[2]-m_p[0];
//   if (m_p[0][0]<0. || m_p[2][0]<0. || m_p[1][0]<0.)
//     THROW(fatal_error,"Invalid kinematics");
//   if (!IsZero(sqr(((p1+p2)-(m_p[0]+m_p[2]+m_p[1]))[0]/(p1+p2)[0])) ||
//       !IsZero(sqr(((p1+p2)-(m_p[0]+m_p[2]+m_p[1])).Abs2()/(p1+p2).Abs2()))) {
//     msg_Error()<<METHOD<<"(): Momentum conservation violated "
//                <<(p1+p2)-(m_p[0]+m_p[2]+m_p[1])<<"\n"<<
//       " Q initial : "<<p1+p2<<" -> "<<(p1+p2)*(p1+p2)<<"\n"<<
//       " Q final   : "<<m_p[0]+m_p[2]+m_p[1]<<" -> "
//                <<(m_p[0]+m_p[2]+m_p[1])*(m_p[0]+m_p[2]+m_p[1])<<"\n";
//   }
//   for (int i(0);i<3;++i) {
//     xrot.RotateBack(m_p[i]);
//     zrot.RotateBack(m_p[i]);
//     cms.BoostBack(m_p[i]);
//   }
//   if (m_p[0][0]<0. || m_p[2][0]<0. || m_p[1][0]<0.)
//     THROW(fatal_error,"Invalid kinematics");
//   li->SetMom(m_p[0]);
//   lj->SetMom(m_p[1]);
//   lk->SetMom(m_p[2]);
}

