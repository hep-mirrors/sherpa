#include "CFPSHOWER++/Shower/Cluster_Definitions.H"
#include "CFPSHOWER++/Shower/Shower.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

const double s_uxeps=1.0e-3;

Cluster_Definitions::Cluster_Definitions(Shower * const shower):
  p_shower(shower), m_mode(0), m_amode(0) {}


Cluster_Param Cluster_Definitions::Cluster(const Cluster_Config & ca)
{
  msg_Out()<<METHOD<<":\n"<<ca<<".\n";
  DEBUG_FUNC(ca);
  p_shower->SetMassSelector(p_ms = ca.p_ms);
  int i(ca.m_i), j(ca.m_j), swap(j<ca.p_ampl->NIn() && j<i);
  if (swap) std::swap<int>(i,j);
  double ws, mu2;
  kernel_type::code type = GetCode(i<ca.p_ampl->NIn(),ca.m_k<ca.p_ampl->NIn());
  Splitting split(KT2(*ca.p_ampl,i,j,ca.m_k,ca.m_mo,
		      ca.m_kin,int(type),(swap?2:0)|(ca.m_mode<<2),ws,mu2));
  bool iss = (i<ca.p_ampl->NIn() || j<ca.p_ampl->NIn());
  if (split.t()>0.0)
    return Cluster_Param(this,ws,split.t(),mu2,0,split.KinScheme(),0,
			 (iss?-1.:1.)*split.GetKinematics()->m_pi,
			 (ca.m_k<ca.p_ampl->NIn()?-1.:1.)*split.GetKinematics()->m_pk,
			 split.GetKinematics()->m_lam);
  if (ca.PureQCD()) return Cluster_Param(this,0.0,0.0,0.0,0);
  return Cluster_Param(this,0.0);
}


Splitting Cluster_Definitions::KT2(const ATOOLS::Cluster_Amplitude &ampl,
				   int i,int j,int k, const ATOOLS::Flavour &mo,
				   const int kin,const int type,const int mode,
				   double &ws,double &mu2)
{
  const ATOOLS::Cluster_Leg * li=ampl.Leg(i), * lj=ampl.Leg(j), * lk=ampl.Leg(k);
  Parton * out1(GetParton(li)), * out2(GetParton(lj)), * spec(GetParton(lk));  // c n s
  Splitting splitting;
  splitting.SetParton(0,out1);
  splitting.SetParton(1,out2);
  splitting.SetSpectator(spec);
  splitting.SetMom(0,li->Mom());
  splitting.SetMom(1,lj->Mom());
  splitting.SetSpecMom(lk->Mom());
  splitting.Set_eta(out1->XB());
  splitting.SetKinScheme(kin);
  splitting.SetClustered(1);
  /*
  splitting.m_kin=kin>=0?kin:p_shower->KinematicsScheme();
  splitting.m_type=type;
  splitting.m_cpl=p_shower->CouplingScheme();
  splitting.m_kfac=(mode&(16<<2))?0:p_shower->KFactorScheme();
  Kernel *sk(p_shower->GetKernel(splitting,(mode&2)?1:0));
  ws=0.0;
  if (sk==NULL) return Splitting(NULL,NULL,-1.0);
  if (!sk->LF()->SetLimits(splitting) ||
      !sk->LF()->Cluster(splitting,1|2)) {
    splitting.m_t=-1.0;
    return splitting;
  }
  msg_Debugging()<<"Splitting: t = "<<splitting.m_t<<" = "<<sqrt(splitting.m_t)
		 <<" ^ 2, z = "<<splitting.m_z<<", phi = "<<splitting.m_phi<<"\n"; 
  ws=sk->Value(splitting);
  msg_Debugging()<<"Kernel: "<<ws<<" ( kfac = "<<splitting.m_kfac
		 <<" )  <-  "<<sk->Class()<<"\n";
  mu2=sk->GF()->Scale(splitting);
  if (p_shower->KFactorScheme() &&
      splitting.m_t>p_shower->TMin(type&1)) {
    splitting.m_kfac=0;
    double K=ws/sk->Value(splitting);
    msg_Debugging()<<"     K: "<<K<<" ( kfac = "<<splitting.m_kfac<<" )\n";
    if (K>0.0 && !IsEqual(K,1.0)) {
      splitting.m_clu=0;
      mu2=sk->GF()->Solve(K*sk->GF()->Coupling(splitting));
    }
  }
  if (ws) ws=ws*splitting.m_Q2/splitting.m_t;
  return splitting;
  */
}


Flavour Cluster_Definitions::ProperFlav(const Flavour & flav) const
{
  Flavour properflav(flav);
  switch (properflav.Kfcode()) {
  case kf_gluon_qgc:
    properflav=Flavour(kf_gluon);
    break;
  default: break;
  }
  return properflav;
}
