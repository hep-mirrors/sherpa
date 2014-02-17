#include "PDF/Main/Jet_Criterion.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

namespace PYTHIA {

  class Pythia_Jet_Criterion: public Jet_Criterion {
  private:

  public:

    Pythia_Jet_Criterion(const std::string &args)
    { }

    bool Jets(Cluster_Amplitude *ampl,int mode)
    {
      DEBUG_FUNC("mode = "<<mode);
      msg_Debugging()<<*ampl<<"\n";
      PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
      double q2cut(jf->Ycut()*sqr(rpa->gen.Ecms()));
      double q2min(std::numeric_limits<double>::max());
      size_t imin(0), jmin(0), kmin(0);
      Flavour mofl;
      for (size_t i(0);i<ampl->Legs().size();++i) {
	Cluster_Leg *li(ampl->Leg(i));
	for (size_t j(Max(i+1,ampl->NIn()));j<ampl->Legs().size();++j) {
	  Cluster_Leg *lj(ampl->Leg(j));
	  if (j<ampl->NIn()) continue;
	  for (size_t k(0);k<ampl->Legs().size();++k) {
	    if (k==i || k==j) continue;
	    Cluster_Leg *lk(ampl->Leg(k));
	    if (i<ampl->NIn() && k>=ampl->NIn()) continue;
	    if (lk->Flav().Strong() &&
		li->Flav().Strong() && lj->Flav().Strong()) {
	      if (i<ampl->NIn()) li->SetMom(-li->Mom());
	      if (k<ampl->NIn()) lk->SetMom(-lk->Mom());
	      double q2ijk(pT2pythia(ampl,*li,*lj,*lk,i<ampl->NIn()?-1:1));
	      msg_Debugging()<<"Q_{"<<ID(li->Id())<<ID(lj->Id())
			     <<","<<ID(lk->Id())<<"} = "<<sqrt(q2ijk)<<"\n";
	      if (i<ampl->NIn()) li->SetMom(-li->Mom());
	      if (k<ampl->NIn()) lk->SetMom(-lk->Mom());
	      if (mode==0) {
		if (q2ijk<q2cut) return false;
	      }
	      else {
		if (q2ijk<q2min) {
		  q2min=q2ijk;
		  mofl=Flavour(kf_gluon);
		  if (li->Flav().IsGluon()) mofl=lj->Flav();
		  if (lj->Flav().IsGluon()) mofl=li->Flav();
		  imin=i;
		  jmin=j;
		  kmin=k;
		}
	      }
	    }
	  }
	}
      }
      if (mode!=0 && imin!=jmin) {
	Vec4D_Vector p=Combine(*ampl,imin,jmin,kmin,mofl);
	if (p.empty()) {
	  msg_Error()<<METHOD<<"(): Combine failed. Use R configuration."<<std::endl;
	  return Jets(ampl,0);
	}
	Cluster_Amplitude *bampl(Cluster_Amplitude::New());
	bampl->SetProc(ampl->Proc<void>());
	bampl->SetMS(ampl->MS());
	bampl->SetNIn(ampl->NIn());
	bampl->SetJF(ampl->JF<void>());
	for (int i(0), j(0);i<ampl->Legs().size();++i) {
	  if (i==jmin) continue;
	  if (i==imin) {
	    bampl->CreateLeg(p[j],mofl,ampl->Leg(i)->Col());
	    bampl->Legs().back()->SetId(ampl->Leg(imin)->Id()|ampl->Leg(jmin)->Id());
	    bampl->Legs().back()->SetK(ampl->Leg(kmin)->Id());	
	  }
	  else {
	    bampl->CreateLeg(p[j],ampl->Leg(i)->Flav(),ampl->Leg(i)->Col());
	  }
	  ++j;
	}
	bool res=Jets(bampl,0);
	bampl->Delete();
	return res;
      }
      msg_Debugging()<<"--- Jet veto ---\n";
      return true;
    }

    double pT2pythia(Cluster_Amplitude* ampl, const Cluster_Leg& RadAfterBranch,
                const Cluster_Leg& EmtAfterBranch,
                const Cluster_Leg& RecAfterBranch, int ShowerType){
      // Save type: 1 = FSR pT definition, else ISR definition
      int Type   = ShowerType;
      // Calculate virtuality of splitting
      int sign = (Type==1) ? 1 : -1;
      Vec4D Q(RadAfterBranch.Mom() + sign*EmtAfterBranch.Mom());
      double Qsq = sign * Q.Abs2();
      // Mass term of radiator
      DEBUG_VAR(ampl->MS());
      double m2Rad = ( abs(RadAfterBranch.Flav().Kfcode()) >= 4
                   && abs(RadAfterBranch.Flav().Kfcode()) < 7)
                   ? ampl->MS()->Mass2(RadAfterBranch.Flav())
                   : 0.;
      // Construct 2->3 variables for FSR
      Vec4D  sum     = RadAfterBranch.Mom() + RecAfterBranch.Mom()
                     + EmtAfterBranch.Mom();
      double m2Dip = sum.Abs2();
      double x1 = 2. * (sum * RadAfterBranch.Mom()) / m2Dip;
      double x3 = 2. * (sum * EmtAfterBranch.Mom()) / m2Dip;
      // Construct momenta of dipole before/after splitting for ISR 
      Vec4D qBR(RadAfterBranch.Mom() - EmtAfterBranch.Mom() + RecAfterBranch.Mom());
      Vec4D qAR(RadAfterBranch.Mom() + RecAfterBranch.Mom()); 
      // Calculate z of splitting, different for FSR and ISR
      double z = (Type==1) ? x1/(x1+x3)
                         : (qBR.Abs2())/( qAR.Abs2());
      // Separation of splitting, different for FSR and ISR
      double pTpyth = (Type==1) ? z*(1.-z) : (1.-z);
      // pT^2 = separation*virtuality
      pTpyth *= (Qsq - sign*m2Rad);
      if(pTpyth < 0.) pTpyth = 0.;
      // Return pT
      return pTpyth;
    }

    ATOOLS::Vec4D_Vector  Combine
    (const Cluster_Amplitude &ampl,int i,int j,int k,const ATOOLS::Flavour &mo)
    {
      Mass_Selector *p_ms=ampl.MS();
      if (i>j) std::swap<int>(i,j);
      Vec4D_Vector after(ampl.Legs().size()-1);
      double mb2(0.0);
      if (i<2) {
	mb2=ampl.Leg(1-i)->Mom().Abs2();
	double mfb2(p_ms->Mass2(ampl.Leg(1-i)->Flav()));
	if ((mfb2==0.0 && IsZero(mb2,1.0e-6)) || IsEqual(mb2,mfb2,1.0e-6)) mb2=mfb2;
      }
      Vec4D pi(ampl.Leg(i)->Mom()), pj(ampl.Leg(j)->Mom());
      Vec4D pk(ampl.Leg(k)->Mom()), pb(i<2?ampl.Leg(1-i)->Mom():Vec4D());
      double mi2=pi.Abs2(), mfi2=p_ms->Mass2(ampl.Leg(i)->Flav());
      double mj2=pj.Abs2(), mfj2=p_ms->Mass2(ampl.Leg(j)->Flav());
      double mk2=pk.Abs2(), mfk2=p_ms->Mass2(ampl.Leg(k)->Flav());
      if ((mfi2==0.0 && IsZero(mi2,1.0e-6)) || IsEqual(mi2,mfi2,1.0e-6)) mi2=mfi2;
      if ((mfj2==0.0 && IsZero(mj2,1.0e-6)) || IsEqual(mj2,mfj2,1.0e-6)) mj2=mfj2;
      if ((mfk2==0.0 && IsZero(mk2,1.0e-6)) || IsEqual(mk2,mfk2,1.0e-6)) mk2=mfk2;
      double mij2=p_ms->Mass2(mo);
      Kin_Args lt;
      if (i>1) {
	if (k>1) lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,2);
	else lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,2);
	if ((k==0 && lt.m_pk[3]<0.0) ||
	    (k==1 && lt.m_pk[3]>0.0) || lt.m_pk[0]<0.0) return Vec4D_Vector();
      }
      else {
	if (k>1) {
	  lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-pb,2);
	}
	else lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,2);
	if ((i==0 && lt.m_pi[3]<0.0) ||
	    (i==1 && lt.m_pi[3]>0.0) || lt.m_pi[0]<0.0) return Vec4D_Vector();
      }
      if (lt.m_stat<0) return Vec4D_Vector();
      for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
	if (m==(size_t)j) continue;
	if (m==(size_t)i) after[l]=i>1?lt.m_pi:-lt.m_pi;
	else if (m==(size_t)k) after[l]=k>1?lt.m_pk:-lt.m_pk;
	else after[l]=lt.m_lam*ampl.Leg(m)->Mom();
	++l;
      }
      return after;
    }

  };// end of class Pythia_Jet_Criterion

}

using namespace PYTHIA;

DECLARE_GETTER(Pythia_Jet_Criterion,"PYTHIA",
	       Jet_Criterion,JetCriterion_Key);
Jet_Criterion *ATOOLS::Getter
<Jet_Criterion,JetCriterion_Key,Pythia_Jet_Criterion>::
operator()(const JetCriterion_Key &args) const
{ return new Pythia_Jet_Criterion(args.m_key); }
void ATOOLS::Getter
<Jet_Criterion,JetCriterion_Key,Pythia_Jet_Criterion>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"The Pythia jet criterion"; }
