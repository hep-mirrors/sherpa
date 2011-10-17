#include "SHERPA/SoftPhysics/Primordial_KPerp.H"

#include "PDF/Remnant/Remnant_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Primordial_KPerp::Primordial_KPerp(std::string _m_path,std::string _m_file):
  p_filled(new std::set<ATOOLS::Particle*>()),
  p_boosted(new std::set<ATOOLS::Blob*>()),
  m_scheme(0)
{
  p_remnants[0] = p_remnants[1] = NULL;
  p_kperp[0] = new std::vector<Vec3D>();
  p_kperp[1] = new std::vector<Vec3D>();
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(_m_path);
  dataread.SetInputFile(_m_file);
  m_scheme        = dataread.GetValue<int>("K_PERP_SCHEME",0);
  // defaults from Z peak
  double defaultmean1(0.0), defaultmean2(0.0),
    defaultsigma1(0.0), defaultsigma2(0.0);
  if (rpa->gen.Beam1().IsHadron() && rpa->gen.Beam2().IsHadron()) {
  if (rpa->gen.Beam1().Kfcode()==kf_p_plus) {
    defaultmean1=0.917;
    defaultsigma1=0.903;
  }
  if (rpa->gen.Beam2().Kfcode()==kf_p_plus) {
    defaultmean2=0.917;
    defaultsigma2=0.903;
  }
  }
  m_kperpmean[0]  = dataread.GetValue<double>("K_PERP_MEAN_1",defaultmean1);
  m_kperpmean[1]  = dataread.GetValue<double>("K_PERP_MEAN_2",defaultmean2);
  m_kperpsigma[0] = dataread.GetValue<double>("K_PERP_SIGMA_1",defaultsigma1);
  m_kperpsigma[1] = dataread.GetValue<double>("K_PERP_SIGMA_2",defaultsigma2);
}

Primordial_KPerp::~Primordial_KPerp()
{
  delete p_boosted;
  delete p_filled;
  delete p_kperp[1];
  delete p_kperp[0];
}

bool Primordial_KPerp::CreateKPerp(ATOOLS::Blob *blob1,ATOOLS::Blob *blob2)
{
  if (blob1==NULL || blob2==NULL) 
    THROW(critical_error,"Called with NULL pointer.");
  double kpm1=m_kperpmean[0], kpm2=m_kperpmean[1];
  double kps1=m_kperpsigma[0], kps2=m_kperpsigma[1];
  size_t m_maxtrials=1000;
  ATOOLS::Blob *blob[2];
  if (blob1->InParticle(0)->Beam()==0) { blob[0]=blob1; blob[1]=blob2; }
  else { blob[0]=blob2; blob[1]=blob1; }
  p_kperp[0]->resize(blob[0]->NOutP()); 
  p_kperp[1]->resize(blob[1]->NOutP());
  p_boosted->clear();
  double Etot=p_remnants[0]->GetBeam()->Energy()
    +p_remnants[1]->GetBeam()->Energy();
  if (m_scheme==0) {
    bool success=true;
    Vec3D sum[2];
    size_t trials=0;
    do {
      size_t min[2];
      do {
	double ran1, ran2, r12, kperp[2], minimum[2];
	min[1]=min[0]=0; minimum[1]=minimum[0]=1.0e37; sum[1]=sum[0]=Vec3D();
	int pairs=ATOOLS::Max(blob1->NOutP(),blob2->NOutP())-1;
	if (blob1->OutParticle(0)->Flav().IsLepton() ||
	    blob2->OutParticle(0)->Flav().IsLepton()) --pairs;
	for (int i=0;i<pairs;++i) {
	  // generate gaussian numbers
	  do {
	    ran1=2.0*ran->Get()-1.0; ran2=2.0*ran->Get()-1.0;
	    r12=ran1*ran1+ran2*ran2;
	  } while (r12>1.0);
	  r12=sqrt(-2.0*log(r12)/r12);
	  // calculate k_\perp
	  kperp[0]=m_kperpmean[0]+Sign(0.5-ran->Get())*m_kperpsigma[0]*ran1*r12;
	  kperp[1]=m_kperpmean[1]+Sign(0.5-ran->Get())*m_kperpsigma[1]*ran2*r12;
	  for (size_t j=0;j<2;++j) {
	    // generate angle
	    do {
	      ran1=2.0*ran->Get()-1.0; ran2=2.0*ran->Get()-1.0;
	      r12=ran1*ran1+ran2*ran2;
	    } while (r12>1.0);
	    if (i<blob[j]->NOutP()-1) { 
	      (*p_kperp[j])[i]=Vec3D(kperp[j]*(ran1*ran1-ran2*ran2)/
				     r12,kperp[j]*2.0*ran1*ran2/r12,0.0);
	      sum[j]=sum[j]-(*p_kperp[j])[i];
	      if (minimum[j]>dabs(kperp[j])) {
		minimum[j]=dabs(kperp[j]); 
		min[j]=i;
	      }
	    }
	  }
	}
	for (size_t j=0;j<2;++j) (*p_kperp[j])[blob[j]->NOutP()-1]=sum[j];
	success=true;
      } while (!success);
      success=true;
      // sort k_\perp values
      for (size_t i=0;i<2;++i) {
	Vec3D copy=(*p_kperp[i])[blob[i]->NOutP()-1];
	(*p_kperp[i])[blob[i]->NOutP()-1]=(*p_kperp[i])[min[i]];
	(*p_kperp[i])[min[i]]=copy;
	for (int j=0;j<blob[i]->NOutP();++j) {
	  double cur=dabs(blob[i]->OutParticle(j)->Momentum()[3]);
	  for (int k=j;k<blob[i]->NOutP();++k) {
	    if (cur>(*p_kperp[i])[k].Abs()) {
	      copy=(*p_kperp[i])[j]; 
	      (*p_kperp[i])[j]=(*p_kperp[i])[k]; 
	      (*p_kperp[i])[k]=copy;
	      break;
	    }
	    else {
	      if (j==blob[i]->NOutP()-1) success=false;
	    }
	  }
	}
      }
      // test whether Energy and momentum of hard scattering can be preserved
      p_filled->clear();
      int tested=0;
      double E=0.0;
      for (int i=0;i<blob[0]->NOutP();++i) {
      	ATOOLS::Particle *cur2, *cur1=blob[0]->OutParticle(i);
	if (FindConnected(cur1,cur2,true,0)) {
	  ++tested;
	  Vec3D kp1=(*p_kperp[0])[i], kp2=(*p_kperp[1])[i];
	  double s, sp, sp1, sp2;
	  Vec4D cms=cur1->Momentum()+cur2->Momentum();
	  s=cms.Abs2();
	  sp=s+sqr((kp1+kp2).Abs());
	  sp1=cur1->Momentum().Abs2()+sqr(kp1.Abs()); 
	  sp2=cur2->Momentum().Abs2()+sqr(kp2.Abs());
	  E+=sqrt(sp/(1.0-sqr(cms[3]/cms[0])));
	  if (((sp-sp1-sp2)*(sp-sp1-sp2)<4.0*sp1*sp2)||
	      (s<sp1)||(s<sp2)) success=false;
	  p_filled->insert(cur2);
	}
	else {
	  E+=cur1->Momentum()[0];
	  if ((*p_kperp[0])[i].Abs()>
	      ATOOLS::dabs(cur1->Momentum()[3])) success=false;
	}
	p_filled->insert(cur1);
      }
      if (E>Etot) success=false;
      for (int i=0;i<blob[1]->NOutP();++i) {
	ATOOLS::Particle *cur=blob[1]->OutParticle(i);
	if (p_filled->find(cur)==p_filled->end()) {
	  if ((*p_kperp[1])[tested++].Abs()>
	      ATOOLS::dabs(cur->Momentum()[3])) success=false;
	}
      }      
      if ((++trials)==m_maxtrials) {
	for(size_t i=0;i<2;++i) {
	  m_kperpmean[i]/=10.0;
	  m_kperpsigma[i]/=10.0;
	  if (ATOOLS::IsZero(m_kperpmean[i])) m_kperpmean[i]=0.0;
	  if (ATOOLS::IsZero(m_kperpsigma[i])) m_kperpsigma[i]=0.0;
	}
	trials=0;
      }
      // accept or reject
      if ((m_kperpmean[0]==0.0)&&(m_kperpmean[1]==0.0)) success=true;
    } while (!success);
  }
  m_current[1]=m_current[0]=-1;
  m_kperpmean[0]=kpm1; 
  m_kperpmean[1]=kpm2;
  m_kperpsigma[0]=kps1; 
  m_kperpsigma[1]=kps2;
  p_filled->clear();
  m_fill=1;
  return true;
}

bool Primordial_KPerp::BoostConnected(ATOOLS::Blob *blob,unsigned int catcher)
{ 
  if (++catcher>100) {
    msg_Error()<<"ERROR in Primordial_KPerp::BoostConnected(..): "
	       <<"   Blob nesting is too deep."<<std::endl;
    return false;
  }
  if ((blob==NULL)||(p_boosted->find(blob)!=p_boosted->end())) return true;
  p_boosted->insert(blob);
  for (int i=0;i<blob->NOutP();++i) {
    Particle *cur=blob->OutParticle(i);
    if (blob->Type()!=btp::Signal_Process && blob->Type()!=btp::Hard_Decay &&
	blob->Type()!=btp::Hard_Collision &&
	(cur->DecayBlob()==NULL || cur->DecayBlob()->Type()!=btp::Signal_Process) &&
	(cur->DecayBlob()==NULL || cur->DecayBlob()->Type()!=btp::Hard_Decay) &&
	(cur->DecayBlob()==NULL || cur->DecayBlob()->Type()!=btp::Hard_Collision)) {
      Vec4D mom=cur->Momentum();
      m_oldcms.Boost(mom);
      m_rotate.Rotate(mom);
      m_newcms.BoostBack(mom);
      cur->SetMomentum(mom);
    }
    if (!BoostConnected(cur->DecayBlob(),catcher)) return false;
  }
  return true;
}

bool Primordial_KPerp::FindConnected(ATOOLS::Particle *particle,ATOOLS::Particle *&connected,
                                     bool forward,unsigned int catcher)
{
  if (++catcher>100) {
    msg_Error()<<"ERROR in Primordial_KPerp::FindConnected(..): "
	       <<"   Blob nesting is too deep."<<std::endl;
    return false;
  }
  if (!forward) {
    Blob *prod=particle->ProductionBlob();
    if (prod!=NULL) {
      if (prod->Type()==btp::Beam) {
        connected=particle;
        return true;
      }
      for (int i=0;i<prod->NInP();++i) {
        if (FindConnected(prod->InParticle(i),connected,false,catcher)) return true;
      }
    }
    else {
      connected=particle;
      return true;
    }
  }
  else {
    Blob *decy=particle->DecayBlob();
    if (decy!=NULL) {
      for (int i=0;i<decy->NInP();++i) {
	Particle *next=decy->InParticle(i);
        if (next!=particle && 
	    next->ProductionBlob()->Type()!=btp::Signal_Process &&
	    next->ProductionBlob()->Type()!=btp::Hard_Decay &&
	    next->ProductionBlob()->Type()!=btp::Hard_Collision) 
	  if (FindConnected(next,connected,false,catcher)) return true;
      }
      THROW(fatal_error,"Inconsistent blob structure");
    }
  }
  return false;
}

double Primordial_KPerp::Lambda2(double sp,double sp1,double sp2) 
{ 
  return (sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2;
}

void Primordial_KPerp::FillKPerp(ATOOLS::Particle *cur1,unsigned int beam)
{
  if (m_kperpmean[0]==0.0 && m_kperpmean[1]==0.0) {
    return;
  }
  if (p_filled->find(cur1)!=p_filled->end()) return;
  ++m_current[beam];
  Vec3D kp1;
  Vec4D mom1, old1=cur1->Momentum();
  kp1=(*p_kperp[beam])[m_current[beam]]+Vec3D(old1[1],old1[2],0.0);
  Particle *cur2;
  if (!FindConnected(cur1,cur2,true,0)) {
    mom1=Vec4D(old1[0],kp1[1],kp1[2],
	       Sign(old1[3])*sqrt(old1[3]*old1[3]-kp1[1]*kp1[1]-kp1[2]*kp1[2])); 
    cur1->SetMomentum(mom1); 
    p_filled->insert(cur1);
    return;
  }
  if (cur1->Flav().IsLepton() || cur2->Flav().IsLepton()) return;
  ++m_current[1-beam];
  Vec4D mom2, old2=cur2->Momentum(), oldcms=old1+old2;
  Vec3D kp2=(*p_kperp[1-beam])[m_current[1-beam]]+Vec3D(old2[1],old2[2],0.0);
  m_oldcms=Poincare(oldcms);
  double sp, sp1, sp2, Enew, pznew, E1, E2, pz1, pz2;
  sp1=old1.Abs2()+sqr(kp1.Abs());
  sp2=old2.Abs2()+sqr(kp2.Abs());
  sp=oldcms.Abs2()+sqr((kp1+kp2).Abs());
  Enew=sqrt(sp/(1.0-sqr(oldcms[3]/oldcms[0])));
  pznew=sqrt(sp*sqr(oldcms[3]/oldcms[0])/(1.0-sqr(oldcms[3]/oldcms[0])));
  double yto=(oldcms[0]+oldcms[3])/(oldcms[0]-oldcms[3]);
  double spo=oldcms.Abs2();
  for (double sign=1.0;sign>=-1.0;sign-=2.0) {
    E1=0.5/sp*((sp+sp1-sp2)*Enew+sign*sqrt(Lambda2(sp,sp1,sp2))*pznew);
    E2=Enew-E1;
    pz1=Sign(old1[3])*sqrt(E1*E1-sp1);
    pz2=Sign(old2[3])*sqrt(E2*E2-sp2);
    double spn1=sqr(E1+E2)-sqr(pz1+pz2)-sqr((kp1+kp2).Abs());
    double spn2=sqr(E1+E2)-sqr(-pz1+pz2)-sqr((kp1+kp2).Abs());
    if (ATOOLS::dabs(spn1-spo)>ATOOLS::dabs(spn2-spo)) pz1*=-1.0;
    spn1=sqr(E1+E2)-sqr(pz1+pz2)-sqr((kp1+kp2).Abs());
    spn2=sqr(E1+E2)-sqr(-pz1+pz2)-sqr((kp1+kp2).Abs());
    double ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2);
    if (ATOOLS::dabs(ytn-yto)>ATOOLS::dabs(1./ytn-yto)) { pz1*=-1.0; pz2*=-1.0; }
    ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2);
    if (dabs(pz1)>dabs(pz2)) { if (Sign(pz1)==Sign(old1[3])) break; }
    else { if (Sign(pz2)==Sign(old2[3])) break; }
  }
  mom1=Vec4D(E1,kp1[1],kp1[2],pz1);
  mom2=Vec4D(E2,kp2[1],kp2[2],pz2);
  m_newcms=Poincare(mom1+mom2);
  cur1->SetMomentum(mom1); 
  cur2->SetMomentum(mom2);
  m_newcms.Boost(mom1); 
  m_newcms.Boost(mom2);
  if (mom2[3]>0.0) m_rotate=Poincare(Vec4D::ZVEC,mom2);
  else m_rotate=Poincare(Vec4D::ZVEC,mom1);
  BoostConnected(cur1->DecayBlob(),0);
  BoostConnected(cur2->DecayBlob(),0);
  p_filled->insert(cur2);
  p_filled->insert(cur1);
  if (!(cur1->Momentum()[0]>0.)) {
    msg_Tracking()<<"Primordial_KPerp::FillKPerp(..): "
			  <<"Parton ("<<cur1->Number()<<") has non-positive energy "
			  <<cur1->Momentum()<<std::endl;
  }
  if (!(cur2->Momentum()[0]>0.)) {
    msg_Tracking()<<"Primordial_KPerp::FillKPerp(..): "
			  <<"Parton ("<<cur2->Number()<<") has non-positive energy "
			  <<cur2->Momentum()<<std::endl;
  }
  return;
}

void Primordial_KPerp::FillKPerp(ATOOLS::Blob *blob)
{
  unsigned int beam=blob->InParticle(0)->Beam();
  if (beam==0) {
    m_current[1]=m_current[0]=-1;
    p_filled->clear();
  }
  if (m_fill) 
    for (int i=0;i<blob->NOutP();++i) 
      FillKPerp(blob->OutParticle(i),beam);
}
