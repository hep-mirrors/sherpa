#include "Primordial_KPerp.H"

#include "Run_Parameter.H"
#include "Data_Read.H"
#include "Random.H"
#include "Vector.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Primordial_KPerp::Primordial_KPerp(std::string _m_path,std::string _m_file):
  p_filled(new std::set<ATOOLS::Particle*>()),
  p_boosted(new std::set<ATOOLS::Blob*>()),
  m_scheme(0)
{
  p_kperp[0] = new std::vector<Vec3D>();
  p_kperp[1] = new std::vector<Vec3D>();
  m_current[1]=m_current[0]=1;
  Data_Read *dataread = new Data_Read(_m_path+_m_file);
  m_scheme=dataread->GetValue<int>("K_PERP_SCHEME",0);
  if (rpa.gen.Beam1().IsHadron() && rpa.gen.Beam2().IsHadron()) { 
    cout<<"Beam1 is hadron"<<endl; 
    m_kperpmean[0]=dataread->GetValue<double>("K_PERP_MEAN_1",0.8);
    m_kperpmean[1]=dataread->GetValue<double>("K_PERP_MEAN_2",0.8);
    m_kperpsigma[0]=dataread->GetValue<double>("K_PERP_SIGMA_1",m_kperpmean[0]);
    m_kperpsigma[1]=dataread->GetValue<double>("K_PERP_SIGMA_2",m_kperpmean[1]);
  }
  else { 
    cout<<"Beam1 is not a hadron"<<endl; 
    m_kperpmean[0]=m_kperpmean[1]=0.0; 
    m_kperpsigma[0]=m_kperpsigma[1]=0.0;
  }
  delete dataread;
}

Primordial_KPerp::~Primordial_KPerp()
{
  delete p_boosted;
  delete p_filled;
  delete p_kperp[1];
  delete p_kperp[0];
}

bool Primordial_KPerp::CreateKPerp(ATOOLS::Blob *blob,unsigned int nparticle)
{
  unsigned int beam=blob->Beam();
  //   if (m_kperpmean[beam]==0.0) return true;
  double m_cut=10.0;
  p_kperp[beam]->resize(nparticle);
  if (beam==0) {
    p_filled->clear();
    p_boosted->clear();
    p_kperp[1-beam]->clear();
  }
  std::vector<std::pair<Vec4D,Vec4D> > connected;
  if (m_scheme==0) {
    bool success;
    connected.clear();
    for (int i=0;i<blob->NOutP();++i) {
      ATOOLS::Particle *cur2, *cur1=blob->OutParticle(i);
      FindConnected(cur1,cur2,true,0);
      connected.push_back(std::pair<Vec4D,Vec4D>(cur1->Momentum(),cur2->Momentum())); 
      if (beam==0) p_kperp[1-beam]->push_back(ATOOLS::Vec3D());
    }
    Vec3D sum, res;
    do {
      success=true;
      do {
	sum=Vec3D();
	double ran1, ran2, ran3, ran4, r122, r342, kperp, next, min=1.0e37;
	for (unsigned int i=0;i<nparticle-1;++i) {
	  next=0.0;
	  do { 
	    if (m_kperpmean[beam]!=0.0) {
	      if (next==0.0) {
		do {
		  ran1=2.0*ran.Get()-1.0; ran2=2.0*ran.Get()-1.0;
		  r122=ran1*ran1+ran2*ran2;
	      } while (r122>1);
		r122=sqrt(-2.0*log(r122)/r122);
		kperp=m_kperpmean[beam]+Sign(0.5-ran.Get())*m_kperpsigma[beam]*ran1*r122;
		next=m_kperpmean[beam]+Sign(0.5-ran.Get())*m_kperpsigma[beam]*ran2*r122;
	      }
	      else {
		kperp=next;
		next=0.0;
	      }
	    }
	    else kperp=0.0;
	    do {
	      ran3=2.0*ran.Get()-1.0; ran4=2.0*ran.Get()-1.0;
	      r342=ran3*ran3+ran4*ran4;
	    } while (r342>1);
	    (*p_kperp[beam])[i]=res=Vec3D(kperp*(ran3*ran3-ran4*ran4)/r342,
					  kperp*2.0*ran3*ran4/r342,0.0);
	    if (min>dabs(kperp)) {
	      min=dabs(kperp);
	      (*p_kperp[beam])[i]=(*p_kperp[beam])[0];
	      (*p_kperp[beam])[0]=res;
	    }
	  } while (kperp>m_cut);
	  sum=sum-res;
	  if ((i<connected.size())&&(m_current[1-beam]==0)) {
	    Vec3D kp1=(*p_kperp[beam])[i], kp2=(*p_kperp[1-beam])[i];
	    double sp, sp1, sp2;
	    sp1=connected[i].first.Abs2()+sqr(kp1.Abs()); 
	    sp2=connected[i].second.Abs2()+sqr(kp2.Abs());
	    sp=(connected[i].first+connected[i].second).Abs2()+sqr((kp1+kp2).Abs());
	    if (((sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2)<0.0) success=false;
	  }
	}
	(*p_kperp[beam])[nparticle-1]=sum;
      } while (exp(-0.5*sqr((m_kperpmean[beam]-sum.Abs())/m_kperpsigma[beam]))<ran.Get());
      m_current[beam]=0;
    } while (!success);
  }
  return true;
}

bool Primordial_KPerp::BoostConnected(ATOOLS::Blob *blob,unsigned int catcher)
{ 
  if (++catcher>100) {
    msg.Error()<<"Primordial_KPerp::BoostConnected(..): "
		       <<"Blob nesting is too deep!"<<std::endl
		       <<"   Cannot boost connected parton."<<std::endl;
    return false;
  }
  if ((blob==NULL)||(p_boosted->find(blob)!=p_boosted->end())) return true;
  p_boosted->insert(blob);
  for (int i=0;i<blob->NOutP();++i) {
    Particle *cur=blob->OutParticle(i);
    Vec4D mom=cur->Momentum();
    m_oldcms.Boost(mom);
    m_rotate.Rotate(mom);
    m_newcms.BoostBack(mom);
    cur->SetMomentum(mom);
    if (!BoostConnected(cur->DecayBlob(),catcher)) return false;
  }
  return true;
}

bool Primordial_KPerp::FindConnected(ATOOLS::Particle *particle,ATOOLS::Particle *&connected,
                                     bool forward,unsigned int catcher)
{
  if (++catcher>100) {
    msg.Error()<<"Primordial_KPerp::FindConnected(..): "
                       <<"Blob nesting is too deep!"<<std::endl
                       <<"   Cannot find connected parton."<<std::endl;
    return false;
  }
  if (!forward) {
    Blob *prod=particle->ProductionBlob();
    if (prod!=NULL) {
      if (prod->Type().find("Beam Remnant")!=std::string::npos) {
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
        if (next!=particle) if (FindConnected(next,connected,false,catcher)) return true;
      }
      for (int i=0;i<decy->NOutP();++i) {
        if (FindConnected(decy->OutParticle(i),connected,true,catcher)) return true;
      }
    }
  }
  return false;
}

bool Primordial_KPerp::FillKPerp(ATOOLS::Particle *cur1,unsigned int beam)
{
  if ((m_kperpmean[0]==0.0)&&(m_kperpmean[1]==0.0)) return true;
  if (p_filled->find(cur1)!=p_filled->end()) return true;
  ++m_current[beam];
  Vec3D kp1; 
  Vec4D mom1, old1=cur1->Momentum();
  size_t i=0;
  if (cur1->Flav().IsDiQuark()) {
    kp1=(*p_kperp[beam])[0];
    if (dabs(old1[3])>kp1.Abs()) --m_current[beam];
    else i=p_kperp[beam]->size();
  }
  else {
    for (i=m_current[beam];i<p_kperp[beam]->size();++i) {
      kp1=(*p_kperp[beam])[i];
      if (dabs(old1[3])>kp1.Abs()) {
	(*p_kperp[beam])[i]=(*p_kperp[beam])[m_current[beam]];
	break; 
      }
    }
  }
  if (i==p_kperp[beam]->size()) {
    ATOOLS::msg.Out()<<"Primordial_KPerp::FillKPerp(..): "
		     <<"Cannot assign k_\\perp to particle! Reject k_\\perp set."<<std::endl;
    return false;
  }
  Particle *cur2;
  if (!FindConnected(cur1,cur2,true,0)) {
    mom1=Vec4D(old1[0],kp1[1],kp1[2],
	       Sign(old1[3])*sqrt(old1[3]*old1[3]-kp1[1]*kp1[1]-kp1[2]*kp1[2])); 
    cur1->SetMomentum(mom1); 
    p_filled->insert(cur1);
    return true;
  }
  ++m_current[1-beam];
  Vec3D kp2;
  Vec4D mom2, old2=cur2->Momentum(), oldcms=old1+old2;
  for (i=m_current[1-beam];i<p_kperp[1-beam]->size();++i) {
    kp2=(*p_kperp[1-beam])[i];
    if (dabs(old2[3])>kp2.Abs()) {
      (*p_kperp[1-beam])[i]=(*p_kperp[1-beam])[m_current[1-beam]];
      break; 
    }
  }
  if (i==p_kperp[1-beam]->size()) {
    ATOOLS::msg.Tracking()<<"Primordial_KPerp::FillKPerp(..): "
			  <<"Cannot assign k_\\perp to particle! Reject k_\\perp set."<<std::endl;
    return false;
  }
  m_oldcms=Poincare(oldcms);
  double sp, sp1, sp2, lam, pf2, E1, E2, pz1, pz2;
  sp1=old1.Abs2()+sqr(kp1.Abs());
  sp2=old2.Abs2()+sqr(kp2.Abs());
  sp=oldcms.Abs2()+sqr((kp1+kp2).Abs());
  pf2=sp/(1.0-sqr(oldcms[3]/oldcms[0]));
  lam=sqrt((sp-sp1-sp2)*(sp-sp1-sp2)/4.0-sp1*sp2);
  double ytn, yto=(oldcms[0]+oldcms[3])/(oldcms[0]-oldcms[3]);
  double spn, spo=oldcms.Abs2();
  for (double sign=1.0;sign>=-1.0;sign-=2.0) {
    E1=(sqrt(pf2)*(sp-sp2+sp1)/2.0+sign*lam*sqrt(pf2-sp))/sp;
    E2=sqrt(pf2)-E1;
    pz1=Sign(old1[3])*sqrt(E1*E1-sp1);
    pz2=Sign(old2[3])*sqrt(E2*E2-sp2);
    spn=pf2-sqr(pz1+pz2)-sqr((kp1+kp2).Abs());
    if (!IsEqual(spn,spo)) pz1*=-1.0;
    ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2);
    if (!IsEqual(ytn,yto)) { pz1*=-1.0; pz2*=-1.0; }
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
  return true;
}

