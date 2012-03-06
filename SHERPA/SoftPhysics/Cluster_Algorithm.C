#include "SHERPA/SoftPhysics/Cluster_Algorithm.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <limits>

using namespace SHERPA;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm():
  p_ampl(NULL), p_clus(NULL), p_jf(NULL) {}

Cluster_Algorithm::~Cluster_Algorithm()
{
  if (p_jf) delete p_jf;
  //if (p_ampl!=NULL) p_ampl->Delete();
}

int Cluster_Algorithm::ColorConnected(const ColorID &i,const ColorID &j) const
{
  return int(i.m_i==j.m_j && i.m_i!=0)+int(i.m_j==j.m_i && i.m_j!=0);
}

void Cluster_Algorithm::
ProjectOnSinglets(Blob * const blob,std::list<ParticleList *> & singlets) {
  //msg_Out()<<METHOD<<" for \n"<<(*blob)<<"\n";
  ParticleList outs, * sing;
  std::list<Particle * >::iterator piter;
  for (int i(0);i<blob->NOutP();++i) outs.push_back(blob->OutParticle(i));
  int col1,col2;
  bool add;
  while (!outs.empty()) {
    sing  = new std::list<Particle *>;
    if (outs.size()==0) break;
    col1  = col2  = 0;
    add   = false;
    msg_Tracking()<<"++++ Start 3 list at begin of out-particle list: "
	     <<outs.size()<<" particles left.\n";
    if (!add && !outs.empty()) {
      piter = outs.begin();
      do {
	msg_Tracking()<<"   Test "<<(*piter)->Flav()<<" "
		 <<"["<<(*piter)->GetFlow(1)<<", "
		 <<(*piter)->GetFlow(2)<<"], "
		 <<"add = "<<add;
	if (sing->empty() && 
	    (*piter)->GetFlow(1)!=0 && (*piter)->GetFlow(2)==0) {
	  col1 = (*piter)->GetFlow(1);
	  col2 = (*piter)->GetFlow(2);
	  sing->push_back((*piter));
	  msg_Tracking()<<" --> start singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else if (col1!=0 && col2==0 && (*piter)->GetFlow(2)==col1) {
	  col1 = (*piter)->GetFlow(1);
	  sing->push_back((*piter));
	  msg_Tracking()<<" --> add to singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else {
	  msg_Tracking()<<" --> ignore.\n";
	  add = false;
	  piter++;
	}
      } while (add && piter!=outs.end() && !outs.empty());
    }
    msg_Tracking()<<"++++ Start anti-3 list at begin of out-particle list: "
	     <<outs.size()<<" particles left, add = "<<add<<".\n";
    if (!add && !outs.empty()) {
      piter = outs.begin();
      do {
	msg_Tracking()<<"   Test "<<(*piter)->Flav()<<" "
		 <<"["<<(*piter)->GetFlow(1)<<", "
		 <<(*piter)->GetFlow(2)<<"], "
		 <<"add = "<<add;
	if (sing->empty() && 
	    (*piter)->GetFlow(1)==0 && (*piter)->GetFlow(2)!=0) {
	  col1 = (*piter)->GetFlow(1);
	  col2 = (*piter)->GetFlow(2);
	  sing->push_front((*piter));
	  msg_Tracking()<<" --> start singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else if (col1==0 && col2!=0 && (*piter)->GetFlow(1)==col2) {
	  col2 = (*piter)->GetFlow(2);
	  sing->push_front((*piter));
	  msg_Tracking()<<" --> add to singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else {
	  msg_Tracking()<<" --> ignore.\n";
	  add = false;
	  piter++;
	}
      } while (add && piter!=outs.end() && !outs.empty());
    }
    msg_Tracking()<<"++++ Start 8 list at begin of out-particle list: "
	     <<outs.size()<<" particles left, add = "<<add<<".\n";
    if (!add && !outs.empty()) {
      piter = outs.begin();
      do {
	msg_Tracking()<<"   Test "<<(*piter)->Flav()<<" "
		 <<"["<<(*piter)->GetFlow(1)<<", "
		 <<(*piter)->GetFlow(2)<<"], "
		 <<"add = "<<add;
	if (sing->empty() && 
	    (*piter)->GetFlow(1)!=0 && (*piter)->GetFlow(2)!=0) {
	  col1 = (*piter)->GetFlow(1);
	  col2 = (*piter)->GetFlow(2);
	  sing->push_front((*piter));
	  msg_Tracking()<<" --> start singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else if (col1!=0 && col2!=0 && (*piter)->GetFlow(1)==col2) {
	  col2 = (*piter)->GetFlow(2);
	  sing->push_front((*piter));
	  msg_Tracking()<<" --> add to singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else if (col1!=0 && col2!=0 && (*piter)->GetFlow(2)==col1) {
	  col1 = (*piter)->GetFlow(1);
	  sing->push_back((*piter));
	  msg_Tracking()<<" --> add to singlet.\n";
	  outs.erase(piter);
	  if (!outs.empty()) piter = outs.begin();
	  add = true;
	}
	else {
	  msg_Tracking()<<" --> ignore.\n";
	  add = false;
	  piter++;
	}
      } while (add && !(piter==outs.end() || outs.empty()));
    }
    if (sing->empty()) {
      delete sing;
      break;
    }
    else singlets.push_back(sing);
    //msg_Out()<<(*sing)<<"\n";
  }
}

double Cluster_Algorithm::
PT2(const Vec4D & pi,const Vec4D & pj,const bool & beam) const 
{
  double pref(1.);
  double pti2(pi.PPerp2()), ptj2(pj.PPerp2());
  double ptij2    = 
    2.*Min(pti2,ptj2)*(cosh(pi.Eta()-pj.Eta())-cos(pi.Phi()-pj.Phi()));
  return (beam?Min(pref*pti2,ptij2):ptij2);
}

bool Cluster_Algorithm::Cluster(Blob *const blob)
{
  //msg_Out()<<"===========================================\n"
  // 	   <<"===========================================\n"
  // 	   <<"===========================================\n"
  //	   <<METHOD<<"(kt^2 = "<<m_minkt2<<")for:\n"<<(*blob)<<"\n";
  std::list<ParticleList * > singlets;
  ProjectOnSinglets(blob,singlets);

  p_ampl=Cluster_Amplitude::New(NULL);
  Vec4D axis1(blob->GetParticle(0)->Momentum());
  Vec4D axis2(blob->GetParticle(1)->Momentum());
  for (int i(0);i<blob->NInP();++i) {
    Particle *const copy(blob->GetParticle(i));
    size_t id(1<<p_ampl->Legs().size());
    ColorID col(copy->GetFlow(1),copy->GetFlow(2));
    col=col.Conj();
    Flavour flav(copy->Flav().Bar());
    Vec4D mom(-copy->Momentum());
    p_ampl->CreateLeg(mom,flav,col,id);
    Cluster_Leg * leg(p_ampl->Legs().back());
    leg->SetNMax(blob->NOutP()+2);
    leg->SetKTStart(0.);
    leg->SetKTMax(0.);
    leg->SetKTVeto(0.);
    leg->GetSpectators().clear();
    leg->SetConnected(false);
  }
  while (!singlets.empty()) {
    ParticleList * sing=singlets.front();
    while (!sing->empty()) {
      Particle *const copy(sing->front());
      size_t id(1<<p_ampl->Legs().size());
      ColorID col(copy->GetFlow(1),copy->GetFlow(2));
      Flavour flav(copy->Flav());
      Vec4D mom(copy->Momentum());
      p_ampl->CreateLeg(mom,flav,col,id);
      Cluster_Leg * leg(p_ampl->Legs().back());
      leg->SetStat(0);
      leg->SetKTStart(m_minkt2);
      leg->SetKTMax(m_minkt2);
      leg->SetKTVeto(m_minkt2);
      leg->SetNMax(blob->NOutP()+3);
      leg->GetSpectators().clear();
      leg->SetConnected(false);
      sing->pop_front();
    }
    delete sing;
    singlets.pop_front();
  }

  ClusterLeg_Vector legs(p_ampl->Legs());
  Cluster_Leg * split, * spect;

  double kt2max, kt2min, kt2FS, sFS, ysplit, totmax(m_minkt2);
  double shat((legs[0]->Mom()+legs[1]->Mom()).Abs2());
  size_t nlegs(legs.size());

  for (size_t i=2;i<nlegs;i++) {
    split   = legs[i];
    ysplit  = dabs(split->Mom().Y());
    kt2max  = m_minkt2; 
    kt2min  = shat; 
    for (size_t j=nlegs;j>0;j--) {
      if (i==j-1) continue;
      spect = legs[j-1];
      int nconn(ColorConnected(split->Col(),spect->Col()));
      if (nconn==0) continue;
      kt2FS = PT2(split->Mom(),spect->Mom(),
		  i==2 || i==nlegs-1 || j==nlegs-1 || j<=3);
      if (j>2) sFS = (split->Mom()+spect->Mom()).Abs2()/4.; 
          else sFS = 0.;
      if (kt2FS>kt2max) kt2max = kt2FS;
      if (kt2FS<kt2min) kt2min = kt2FS;
      if (j>2) {
	do {
	  split->AddToSpectators(spect);
	  nconn--;
	} while (nconn>0);
	split->SetConnected(true);
      }
    }
    if (!split->Connected()) {
      int j(i);
      while (j==i) { j = 2+int(ATOOLS::ran->Get()*(nlegs-2)); }
      split->AddToSpectators(legs[j]);
      split->SetConnected(false);
      kt2max = split->Mom().PPerp2(); //m_minkt2; 
      kt2min = split->Mom().PPerp2();
    }
    if (kt2max>totmax) totmax = kt2max;
    split->SetKTStart(kt2max);
    split->SetKTVeto(kt2min);
    split->SetKTMax(kt2min);
  }
  p_ampl->SetNIn(blob->NInP());
  p_ampl->SetOrderEW(0);
  p_ampl->SetOrderQCD(blob->NOutP());

  p_ampl->SetMS(this);
  p_ampl->SetKT2(totmax);
  p_ampl->SetMuR2(totmax);
  p_ampl->SetMuF2(totmax);
  p_ampl->SetMu2(totmax);

  return true;
}

double Cluster_Algorithm::Mass(const Flavour &fl) const
{
  return fl.Mass();
}

