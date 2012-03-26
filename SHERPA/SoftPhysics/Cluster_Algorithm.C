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
  p_ampl(NULL), p_clus(NULL), p_jf(NULL),m_showerfac(1.) {}

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
  double pref(4.);
  double pti2(pi.PPerp2()), ptj2(pj.PPerp2());
  double ptij2    = 
    2.*Min(pti2,ptj2)*(cosh(pi.Eta()-pj.Eta())-cos(pi.Phi()-pj.Phi()));
  return m_showerfac*(beam?Min(pref*pti2,ptij2):ptij2);
}

bool Cluster_Algorithm::Cluster(Blob *const blob)
{
  std::list<ParticleList * > singlets;
  ProjectOnSinglets(blob,singlets);

  double scale((m_mode>2?m_minkt2:0.)/1.);
  double ymin(10000.),ymax(-10000.);
  int iymin(-1),iymax(-1),n(1);

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
      n++;
      Particle *const copy(sing->front());
      size_t id(1<<p_ampl->Legs().size());
      ColorID col(copy->GetFlow(1),copy->GetFlow(2));
      Flavour flav(copy->Flav());
      Vec4D mom(copy->Momentum());
      if(mom.Y()<ymin){
	ymin = mom.Y();
	iymin = n;
      }
      if(mom.Y()>ymax){
	ymax = mom.Y();
	iymax = n;
      }
      p_ampl->CreateLeg(mom,flav,col,id);
      Cluster_Leg * leg(p_ampl->Legs().back());
      leg->SetStat(0);
      leg->SetKTStart(scale);
      leg->SetKTMax(scale);
      leg->SetKTVeto(scale);
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
    kt2max  = 0.;
    kt2min  = (m_mode>2)?scale:shat; 
    for (size_t j=nlegs;j>0;j--) {
      if (i==j-1) continue;
      spect = legs[j-1];
      int nconn(ColorConnected(split->Col(),spect->Col()));
      if (nconn==0) continue;
      kt2FS = PT2(split->Mom(),spect->Mom(),!m_resc && (i==iymin||i==iymax));
      if (j>2) sFS = (split->Mom()+spect->Mom()).Abs2(); 
          else sFS = 0.;
      if (kt2FS>kt2max) kt2max = kt2FS;
      //if (sFS>kt2max) kt2max = sFS;
      switch (m_mode) {
      case 7:
      case 6:
      case 5:
      case 4:
	if (kt2FS>kt2min) kt2min = kt2FS;
	break;
      case 3:
      case 2:
      case 1:
      case 0:
      default:
	if (kt2FS<kt2min) kt2min = kt2FS;
	break;
      }
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
      switch (m_mode) {
      case 7:
	kt2max = kt2min = scale;
	break;
      case 6:
	kt2max = kt2min = scale/2.;
	break;
      case 5:
	kt2max = kt2min = sqrt(split->Mom().PPerp2()*scale);
	break;
      case 4:
      case 3:
	kt2max = PT2(split->Mom(),legs[j]->Mom(),true);
	break;
      case 2:
      case 1:
      case 0:
      default:
	kt2max = kt2min = Min(scale/16.,split->Mom().PPerp2());
	break;
      }
    }
    if (kt2max>totmax) totmax = kt2max;
    split->SetKTStart(kt2max);
    split->SetKTMax(kt2max);
    split->SetKTVeto(kt2min);
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

