#include "AHADIC++/Tools/Cluster_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Splitter::Cluster_Splitter() : 
  Splitter_Base(), 
  m_nmax(size_t(hadpars->Get(string("MaxNumberOfPairs")))),
  m_etax(hadpars->Get(string("SplitExponent"))),
  m_etax_lead(hadpars->Get(string("SplitLeadExponent"))),
  m_etay(hadpars->Get(string("SpectExponent"))),
  m_etay_lead(hadpars->Get(string("SpectLeadExponent"))),
  p_trip(NULL), p_anti(NULL)
{
  m_anapath = string("cluster");
}

Cluster_Splitter::~Cluster_Splitter() {}

bool Cluster_Splitter::operator()(Cluster * cluster) {
  Reset();
  if ((cluster->GetTrip()->m_flav.HadMass()+
       cluster->GetAnti()->m_flav.HadMass()+2.*m_mmin)>cluster->Mass()) {
    return false;
  }
  if (!SelectSplitter(cluster->GetTrip(),cluster->GetAnti())) abort();
  DefineTags();
  ConstructTrafos();
  if (ConstructLightC() && ConstructSystem(cluster)) {
    if (m_ana) Analysis();
    Reset();
    if (!cluster->EnsureMomentum() && !EnforceMomentum(cluster)) {
      return false;
    }
    return true;
  }
  UndoTrafos();
  Reset();
  Cluster_List * clist(cluster->GetClusters());
  while (!clist->empty()) {
    delete clist->front();
    clist->pop_front();
  }
  return false;
}

bool Cluster_Splitter::  
SelectSplitter(Proto_Particle * part1,Proto_Particle * part2) {
  Flavour tflav(part1->m_flav), aflav(part2->m_flav);
  bool q1(tflav.IsQuark() || tflav.IsDiQuark());
  bool q2(aflav.IsQuark() || aflav.IsDiQuark());
  if (!q1 || !q2) return false;
  bool hit1(part1->m_info=='L' || part1->m_info=='B');
  bool hit2(part2->m_info=='L' || part2->m_info=='B');
  m_swap = ((hit2 && !hit1) ||
  	    (((hit1 && hit2) || (!hit1 && !hit2)) && ran->Get()<0.5));

  p_split = m_swap?part2:part1;
  p_spect = m_swap?part1:part2;
  return true;
}

bool Cluster_Splitter::ConstructSystem(Cluster * cluster) {
  m_cms = p_split->m_mom + p_spect->m_mom;
  pair<double,double> exponents(FixExponents());
  bool hit;
  double pt2max(m_pt2max);
  if (m_leadsplit) pt2max *= m_pt2max/Max(m_pt2max,m_LC.m_msplit2);
  if (m_leadspect) pt2max *= m_pt2max/Max(m_pt2max,m_LC.m_mspect2);
  m_pt2min = m_pt02 * 
    m_pt02/Max(m_pt02,m_LC.m_msplit2) * 
    m_pt02/Max(m_pt02,m_LC.m_mspect2);
  m_added  = true; 
  //!(p_spect->m_flav==Flavour(kf_b)||
  //	       p_spect->m_flav==Flavour(kf_b).Bar()||
  //	       p_split->m_flav==Flavour(kf_b)||
  //	       p_split->m_flav==Flavour(kf_b).Bar()); //false; //true; 
  for (size_t i=0;i<10;i++) {
    m_pairs = m_isbeam?1:SelectNumberOfPairs(m_nmax);
    for (size_t i=0;i<m_pairs;i++) {
      long int calls(0);
      m_popped.push_back(new PoppedPair);
      do {
	ConstructKinematics(exponents.first,exponents.second);
	hit = (SelectFlavour(m_popped.back()->m_sqq-(m_added?m_pt2min:0.)) &&  
	       AcceptSystem(pt2max));
      } while (!hit && calls++<=1000);
      if (hit) {
	m_sumx += m_popped.back()->m_x; 
	m_sumy += m_popped.back()->m_y;
      }
      if (!hit && calls>1000) Reset();
    }
    if (hit) break;
  }
  if (!hit) {
    return false;
  }
  MakeKinematics();
  MakeClusters(cluster);
  return true;
}

void Cluster_Splitter::
ConstructKinematics(const double & etax,const double & etay) {
  bool   spectHF(p_spect->m_flav==Flavour(kf_b)||
		 p_spect->m_flav==Flavour(kf_b).Bar());
  bool   splitHF(p_split->m_flav==Flavour(kf_b)||
		 p_split->m_flav==Flavour(kf_b).Bar());
  double sqqmin(4.*m_mmin2+(m_added?m_pt2min:0.));
  double msplit2hat(m_LC.m_msplit2/m_LC.m_smandel); 
  double mspect2hat(m_LC.m_mspect2/m_LC.m_smandel); 
  double xarg(sqqmin/m_LC.m_smandel);
  double xmin(xarg);//?sqrt(xarg):xarg));
  double xmax(1.-m_LC.m_msplit2/m_LC.m_smandel-m_sumx);
  double disc(Max(xmax/2.,4.*sqrt(xarg)));
  if (xmax>disc & xmin<disc) xmax = disc;
  double offsetx(m_pt2min/m_LC.m_smandel),offsety;
  double ymin,ymax,x,y,z,sqq,weight;
  if (spectHF && !splitHF) {
    ymin    = xarg;
    ymax    = 1.-m_LC.m_mspect2/m_LC.m_smandel-m_sumx;
    offsety = m_pt2min/m_LC.m_smandel;
  }
  long int calls(0);
  weight=0;
  do {
    if (spectHF && !splitHF) {
      y       = SelectY(ymin,ymax,etay,offsety);
      xmin    = sqqmin/(y*m_LC.m_smandel); 
      xmax    = 1.-msplit2hat-m_sumx;
      if (xmax>disc & xmin<disc) xmax = disc;
      offsetx = offsety/y;
      x       = SelectY(xmin,xmax,etax,offsetx);
    }
    else {
      x       = SelectY(xmin,xmax,etax,offsetx);
      ymin    = sqqmin/(x*m_LC.m_smandel); 
      //if (spectHF && ymin>disc) continue;
      ymax    = 1.-mspect2hat-m_sumy;
      if (ymax>disc & ymin<disc) ymax = disc;
      offsety = offsetx/x;
      y       = SelectY(ymin,ymax,etay,offsety);
    }
    sqq     = x*y*m_LC.m_smandel;
    if (sqq<sqqmin || 
	1.-(m_sumx+x)<m_LC.m_msplit2/m_LC.m_smandel ||
	1.-(m_sumy+y)<m_LC.m_mspect2/m_LC.m_smandel ||
	((1.-m_sumx-x)*(1.-m_sumy-y)*m_LC.m_smandel<
	 sqr(m_LC.m_mspect+m_LC.m_msplit))) {
      weight = 0.;
    }
    else { 
      z      = SelectZ(4.*m_mmin2/sqq,m_leadspect || m_leadsplit);
      weight = exp(-(sqq-sqqmin)/(4.*m_pt02));
    }
    calls++;
  } while (weight<ran->Get() && calls<=1000);
  PoppedPair * pop(m_popped.back());
  if (calls<=1000) {
    //if (p_spect->m_flav==Flavour(kf_b)||p_spect->m_flav==Flavour(kf_b).Bar()||
    //	p_split->m_flav==Flavour(kf_b)||p_split->m_flav==Flavour(kf_b).Bar()) {
    //  msg_Out()<<METHOD<<" for x = "<<x<<" in ["<<xmin<<", "<<xmax<<"] and "
    //	       <<"y = "<<y<<" in ["<<ymin<<", "<<ymax<<"] with offsets "
    //	       <<offsetx<<" & "<<offsety<<".\n";
    //}
    pop->m_x = x; pop->m_y = y; pop->m_z = z; pop->m_sqq = sqq;
  }
  else {
    pop->m_x = 0.; pop->m_y = 0.; pop->m_z = 0.; pop->m_sqq = 0.;
  }
  pop->m_kt2 = 0.; pop->m_flav = Flavour(kf_none); pop->m_mpop2 = 0.;
  for (size_t i(0);i<2;i++) pop->m_outmom[i] = ATOOLS::Vec4D(0.,0.,0.,0.);
}

bool Cluster_Splitter::AcceptSystem(const double & pt2max) {
  PoppedPair * pop(m_popped.back());
  pop->m_kt2 = pop->m_z*(1.-pop->m_z)*pop->m_sqq-pop->m_mpop2;
  if (pop->m_kt2 < 0. || pop->m_kt2 > pt2max) {
    return false;
  }
  return (((*p_as)(pop->m_sqq,false)*(*p_as)(pop->m_kt2,false))/
	  sqr(p_as->MaxValue())) > ran->Get();
}

pair<double,double> Cluster_Splitter::FixExponents() {
  pair<double,double> exponents(m_leadsplit?m_etax_lead:m_etax,
				m_leadspect?m_etay_lead:m_etay);
  if (m_isbeam) 
    exponents.first = exponents.second = m_etax_lead+m_etay_lead+2;
  return exponents;
}

void Cluster_Splitter::MakeKinematics() {
  Vec4D test(0.,0.,0.,0.),test1(0.,0.,0.,0.);
  m_sumx = m_sumy = 0.;
  for (list<PoppedPair *>::iterator pit=m_popped.begin();
       pit!=m_popped.end();pit++) 
    MakePairKinematics((*pit),test,test1);
  MakeSplitterAndSpectatorMoms(test,test1);
}

void Cluster_Splitter::
MakeSplitterAndSpectatorMoms(Vec4D & test,Vec4D & test1) {
  AlphaBeta((1.-m_sumx)*(1.-m_sumy)*m_LC.m_smandel,m_alphanew,m_betanew);
  test += m_splitmom = 
    (1.-m_sumx)*(1.-m_alphanew)*m_LC.m_pA+
    (1.-m_sumy)*m_betanew*m_LC.m_pB;
  test += m_spectmom = 
    (1.-m_sumx)*m_alphanew*m_LC.m_pA+
    (1.-m_sumy)*(1.-m_betanew)*m_LC.m_pB;
  // if (((p_spect->m_flav==Flavour(kf_b)||
  // 	p_spect->m_flav==Flavour(kf_b).Bar())&&
  //      dabs(m_spectmom[3]/p_spect->m_mom[3])<0.75) ||
  //     ((p_split->m_flav==Flavour(kf_b)||
  // 	p_split->m_flav==Flavour(kf_b).Bar())&&
  //      dabs(m_splitmom[3]/p_split->m_mom[3])<0.75)) {
  //   msg_Out()<<"\n\n\n"
  // 	     <<"GOTCHA!!! ================================================\n"
  // 	     <<METHOD<<" for "<<(*m_popped.begin())->m_flav<<", "
  // 	     <<"x = "<<(*m_popped.begin())->m_x<<", "
  // 	     <<"y = "<<(*m_popped.begin())->m_y<<", "
  // 	     <<"z = "<<(*m_popped.begin())->m_z<<", "
  // 	     <<" and kt = "<<sqrt((*m_popped.begin())->m_kt2)<<":\n"
  // 	     <<"spect ("<<p_spect->m_flav<<"): "
  // 	     <<p_spect->m_mom<<" ---> "<<m_spectmom<<";\n"
  // 	     <<"split ("<<p_split->m_flav<<"): "
  // 	     <<p_split->m_mom<<" ---> "<<m_splitmom<<"\n"
  // 	     <<"   from light cone = "<<m_LC.m_pA<<"/"<<m_LC.m_pB<<".\n"
  // 	     <<"==========================================================\n";
  // }
  m_rotat.RotateBack(m_splitmom);
  m_rotat.RotateBack(m_spectmom);
  m_boost.BoostBack(m_splitmom);
  m_boost.BoostBack(m_spectmom);
  test1+=m_splitmom+m_spectmom;
}

void Cluster_Splitter::
MakePairKinematics(PoppedPair * pop,Vec4D & test,Vec4D & test1) {
  double kt(sqrt(pop->m_kt2));
  double phi(2.*M_PI*ran->Get());
  Vec4D  kperp(0.,kt*cos(phi),kt*sin(phi),0.);
  test += pop->m_outmom[0]  = 
    pop->m_x*pop->m_z*m_LC.m_pA + 
    pop->m_y*(1.-pop->m_z)*m_LC.m_pB + kperp;
  test += pop->m_outmom[1]  = 
    pop->m_x*(1.-pop->m_z)*m_LC.m_pA + 
    pop->m_y*pop->m_z*m_LC.m_pB - kperp;
  m_sumx += pop->m_x;
  m_sumy += pop->m_y;
  for (size_t i(0);i<2;i++) {
    m_rotat.RotateBack(pop->m_outmom[i]);
    m_boost.BoostBack(pop->m_outmom[i]);
  }
  test1+=pop->m_outmom[0]+pop->m_outmom[1];
}

void Cluster_Splitter::MakeClusters(Cluster * cluster) {
  SelectPartners();
  MakeSplitterAndSpectatorClusters(cluster);
  MakeOtherClusters(cluster);
}

void Cluster_Splitter::MakeOtherClusters(Cluster * cluster) {
  if (m_popped.size()==1) return;
  if (m_popped.size()==2) {
    if (p_trip && p_anti) {
      Cluster * newcluster(new Cluster(p_trip,p_anti));
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      return;
    }
    else abort();
  }
  Proto_Particle * trip(p_trip),* anti(p_anti), * part;
  size_t winmom;
  bool   lastistrip;
  for (list<PoppedPair *>::iterator pit=m_popped.begin();
       pit!=m_popped.end();pit++) {
    winmom = ran->Get()<0.5?0:1;
    if (ran->Get()<0.5) {
      part = new Proto_Particle((*pit)->m_flav.Bar(),
				(*pit)->m_outmom[winmom],'l');
      Cluster * newcluster(new Cluster(trip,part));
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      trip = new Proto_Particle((*pit)->m_flav,
				(*pit)->m_outmom[1-winmom],'l');
    }
    else {
      part = new Proto_Particle((*pit)->m_flav,
				(*pit)->m_outmom[winmom],'l');
      Cluster * newcluster(new Cluster(part,anti));
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      anti = new Proto_Particle((*pit)->m_flav.Bar(),
				(*pit)->m_outmom[1-winmom],'l');
    }
  }
  Cluster * newcluster(new Cluster(trip,anti));
  newcluster->SetPrev(cluster);
  cluster->push_back(newcluster);
}

void Cluster_Splitter::MakeSplitterAndSpectatorClusters(Cluster * cluster) {
  Proto_Particle * newspect(new Proto_Particle(*p_spect));
  Proto_Particle * newsplit(new Proto_Particle(*p_split));
  newspect->m_mom = m_spectmom;
  newsplit->m_mom = m_splitmom; 
  Proto_Particle * splitp, * spectp;
  Cluster * newcluster;
  p_trip = p_anti = NULL;
  char info((cluster->GetTrip()->m_info=='B' || 
	     cluster->GetAnti()->m_info=='B') ? 'B':'l');
  if (!m_swap) {
    // split = trip, spect = anti
    splitp = new Proto_Particle((*m_popsplit)->m_flav.Bar(),
				(*m_popsplit)->m_outmom[m_popspliti],info);
    spectp = new Proto_Particle((*m_popspect)->m_flav,
				(*m_popspect)->m_outmom[m_popspecti],info);
    if (!Rearrange()) {
      newcluster = new Cluster(newsplit,splitp);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      newcluster = new Cluster(spectp,newspect);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
    }
    else {
      newcluster = new Cluster(newsplit,newspect);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      newcluster = new Cluster(spectp,splitp);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
    }
    if (m_popsplit!=m_popspect) {
      p_trip = new Proto_Particle((*m_popsplit)->m_flav,
				  (*m_popsplit)->m_outmom[1-m_popspliti],info);
      p_anti = new Proto_Particle((*m_popspect)->m_flav.Bar(),
				  (*m_popspect)->m_outmom[1-m_popspecti],info);
    }
  }
  else {
    // split = anti, spect = trip
    splitp = new Proto_Particle((*m_popsplit)->m_flav,
				(*m_popsplit)->m_outmom[m_popspliti],info);
    spectp = new Proto_Particle((*m_popspect)->m_flav.Bar(),
				(*m_popspect)->m_outmom[m_popspecti],info);
    if (!Rearrange()) {
      newcluster = new Cluster(splitp,newsplit);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      newcluster = new Cluster(newspect,spectp);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
    }
    else {
      newcluster = new Cluster(splitp,spectp);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
      newcluster = new Cluster(newspect,newsplit);
      newcluster->SetPrev(cluster);
      cluster->push_back(newcluster);
    }
    if (m_popsplit!=m_popspect) {
      p_trip = new Proto_Particle((*m_popspect)->m_flav,
				  (*m_popspect)->m_outmom[1-m_popspecti],info);
      p_anti = new Proto_Particle((*m_popsplit)->m_flav.Bar(),
				  (*m_popsplit)->m_outmom[1-m_popspliti],info);
    }
  }
}  

bool Cluster_Splitter::Rearrange() {
  if (m_pairs!=1 || m_isbeam) return false;
  PoppedPair * pop(*m_popped.begin());
  double y1(log(((1.-m_sumx)*(1.-m_alphanew))/((1.-m_sumy)*m_betanew)));
  double y2(log(((1.-m_sumx)*m_alphanew)/((1.-m_sumy)*(1.-m_betanew))));
  double y3(log((m_sumx*pop->m_z)/(m_sumy*(1.-pop->m_z))));
  double y4(log((m_sumx*(1.-pop->m_z))/(m_sumy*pop->m_z)));
  return (0.5*dabs((y1-y2)*(y3-y4))/
	  (dabs((y1-y3)*(y2-y4))+dabs((y1-y4)*(y2-y3))))>ran->Get();
}

void Cluster_Splitter::SelectPartners() {
  double masssplit(1.e12), massspect(1.e12), testsplit, testspect;
  m_popsplit  = m_popspect  = m_popped.end();
  m_popspliti = m_popspecti = 2;
  for (list<PoppedPair *>::iterator pit=m_popped.begin();
       pit!=m_popped.end();pit++) {
    for (size_t i(0);i<2;i++) {
      testsplit = (m_splitmom+(*pit)->m_outmom[i]).Abs2();
      if (testsplit<masssplit && 
	  ((m_popped.size()==1 && i!=m_popspecti) ||
	   pit!=m_popspect)) {
	masssplit   = testsplit;
	m_popsplit  = pit;
	m_popspliti = i;
	break;
      }
    }
  }
  for (list<PoppedPair *>::iterator pit=m_popped.begin();
       pit!=m_popped.end();pit++) {
    for (size_t i(0);i<2;i++) {
      testspect = (m_spectmom+(*pit)->m_outmom[i]).Abs2();
      if (testspect<massspect && 
	  ((m_popped.size()==1 && i!=m_popspliti) ||
	   pit!=m_popsplit)) {
	massspect   = testspect;
	m_popspect  = pit;
	m_popspecti = i;
	break;
      }
    }
  }
}


bool Cluster_Splitter::PoppedMassPossible(const double & m2) {
  PoppedPair * pop(m_popped.back());
  pop->m_kt2 = pop->m_z*(1.-pop->m_z)*pop->m_sqq-m2;
  if (pop->m_kt2<0.) {
    //if (m2<=1.001*m_mmin2)
    //msg_Out()<<"      --> "<<METHOD<<"(sqq = "<<pop->m_sqq<<", "
    //	       <<"z = "<<pop->m_z<<" --> "<<pop->m_kt2<<").\n";
    return false;
  }
  double sumx = m_sumx+pop->m_x;
  double sumy = m_sumy+pop->m_y;
  double snew((1.-sumx)*(1.-sumy)*m_LC.m_smandel);
  if (snew<sqr(m_LC.m_mspect+m_LC.m_msplit) ||
      sumx>1. || sumy>1.||
      !AlphaBeta(snew,m_alphanew,m_betanew)) {
    //msg_Out()<<"      --> "<<METHOD<<" testing LC's failed.\n";
    return false;
  }
  Vec4D tsplit((1.-sumx)*(1.-m_alphanew)*m_LC.m_pA+
	       (1.-sumy)*m_betanew*m_LC.m_pB);
  Vec4D tspect((1.-sumx)*m_alphanew*m_LC.m_pA+
	       (1.-sumy)*(1.-m_betanew)*m_LC.m_pB);
  if (dabs(tsplit.Abs2()/m_LC.m_msplit2-1.)>1.e-6 ||
      dabs(tspect.Abs2()/m_LC.m_mspect2-1.)>1.e-6) {
    //msg_Out()<<"      --> "<<METHOD<<" failed: "
    //	     <<dabs(tsplit.Abs2())<<"/"<<m_LC.m_msplit2
    //	     <<" -> "<<dabs(tsplit.Abs2()/m_LC.m_msplit2-1.)<<" and "
    //	     <<dabs(tspect.Abs2())<<"/"<<m_LC.m_mspect2
    //	     <<" -> "<<dabs(tspect.Abs2()/m_LC.m_mspect2-1.)<<".\n";
    return false;
  }
  return true;
}

size_t Cluster_Splitter::SelectNumberOfPairs(const size_t & nmax) {
  double deltaY(0.5*(dabs(log((1.-m_LC.m_beta)/m_LC.m_alpha)-
			  log(m_LC.m_beta/(1.-m_LC.m_alpha)))));
  size_t number(0);
  while (number<1 || number>nmax) number = ran->Poissonian(deltaY/(2.*M_PI));
  return number;
}

bool Cluster_Splitter::EnforceMomentum(Cluster * cluster) {
  if (cluster->GetClusters()->empty()) abort();
  Vec4D summom(0.,0.,0.,0.);
  for (Cluster_Iterator cit(cluster->GetClusters()->begin());
       cit!=cluster->GetClusters()->end();cit++) {
    summom += (*cit)->Momentum();
  }
  //msg_Out()<<"      --> "<<METHOD<<" for \n"
  //	   <<summom<<" vs. "<<cluster->Momentum()<<" --> "
  //	   <<(summom-cluster->Momentum())<<".\n";
  //for (Cluster_Iterator cit(cluster->GetClusters()->begin());
  //   cit!=cluster->GetClusters()->end();cit++) {
  //msg_Out()<<(**cit);
  //}
  Poincare rest(summom);
  Poincare back(cluster->Momentum());
  for (Cluster_Iterator cit(cluster->GetClusters()->begin());
       cit!=cluster->GetClusters()->end();cit++) {
    (*cit)->Boost(rest);
    (*cit)->BoostBack(back);
  }
  if (!cluster->EnsureMomentum()) abort();
  return true;
}
