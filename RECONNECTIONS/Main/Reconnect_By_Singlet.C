#include "RECONNECTIONS/Main/Reconnect_By_Singlet.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnect_By_Singlet::Reconnect_By_Singlet() :
  Reconnection_Base(), m_weights(Reconnection_Weights(this)) {}

Reconnect_By_Singlet::~Reconnect_By_Singlet() {}

int Reconnect_By_Singlet::operator()(Blob_List *const blobs) {
  if (!HarvestParticles(blobs))               return -1;
  if (m_cols[0].empty() && m_cols[1].empty()) return 0;
  MakeSinglets();
  //if (m_analysis) FillMassesInHistogram(m_histomap["Reconn_MassBefore"]);
  m_weights.FillTables();
  ReshuffleSinglets();
  ReconnectSinglets();
  //if (m_analysis) FillMassesInHistogram(m_histomap["Reconn_MassAfter"]);
  FlattenSinglets();
  return true;
}

void Reconnect_By_Singlet::Reset() {
  m_singlets.clear();
  m_weights.Reset();
  Reconnection_Base::Reset();
}

void Reconnect_By_Singlet::SetParameters() {
  // Pmode is the mode for the distance measure in momentum space.
  // 0 - mode is "linear":    dist = log(1+sij/Q0^2)
  // 1 - mode is "power law": dist = exp[eta * log(1+sij/Q0^2) ] 
  auto  s = Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  m_Pmode     = s["PMODE"].SetDefault(0).Get<int>();
  m_Q02       = sqr(s["Q_0"].SetDefault(1.00).Get<double>());
  m_etaQ      = sqr(s["etaQ"].SetDefault(0.16).Get<double>());
  m_R02       = sqr(s["R_0"].SetDefault(1.00).Get<double>());
  m_etaR      = sqr(s["etaR"].SetDefault(0.16).Get<double>());
  m_reshuffle = 1./(s["Reshuffle"].SetDefault(1./3.).Get<double>());
  m_restring  = 1./(s["Restring"].SetDefault(1./3.).Get<double>());
}

void Reconnect_By_Singlet::MakeSinglets() {
  // Start singlets with triplet particles only (quarks or anti-diquarks, with GetFlow(2)==0,
  // in FindStart), and look for fitting anti-triplet colours (in FindNext), until we hit an
  // anti-triplet particle and the singlet is finished.  Once triplet particles are exhausted,
  // repeat with gluons (again, in FindStart), and make sure we book the starting gluon only
  // once.
  Particle * part, * start = FindStart();
  size_t col = start->GetFlow(1);
  Part_List * sing = new Part_List;
  m_singlets.push_back(sing);
  sing->push_back(start);
  while (start) {
    part = FindNext(col);
    if (part!=start) sing->push_back(part);
    col  = part->GetFlow(1);
    if (col==0 || part==start) {
      start = FindStart();
      if (start) {
	col  = start->GetFlow(1);
	sing = new Part_List;
	m_singlets.push_back(sing);
	sing->push_back(start);
      }
    }
    else m_cols[0].erase(col);
  }
}

void Reconnect_By_Singlet::FlattenSinglets() {
  m_particles.clear();
  while (!m_singlets.empty()) {
    Part_List * singlet = m_singlets.front();
    while (!singlet->empty()) {
      m_particles.push_back(singlet->front());
      singlet->pop_front();
    }
    m_singlets.pop_front();
  }
}


Particle * Reconnect_By_Singlet::FindStart() {
  Particle * start(NULL);
  // Look for triplet-only particles.
  for (std::map<unsigned int, ATOOLS::Particle *>::iterator cpit = m_cols[0].begin();
       cpit!=m_cols[0].end();cpit++) {
    if (cpit->second->GetFlow(2)==0) {
      start = cpit->second;
      break;
    }
  }
  // Didn't find a triplet particle, look for any particle with triplet colour
  if (start==NULL && m_cols[0].size()>0) start = m_cols[0].begin()->second;
  if (start!=0) {
    // Make sure the "start" parton is taken out of potential singlet starters.
    m_cols[0].erase(start->GetFlow(1));
    return start;
  }
  // Nothing found, cannot start other singlets.
  return NULL;
}

Particle * Reconnect_By_Singlet::FindNext(const size_t & col) {
  std::map<unsigned int, ATOOLS::Particle *>::iterator cpit = m_cols[1].find(col);
  if (cpit==m_cols[1].end()) {
    msg_Error()<<"Error in "<<METHOD<<" did not find particle with colour[1] = "<<col<<"\n";
    exit(1);
  }
  Particle * part = cpit->second;
  m_cols[1].erase(cpit);
  return part;
}



void Reconnect_By_Singlet::ReshuffleSinglets() {
  // Go through the singlets.  Reshuffling makes sense only if you have more than 4
  // particles in the singlet.  Keep in mind: no splitting of singlets into two as
  // result of reshuffling.  
  for (list<Part_List *>::iterator sit=m_singlets.begin();sit!=m_singlets.end();sit++) {
    if ((*sit)->size()<4) continue;
    bool hit=true;
    while(hit) hit = ReshuffleSinglet((*sit));
  }
}

bool Reconnect_By_Singlet::ReshuffleSinglet(Part_List * singlet) {
  Part_Iterator pit[6], stopit=singlet->end(); ((stopit--)--)--;
  // Logic: check if one particle should be reinserted at another position.
  // Assume particle 4 should be inserted between 0 and 1.  Then we need to compare
  // <01> * <34> * <45> with <04> * <41> * <35>
  // Similarly, if we insert 1 between 4 and 5, we need to compare
  // <01> * <12> * <45> with <02> * <41> * <15>
  pit[0] = singlet->begin(); 
  for (size_t i=1;i<3;i++) { pit[i] = pit[i-1]; pit[i]++; }  
  do {
    if (pit[2]==singlet->end()) return false;
    pit[3] = pit[1];
    for (size_t i=4;i<6;i++) {
      pit[i] = pit[i-1]; pit[i]++;
      if (pit[i]==singlet->end() || pit[i-1]==singlet->end()) break;
    }  
    while (pit[4]!=singlet->end() && pit[5]!=singlet->end()) {
      double dist01345 = (m_weights(*pit[0],*pit[1]) *
			  m_weights(*pit[3],*pit[4]) *
			  m_weights(*pit[4],*pit[5]));
      double dist04135 = (m_weights(*pit[0],*pit[4]) *
			  m_weights(*pit[4],*pit[1]) *
			  m_weights(*pit[3],*pit[5]));
      if (dist04135!=0 && dist04135 < ran->Get() * dist01345) {
	m_weights.SetWeight(*pit[0],*pit[1],0.);
	m_weights.SetWeight(*pit[3],*pit[4],0.);
	m_weights.SetWeight(*pit[4],*pit[5],0.);
	(*pit[4])->SetFlow(2,(*pit[0])->GetFlow(1));
	(*pit[1])->SetFlow(2,(*pit[4])->GetFlow(1));
	(*pit[5])->SetFlow(2,(*pit[3])->GetFlow(1));
	//msg_Out()<<"Shuffle "<<(**pit[4])<<"\n";
	singlet->insert(pit[1],*pit[4]);
	singlet->erase(pit[4]);
	return true;
      }
      double dist01245 = (m_weights(*pit[0],*pit[1]) *
			  m_weights(*pit[1],*pit[2]) *
			  m_weights(*pit[4],*pit[5]));
      double dist02415 = (m_weights(*pit[0],*pit[2]) *
			  m_weights(*pit[4],*pit[1]) *
			  m_weights(*pit[1],*pit[5]));
      if (dist02415!=0 && dist02415 < ran->Get() * dist01245) {
	m_weights.SetWeight(*pit[0],*pit[1],0.);
	m_weights.SetWeight(*pit[1],*pit[2],0.);
	m_weights.SetWeight(*pit[4],*pit[5],0.);
	(*pit[2])->SetFlow(2,(*pit[0])->GetFlow(1));
	(*pit[1])->SetFlow(2,(*pit[4])->GetFlow(1));
	(*pit[5])->SetFlow(2,(*pit[1])->GetFlow(1));
	//msg_Out()<<"Shuffle "<<(**pit[1])<<"\n";
	singlet->insert(pit[5],*pit[1]);
	singlet->erase(pit[1]);
	return true;
      }
      pit[3]++;
      for (size_t i=4;i<6;i++) {pit[i] = pit[i-1]; pit[i]++; }  
    }
    for (size_t i=0;i<3;i++) pit[i]++;
  } while (pit[2]!=singlet->end());
  return false;
}

void Reconnect_By_Singlet::ReconnectSinglets() {
  list<Part_List *>::iterator sit1, sit2;
  Part_Iterator pit11, pit12, pit21, pit22;
  double dist11, dist12, dist21, dist22;
  bool hit;
  if (m_singlets.size()<2) return;
  do {
    hit = false;
    sit1=m_singlets.begin();
    sit2=m_singlets.begin();sit2++;
    while (sit2!=m_singlets.end() && !hit) {
      // Check for pairs of consecutive partons in both singlets
      // (pit11 & pit12 in singlet 1 and pit21 & pit22 in singlet 2)
      // if they should be reconnected cross-wise, i.e., if singlet 1 contains of
      // parton up to and including pit11 and continues with pit22 to the end of the
      // original singlet 2 and vice versa.  The decision is based on distances between
      // pit11/pit12 and pit21/pit22 vs. pit11/pit22 and pit21/pit12.
      pit11  = (*sit1)->begin(); pit12=pit11; pit12++;
      dist11 = m_weights((*pit11),(*pit12));
      while (pit12!=(*sit1)->end() && !hit) {
	pit21=(*sit2)->begin(); pit22=pit21; pit22++;
	while (pit22!=(*sit2)->end() && !hit) {
	  //msg_Out()<<METHOD<<" tests to shuffle particles: "
	  //	   <<(*pit11)->Number()<<" & "<<(*pit12)->Number()<<" vs "
	  //	   <<(*pit21)->Number()<<" & "<<(*pit22)->Number()<<"\n";
	  dist12 = m_weights((*pit11),(*pit22));
	  dist21 = m_weights((*pit21),(*pit12));
	  dist22 = m_weights((*pit21),(*pit22));
	  double rand = ran->Get();
	  if ((dist11*dist22)>(dist21*dist12)*rand) {
	    //msg_Out()<<"   ("<<dist11<<" * "<<dist22<<" = "<<(dist11*dist22)<<") "
	    //	     <<"vs. ("<<dist21<<" * "<<dist12<<" = "<<(dist12*dist21)<<" * "<<rand<<").\n";
	    hit = true;
	    SpliceSinglets((*sit1),(*sit2),pit12,pit22);
	    AftermathOfSlicing((*pit11),(*pit12),(*pit21),(*pit22));
	    break;
	  }
	  if (!hit) { pit21++; pit22++; }
	}
	if (!hit) { pit11++; pit12++; }
      }
      sit2++;
    }
    sit1++;
  } while (hit);
}

void Reconnect_By_Singlet::SpliceSinglets(Part_List * sing1,Part_List * sing2,
					  Part_Iterator & pit1,Part_Iterator & pit2) {
  Part_List help;
  help.splice(help.begin(),*sing1,pit1,sing1->end());
  sing1->splice(sing1->end(),*sing2,pit2,sing2->end());
  for (Part_Iterator pit=help.begin();pit!=help.end();pit++) sing2->push_back(*pit);
  // msg_Out()<<"--------- Singlet 2 with "<<sing2->size()<<" particles.\n";
  // for (Part_Iterator pit=sing2->begin();pit!=sing2->end();pit++) {
  //   msg_Out()<<"  "<<(*pit)->Number()<<" ["<<(*pit)->GetFlow(1)<<", "<<(*pit)->GetFlow(2)<<"]\n";
  // }  
}

void Reconnect_By_Singlet::AftermathOfSlicing(Particle * part11,Particle * part12,
					      Particle * part21,Particle * part22) {
  // After a slicing has happened, we set the weights to zero to ensure that the slicing is not
  // unmade in subsequent sweeps.
  m_weights.SetWeight(part11,part12,0.);
  m_weights.SetWeight(part11,part22,0.);
  m_weights.SetWeight(part21,part22,0.);
  m_weights.SetWeight(part21,part12,0.);
  part22->SetFlow(2,part11->GetFlow(1));
  part12->SetFlow(2,part21->GetFlow(1));
}

double Reconnect_By_Singlet::Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) {
  return (MomDistance(trip,anti) *
	  PosDistance(trip,anti) *
	  ColDistance(trip,anti));
}

double Reconnect_By_Singlet::MomDistance(Particle * part1,Particle * part2) {
  // Here we take a variant of the Lund lambda measure for the distance in momentum space
  double p1p2 = ((part1->Momentum()+part2->Momentum()).Abs2() -
		 (part1->Momentum().Abs2()+part2->Momentum().Abs2()));
  return m_Pmode==0 ? log(1.+p1p2/m_Q02) : pow(1.+p1p2/m_Q02,m_etaQ);
}

double Reconnect_By_Singlet::PosDistance(Particle * part1,Particle * part2) {
  double xdist2 = dabs((part1->XProd().Perp()-part2->XProd().Perp()).Abs2());
  return xdist2<1.e-6? 1. : pow(xdist2/m_R02, m_etaR);
}

double Reconnect_By_Singlet::ColDistance(Particle * part1,Particle * part2) {
  // For colour connected partons there is no extra colour suppression, which would increase
  // the distance in our combined space.
  if ((part1->GetFlow(1)==part2->GetFlow(2) && part1->GetFlow(1)!=0) ||
      (part1->GetFlow(2)==part2->GetFlow(1) && part1->GetFlow(2)!=0)) return 1.;
  // We find out if the two partons are in the smae singlet or not
  Part_List * singlet1(NULL), * singlet2(NULL);
  Part_Iterator pit1, pit2;
  for (list<Part_List *>::iterator sit=m_singlets.begin();sit!=m_singlets.end();sit++) {
    pit1=find((*sit)->begin(),(*sit)->end(),part1);
    if (pit1!=(*sit)->end()) singlet1 = (*sit);
    pit2=find((*sit)->begin(),(*sit)->end(),part2);
    if (pit2!=(*sit)->end()) singlet2 = (*sit);
    if (singlet1!=NULL && singlet2!=0) break;
  }
  // If they are not in the same singlet, the distance in colour space is given by the
  // "flat" re-stringing probability m_restring;
  if (singlet1!=singlet2) return m_restring;
  // If they are in the same singlet, we first figure out how many steps in colour space
  // the two partons are apart from each other.  The reshuffling probability is given by
  // a weight m_reshuffle ^ (steps-1)
  // TODO: Deal with gluon rings here.
  int steps=0;
  pit2=pit1;
  while (pit2!=singlet1->end() && (*pit2)!=part2) {
    pit2++;
    steps++;
  }
  int minsteps=steps;
  steps = 0;
  pit2=pit1;
  while (pit2!=singlet1->end() && (*pit2)!=part2) {
    pit2--;
    steps++;
  }
  if (steps!=0 && steps<minsteps) minsteps = steps;
  return pow(m_reshuffle,Max(2,minsteps-1));
}

