#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnect_Statistical::Reconnect_Statistical() : Reconnection_Base() {}

Reconnect_Statistical::~Reconnect_Statistical() {}

void Reconnect_Statistical::ReadWeightParameters(Default_Reader *const defaultreader) {
  // Pmode is the mode for the distance measure in momentum space.
  // 0 - mode is "linear":    dist = log(1+sij/Q0^2)
  // 1 - mode is "power law": dist = exp[eta * log(1+sij/Q0^2) ] 
  m_Pmode     = defaultreader->GetValue<double>("RECONNECTIONS::PMODE",0);
  m_Q02       = sqr(defaultreader->GetValue<double>("RECONNECTIONS::Q_0",0.25));
  m_eta       = sqr(defaultreader->GetValue<double>("RECONNECTIONS::eta",0.16));
  m_R02       = sqr(defaultreader->GetValue<double>("RECONNECTIONS::R_0",1.));
  m_reshuffle = 1./defaultreader->GetValue<double>("RECONNECTIONS::RESHUFFLE",1./3.);
  m_kappa     = defaultreader->GetValue<double>("RECONNECTIONS::kappa",2.);
}

void Reconnect_Statistical::Reset() {
  m_collist.clear();
  Reconnection_Base::Reset();
}

bool Reconnect_Statistical::operator()(Blob_List *const blobs) {
  if (!HarvestParticles(blobs)) return false;
  m_norm = TotalLength()/pow(m_parts[0].size(),m_kappa);
  //msg_Out()<<METHOD<<" with "<<m_particles.size()<<" "
  //	   <<"["<<m_parts[0].size()<<"/"<<m_parts[1].size()<<"] particles: "
  //	   <<"norm = "<<m_norm<<"\n";
  for (map<unsigned int, Particle *>::iterator cit=m_cols[0].begin();
       cit!=m_cols[0].end();cit++) m_collist.push_back(cit->first);
  size_t N = m_collist.size();
  //msg_Out()<<"   engaging in "<<sqr(N)<<" reshuffling attempts.\n";
  unsigned int col[2];
  for (size_t i=0;i<sqr(N);i++) {
    if (!SelectColourPair(N,col[0],col[1])) break;
    if (!AttemptSwap(col)) return false;;
  }
  UpdateColours();
  return true;
}

bool Reconnect_Statistical::SelectColourPair(const size_t & N,
					     unsigned int & col1, unsigned int & col2) {
  unsigned int trials=0;
  do {
    col1 = m_collist[int(ran->Get()*N)];
    col2 = m_collist[int(ran->Get()*N)];
    if ((trials++)==100) { col1 = col2 = 0; return false; }
  } while (col1 == col2 ||
	   m_cols[0][col1]==m_cols[1][col2] ||
	   m_cols[0][col2]==m_cols[1][col1]);
  return true;
}

bool Reconnect_Statistical::AttemptSwap(const unsigned int col[2]) {
  if (m_cols[0].find(col[0])==m_cols[0].end() ||
      m_cols[0].find(col[1])==m_cols[0].end() ||
      m_cols[1].find(col[0])==m_cols[1].end() ||
      m_cols[1].find(col[1])==m_cols[1].end()) {
    msg_Error()<<"Error in "<<METHOD<<": ill-defined colours.\n";
    return false;
  }
  Particle * part[4];
  for (size_t i=0;i<2;i++) {
    for (size_t j=0;j<2;j++) part[2*i+j] = m_cols[i][col[j]];
  }
  //msg_Out()<<METHOD<<"("<<col[0]<<", "<<col[1]<<"): "
  //	   <<"["<<part[0]->Number()<<"/"<<part[1]->Number()<<"] & "
  //	   <<"["<<part[2]->Number()<<"/"<<part[3]->Number()<<"]\n";
  double dist0  = Distance(part[0],part[2]), dist1  = Distance(part[1],part[3]);
  double ndist0 = Distance(part[0],part[3]), ndist1 = Distance(part[1],part[2]);
  double prob   = exp(-((ndist0+ndist1)-(dist0+dist1))/m_norm);
  //msg_Out()<<"Swap in "<<col[0]<<" <---> "<<col[1]<<": "
  //	   <<"["<<part[0]->Number()<<"|"<<part[2]->Number()<<"]"
  //	   <<"["<<part[1]->Number()<<"|"<<part[3]->Number()<<"] = "
  //	   <<(dist0+dist1)<<" vs. "
  //	   <<"["<<part[0]->Number()<<"|"<<part[3]->Number()<<"]"
  //	   <<"["<<part[1]->Number()<<"|"<<part[2]->Number()<<"] = "
  //	   <<(ndist0+ndist1)<<": prob = "<<prob<<".\n";
  if (  prob>ran->Get() ) {
    //msg_Out()<<" --> Yay! Swap entries: "
    //	     <<col[0]<<" = "<<part[2]->Number()<<" ("<<m_cols[1][col[0]]->Number()<<") and "
    //	     <<col[1]<<" = "<<part[3]->Number()<<" ("<<m_cols[1][col[1]]->Number()<<")\n";
    m_cols[1][col[0]] = part[3];
    m_cols[1][col[1]] = part[2];
    //msg_Out()<<" --> becomes: "
    //	     <<col[0]<<" = "<<part[3]->Number()<<" ("<<m_cols[1][col[0]]->Number()<<") and "
    //	     <<col[1]<<" = "<<part[2]->Number()<<" ("<<m_cols[1][col[1]]->Number()<<")\n";
  }
  return true;
}

void Reconnect_Statistical::UpdateColours() {
  for (size_t i=0;i<2;i++) {
    for (map<unsigned int,Particle *>::iterator cit=m_cols[i].begin();
	 cit!=m_cols[i].end();cit++) {
      //msg_Out()<<"   ("<<i<<") "<<cit->second->Number()<<": "<<cit->second->GetFlow(i+1)<<" --> ";
      cit->second->SetFlow(i+1,cit->first);
      //msg_Out()<<cit->second->GetFlow(i+1)<<"\n";
    }
  }
}

double Reconnect_Statistical::Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) {
  return (MomDistance(trip,anti) * PosDistance(trip,anti) * ColDistance(trip,anti));
}

double Reconnect_Statistical::MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) {
  double p1p2 = ( ((trip->Flav().IsGluon() ? 0.5 : 1.) * trip->Momentum()+
		   (anti->Flav().IsGluon() ? 0.5 : 1.) * anti->Momentum()).Abs2() -
		  (trip->Momentum().Abs2()+anti->Momentum().Abs2()) );
  return log(1.+p1p2/m_Q02);
}

double Reconnect_Statistical::PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) {
  double xdist2 = dabs((trip->XProd()-anti->XProd()).Abs2());
  return xdist2<m_R02 ? 1. : sqrt(xdist2/m_R02);
}

double Reconnect_Statistical::ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) {
  return trip->GetFlow(1)==anti->GetFlow(2) ? 1. : m_reshuffle;
}

double Reconnect_Statistical::TotalLength() {
  double total = 0.;
  Particle * part1, * part2;
  for (map<unsigned int, Particle *>::iterator cit=m_cols[0].begin();
       cit!=m_cols[0].end();cit++) {
    part1 = cit->second;
    part2 = m_cols[1].find(cit->first)->second;
    total += Distance(part1,part2);
  }
  //msg_Out()<<METHOD<<" = "<<total<<" --------------------------------------\n";
  return total;
}

