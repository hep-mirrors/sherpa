#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnect_Statistical::Reconnect_Statistical() : Reconnection_Base() {}

Reconnect_Statistical::~Reconnect_Statistical() {}

void Reconnect_Statistical::SetParameters() {
  // Pmode is the mode for the distance measure in momentum space.
  // 0 - mode is "linear":    dist = log(1+sij/Q0^2)
  // 1 - mode is "power law": dist = exp[eta * log(1+sij/Q0^2) ] 
  auto  s = Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  m_Pmode     = s["RECONNECTIONS::PMODE"].SetDefault(0).Get<int>();
  m_Q02       = sqr(s["RECONNECTIONS::Q_0"].SetDefault(1.00).Get<double>());
  m_etaQ      = sqr(s["RECONNECTIONS::etaQ"].SetDefault(0.16).Get<double>());
  m_R02       = sqr(s["RECONNECTIONS::R_0"].SetDefault(1.00).Get<double>());
  m_etaR      = sqr(s["RECONNECTIONS::etaR"].SetDefault(0.16).Get<double>());
  m_reshuffle = 1./(s["RECONNECTIONS::RESHUFFLE"].SetDefault(1./3.).Get<double>());
  m_kappa     = 1./(s["RECONNECTIONS::KAPPA"].SetDefault(2.).Get<double>());
}

void Reconnect_Statistical::Reset() {
  m_collist.clear();
  Reconnection_Base::Reset();
}

int Reconnect_Statistical::operator()(Blob_List *const blobs) {
  if (!HarvestParticles(blobs))               return -1;
  if (m_cols[0].empty() && m_cols[1].empty()) return 0;
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

bool Reconnect_Statistical::
SelectColourPair(const size_t & N,unsigned int & col1, unsigned int & col2) {
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

double Reconnect_Statistical::Distance(Particle * trip,Particle * anti) {
  return (MomDistance(trip,anti) *
	  PosDistance(trip,anti) *
	  ColDistance(trip,anti));
}

double Reconnect_Statistical::MomDistance(Particle * trip,Particle * anti) {
  double p1p2 = ( ((trip->Flav().IsGluon() ? 0.5 : 1.) * trip->Momentum()+
		   (anti->Flav().IsGluon() ? 0.5 : 1.) * anti->Momentum()).
		  Abs2() -
		  (trip->Momentum().Abs2()+anti->Momentum().Abs2()) );
  return m_Pmode==0 ? log(1.+p1p2/m_Q02) : pow(1.+p1p2/m_Q02,m_etaQ);
}

double Reconnect_Statistical::PosDistance(Particle * trip,Particle * anti) {
  double xdist2 = dabs((trip->XProd()-anti->XProd()).Abs2());
  return xdist2<m_R02 ? 1. : pow(xdist2/m_R02, m_etaR);
}

double Reconnect_Statistical::ColDistance(Particle * trip,Particle * anti) {
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

