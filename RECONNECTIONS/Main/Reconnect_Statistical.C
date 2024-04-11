#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnect_Statistical::Reconnect_Statistical() :
  Reconnection_Base(),
  m_Pmode(0), m_Q02(1.), m_etaQ(0.16), m_reshuffle(1./9.), m_kappa(1.) {}


Reconnect_Statistical::~Reconnect_Statistical() {
  m_collist.clear();
}

void Reconnect_Statistical::SetParameters() {
  // Pmode is the mode for the distance measure in momentum space,
  // based on the notion of the string area law, cf. hep-ph/9812423, where the
  // area of a "string" made up of two coloured particles i and j is given by
  // pi*pj-mi*mj (note we assume the gluons to distribute their momentum equally
  // between the two colours).
  // 1 - mode is "power law":   dist = ((pipj-mimj)/norm)^kappa
  // 0 - mode is "logarithmic": dist = log[1+(pipj-mimj)/m_Q02]
  // where norm is the total area of the ordered ensemble
  auto s = Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  string pm   = s["PMODE"].SetDefault("Log").Get<string>(); 
  m_Pmode     = (pm==("Power")?1 : 0); 
  m_Q02       = sqr(s["Q_0"].SetDefault(1.00).Get<double>());
  m_etaQ      = sqr(s["ETA_Q"].SetDefault(0.1).Get<double>());
  m_R02       = sqr(s["R_0"].SetDefault(100.).Get<double>());
  m_etaR      = sqr(s["ETA_R"].SetDefault(0.16).Get<double>());
  m_reshuffle = s["RESHUFFLE"].SetDefault(1./9.).Get<double>();
  m_kappa     = s["KAPPA"].SetDefault(1.).Get<double>();
}

void Reconnect_Statistical::Reset() {
  m_collist.clear();
  Reconnection_Base::Reset();
}

void Reconnect_Statistical::FixPMode(const string & pm) {
  if      (pm=="power" || pm=="Power") m_Pmode=1;
  else if (pm=="log"   || pm=="Log")   m_Pmode=0;
  else {
    msg_Error()<<"Error in "<<METHOD<<"("<<pm<<") is unknown tag.\n"
	       <<"   Will use log-scaling and hope for the best.\n";
    m_Pmode = 0;
  }
}

int Reconnect_Statistical::operator()(Blob_List *const blobs) {
  if (!HarvestParticles(blobs))               return -1;
  if (m_cols[0].empty() && m_cols[1].empty()) return 0;
  m_norm = m_Pmode==1 ? TotalLength() : 0.;
  for (map<unsigned int, Particle *>::iterator cit=m_cols[0].begin();
       cit!=m_cols[0].end();cit++) m_collist.push_back(cit->first);
  size_t N = m_collist.size();
  unsigned int col[2];
  for (size_t i=0;i<sqr(N);i++) {
    if (!SelectColourPair(N,col[0],col[1])) break;
    if (!AttemptSwap(col)) return 0;
  }
  UpdateColours();
  m_collist.clear();
  return 1;
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
  double dist0  = Distance(part[0],part[2]), dist1  = Distance(part[1],part[3]);
  double ndist0 = Distance(part[0],part[3]), ndist1 = Distance(part[1],part[2]);
  double prob   = m_reshuffle*(1.-exp(-m_etaQ*((dist0+dist1)-(ndist0+ndist1))));
  if (ran->Get()<prob) {
    m_cols[1][col[0]] = part[3];
    m_cols[1][col[1]] = part[2];
  }
  return true;
}

void Reconnect_Statistical::UpdateColours() {
  for (size_t i=0;i<2;i++) {
    for (map<unsigned int,Particle *>::iterator cit=m_cols[i].begin();
	 cit!=m_cols[i].end();cit++) {
      cit->second->SetFlow(i+1,cit->first);
    }
  }
}

double Reconnect_Statistical::Distance(Particle * trip,Particle * anti) {
  return MomDistance(trip,anti);
}

double Reconnect_Statistical::MomDistance(Particle * trip,Particle * anti) {
  double p1p2 = trip->Momentum() * anti->Momentum();
  if (trip->Flav().IsGluon()) p1p2 /= 2.;
  if (anti->Flav().IsGluon()) p1p2 /= 2.;
  double m1m2 = trip->Flav().HadMass() * anti->Flav().HadMass();
  return (m_Pmode==1?
	  pow((p1p2-m1m2),m_kappa)/m_norm :
	  log(1+(p1p2-m1m2)/m_Q02) );
}

double Reconnect_Statistical::PosDistance(Particle * trip,Particle * anti) {
  double xdist2 = dabs((trip->XProd()-anti->XProd()).Abs2());
  return xdist2<m_R02 ? 1. : pow(xdist2/m_R02, m_etaR);
}

double Reconnect_Statistical::ColDistance(Particle * trip,Particle * anti) {
  return trip->GetFlow(1)==anti->GetFlow(2) ? 1. : m_reshuffle;
}

double Reconnect_Statistical::TotalLength() {
  double total = 1.;
  Particle * part1, * part2;
  for (map<unsigned int, Particle *>::iterator cit=m_cols[0].begin();
       cit!=m_cols[0].end();cit++) {
    part1  = cit->second;
    part2  = m_cols[1].find(cit->first)->second;
    total += Distance(part1,part2);
  }
  return total/m_parts[0].size();
}

