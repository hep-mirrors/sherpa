#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Histogram.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnect_Statistical::Reconnect_Statistical() : Reconnection_Base() 
{
	// define histogram for the total stringlength and add to histogram map
	// (no idea what appropriate maximum / binning is yet)
	m_histomap[string("total_stringlength")] = new Histogram(0, 100.0, 10.0e18, 200.);
	this->m_stringlength_totals.clear();
}

Reconnect_Statistical::~Reconnect_Statistical() {
  //PlotTotalLength(TotalLength());
  m_collist.clear();
  //msg_Info() << std::endl << "Min. total stringlength: " << this->find_min_totalLength(this->m_stringlength_totals) << std::endl;
  //msg_Info() << std::endl << "Max. total stringlength: " << this->find_max_totalLength(this->m_stringlength_totals) << std::endl;
  write_stringLengths_to_file(this->m_stringlength_totals);
}

void Reconnect_Statistical::SetParameters() {
  // Pmode is the mode for the distance measure in momentum space,
  // based on the notion of the string are law, cf. hep-ph/9812423, where the
  // area of a "string" made up of two coloured particles i and j is given by
  // pi*pj-mi*mj (note we assume the gluons to distribute their momentum equally
  // between the two colours).
  // 0 - mode is "linear":    dist = etaQ * (pipj-mimj)/norm
  // 1 - mode is "power law": dist = log[1+etaQ * (pipj-mimj)/norm]
  // where norm is the total area of the ordered ensemble
  auto s = Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  m_Pmode     = s["PMODE"].SetDefault(0).Get<int>();
  m_Q02       = sqr(s["Q_0"].SetDefault(1.00).Get<double>());
  m_etaQ      = sqr(s["etaQ"].SetDefault(0.1).Get<double>());
  m_R02       = sqr(s["R_0"].SetDefault(100.).Get<double>());
  m_etaR      = sqr(s["etaR"].SetDefault(0.16).Get<double>());
  m_reshuffle = s["Reshuffle"].SetDefault(1./9.).Get<double>();
  m_kappa     = s["kappa"].SetDefault(2.).Get<double>();
}

void Reconnect_Statistical::Reset() {
  m_collist.clear();
  Reconnection_Base::Reset();
}

int Reconnect_Statistical::operator()(Blob_List *const blobs) {
  if (!HarvestParticles(blobs))               return -1;
  if (m_cols[0].empty() && m_cols[1].empty()) return 0;
  m_norm = TotalLength();
  PlotTotalLength(TotalLength()); // don't call twice, put in var ...
  for (map<unsigned int, Particle *>::iterator cit=m_cols[0].begin();
       cit!=m_cols[0].end();cit++) { m_collist.push_back(cit->first); msg_Info() << std::endl << "Particle->Number(): " << cit->second->Number() << std::endl; }
  size_t N = m_collist.size();
  unsigned int col[2];
  for (size_t i=0;i<sqr(N);i++) {
    if (!SelectColourPair(N,col[0],col[1])) break;
    if (!AttemptSwap(col)) return false;;
  }
  UpdateColours();
  m_collist.clear();
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
  double dist0  = Distance(part[0],part[2]), dist1  = Distance(part[1],part[3]);
  double ndist0 = Distance(part[0],part[3]), ndist1 = Distance(part[1],part[2]);
  double prob   = m_reshuffle * exp(-m_etaQ*((ndist0+ndist1)-(dist0+dist1)));
  if (prob>ran->Get()) {
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
  return p1p2-m1m2;
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
    total *= Distance(part1,part2);
    msg_Info() << "Distance: " << Distance(part1, part2) << std::endl;
    msg_Info() << "Total" << total << std::endl;
  }
  return total/pow(m_parts[0].size(),m_kappa);
  //const double norm_total = total/pow(m_parts[0].size(),m_kappa);
  //this.m_stringlength_totals.push_back(norm_total);
  //return norm_total;
}


void Reconnect_Statistical::PlotTotalLength(double total_length)
{
	// define and fill histogram of total string length ...
	//std::cout << std::endl << "total string length: " << total_length << std::endl;
	msg_Info() << std::endl << "total string length: " << total_length << std::endl;
	m_histomap[string("total_stringlength")]->Insert(total_length);
	this->m_stringlength_totals.push_back(total_length);
}

/*double Reconnect_Statistical::find_min_totalLength(std::vector<double> total_vec)
{
	double min_value = total_vec[0];
	for(std::vector<double>::iterator it = total_vec.begin()+1; it != total_vec.end(); it++)
		// check that total_vec.begin()+1 adds 1 to index, not to value at index=0
		if(*it < min_value) min_value = *it;
	return min_value;
}

double Reconnect_Statistical::find_max_totalLength(std::vector<double> total_vec)
{
	double max_value = total_vec[0];
	for(std::vector<double>::iterator it = total_vec.begin()+1; it != total_vec.end(); it++)
		if(*it > max_value) max_value = *it;
	return max_value;
}*/

void Reconnect_Statistical::write_stringLengths_to_file(std::vector<double> total_vec)
{
	std::ofstream output_file("./Reconnection_Analysis/stringlength_vec.txt");
	std::ostream_iterator<double> output_iterator(output_file, "\n");
	std::copy(total_vec.begin(), total_vec.end(), output_iterator);
}

