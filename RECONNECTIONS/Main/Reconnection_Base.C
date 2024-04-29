#include "RECONNECTIONS/Main/Reconnection_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Base::Reconnection_Base() :
  m_analysis(true)
{ }

Reconnection_Base::~Reconnection_Base() {
  if (m_analysis) {
    for (map<string,Histogram *>::iterator hit=m_histomap.begin();
	 hit!=m_histomap.end();hit++) {
      Histogram * histo = hit->second;
      string name  = string("Reconnection_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}

void Reconnection_Base::Initialize() {
  SetParameters(); 
  if (m_analysis) {
    m_histomap[string("Reconn_MassBefore")] = new Histogram(0,0.0,100.0,200);
    m_histomap[string("Reconn_MassAfter")]  = new Histogram(0,0.0,100.0,200);
  }
}

void Reconnection_Base::Reset() {
  for (size_t pos=0;pos<2;pos++) { m_cols[pos].clear(); m_parts[pos].clear(); }
  m_particles.clear();
}

bool Reconnection_Base::HarvestParticles(Blob_List * blobs) {
  // Extract all coloured particles from hadronization blob(s) and fill them into
  // relevant lists.  We assume that there already exists one (or more) active
  // fragmentation blob(s), with incoming coloured particles.
  // TODO: Will have to check how this works with two hadronization blobs, for
  //       example in e+e- -> W+W- -> 4 quarks or, more tricky, pp -> WW -> 4 quarks
  m_found=false;
  Blob * blob;
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    blob = (*bit);
    if (!blob->Has(blob_status::needs_reconnections)) continue;
    m_found = true;
    blob->SetTypeSpec("Colour Reconnections");
    for (int i=0;i<blob->NInP();i++) HarvestParticleInfo(blob->InParticle(i));
    blob->UnsetStatus(blob_status::needs_reconnections |
		      blob_status::needs_hadronization);
  }
  //msg_Out()<<METHOD<<" harvested "
  //	   <<m_cols[0].size()<<" colours and "<<m_cols[1].size()<<" anti-colours.\n";
  return (m_cols[0].size()==m_cols[1].size());
}

void Reconnection_Base::HarvestParticleInfo(ATOOLS::Particle * part) {
  // Only add strong particles with colours different from zero.
  unsigned int col[2];
  for (size_t pos=0;pos<2;pos++) col[pos] = part->GetFlow(pos+1);
  if (col[0]==0 && col[1]==0) return;
  Particle * copy = new Particle(*part);
  colpair cols = colpair(col[0],col[1]);
  // We work with colour pairs <triplet, anti-triplet>.
  // Filling maps of triplet/anti-triplet colour indices to the particles;
  // constructing the singlets will use this, the two sets m_parts[] will be
  // used for the construction of the "distances" between particles.
  for (size_t pos=0;pos<2;pos++) {
    if (col[pos]!=0) {
      m_cols[pos][col[pos]] = copy;
      m_parts[pos].insert(copy);
    }
  }
  m_particles.push_back(copy);
  // Make sure we know the decay blob of the original particle - we will add the
  // copy as ingoing particle to the same blob, but with possibly reshuffled colours;
  // the originl particles will end in a reconnection blob, and their copies will
  // form its outgoing particles.
  copy->SetDecayBlob(part->DecayBlob());
  copy->SetProductionBlob(NULL);
  //msg_Out()<<"* |"<<copy<<"| ("<<copy->Number()<<"): "<<copy->Momentum()<<" @ "<<copy->XProd()
  //	   <<" ["<<copy->ProductionBlob()<<"]\n";
}

double Reconnection_Base::get_total_hadMass(std::map<unsigned int, ATOOLS::Particle*> m_cols[2])
{
	double total_hadMass = 0;
	for(std::map<unsigned int, ATOOLS::Particle*>::iterator cit = m_cols[0].begin(); cit != m_cols[0].end(); cit++)
	{
		double m1 = cit->second->Flav().HadMass();
		double m2 = m_cols[1].find(cit->first)->second->Flav().HadMass();
		total_hadMass += (m1*m2);
	}

	return total_hadMass;
}

void Reconnection_Base::FillMassesInHistogram(std::map<unsigned int, ATOOLS::Particle *> col_pairs[2], std::string hist_name)
{
	std::cout << std::endl << "Filling mass histograms..." << std::endl;

	ATOOLS::Histogram *mass_hist = new ATOOLS::Histogram(0, 0.0, 100.0, 200); // maybe pass min/maxd as arg
	for(std::map<unsigned int, ATOOLS::Particle*>::iterator cit = m_cols[0].begin(); cit != m_cols[0].end(); cit++)
	{
		double m = cit->second->Flav().HadMass() * m_cols[1].find(cit->first)->second->Flav().HadMass();
		mass_hist->Insert(m);
		std::cout << "\n fill value: " << m << std::endl;
	}
	
	m_histomap[hist_name] = mass_hist;

}

void Reconnection_Base::FillLengthInHistogram(float total_length, std::string hist_name)
{
	std::cout << std::endl << "Filling string-length histograms..." << std::endl; // replace eith proper messaging / logger function
	ATOOLS::Histogram *length_hist = new ATOOLS::Histogram(0, 0.0, 100.0, 200); // what should correct binning be?
	// hist with max 20, max 100, max 1,000, then maybe one between say 900 and 1000 (i.e. not starting from zero)
	ATOOLS::Histogram *length_hist_20   = new ATOOLS::Histogram(0,0.0,20.0,10);
	ATOOLS::Histogram *length_hist_100  = new ATOOLS::Histogram(0,0.0,100.0,10);
	ATOOLS::Histogram *length_hist_1000 = new ATOOLS::Histogram(0,0.0,1000.0,10);
	length_hist->Insert(total_length);
	length_hist_20->Insert(total_length);
	length_hist_100->Insert(total_length);
	length_hist_1000->Insert(total_length);
	std::cout << "\n fill value: " << total_length << std::endl;
	m_histomap[hist_name] = length_hist;
	std::string name = hist_name; name.append("_20");
	m_histomap[name] = length_hist_20;
	name = hist_name; name.append("_100");
	m_histomap[name] = length_hist_100;
	name = hist_name; name.append("_1000");
	m_histomap[name] = length_hist_1000;
}

