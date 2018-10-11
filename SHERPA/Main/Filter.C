#include "SHERPA/Main/Filter.H"
#include "SHERPA/Tools/Output_Base.H"
#include "ATOOLS/Org/Default_Reader.H"


using namespace SHERPA;
using namespace ATOOLS;


Filter::Filter() :
  m_on(false) { }

Filter::~Filter() {
  for (std::map<Flavour, FilterCriterion *>::iterator fit=m_filters.begin();
       fit!=m_filters.end();fit++) delete fit->second;
  m_filters.clear();
}

bool Filter::Init(const std::string & dir,const std::string & file) {
  Default_Reader reader;
  reader.AddIgnore("[");
  reader.AddIgnore("]");
  reader.SetAddCommandLine(false);
  reader.SetInputPath(dir);
  reader.SetInputFile(file);
  reader.SetMatrixType(mtc::transposed);
  std::vector<std::vector<std::string> > entries;
  reader.MatrixFromFile(entries,"");
  for (size_t i=0;i<entries.size();i++) Add(entries[i]);
  //msg_Out()<<METHOD<<" has "<<m_filters.size()<<" individual filters.\n";
  return (m_on = m_filters.size()>0);
}

void Filter::Add(const std::vector<std::string> & entry) {
  if (entry.size()!=7) return;
  FilterCriterion * crit = new FilterCriterion;
  crit->m_flav   = Flavour(stoi(entry[0]));
  crit->m_etamin = stof(entry[1]);
  crit->m_etamax = stof(entry[2]);
  crit->m_pTmin  = stof(entry[3]);
  crit->m_pTmax  = stof(entry[4]);
  crit->m_Nmin   = stoi(entry[5]);
  crit->m_Nmax   = stoi(entry[6]);
  m_filters[crit->m_flav] = crit;
}

bool Filter::operator()(Blob_List * blobs) {
  if (!m_on) return true;
  Reset();
  HarvestActiveParticles(blobs);
  FilterAccepted();
  return Check();
}

bool Filter::Check() {
  for (std::map<Flavour, FilterCriterion *>::iterator flit=m_filters.begin();
       flit!=m_filters.end();flit++) {
    //msg_Out()<<METHOD<<" looks for "<<flit->first;
    std::map<Flavour, int>::iterator accit = m_accepted.find(flit->first);
    if (accit==m_accepted.end()) {
      //msg_Out()<<" --> not found at all.\n";
      return false;
    }
    if (accit->second < flit->second->m_Nmin ||
	accit->second > flit->second->m_Nmax) {
      //msg_Out()<<" --> wrong multiplicity.\n";
      return false;
    }
    //msg_Out()<<" --> ok.\n";
  }
  return true;
}

void Filter::FilterAccepted() {
  for (std::list<Particle *>::iterator pit=m_particles.begin();
       pit!=m_particles.end();pit++) {
    Flavour flav = (*pit)->Flav();
    std::map<Flavour, FilterCriterion *>::iterator flit = m_filters.find(flav); 
    if (flit==m_filters.end()) continue;
    Vec4D   mom  = (*pit)->Momentum();
    double  eta  = mom.Eta(), pT = mom.PPerp();
    FilterCriterion * crit = flit->second;
    if (eta>=crit->m_etamin && eta<=crit->m_etamax &&
	pT>=crit->m_pTmin   && pT<=crit->m_pTmax) {
      std::map<Flavour, int>::iterator accit = m_accepted.find(flav);
      if (accit==m_accepted.end()) m_accepted[flav]=1;
      else accit->second++;
    }
  }
  //for (std::map<Flavour,int>::iterator accit=m_accepted.begin();
  //     accit!=m_accepted.end();accit++) {
  //  msg_Out()<<"  --> found "<<accit->first<<": "<<accit->second<<"\n";
  //}
}

void Filter::HarvestActiveParticles(Blob_List * blobs) {
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    const Particle_Vector & particles = (*bit)->GetOutParticles();
    for (Particle_Vector::const_iterator pit=particles.begin();
	 pit!=particles.end();pit++) {
      if (!(*pit)->DecayBlob()) m_particles.push_back((*pit));
    }
  }
}

void Filter::Reset() {
  m_particles.clear();
  m_accepted.clear();
}
