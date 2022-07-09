#include "SHERPA/Main/Filter.H"
#include "SHERPA/Tools/Output_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"


using namespace SHERPA;
using namespace ATOOLS;


Filter::Filter() :
  m_on(false) { }

Filter::~Filter() {
  for (std::map<Flavour, FilterCriterion *>::iterator fit=m_filters.begin();
       fit!=m_filters.end();fit++) delete fit->second;
  m_filters.clear();
}

bool Filter::Init() {
  for (auto s : Settings::GetMainSettings()["FILTERS"].GetItems()) {
    Add(s);
  }
  return (m_on = m_filters.size()>0);
}

void Filter::Add(Scoped_Settings& s) {
  if (!s["Kf"].IsSetExplicitly())
    THROW(fatal_error, "Every filter needs to specify \"Kf: <kf_code>\".");
  FilterCriterion * crit = new FilterCriterion;
  crit->m_flav   = Flavour(s["Kf"].SetDefault(0).Get<int>());
  crit->m_etamin = s["EtaMin"].SetDefault(-1e20).Get<double>();
  crit->m_etamax = s["EtaMax"].SetDefault( 1e20).Get<double>();
  crit->m_pTmin  = s["PTMin"] .SetDefault(    0).Get<double>();
  crit->m_pTmax  = s["PTMax"] .SetDefault( 1e20).Get<double>();
  crit->m_Nmin   = s["NMin"]  .SetDefault(    0).Get<int>();
  crit->m_Nmax   = s["NMax"]  .SetDefault( 1e20).Get<int>();
  m_filters[crit->m_flav] = crit;
}

void Filter::ShowSyntax(int mode)
{
  if (!msg_LevelIsInfo() || mode == 0) return;
  msg_Out() << METHOD << "(): {\n\n";
  msg_Out() << "  FILTERS:\n";
  msg_Out() << "  - Kf: <kf_code>\n";
  msg_Out() << "    # optional:\n";
  msg_Out() << "    EtaMin: <eta_min>\n";
  msg_Out() << "    EtaMax: <eta_max>\n";
  msg_Out() << "    PTMin:  <pt_min>\n";
  msg_Out() << "    PTMax:  <pt_max>\n";
  msg_Out() << "    NMin:  <pt_min>\n";
  msg_Out() << "    NMax:  <pt_max>\n";
  msg_Out() << "  - ...\n";
  msg_Out() << "\n}" << std::endl;
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
