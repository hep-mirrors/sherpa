#include "PHASIC++/Main/Phase_Space_Enhance.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHASIC;
using namespace ATOOLS;

Phase_Space_Enhance::Phase_Space_Enhance()
    : p_obs(nullptr), p_func(nullptr), p_histo(nullptr), p_histo_current(nullptr),
      m_func_min(std::numeric_limits<double>::lowest()),
      m_func_max(std::numeric_limits<double>::max()), m_xs(1.), m_factor(1.),
      p_moms(nullptr), p_flavs(nullptr), m_nflavs(0)
{
  RegisterDefaults();
  Settings& s = Settings::GetMainSettings();
  m_xs = s["ENHANCE_XS"].Get<int>();
}

void Phase_Space_Enhance::Init(Phase_Space_Handler * psh) {
  p_moms   = &psh->Momenta().front();
  p_flavs  = &psh->Flavs().front();
  m_nflavs = psh->Process()->Process()->NIn()+psh->Process()->Process()->NOut();
}

Phase_Space_Enhance::~Phase_Space_Enhance() {
  if (p_obs)           delete p_obs;
  if (p_func)          delete p_func;
  if (p_histo)         delete p_histo;
  if (p_histo_current) delete p_histo_current;
}

double Phase_Space_Enhance::operator()() {
  if (p_func==nullptr) return 1.0;
  else return (*p_func)(p_moms,p_flavs,m_nflavs);
}

double Phase_Space_Enhance::Factor(double totalxs)
{
  m_factor=1.0;
  if (p_obs) {
    double obs=(*p_obs)(p_moms,p_flavs,m_nflavs);
    if (obs>=p_histo->Xmax()) obs=p_histo->Xmax()-1e-30;
    if (obs<=p_histo->Xmin()) obs=p_histo->Xmin()+1e-30;
    double dsigma=p_histo->Bin(obs);
    if (dsigma>0.0) m_factor *= 1.0/dsigma;
  }
  if (p_func) {
    double obs=(*p_func)(p_moms,p_flavs,m_nflavs);
    if (obs<m_func_min) obs=m_func_min;
    if (obs>m_func_max) obs=m_func_max;
    m_factor *= obs;
  }
  if (m_xs && totalxs>0.0) m_factor /= totalxs;
  return m_factor;
}

void Phase_Space_Enhance::RegisterDefaults() {
  Settings& settings = Settings::GetMainSettings();
  settings["ENHANCE_XS"].SetDefault(0);
}

void Phase_Space_Enhance::SetObservable(const std::string &enhanceobs,
					Process_Base * const process)
{
  if (enhanceobs.empty() || enhanceobs=="1")
    return;
  if (p_obs)
    THROW(fatal_error, "Overwriting ME enhance observable.");
  std::vector<std::string> parts;
  std::stringstream ss(enhanceobs);
  std::string item;
  while(std::getline(ss, item, '|'))
    parts.push_back(item);
  if (parts.size()<3)
    THROW(fatal_error,"Wrong syntax in enhance observable.");
  p_obs = Enhance_Observable_Base::Getter_Function::GetObject
    (parts[0],Enhance_Arguments(process,parts[0]));
  if (p_obs==nullptr) {
    msg_Error()<<METHOD<<"(): Enhance observable not found. Try 'VAR{..}'.\n";
    THROW(fatal_error,"Invalid enhance observable");
  }
  double enhancemin=ToType<double>(parts[1]);
  double enhancemax=ToType<double>(parts[2]);
  int nbins=parts.size()>3?ToType<int>(parts[3]):100;

  p_histo = new Histogram(1,enhancemin,enhancemax,nbins,"enhancehisto");
  p_histo->InsertRange(enhancemin, enhancemax, 1.0);
  p_histo->MPISync();
  p_histo->Scale(1.0/p_histo->Integral());
  p_histo_current = new Histogram(p_histo->Type(),
                                  p_histo->Xmin(),
                                  p_histo->Xmax(),
                                  p_histo->Nbin(),
                                  "enhancehisto_current");
}

void Phase_Space_Enhance::SetFunction(const std::string &enhancefunc,
				      Process_Base * const process)
{
  if (enhancefunc.empty() || enhancefunc=="1") return;
  if (p_func) THROW(fatal_error, "Overwriting ME enhance function.");
  std::vector<std::string> parts;
  std::stringstream ss(enhancefunc);
  std::string item;
  while(std::getline(ss, item, '|')) parts.push_back(item);
  if (parts.empty()) THROW(fatal_error,"Wrong syntax in enhance function.");
  p_func = Enhance_Observable_Base::Getter_Function::GetObject
    (parts[0],Enhance_Arguments(process,parts[0]));
  if (p_func==nullptr) {
    msg_Error()<<METHOD<<"(): Enhance function not found. Try 'VAR{...}'.\n";
    THROW(fatal_error,"Invalid enhance function");
  }
  if (parts.size()>2) {
    m_func_min=ToType<double>(parts[1]);
    m_func_max=ToType<double>(parts[2]);
  }
}

void Phase_Space_Enhance::AddPoint(double xs) {
  if (p_histo) {
    double obs((*p_obs)(p_moms,p_flavs,m_nflavs));
    p_histo_current->Insert(obs,xs/m_factor);
  }
}

void Phase_Space_Enhance::Optimize() {
  if (!p_histo) return;
  p_histo_current->MPISync();
  for (int i(0);i<p_histo_current->Nbin()+2;++i)
    p_histo_current->SetBin(i,dabs(p_histo_current->Bin(i)));
  p_histo_current->Scale(1.0/p_histo_current->Integral());
  p_histo->AddGeometric(p_histo_current);
  p_histo->Scale(1.0/p_histo->Integral());
  p_histo_current->Reset();
}

void Phase_Space_Enhance::ReadIn(const std::string &pID) {
  if (!p_histo) return;
  delete p_histo;
  p_histo = new ATOOLS::Histogram(pID+"/MC_Enhance.histo");
  delete p_histo_current;
  p_histo_current =
    new Histogram(p_histo->Type(),p_histo->Xmin(),p_histo->Xmax(),p_histo->Nbin(),
		  "enhancehisto_current");
}
