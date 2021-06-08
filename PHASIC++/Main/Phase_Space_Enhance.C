#include "PHASIC++/Main/Phase_Space_Enhance.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHASIC;
using namespace ATOOLS;

Phase_Space_Enhance::Phase_Space_Enhance() :
  p_obs(NULL), p_func(NULL), p_histo(NULL), p_histo_current(NULL),
  m_xs(1.), m_factor(1.)
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
  if (p_func==NULL) return 1.0;
  else return (*p_func)(p_moms,p_flavs,m_nflavs);
}

double Phase_Space_Enhance::Factor(Process_Base *const process,const double & totalxs)
{
  if (p_obs==NULL) return 1.0;
  double obs=p_histo?p_histo->Xmin():0.0;
  if (!process->Info().Has(nlo_type::rsub)) obs=(*p_obs)(p_moms,p_flavs,m_nflavs);
  else {
    // fixed-order RS, read out with R kinematics
    if (process->Info().m_nlomode==1) {
      obs=(*p_obs)(p_moms,p_flavs,m_nflavs);
    }
    // MC@NLO H, read out with H kinematics
    else {
      obs=(*p_obs)(p_moms,p_flavs,m_nflavs);
    }
  }
  if (p_histo==NULL) return obs;
  if (obs>=p_histo->Xmax()) obs=p_histo->Xmax()-1e-12;
  if (obs<=p_histo->Xmin()) obs=p_histo->Xmin()+1e-12;
  double dsigma=p_histo->Bin(obs);
  if (dsigma<=0.0) {
    PRINT_INFO("Warning: Tried enhancement with dsigma/dobs("<<obs<<")="<<dsigma<<".");
    dsigma=1.0;
  }
  if (m_xs && totalxs>0.0) return 1.0/dsigma/totalxs;
  else return 1.0/dsigma;
}

void Phase_Space_Enhance::RegisterDefaults() {
  Settings& settings = Settings::GetMainSettings();
  settings["ENHANCE_XS"].SetDefault(0);
}

void Phase_Space_Enhance::SetObservable(const std::string &enhanceobs,
					Process_Base * const process)
{
  if (enhanceobs!="1") {
    if (p_obs)
      THROW(fatal_error, "Overwriting ME enhance observable.");
    std::vector<std::string> parts;
    std::stringstream ss(enhanceobs);
    std::string item;
    while(std::getline(ss, item, '|')) {
      parts.push_back(item);
    }
    if (parts.size()<3 || parts.size()>4)
      THROW(fatal_error,"Wrong syntax in enhance observable.");
    p_obs = Enhance_Observable_Base::Getter_Function::GetObject
      (parts[0],Enhance_Arguments(process,parts[0]));
    if (p_obs==NULL) {
      msg_Error()<<METHOD<<"(): Enhance observable not found. Try 'VAR{..}'.\n";
      THROW(fatal_error,"Invalid enhance observable");
    }
    double enhancemin=ToType<double>(parts[1]);
    double enhancemax=ToType<double>(parts[2]);
    double nbins=parts.size()>3?ToType<size_t>(parts[3]):100;
    
    p_histo = new Histogram(1,enhancemin,enhancemax,nbins,"enhancehisto");
    p_histo->InsertRange(enhancemin, enhancemax, 1.0);
    p_histo->MPISync();
    p_histo->Scale(1.0/p_histo->Integral());
    p_histo_current =
      new Histogram(p_histo->Type(),p_histo->Xmin(),p_histo->Xmax(),p_histo->Nbin(),
		    "enhancehisto_current");
  }
}

void Phase_Space_Enhance::SetFunction(const std::string &enhancefunc,
				      Process_Base * const process)
{
  if (enhancefunc!="1") {
    if (p_func)
      THROW(fatal_error,"Attempting to overwrite enhance function");
    p_func = Enhance_Observable_Base::Getter_Function::GetObject
      (enhancefunc,Enhance_Arguments(process,enhancefunc));
    if (p_func==NULL) {
      msg_Error()<<METHOD<<"(): Enhance function not found. Try 'VAR{..}'.\n";
      THROW(fatal_error,"Invalid enhance function.");
    }
  }
}

void Phase_Space_Enhance::AddPoint(const double xs,Process_Base * process) {
  if (p_histo) {
    if (!process->Info().Has(nlo_type::rsub)) {
      double obs((*p_obs)(p_moms,p_flavs,m_nflavs));
      p_histo_current->Insert(obs,xs/m_factor);
    }
    else {
      // fixed-order RS, fill with RS weight and R kinematics
      if (process->Info().m_nlomode==1) {
        double obs((*p_obs)(p_moms,p_flavs,m_nflavs));
        p_histo_current->Insert(obs,xs/m_factor);
      }
      // MC@NLO H, read out with H kinematics
      else {
        double obs((*p_obs)(p_moms,p_flavs,m_nflavs));
        p_histo_current->Insert(obs,xs/m_factor);
      }
    }
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
