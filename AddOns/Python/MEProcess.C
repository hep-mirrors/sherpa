#include "MEProcess.H"

#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "PHASIC++/Channels/Rambo.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <sstream>
#include <algorithm>
#include <string>

MEProcess::MEProcess(SHERPA::Sherpa *a_Generator) :
  m_nlotype(ATOOLS::nlo_type::lo),
  p_amp(ATOOLS::Cluster_Amplitude::New()),
  p_gen(a_Generator), p_proc(NULL), p_rambo(NULL),
  m_ncolinds(0), m_npsp(0), m_nin(0), m_nout(0), p_colint(NULL)
{
  m_kpz[0]=m_kpz[1]=0.;
  ATOOLS::Settings::GetMainSettings()
    .DeclareMatrixSettingsWithEmptyDefault({"MOMENTA"});
}

MEProcess::~MEProcess()
{
  if (p_rambo)         { delete p_rambo; p_rambo=NULL; }
}

std::string MEProcess::Name() const
{
  if(!p_proc) THROW(fatal_error, "Process not initialized");
  return p_proc->Name();
}

void MEProcess::SetMomentumIndices(const std::vector<int> &pdgs)
{
  // fill vector m_mom_inds, such that there is a correspondence
  // p_ampl->Leg(m_mom_inds[i]) <--> pdgs[i]
  DEBUG_FUNC(m_nin<<"->"<<m_nout<<": "<<pdgs);
  if(pdgs.size()<m_nin+m_nout) 
    THROW(fatal_error, "Wrong number of pdg codes given.");
  for (size_t i(0); i<m_nin+m_nout; i++) {
    // find first occurence of flavour 'pdgs[i]' among external legs
    ATOOLS::Flavour flav(abs(pdgs[i]),pdgs[i]<0?true:false);
    bool found = false;
    for (size_t j(0); j<m_nin+m_nout; j++) {
      ATOOLS::Flavour thisflav(j<m_nin?p_amp->Leg(j)->Flav().Bar():p_amp->Leg(j)->Flav());
      if (thisflav==flav) {
	msg_Debugging()<< flav <<" <-> "<< thisflav <<std::endl;
        // if the index j is already assigned, continue
        if (std::find(m_mom_inds.begin(),m_mom_inds.end(),j)!=m_mom_inds.end())
          continue;
        m_mom_inds.push_back(j);
        found=true;
        break;
      }
    }
    if(!found) THROW(fatal_error, "Could not map pdg codes.");
  }
  msg_Debugging()<<m_mom_inds<<std::endl;
}

size_t MEProcess::NumberOfPoints()
{
  if (m_npsp>0) return m_npsp;
  m_npsp=ATOOLS::Settings::GetMainSettings()["MOMENTA"].GetItemsCount();
  return m_npsp;
}

void MEProcess::ReadProcess(size_t n)
{
  DEBUG_FUNC("n="<<n);
  size_t id{ 0 };

  if (NumberOfPoints() == 0)
    THROW(missing_input,"Define momenta using the MOMENTA settings.");

  ATOOLS::Settings& main_settings = ATOOLS::Settings::GetMainSettings();
  main_settings.DeclareMatrixSettingsWithEmptyDefault({ "MOMENTA" });
  ATOOLS::Scoped_Settings s{ main_settings["MOMENTA"].GetItemAtIndex(n - 1) };
  for (const auto& row : s.GetMatrix<std::string>()) {
    msg_Debugging()<<row<<std::endl;
    // either "flav mom" or "flav mom col"
    if (row.size()==2 && row[0]=="End" && row[1]=="point") break;
    if (row.size()!=2 && row.size()!=5 && row.size()!=7) continue;
    if (row[0]=="KP_z_0") {
      msg_Debugging()<<"Set KP-eta values for Beam 0."<<std::endl;
      m_kpz[0]=ATOOLS::ToType<double>(row[1]);
    }
    else if (row[0]=="KP_z_1") {
      msg_Debugging()<<"Set KP-eta values for Beam 1."<<std::endl;
      m_kpz[1]=ATOOLS::ToType<double>(row[1]);
    }
    else if (row[0]=="NLOType") {
      m_nlotype=ATOOLS::ToType<ATOOLS::nlo_type::code>(row[1]);
    }
    else {
      int kf(ATOOLS::ToType<int>(row[0]));
      ATOOLS::Vec4D p(ATOOLS::ToType<double>(row[1]),
                      ATOOLS::ToType<double>(row[2]),
                      ATOOLS::ToType<double>(row[3]),
                      ATOOLS::ToType<double>(row[4]));
      if (kf!=(long int)p_proc->Flavours()[id]){
        std::stringstream err;
        err << "Momenta must be listed flavour-ordered in run card: " << p_proc->Flavours();
        THROW(fatal_error, err.str());
      }
      ATOOLS::ColorID col(0,0);
      if (row.size()==7) col=ATOOLS::ColorID(ATOOLS::ToType<size_t>(row[5]),
                                             ATOOLS::ToType<size_t>(row[6]));
      SetMomentum(id,p);
      SetColor(id,col);
      id++;
    }
  }
  msg_Debugging()<<*p_amp<<std::endl
                 <<"kpz0="<<m_kpz[0]<<", kpz1="<<m_kpz[1]
                 <<", nlo-type="<<m_nlotype<<std::endl;
}

void MEProcess::SetMomenta(const std::vector<double*> &p)
{
  for (unsigned int i(0); i<m_nin; i++)
    p_amp->Leg(m_mom_inds[i])->SetMom(ATOOLS::Vec4D(-p[i][0], -p[i][1],
                                                    -p[i][2], -p[i][3]));
  for (unsigned int i(m_nin); i<p.size(); i++)
    p_amp->Leg(m_mom_inds[i])->SetMom(ATOOLS::Vec4D( p[i][0],  p[i][1],
                                                     p[i][2],  p[i][3]));
}

void MEProcess::SetMomenta(const ATOOLS::Vec4D_Vector &p)
{
  for (size_t i(0);i<m_nin;i++)        p_amp->Leg(m_mom_inds[i])->SetMom(-p[i]);
  for (size_t i(m_nin);i<p.size();i++) p_amp->Leg(m_mom_inds[i])->SetMom( p[i]);
}

ATOOLS::Vec4D_Vector MEProcess::GetMomenta()
{
  ATOOLS::Vec4D_Vector mom;
  for (unsigned int i(0); i<m_nin; i++)
    mom.push_back( -p_amp->Leg(m_mom_inds[i])->Mom());
  for (unsigned int i(m_nin); i<m_nin+m_nout; i++)
    mom.push_back(  p_amp->Leg(m_mom_inds[i])->Mom());
  return mom;
}

void MEProcess::SetMomentum(const size_t &index, const double &e,
                            const double &px, const double &py,
                            const double &pz)
{
  if (index<m_nin)
    p_amp->Leg(m_mom_inds[index])->SetMom(ATOOLS::Vec4D(-e, -px, -py, -pz));
  else
    p_amp->Leg(m_mom_inds[index])->SetMom(ATOOLS::Vec4D(+e, +px, +py, +pz));
}

void MEProcess::SetMomentum(const size_t &index, const ATOOLS::Vec4D &p)
{
  if (index<m_nin) p_amp->Leg(m_mom_inds[index])->SetMom(-p);
  else             p_amp->Leg(m_mom_inds[index])->SetMom(p);
}

void MEProcess::SetColor(const size_t &index, const ATOOLS::ColorID& col)
{
  if (index<m_nin) p_amp->Leg(m_mom_inds[index])->SetCol(col.Conj());
  else             p_amp->Leg(m_mom_inds[index])->SetCol(col);
}

void MEProcess::AddInFlav(const int &id)
{
  msg_Debugging()<<METHOD<<"(): "<<id<<std::endl;
  ATOOLS::Flavour flav(id>0?id:-id, id>0 ? true : false);
  p_amp->CreateLeg(ATOOLS::Vec4D(), flav);
  p_amp->SetNIn(p_amp->NIn()+1);
  // PHASIC::Process_Base::SortFlavours(p_amp);
  m_inpdgs.push_back(id);
  m_flavs.push_back(flav);
  m_nin+=1;
}

void MEProcess::AddOutFlav(const int &id)
{
  msg_Debugging()<<METHOD<<"(): "<<id<<std::endl;
  ATOOLS::Flavour flav(id>0?id:-id, id>0 ? false : true);
  p_amp->CreateLeg(ATOOLS::Vec4D(), flav);
  // PHASIC::Process_Base::SortFlavours(p_amp);
  m_outpdgs.push_back(id);
  m_flavs.push_back(flav);
  m_nout+=1;
}

void MEProcess::AddInFlav(const int &id, const int &col1, const int &col2)
{
  msg_Debugging()<<METHOD<<"(): "<<id<<" ("<<col1<<","<<col2<<")"<<std::endl;
  ATOOLS::Flavour flav(id>0?id:-id, id>0 ? false : true);
  p_amp->CreateLeg(ATOOLS::Vec4D(), flav,
                   ATOOLS::ColorID(col1, col2));
  p_amp->SetNIn(p_amp->NIn()+1);
  // PHASIC::Process_Base::SortFlavours(p_amp);
  m_inpdgs.push_back(id);
  m_flavs.push_back(flav);
  m_nin+=1;
}

void MEProcess::AddOutFlav(const int &id, const int &col1, const int &col2)
{
  msg_Debugging()<<METHOD<<"(): "<<id<<" ("<<col1<<","<<col2<<")"<<std::endl;
  ATOOLS::Flavour flav(id>0?id:-id, id>0 ? false : true);
  p_amp->CreateLeg(ATOOLS::Vec4D(), flav,
                   ATOOLS::ColorID(col1, col2));
  // PHASIC::Process_Base::SortFlavours(p_amp);
  m_outpdgs.push_back(id);
  m_flavs.push_back(flav);
  m_nout+=1;
}

double MEProcess::GenerateColorPoint()
{
  if (p_colint==0) THROW(fatal_error, "No color integrator. Make sure Comix is used.");
  p_colint->GeneratePoint();
  for (size_t i=0; i<p_amp->Legs().size(); ++i)
    p_amp->Leg(i)->SetCol(ATOOLS::ColorID(p_colint->I()[i],p_colint->J()[i]));
  SetColors();
  return p_colint->GlobalWeight();
}

void MEProcess::SetColors()
{ 
  if (p_colint==0) THROW(fatal_error, "No color integrator. Make sure Comix is used.");
  PHASIC::Int_Vector ci(p_amp->Legs().size());
  PHASIC::Int_Vector cj(p_amp->Legs().size());
  for (size_t i=0; i<p_amp->Legs().size(); ++i){
      ci[i] = p_amp->Leg(i)->Col().m_i;
      cj[i] = p_amp->Leg(i)->Col().m_j;
    }
  p_colint->SetI(ci);
  p_colint->SetJ(cj);
}

void MEProcess::PrintProcesses() const
{
  SHERPA::Matrix_Element_Handler* me_handler = p_gen->GetInitHandler()
    ->GetMatrixElementHandler();
  msg_Info()<<"Available processes:"<<std::endl;
  for (size_t i(0); i<me_handler->ProcMaps().size(); i++) {
    for (PHASIC::NLOTypeStringProcessMap_Map::const_iterator
         it=(*me_handler->ProcMaps()[i]).begin();
         it!=(*me_handler->ProcMaps()[i]).end();++it) {
      for (PHASIC::StringProcess_Map::const_iterator
           sit=it->second->begin();sit!=it->second->end();++sit) {
        msg_Info()<<sit->first<<" : "<<sit->second->Name()<<std::endl;
      }
    }
  }
}

PHASIC::Process_Base* MEProcess::FindProcess()
{
  SHERPA::Matrix_Element_Handler* me_handler = p_gen->GetInitHandler()
    ->GetMatrixElementHandler();
  PHASIC::Process_Vector procs = me_handler->AllProcesses();
  if (procs.size()>1) THROW(fatal_error,"More than one process initialised.");
  return procs[0];
}

PHASIC::Process_Base* MEProcess::FindProcess(const ATOOLS::Cluster_Amplitude* amp)
{
  DEBUG_FUNC("");
  SHERPA::Matrix_Element_Handler*
      me_handler(p_gen->GetInitHandler()->GetMatrixElementHandler());
  std::string name = PHASIC::Process_Base::GenerateName(amp);
  msg_Debugging()<<"Looking for "<<name<<std::endl;
  for (size_t i(0); i<me_handler->ProcMaps().size(); i++)
    {
      if(me_handler->ProcMaps()[i]->find(m_nlotype)==
	 me_handler->ProcMaps()[i]->end()) continue;
      PHASIC::StringProcess_Map::const_iterator
        pit(me_handler->ProcMaps()[i]->find(m_nlotype)
            ->second->find(name));
      if (pit == me_handler->ProcMaps()[i]->find(m_nlotype)
	                    ->second->end()) continue;
      return pit->second;
    }
  return NULL;
}

void MEProcess::Initialize()
{
  DEBUG_FUNC((p_proc?p_proc->Name():"no process set yet"));

  if (msg_LevelIsDebugging()) PrintProcesses();

  // Try to find process by amplitude (assumes it has been filled
  // through AddInFlav methods etc)
  p_proc=FindProcess(p_amp);

  if(!p_proc)
    {
      // If not found, assume no flavours have been added and we just want
      // the first (and only) process initialized through the run card. Check
      // first if amplitude is really empty to avoid returning a wrong proc.
      if(p_amp->Legs().size()!=0) THROW(fatal_error, "Requested process not found");
      
      // Retrieve proc initialized through run card and fill amplitude
      p_proc = FindProcess();
      for (size_t i(0);i<p_proc->Flavours().size();++i)
	{
	  if(i<p_proc->NIn()) AddInFlav((long int)p_proc->Flavours()[i]);
	  else                AddOutFlav((long int)p_proc->Flavours()[i]);
	}
    }

  // Check if any of the flavours are 'containters'.
  // In this case we can't assign any color structure
  bool container(false);
  for (size_t i(0);i<p_proc->Flavours().size();++i)
    if (p_proc->Flavours()[i].Size() > 1)
      container=true;
  if(!container) SetUpColorStructure();

  std::vector<int> allpdgs;
  for (std::vector<int>::const_iterator it=m_inpdgs.begin();
       it!=m_inpdgs.end(); it++)  allpdgs.push_back(*it);
  for (std::vector<int>::const_iterator it=m_outpdgs.begin();
       it!=m_outpdgs.end(); it++) allpdgs.push_back(*it);
  SetMomentumIndices(allpdgs);
  
  p_rambo = new PHASIC::Rambo(m_nin,m_nout,
			      &m_flavs.front(),
			      p_proc->Generator());
}

void MEProcess::SetUpColorStructure()
{

  for (unsigned int i = 0; i<p_amp->Legs().size(); i++) {
    if (p_amp->Leg(i)->Flav().Strong()) {
      int scharge = p_amp->Leg(i)->Flav().StrongCharge();
      if (scharge == 8)
        m_gluinds.push_back(i);
      else if (scharge == -3)
        m_quabarinds.push_back(i);
      else if (scharge == 3)
        m_quainds.push_back(i);
      else {
	std::stringstream msg; msg << p_amp->Leg(i)->Flav();
        THROW(fatal_error, "External leg with unknown strong charge detected: "+msg.str());
      }
    }
  }

  m_ncolinds = 2*m_gluinds.size() + m_quabarinds.size() + m_quainds.size();
  if (m_ncolinds%2) THROW(fatal_error, "Odd number of color indices");

  for (int i(0); i<pow(3, m_ncolinds); i++) {
    int k(i);
    int mod(0);
    int r(0), g(0), b(0);
    int rb(0), gb(0), bb(0);
    std::vector<int> combination;
    for (int m(0); m<m_ncolinds/2; m++) {
      mod  = k%3;
      switch(mod) {
      case 0: r+=1;
      case 1: g+=1;
      case 2: b+=1;
      }
      combination.push_back(mod+1);
      k = (k-mod)/3;
    }
    for (int m(m_ncolinds/2); m<m_ncolinds; m++) {
      mod  = k%3;
      switch(mod) {
      case 0: rb+=1;
      case 1: gb+=1;
      case 2: bb+=1;
      }
      combination.push_back(mod+1);
      k = (k-mod)/3;
    }
    if (rb==r&&gb==g&&bb==b) m_colcombinations.push_back(combination);
  }

  if(p_proc->Integrator()->ColorIntegrator()!=NULL)
    p_colint = p_proc->Integrator()->ColorIntegrator();
}

double MEProcess::TestPoint(const double& E){
  ATOOLS::Vec4D_Vector p = p_rambo->GeneratePoint(E);
  SetMomenta(p);
  if(p_colint!=NULL) GenerateColorPoint();
  p_rambo->GenerateWeight(&p[0]);
  return p_rambo->Weight();
}

double MEProcess::MatrixElement()
{
  if(p_colint!=NULL) p_colint->SetWOn(false);
  double res(p_proc->Differential(*p_amp, ATOOLS::Variations_Mode::nominal_only, 1|4));
  if(p_colint!=NULL) p_colint->SetWOn(true);
  // Cancel out initial state swap factor
  // which can be accessed through
  // PHASIC::Process_Base::ISSymFac()
  res *= p_proc->ISSymFac();
  return res;
}

double MEProcess::CSMatrixElement()
{
  if (p_colint==NULL) return MatrixElement();
  GenerateColorPoint();
  double r_csme(0.);
  std::vector<std::vector<int> >::const_iterator it;
  std::vector<int>::const_iterator jt;
  for(it=m_colcombinations.begin(); it!=m_colcombinations.end(); ++it) {
    int ind(0);
    int indbar(m_ncolinds/2);
    for(jt=m_gluinds.begin(); jt!=m_gluinds.end(); ++jt) {
      p_amp->Leg(*jt)->SetCol(ATOOLS::ColorID((*it)[ind], (*it)[indbar]));
      ind+=1;
      indbar+=1;
    }
    for(jt=m_quainds.begin(); jt!=m_quainds.end(); ++jt) {
      p_amp->Leg(*jt)->SetCol(ATOOLS::ColorID((*it)[ind], 0));
      ind+=1;
    }
    for(jt=m_quabarinds.begin(); jt!=m_quabarinds.end(); ++jt) {
      p_amp->Leg(*jt)->SetCol(ATOOLS::ColorID(0,(*it)[indbar] ));
      indbar+=1;
    }
    if(ind!=m_ncolinds/2)  THROW(fatal_error, "Internal Error");
    if(indbar!=m_ncolinds) THROW(fatal_error, "Internal Error");
    SetColors();
    r_csme+=MatrixElement();
  }
  // Cancel out initial state swap factor
  // which can be accessed through
  // PHASIC::Process_Base::ISSymFac()
  r_csme *= p_proc->ISSymFac();
  return r_csme;
}

double MEProcess::GetFlux()
{
  ATOOLS::Vec4D p0(-p_amp->Leg(0)->Mom());
  ATOOLS::Vec4D p1(-p_amp->Leg(1)->Mom());
  return 0.25/sqrt(ATOOLS::sqr(p0*p1)-p0.Abs2()*p1.Abs2());
}

std::string MEProcess::GeneratorName()
{
  std::string loopgen("");
  if (p_proc->Info().m_fi.m_nlotype&ATOOLS::nlo_type::loop)
    loopgen="+"+p_proc->Info().m_loopgenerator;
  return p_proc->Generator()->Name()+loopgen;
}

ATOOLS::Flavour MEProcess::GetFlav(size_t i)
{
  if (i>=m_nin+m_nout) THROW(fatal_error,"Index out of bounds.");
  ATOOLS::Flavour fl=p_amp->Leg(i)->Flav();
  return (i<m_nin?fl.Bar():fl);
}

std::vector<double> MEProcess::NLOSubContributions()
  {
    if (p_proc->IsGroup() && p_proc->Size()>1)
      THROW(not_implemented, "Not implemented for process groups");
    
    PHASIC::Process_Base* proc = 
      p_proc->IsGroup() ? (*p_proc)[0] : p_proc;

    std::vector<double> ret;
    if(proc->GetRSSubevtList())
      for(auto& sub : *(proc->GetRSSubevtList()))
	ret.push_back(sub->m_result);

    return ret;
  }
