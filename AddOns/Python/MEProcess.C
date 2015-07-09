#include "MEProcess.H"

#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <sstream>
#include <algorithm>

MEProcess::MEProcess(SHERPA::Sherpa *a_Generator) :
  m_name(""), p_amp(ATOOLS::Cluster_Amplitude::New()),
  p_gen(a_Generator), p_proc(NULL), m_ncolinds(0),
  m_npsp(0), m_nin(0), m_nout(0)
{
}

MEProcess::~MEProcess()
{
}

bool MEProcess::HasColorIntegrator()
{
  return (p_proc->Integrator()->ColorIntegrator() != 0);
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
  ATOOLS::Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(ATOOLS::rpa->GetPath());
  reader.SetInputFile(ATOOLS::rpa->gen.Variable("MOMENTA_DATA_FILE"));
  m_npsp=reader.GetValue<size_t>("NUMBER_OF_POINTS",1);
  return m_npsp;
}

void MEProcess::SetMomenta(size_t n)
{
  ATOOLS::Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(ATOOLS::rpa->GetPath());
  reader.SetInputFile(ATOOLS::rpa->gen.Variable("MOMENTA_DATA_FILE"));
  std::vector<std::vector<std::string> > momdata;
  if (!reader.MatrixFromFile(momdata,""))
    THROW(missing_input,"No data in "+ATOOLS::rpa->GetPath()
                        +ATOOLS::rpa->gen.Variable("MOMENTA_DATA_FILE")+"'.");
  size_t begin(0),id(0);
  for (size_t nf(0);nf<momdata.size();++nf) {
    std::vector<std::string> &cur(momdata[nf]);
    if (cur.size()==2 && cur[0]=="Point" &&
        ATOOLS::ToType<int>(cur[1])==n) { begin=nf+1; break; }
  }
  for (size_t nf(begin);nf<momdata.size();++nf) {
    std::vector<std::string> &cur(momdata[nf]);
    // either "flav mom" or "flav mom col"
    if (cur.size()==2 && cur[0]=="End" && cur[1]=="point") break;
    if (cur.size()!=5 && cur.size()!=7) continue;
    int kf(ATOOLS::ToType<int>(cur[0]));
    ATOOLS::Vec4D p(ATOOLS::ToType<double>(cur[1]),
                    ATOOLS::ToType<double>(cur[2]),
                    ATOOLS::ToType<double>(cur[3]),
                    ATOOLS::ToType<double>(cur[4]));
    ATOOLS::ColorID col(0,0);
    if (cur.size()==7) col=ATOOLS::ColorID(ATOOLS::ToType<size_t>(cur[5]),
                                           ATOOLS::ToType<size_t>(cur[6]));
    int kfamp(p_amp->Leg(id)->Flav().Kfcode());
    if (id<m_nin) kfamp=-kfamp;
    if (p_amp->Leg(id)->Flav().IsAnti()) kfamp=-kfamp;
    if (kf!=kfamp) THROW(fatal_error,"Wrong momentum ordering.");
    if (id<m_nin) p_amp->Leg(id)->SetMom(-p);
    else          p_amp->Leg(id)->SetMom(p);
    if (id<m_nin) p_amp->Leg(id)->SetCol(col.Conj());
    else          p_amp->Leg(id)->SetCol(col);
    id++;
  }
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

void MEProcess::AddInFlav(const int &id)
{
  DEBUG_FUNC(id);
  p_amp->CreateLeg(ATOOLS::Vec4D(),
                   ATOOLS::Flavour(id>0?id:-id, id>0 ? true : false));
  p_amp->SetNIn(p_amp->NIn()+1);
  m_inpdgs.push_back(id);
  m_nin+=1;
}

void MEProcess::AddOutFlav(const int &id)
{
  DEBUG_FUNC(id);
  p_amp->CreateLeg(ATOOLS::Vec4D(),
                   ATOOLS::Flavour(id>0?id:-id, id>0 ? false : true));
  m_outpdgs.push_back(id);
  m_nout+=1;
}

void MEProcess::AddInFlav(const int &id, const int &col1, const int &col2)
{
  DEBUG_FUNC(id<<" ("<<col1<<","<<col2<<")");
  p_amp->CreateLeg(ATOOLS::Vec4D(),
                   ATOOLS::Flavour(id>0?id:-id, id>0 ? false : true),
                   ATOOLS::ColorID(col1, col2));
  p_amp->SetNIn(p_amp->NIn()+1);
  m_inpdgs.push_back(id);
  m_nin+=1;
}

void MEProcess::AddOutFlav(const int &id, const int &col1, const int &col2)
{
  DEBUG_FUNC(id<<" ("<<col1<<","<<col2<<")");
  p_amp->CreateLeg(ATOOLS::Vec4D(),
                   ATOOLS::Flavour(id>0?id:-id, id>0 ? false : true),
                   ATOOLS::ColorID(col1, col2));
  m_outpdgs.push_back(id);
  m_nout+=1;
}

double MEProcess::GenerateColorPoint()
{
  SP(PHASIC::Color_Integrator) CI = (p_proc->Integrator()->ColorIntegrator());
  if (CI == 0)
    THROW(fatal_error, "No color integrator. Make sure Comix is used.");
  CI->GeneratePoint();
  for (size_t i=0; i<p_amp->Legs().size(); ++i)
    p_amp->Leg(i)->SetCol(ATOOLS::ColorID(CI->I()[i],CI->J()[i]));
  return CI->GlobalWeight();
}

void MEProcess::SetColors()
{ 
  PHASIC::Int_Vector ci(p_amp->Legs().size());
  PHASIC::Int_Vector cj(p_amp->Legs().size());
  SP(PHASIC::Color_Integrator) CI = (p_proc->Integrator()->ColorIntegrator());
  if (CI==0)
    THROW(fatal_error, "No color integrator. Make sure Comix is used.");
  CI->GeneratePoint();
  for (size_t i=0; i<p_amp->Legs().size(); ++i)
    {
      ci[i] = p_amp->Leg(i)->Col().m_i;
      cj[i] = p_amp->Leg(i)->Col().m_j;
    }
  CI->SetI(ci);
  CI->SetJ(cj);
}

PHASIC::Process_Base* MEProcess::FindProcess()
{
  SHERPA::Matrix_Element_Handler* me_handler = p_gen->GetInitHandler()->GetMatrixElementHandler();
  PHASIC::Process_Base::SortFlavours(p_amp);
  m_name = PHASIC::Process_Base::GenerateName(p_amp);
  for (unsigned int i(0); i<me_handler->ProcMaps().size(); i++)
    for (PHASIC::NLOTypeStringProcessMap_Map::const_iterator sit(me_handler->ProcMaps()[i]->begin());
	 sit!=me_handler->ProcMaps()[i]->end();++sit) {
      PHASIC::StringProcess_Map::const_iterator pit(sit->second->find(m_name));
      if (pit!=sit->second->end()) return pit->second;
    }
  return NULL;
}

void MEProcess::Initialize()
{
  p_proc = FindProcess();
  // if no process was found, assume there is only
  // one initialized in the run card and take that one
  DEBUG_FUNC((p_proc?p_proc->Name():"no process set yet"));
  if(!p_proc){
    SHERPA::Matrix_Element_Handler* me_handler = p_gen->GetInitHandler()
      ->GetMatrixElementHandler();
    PHASIC::Process_Vector procs = me_handler->AllProcesses();
    if (procs.size()>1) THROW(fatal_error,"More than one process initialised.");
    p_proc=procs[0];
    msg_Debugging()<<"Process: "<<p_proc->Name()<<std::endl;
    // fill cluster amplitude according to process
    m_nin=p_proc->NIn();
    m_nout=p_proc->NOut();
    for (size_t i(0);i<p_proc->Flavours().size();++i) {
      ATOOLS::Flavour fl=p_proc->Flavours()[i];
      if (i<m_nin) fl=fl.Bar();
      p_amp->CreateLeg(ATOOLS::Vec4D(),fl);
      m_inpdgs.push_back(fl.IsAnti()?-fl.Kfcode():fl.Kfcode());
    }
    p_amp->SetNIn(m_nin);
  }
  m_name=p_proc->Name();
  for (unsigned int i = 0; i<p_amp->Legs().size(); i++) {
    if (p_amp->Leg(i)->Flav().Strong()) {
      int scharge = p_amp->Leg(i)->Flav().StrongCharge();
      if (scharge == 8)
        m_gluinds.push_back(i);
      else if (scharge == -3)
        m_quabarinds.push_back(i);
      else if (scharge == 3)
        m_quainds.push_back(i);
      else
        THROW(fatal_error, "External leg with unknown strong charge detected.");
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

  std::vector<int> allpdgs;
  for (std::vector<int>::const_iterator it=m_inpdgs.begin();
       it!=m_inpdgs.end(); it++)  allpdgs.push_back(*it);
  for (std::vector<int>::const_iterator it=m_outpdgs.begin();
       it!=m_outpdgs.end(); it++) allpdgs.push_back(*it);
  SetMomentumIndices(allpdgs);
  if (p_proc->Integrator()->ColorIntegrator()!=NULL)
    p_proc->Integrator()->ColorIntegrator()->GeneratePoint();
}

double MEProcess::MatrixElement()
{
  if (!HasColorIntegrator()) return p_proc->Differential(*p_amp);
  SP(PHASIC::Color_Integrator) ci(p_proc->Integrator()->ColorIntegrator());
  ci->SetWOn(false);
  double res(p_proc->Differential(*p_amp));
  ci->SetWOn(true);
  return res;
}

double MEProcess::CSMatrixElement()
{
  if (!HasColorIntegrator()) return p_proc->Differential(*p_amp,1|4);
  SP(PHASIC::Color_Integrator) ci(p_proc->Integrator()->ColorIntegrator());
  ci->SetWOn(false);
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
    r_csme+=p_proc->Differential(*p_amp,1|4);
  }
  ci->SetWOn(true);
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
  if (p_proc->Info().m_fi.m_nloqcdtype&PHASIC::nlo_type::loop ||
      p_proc->Info().m_fi.m_nloewtype&PHASIC::nlo_type::loop)
    loopgen="+"+p_proc->Info().m_loopgenerator;
  return p_proc->Generator()->Name()+loopgen;
}

ATOOLS::Flavour MEProcess::GetFlav(size_t i)
{
  if (i>=m_nin+m_nout) THROW(fatal_error,"Index out of bounds.");
  ATOOLS::Flavour fl=p_amp->Leg(i)->Flav();
  return (i<m_nin?fl.Bar():fl);
}

