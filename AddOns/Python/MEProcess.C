#include "MEProcess.H"

#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"

#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <sstream>
#include <algorithm>

MEProcess::MEProcess(SHERPA::Sherpa *a_Generator){
  m_amp = ATOOLS::Cluster_Amplitude::New();
  m_gen = a_Generator;
  m_proc = NULL;
  m_nin = 0;
  m_nout = 0;
  m_ncolinds = 0;
}

MEProcess::~MEProcess(){
}

// TO BE REMOVED IN FINAL VERSION OR TO BE REPLACED BY PROPER METHOD TO DETERMINE WHICH 
// ME GENERATOR IS USED
bool MEProcess::HasColorIntegrator(){
  return  (m_proc->Integrator()->ColorIntegrator() != 0);
}

void MEProcess::SetMomentumIndices(const std::vector<int> &pdgs){
  if(pdgs.size()<m_nin+m_nout) 
    THROW(fatal_error, "Number of pdg codes in vector does not match total number of external particles");
  for (unsigned int i(0); i<m_nin; i++)
    {
      // find the first occurence of a flavour of type pdgs[i] among the external legs
      bool found = false;
      //std::cout << "Looking for pdgid " << (long int)(ATOOLS::Flavour(abs(pdgs[i]), pdgs[i]<0?false:true))<< std::endl;
      for (unsigned int j(0); j < m_nin; j++)
	{
	  //std::cout << "Checking " << (long int)(m_amp->Leg(j)->Flav()) << std::endl;
	  if(((long int)(m_amp->Leg(j)->Flav()))==((long int)(ATOOLS::Flavour(abs(pdgs[i]), pdgs[i]<0?false:true))))
	    {
	      //std::cout << "Found match, checking if already assigned" << std::endl;
	      // if the index j is already assigned, continue searching
	      if(std::find(m_mom_inds.begin(), m_mom_inds.end(), j) !=m_mom_inds.end())
		continue;
	      m_mom_inds.push_back(j);
	      found=true;
	      break;
	    }
	}
      if(!found)
	THROW(fatal_error, "Could not assign all pdg-ids in vector to external legs of cluster amplitude");
    }
  for (unsigned int i(m_nin); i<m_nin+m_nout; i++)
    {
      // find the first occurence of a flavour of type pdgs[i] among the external legs
      bool found = false;
      //std::cout << "Looking for pdgid " << (long int)(ATOOLS::Flavour(abs(pdgs[i]), pdgs[i]<0?true:false))<< std::endl;
      for (unsigned int j(m_nin); j < m_nin+m_nout; j++)
	{
	  //std::cout << "Checking " << (long int)(m_amp->Leg(j)->Flav()) << std::endl;
	  if(((long int)(m_amp->Leg(j)->Flav()))==((long int)(ATOOLS::Flavour(abs(pdgs[i]), pdgs[i]<0?true:false))))
	    {
	      //std::cout << "Found match, checking if already assigned" << std::endl;
	      // if the index j is already assigned, continue searching
	      if(std::find(m_mom_inds.begin(), m_mom_inds.end(), j) !=m_mom_inds.end())
		continue;
	      m_mom_inds.push_back(j);
	      found=true;
	      break;
	    }
	}
      if(!found)
	THROW(fatal_error, "Could not assign all pdg-ids in vector to external legs of cluster amplitude");
    }  
}

void MEProcess::SetMomenta(double **p){
  for (unsigned int i(0); i<m_nin; i++) 
    {
      m_amp->Leg(m_mom_inds[i])->SetMom(ATOOLS::Vec4D(-p[i][0], -p[i][1], -p[i][2], -p[i][3]));
      //double mass =  sqrt(p[i][0]*p[i][0] - p[i][1]*p[i][1] - p[i][2]*p[i][2] - p[i][3]*p[i][3]);
      //std:: cout << "Flavour " << m_amp->Leg(m_mom_inds[i])->Flav() << " Mass " << mass << std::endl;
    }
  for (unsigned int i(m_nin); m_nin+m_nout; i++)
    {
      m_amp->Leg(m_mom_inds[i])->SetMom(ATOOLS::Vec4D(p[i][0], p[i][1], p[i][2], p[i][3]));  
      //double mass =  sqrt(p[i][0]*p[i][0] - p[i][1]*p[i][1] - p[i][2]*p[i][2] - p[i][3]*p[i][3]);
      //std:: cout << "Flavour " << m_amp->Leg(m_mom_inds[i])->Flav() << " Mass " << mass << std::endl;
    }
}

void MEProcess::SetMomentum(int index, double e, double px, double py, double pz){
  if (index<m_nin)
    m_amp->Leg(m_mom_inds[index])->SetMom(ATOOLS::Vec4D(-e, -px, -py, -pz));
  else
    m_amp->Leg(m_mom_inds[index])->SetMom(ATOOLS::Vec4D(+e, +px, +py, +pz));
}

void MEProcess::AddInFlav(const int &id){
  m_amp->CreateLeg(ATOOLS::Vec4D(),   ATOOLS::Flavour(id>0?id:-id, id>0 ? true : false));
  m_amp->SetNIn(m_amp->NIn()+1);
  m_inpdgs.push_back(id);
  m_nin+=1;
}

void MEProcess::AddOutFlav(const int &id){
  m_amp->CreateLeg(ATOOLS::Vec4D(), ATOOLS::Flavour(id>0?id:-id, id>0 ? false : true));
  m_outpdgs.push_back(id);
  m_nout+=1;
}

void MEProcess::AddInFlav(const int &id, const int &col1, const int &col2){
  m_amp->CreateLeg(ATOOLS::Vec4D(), ATOOLS::Flavour(id>0?id:-id, id>0 ? false : true), ATOOLS::ColorID(col1, col2));
  m_amp->SetNIn(m_amp->NIn()+1);
  m_inpdgs.push_back(id);
  m_nin+=1;
}

void MEProcess::AddOutFlav(const int &id, const int &col1, const int &col2){
  m_amp->CreateLeg(ATOOLS::Vec4D(), ATOOLS::Flavour(id>0?id:-id, id>0 ? false : true), ATOOLS::ColorID(col1, col2));
  m_outpdgs.push_back(id);
  m_nout+=1;
}

double MEProcess::GenerateColorPoint(){
  SP(PHASIC::Color_Integrator) CI = (m_proc->Integrator()->ColorIntegrator());
  if (CI == 0)
    THROW(fatal_error, "No color integrator. Make sure Comix is used and that the initialization method has already been called prior to calling this function");
  CI->GeneratePoint();
  for (size_t i=0; i<m_amp->Legs().size(); ++i)
    m_amp->Leg(i)->SetCol(ATOOLS::ColorID(CI->I()[i],CI->J()[i]));
  return CI->GlobalWeight();
}

void MEProcess::SetColors()
{ 
  PHASIC::Int_Vector ci(m_amp->Legs().size());
  PHASIC::Int_Vector cj(m_amp->Legs().size());
  SP(PHASIC::Color_Integrator) CI = (m_proc->Integrator()->ColorIntegrator());
  if (CI==0)
    THROW(fatal_error, "No color integrator. Make sure Comix is used and that the initialization method has already been called prior to calling this function");
  CI->GeneratePoint();
  for (size_t i=0; i<m_amp->Legs().size(); ++i)
    {
      ci[i] = m_amp->Leg(i)->Col().m_i;
      cj[i] = m_amp->Leg(i)->Col().m_j;
    }
  CI->SetI(ci);
  CI->SetJ(cj);
}

void MEProcess::Initialize(){
  SHERPA::Matrix_Element_Handler* me_handler = m_gen->GetInitHandler()->GetMatrixElementHandler();
  PHASIC::Process_Base::SortFlavours(m_amp);
  m_name = PHASIC::Process_Base::GenerateName(m_amp);
  for (unsigned int i(0); i<me_handler->ProcMaps().size(); i++)
    {
      PHASIC::StringProcess_Map::const_iterator pit(me_handler->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->find(m_name));
    
      //FOR DEBUGGING PURPOSES
      //   std::cout << "Initialized Processes: " << std::endl;
      //   for (PHASIC::StringProcess_Map::const_iterator 
      //   	   it(me_handler->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->begin()); 
      //   	 it !=me_handler->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->end();
      //   	 ++it)
      //     {
      //   	std::cout << "Process " << (it->first)<< std::endl;
      //     }

        if(pit == me_handler->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->end())
          continue;
        else{
          m_proc = pit->second;
          break;
        }
    }
  if (!m_proc)
    THROW(fatal_error, std::string("Process ").append(m_name).append(" not found"));
  for(unsigned int i = 0; i<m_amp->Legs().size(); i++)
    {
      if(m_amp->Leg(i)->Flav().Strong()) 
	{
	  int scharge = m_amp->Leg(i)->Flav().StrongCharge();
	  //std::cout << "Strong flavour " << m_amp->Leg(i)->Flav().Kfcode() << " with charge " << m_amp->Leg(i)->Flav().StrongCharge() << std::endl;
	  if (scharge == 8)
	    m_gluinds.push_back(i);
	  else if (scharge == -3)
	    m_quabarinds.push_back(i);
	  else if (scharge == 3)
	    m_quainds.push_back(i);
	  else
	    THROW(fatal_error, "External leg with strong charge other than 3, -3 or 8 detected, dunno what to do.");
	}
    }
  m_ncolinds = 2*m_gluinds.size() + m_quabarinds.size() + m_quainds.size();
  if(m_ncolinds%2)
    THROW(fatal_error, "Odd number of color indices");
  for(int i(0); i<pow(3, m_ncolinds); i++)
    {
      int k(i);
      int mod(0);
      int r(0), g(0), b(0);
      int rb(0), gb(0), bb(0);
      std::vector<int> combination;
      for(int m(0); m<m_ncolinds/2; m++)
	{
	  mod  = k%3;
	  switch(mod)
	    {
	    case 0: r+=1;
	    case 1: g+=1;
	    case 2: b+=1;
	    }
	  combination.push_back(mod+1);
	  k = (k-mod)/3;
	}
      for(int m(m_ncolinds/2); m<m_ncolinds; m++)
	{
	  mod  = k%3;
	  switch(mod)
	    {
	    case 0: rb+=1;
	    case 1: gb+=1;
	    case 2: bb+=1;
	    }
	  combination.push_back(mod+1);
	  k = (k-mod)/3;
	}
      if(rb==r&&gb==g&&bb==b)
	m_colcombinations.push_back(combination);
    }

  // for(int i(0); i<m_colcombinations.size(); i++)
  //   {
  //     std::cout << "One Combination:" <<std::endl;
  //     for (int j(0); j<m_colcombinations[i].size()/2; j++)
  // 	std::cout << m_colcombinations[i][j];
  //     std::cout << std::endl;
  //     for (int j(m_colcombinations[i].size()/2); j<m_colcombinations[i].size(); j++)
  //     	std::cout << m_colcombinations[i][j];
  //     std::cout <<std::endl;
  //   }

  std::vector<int> allpdgs;
  for(std::vector<int>::const_iterator it=m_inpdgs.begin(); it!=m_inpdgs.end(); it++)
    allpdgs.push_back(*it);
  for(std::vector<int>::const_iterator it=m_outpdgs.begin(); it!=m_outpdgs.end(); it++)
    allpdgs.push_back(*it);
  for(std::vector<int>::const_iterator it=allpdgs.begin(); it!=allpdgs.end(); it++)
    std::cout << *it << std::endl;
  SetMomentumIndices(allpdgs);
}

double MEProcess::MatrixElement(){
  return m_proc->Differential(*m_amp);
}

double MEProcess::CSMatrixElement(){
  SP(PHASIC::Color_Integrator) ci(m_proc->Integrator()->ColorIntegrator());
  if(ci==0) // assume Amegic
    return m_proc->Differential(*m_amp);
  //ci->SetWOn(false);
  double r_csme(0.);
  for(std::vector<std::vector<int> >::const_iterator it=m_colcombinations.begin(); it!=m_colcombinations.end(); ++it)
    {
      int ind(0);
      int indbar(m_ncolinds/2);
      for(std::vector<int>::const_iterator jt=m_gluinds.begin(); jt!=m_gluinds.end(); ++jt)
	{
	  m_amp->Leg(*jt)->SetCol(ATOOLS::ColorID((*it)[ind], (*it)[indbar]));
	  ind+=1;
	  indbar+=1;
	}
      for(std::vector<int>::const_iterator jt=m_quainds.begin(); jt!=m_quainds.end(); ++jt)
	{
	  m_amp->Leg(*jt)->SetCol(ATOOLS::ColorID((*it)[ind], 0));
	  ind+=1;
	} 
      if(ind!=m_ncolinds/2)
	THROW(fatal_error, "Fucked shit up");
      for(std::vector<int>::const_iterator jt=m_quabarinds.begin(); jt!=m_quabarinds.end(); ++jt)
	{
	  m_amp->Leg(*jt)->SetCol(ATOOLS::ColorID(0,(*it)[indbar] ));
	  indbar+=1;
	}
      if(indbar!=m_ncolinds)
	THROW(fatal_error, "Fucked shit up");
      SetColors();
      r_csme+=m_proc->Differential(*m_amp);
    }
  //ci->SetWOn(true);
  return r_csme;
}
