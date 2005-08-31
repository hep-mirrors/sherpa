#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "All_Processes.H"
#include "Random.H"
#include "Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;


/*----------------------------------------------------------------------------------
  
  Management of the processes included in All_Processes

  ----------------------------------------------------------------------------------*/

#include "Debugger.H"

void All_Processes::Add(Process_Base * _proc)
{
  _proc->SetParent(this);
  m_procs.push_back(_proc);
}

bool  All_Processes::SelectOne() {
  if (m_weventmode<0) return SelectOneFromList(); 
  DeSelect();
  if (m_totalxs==0) p_selected = m_procs[int(ran.Get()*m_procs.size())];
  else {
    double disc = m_totalxs * ran.Get(); 
    for (size_t i=0;i<m_procs.size();i++) {
      disc -= m_procs[i]->TotalXS();
      if (disc<0.) {
	p_selected = m_procs[i];
	p_selected->DeSelect();
	return true;
      }
    }
    if (disc>0.) { 
      msg.Error()<<"Error in All_Processes::SelectOne() : "<<std::endl
		 <<"   Total xsec = "<<m_totalxs<<", continue and hope for the best. "<<std::endl;
      return false;
    }
  }
  return true;
}

void All_Processes::RescaleXSec(double) {
  double sumxs=0., summax=0.;
  for (size_t i=0;i<m_procs.size();++i) {
    sumxs +=m_procs[i]->TotalXS();
    summax+=m_procs[i]->Max();
  }
  m_totalxs = sumxs;
  m_max     = summax;
}

void All_Processes::SetupEnhance() {
  double sum = 0.;
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->SetupEnhance();
    sum += m_procs[i]->TotalXS();
  }
  SetTotalXS(sum);
}


/*----------------------------------------------------------------------------------
  
  Initialization of the process included in All_Processes

  ----------------------------------------------------------------------------------*/


int All_Processes::InitAllProcesses(Interaction_Model_Base * model,Topology * top,Vec4D *& moms)
{
  vector<Process_Base *> links,errs;
  bool okay     = 1;
  int totalsize = 0;
  int procs     = 0;
  int current_atom = 0;
  for (size_t i=0;i<m_procs.size();i++) {
    msg_Tracking()<<"========================================================="<<endl
		  <<"========================================================="<<endl
		  <<"Process_Group::InitAmplitude for "<<m_procs[i]->Name()<<endl;
    if (moms) { delete [] moms; moms = 0; }
    if (m_atoms) current_atom = links.size();

    switch(m_procs[i]->InitAmplitude(model,top,moms,links,errs,totalsize,procs,current_atom)) {
    case 1:break;
    case 0:okay = 0;break;
    default:
      msg.Error()<<"Error in All_Processes::InitAllProcesses : failed"<<endl;
      return -1;
    }
  }

  if (moms) { delete [] moms; moms = 0; }
  if (totalsize==0) {
    msg.Error()<<"Error in All_Processes::InitAllProcesses : "<<endl
	       <<" No amplitude constructed for any process. Nothing to integrate."<<endl;
    return -1;
  }
  if (okay) {
    for (size_t i=0;i<links.size();i++) {
      msg_Tracking()<<"========================================================="<<endl
		    <<"========================================================="<<endl
		    <<"All_Processes::SetUpIntegrator for "<<links[i]->Name()<<endl;
      if (!(links[i]->SetUpIntegrator()))       okay = 0;
    }
  }
  
  if (okay) {
    for (size_t i=0;i<m_procs.size();i++) {
      msg_Tracking()<<"========================================================="<<endl
		    <<"========================================================="<<endl
		    <<"Process_Group::SetUpIntegrator for "<<m_procs[i]->Name()<<endl;
      if (m_procs[i]->Partner()==NULL) {      
	if (!(m_procs[i]->SetUpIntegrator())) okay = 0;
      }
    }
  }

  msg_Tracking()<<"All_Processes::InitAllProcesses ";
  if (okay) msg_Tracking()<<" successful."<<endl
			  <<"  "<<procs<<" processes using "<<totalsize<<" libraries."<<endl;
       else msg_Tracking()<<" failed."<<endl
			  <<"  "<<totalsize<<" libraries "<<" for "
			  <<procs<<" processes created."<<endl;
  return okay;
}



/*----------------------------------------------------------------------------------
  
  Evaluation of the processes

  ----------------------------------------------------------------------------------*/

bool All_Processes::CalculateTotalXSec(string _resdir)
{
  bool okay = 1;
  for (size_t i=0;i<m_procs.size();i++) {
    msg_Info()<<"All_Processes::CalculateTotalXSec for "<<m_procs[i]->Name()<<endl;
    if (!m_procs[i]->PSHandler()->CreateIntegrators()) return 0;
    if (!(m_procs[i]->CalculateTotalXSec(_resdir))) okay = 0;
    else m_totalxs += m_procs[i]->TotalXS();
  }
  if (m_totalxs==0.) okay=0;

  if (m_procs.size()>0) {
    for (size_t i=0;i<m_procs.size();i++)
    m_procs[i]->AddToDataCollector(i);
  }


  return okay;
}

bool All_Processes::OneEvent(double _mass) {
  SelectOne();
  return dynamic_cast<Process_Base*>(p_selected)->OneEvent(_mass);
}

bool All_Processes::SameEvent() {
  if (p_selected) 
    return p_selected->SameEvent();
  msg.Error()<<"ERROR in All_Processes::SameEvent() : continue and hope for the best. "<<endl;
  return 0;
}
