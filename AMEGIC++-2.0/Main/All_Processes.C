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

void All_Processes::Add(Process_Base * _proc)
{
  m_procs.push_back(_proc);
}

void  All_Processes::SelectOne() {
  DeSelect();
  if (m_totalxs==0) p_selected = m_procs[int(ran.Get()*m_procs.size())];
  else {
    double disc = m_totalxs * ran.Get(); 
    for (int i=0;i<m_procs.size();i++) {
      disc -= m_procs[i]->Total();
      if (disc<0.) {
	p_selected = m_procs[i];
	p_selected->DeSelect();
	return;
      }
    }
    if (disc>0.) { 
      msg.Error()<<"Error in Process_Group::SelectOne() : "<<"Total xsec = "<<m_totalxs<<std::endl;
      return;
    }
  }
}

void All_Processes::RescaleXSec(double) {
  double sumxs=0., summax=0.;
  for (int i=0;i<m_procs.size();++i) {
    sumxs +=m_procs[i]->Total();
    summax+=m_procs[i]->Max();
  }
  m_totalxs = sumxs;
  m_max     = summax;
}


/*----------------------------------------------------------------------------------
  
  Initialization of the process included in All_Processes

  ----------------------------------------------------------------------------------*/


int All_Processes::InitAllProcesses(Interaction_Model_Base * model,Topology * top,Vec4D *& moms,
				    vector<double> & results,vector<Single_Process *> & links)
{
  bool okay = 1;
  for (int i=0;i<m_procs.size();i++) {
    msg.Debugging()<<"========================================================="<<endl
		   <<"========================================================="<<endl
		   <<"Process_Group::InitAmplitude for "<<m_procs[i]->Name()<<endl;
    if (moms) { delete [] moms; moms = 0; }
    if (!(m_procs[i]->InitAmplitude(model,top,moms,results,links))) okay = 0;
  }

  if (moms) { delete [] moms; moms = 0; }
  if (results.size()==0) {
    msg.Error()<<"Error in All_Processes::InitAllProcesses : "<<endl
	       <<" No amplitude constructed for any process. Nothing to integrate."<<endl;
    return -1;
  }

  for (int i=0;i<results.size();i++) {
    msg.Debugging()<<"========================================================="<<endl
		   <<"========================================================="<<endl
		   <<"All_Processes::SetUpIntegrator for "<<links[i]->Name()<<endl;
    if (!(links[i]->SetUpIntegrator()))       okay = 0;
  }
  if (okay) {
    for (int i=0;i<m_procs.size();i++) {
      msg.Debugging()<<"========================================================="<<endl
		     <<"========================================================="<<endl
		     <<"Process_Group::SetUpIntegrator for "<<m_procs[i]->Name()<<endl;
      moms = 0;
      if (m_procs[i]->Partner()==NULL) {      
	if (!(m_procs[i]->SetUpIntegrator())) okay = 0;
      }
    }
  }

  msg.Tracking()<<"All_Processes::InitAllProcesses ";
  if (okay) msg.Tracking()<<" successful."<<endl;
       else msg.Tracking()<<" failed."<<endl;
  return okay;
}

/*----------------------------------------------------------------------------------
  
  Evaluation of the processes

  ----------------------------------------------------------------------------------*/

bool All_Processes::CalculateTotalXSec(string _resdir)
{
  bool okay = 1;
  for (int i=0;i<m_procs.size();i++) {
    msg.Tracking()<<"All_Processes::CalculateTotalXSec for "<<m_procs[i]->Name()<<endl;
    if (!(m_procs[i]->CalculateTotalXSec(_resdir))) okay = 0;
                                               else m_totalxs += m_procs[i]->Total();
  }
  if (m_totalxs<=0.) okay=0;
  return okay;
}

bool All_Processes::OneEvent() {
  SelectOne();
  return p_selected->OneEvent();
}

bool All_Processes::SameEvent() {
  if (p_selected) 
    return p_selected->SameEvent();
  msg.Error()<<" ERROR: in bool All_Processes::SameEvent() "<<endl;
  return 0;
}
