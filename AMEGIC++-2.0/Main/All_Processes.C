#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "All_Processes.H"
#include "Random.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace std;

/*----------------------------------------------------------------------------------
  
  Management of the processes included in All_Processes

  ----------------------------------------------------------------------------------*/

void All_Processes::Add(Process_Base * _proc)
{
  msg.Tracking()<<"Add process group "<<_proc->Name()<<" to All_Processes."<<endl;
  procs.push_back(_proc);
}

void  All_Processes::SelectOne() {
  //  msg.Out()<<"AllProcesses::SelectOne : totalxs, max = "<<totalxs<<", "<<max<<endl;
  DeSelect();
  if (totalxs==0) selected = procs[int(ran.Get()*procs.size())];
  else {
    double disc = totalxs * ran.Get(); 
    //    cout<<" disc="<<disc<<endl;
    for (int i=0;i<procs.size();i++) {
      disc -= procs[i]->Total();
      if (disc<0.) {
	selected = procs[i];
	selected->DeSelect();
	return;
      }
      //    cout<<" disc="<<disc<<endl;
    }
    //    cout<<" selected="<<selected->Name()<<endl;
    if (disc>0.) { 
      msg.Error()<<"Error in Process_Group::SelectOne() : ";
      msg.Error()<<"Total xsec = "<<totalxs<<endl;
      return;
    }
  }
}

void All_Processes::RescaleXSec(double) {
  double sumxs=0., summax=0.;
  for (int i=0;i<procs.size();++i) {
    sumxs +=procs[i]->Total();
    summax+=procs[i]->Max();
  }

  cout<<" in RescaleXSec of "<<name<<endl;
  cout<<"  "<<rfactor<<" "<<totalxs<<" "<<max<<" "<<endl;
  totalxs=sumxs;
  max    =summax;
  cout<<"  "<<rfactor<<" "<<totalxs<<" "<<max<<" "<<endl;

}


/*----------------------------------------------------------------------------------
  
  Initialization of the process included in All_Processes

  ----------------------------------------------------------------------------------*/


int All_Processes::InitAllProcesses(Topology* top,Vec4D *& moms,
				    vector<double> & results,vector<Single_Process *> & links)
{
  bool okay = 1;
  for (int i=0;i<procs.size();i++) {
    msg.Tracking()<<"========================================================="<<endl;
    msg.Tracking()<<"========================================================="<<endl;
    msg.Tracking()<<"Process_Group::InitAmplitude for "<<procs[i]->Name()<<endl;
    if (moms) { delete [] moms; moms = 0; }
    if (!(procs[i]->InitAmplitude(top,moms,results,links))) okay = 0;
  }

  if (moms) { delete [] moms; moms = 0; }
  if (results.size() == 0) {
    msg.Error()<<"Error in All_Processes::InitAllProcesses : "<<endl
	       <<" No amplitude constructed for any process. Nothing to integrate."<<endl;
    return -1;
  }

  msg.Tracking()<<"Set up "<<results.size()<<" integrators : "<<endl;
  for (int i=0;i<results.size();i++) {
    msg.Tracking()<<"========================================================="<<endl;
    msg.Tracking()<<"========================================================="<<endl;
    msg.Tracking()<<"All_Processes::SetUpIntegrator for "<<links[i]->Name()<<endl;
    if (!(links[i]->SetUpIntegrator()))       okay = 0;
  }
  if (okay) {
    for (int i=0;i<procs.size();i++) {
      msg.Tracking()<<"========================================================="<<endl;
      msg.Tracking()<<"========================================================="<<endl;
      msg.Tracking()<<"Process_Group::SetUpIntegrator for "<<procs[i]->Name()<<endl;
      moms = 0;
      if (!(procs[i]->SetUpIntegrator())) okay = 0;
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

bool All_Processes::CalculateTotalXSec()
{
  bool okay = 1;
  for (int i=0;i<procs.size();i++) {
    msg.Tracking()<<"All_Processes::CalculateTotalXSec for "<<procs[i]->Name()<<endl;
    if (!(procs[i]->CalculateTotalXSec())) okay = 0;
                                      else totalxs += procs[i]->Total();
  }
  return okay;
}

bool All_Processes::OneEvent() {
  SelectOne();
  msg.Debugging()<<"Selected Process_Group : "<<selected->Name()<<endl;
  return selected->OneEvent();
}

bool All_Processes::SameEvent() {
  if (selected) 
    return selected->SameEvent();
  msg.Error()<<" ERROR: in bool All_Processes::SameEvent() "<<endl;
  return 0;
}
