#include "ALPACA/EventGeneration/Event_Generator.H"

#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"

#include <string>
#include <vector>
#include <cmath>
#include <complex>

#include<stdio.h>

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Event_Generator::Event_Generator(shared_ptr<list<shared_ptr<Parton>>> ptr_parton_list,
                                 shared_ptr<list<pair<int,int>>> ptr_flow_backtrack):
             m_tau(0.), p_partons(ptr_parton_list),
             p_flow_backtrack(ptr_flow_backtrack),
             m_taumax(HIPars.TauMax()),
             m_only_gluons(HIPars.OnlyGluons()),
             m_tsample_min(HIPars.tsampleMin()), m_tsample_max(HIPars.tsampleMax()),
             m_show_progress_bar(HIPars.ShowProgressBar()),
             m_p_min(HIPars.pMin()), 
             m_extra_tau(1.), m_taumax_scaling_limit(HIPars.TauMaxScalingLimit()), 
             m_no_process_end(false),
             m_test_double(HIPars.TestDouble()), m_test_bool(HIPars.TestBool())
{
  //Initialize pointers
  //Lists of partons
  //p_partons = make_shared<list<shared_ptr<Parton>>>();
  p_removed_partons = make_shared<list<shared_ptr<Parton>>>();


  //Lists of taus
  p_taubars = make_shared<vector<pair<double,taubarpair>>>();
  p_tausplit = make_shared<vector<pair<double,splitpair>>>();

  //Misc.
  p_add_histo = make_shared<bool>(true);
  p_tau_restart_list = make_shared<pair<double, int>>(HIPars.TauRestart(), 1);

  //Handlers
  p_colour_handler = make_shared<Colour_Handler>(p_partons, p_flow_backtrack);
  p_analysis_handler = make_shared<Analysis_Handler>(p_partons, p_removed_partons, p_tau_restart_list, 
                                                     p_add_histo, p_taubars, p_tausplit);
  p_dyn_quant_handler = make_shared<Dynamic_Quantities_Handler>(p_partons, p_removed_partons);
  p_dyn_quant_handler->SetAnalysisHandler(p_analysis_handler);
  p_inelastic_kernel = make_shared<Inelastic_Kernel>(p_dyn_quant_handler, p_analysis_handler);
  p_elastic_kernel = make_shared<Elastic_Kernel>(p_dyn_quant_handler, p_analysis_handler);

  kinematicspair temp_placeholder;
  pair<double, double> temp_placeholder_2;
  pair<kinematicspair, pair<double, double>> new_temp_pair(temp_placeholder, temp_placeholder_2);
  p_save_merging_kinematics = make_shared<pair<shared_ptr<Parton>, kinematicspair>>(nullptr, temp_placeholder);
  p_save_split_kinematics = make_shared<pair<shared_ptr<Parton>, kinematicspair>>(nullptr, temp_placeholder);
  p_save_scatter_kinematics = make_shared<pair<kinematicspair, pair<double, double>>>(new_temp_pair);

  p_scatter_merge_handler = make_shared<Scatter_Merge_Handler>(p_partons, p_removed_partons, p_dyn_quant_handler, 
                                                               p_analysis_handler, p_inelastic_kernel, p_elastic_kernel, 
                                                               p_tau_restart_list, p_save_merging_kinematics, 
                                                               p_save_scatter_kinematics, p_taubars, p_tausplit,
                                                               p_colour_handler);
  p_split_handler = make_shared<Split_Handler>(p_partons, p_removed_partons, p_dyn_quant_handler, 
                                               p_analysis_handler, p_inelastic_kernel, p_tau_restart_list, 
                                               p_save_split_kinematics, p_taubars, p_tausplit,
                                               p_colour_handler);

  //Misc.
  if(p_tau_restart_list->first > m_taumax*m_taumax_scaling_limit){
    p_tau_restart_list->first = m_taumax*m_taumax_scaling_limit + 1.;
  }
  *p_add_histo = true;
}

Event_Generator::~Event_Generator(){
}


int Event_Generator::GenerateEvent() {

  //Reset variables from previous event
  Reset();

  //Setup new parton list
  //p_partons->insert(p_partons->begin(), ptr_parton_list->begin(), ptr_parton_list->end());

  //If timekeeper is chosen as labframe, set for all partons
  if(HIPars.Timekeeper() == 1){
    SetTimekeeperLabFrame();
  }

  //Fill parton list and initial tau vectors
  msg_Out() << METHOD << ": FillTautables" << endl;
  FillTauTables();

  //Print event info pre event
  msg_Out() << METHOD << ": PrintEventInfo" << endl;
  p_analysis_handler->PrintEventInfo(true, 0.);

  //Histogram pre event data
  msg_Out() << METHOD << ": PrePostSampling" << endl;
  p_analysis_handler->PrePostSampling(true, m_tau, p_dyn_quant_handler);

  //Run evolution of event
  msg_Out() << METHOD << ": EvolveInTau" << endl;
  double nscat = EvolveInTau();

  //Print event info post event
  msg_Out() << METHOD << ": PrintEventInfo" << endl;
  p_analysis_handler->PrintEventInfo(false, m_tau);

  //Histogram post event data
  msg_Out() << METHOD << ": PrePostSampling" << endl;
  p_analysis_handler->PrePostSampling(false, m_tau, p_dyn_quant_handler);

  return true;

}

void Event_Generator::FillTauTables() {
  if(p_partons->size() == 0){
    msg_Out() << METHOD << ": WARNING: No partons initialized, p_partons->size() = 0. Not allowed, exit." << endl;
    exit(1.);
  }

  // #### Fill table of taubars (tau for 2->2 collisions and 2->1 mergings) ####
  double taubar, dist;
  for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
    Vec4D x = (*iter1)->Position();
    Vec4D p = (*iter1)->Momentum();
    msg_Out() << "Checking parton " << (*iter1)->Number() << ", " << (*iter1)->Flav() << ", [" << (*iter1)->GetFlow(1) << ", " << (*iter1)->GetFlow(2) << "]" <<endl;
    msg_Out() << "x = [" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << "]" << endl;
    msg_Out() << "p = [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "], p^2 = " << p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3] << endl;
    for (list<shared_ptr<Parton>>::iterator iter2=iter1; iter2!=p_partons->end(); iter2++) {
      if ((*iter1)==(*iter2)) continue;
   
      std::pair<double,taubarpair> temp_pair = p_scatter_merge_handler->FindTauScatterMerge(*iter1, *iter2, m_tau);
      if (temp_pair.first > m_tau && temp_pair.first <= p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1)) {
        p_taubars->push_back(temp_pair);
      }
    }
  }
  sort(p_taubars->begin(), p_taubars->end(), sortTaubars);

  // #### Fill map with tausplit (tau for particle splitting) ####
  pair<double, splitpair> temp_splitpair;
  for (list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
    std::shared_ptr<Split_Variables> temp_split_variables = make_shared<Split_Variables>(*iter);
    temp_splitpair = p_split_handler->FindTauSplit(*iter, m_tau);
    p_tausplit->push_back(temp_splitpair);
  }
  sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);
}

void Event_Generator::ResetScatterMergeSplitList() {
  //Clear old list
  p_taubars->clear();

  // #### Fill table of taubars (tau for 2->2 collisions and 2->1 mergings) ####
  double taubar, dist;
  for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
    for (list<shared_ptr<Parton>>::iterator iter2=iter1; iter2!=p_partons->end(); iter2++) {
      if ((*iter1)==(*iter2)) continue;
   
      std::pair<double,taubarpair> temp_pair = p_scatter_merge_handler->FindTauScatterMerge(*iter1, *iter2, m_tau);
      if (temp_pair.first > m_tau && temp_pair.first <= p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1)) {
        p_taubars->push_back(temp_pair);
      }
    }
  }
  sort(p_taubars->begin(), p_taubars->end(), sortTaubars);

  //Clear old list
  p_tausplit->clear();

  pair<double, splitpair> temp_splitpair;
  for (list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
    temp_splitpair = p_split_handler->FindTauSplit(*iter, m_tau);
    p_tausplit->push_back(temp_splitpair);
  }
  sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);
}

double Event_Generator::EvolveInTau() {
  //Print progress
  float progress = 0.0;
  if(m_show_progress_bar){
    msg_Out() << "\n\n  ## Event " << rpa->gen.NumberOfGeneratedEvents()+1 << " scattering progress ##" << endl;
    ProgressBar(progress);
  }
  double even_progressbar_update = 0.5;

  double nscat(0.);
  Vec4D shift;
  bool sampled_f_midrun = false;

  //Main loop in the evolution of the parton cascade, evolves in tau
  while (m_tau < m_taumax_scaling_limit*m_taumax) {
    //Reset tau for scattering, merge and split
    double tauscat(-1.), tausplit(-1.), taumerge(-1.), dist2(-1.), taubar(-1.);

    if(m_tau >= even_progressbar_update){
      if(m_show_progress_bar){
        if(m_tau > 0.001){
          progress = m_tau/(m_taumax_scaling_limit*m_taumax);
          ProgressBar(progress);
        }
      }
      even_progressbar_update += 0.5;
    }

    int process_type = 0; //0 for no process found, 1: 2->2, 2: 2->1, 3: 1->2

    pair<double,taubarpair> taubar_pair;
    pair<double,splitpair> split_pair;
    pair<double,shared_ptr<Parton>> bounds_pair;
    pair<int, pair<double, double>> taubar_check;
    pair<bool, splitpair> tausplit_check;

    vector<pair<double,taubarpair>>::iterator scatit;
    vector<pair<double,taubarpair>>::iterator mergeit;
    vector<pair<double,splitpair>>::iterator splitit;
    vector<pair<double,shared_ptr<Parton>>>::iterator biter;

    bool process_found = false;
    while(!process_found){
      if(p_taubars->size() > 0){
        taubar_pair = (*p_taubars)[0];
        taubar = taubar_pair.first;
      } else{
        taubar = 2.*m_taumax_scaling_limit*m_taumax;
      }
      
      if(p_tausplit->size() > 0){
        split_pair = (*p_tausplit)[0];
        tausplit = split_pair.first;
      } else{
        tausplit = 2.*m_taumax_scaling_limit*m_taumax;
      }

      //No processes left before max tau, ending
      if(taubar > m_taumax_scaling_limit*m_taumax && tausplit > m_taumax_scaling_limit*m_taumax){
        m_no_process_end = true;
        return nscat;
      }

      if(taubar < m_tau || tausplit < m_tau){
        msg_Out() << METHOD << ": WARNING: taubar = " << taubar << ", tausplit = " << tausplit  << ", some time less than m_tau = " << m_tau << endl; 
        exit(1);
      }
      
      bool reset_skip = false;
      double first_tau = taubar;
      if(tausplit < first_tau) first_tau = tausplit;

      if(p_tau_restart_list->first > 0){ //If <0 then reset later after each process occurs
        if(first_tau >= p_tau_restart_list->first*double(p_tau_restart_list->second)){
          m_tau = p_tau_restart_list->first*double(p_tau_restart_list->second);
          p_tau_restart_list->second = p_tau_restart_list->second + 1;
          ResetScatterMergeSplitList();
          reset_skip = true;
        }
      }

      if(!reset_skip){
        if(taubar < tausplit){ 
          // #### Case of 2->2 or 2->1 ####
          m_tau = taubar;
          taubar_check = p_scatter_merge_handler->DoesScatterOrMerge(taubar_pair.first, taubar_pair.second);
          //msg_Out() << "taubar_check.first = " << taubar_check.first << endl;
          if (taubar_check.first != 0){ //Accepted
            if(taubar_check.first == 1){ //Case of 2->2 scattering
              //msg_Out() << "scatter" << endl;
              tauscat = taubar_pair.first;
              scatit = p_taubars->begin();
              taumerge = (m_taumax_scaling_limit+1.)*m_taumax;
              process_type = 1;
            } else{ //Case of 2->1 merging
              //msg_Out() << "\n\n\n\n #### MERGING #### \n\n\n\n" << endl;
              taumerge = taubar_pair.first;
              mergeit = p_taubars->begin();
              tauscat = (m_taumax_scaling_limit+1)*m_taumax;
              process_type = 2;
            }
            process_found = true;
          } else{ //Rejected
            std::pair<double,taubarpair> temp_taubarpair = p_scatter_merge_handler->FindTauScatterMerge(taubar_pair.second.first[0],taubar_pair.second.first[1], m_tau);
            //if (temp_taubarpair.first > m_tau && temp_taubarpair.first <= m_taumax_scaling_limit*m_taumax) {
            if (temp_taubarpair.first > m_tau && temp_taubarpair.first <= p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1)) {
              p_taubars->push_back(temp_taubarpair);
            }
            p_taubars->erase(p_taubars->begin());
            sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
          }

        } else if(tausplit < taubar){ 
          // #### Case of 1->2 ####
          m_tau = tausplit;
          tausplit_check = p_split_handler->DoesSplit(split_pair.first, split_pair.second);
          if(tausplit_check.first){ //Accepted
            process_type = 3;
            splitit = p_tausplit->begin();
            process_found = true;
          } else{ //Rejected
            pair<double, splitpair> temp_splitpair = p_split_handler->FindTauSplit(split_pair.second.first.first, m_tau);
            (*p_tausplit)[0] = temp_splitpair;
            sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);
          }

        } else{
          msg_Out() << METHOD << ": WARNING: no process found at m_tau = " << m_tau << endl;
          process_type = 0;
          process_found = true;
        }
      }
    }

    if (process_type == 0){
      //msg_Out() << METHOD << ": Some taubar for a process is not found, i.e. either 2->2, 1->2, 2->1 or boundary violation, this should not occur.\n";
      //msg_Out() << METHOD <<  "m_tau = " << m_tau << endl;
      //msg_Out() << METHOD << "tauscat = " << tauscat << ", tausplit = " << tausplit << ", taumerge = " << taumerge << endl;
      return nscat;
    }
    
    if (process_type == 3) { //Case for 1->2 collinear split
      //Do splitting
      shared_ptr<Parton> part_in = tausplit_check.second.first.first;
      shared_ptr<Parton> part_recoil = tausplit_check.second.first.second; //Ignore distance and shift for recoil, only momentum is relevant
      flavourpair out_flavours = tausplit_check.second.second.second;
      double t_split = (part_in->Position(tausplit))[0];
      
      if(!m_do_kinematics){ //Just bookkeep splittings, no kinematic update
        p_analysis_handler->BookkeepSplit(part_in, out_flavours, t_split);

        m_tau = tausplit;

        UpdateTauTables(part_in, 1);
        sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
        sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);
      } else{
        pair<double, partpair> new_partons = p_split_handler->DoSplittingKinematics(tausplit, tausplit_check.second);

        if(new_partons.first){ //Kinematics still allowed or new kinematics allowed
          p_analysis_handler->BookkeepSplit(part_in, out_flavours, t_split);

          shared_ptr<Parton> part_out_1 = new_partons.second.first;
          shared_ptr<Parton> part_out_2 = new_partons.second.second;

          m_tau = tausplit;
          if(m_show_progress_bar){
            if(m_tau > 0.001){
              progress = m_tau/(m_taumax_scaling_limit*m_taumax);
              ProgressBar(progress);
            }
          }

          UpdateTauTables(part_out_1, 1);
          UpdateTauTables(part_out_2, 1);
          UpdateTauTables(part_recoil, 1);
          sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
          sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);
        } else{ //Not possible to find new kinematics for splitting, does not split
          msg_Out() << METHOD << ": Split failed at m_tau = " << m_tau << " - for " << part_in->Number()  << " and " << part_recoil->Number() << endl;
          pair<double, splitpair> temp_splitpair = p_split_handler->FindTauSplit(split_pair.second.first.first, m_tau);
          (*p_tausplit)[0] = temp_splitpair;
          sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);
        }
      }
    }

    if (process_type == 2) { //Case for 2->1 merging
      //Do merging
      shared_ptr<Parton> part_in_1 = mergeit->second.first[0];
      shared_ptr<Parton> part_in_2 = mergeit->second.first[1];
      shared_ptr<Parton> part_recoil = p_save_merging_kinematics->first;
      double t_merge = ((part_in_1->Position(tausplit))[0] + (part_in_2->Position(tausplit))[0])/2.;

      if(!m_do_kinematics){ //Just bookkeep mergings, do not do kinematics
        p_analysis_handler->BookkeepMerge(part_in_1, part_in_2, t_merge);

        m_tau = taumerge;

        UpdateTauTables(part_in_1, 0);
        UpdateTauTables(part_in_2, 0);
        sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
      } else{
        pair<bool,shared_ptr<Parton>> merge_kinematics = p_scatter_merge_handler->DoMergingKinematics(part_in_1, part_in_2, taumerge, taubar_check.second.second);
        if(merge_kinematics.first){ //Merging still allowed
          p_analysis_handler->BookkeepMerge(part_in_1, part_in_2, t_merge);
          shared_ptr<Parton> part_out = merge_kinematics.second;

          m_tau = taumerge;

          if(m_show_progress_bar){
            if(m_tau > 0.001){
              progress = m_tau/(m_taumax_scaling_limit*m_taumax);
              ProgressBar(progress);
            }
          }

          UpdateTauTables(part_out, 1);
          UpdateTauTables(part_recoil, 1);
          sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
          sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);

        } else{ //Rejected due to kinematics
          msg_Out() << METHOD << ": ERROR: Merging kinematics failed at m_tau = " << m_tau << endl;
          std::pair<double,taubarpair> temp_taubarpair = p_scatter_merge_handler->FindTauScatterMerge(part_in_1, part_in_2, m_tau);
          //if (temp_taubarpair.first > m_tau && temp_taubarpair.first <= m_taumax_scaling_limit*m_taumax) {
          if (temp_taubarpair.first > m_tau && temp_taubarpair.first <= p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1)) {
            p_taubars->push_back(temp_taubarpair);
          }
          p_taubars->erase(p_taubars->begin());
          sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
        }
      }

    }
    
 
    if(process_type == 1){ //Case for 2->2 scattering
      //do scattering
      shared_ptr<Parton> part_in_1 = scatit->second.first[0];
      shared_ptr<Parton> part_in_2 = scatit->second.first[1];
      Flavour flav_out_1 = p_save_scatter_kinematics->first.second.first;
      Flavour flav_out_2 = p_save_scatter_kinematics->first.second.second;
      flavourpair out_flavours = make_pair(flav_out_1, flav_out_2);
      double t_scatter = ((part_in_1->Position(tausplit))[0] + (part_in_2->Position(tausplit))[0])/2.;

      if(p_scatter_merge_handler->DoScatteringKinematics(part_in_1, part_in_2, tauscat, shift)){
        p_analysis_handler->BookkeepScatter(part_in_1, part_in_2, out_flavours, t_scatter);

        if(m_show_progress_bar){
          if(m_tau > 0.001){
            progress = m_tau/(m_taumax_scaling_limit*m_taumax);
            ProgressBar(progress);
          }
        }

        UpdateTauTables(part_in_1, 1);
        UpdateTauTables(part_in_2, 1);

        sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
        sort(p_tausplit->begin(),p_tausplit->end(), sortTausplits);

      } else { //Rejected due to kinematics or other problem in p_scatter_merge_handler->DoScatteringKinematics
        msg_Out() << METHOD << ": Failed scatter at m_tau = " << m_tau << " - old particles: " << part_in_1->Number()  << ", " << part_in_1->Number() << endl;
        std::pair<double,taubarpair> temp_taubarpair = p_scatter_merge_handler->FindTauScatterMerge(part_in_1, part_in_2, m_tau);
        //if (temp_taubarpair.first > m_tau && temp_taubarpair.first <= m_taumax_scaling_limit*m_taumax) {
        if (temp_taubarpair.first > m_tau && temp_taubarpair.first <= p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1)) {
          p_taubars->push_back(temp_taubarpair);
        }
        p_taubars->erase(p_taubars->begin());
        sort(p_taubars->begin(), p_taubars->end(), sortTaubars);
      }
    }

    //Check if all particles have t > m_taumax, if not keep looping (until tau > 5.*m_taumax)
    if(m_tau >= m_extra_tau*m_taumax){
      int particles_left = 0;
      Vec4D x_temp;
      for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
        x_temp = (*iter1)->Position(m_tau);
        if(x_temp[0] < m_taumax){
          msg_Out() << "ID = " << (*iter1)->Number() << ", t = " << x_temp[0] << ", tau = " << m_tau << endl; 
          particles_left++;
        }
      }

      if(particles_left == 0){ //No particles left with t < m_taumax, end event
        m_tau = (m_taumax_scaling_limit+1)*m_taumax;
      } else{
        m_extra_tau = m_extra_tau + 0.1;
      }
    }


    if(m_tau >= m_taumax_scaling_limit*m_taumax){
      int particles_left = 0;
      Vec4D x_temp;
      for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
        x_temp = (*iter1)->Position(m_tau);
        if(x_temp[0] < m_taumax){
          particles_left++;
        }
      }
    }

    if(p_tau_restart_list->first < 0){ //If <0 then reset after each process occurs
      ResetScatterMergeSplitList();
    }

  }
  return nscat;
}


void Event_Generator::UpdateTauTables(shared_ptr<Parton> part1, int b_type) {
  /*
    Updates the three tau tables, scattering, splitting and boundary violations for part1.
    b_type indicates which kind of update should be done:
      0: Updates all lists except tausplit (case for boundary violation, 2->2 scattering, 1->2 splitting recoiler and 2->1 merging recoiler)
      1: Updates all lists including new tausplit for part1 (Case for new parton in 1->2 splitting and new parton in 2->1 merging)
      2: Updates only tau_split list for part1 (Case for rejected splitting)
  */

  if(part1->GetSplitMergeStatus()){
    msg_Out() << METHOD << ": WARNING: updating tables for particle " << part1->Number() << " with SplitMergeStatus = true" << endl;
  }

  if(b_type != 2){
    //Update table with scattering taus, only entries with part1 present
    for (vector<pair<double,taubarpair>>::iterator mapit=p_taubars->begin(); mapit!=p_taubars->end(); mapit++){
      //Delete all instances of part1
      if ((mapit->second.first[0]->Number()==part1->Number()) || (mapit->second.first[1]->Number()==part1->Number())) {
        p_taubars->erase(mapit--);
      } else if(mapit->first < m_tau){
        p_taubars->erase(mapit--);
      }
    }
    //Find and add all taubars for part1 vs all other particles in p_partons
    for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
      if ((*iter1)==part1) continue; //Do not include closest approach to itself
      if ((*iter1)==part1->GetLastScatter()) continue; //Do not include closest approach to last elastic scattering particle
      //if ((*iter1)==part1->GetSecondLastScatter()) continue; //Do not include closest approach to second last elastic scattering particle
      //if ((*iter1)==part1->GetThirdLastScatter()) continue; //Do not include closest approach to third last elastic scattering particle
      //if ((*iter1)==part1->GetSplitPartner()) continue; //Do not include closest approach to splitting partner

      std::pair<double,taubarpair> temp_pair = p_scatter_merge_handler->FindTauScatterMerge(part1, *iter1, m_tau);

      if (temp_pair.first > m_tau && temp_pair.first <= p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1) ) {
        p_taubars->push_back(temp_pair);
      }
    }
  }
  
  //Update tausplit in p_tausplit for part1
  if(b_type == 1 || b_type == 2){ //Only update splitting for a parton when that parton splits
    pair<double, splitpair> temp_splitpair;
    for (vector<std::pair<double,splitpair>>::iterator mapit2=p_tausplit->begin(); mapit2!=p_tausplit->end(); mapit2++){
      if ((mapit2->second.first.first)==part1){
        temp_splitpair = p_split_handler->FindTauSplit(mapit2->second.first.first, m_tau);
        *mapit2 = temp_splitpair;
      }
    }
  }
}

void Event_Generator::SetTimekeeperLabFrame(){
  Vec4D a_hat(1.,0.,0.,0.);
  /*
  if(m_timekeeper == 2){
    a_hat = m_P_global/(sqrt(m_P_global*m_P_global));
  } else if (m_timekeeper != 1){
    msg_Out() << METHOD << ": WARNING: m_timekeeper != 1 or 2 but trying to set a_hat" << endl;
  }
  */
 for (list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
    (*iter)->SetTimekeeper(1, a_hat);
  }
}



void Event_Generator::Reset() {
  //p_partons->clear();
  p_removed_partons->clear();
  p_taubars->clear();
  p_tausplit->clear();
  m_tau = 0.;
  m_extra_tau = 1.;
  p_tau_restart_list->second = 1;
  m_no_process_end = false;

  p_analysis_handler->ResetEventVariables();
}


void Event_Generator::CleanUp(const size_t & mode) {
}



void Event_Generator::ProgressBar(float progress){
  /*
    Function used to update the progress bar for each event
  */
  int pos, barWidth = 70;
  cout << "  [";
  pos = barWidth * progress;
  for (int ibar = 0; ibar < barWidth; ++ibar) {
      if (ibar < pos) cout << "=";
      else if (ibar == pos) cout << "|";
      else cout << " ";
  }

  cout << "] " << int(progress * 100.0) << " %, tau = " << m_tau << ", p_parton = " << p_partons->size() << ", p_taubar = " << p_taubars->size();
  cout << ", p_tausplit = " << p_tausplit->size() << ", p_removed_partons = " << p_removed_partons->size() <<  " \r";
  cout.flush();
}