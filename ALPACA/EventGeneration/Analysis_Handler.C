#include "ALPACA/EventGeneration/Analysis_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Model_Base.H"

#include <string>
#include <vector>
#include <cmath>
#include <stdio.h>


using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Analysis_Handler::Analysis_Handler(shared_ptr<list<shared_ptr<Parton>>> ptr_partons, 
                                   shared_ptr<list<shared_ptr<Parton>>> ptr_removed_partons,
                                   shared_ptr<pair<double, int>> ptr_tau_restart_list,
                                   shared_ptr<bool> ptr_add_histo,
                                   shared_ptr<vector<pair<double,taubarpair>>> ptr_taubars,
                                   shared_ptr<vector<pair<double,splitpair>>> ptr_tausplits):
  p_partons(ptr_partons), p_removed_partons(ptr_removed_partons), 
  p_tau_restart_list(ptr_tau_restart_list),
  p_taubars(ptr_taubars), p_tausplits(ptr_tausplits),
  p_add_histo(ptr_add_histo),
  m_show_event_information(HIPars.ShowEventInformation()),
  m_seed_str("_s" + to_string(ran->GetSeed())), m_p_min(HIPars.pMin()),
  m_tsample_min(HIPars.tsampleMin()), m_tsample_max(HIPars.tsampleMax()),
  m_N_scatter_event(0), m_N_split_event(0), m_N_merge_event(0),
  m_v2_w_total(0.), m_v2_total(0.)
{
  for (int i = 0; i < 3; ++i) {
    m_N_split_process_event[i] = 0;
    m_N_merge_process_event[i] = 0;
  }
  for (int i = 0; i < 7; ++i) {
    m_N_scatter_process_event[i] = 0;
  }

  //Initialize histograms
  p_histomap = make_shared<map<string, shared_ptr<ATOOLS::Histogram>>>();
  p_histomap2D = make_shared<map<string, shared_ptr<ATOOLS::Histogram_2D>>>();

  msg_Out() << METHOD << ": Initializing analysis handler" << endl;
  msg_Out() << "#### HIPars test = " << HIPars.ShowEventInformation() << endl;

  //1D Histograms
  int N_process_max = 50;
  int N_process_bin = 50;
  (*p_histomap)[string("N_scatter" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,N_process_max,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_gg_gg" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_qqb_gg" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_gq_gq" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_qiqi_qiqi" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_qiqib_qiqib" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_qiqib_qjqjb" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_scatter_ratio_qiqj_qiqj" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_split" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,N_process_max,N_process_bin);
  (*p_histomap)[string("N_split_ratio_ggg" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_split_ratio_gqq" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_split_ratio_qgq" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_merge" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,N_process_max,N_process_bin);
  (*p_histomap)[string("N_merge_ratio_ggg" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_merge_ratio_gqq" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);
  (*p_histomap)[string("N_merge_ratio_qgq" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,1.,N_process_bin);

  (*p_histomap)[string("t_scatter" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,m_tsample_max,50);
  (*p_histomap)[string("t_split" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,m_tsample_max,50);
  (*p_histomap)[string("t_merge" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,m_tsample_max,50);

  (*p_histomap)[string("N_overdraw_split" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,100,100);
  (*p_histomap)[string("N_overdraw_scatter_merge" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,100,100);
  (*p_histomap)[string("val_overdraw_split" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,100,100);
  (*p_histomap)[string("val_overdraw_scatter_merge" + m_seed_str)] = make_shared<ATOOLS::Histogram>(0,0.,100,100);

  //2D histograms
  (*p_histomap2D)[string("t+mg2" + m_seed_str)] = std::make_shared<ATOOLS::Histogram_2D>(0,0.,m_tsample_max,25,0.,1.,100);
  (*p_histomap2D)[string("t+mq2" + m_seed_str)] = std::make_shared<ATOOLS::Histogram_2D>(0,0.,m_tsample_max,25,0.,1.,100);
  (*p_histomap2D)[string("t+Tstar" + m_seed_str)] = std::make_shared<ATOOLS::Histogram_2D>(0,0.,m_tsample_max,25,0.,1.,100);
};

Analysis_Handler::~Analysis_Handler() {
  //Finalize histograms before removing Analysis_Handler
  if (!(*p_histomap).empty()) {
    shared_ptr<Histogram> histo;
    string name;
    for (map<string,shared_ptr<ATOOLS::Histogram>>::iterator hit = (*p_histomap).begin(); hit != (*p_histomap).end(); hit++) {
      histo = hit->second;
      name = string("HI_Analysis/")+hit->first+string(".dat");
      //histo->Finalize();
      histo->Output(name);
    }
    (*p_histomap).clear();
  } else{
    msg_Out() << "Empty histogram map." << endl;
  }
  if (!(*p_histomap2D).empty()) {
    shared_ptr<Histogram_2D> histo;
    string name;
    for (map<string,shared_ptr<ATOOLS::Histogram_2D>>::iterator hit = (*p_histomap2D).begin(); hit != (*p_histomap2D).end(); hit++) {
      histo = hit->second;
      name = string("HI_Analysis2D/")+hit->first+string(".dat");
      //histo->Finalize();
      histo->Output(name);
    }
    (*p_histomap2D).clear();
  }
  else{
    msg_Out() << "Empty histogram2D map." << endl;
  }
};

void Analysis_Handler::Histo(std::string histo_name, double histo_val){
  //Check if histogram exists
  if(p_histomap->find(histo_name + m_seed_str) != p_histomap->end()){
    if(*p_add_histo){
      //Check if value is outside of bin max/min for histogram. If so, set it to corresponding boundary.
      if(histo_val <= (*p_histomap)[string(histo_name + m_seed_str)]->Xmin()) histo_val = 1.01*((*p_histomap)[string(histo_name + m_seed_str)]->Xmin());
      if(histo_val >= (*p_histomap)[string(histo_name + m_seed_str)]->Xmax()) histo_val = 0.99*((*p_histomap)[string(histo_name + m_seed_str)]->Xmax());
      //Add entry to histogram
      (*p_histomap)[string(histo_name + m_seed_str)]->Insert(histo_val);
    }
  } else{
    msg_Out() << METHOD << ": WARNING: histogram with name " << histo_name << " does not exist, will exit()." << endl;
    exit(1.);
  }
};

void Analysis_Handler::Histo2D(std::string histo_name, double histo_val_1, double histo_val_2){
  //Check if histogram exists
  if(p_histomap2D->find(histo_name + m_seed_str) != p_histomap2D->end()){
    if(*p_add_histo){
      //Check if value is outside of bin max/min for histogram. If so, set it to corresponding boundary.
      if(histo_val_1 <= (*p_histomap2D)[string(histo_name + m_seed_str)]->Xmin()) histo_val_1 = 1.01*((*p_histomap2D)[string(histo_name + m_seed_str)]->Xmin());
      if(histo_val_1 >= (*p_histomap2D)[string(histo_name + m_seed_str)]->Xmax()) histo_val_1 = 0.99*((*p_histomap2D)[string(histo_name + m_seed_str)]->Xmax());
      if(histo_val_2 <= (*p_histomap2D)[string(histo_name + m_seed_str)]->Ymin()) histo_val_2 = 1.01*((*p_histomap2D)[string(histo_name + m_seed_str)]->Ymin());
      if(histo_val_2 >= (*p_histomap2D)[string(histo_name + m_seed_str)]->Ymax()) histo_val_2 = 0.99*((*p_histomap2D)[string(histo_name + m_seed_str)]->Ymax());
      //Add entry to histogram
      (*p_histomap2D)[string(histo_name + m_seed_str)]->Insert(histo_val_1, histo_val_2);
    }
  } else{
    msg_Out() << METHOD << ": WARNING: histogram2D with name " << histo_name << " does not exist, will exit()." << endl;
    exit(1.);
  }
};

void Analysis_Handler::BookkeepScatter(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, flavourpair out_flavours, double t){
  //Function to bookkeep each elastic scatter in the Event_Generator, including which type it is
  //for histograming once the event is done

  Flavour flav_in_1 = part_in_1->Flav();
  Flavour flav_in_2 = part_in_2->Flav();
  Flavour flav_out_1 = out_flavours.first;
  Flavour flav_out_2 = out_flavours.second;
  
  //Find which type of elastic scatter it is
  int process_type_scatter = -1;
  if((flav_in_1.IsGluon() && flav_in_2.IsGluon()) || (flav_out_1.IsGluon() && flav_out_2.IsGluon())){ 
    if(flav_in_1.IsGluon() && flav_in_2.IsGluon() && flav_out_1.IsGluon() && flav_out_2.IsGluon()){ 
      process_type_scatter = 0; // gg->gg
    } else{ 
      process_type_scatter = 1; //qqb->gg
    }  
  } else if((flav_in_1.IsGluon() && flav_in_2.IsFermion()) || (flav_in_1.IsFermion() && flav_in_2.IsGluon())){ 
    process_type_scatter = 2; //gq->gq
  } else{
    if(flav_in_1 == flav_in_2 && flav_in_1 == flav_out_1){
      process_type_scatter = 3; //qiqi->qiqi
    } else if(flav_in_1 == flav_in_2.Bar()){
      if(flav_in_1 == flav_out_1 || flav_in_1 == flav_out_2){
        process_type_scatter = 4; //qiqib->qiqib
      } else{
        process_type_scatter = 5; //qiqib->qjqjb
      }
    } else{
      process_type_scatter = 6; //qiqj->qiqj
    }
  }

  //Bookkeep scatter
  if(process_type_scatter >= 0){
    m_N_scatter_event++;
    m_N_scatter_process_event[process_type_scatter]++;
    Histo("t_scatter", t);
  } else{
    msg_Out() << METHOD << ": WARNING: could not find process type." << endl;
  }
};

void Analysis_Handler::BookkeepMerge(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double t){
  //Function to bookkeep each inelastic merging in the Event_Generator, including which type it is
  //for histograming once the event is done
  Flavour flav_in_1 = part_in_1->Flav();
  Flavour flav_in_2 = part_in_2->Flav();
  
  int process_type_merge;
  if(flav_in_1.IsGluon() && flav_in_2.IsGluon()){ // gg->g
    process_type_merge = 0;
  } else if(flav_in_1.IsQuark() && flav_in_2.IsQuark()){ // qqb -> g
    process_type_merge = 1;
  } else{ // gq->q
    process_type_merge = 2;
  }
  m_N_merge_event++;
  m_N_merge_process_event[process_type_merge]++;
  Histo("t_merge", t);
};

void Analysis_Handler::BookkeepSplit(std::shared_ptr<Parton> part_in, flavourpair out_flavours, double t){
  //Function to bookkeep each inelastic splitting in the Event_Generator, including which type it is
  //for histograming once the event is done
  Flavour flav_in = part_in->Flav();
  int process_type_split;
  if(flav_in.IsGluon()){ 
    if(out_flavours.first.IsGluon()){ // g->gg
      process_type_split = 0;
    } else{ // g->qqb
      process_type_split = 1;
    }
  } else{ // q->gq
    process_type_split = 2;
  }

  m_N_split_event++;
  m_N_split_process_event[process_type_split]++;
  Histo("t_split", t);
};

void Analysis_Handler::BookkeepOverdraw(double val, bool is_split){
  if(is_split){
    m_N_split_overdraw_event.push_back(val);
    Histo("val_overdraw_split", val);
  } else{
    m_N_scatter_merge_overdraw_event.push_back(val);
    Histo("val_overdraw_scatter_merge", val);
  }
};

void Analysis_Handler::ResetEventVariables() {
  //Called once an event is done to reset variables that just
  //bookkeep on an event-to-event basis
  m_N_scatter_event = 0;
  m_N_split_event = 0;
  m_N_merge_event = 0;
  for (int i = 0; i < 3; ++i) {
    m_N_split_process_event[i] = 0;
    m_N_merge_process_event[i] = 0;
  }
  for (int i = 0; i < 7; ++i) {
    m_N_scatter_process_event[i] = 0;
  }
  m_N_split_overdraw_event.clear();
  m_N_scatter_merge_overdraw_event.clear();
};


void Analysis_Handler::PrintEventInfo(bool pre_event, double tau){
  //Print general run information pre and post event

  if(pre_event){
    // #### Pre Event ####
    //if(m_show_event_information){
    if(true){
      msg_Out() << "\n\n\n#### New event (" << int(rpa->gen.NumberOfGeneratedEvents())+1 << " of " << rpa->gen.NumberOfEvents() << ") ####\n" << endl;
      msg_Out() << "  ## General info ##" << endl;
      msg_Out() << "  Seed = " << to_string(ran->GetSeed()) << endl;
      msg_Out() << "  Seed_str = " << m_seed_str << endl;
      msg_Out() << "  add_histo = " << *p_add_histo << endl;
      msg_Out() << "  timekeeper = " << HIPars.Timekeeper() << endl;
      if(HIPars.Timekeeper() == 1){
        msg_Out() << "    Lambda = " << HIPars.Lambda() << endl;
      }
      msg_Out() << "  taumax = " << HIPars.TauMax() << endl;
      msg_Out() << "  alphaS = " << HIPars.AlphaS() << endl;
      if(HIPars.XSec_Form() == xsec_form::BlackDisc) {
        msg_Out() << "  Xsec_Form = BlackDisc" << endl;
      } else {
        msg_Out() << "  Xsec_Form = Gauss" << endl;
      }
      msg_Out() << "  fixed_sigma = [" << HIPars.FixedSigma().first << ", " << HIPars.FixedSigma().second << "]" << endl;
      msg_Out() << "  splitting_merging = " << HIPars.SplittingMerging() << endl;
      if(HIPars.SplittingMerging()){
        msg_Out() << "    formation_time = " << HIPars.FormationTime() << endl;
      }
      msg_Out() << "  OE_mult_scatter = " << HIPars.OEMultScatter() << endl;
      msg_Out() << "  OE_mult_merge = " << HIPars.OEMultMerge() << endl;
      msg_Out() << "  OE_mult_split = " << HIPars.OEMultSplit() << endl;
      msg_Out() << "  gaussian_kT2 = " << HIPars.GaussiankT2() << endl;
      msg_Out() << "  kT2_reg = " << HIPars.kT2Reg() << endl;
      msg_Out() << "  fixed_gamma = [" << HIPars.FixedGamma().first << ", " << HIPars.FixedGamma().second << "]" << endl;
      msg_Out() << "  p_min = " << HIPars.pMin() << endl;
      msg_Out() << "  only gluons = " << HIPars.OnlyGluons() << endl;
      msg_Out() << "  tau_reset = [" << p_tau_restart_list->first << ", " << p_tau_restart_list->second << "]" << endl;
      msg_Out() << "  include_bose_factors = " << HIPars.IncludeBoseFactors() << endl;
      msg_Out() << "  f_r = " << HIPars.fr() << endl;
      msg_Out() << "  f_delta_p = " << HIPars.fDeltap() << endl;
      msg_Out() << "  f_shell = " << HIPars.fShell() << endl;
      msg_Out() << "  f_N_max = " << HIPars.fNMax() << endl; 
      msg_Out() << "  m2_N_include = " << HIPars.NInclude() << endl;
      msg_Out() << "  m2_min_scale = " << HIPars.M2MinScale() << endl;

      msg_Out() << "  test_double = " << HIPars.TestDouble() << endl;
      msg_Out() << "  test_bool = " << HIPars.TestBool() << endl;


      msg_Out() << "\n  ## Pre collsions ##" << endl;
      msg_Out() << "  p_partons->size() = " << p_partons->size() << endl;
      msg_Out() << "  p_taubars->size() = " << p_taubars->size() << endl;
      msg_Out() << "  p_tausplits->size() = " << p_tausplits->size() << endl;
    }


    double px_pre = 0, py_pre = 0, pz_pre = 0;
    double total_energy_pre = 0.;
    bool no_error = true;
    int N_f = 0, N_b = 0;
    for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
      Vec4D p = (*iter1)->Momentum();
      px_pre += p[1];
      py_pre += p[2];
      pz_pre += p[3];
      total_energy_pre = total_energy_pre + p[0];
      if(p.Abs2() > 1.e-6){
        msg_Out() << "  Error, parton not massless, p.abs2() = " << p.Abs2() << endl;
        no_error = false;
      }

      //Count fermion, boson
      if((*iter1)->Flav().IsFermion()){
        N_f++;
      } else{
        N_b++;
      }
    }
    if(m_show_event_information){
      msg_Out() << "  Number of fermions = " << N_f << endl;
      msg_Out() << "  Number of bosons = " << N_b << endl;
      msg_Out() << "  Total energy of the system = " << total_energy_pre << endl;
      msg_Out() << "  Total 3 momentum of the system = (" << px_pre << ", " << py_pre << ", " << pz_pre << ")" << endl;
      if(no_error){
        msg_Out() << "  Particle 4 momentum check ok, pp=0 for all particles." << endl;
      } else{
        msg_Out() << "  WARNING: not all particles massless, tolerance 1.e-6 " << endl;
      }
      CheckIfSinglet();
    }

    

  } else{
    // #### Post Event ####
    if(m_show_event_information){
      msg_Out() << "\n\n  ## Post collisions ##" << endl;
      if(HIPars.SplittingMerging()){
        msg_Out() << "  Splitting/merging statistics between between tmin = " << m_tsample_min << " and tmax = " << m_tsample_max << endl;
        msg_Out() << "  N_split = " << m_N_split_event  << endl;
        msg_Out() << "    g->gg = " << m_N_split_process_event[0] << endl;
        msg_Out() << "    g->qqb = " << m_N_split_process_event[1] << endl;
        msg_Out() << "    q->gq = " << m_N_split_process_event[2] << endl;
        
        msg_Out() << "  N_merge = " << m_N_merge_event  << endl;
        msg_Out() << "    gg->g = " << m_N_merge_process_event[0] << endl;
        msg_Out() << "    qqb->g = " << m_N_merge_process_event[1] << endl;
        msg_Out() << "    qg->q = " << m_N_merge_process_event[2] << endl;
        
      }
      if(HIPars.ElasticScattering()){
        msg_Out() << "  N_scatter = " << m_N_scatter_event << endl;
        msg_Out() << "    gg<->gg = " << m_N_scatter_process_event[0] << endl;
        msg_Out() << "    gg<->qqb = " << m_N_scatter_process_event[1] << endl;
        msg_Out() << "    gq<->gq = " << m_N_scatter_process_event[2] << endl;
        msg_Out() << "    qiqi<->qiqi = " << m_N_scatter_process_event[3] << endl;
        msg_Out() << "    qiqib<->qiqib = " << m_N_scatter_process_event[4] << endl;
        msg_Out() << "    qiqib<->qjqjb = " << m_N_scatter_process_event[5] << endl;
        msg_Out() << "    qiqj<->qiqj = " << m_N_scatter_process_event[6] << endl;
      }
      msg_Out() << "  N_scatter_merge_overdraw = " << m_N_scatter_merge_overdraw_event.size()  << endl;
      msg_Out() << "  N_split_overdraw = " << m_N_split_overdraw_event.size()  << endl;
      msg_Out() << "  p_partons->size() = " << p_partons->size() << endl;
      msg_Out() << "  p_taubars->size() = " << p_taubars->size() << endl;
      msg_Out() << "  p_tausplits->size() = " << p_tausplits->size() << endl;
      msg_Out() << "  tau_reset = [" << p_tau_restart_list->first << ", " << p_tau_restart_list->second << "]" << endl;
    }

    double px = 0, py = 0, pz = 0;
    double total_energy_post = 0.;
    bool no_error = true;
    int N_f = 0;
    int N_b = 0;
    for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
      Vec4D p = (*iter1)->Momentum();
      Vec3D p3(p[1],p[2],p[3]);
      px += p[1];
      py += p[2];
      pz += p[3];
      total_energy_post = total_energy_post + p[0];
      if(p.Abs2() > 1.e-6){
        msg_Out() << "  Error, parton not massless, p.abs2() = " << p.Abs2() << endl;
        no_error = false;
      }
      //Count fermion, boson
      if((*iter1)->Flav().IsFermion()){
        N_f++;
      } else if((*iter1)->Flav().IsBoson()){
        N_b++;
      }
    }

    if(m_show_event_information){
      msg_Out() << "  Number of fermions = " << N_f << endl;
      msg_Out() << "  Number of bosons = " << N_b << endl;
      msg_Out() << "  Total energy of the system = " << total_energy_post << endl;
      msg_Out() << "  Total 3 momentum of the system = (" << px << ", " << py << ", " << pz << ")" << endl;
      if(no_error){
        msg_Out() << "  Particle 4 momentum check ok, pp=0 for all particles." << endl;
      } else{
        msg_Out() << "  WARNING: not all particles massless, tolerance 1.e-6 " << endl;
      }
    }
    
    CheckIfSinglet();
  }
}

void Analysis_Handler::CheckIfSinglet(){
  int flow_id_1, flow_id_2;
  int found_id_1, found_id_2;
  int found_err_1, found_err_2;
  bool all_found = true;
  for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
    flow_id_1 = (*iter1)->GetFlow(1);
    flow_id_2 = (*iter1)->GetFlow(2);
    found_id_1 = 0;
    found_id_2 = 0;
    found_err_1 = 0;
    found_err_2 = 0;
    for (list<shared_ptr<Parton>>::iterator iter2=p_partons->begin(); iter2!=p_partons->end(); iter2++) {
      if((*iter1) == (*iter2)) continue;

      if(flow_id_1 == (*iter2)->GetFlow(2) && flow_id_1 != 0) found_id_1++;
      if(flow_id_2 == (*iter2)->GetFlow(1) && flow_id_2 != 0) found_id_2++;
      if(flow_id_1 == (*iter2)->GetFlow(1) && flow_id_1 != 0) found_err_1++;
      if(flow_id_2 == (*iter2)->GetFlow(2) && flow_id_2 != 0) found_err_2++;
    }

    if(found_id_1 != 1 && flow_id_1 != 0){
      msg_Out() << "  WARNING: parton " << (*iter1)->Number() << ", flow_id = " << flow_id_1 << " has " << found_id_1 << " anti-partners, logged as external flow, [" << flow_id_1 << ", " << flow_id_2 << "]" << endl;
      all_found = false;
    }
    if(found_id_2 != 1 && flow_id_2 != 0){
      msg_Out() << "  WARNING: parton " << (*iter1)->Number() << ", flow_id = " << flow_id_2 << " has " << found_id_2 << " anti-partners, logged as external flow, [" << flow_id_1 << ", " << flow_id_2 << "]" << endl;
      all_found = false;
    }
    if(found_err_1 != 0 && flow_id_1 != 0){
      msg_Out() << "  WARNING: flow_id = " << flow_id_1 << " has " << found_err_1 << " copies" << endl;
      all_found = false;
    }
    if(found_err_2 != 0 && flow_id_2 != 0){
      msg_Out() << "  WARNING: flow_id = " << flow_id_2 << " has " << found_err_2 << " copies" << endl;
      all_found = false;
    }
  }

  if(all_found){
    msg_Out() << "  All colour indices match, total is singlet" << endl;
  } else{
    msg_Out() << "  Not a singlet in total" << endl;
    //exit(1.);
  }
}

void Analysis_Handler::PrePostSampling(bool pre_event, double tau, shared_ptr<Dynamic_Quantities_Handler> dyn_quant_handler){
  //Sample certain quantitites from the parton ensamble, pre and post each event.
  double f_psd, px_psd, py_psd, pz_psd, x_psd, y_psd, z_psd;
  Vec4D temp_p_vec;

  if(pre_event){
    // #### Pre Event ####
    
  } else{
    // #### Post Event ####

    //Bookkeep event scatterings, splitting and mergings
    Histo("N_scatter", m_N_scatter_event);
    Histo("N_scatter_ratio_gg_gg", double(m_N_scatter_process_event[0])/double(m_N_scatter_event));
    Histo("N_scatter_ratio_qqb_gg", double(m_N_scatter_process_event[1])/double(m_N_scatter_event));
    Histo("N_scatter_ratio_gq_gq", double(m_N_scatter_process_event[2])/double(m_N_scatter_event));
    Histo("N_scatter_ratio_qiqi_qiqi", double(m_N_scatter_process_event[3])/double(m_N_scatter_event));
    Histo("N_scatter_ratio_qiqib_qiqib", double(m_N_scatter_process_event[4])/double(m_N_scatter_event));
    Histo("N_scatter_ratio_qiqib_qjqjb", double(m_N_scatter_process_event[5])/double(m_N_scatter_event));
    Histo("N_scatter_ratio_qiqj_qiqj", double(m_N_scatter_process_event[6])/double(m_N_scatter_event));

    Histo("N_merge", m_N_merge_event);
    Histo("N_merge_ratio_ggg", double(m_N_merge_process_event[0])/double(m_N_merge_event));
    Histo("N_merge_ratio_gqq", double(m_N_merge_process_event[1])/double(m_N_merge_event));
    Histo("N_merge_ratio_qgq", double(m_N_merge_process_event[2])/double(m_N_merge_event));

    Histo("N_split", m_N_split_event);
    Histo("N_split_ratio_ggg", double(m_N_split_process_event[0])/double(m_N_merge_event));
    Histo("N_split_ratio_gqq", double(m_N_split_process_event[1])/double(m_N_merge_event));
    Histo("N_split_ratio_qgq", double(m_N_split_process_event[2])/double(m_N_merge_event));

    //Bookkeep number of times accept/reject ratio was > 1 in the event
    Histo("N_overdraw_split", m_N_split_overdraw_event.size());
    Histo("N_overdraw_scatter_merge", m_N_scatter_merge_overdraw_event.size());
    
    list<shared_ptr<Parton>> all_partons = *p_partons;
    all_partons.insert(all_partons.end(), p_removed_partons->begin(), p_removed_partons->end() );

    double pT_sum = 0.;
    double v2_sum = 0.;
    double v2_sum_w = 0.;
    int N_v2_added = 0;
    int N_scattering_rate_added = 0;
    double scattering_rate_tot = 0.;
    double scattering_rate_tot_v2 = 0.;
    pair<bool, pair<pair<Vec4D, Vec4D>, Flavour>> xp_history_v2;
    std::vector<std::pair<int, std::pair<std::pair<ATOOLS::Vec4D,ATOOLS::Vec4D>, ATOOLS::Flavour>>> xp_history_full;
    for (list<shared_ptr<Parton>>::iterator iter1=all_partons.begin(); iter1!=all_partons.end(); iter1++) {
      //Histogram mean free path
      xp_history_full = (*iter1)->GetAllXPHistory();
      double N = xp_history_full.size() - 2.;
      double delta_t_scatter = 0.;
      N_scattering_rate_added++;
      if(N > 0.5){
        delta_t_scatter = xp_history_full[xp_history_full.size()-1].second.first.first[0] - xp_history_full[1].second.first.first[0];
        scattering_rate_tot = scattering_rate_tot + 1./(N/delta_t_scatter);
        scattering_rate_tot_v2 = scattering_rate_tot_v2 + N/delta_t_scatter;
        //if(*p_add_histo) (*p_histomap)[string("mfp_v2" + m_seed_str)]->Insert(1./(N/delta_t_scatter));
      }
      
      xp_history_v2 = (*iter1)->GetXP(m_tsample_max);
      if(xp_history_v2.first){
        Vec4D temp_p = xp_history_v2.second.first.second;
        double pT = sqrt(temp_p[1]*temp_p[1] + temp_p[2]*temp_p[2]);
        double phi_p = atan2(temp_p[2], temp_p[1]);
        v2_sum = v2_sum + cos(2.*phi_p);
        v2_sum_w = v2_sum_w + pT*cos(2.*phi_p);
        pT_sum = pT_sum + pT;
        N_v2_added++;
      }
    }

    scattering_rate_tot = scattering_rate_tot/double(N_scattering_rate_added);
    scattering_rate_tot_v2 = scattering_rate_tot_v2/double(N_scattering_rate_added);
    //if(*p_add_histo) (*p_histomap)[string("mfp_v3" + m_seed_str)]->Insert(scattering_rate_tot);
    //if(*p_add_histo) (*p_histomap)[string("mfp_v4" + m_seed_str)]->Insert(1./scattering_rate_tot_v2);

    if(N_v2_added > 0){
      v2_sum = v2_sum/double(N_v2_added);
    }
    if(pT_sum > 0.){
      v2_sum_w = v2_sum_w/pT_sum;
    }


    double N_tot_events = rpa->gen.NumberOfEvents();

    if(*p_add_histo){
      m_v2_total = m_v2_total + v2_sum/N_tot_events;
      m_v2_w_total = m_v2_w_total + v2_sum_w/N_tot_events;
    }

    if(m_show_event_information) msg_Out() << "  v2_w = " << v2_sum_w << endl;

    if(int(rpa->gen.NumberOfGeneratedEvents())+1 == rpa->gen.NumberOfEvents()){ //If last event, output v2 data into separate files
      msg_Out() << "  Last event, printing v2_w and v2 values to separate file." << endl;

      ofstream outfile_g;
      int r_1 = std::remove(("v2/v2" + m_seed_str + ".dat").c_str());
      outfile_g.open(string("v2/v2" + m_seed_str + ".dat"), std::ios::app);
      outfile_g << m_v2_w_total << " " << m_v2_total;
      outfile_g.close();
    }
  }
}