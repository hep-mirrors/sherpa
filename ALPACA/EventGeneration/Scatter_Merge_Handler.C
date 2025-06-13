#include "ALPACA/EventGeneration/Scatter_Merge_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"


#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <stdio.h>

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Scatter_Merge_Handler::Scatter_Merge_Handler(shared_ptr<list<shared_ptr<Parton>>> ptr_partons,
                                             shared_ptr<list<shared_ptr<Parton>>> ptr_removed_partons,
                                             shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler,
                                             shared_ptr<Analysis_Handler> ptr_analysis_handler,
                                             shared_ptr<Inelastic_Kernel> ptr_inelastic_kernel,
                                             shared_ptr<Elastic_Kernel> ptr_elastic_kernel,
                                             shared_ptr<pair<double, int>> ptr_tau_restart_list,
                                             shared_ptr<pair<shared_ptr<Parton>, kinematicspair>>  ptr_save_merging_kinematics,
                                             shared_ptr<pair<kinematicspair, pair<double,double>>>  ptr_save_scatter_kinematics,
                                             shared_ptr<vector<pair<double,taubarpair>>> ptr_taubars,
                                             shared_ptr<vector<pair<double,splitpair>>> ptr_tausplits,
                                             shared_ptr<Colour_Handler> ptr_colour_handler):
  p_partons(ptr_partons), p_removed_partons(ptr_removed_partons), 
  p_dyn_quant_handler(ptr_dyn_quant_handler), 
  p_colour_handler(ptr_colour_handler),
  p_inelastic_kernel(ptr_inelastic_kernel), 
  p_elastic_kernel(ptr_elastic_kernel),
  p_analysis_handler(ptr_analysis_handler),
  p_taubars(ptr_taubars), p_tausplits(ptr_tausplits),
  p_tau_restart_list(ptr_tau_restart_list),
  m_splitting_merging(HIPars.SplittingMerging()), 
  m_taumax(HIPars.TauMax()), 
  m_only_gluons(HIPars.OnlyGluons()), m_p_min(HIPars.pMin()), 
  m_fixed_sigma(HIPars.FixedSigma()),
  m_OE_mult_scatter(HIPars.OEMultScatter()), m_OE_mult_merge(HIPars.OEMultMerge()),
  m_timekeeper(HIPars.Timekeeper()), m_include_bose_factors(HIPars.IncludeBoseFactors()),
  m_elastic_scattering(HIPars.ElasticScattering()), m_taumax_scaling_limit(HIPars.TauMaxScalingLimit()),
  p_save_merging_kinematics(ptr_save_merging_kinematics), p_save_scatter_kinematics(ptr_save_scatter_kinematics)
{};


Scatter_Merge_Handler::~Scatter_Merge_Handler() {
};


pair<double,taubarpair> Scatter_Merge_Handler::FindTauScatterMerge(shared_ptr<Parton> part_in_1, shared_ptr<Parton> part_in_2, double tau){
  /*
    Function to find taubar for closest approach of a particle pair where they either scatter or merge.
    Cross sections use are overestimates in terms of Bose/Pauli factors, and the variables m_g^2, m_f^2 and T^*
    
    Returns:
      [taubar, [xsec used, [type of process(0: scatter, 1: merge), formation tau (0 if scatter)]]]
  */
  //msg_Out() << "Enter FindTauScatterMerge" << endl;
  pair<double,taubarpair> retValue;
  vector<shared_ptr<Parton>> partvec = {part_in_1, part_in_2};
  //Remove any bias that comes from initializing the parton list in a certain way, e.g. gluons first. Should be random which is which.
  if(ran->Get() < 0.5){
    partvec[0] = part_in_2;
    partvec[1] = part_in_1;
  }
  partpair partpair_temp;

  double shat = (partvec[0]->Momentum()+partvec[1]->Momentum()).Abs2();

  bool can_merge = m_splitting_merging;


  //Set starting tau for looking for closest approach
  pair<double, double> formation_tau_1 = part_in_1->GetFormationTau();
  pair<double, double> formation_tau_2 = part_in_2->GetFormationTau();
  
  double tau_start = tau;
  /* Temporary: hc while formation time*/
  if(formation_tau_1.second > tau && formation_tau_2.second > tau){
    if(formation_tau_1.second > formation_tau_2.second){
      tau_start = formation_tau_1.second;
    } else{
      tau_start = formation_tau_2.second;
    }
  }
  
  shared_ptr<Parton> part_merge_recoil;
  if(can_merge){
    part_merge_recoil = (p_dyn_quant_handler->FindRecoilPartner(tau_start, partvec[0], partvec[1]))[0].second;
    partvec.push_back(part_merge_recoil);
  }
  
  
  double T_star, mg2, mf2;
  if(m_splitting_merging){
    vector<double> temp_m2 = p_dyn_quant_handler->Getm2Tstar(partvec[0], partvec[1], tau_start);
    mg2 = temp_m2[0];
    mf2 = temp_m2[1];
    T_star = temp_m2[2];

    vector<double> mu2_temp;
    mu2_temp.push_back(mf2);
    mu2_temp.push_back(mg2);
    p_elastic_kernel->SetMu2Dyn(mu2_temp);
  }
  

  double xsec_scatter = 0.;
  double oe_1, oe_2, oe_3;
  Flavour flav_in_1 = partvec[0]->Flav();
  Flavour flav_in_2 = partvec[1]->Flav();
  //Temporary disabled
  if(false){
    //Setup overestimates for Bose/Pauli factors, if enabled
    if(m_include_bose_factors){
      //Temporary, if to be used setup dynamic OE for Bose factors
      oe_1 = 5.;
      oe_2 = 5.;
    } else{
      oe_1 = 1.;
      oe_2 = 1.;
    }

    //Find xsec for 2->2 scattering
    
    if(m_elastic_scattering){
      if(m_fixed_sigma.first){
        xsec_scatter = m_fixed_sigma.second;
      } else{
        xsec_scatter = p_elastic_kernel->FindScatteringXsec(flav_in_1,flav_in_2,shat); //NOTE: this should be an overestimate of xsec_scatter wrt. mg^2, using analytical value for now
      }
      xsec_scatter = m_OE_mult_scatter*xsec_scatter*oe_1*oe_2;
    }
  }


  //Find xsec for 2->1 merging
  double xsec_merge = 0.;
  double x_merge;
  
  double x;
  if(can_merge){
    Vec4D Pi = partvec[0]->Momentum();
    Vec4D Pj = partvec[1]->Momentum();
    Vec4D Pk = part_merge_recoil->Momentum();
    x = Pi*Pk/(Pi*Pk + Pj*Pk);
    double Q2 = (Pi + Pj + Pk)*(Pi + Pj + Pk);

    int N_kinematic_reject = 0;
    vector<pair<double, shared_ptr<Parton>>> include_map;
    pair<bool,vector<double>> merging_kinematics = FindMergingKinematics(partvec[0], partvec[1], part_merge_recoil);

    while(!merging_kinematics.first){
      N_kinematic_reject++;
      if(N_kinematic_reject < 2){
        include_map = p_dyn_quant_handler->FindRecoilPartner(tau_start, partvec[0], partvec[1]);
      }

      if(include_map.size() >= N_kinematic_reject+3){
        part_merge_recoil = include_map[N_kinematic_reject+2].second;
      } else{ //Not enough recoilers to pick from, merging false
        merging_kinematics.first = false;
        break;
      }

      if(part_merge_recoil == partvec[0] || part_merge_recoil == partvec[1]){
        //Not allowed as recoil, skip to next item in list
        continue;
      }

      merging_kinematics = FindMergingKinematics(partvec[0], partvec[1], part_merge_recoil);

      if(N_kinematic_reject > 20.){
        merging_kinematics.first = false;
        break;
      }
    }

    if(merging_kinematics.first){
      //Boost-no-boost
      double p_out = merging_kinematics.second[5];
      double kT2 = merging_kinematics.second[0];
      //double p_out = merging_kinematics.second[3]; //3 for boosted, 5 for non-boosted

      if(m_include_bose_factors){
        //Temporary, if to be used later include dynamic OE
        oe_3 = 5.;
      } else{
        oe_3 = 1.;
      }
      //Boost-no-boost
      xsec_merge = p_inelastic_kernel->FindMergingXsec(partvec[0], partvec[1], T_star, mg2, mf2, x, kT2, Pi[0], Pj[0], p_out, Q2*x*(1.-x));
      //xsec_merge = FindMergingXsec(partvec[0], partvec[1], T_star, mg2, mf2, x, kT2, merging_kinematics.second[1], merging_kinematics.second[3], p_out, Q2*x*(1.-x));
      //ThermalUO: add extra factors for oe
      xsec_merge = m_OE_mult_merge*xsec_merge*oe_3;
    } else{
      //msg_Out() << METHOD << ": WARNING: mergning kinematics rejected in finding tau scatter/merge, xsec_merge=0" << endl;
      xsec_merge = 0.;
    }
  }

  double taubar = partvec[0]->TauBar(partvec[1]);
  //msg_Out() << "#### Taubar found = " << taubar << endl;
  //msg_Out() << "Limit lower = " << tau_start << endl;
  //msg_Out() << "Limit upper 1 = " << p_tau_restart_list->first << endl;
  //msg_Out() << "Limit upper 2 = " << p_tau_restart_list->second << endl;

  double delta_tau = 0.;
  int process_type = 0; //Case of scattering
  double xsec_used;

  xsec_used = 1.;


  if(taubar > tau_start && taubar < p_tau_restart_list->first*(double(p_tau_restart_list->second) + 0.1)){
    /* Temporary disabled
    */
  } else{ //If t_split is later than m_tau_max, set a taubar that can never be reached
    taubar = (m_taumax_scaling_limit+1.)*m_taumax;
  }
  
  Vec4D shift(0.,0.,0.,0.);
  processpair processpair_temp = make_pair(make_pair(xsec_used, make_pair(process_type, delta_tau)),shift);
  taubarpair taubarpair_temp = make_pair(partvec, processpair_temp);
  retValue = make_pair(taubar, taubarpair_temp);

  if(taubar == tau){
    msg_Out() << METHOD <<  "Warning, taubar = m_tau = " << taubar << ", for partons " << partvec[0]->Number() << ", " << partvec[1]->Number() << endl;
  }
  //msg_Out() << "Sending again taubar = " << retValue.first << " with m_taumax = " << m_taumax <<"\n\n" << endl;
  return retValue;
}


pair<bool, vector<double>> Scatter_Merge_Handler::FindMergingKinematics(shared_ptr<Parton> part_in_1, shared_ptr<Parton> part_in_2, shared_ptr<Parton> part_recoil){
  /*
    Function to find the 4-momentum Pp of the outgoing parton in a 2->1 merging process.
    Return value is a pair consisting of a bool to indicate the calculation of Pp is possible, and if true then also the 4-momentum Pp.
  */

  bool merging_kinematics = true;
  Vec4D Pijt, Pkt;
  Vec4D Pk = part_recoil->Momentum(); 
  Vec4D Pi = part_in_1->Momentum();
  Vec4D Pj = part_in_2->Momentum();
  Vec4D Pk_save = Pk;
  Flavour flav_in_1 = part_in_1->Flav();
  Flavour flav_in_2 = part_in_2->Flav();
  vector<double> retVec = {0., 0., 0., 0., 0., 0., 0.};

  if(part_in_1->GetSplitMergeStatus() || part_in_2->GetSplitMergeStatus() || part_recoil->GetSplitMergeStatus()){
    msg_Out() << METHOD << ": WARNING: " << part_in_1->Number() << " " << part_in_2->Number() << " " << part_recoil->Number() << ", split status = " << part_in_1->GetSplitMergeStatus() << part_in_2->GetSplitMergeStatus() << part_recoil->GetSplitMergeStatus() << endl;
    //exit(1);
  }

  Flavour flav_out;
  if(flav_in_1.IsGluon() && flav_in_2.IsGluon()){ // gg -> g
    flav_out = kf_gluon;
  } else if( (flav_in_1.IsFermion() && flav_in_2.IsGluon()) || (flav_in_1.IsGluon() && flav_in_2.IsFermion()) ){ // gq -> q
    if(flav_in_1.IsFermion()){
      flav_out = flav_in_1;
    } else{
      flav_out = flav_in_2;
    }
  } else if( (flav_in_1.IsAnti() && !flav_in_2.IsAnti()) || (!flav_in_1.IsAnti() && flav_in_2.IsAnti()) ){ // qqb -> g
    flav_out = kf_gluon;
  }
  

  double x = Pi*Pk/(Pi*Pk + Pj*Pk);
  double Q2 = (Pi + Pj + Pk)*(Pi + Pj + Pk);

  //Calculate kinematics
  Vec4D cms(Pi+Pj+Pk);
  Poincare boost(cms);
  boost.Boost(Pi);
  boost.Boost(Pj);
  boost.Boost(Pk);
  Poincare rot(Pk,Vec4D(0.,0.,0.,1.));
  rot.Rotate(Pk);
  rot.Rotate(Pi);
  rot.Rotate(Pj);
  

  double yijk = Pi*Pj/(Pi*Pj + Pi*Pk + Pj*Pk);
  double kT2 = Q2*yijk*x*(1.-x);

  //First check to see if kinematics are allowed
  if(yijk >= 1 || yijk <= 0){
    msg_Out() <<" \n\n" << METHOD << ": ERROR: yijk = " << yijk << " >= 1. Q2 = " << Q2 << ", kT2 = " << kT2 << ", x = " << x << endl;
    merging_kinematics = false;
  }
  
  if(Pi[3] > 0 || Pj[3] > 0){
    merging_kinematics = false;
  }

  Pijt = Pi + Pj - yijk*Pk/(1.-yijk);
  Pkt = Pk/(1.-yijk);

  retVec[0] = kT2;
  retVec[1] = Pi[0];
  retVec[2] = Pj[0];
  retVec[3] = Pijt[0];
  retVec[4] = 2.*Pijt*Pkt*x*(1.-x);

  rot.RotateBack(Pijt);
  rot.RotateBack(Pkt);
  boost.BoostBack(Pijt);
  boost.BoostBack(Pkt);

  retVec[5] = Pijt[0];

  //Second check to see if kinematics are allowed
  if(merging_kinematics){
    double tol = 1.e-7;
    if(abs(Pijt.Abs2()) > tol || abs(Pkt.Abs2()) > tol){
      //msg_Out() << METHOD << ": ERROR, some parton not massless for parton " << part_in_1->Number() << " merging with parton " << part_in_2->Number() << " with recoil parton " << part_recoil->Number() << endl;
      //msg_Out() << METHOD << ": Pijt.abs2() = " << Pijt.Abs2() << ", Pkt.abs2() = " << Pkt.Abs2() << endl;
      merging_kinematics = false;
    }

    Vec4D p_compare = part_in_1->Momentum() + part_in_2->Momentum() + part_recoil->Momentum() - Pijt - Pkt;
    if(p_compare[0] > tol ||p_compare[1] > tol || p_compare[2] > tol || p_compare[3] > tol){
      //msg_Out() << METHOD << ": ERROR, four momentum not conserved for parton " << part_in_1->Number() << " merging with parton " << part_in_2->Number() << " with recoil parton " << part_recoil->Number() << endl;
      //msg_Out() << "p_in + p_recoil - Pi - Pj - Pk = [" << p_compare[0] << ", " << p_compare[1] << ", " << p_compare[2] << ", " << p_compare[3]  << "]" << endl;
      merging_kinematics = false;
    }

    if(Pijt[0] < m_p_min || Pkt[0] < m_p_min){
      //msg_Out() << METHOD << ": ERROR, some parton with negative energy for parton " << part_in_1->Number() << " merging with parton " << part_in_2->Number() << " with recoil parton " << part_recoil->Number() << endl;
      merging_kinematics = false;
    }
  }
  
  if(merging_kinematics){
    vec4dpair Pijt_Pkt = make_pair(Pijt,Pkt);
    kinematicspair temp_kinematics_pair = make_pair(Pijt_Pkt, make_pair(flav_out,flav_out));
    *p_save_merging_kinematics = make_pair(part_recoil, temp_kinematics_pair);
  }
  

  return make_pair(merging_kinematics, retVec);
}

pair<int, pair<double, double>> Scatter_Merge_Handler::DoesScatterOrMerge(const double taubar, taubarpair taubar_pair) {
  /*
    Function to check if the pair part1 and part2 either scatter elastically (2->2) or merge (2->1).
    Returns a int which indicates if any of the processes takes place (0 failed, 1 scatter, 2 merge)
    and a double with the accept/reject ratio

  */
  int retValue = false;
  double accept_reject = 0.;
  double formation_time = 0.;
  double xsec_oe = taubar_pair.second.first.first;
  int process_type = taubar_pair.second.first.second.first;
  shared_ptr<Parton> part_in_1 = taubar_pair.first[0];
  shared_ptr<Parton> part_in_2 = taubar_pair.first[1];
  double shat((part_in_1->Momentum()+part_in_2->Momentum()).Abs2());
  if (shat < 0.) return make_pair(retValue, make_pair(accept_reject, formation_time));

  Vec4D shift = taubar_pair.second.second;

  pair<double, double> formation_tau_1 = part_in_1->GetFormationTau();
  pair<double, double> formation_tau_2 = part_in_2->GetFormationTau();
  bool inside_formation_tau_1 = (taubar >= formation_tau_1.first) && (taubar <= formation_tau_1.second);
  bool inside_formation_tau_2 = (taubar >= formation_tau_2.first) && (taubar <= formation_tau_2.second);

  /* Temporary: hc while formation time*/
  if(!inside_formation_tau_1 && !inside_formation_tau_2){ //If current tau is within the formation tau for one of the particles, elastic scattering or merging is not allowed
    shared_ptr<Parton> part_recoil;    
    
    double T_star, mg2, mf2;
    
    vector<double> temp_m2 = p_dyn_quant_handler->Getm2Tstar(part_in_1, part_in_2, taubar);
    mg2 = temp_m2[0];
    mf2 = temp_m2[1];
    T_star = temp_m2[2];

    
    double xsec_scatter = 0.;
    double xsec_merge = 0.;


    if(m_elastic_scattering){ //Case of 2->2 scattering   
      Flavour flav_in_1 = part_in_1->Flav();
      Flavour flav_in_2 = part_in_2->Flav();

      vector<double> m2_temp;
      m2_temp.push_back(mf2);
      m2_temp.push_back(mg2);
      p_elastic_kernel->SetMu2Dyn(m2_temp);


      //Find 2->2 kinematics
      bool scatter_kinematics = FindScatteringKinematics(part_in_1, part_in_2, taubar, mg2, mf2);
      //msg_Out() << "scatter_kinematics = " << scatter_kinematics << endl;
      if(!scatter_kinematics){
        //msg_Out() << METHOD << ": Ending because scattering_kinematics = " << scatter_kinematics << " at m_tau = " << m_tau << endl;
        xsec_scatter = 0.;
      } else{
        //Include Bose/Pauli factors
        double bose_out_1, bose_out_2;
        if(m_include_bose_factors){
          if(p_save_scatter_kinematics->first.second.first.IsFermion()){
            double fp = p_dyn_quant_handler->GetPSD(part_in_1, part_in_2, taubar, p_save_scatter_kinematics->first.first.first, true, T_star);
            bose_out_1 = 1.-fp;
          } else{
            double fp = p_dyn_quant_handler->GetPSD(part_in_1, part_in_2, taubar, p_save_scatter_kinematics->first.first.first, false, T_star);
            bose_out_1 = 1.+fp;
          }
          if(p_save_scatter_kinematics->first.second.second.IsFermion()){
            double fk = p_dyn_quant_handler->GetPSD(part_in_2, part_in_1, taubar, p_save_scatter_kinematics->first.first.second, true, T_star);
            bose_out_2 = 1.-fk;
          } else{
            double fk = p_dyn_quant_handler->GetPSD(part_in_2, part_in_1, taubar, p_save_scatter_kinematics->first.first.second, false, T_star);
            bose_out_2 = 1.+fk;
          }

          if(bose_out_1 < 0){
            bose_out_1 = 0.;
          }
          if(bose_out_2 < 0){
            bose_out_2 = 0.;
          }
        } else{
          bose_out_1 = 1.;
          bose_out_2 = 1.;
        }

        if(m_fixed_sigma.first){
          xsec_scatter = m_fixed_sigma.second;
        } else{
          xsec_scatter = p_elastic_kernel->FindScatteringXsec(flav_in_1,flav_in_2,shat);
          //double anisotropic_factor = shat/(2*(part_in_1->Momentum())[0]*(part_in_2->Momentum())[0]);
          //xsec_scatter = anisotropic_factor*xsec_scatter;
        }

        xsec_scatter = xsec_scatter*bose_out_1*bose_out_2;
      }
    }
    
    if(m_splitting_merging){ //Case of 2->1 merging
      //Check if recoil partner has changed due to splitting or being recoil for another process
      //If so, find new recoil partner
      part_recoil = taubar_pair.first[2];

      //If recoil partner has already split, find new recoil partner.
      bool recoil_split = false;
      if(part_recoil->GetSplitMergeStatus()){
        part_recoil = (p_dyn_quant_handler->FindRecoilPartner(taubar, part_in_1, part_in_2))[0].second;
        recoil_split = true;
      }

      //Find 2->1 xsec
      Vec4D Pi = part_in_1->Momentum();
      Vec4D Pj = part_in_2->Momentum();
      Vec4D Pk = part_recoil->Momentum();
      double x = Pi*Pk/(Pi*Pk + Pj*Pk);

      if( !(x <= 0. || x >= 1.) ){
        //Find 2->1 kinematics
        int N_kinematic_reject = 0;
        vector<pair<double, shared_ptr<Parton>>> include_map;

        pair<bool,vector<double>> merging_kinematics = FindMergingKinematics(part_in_1, part_in_2, part_recoil);

        while(!merging_kinematics.first){
          N_kinematic_reject++;
          if(N_kinematic_reject < 2){
            include_map = p_dyn_quant_handler->FindRecoilPartner(taubar, part_in_1, part_in_2);
          }

          if(include_map.size() >= N_kinematic_reject+3){
            part_recoil = include_map[N_kinematic_reject+2].second;
          } else{ //Not enough recoilers to pick from, merging false
            merging_kinematics.first = false;
            break;
          }

          if(part_recoil == part_in_1 || part_recoil == part_in_2){
            //Not allowed as recoil, skip to next item in list
            continue;
          }

          merging_kinematics = FindMergingKinematics(part_in_1, part_in_2, part_recoil);

          if(N_kinematic_reject > 20){
            merging_kinematics.first = false;
            break;
          }
        }


        if(!merging_kinematics.first){
          //msg_Out() << METHOD << ": Ending because merging_kinematics = " << merging_kinematics << " at m_tau = " << m_tau << endl;
          xsec_merge = 0.;
        } else{
          //If kinematics allowed, evaluate xsec_merge
          double p_out_temp = Abs3(p_save_merging_kinematics->second.first.first);

          double kT2 = merging_kinematics.second[0];

          xsec_merge = p_inelastic_kernel->FindMergingXsec(part_in_1, part_in_2, T_star, mg2, mf2, x, kT2, Pi[0], Pj[0], p_out_temp, merging_kinematics.second[4]);

          //Include Bose/Pauli factors
          double bose_out_3, bose_out_3_oe;

          if(m_include_bose_factors){
            if(p_save_merging_kinematics->second.second.first.IsFermion()){
              double fp = p_dyn_quant_handler->GetPSD(part_in_1, part_in_2, taubar, p_save_merging_kinematics->second.first.first, true, T_star);
              bose_out_3 = 1.-fp;
            } else{
              double fp = p_dyn_quant_handler->GetPSD(part_in_1, part_in_2, taubar, p_save_merging_kinematics->second.first.first, false, T_star);
              bose_out_3 = 1.+fp;
            }
            if(bose_out_3 < 0){
              bose_out_3 = 0.;
            }
          } else{
            bose_out_3 = 1.;
          }

          xsec_merge = xsec_merge*bose_out_3;

          formation_time = p_inelastic_kernel->FindFormationTime(Pi[0]+Pj[0], T_star, mg2, -1.);
        }

      } else{ //Reject merger due to x being too small or large w.r.t. m_x_min
        xsec_merge = 0.;
      }
    }


    accept_reject = (xsec_scatter+xsec_merge)/xsec_oe;


    bool accept_process = false;
    double dij2 = part_in_1->Dist2(part_in_2, taubar);
    if(sqrt(dij2) <= sqrt((xsec_scatter+xsec_merge)/M_PI)){
      accept_process = true;
    }

    if(accept_process){
      //Accepted scatter/merge, decide which based on xsec ratio
      if(xsec_scatter == 0. && xsec_merge == 0.){
        //No scatter or merge
        retValue = 0;
      } else if(xsec_scatter == 0. && xsec_merge != 0.){
        //Merge
        retValue = 2;
      } else if(xsec_scatter != 0. && xsec_merge == 0.){
        //Scatter
        retValue = 1;
      } else{
        //Decide based on ratio
        if(ran->Get() <= xsec_scatter/(xsec_merge+xsec_scatter)){
          //Scatter
          retValue = 1;
        } else{
          //Merge
          retValue = 2;
        }
      }
    } else{
      retValue = 0;
    }


    return make_pair(retValue, make_pair(accept_reject, formation_time));
  } else{ //If taubar is within one particles formation time no scattering or merging is allowed
    retValue = 0;
    return make_pair(retValue, make_pair(accept_reject, formation_time));
  }
}


int Scatter_Merge_Handler::DoScatteringKinematics(shared_ptr<Parton> part1, shared_ptr<Parton> part2, const double taubar, Vec4D shift) {
  part1->SetLastScatter(part2);
  part2->SetLastScatter(part1);

  Vec4D x1, x2, x2_shifted;
  x1 = part1->Position(taubar);
  part1->SetPosition(x1);
  x2 = part2->Position(taubar);
  part2->SetPosition(x2);
  x2_shifted = x2 + shift;

  Vec4D p1(part1->Momentum()), p2(part2->Momentum());
  double shat = (p1+p2).Abs2();
  int oldindex(0),newindex(0);
  bool replace;
  if (shat < 0.) PRINT_VAR(shat); 
  if (shat < 0.) return 0;
  Vec4D p3 = p_save_scatter_kinematics->first.first.first;
  Vec4D p4 = p_save_scatter_kinematics->first.first.second;
  double that = (p1-p3).Abs2();
  Flavour in1 = part1->Flav();
  Flavour in2 = part2->Flav();
  Flavour out1 = p_save_scatter_kinematics->first.second.first;
  Flavour out2 = p_save_scatter_kinematics->first.second.second;
  double mg2 = p_save_scatter_kinematics->second.first;
  double mf2 = p_save_scatter_kinematics->second.second;

  //msg_Out() << "\nPre scatter: " << part1->Flav() << " " << part2->Flav() << " -> " << out1 << " " << out2 << endl;
  //msg_Out() << "Pre scatter:  [" << part1->GetFlow(1) << ", " << part1->GetFlow(2) << "] [" << part2->GetFlow(1) << ", " << part2->GetFlow(2) << "]" << endl;

  p_colour_handler->UpdateColoursScatter(part1, part2, out1, out2, that, mg2, mf2);

  //msg_Out() << "Post scatter:  [" << part1->GetFlow(1) << ", " << part1->GetFlow(2) << "] [" << part2->GetFlow(1) << ", " << part2->GetFlow(2) << "]" << endl;
  

  if(part1->GetFlow(1) == part1->GetFlow(2) || part2->GetFlow(1) == part2->GetFlow(2)){
    msg_Out() << METHOD << " WARNING: Post scatter:  [" << part1->GetFlow(1) << ", " << part1->GetFlow(2) << "] [" << part2->GetFlow(1) << ", " << part2->GetFlow(2) << "]" << endl;
  }
  // Histogram mean free path and interaction time
  double delta_t1 = x1[0] - part1->GetLastXPt();
  double delta_t2 = x2[0] - part2->GetLastXPt();

  part1->AddXPHistory(1, part1->Position(taubar), p3, out1, taubar);
  part2->AddXPHistory(1, part2->Position(taubar), p4, out2, taubar);

  part1->SetMomentum(p3);
  part1->SetFlav(out1);
  part2->SetMomentum(p4);
  part2->SetFlav(out2);
  part1->AddNScatter();
  part2->AddNScatter();

  if(m_timekeeper == 1){
    Vec4D a_hat(1.,0.,0.,0.);
    part1->SetTimekeeper(1, a_hat);
    part2->SetTimekeeper(1, a_hat);
  }
  
  part1->SetPosition(x1);
  part1->SetInittau(taubar);

  part2->SetPosition(x2);
  part2->SetInittau(taubar);

  return 1;
}


bool Scatter_Merge_Handler::FindScatteringKinematics(std::shared_ptr<Parton> part1, std::shared_ptr<Parton> part2, const double taubar, double mg2, double mf2){
  /*
    Function to find the effective mass of the system of based on partons part1 and part2, 
    for a fermion and a boson, summing others particles 3-momentum.

    Input:
      taubar      Tau at which mu2 is calculated (tau for closest approach of part1 and part2).
      part1       First particle in the collsion.
      part2       Second particle in the collsion.

    Output:
      val         pair with [{p3,p4},{that,phi}]
  */
  
  Vec4D p3(0.,0.,0.,0.), p4(0.,0.,0.,0.);

  Vec4D p1(part1->Momentum()), p2(part2->Momentum());
  Vec4D p1_save = p1;
  double E1 = p1[0];
  double E2 = p2[0];
  double shat((p1+p2).Abs2());
  if (shat < 0.) PRINT_VAR(shat); 
  if (shat < 0.) return false;
  Vec4D cms(p1+p2);
  Poincare boost(cms);
  boost.Boost(p1);
  boost.Boost(p2);
  Poincare rot(p1,Vec4D(0.,0.,0.,1.));
  rot.Rotate(p1);
  rot.Rotate(p2);

  double that,tmin,tmax;
  double pt, phi;
  
  tmax = p1.Abs2()-sqrt(shat)*(p1[0]-p1[3]);
  tmin = p1.Abs2()-sqrt(shat)*(p1[0]+p1[3]);

  if (tmin > tmax) return false;
  Flavour in1(part1->Flav()), in2(part2->Flav()), out1(kf_none), out2(kf_none);
  bool found_p = false;
  int N_iter = 0;
  while(!found_p){
    out1 = Flavour(kf_none);
    out2 = Flavour(kf_none);
    do{
      that = p_elastic_kernel->SelectThat(in1,in2,shat,out1,out2); //If fixed sigma, sets correct flavours and returns that=1
      if(that > 0. && m_fixed_sigma.first){
        that = -(tmax-tmin)*ran->Get();
      } else if (that>0.) {
        msg_Out() << METHOD << ": ERROR: that = " << that << " > 0 at taubar = " << taubar << endl;
        msg_Out() << METHOD << ": for << " << in1 << in2 << " to " << out1 << out2 << endl;
        msg_Out() << METHOD << ": particle number << " << part1->Number() << " and " << part2->Number() << endl;
        vector<double> m2_temp = p_elastic_kernel->GetMu2Dyn();
        msg_Out() << METHOD << ": m2 in xsec =  [" << m2_temp[0] << ", " << m2_temp[1] << "]\n" << endl;
        //m_tau = taubar;
        //part1->SetInittau(taubar);
        //part2->SetInittau(taubar);
        return false;
      }
    } while (that < tmin || that > tmax);
    
    pt = sqrt(shat/4. - sqr((that-p1.Abs2()+p1[0]*sqrt(shat))/(2.*p1[3])));
    phi = 2.*M_PI*ran->Get();
    p3 = Vec4D(sqrt(shat)/2.,pt*cos(phi),pt*sin(phi),sqrt(shat/4.-pt*pt));
    p4 = p1+p2-p3;

    if(p3[0] < 0. || p4[0] < 0.){
      msg_Out() << METHOD << ": ERROR: In COM: E_3 = " << p3[0] << ", E_4" << p4[0] << " for t_hat = " << that << ", that_min = " << tmin << endl;
      msg_Out() << METHOD << ": shat = " << shat << " ,shat/2 = " << shat/2. << ", shat/4 = " << shat/4. << endl;
    }

    rot.RotateBack(p3);
    rot.RotateBack(p4);
    boost.BoostBack(p3);
    boost.BoostBack(p4);

    if(p3[0] < m_p_min || p4[0] < m_p_min){
      //p_save_scatter_kinematics->first = make_pair(p3,p4);
      //p_save_scatter_kinematics->second = make_pair(out1,out2);
      //return false;
      //msg_Out() << METHOD << ": resampling scattering due to p3[0] = " << p3[0] << ", p4[0] = " << p4[0] << ", p_min = " << m_p_min << endl;
    } else{
      found_p = true;
    }

    N_iter++;

    //ThermalUO: increase from 20 to 50
    if(N_iter >= 50){
      msg_Out() << METHOD << ": rejecting scattering due to p3[0] = " << p3[0] << ", p4[0] = " << p4[0] << ", p_min = " << m_p_min << " after 50 iterations" << endl;
      p_save_scatter_kinematics->first.first = make_pair(p3,p4);
      p_save_scatter_kinematics->first.second = make_pair(out1,out2);
      p_save_scatter_kinematics->second = make_pair(mg2,mf2);
      return false;
    }
    //found_p = true;
  }
  /**/

  p_save_scatter_kinematics->first.first = make_pair(p3,p4);
  p_save_scatter_kinematics->first.second = make_pair(out1,out2);
  p_save_scatter_kinematics->second = make_pair(mg2,mf2);
  double delta_E;
  //msg_Out() << "dE_1 = " << abs(E1 - p3[0]) << ", dE_2 = " << abs(E2 - p4[0]) << ", (p_1-p_3)^2 = " << (p1_save-p3)*(p1_save-p3) << ", (p_1-4)^2 = " << (p1_save-p4)*(p1_save-p4) << ", t = " << that << endl;
  
  double tol = pow(10,-6);
  if(abs((p1_save-p3)*(p1_save-p3) - that) < tol){
    delta_E = abs(E1 - p3[0]);
  } else{
    delta_E = abs(E1 - p4[0]);
  }
  
  //m_save_elastic_that_dE = make_pair(that, delta_E);
  
  if(p3[0] < 0. || p4[0] < 0.){
    msg_Out() << METHOD << ": ERROR: E_3 = " << p3[0] << ", E_4" << p4[0] << " for t_hat = " << that << ", that_min = " << tmin << endl;
    msg_Out() << METHOD << ": shat = " << shat << " ,shat/2 = " << shat/2. << ", shat/4 = " << shat/4. << endl;
  }

  if(p3[0] < m_p_min || p4[0] < m_p_min){
    msg_Out() << METHOD << ": rejecting scattering due to p3[0] = " << p3[0] << ", p4[0] = " << p4[0] << ", p_min = " << m_p_min << endl;
    return false;
  } else{
    return true;
  }
}

std::pair<bool,std::shared_ptr<Parton>> Scatter_Merge_Handler::DoMergingKinematics(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double taumerge, double t_form){
  /*
    Function to update the kinematics for a 2->1 merging process after it has been accepted.
    Removes part_in_1 and part_in_2 from all relevant lists and creates a new parton with 4-momentum Pp.
    Returns a bool which indicates if the kinematics update was possible (and hence performed) or not,
    and a pointer to the new parton which is a result of the merging
  */
  shared_ptr<Parton> part_recoil = p_save_merging_kinematics->first;
  Vec4D Pp = p_save_merging_kinematics->second.first.first;
  Vec4D Rp = p_save_merging_kinematics->second.first.second;
  Vec4D x_temp = (part_in_1->Position(taumerge) + part_in_2->Position(taumerge))/2.;
  Vec4D x;

  if(ran->Get() <= 0.5){
    x = part_in_1->Position(taumerge);
  } else{
    x = part_in_2->Position(taumerge);
  }

  
  Vec4D temp_x1 = part_in_1->Position(taumerge);
  Vec4D temp_x2 = part_in_2->Position(taumerge);

  Flavour flav_out = p_save_merging_kinematics->second.second.first;

  //Bookkeep if the merging particles are previously split from the same particle
  double prev_split_partner = 0.;
  if(part_in_1->GetSplitPartner() == part_in_2->GetSplitPartner()){
    prev_split_partner = 1.;
  }

  double last_t_1 = part_in_1->GetLastXPt();
  double last_t_2 = part_in_2->GetLastXPt();
  double last_t = last_t_1;
  if(last_t_2 > last_t_1){
    last_t = last_t_2;
  }
  bool before_ini = false;
  if(last_t <= 0){
    before_ini = true;
  }

  double current_t = x[0];
  if(last_t < current_t || before_ini){
    if(abs(current_t - last_t) > t_form || before_ini){
      t_form = ran->Get()*t_form;
    } else{
      t_form = ran->Get()*(t_form - (current_t - last_t) );
    }
  }
  if(t_form < 0.){
    msg_Out() << "WARNING: t_form = " << t_form << endl;
    t_form = 0.;
  }

  // #### Create new particle with Pp ####
  //Create new partons
  Parton part_out(flav_out,Pp,x,taumerge);
  part_out.SetNumber();
  part_out.AddXPHistory(5, x,Pp,flav_out,taumerge);
  if(m_timekeeper == 1){
    Vec4D a_hat(1.,0.,0.,0.);
    part_out.SetTimekeeper(1, a_hat);
  }
  part_out.SetFormationTime(x[0], t_form);
  shared_ptr<Parton> part_out_ptr = make_shared<Parton>(part_out);
  part_out_ptr->SetLastScatter(part_recoil);

  //msg_Out() << "\nPre merge: " << part_in_1->Flav() << " " << part_in_2->Flav() << " -> " << part_out_ptr->Flav() << endl;
  //msg_Out() << "Pre merge:  [" << part_in_1->GetFlow(1) << ", " << part_in_1->GetFlow(2) << "] [" << part_in_2->GetFlow(1) << ", " << part_in_2->GetFlow(2) << "] -> [" << part_out_ptr->GetFlow(1) << ", " << part_out_ptr->GetFlow(2) << "]" << endl;

  p_colour_handler->UpdateColoursMerge(part_in_1, part_in_2, part_out_ptr);

  //msg_Out() << "Post merge: " << part_in_1->Number() << "  [" << part_in_1->GetFlow(1) << ", " << part_in_1->GetFlow(2) << "], " << part_in_2->Number() << " [" << part_in_2->GetFlow(1) << ", " << part_in_2->GetFlow(2) << "] -> " << part_out_ptr->Number() << " [" << part_out_ptr->GetFlow(1) << ", " << part_out_ptr->GetFlow(2) << "]" << endl;
  

  if(part_out_ptr->GetFlow(1) == part_out_ptr->GetFlow(2)){
    msg_Out() << METHOD << " WARNING: Post merge: " << part_in_1->Number() << "  [" << part_in_1->GetFlow(1) << ", " << part_in_1->GetFlow(2) << "], " << part_in_2->Number() << " [" << part_in_2->GetFlow(1) << ", " << part_in_2->GetFlow(2) << "] -> " << part_out_ptr->Number() << " [" << part_out_ptr->GetFlow(1) << ", " << part_out_ptr->GetFlow(2) << "]" << endl;
  }
  
  //Add new partons to the total list of partons (this also covers the taubar list when all lists are updated for the outgoing) 
  p_partons->push_back(part_out_ptr);

  //Add new partons to the tausplit list
  partpair temp_partpair = make_pair(part_out_ptr, part_out_ptr);
  kinematicspair_2 temp_kinematicspair;
  splitpair temp_splitpair = make_pair(temp_partpair, temp_kinematicspair);

  p_tausplits->push_back(make_pair((m_taumax_scaling_limit+1.)*m_taumax, temp_splitpair)); //Only relevant information is .second.first.first, i.e. part_out. All other information will be filled in later.

  // #### Update kinematics for recoil parton ####
  part_recoil->AddXPHistory(6, part_recoil->Position(taumerge),Rp,part_recoil->Flav(),taumerge);
  part_recoil->SetPosition(part_recoil->Position(taumerge));
  part_recoil->SetInittau(taumerge);
  part_recoil->SetMomentum(Rp);
  if(m_timekeeper == 1){
    Vec4D a_hat(1.,0.,0.,0.);
    part_recoil->SetTimekeeper(1, a_hat);
  }

  // #### Remove incoming partons ####

  //Set split status indicating that it has merged (needed if it is referenced as a recoil partner for cases of future splitting or merging even though it is gone from all the lists)
  part_in_1->SetSplitMergeStatus(true);
  part_in_1->AddXPHistory(7, part_in_1->Position(taumerge),part_in_1->Momentum(),part_in_1->Flav(),taumerge);
  part_in_2->SetSplitMergeStatus(true);
  part_in_2->AddXPHistory(7, part_in_2->Position(taumerge),part_in_2->Momentum(),part_in_2->Flav(),taumerge);

  p_removed_partons->push_back(part_in_1);
  p_removed_partons->push_back(part_in_2);


  //From list of all partons (pointers)
  for (list<shared_ptr<Parton>>::iterator mapit0 = p_partons->begin(); mapit0 != p_partons->end(); ) {
    if ((*mapit0)->Number() == part_in_1->Number() || (*mapit0)->Number() == part_in_2->Number()) {
      mapit0 = p_partons->erase(mapit0);  // erase returns the next valid iterator
    } else {
      ++mapit0;  // only increment if we didn't erase
    }
  }

  //From scattering (closest approach) list
  bool er = false;
  for (vector<pair<double,taubarpair>>::iterator mapit1=p_taubars->begin(); mapit1!=p_taubars->end(); mapit1++){
    if ((mapit1->second.first[0]->Number()==part_in_1->Number()) || (mapit1->second.first[1]->Number()==part_in_1->Number())) {
      er = true;
    }
    if ((mapit1->second.first[0]->Number()==part_in_2->Number()) || (mapit1->second.first[1]->Number()==part_in_2->Number())) {
      er = true;
    }
    if(er){
      p_taubars->erase(mapit1--);
      er = false;
    }
  }

  //From splitting list
  for (vector<std::pair<double,splitpair>>::iterator mapit2=p_tausplits->begin(); mapit2!=p_tausplits->end(); mapit2++){
    if (mapit2->second.first.first->Number()==part_in_1->Number() || mapit2->second.first.first->Number()==part_in_2->Number()){
      p_tausplits->erase(mapit2--);
    }
  }

  return make_pair(true, part_out_ptr);
}

double Scatter_Merge_Handler::Abs3(Vec4D v){
  /*
    Function to find the Euclidian spatial distance of a 4-vector

    Input:
      v     4D vector for which to find the absolute spatial distance from [0,0,0,0] for.
  */
  return sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}