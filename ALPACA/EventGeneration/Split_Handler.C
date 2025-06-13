#include "ALPACA/EventGeneration/Split_Handler.H"

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

using namespace ALPACA;

Split_Handler::Split_Handler(shared_ptr<list<shared_ptr<Parton>>> ptr_partons,
                             shared_ptr<list<shared_ptr<Parton>>> ptr_removed_partons,
                             shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler,
                             shared_ptr<Analysis_Handler> ptr_analysis_handler,
                             shared_ptr<Inelastic_Kernel> ptr_inelastic_kernel,
                             shared_ptr<pair<double, int>> ptr_tau_restart_list,
                             shared_ptr<pair<shared_ptr<Parton>, kinematicspair>>  ptr_save_split_kinematics,
                             shared_ptr<vector<pair<double,taubarpair>>> ptr_taubars,
                             shared_ptr<vector<pair<double,splitpair>>> ptr_tausplits,
                             shared_ptr<Colour_Handler> ptr_colour_handler):
  p_partons(ptr_partons), p_removed_partons(ptr_removed_partons), p_dyn_quant_handler(ptr_dyn_quant_handler), 
  p_analysis_handler(ptr_analysis_handler), p_inelastic_kernel(ptr_inelastic_kernel),
  p_colour_handler(ptr_colour_handler),
  p_taubars(ptr_taubars), p_tausplits(ptr_tausplits),
  p_tau_restart_list(ptr_tau_restart_list),
  m_splitting_merging(HIPars.SplittingMerging()), 
  m_taumax(HIPars.TauMax()), 
  m_only_gluons(HIPars.OnlyGluons()), m_p_min(HIPars.pMin()), 
  m_OE_mult_split(HIPars.OEMultSplit()),
  m_timekeeper(HIPars.Timekeeper()), m_include_bose_factors(HIPars.IncludeBoseFactors()), 
  m_taumax_scaling_limit(HIPars.TauMaxScalingLimit()),
  m_tsample_min(HIPars.tsampleMin()), m_tsample_max(HIPars.tsampleMax()),
  p_save_split_kinematics(ptr_save_split_kinematics)
{};


Split_Handler::~Split_Handler() {
};


pair<double,splitpair> Split_Handler::FindTauSplit(shared_ptr<Parton> part_in, double tau){
  /*
  Find tau (and recoil partner, branching energy fraction x, and kinematics) for part_in to split using the Sudakov Veto algorithm.
  This function only gives the t sampled from overestimate.
  */

  double t_split, tau_split;
  double gamma;
  double x;
  Flavour flavour_in = part_in->Flav();
  Flavour flavour_out_1, flavour_out_2;
  flavourpair flav_pair_out;
  shared_ptr<Parton> part_recoil, part_recoil_1, part_recoil_2;
  vector<Vec4D> Pp_Kp_Rp_save;
  vector<double> Tstar_mg2_mf2_x_kT2;
  double kT2, kT2_2;

  // #### Sudakov Veto algorithm ####
  double x_min = m_p_min/part_in->Abs3Momentum();
  double x_max = 1.-x_min;

  // No allowed splitting processes if there are only quarks in the system, if splitting is turned off or if x_min >= 0.5 meaning p is too low
  bool can_split = true;
  if(!m_splitting_merging || x_min >= 0.5){
    can_split = false;
  }

  if(can_split){
    Vec4D P_in = part_in->Momentum();
    double t_split_check = (part_in->Position(tau))[0];

    //Initialize variables for veto algorithm
    double R1, R2;
    vector<Vec4D> Pp_Kp_Rp, Pp_Kp_Rp_2;
    double bose_out_1, bose_out_2; // (1+-f) terms
    double C_g = 16.;
    double nu_s;
    if(flavour_in.IsQuark()){
      nu_s = 6.;
    } else{
      nu_s = 16.;
    }
    
    int N_count = 0; //Keep track of number of iterations in the veto algorithm
    int N_count_1 = 0; //Number of iterations for channel 1
    int N_count_2 = 0; //Number of iterations for channel 2
    int N_cutoff = 100000; //Cutoff for number of iterations if no solution is found in the veto algorithm
    int N_inc_limit;


    //Set initial parameter values for tau. Assumed not to change over these timescales.
    if(part_in->GetSplitRecoilPartner() == nullptr){
      part_recoil = (p_dyn_quant_handler->FindRecoilPartner(tau, part_in, part_in))[0].second;
    } else{
      part_recoil = part_in->GetSplitRecoilPartner();
    }
    part_recoil_1 = part_recoil;
    part_recoil_2 = part_recoil;

    double T_star, mg2, mf2;
    vector<double> temp_m2 = p_dyn_quant_handler->Getm2Tstar(part_in, part_in, tau);

    mg2 = temp_m2[0];
    mf2 = temp_m2[1];
    T_star = temp_m2[2];
    
    
    double Q2 = 2.*(part_in->Momentum())*(part_recoil->Momentum());
    double oe_1, oe_2; //Overestimates for Bose factors (1+-f) for main channel (q->gq or g->gg)
    if(m_include_bose_factors){
      double fp, fk;
      if(flavour_in.IsQuark()){
        fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tau, P_in, true, T_star);
        fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tau, P_in, false, T_star);
        oe_1 = 1.-fp;
        oe_2 = 1.+fk;
      } else{
        fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tau, P_in, false, T_star);
        fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tau, P_in, false, T_star);
        oe_1 = 1.+fp;
        oe_2 = 1.+fk;
      }
      if(oe_1 < 0){
        oe_1 = 0.;
      }
      if(oe_2 < 0){
        oe_2 = 0.;
      }
    } else{
      oe_1 = 1.;
      oe_2 = 1.;
    }

    double oe_1_2 = 0.;
    double oe_2_2 = 0.; //Overestimates for Bose factors (1+-f) for second channel if g->qqb possible/allowed
    if(flavour_in.IsGluon() && !m_only_gluons){
      if(m_include_bose_factors){
        double fp, fk;
        fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tau, P_in, true, T_star);
        fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tau, P_in, true, T_star);
        oe_1_2 = 1.-fp;
        oe_2_2 = 1.-fk;
        if(oe_1_2 < 0){
          oe_1_2 = 0.;
        }
        if(oe_2_2 < 0){
          oe_2_2 = 0.;
        }
      } else{
        oe_1_2 = 1.;
        oe_2_2 = 1.;
      }
    }

    double f; //Instantaneous splitting probability f(t,x) (=f(x)) from AMY
    t_split = (part_in->Position(tau))[0]; //Set initial t
    double t_max = (part_in->Position(m_taumax_scaling_limit*m_taumax))[0]; //Set max t

    pair<double,double> formation_time = part_in->GetFormationTime();
    double t_form_min = formation_time.first; 
    double t_form_max = formation_time.second;

    if(t_split >= t_form_min && t_split <= t_form_max){
      //If particle is recoil partner or new particle from splitting/merging it can already have formation time
      //if t_split initially is within formation time, jump to t_formation_max as start
      t_split = t_form_max;
    }

    double t_split_save = t_split; //Save previous t value in case draw needs to be redone

    vector<double> g; //Overestimate of gamma for veto algorithm, [0] for main channel, [1] for extra channel g->qqb if allowed

    //Boost-no-boost extra factor of *4
    if(!m_fixed_gamma.first){
      g = FindSplittingOE(part_in, T_star, mg2);
      g[0] = g[0]*oe_1*oe_2;
      g[1] = g[1]*oe_1_2*oe_2_2;
      //g[0] = 4.*g[0]*oe_1*oe_2;
      //g[1] = 4.*g[1]*oe_1_2*oe_2_2;
    } else{
      double g_fixed_1 = pow(2.*M_PI,3.)*m_fixed_gamma.second*oe_1*oe_2/(2.*nu_s*part_in->Abs3Momentum());
      double g_fixed_2 = 0.;

      if(flavour_in.IsGluon() && !m_only_gluons){
        g_fixed_2 = pow(2.*M_PI,3.)*m_fixed_gamma.second*oe_1_2*oe_2_2/(2.*nu_s*part_in->Abs3Momentum());
      }
      
      g = {g_fixed_1, 2.*3.*g_fixed_2}; //Overestimate with fixed gamma and a factor of 4 for the Bose terms, extra factor of 3*2 in g[1] for flavours and double counting g->qqb and g->qbq
      if(flavour_in.IsQuark()){
        g[0] = 2.*g[0]; //Needed since double counting q->gq and q->qg are distinct accounted for by 1/2 in cross section
      }
    }

    //ThermalUO: add extra factors for oe
    g[0] = m_OE_mult_split*g[0];
    g[1] = m_OE_mult_split*g[1];
    //Veto algorithm to draw tau
    bool cond_1 = false; //Separate conditions for when Veto Algorithm is done when multichannel is an option
    int N_kinematic_reject_1 = 0;
    // #### Sample first channel (q->gq or g->gg) and second if allowed (g->qqb) ####
    vector<pair<double, shared_ptr<Parton>>> include_map_1;
    while(!cond_1){ //Sample t_split
      if(t_split >= t_form_min && t_split <= t_form_max){
        //If particle is recoil partner or new particle from splitting it can already have formation time
        //if t_split ends up is within formation time, jump to t_formation_max as new start
        t_split = t_form_max;
        t_split_save = t_form_max;
      }

      R1 = ran->Get();
      R2 = ran->Get();

      t_split = t_split - log(R1)/((g[0]+g[1])*(x_max-x_min)); //Draw tau from g
      x = x_min + R2*(x_max - x_min); //Draw x from g

      //Find gamma to evaluate actual splitting probability for main channel
      gamma = 0.;

      //Find kinematics for sampled x to evaluate Bose terms
      pair<vector<Vec4D>, double> splitting_kinematics = FindSplittingKinematics(part_in, part_recoil_1, x, mg2, -1.);
      Pp_Kp_Rp = splitting_kinematics.first;
      kT2 = splitting_kinematics.second;

      if(kT2 < 0){ //Sampling of kinematics failed, ignore this draw, change recoil partner and redo
        t_split = t_split_save;
        N_kinematic_reject_1++;
        if(N_kinematic_reject_1 < 2){
          include_map_1 = p_dyn_quant_handler->FindRecoilPartner(tau, part_in, part_in);
        }

        bool found_recoil_part = false;

        while(!found_recoil_part){
          if(include_map_1.size() >= N_kinematic_reject_1+2){
            part_recoil_1 = include_map_1[N_kinematic_reject_1+1].second;
          } else{//No recoil parton found, list not long enough
            found_recoil_part = true;
            cond_1 = true;
            t_split = -1.;
          }

          if(part_recoil_1 == part_in){
            N_kinematic_reject_1++;
            continue;
          } else{
            found_recoil_part = true;
          }
        }
      } else{
        cond_1 = true;
        t_split_save = t_split;
      }

      N_count_1++;
      if(N_count_1 > N_cutoff){
        msg_Out() << "#### ERROR: N_count_1 = " << N_count_1 << " in Veto algorithm for main channel, ending sampling" << endl;
        cond_1 = true;
        t_split = -1.;
      }

      if(N_kinematic_reject_1 > 20 && !cond_1){
        //msg_Out() << "#### ERROR: N_kinematic_reject_1 = " << N_kinematic_reject_1 << " in Veto algorithm for main channel, ending sampling" << endl;
        cond_1 = true;
        t_split = -1.;
      }
    }

    //Evaluate if veto algorithm is done (and for which channel)
    bool cond = false; //Condition for if a t_split has been found
    if(t_split > 0){
      cond = true;
    }

    //If veto algorithm is done, set paramaters to save for output
    if(cond){
      Pp_Kp_Rp_save = Pp_Kp_Rp;
      Tstar_mg2_mf2_x_kT2 = {T_star, mg2, mf2, x, kT2, g[0]+g[1]};
      part_in->SetSplitRecoilPartner(part_recoil_1);
      part_recoil = part_recoil_1;

      tau_split = part_in->Tau(t_split);

      if(ran->Get() <= 0.5){
        flav_pair_out = make_pair(kf_gluon, flavour_in);
      } else{
        flav_pair_out = make_pair(flavour_in, kf_gluon);
      }

    } else{ //If no t has been found
      gamma = 0.;
      tau_split = (m_taumax_scaling_limit+1.)*m_taumax;
      x = 0.5;
    }

  } else{ //Case of no splittings allowed
    gamma = 0.;
    tau_split = (m_taumax_scaling_limit+1.)*m_taumax;
    x = 0.5;
  }
  
  
  partpair part_pair = make_pair(part_in, part_recoil);
  energypair energy_pair = make_pair(Pp_Kp_Rp_save, Tstar_mg2_mf2_x_kT2);
  kinematicspair_2 kinematics_pair = make_pair(energy_pair, flav_pair_out);
  splitpair split_pair = make_pair(part_pair, kinematics_pair);
  return make_pair(tau_split, split_pair);
}


vector<double> Split_Handler::FindSplittingOE(shared_ptr<Parton> part_in, double T_star, double mg2){
  /*
    Find the (constant )overestimate g of the instantaneous splitting probability f(t,x) = f(x) < g(t,x) = g
      g = P_1(p) * P_2(x_min) * mu2_LPM_oe(T_star,p,x_max,x_min)
  */
  vector<double> g = {0., 0.};
  Flavour flavour_in = part_in->Flav();
  double p = (part_in->Momentum())[0];
  double nu_A, d_A, C_A, nu_F, d_F, C_F;
  nu_F = 6.;
  d_F = 3.;
  C_F = 4./3.;
  nu_A = 16.;
  d_A = 8.;
  C_A = 3.;
  double xi = 9.09916;
  double N_c = 3.; //Temporary, just for gluons
  double lambda = 4.*M_PI*(p_dyn_quant_handler->GetAlphaS())*N_c;
  double x_min = m_p_min/p;
  double x_max = 1.-x_min;
  double C_sum, P_1, P_2;
  P_1 = sqrt(2.)*(p_dyn_quant_handler->GetAlphaS())*lambda*T_star*xi/(2.*pow(2.*M_PI,2.));

  if(flavour_in.IsQuark()){ //Case for just q->gq
    C_sum = (17./3.)/C_A;
    P_2 = C_sum*(d_F*C_F/nu_F)*(1+pow(x_max,2.))/x_min;
    g[0] = 2.*P_1*P_2; //*2 for double counting q->gq and q->qg (1/2 included in kernel for this)
  } else{                   //Case for g->gg
    C_sum = 9./C_A;
    P_2 = C_sum*(d_A*C_A/nu_A)*(1+pow(x_min,4.)+pow(x_max,4.))/(x_min*x_max);
    g[0] = P_1*P_2;
    if(!m_only_gluons){     //Case for g->qqb
      C_sum = (17./3)/C_A;
      P_2 = C_sum*(d_F*C_F/nu_F)*(x_min*x_min+x_max*x_max);
      g[1] = 2.*3.*P_1*P_2; //*2 for double counting g->qqb and g->qbq (1/2 included in kernel for this), *3 for flavours
    }
  }

  return g;
}

pair<vector<ATOOLS::Vec4D>, double> Split_Handler::FindSplittingKinematics(shared_ptr<Parton> part_in, shared_ptr<Parton> part_recoil, double x, double mg2, double kT2){
  /*
    Function to find the kinematics for a particle splitting ~collinearly.
    Only finds pp and kp, does not crate/delete particles and update lists.
    This is done to then evaluate probability of splitting into this pp and kp.
    Updates of list etc. is done in DoSplittingKinematics after the probability check.
  */
  vector<Vec4D> vec_Pi_Pj_Pk;
  Vec4D Pi, Pj, Pk;
  Vec4D Pijt = part_in->Momentum();
  Vec4D Pkt = part_recoil->Momentum();
  Vec4D Pijt_save = Pijt;
  Vec4D Pkt_save = Pkt;
  double Q2 = 2.*Pijt*Pkt;
  double kT2_min = 0.;
  double kT2_max = Q2*x*(1.-x);
  double kT;

  if(kT2 < 0){ //Set as input if it is the first calculation of kinematics to find tausplit
    kT = sqrt(p_inelastic_kernel->SamplekT2(kT2_min, kT2_max));
  } else{      //>0 if new kinematics has to be calculated at tausplit (due to e.g. elastic scattering of primary or recoil parton in the time leading up to tausplit)
    if(kT2 <= kT2_max && kT2 > kT2_min){ //If old kT is still allowed, keep that draw, otherwise redraw
      kT = sqrt(kT2);
    } else{
      kT = sqrt(p_inelastic_kernel->SamplekT2(kT2_min, kT2_max));
    } 
  }
  double phi = ran->Get()*2.*M_PI;
  Vec4D kT_vec(0., kT*cos(phi), kT*sin(phi), 0.);

  
  //Calculate kinematics, Pp and Kp, from z
  Vec4D cms(Pijt+Pkt);
  Poincare boost(cms);
  boost.Boost(Pijt);
  boost.Boost(Pkt);
  Poincare rot(Pijt,Vec4D(0.,0.,0.,1.));
  rot.Rotate(Pijt);
  rot.Rotate(Pkt);

  double yijk = kT*kT/(Q2*x*(1-x));

  Pi = x*Pijt + kT*kT*Pkt/(x*Q2) + kT_vec;
  Pj = (1.-x)*Pijt + kT*kT*Pkt/((1.-x)*Q2) - kT_vec;
  Pk = (1.-yijk)*Pkt;

  Vec4D Pijt_boost = Pijt;
  Vec4D Pi_boost = Pi;
  Vec4D Pj_boost = Pj;

  rot.RotateBack(Pi);
  rot.RotateBack(Pj);
  rot.RotateBack(Pk);
  boost.BoostBack(Pi);
  boost.BoostBack(Pj);
  boost.BoostBack(Pk);
  
  vec_Pi_Pj_Pk.push_back(Pi);
  vec_Pi_Pj_Pk.push_back(Pj);
  vec_Pi_Pj_Pk.push_back(Pk);
  vec_Pi_Pj_Pk.push_back(Pijt_save); //Also save original momentum of primary and recoil to check if it has changed at tausplit in DoSplittingKinematics()
  vec_Pi_Pj_Pk.push_back(Pkt_save);
  vec_Pi_Pj_Pk.push_back(Pijt_boost);
  vec_Pi_Pj_Pk.push_back(Pi_boost);
  vec_Pi_Pj_Pk.push_back(Pj_boost);

  kT2 = kT*kT;

  //Check to see if kinematics are allowed, otherwise return kT2 = -1
  double tol = 1.e-6;
  if(Pi.Abs2() > tol || Pj.Abs2() > tol ||  Pk.Abs2() > tol){
    msg_Out() << METHOD << ": ERROR, some parton not massless for parton " << part_in->Number() << " splitting with recoil parton " << part_recoil->Number() << endl;
    msg_Out() << METHOD << ": Pi.abs2() = " << Pi.Abs2() << ", Pj.abs2() = " << Pj.Abs2() << ", Pk.abs2() = " << Pk.Abs2() << endl;
    kT2 = -1.;
  }

  Vec4D p_compare = part_in->Momentum() + part_recoil->Momentum() - Pi - Pj - Pk;
  if(p_compare[0] > tol ||p_compare[1] > tol || p_compare[2] > tol || p_compare[3] > tol){
    msg_Out() << METHOD << ": ERROR, four momentum not conserved for parton " << part_in->Number() << " splitting with recoil parton " << part_recoil->Number() << endl;
    msg_Out() << "p_in + p_recoil - Pi - Pj - Pk = [" << p_compare[0] << ", " << p_compare[1] << ", " << p_compare[2] << ", " << p_compare[3]  << "]" << endl;
    kT2 = -1.;
  }

  if(Pi[0] < m_p_min || Pj[0] < m_p_min || Pk[0] < m_p_min){
    //msg_Out() << METHOD << ": ERROR, some parton with energy = " << Pi[0] << ", " << Pj[0] << ", " << Pk[0] << " < p_min for parton " << part_in->Number() << " splitting with recoil parton " << part_recoil->Number() << endl;
    kT2 = -1.;
  }


  return make_pair(vec_Pi_Pj_Pk, kT2);
}


pair<double, partpair> Split_Handler::DoSplittingKinematics(double tausplit, splitpair split_pair){
    /*
      Function to preform the kinematics for a particle splitting ~collinearly. 
        - Check if recoil partner still exists and if momentum has changed along the way
        - create the 2 new particles with pp and kp, and add them to all 3+1 lists (taubar, tausplit, taubound + p_partons)
        - remove the incoming particles from all 3+1 lists
        - update kinematics for recoil particle
    */

    shared_ptr<Parton> part_in = split_pair.first.first;
    shared_ptr<Parton> part_recoil = split_pair.first.second;
    vector<Vec4D> Pp_Kp_Rp = split_pair.second.first.first;
    vector<double> Tstar_mg2_mf2_x_kT2 = split_pair.second.first.second;
    flavourpair out_flavours = split_pair.second.second;


    Vec4D Pp = Pp_Kp_Rp[0];
    Vec4D Kp = Pp_Kp_Rp[1];
    Vec4D Rp = Pp_Kp_Rp[2];

    Vec4D x = part_in->Position(tausplit);
    Vec4D x2 = part_in->Position(tausplit);

    x[0] = x[0]-0.0000001;

    //Find formation time for the new pair
    double E_in = (part_in->Momentum())[0];
    double T_star = Tstar_mg2_mf2_x_kT2[0];
    double mg2 = Tstar_mg2_mf2_x_kT2[1];
    double kT2 = Tstar_mg2_mf2_x_kT2[4];
    double t_form;
    t_form = p_inelastic_kernel->FindFormationTime(E_in, T_star, mg2, -1.);

    double last_t = part_in->GetLastXPt();
    bool before_ini = false;
    if(last_t <= 0){
      before_ini = true;
    }
    double current_t = x[0];
    if(abs(current_t - last_t) > t_form || before_ini){
      t_form = ran->Get()*t_form;
    } else{
      t_form = ran->Get()*(t_form - (current_t - last_t) );
    }
    if(t_form < 0.){
      msg_Out() << "WARNING: t_form = " << t_form << endl;
      t_form = 0.;
    }

    // #### Create 2 new particles with Pp and Kp #### 

    //Create new partons
    Parton part_out_1(out_flavours.first,Pp,x,tausplit);
    Parton part_out_2(out_flavours.second,Kp,x,tausplit);
    part_out_1.SetNumber();
    part_out_1.AddXPHistory(2,x,Pp,out_flavours.first,tausplit);
    if(m_timekeeper == 1){
      Vec4D a_hat(1.,0.,0.,0.);
      part_out_1.SetTimekeeper(1, a_hat);
    }
    part_out_2.SetNumber();
    part_out_2.AddXPHistory(2,x,Kp,out_flavours.second,tausplit);
    if(m_timekeeper == 1){
      Vec4D a_hat(1.,0.,0.,0.);
      part_out_2.SetTimekeeper(1, a_hat);
    }
    part_out_1.SetFormationTime(x[0], t_form);
    part_out_2.SetFormationTime(x[0], t_form);
    shared_ptr<Parton> part_out_1_ptr = make_shared<Parton>(part_out_1);
    shared_ptr<Parton> part_out_2_ptr = make_shared<Parton>(part_out_2);
    part_out_1_ptr->SetSplitPartner(part_out_2_ptr);
    part_out_2_ptr->SetSplitPartner(part_out_1_ptr);

    //msg_Out() << "\nPre split: " << part_in->Flav() << " -> " << part_out_1_ptr->Flav() << " " << part_out_2_ptr->Flav() << endl;
    //msg_Out() << "Pre split:  [" << part_in->GetFlow(1) << ", " << part_in->GetFlow(2) << "] -> [" << part_out_1_ptr->GetFlow(1) << ", " << part_out_1_ptr->GetFlow(2) << "] [" << part_out_2_ptr->GetFlow(1) << ", " << part_out_2_ptr->GetFlow(2) << "]" << endl;

    p_colour_handler->UpdateColoursSplit(part_in, part_out_1_ptr, part_out_2_ptr);

    //msg_Out() << "Post split: " << part_in->Number() << " [" << part_in->GetFlow(1) << ", " << part_in->GetFlow(2) << "] -> " << part_out_1_ptr->Number() << " [" << part_out_1_ptr->GetFlow(1) << ", " << part_out_1_ptr->GetFlow(2) << "], " << part_out_2_ptr->Number() << " [" << part_out_2_ptr->GetFlow(1) << ", " << part_out_2_ptr->GetFlow(2) << "]" << endl;
    

    if(part_out_1_ptr->GetFlow(1) == part_out_1_ptr->GetFlow(2) || part_out_2_ptr->GetFlow(1) == part_out_2_ptr->GetFlow(2)){
      msg_Out() << METHOD << "WARNING: Post split: " << part_in->Number() << " [" << part_in->GetFlow(1) << ", " << part_in->GetFlow(2) << "] -> " << part_out_1_ptr->Number() << " [" << part_out_1_ptr->GetFlow(1) << ", " << part_out_1_ptr->GetFlow(2) << "], " << part_out_2_ptr->Number() << " [" << part_out_2_ptr->GetFlow(1) << ", " << part_out_2_ptr->GetFlow(2) << "]" << endl;
    }
    partpair new_partons = make_pair(part_out_1_ptr, part_out_2_ptr);

    //Add new partons to the total list of partons (this also covers the taubar list when all lists are updated for the outgoing) 
    p_partons->push_back(part_out_1_ptr);
    p_partons->push_back(part_out_2_ptr);

    //Add new partons to the tausplit list
    partpair partpair_1 = make_pair(part_out_1_ptr, part_out_1_ptr);
    partpair partpair_2 = make_pair(part_out_2_ptr, part_out_2_ptr);
    kinematicspair_2 temp_kinematicspair;
    splitpair splitpair_1 = make_pair(partpair_1, temp_kinematicspair);
    splitpair splitpair_2 = make_pair(partpair_2, temp_kinematicspair);

    p_tausplits->push_back(make_pair((m_taumax_scaling_limit+1.)*m_taumax, splitpair_1)); //Only relevant information is .second.first.first, i.e. part_in. All other information will be filled in later.
    p_tausplits->push_back(make_pair((m_taumax_scaling_limit+1.)*m_taumax, splitpair_2));

    // #### Update kinematics for recoil parton ####
    part_recoil->AddXPHistory(3, part_recoil->Position(tausplit),Rp,part_recoil->Flav(),tausplit);
    part_recoil->SetPosition(part_recoil->Position(tausplit));
    part_recoil->SetInittau(tausplit);
    part_recoil->SetMomentum(Rp);
    if(m_timekeeper == 1){
      Vec4D a_hat(1.,0.,0.,0.);
      part_recoil->SetTimekeeper(1, a_hat);
    }

    // #### Remove incoming parton ####

    //Set split status indicating that it has splitted (needed if it is referenced as a recoil partner for cases of future splitting)
    part_in->SetSplitMergeStatus(true);
    part_in->AddXPHistory(4, part_in->Position(tausplit),part_in->Momentum(),part_in->Flav(),tausplit);

    p_removed_partons->push_back(part_in);

    //From list of all partons (pointers)
    for (list<shared_ptr<Parton>>::iterator mapit0=p_partons->begin(); mapit0!=p_partons->end();){
      if ((*mapit0)->Number()==part_in->Number()) {
        mapit0 = p_partons->erase(mapit0);
      } else{
        ++mapit0;
      }
    }

    //From scattering (closest approach) list
    for (vector<pair<double,taubarpair>>::iterator mapit1=p_taubars->begin(); mapit1!=p_taubars->end(); mapit1++){
      if ((mapit1->second.first[0]->Number()==part_in->Number()) || (mapit1->second.first[1]->Number()==part_in->Number())) {
        p_taubars->erase(mapit1--);
      }
    }

    //From splitting list
    for (vector<pair<double,splitpair>>::iterator mapit2=p_tausplits->begin(); mapit2!=p_tausplits->end(); mapit2++){
      if (mapit2->second.first.first->Number()==part_in->Number()) {
        p_tausplits->erase(mapit2--);
      }
    }

  return make_pair(true, new_partons);
}

pair<bool, splitpair> Split_Handler::DoesSplit(double tausplit, splitpair split_pair){
  /*
    Evaluate if particle will split
  */

  shared_ptr<Parton> part_in = split_pair.first.first;
  shared_ptr<Parton> part_recoil = split_pair.first.second;
  vector<Vec4D> Pp_Kp_Rp = split_pair.second.first.first;
  vector<double> Tstar_mg2_mf2_x_kT2 = split_pair.second.first.second;
  flavourpair out_flavours = split_pair.second.second; //Not used
  Flavour flavour_in = part_in->Flav();
  Flavour flavour_out_1, flavour_out_2;
  if(flavour_in.IsFermion()){ //q->gq
    if(ran->Get() <= 0.5){
      flavour_out_1 = kf_gluon;
      flavour_out_2 = flavour_in;
    } else{
      flavour_out_1 = flavour_in;
      flavour_out_2 = kf_gluon;
    }
  } else{ //g->gg (g->qqb treated separately later for FindGamma())
    flavour_out_1 = kf_gluon;
    flavour_out_2 = kf_gluon;
  }
  
  double kT2 = Tstar_mg2_mf2_x_kT2[4];

  double T_star, mg2, mf2;
  vector<double> temp_m2 = p_dyn_quant_handler->Getm2Tstar(part_in, part_in, tausplit);

  mg2 = temp_m2[0];
  mf2 = temp_m2[1];
  T_star = temp_m2[2];


  bool allow_splitting = true;


  //If recoil partner has already split, find new recoil partner.
  bool change_in = false, change_recoil = false;
  Vec4D x = part_in->Position(tausplit);
  Vec4D x2 = part_in->Position(tausplit);
  
  bool new_recoil = false;
  if(part_recoil->GetSplitMergeStatus()){
    change_recoil = true;
    new_recoil = true;
  }

  //Check if momentum has changed since tausplit was calculated (due to elastic scatterings for either primary or recoil parton, or recoil acting as recoil for other split)
  //If so, recalculate new kinematics keeping the same k_T draw (unless it is now outside of phase space, then redraw k_T)
  partpair new_partons;
  Vec4D Pp, Kp, Rp, Pijt_boost, Pp_boost, Kp_boost;
  Vec4D Pin = part_in->Momentum();
  Vec4D Precoil = part_recoil->Momentum();
  Vec4D Pin_save = Pp_Kp_Rp[3];
  Vec4D Precoil_save = Pp_Kp_Rp[4];
  if(Pin[0] != Pin_save[0] || Pin[1] != Pin_save[1] || Pin[2] != Pin_save[2] || Pin[3] != Pin_save[3]){
    change_in = true;
    //msg_Out() << "Change in primary parton when doing splitting kinematics" << endl;
  }

  if(!change_recoil){
    if(Precoil[0] != Precoil_save[0] || Precoil[1] != Precoil_save[1] || Precoil[2] != Precoil_save[2] || Precoil[3] != Precoil_save[3]){
      //msg_Out() << "Change in recoil parton when doing splitting kinematics" << endl;
      change_recoil = true;
    }
  }

  if(change_in || change_recoil){
    pair<vector<Vec4D>,double> Pp_Kp_Rp_new;
    vector<pair<double, shared_ptr<Parton>>> include_map = p_dyn_quant_handler->FindRecoilPartner(tausplit, part_in, part_in);
    int N_kinematic_reject = 0;
    double kT2_new = -1;
    while(kT2_new < 0){
      if(include_map.size() >= N_kinematic_reject+1){
        part_recoil = include_map[N_kinematic_reject].second;
      } else{
        kT2_new = -1.;
        break;
      }

      if(part_recoil == part_in){
        //Not allowed as recoil, skip to next item in list
        N_kinematic_reject++;
        continue;
      }

      Pp_Kp_Rp_new = FindSplittingKinematics(part_in, part_recoil, Tstar_mg2_mf2_x_kT2[3], Tstar_mg2_mf2_x_kT2[1], Tstar_mg2_mf2_x_kT2[4]);
      kT2_new = Pp_Kp_Rp_new.second;
      N_kinematic_reject++;

      if(N_kinematic_reject > 20){
        kT2_new = -1.;
        break;
      }
    }

    if(kT2_new < 0){ //Redraw of kinematics failed, do not allow splitting
      //msg_Out() << METHOD << ": ERROR: redraw of kinematics failed, kT2_new = " << kT2_new << ", does not allow splitting" << endl;
      allow_splitting = false;
      part_in->SetSplitRecoilPartner(nullptr);
    } else if((Pp_Kp_Rp_new.first[0])[0] <= m_p_min || (Pp_Kp_Rp_new.first[1])[0] <= m_p_min){
      //msg_Out() << METHOD << ": ERROR: redraw of kinematics failed, E_p or E_k < m_p_min, does not allow splitting" << endl;
      allow_splitting = false;
      part_in->SetSplitRecoilPartner(nullptr);
    } else{
      Pp = Pp_Kp_Rp_new.first[0];
      Kp = Pp_Kp_Rp_new.first[1];
      Rp = Pp_Kp_Rp_new.first[2];
      Pijt_boost = Pp_Kp_Rp_new.first[5];
      Pp_boost = Pp_Kp_Rp_new.first[6];
      Kp_boost = Pp_Kp_Rp_new.first[7];
      kT2 = kT2_new;
      part_in->SetSplitRecoilPartner(part_recoil);
    }
  } else{
    Pp = Pp_Kp_Rp[0];
    Kp = Pp_Kp_Rp[1];
    Rp = Pp_Kp_Rp[2];
    Pijt_boost = Pp_Kp_Rp[5];
    Pp_boost = Pp_Kp_Rp[6];
    Kp_boost = Pp_Kp_Rp[7];
  }

  if(allow_splitting){
    //Accept/reject
    double bose_out_1, bose_out_2, nu_s;
    double bose_out_1_2 = 0.;
    double bose_out_2_2 = 0.;
    
    //boost-no-boost
    if(m_include_bose_factors){
      double fp, fk;
      if(flavour_in.IsQuark()){
        if(flavour_out_1.IsQuark()){
          fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Pp, true, T_star);
          fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Kp, false, T_star);
          bose_out_1 = 1.-fp;
          bose_out_2 = 1.+fk;
        } else{
          fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Pp, false, T_star);
          fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Kp, true, T_star);
          bose_out_1 = 1.+fp;
          bose_out_2 = 1.-fk;
        }
      } else{
        fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Pp, false, T_star);
        fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Kp, false, T_star);
        bose_out_1 = 1.+fp;
        bose_out_2 = 1.+fk;
      }

      if(!m_only_gluons){
        fp = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Pp, true, T_star);
        fk = p_dyn_quant_handler->GetPSD(part_in, part_in, tausplit, Kp, true, T_star);
        bose_out_1_2 = 1.-fp;
        bose_out_2_2 = 1.-fk;
      }
      
      if(bose_out_1 < 0){
        bose_out_1 = 0.;
      }
      if(bose_out_2 < 0){
        bose_out_2 = 0.;
      }
      if(bose_out_1_2 < 0){
        bose_out_1_2 = 0.;
      }
      if(bose_out_2_2 < 0){
        bose_out_2_2 = 0.;
      }

    } else{
      bose_out_1 = 1.;
      bose_out_2 = 1.;
      if(!m_only_gluons){
        bose_out_1_2 = 1.;
        bose_out_2_2 = 1.;
      }
    }


    if(flavour_in.IsFermion()){
      nu_s = 6;
    } else{
      nu_s = 16; 
    }

    double gamma;
    //Boost-no-boost
    double p_in = part_in->Abs3Momentum();

    if(m_fixed_gamma.first){
      gamma = m_fixed_gamma.second;
    } else{
      gamma = p_inelastic_kernel->FindGamma(p_in, flavour_in, flavour_out_1, flavour_out_2, Tstar_mg2_mf2_x_kT2[3], T_star, Tstar_mg2_mf2_x_kT2[1], Tstar_mg2_mf2_x_kT2[2]);
    }
    double f_1 =  pow(2.*M_PI,3.)*gamma*bose_out_1*bose_out_2/(2.*nu_s*p_in);
    if(flavour_in.IsFermion()){
      f_1 *= 2.; //*2 for double counting q->gq and q->qg (1/2 included in kernel for this)
    }
    double f = f_1;

    double f_2 = 0.;
    if(flavour_in.IsGluon() && !m_only_gluons){ //g->qqb
      flavour_out_1 = Flavour(kf_d);
      flavour_out_2 = Flavour(kf_d).Bar();
      if(m_fixed_gamma.first){
        gamma = m_fixed_gamma.second;
      } else{
        gamma = p_inelastic_kernel->FindGamma(p_in, flavour_in, flavour_out_1, flavour_out_2, Tstar_mg2_mf2_x_kT2[3], T_star, Tstar_mg2_mf2_x_kT2[1], Tstar_mg2_mf2_x_kT2[2]);
      }
      f_2 = pow(2.*M_PI,3.)*gamma*bose_out_1_2*bose_out_2_2/(2.*nu_s*p_in);
      f_2 *= 2.*3.; //*2 for double counting g->qqb and g>qbq (1/2 included in kernel for this), *3 for flavours
      f += f_2;
    }
    double f_oe = Tstar_mg2_mf2_x_kT2[5];

    double accept_reject = f/f_oe;

    if(accept_reject > 1.001){
      msg_Out() << METHOD << ": WARNING: accept_reject = " << accept_reject << " > 1 for splittings, T_star = " << T_star << " change_in, change_recoil, new_recoil = " << change_in << " " << change_recoil << " " << new_recoil << endl;
      p_analysis_handler->BookkeepOverdraw(accept_reject, true);
    }

    if(ran->Get() <= accept_reject){
      allow_splitting = true;

      flavourpair out_flavours;

      if(flavour_in.IsQuark()){ //Case for q->gq
        out_flavours = make_pair(flavour_out_1, flavour_out_2);
      } else{
        //Select g->gg or g->qqb (if allowed, otherwise f_2 = 0.)
        if(ran->Get() < f_1/f){ //Case for g->gg
          out_flavours = make_pair(kf_gluon, kf_gluon);
        } else{ //Case for g->qqb
          double R4 = ran->Get();
          if(R4 <= 1./3.){
            Flavour out_u(kf_u);
            out_flavours = make_pair(out_u, out_u.Bar());
          } else if (R4 <= 2./3.){
            Flavour out_d(kf_d);
            out_flavours = make_pair(out_d, out_d.Bar());
          } else{
            Flavour out_s(kf_s);
            out_flavours = make_pair(out_s, out_s.Bar());
          }
        }
      }

      split_pair.second.second = out_flavours;

      split_pair.first.second = part_recoil;
      
      Pp_Kp_Rp[0] = Pp;
      Pp_Kp_Rp[1] = Kp;
      Pp_Kp_Rp[2] = Rp;
      split_pair.second.first.first = Pp_Kp_Rp;

      Tstar_mg2_mf2_x_kT2[0] = T_star;
      Tstar_mg2_mf2_x_kT2[4] = kT2;
      split_pair.second.first.second = Tstar_mg2_mf2_x_kT2;

    } else{
      allow_splitting = false;
    }
  }
  return make_pair(allow_splitting, split_pair);
}

double Split_Handler::Abs3(Vec4D v){
  /*
    Function to find the Euclidian spatial distance of a 4-vector

    Input:
      v     4D vector for which to find the absolute spatial distance from [0,0,0,0] for.
  */
  return sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}