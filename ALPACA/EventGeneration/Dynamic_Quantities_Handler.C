#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
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

Dynamic_Quantities_Handler::Dynamic_Quantities_Handler(shared_ptr<list<shared_ptr<Parton>>> ptr_partons, 
                                                       shared_ptr<list<shared_ptr<Parton>>> ptr_removed_partons):
  p_partons(ptr_partons), p_removed_partons(ptr_removed_partons),
  m_N_inc(HIPars.NInclude()),
  m_alphaS(HIPars.AlphaS()),
  p_alphaS(static_cast<MODEL::Strong_Coupling *> (MODEL::s_model->GetScalarFunction(string("strong_cpl")))),
  m_timekeeper(HIPars.Timekeeper()), m_only_gluons(HIPars.OnlyGluons()),
  m_f_shell(HIPars.fShell()), m_f_r(HIPars.fr()), m_m2_min_scale(HIPars.M2MinScale()),
  m_f_delta_p(HIPars.fDeltap()), m_f_N_max(HIPars.fNMax()),
  p_analysis_handler(nullptr)
{};

Dynamic_Quantities_Handler::~Dynamic_Quantities_Handler() {
};

void Dynamic_Quantities_Handler::SetAnalysisHandler(shared_ptr<Analysis_Handler> ptr_analysis_handler) {
  p_analysis_handler = ptr_analysis_handler;
}

double Dynamic_Quantities_Handler::GetRunningAlphaS(double shat){
    return (*p_alphaS)(shat);
}

double Dynamic_Quantities_Handler::GetAlphaS(){
    return m_alphaS;
}

vector<double> Dynamic_Quantities_Handler::Getm2Tstar(shared_ptr<Parton> part1, shared_ptr<Parton> part2, double tau){
  /*
    Function to find m_g^2 and T_* at the average time and position of part1 and part2 at tau.
    Including m_N_inc partons in the calculations.
  */

  double mg2 = 0.;
  double mq2 = 0.;
  double T_star = 0.;

  double  alpha_s = GetAlphaS();
  double  C_A = 3.;
  double  C_F = 4./3.;
  double  d_A = 8.;
  double  sigma_px = 1./2.;
  double  sigma_x = sqrt(sigma_px);
  double  sigma_p = sqrt(sigma_px);
  
  sigma_p = sqrt(sigma_px);
  sigma_x = sqrt(sigma_px);
  double  Tsp_regular_term_gluon = 3.*M_PI*(GetAlphaS())*( 1+1./(16.*pow(2.*sigma_p*sigma_x,3)) )/4.;
  double  Tsp_regular_term_fermion = (8./48.)*(16./6.)*3.*M_PI*(GetAlphaS())*( 1+1./(6.*pow(2.*sigma_p*sigma_x,3)) )/4.;
  double  Tsp_ij_term_gluon = M_PI*alpha_s*C_A/(16.*pow(sigma_x*sigma_p,3.)*16.);
  double  Tsp_ij_term_fermion = M_PI*alpha_s*C_F/(16.*pow(sigma_x*sigma_p,3.)*36.);

  Vec4D     x1_ini = part1->Position(tau); 
  Vec4D     x2_ini = part2->Position(tau);
  Vec4D     x_ini = (x1_ini + x2_ini)/2.;
  double    t = x_ini[0];
  Flavour   flav_1 = part1->Flav();
  Flavour   flav_2 = part2->Flav();

  //Sort all particles w.r.t. distance to midpoint
  double dist;
  Vec4D x_temp, p_temp;
  vector<pair<double, pair<pair<Vec4D, Vec4D>, Flavour>>> rad_values, rad_values_2, rad_values_3, rad_values_Tstar;


  bool include_all = true;
  int N_existing = 0;
  Flavour flav_temp;

  if(m_timekeeper == 1){
    for(list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
      //if((*iter) == part1) continue;
      //if((*iter) == part2) continue;

      x_temp = (*iter)->Position(tau);
      p_temp = (*iter)->Momentum();
      flav_temp = (*iter)->Flav();
      dist = FindClosestR(x_ini, x_temp);
      rad_values.push_back(make_pair(dist, make_pair(make_pair(x_temp, p_temp), flav_temp)));
      N_existing++;
      if(!include_all){ //If m_N_inc < 0 use all partons
        if(N_existing >= m_N_inc+1) break; //If fermions are to be added then this break needs to be moved
      }
    }

  } else{
    pair<bool, pair<pair<Vec4D, Vec4D>, Flavour>> xp_history;

    list<shared_ptr<Parton>> all_partons;
    all_partons = *p_partons;
    if(p_removed_partons->size() > 0){
      all_partons.insert(all_partons.end(), p_removed_partons->begin(), p_removed_partons->end() );
    }
    
    for(list<shared_ptr<Parton>>::iterator iter=all_partons.begin(); iter!=all_partons.end(); iter++) {
      if((*iter) == part1) continue;
      if((*iter) == part2) continue;
      xp_history = (*iter)->GetXP(t);
      if(xp_history.first){ //If the parton exists at time t
        x_temp = xp_history.second.first.first;
        //msg_Out() << "x_temp = " << x_temp[0] << ", " << x_temp[1] << ", " << x_temp[2] << ", " << x_temp[3] << endl;
        p_temp = xp_history.second.first.second;
        
        flav_temp = xp_history.second.second;
        dist = FindClosestR(x_ini, x_temp);
        //msg_Out() << "t = " << t << ", dist = " << dist << ", p_temp[0] = " << p_temp[0] << ", x_temp[0] = " << x_temp[0] << endl;
        rad_values.push_back(make_pair(dist, make_pair(make_pair(x_temp, p_temp), flav_temp)));
        N_existing++;
        if(!include_all){ //If m_N_inc < 0 use all partons
          if(N_existing >= m_N_inc+1) break;
        }
      }
    }
  }

  if(rad_values.size() > 0){
    sort(rad_values.begin(), rad_values.end(), sortM2);
  }

  rad_values_2 = rad_values;

  Vec4D cms(0.,0.,0.,0.);

  int N_com_added = 0;
  for(vector<pair<double, pair<pair<Vec4D, Vec4D>, Flavour>>>::iterator iter=rad_values_2.begin(); iter!=rad_values_2.end(); iter++){
    p_temp = (*iter).second.first.second;
    cms = cms + p_temp;
    N_com_added++;
    if(!include_all){ //If m_N_inc < 0 use all partons
      if(N_com_added >= m_N_inc+1) break;
    }
  }

  //If everything should be calculated in COM
  Poincare boost(cms);
  boost.Boost(x_ini);
  N_com_added = 0;
  for(vector<pair<double, pair<pair<Vec4D, Vec4D>, Flavour>>>::iterator iter=rad_values_2.begin(); iter!=rad_values_2.end(); iter++){
    x_temp = (*iter).second.first.first;
    p_temp = (*iter).second.first.second;
    boost.Boost(x_temp);
    boost.Boost(p_temp);
    dist = FindClosestR(x_ini, x_temp);
    (*iter).second.first.first = x_temp;
    (*iter).second.first.second = p_temp;
    (*iter).first = dist;
    rad_values_3.push_back(make_pair(dist, make_pair(make_pair(x_temp, p_temp), (*iter).second.second)));
    N_com_added++;
    if(!include_all){ //If m_N_inc < 0 use all partons
      if(N_com_added >= m_N_inc+1) break;
    }
  }

  if(rad_values_3.size() > 0){
    sort(rad_values_3.begin(), rad_values_3.end(), sortM2);
  }

  //If not enough values found, add values far away
  if(rad_values.size() < m_N_inc){
    msg_Out() << METHOD << ": WARNING: not enough particles found in GetM2" << endl;
    Vec4D temp_vec(0.,0.,0.,0.);
    for(int j = 0; j <= m_N_inc; j++){
      rad_values.push_back(make_pair(pow(10.,10.), make_pair(make_pair(temp_vec,temp_vec),flav_1))); 
    }
  }

  //If not enough values found, add values far away
  if(rad_values_3.size() < m_N_inc){
    msg_Out() << METHOD << ": WARNING: not enough particles found in GetM2 at tau = " << tau << ", p_partons->size() = " << p_partons->size() << endl;
    Vec4D temp_vec(0.,0.,0.,0.);
    for(int j = 0; j <= m_N_inc; j++){
      rad_values_2.push_back(make_pair(pow(10.,10.), make_pair(make_pair(temp_vec,temp_vec),flav_1))); 
    }
  }

  if(N_existing < m_N_inc){
    m_N_inc = N_existing-1.;
  }

  //Loop through closest particles and add to sum
  double mg2_v2 = 0.;
  int N_added = 0;
  double nu_g = 16.;
  double nu_q = 18.;
  Flavour flav;
  Vec4D p, x, x_rel;
  vector<pair<pair<Vec4D,Vec4D>,Flavour>> ij_vectors;

  double fg_bar;

  double R_max_1 = 0.;
  double R_max_2 = 0.;

  double R_temp;

  // #### Find m^2 ####
  for(vector<pair<double, pair<pair<Vec4D, Vec4D>, Flavour>>>::iterator iter=rad_values.begin(); iter!=rad_values.end(); iter++){
    x = (*iter).second.first.first;
    p = (*iter).second.first.second;
    flav = (*iter).second.second;
    x_rel = x_ini - x;
    //Add regular terms
    if(flav.IsFermion()){
      mg2 += alpha_s*C_F*M_PI/(2.*Abs3(p));
      mq2 += 16.*M_PI*alpha_s/(3.*nu_q*Abs3(p));
    } else{
      mg2 += alpha_s*C_A*M_PI/(2.*Abs3(p));
      mq2 += 32.*M_PI*alpha_s/(3.*nu_g*Abs3(p));
    }

    N_added += 1;

    R_temp = (*iter).first;
    if(R_temp > R_max_1) R_max_1 = R_temp;

    if(!include_all){ //If m_N_inc < 0 use all partons
      if(N_added >= m_N_inc) break;
    }

  }


  N_added = 0;
  double histo_ij_terms = 0.;
  rad_values_Tstar = rad_values_3;
  for(vector<pair<double, pair<pair<Vec4D, Vec4D>, Flavour>>>::iterator iter=rad_values_Tstar.begin(); iter!=rad_values_Tstar.end(); iter++){
    x = (*iter).second.first.first;
    p = (*iter).second.first.second;
    flav = (*iter).second.second;
    //Add regular terms
    if(flav.IsFermion()){
      T_star += M_PI*alpha_s*C_F*(1. + 1./(6.*pow(2.*sigma_x*sigma_p, 3.)))/4.;
      mg2_v2 += alpha_s*C_F*M_PI/(2.*Abs3(p));
    } else{
      // Original method
      T_star += M_PI*alpha_s*C_A*(1. + 1./(16.*pow(2.*sigma_x*sigma_p, 3.)))/4.;
      
      /*
      //New method using GetPSD
      Parton parton_temp_g(Flavour(kf_gluon),Vec4D(0., 0., 0., 0.),Vec4D(t, x[1], x[2], x[3]),tau);
      fg_bar = GetPSD(make_shared<Parton>(parton_temp_g), make_shared<Parton>(parton_temp_g), tau, p, false, 1000.);
      T_star += 3.*M_PI*alpha_s*(1.+fg_bar)/4.;
      */

      mg2_v2 += alpha_s*C_A*M_PI/(2.*Abs3(p));
    }


    R_temp = (*iter).first;
    if(R_temp > R_max_2) R_max_2 = R_temp;

    
    /* Original method */
    //For Tstar: Add mixed term ij, (i=/=j)
    double T_star_add, const_temp;
    for(vector<pair<pair<Vec4D,Vec4D>,Flavour>>::iterator iter3=ij_vectors.begin(); iter3!=ij_vectors.end(); iter3++){
      //if(flav == iter3->second){ //Only same flavour mixed terms contributes
      //Test boost to COM
      if((flav.IsGluon() && iter3->second.IsGluon()) || (flav.IsQuark() && iter3->second.IsQuark())){
        if(flav.IsQuark()){
          T_star += Tsp_ij_term_fermion*exp(-pow(Abs3(x-iter3->first.first),2)/(4.*sigma_x*sigma_x))*exp(-pow(Abs3(p-iter3->first.second),2)/(4.*sigma_p*sigma_p));
        } else{
          //const_temp = Abs3(x-iter3->first.first)/Abs3(p-iter3->first.second);
          //sigma_x = sqrt(sigma_px*const_temp);
          //sigma_p = sqrt(sigma_px/const_temp);
          T_star_add = Tsp_ij_term_gluon*exp(-pow(Abs3(x-iter3->first.first),2)/(4.*sigma_x*sigma_x))*exp(-pow(Abs3(p-iter3->first.second),2)/(4.*sigma_p*sigma_p));
          T_star += T_star_add;
          histo_ij_terms += T_star_add;
        }
      }
    }
    //sigma_x = sqrt(sigma_px);
    //sigma_p = sqrt(sigma_px);
    ij_vectors.push_back(make_pair(make_pair(x,p),flav));
    
    
    N_added += 1;
    if(!include_all){ //If m_N_inc < 0 use all partons
      if(N_added >= m_N_inc) break;
    }
  }

  //Normalize
  double V_1, V_2;

  V_1 = 4.*M_PI*pow(R_max_1, 3.)/3.;
  V_2 = 4.*M_PI*pow(R_max_2, 3.)/3.;
  //V_1 = 4.*M_PI*pow((rad_values[m_N_inc].first+rad_values[m_N_inc-1].first)/2., 3.)/3.;
  //V_2 = 4.*M_PI*pow((rad_values_Tstar[m_N_inc].first+rad_values_Tstar[m_N_inc-1].first)/2., 3.)/3.; //W.r.t. m_N_inc

  double Lambda_QCD = 0.2;

  


  if(V_1 > 0. && V_2 > 0.){
    mg2 = mg2/V_1;
    mq2 = mq2/V_1;
    mg2_v2 = mg2_v2/V_2;

    if(mg2 <= m_m2_min_scale*Lambda_QCD*Lambda_QCD  || isinf(mg2) || isnan(mg2)){
      //msg_Out() << METHOD << ": WARNING: mg2 = " << mg2 << ", setting to " << m_m2_min_scale*Lambda_QCD*Lambda_QCD << endl;
      mg2 = m_m2_min_scale*Lambda_QCD*Lambda_QCD;
    }

    if(mq2 <= m_m2_min_scale*Lambda_QCD*Lambda_QCD  || isinf(mq2) || isnan(mq2)){
      //msg_Out() << METHOD << ": WARNING: mq2 = " << mq2 << ", setting to " << m_m2_min_scale*Lambda_QCD*Lambda_QCD << endl;
      mq2 = m_m2_min_scale*Lambda_QCD*Lambda_QCD;
    }
  
    if(mg2_v2 <= m_m2_min_scale*Lambda_QCD*Lambda_QCD  || isinf(mg2_v2) || isnan(mg2_v2)){
      //msg_Out() << METHOD << ": WARNING: mg2_v2 = " << mg2_v2 << ", setting to " << m_m2_min_scale*Lambda_QCD*Lambda_QCD << endl;
      mg2_v2 = m_m2_min_scale*Lambda_QCD*Lambda_QCD;
    }

    T_star = T_star/(V_2*mg2_v2);
  } else{
    mg2 = m_m2_min_scale*Lambda_QCD*Lambda_QCD;
    mq2 = m_m2_min_scale*Lambda_QCD*Lambda_QCD;
    T_star = 1.;
  }


  if(isinf(T_star) || isnan(T_star) || T_star < 0.000001){
    msg_Out() << "Warning: T_star = " << T_star << ", mg2 = " << mg2 << ", mq2 = " << mq2 << "mg2_v2 = " << mg2_v2 << endl;
    msg_Out() << part1 << ", " << part2 << ", tau = " << tau << endl;
    msg_Out() << "V_1 = " << V_1 << ", V_2 = " << V_2 << ", R_max_1 = " << R_max_1 << " - " << R_max_2 << endl;
    T_star = 1.; 
  }

  p_analysis_handler->Histo2D("t+mg2", t, mg2);
  p_analysis_handler->Histo2D("t+mq2", t, mq2);
  p_analysis_handler->Histo2D("t+Tstar", t, T_star);
    
  vector<double> retVal = {mg2, mq2, T_star};

  return retVal;
}


double Dynamic_Quantities_Handler::GetPSD(shared_ptr<Parton> part1, shared_ptr<Parton> part2, double tau, Vec4D p, bool isQuark, double Tstar){
  /*
    Function to find the PSD values f(p') and f(k') for outgoing momenta for xsec evaluation.
    Uses a fixed number of particles dN to determine dV, and then 
      f = (2\pi)^3*dN/(C_deg*dV)
    
    NOTE: should not rely on part2
  */

  double f_p;

  Vec4D x = part1->Position(tau);
  double t = x[0];

  double dist;
  Vec4D deltaX, deltaP, deltaP_2, x_temp_ini, x_temp, p_temp_ini, p_temp;
  bool isFerm = part1->Flav().IsFermion();
  Flavour flav = part1->Flav();
  vector<pair<double, pair<double, double>>> rad_values, rad_values_2;
  double absDeltaP, absDeltaX, absDeltaX2, absDeltaXP;

  double delta_p = m_f_delta_p;
  double delta_p_original = m_f_delta_p;
  double Vp, r_1, r_2;
  
  if(m_f_shell){
    r_1 = p[0] - delta_p/2;
    r_2 = p[0] + delta_p/2;
    if(r_1 < 0.){
      r_1 = 0.;
    }
    Vp = 4.*M_PI*(r_2*r_2*r_2 - r_1*r_1*r_1)/3.;
  } else{
    Vp = 4.*M_PI*pow(delta_p, 3.)/3.;
  }

  double Vx;
  double r_max = m_f_r;

  pair<bool, pair<pair<Vec4D, Vec4D>, Flavour>> xp_history;

  list<shared_ptr<Parton>> all_partons;
  all_partons = *p_partons;
  if(p_removed_partons->size() > 0){
    all_partons.insert(all_partons.end(), p_removed_partons->begin(), p_removed_partons->end() );
  }

  Flavour temp_flav;
  double invariant_mass, invariant_dist;
  int N = 0;
  bool found_1 = false;
  double N_added_extra = 1.;

  double R_max_inc = 0.;

  double R_temp;


  while(!found_1){
    for (list<shared_ptr<Parton>>::iterator iter=all_partons.begin(); iter!=all_partons.end(); iter++) {
      if((*iter) == part1) continue;
      if((*iter) == part2) continue;
      xp_history = (*iter)->GetXP(t);
      if(xp_history.first){ //If the parton exists at time t
        temp_flav = xp_history.second.second;
        if((!isQuark && temp_flav.IsGluon()) || (isQuark && temp_flav.IsFermion())){
          x_temp = xp_history.second.first.first;
          x_temp_ini = x;
          p_temp = xp_history.second.first.second;
          p_temp_ini = p;

          deltaX = x_temp_ini - x_temp;
          deltaP = p_temp_ini - p_temp;
          absDeltaP = Abs3(deltaP);
          absDeltaX = FindClosestR(x, x_temp); //include all

          R_temp = absDeltaX;
     
          if(absDeltaX < r_max || r_max <= 0){
            if(m_f_shell){
              if(p_temp[0] >= p[0] - delta_p/2. && p_temp[0] <= p[0] + delta_p/2. ){
                N += 1;
                if(R_temp > R_max_inc) R_max_inc = R_temp;
              }
            } else{
              if(absDeltaP <= delta_p){
                N += 1;
                if(R_temp > R_max_inc) R_max_inc = R_temp;
              }
            }
          }

        }
      }
    }

    if(N > 0){
      found_1 = true;
    } else{
      N = 0.;
      N_added_extra = N_added_extra + 1.;
      delta_p = delta_p_original*N_added_extra;
      if(m_f_shell){
        r_1 = p[0] - delta_p/2;
        r_2 = p[0] + delta_p/2;
        if(r_1 < 0.){
          r_1 = 0.;
        }
        Vp = 4.*M_PI*(r_2*r_2*r_2 - r_1*r_1*r_1)/3.;
      } else{
        Vp = 4.*M_PI*pow(delta_p, 3.)/3.;
      }
    }
    
    if(N_added_extra > m_f_N_max){
      found_1 = true;
    }
  }

  double nu_s;
  if(isQuark == true){
    nu_s = 6.*6.;
  } else{
    nu_s = 16.;
  }

  double p2 = Abs3(p)*Abs3(p);

  if(r_max <= 0){
    Vx = 4.*M_PI*pow(R_max_inc, 3.)/3.;
  } else{
    Vx = 4.*M_PI*pow(r_max, 3.)/3.;
  }

  double vol = Vx*Vp;
  if(vol > 0.){
    f_p = pow(2.*M_PI,3.)*double(N)/(nu_s*vol);
  } else{
    f_p = 0.;
  }

  if(isnan(f_p) || isinf(f_p) || f_p < 0.){
    msg_Out() << METHOD << ": WARNING: f_p = " << f_p << ", setting it to 0" << endl;
  }
  
  if(f_p > Tstar/p[0]){
    f_p = Tstar/p[0];
  }

  return f_p;
}


vector<pair<double, shared_ptr<Parton>>> Dynamic_Quantities_Handler::FindRecoilPartner(double tau, shared_ptr<Parton> part_1, shared_ptr<Parton> part_2){
  /*
    Function  used to pick recoil partner, returns vector of pair<distance, recoil parton> ordered by distance to midpoint of part_1 and part_2 at t_coll.
    For splitting, set part_1 = part_2
  */

  Vec4D   x_ini_1 = part_1->Position(tau);
  Vec4D   x_ini_2 = part_2->Position(tau);
  Vec4D   x_coll = (x_ini_1 + x_ini_2)/2.;
  double  t_coll = x_coll[0];

 Vec4D x_temp;
 double dist;
  vector<pair<double, shared_ptr<Parton>>> dist_list;
  for(list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
    if((*iter) == part_1 || (*iter) == part_2) continue;

    x_temp = (*iter)->PositionFromt(t_coll);
    dist = FindClosestR(x_coll, x_temp);
    
    dist_list.push_back(make_pair(dist, *iter));
  } 
  
  sort(dist_list.begin(), dist_list.end(), sortDistList);

  if(dist_list.size() == 0){
    msg_Out() << METHOD << ": WARNING: no recoil partner found" << endl;
    dist_list.push_back(make_pair(0., nullptr));
  }

  return dist_list;
}

double Dynamic_Quantities_Handler::FindClosestR(Vec4D x1, Vec4D x2){
  double R = Abs3(x1 - x2);
  return R;
}

double Dynamic_Quantities_Handler::Abs3(Vec4D v){
  /*
    Function to find the Euclidian spatial distance of a 4-vector

    Input:
      v     4D vector for which to find the absolute spatial distance from [0,0,0,0] for.
  */
  return sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}