#include "YFS/Main/Define_Dipoles.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__LOOPTOOLS
  #include "clooptools.h"
#endif

using namespace YFS;
using namespace ATOOLS;
using namespace std;

// #define C0(A,B,C,D,E,F) Master_Triangle(A,B,C,D,E,F)


Define_Dipoles::Define_Dipoles() {
  m_in = 2; // This is fine in YFS. It will not work for any other inital state multiplicity
  m_softphotonSum = {0., 0., 0., 0.};
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };;
  s.DeclareVectorSettingsWithEmptyDefault({"DIPOLE_FLAV"});
  m_N = s["DIPOLE_N"].SetDefault(0).Get<int>();
  std::vector<std::string> dipoles = s["DIPOLE_FLAV"].GetVector<std::string>();
  for (auto d : dipoles) {
    std::vector<int> tmp;
    size_t pos(d.find(" "));
    tmp.push_back(std::stoi(d.substr(0, pos)));
    tmp.push_back(std::stoi(d.substr(pos + 1)));
    m_dip.push_back(tmp);
  }
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  #ifdef USING__LOOPTOOLS
    FORTRAN(ltini)();
    Setlambda(0);
  #endif
}

Define_Dipoles::~Define_Dipoles() {
  if(p_yfsFormFact) delete p_yfsFormFact;
}


void Define_Dipoles::MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  if(!HasISR()) return;
  if ((mom.size() < 2 || fl.size() < 2) ) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::initial<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for dipoletype");
  }
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_test_dip.clear();
  m_flav_label.clear();
  m_softphotonSum *= 0;
  m_out = fl.size() - m_in;
  m_olddipoles.clear();
  m_dipolesII.clear();
  m_bornmomenta = born;
  Dipole_II(fl, mom);
}


void Define_Dipoles::MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born) {
  if(m_mode==yfsmode::fsr) return;
  // if(m_ifisub==0) return;
  if ((mom.size() != fl.size())) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::ifi<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for dipoletype");
  }
  if (!HasFSR() ) return;
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  // m_test_dip.clear();
  // m_flav_label.clear();
  // m_softphotonSum *= 0;
  m_out = fl.size() - m_in;
  // m_olddipoles.clear();
  m_dipolesIF.clear();
  Dipole_IF(fl, mom, born);
}

void Define_Dipoles::MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  if ((mom.size() != fl.size())) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::ifi<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for dipoletype");
  }
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_test_dip.clear();
  m_flav_label.clear();
  m_softphotonSum *= 0;
  m_out = fl.size() - m_in;
  m_olddipoles.clear();
  m_dipolesFF.clear();
  m_bornmomenta = born;
  Dipole_FF(fl, mom);
}

void Define_Dipoles::MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born ) {
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_test_dip.clear();
  m_flav_label.clear();
  m_softphotonSum *= 0;
  m_bornmomenta = born;
  m_out = fl.size() - m_in;
  m_olddipoles.clear();
  m_dipolesFF.clear();
  int fsrc=0;
  for (int i = 2; i < fl.size(); ++i) if(fl[i].Charge()!=0) fsrc++;
  // m_dipolesIF.clear();
  for(size_t i = 0; i < fl.size(); ++i)
  {
    m_flav_label[fl[i]] = i;
  }
  if (!HasFSR() ) return;
  if (fl.size() == 4) {
    Flavour_Vector ff;
    Vec4D_Vector mm, bm;
    m_flav_label[fl[2]] = 2;
    m_flav_label[fl[3]] = 3;
    for(size_t i = 2; i < fl.size(); i++) {
      // if(fl[i].IntCharge()!=0 && !fl[i].IsQCD()){
      if(fl[i].IntCharge()!=0){
        ff.push_back(fl[i]);
        mm.push_back(mom[i]);
        bm.push_back(m_bornmomenta[i]);
      }
    }
    if(ff.size()==0) return;
    Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
    D.SetResonance(true);
    D.SetFlavLab(2,3);
    // IsResonant(D);
    Dipole_FF(ff, mm);
    m_dipolesFF.push_back(D);
    return;
  }
  map<ATOOLS::Flavour, ATOOLS::Vec4D>::iterator itr;
  if (m_dip.size() != 0) {
    for (auto a : m_dip) {
      Get4Mom(fl, mom); // makes map for flavour momentum
      // else Get4Mom(fl,mom);
      Flavour_Vector ff;
      Flavour f;
      Vec4D_Vector mm, bm;
      // for(auto f: fl){
      for(size_t i = 2; i < fl.size(); ++i)
      {
        f = fl[i];
        if (f.Charge()!=0) {
          // msg_Error()<<"YFS FSR only defined for Final state leptons and not for "<<f<<std::endl;
          // }
          ff.push_back(f);
          mm.push_back(m_test_dip[f]);
          bm.push_back(m_born_dip[f]);
          m_flav_label[f] = i;
          if (!IsEqual(f.Mass(), m_test_dip[f].Mass(), 1e-5)) {
            msg_Error() << "Incorrect mass mapping in dipole" << std::endl
                        << "Flavour mass is " << f.Mass() << std::endl
                        << "Four-Momentum mass is " << m_test_dip[f].Mass() << std::endl;
          }
          // m_mom_label[m_test_dip[f]] = i;
          if (ff.size() == 2) break;
        }
        // i+=1;
      }
      Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
      Dipole_FF(ff, mm);
      if(fsrc==2)  {
        D.SetResonance(true);
        D.SetFlavLab(2,3);
      }
      else IsResonant(D);
      m_dipolesFF.push_back(D);
      msg_Debugging() << "Added " << ff << " to dipole " << a << std::endl;
    }
  }
  else {
    Get4Mom(fl, mom);
    Flavour_Vector ff;
    Vec4D_Vector mm, bm;
    int N = 0; // number of leptons minus the inital state
    // for (auto f : fl) {
    for(int i=2; i < fl.size(); i++){
      if (fl[i].Charge()!=0) N += 1;
    }
    if (N == 2) {
      //only two leptons in final state
      // one dipole
      Flavour_Vector ff;
      Vec4D_Vector mm, bm;
      std::vector<int> id;
      // m_flav_label[fl[2]] = 2;
      // m_flav_label[fl[3]] = 3;
      for(size_t i = 2; i < fl.size(); i++) {
        if (fl[i].Charge()!=0) {
          ff.push_back(fl[i]);
          mm.push_back(mom[i]);
          bm.push_back(m_bornmomenta[i]);
          id.push_back(i);
        }
      }
      Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
      if(fsrc==2)  D.SetResonance(true);
      else IsResonant(D);
      Dipole_FF(ff, mm);
      D.SetFlavLab(id[0],id[1]);
      m_dipolesFF.push_back(D);
      return;
    }
    vector<vector<int>> pairings;
    vector<int> curr_pairing, available_nums;
    for(int i = 1; i <= N; i++) {
      available_nums.push_back(i);
    }
    generate_pairings(pairings, curr_pairing, available_nums);
    int k = 0;
    int d1, d2;
    for(size_t i = 0; i < pairings.size(); i++) {
      for(size_t j = 0; j < pairings[i].size(); j++) {
        if (k == 0) d1 = pairings[i][j] + 1;
        else if (k == 1) d2 = pairings[i][j] + 1;
        // PRINT_VAR(pairings[i][j] + 1);
        k += 1;
        if (k == 2) {
          Flavour f1 = fl[d1];
          Flavour f2 = fl[d2];
          if(f1.Charge()!=0 && f2.Charge()!=0){
            ff.push_back(f1);
            ff.push_back(f2);
            mm.push_back(mom[d1]);
            mm.push_back(mom[d2]);
            bm.push_back(m_bornmomenta[d1]);
            bm.push_back(m_bornmomenta[d2]);
            Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
            Dipole_FF(ff, mm);
            IsResonant(D);
            // if(D.IsResonance()) PRINT_VAR(D);
            D.SetFlavLab(d1,d2);
            m_dipolesFF.push_back(D);
            ff.clear();
            mm.clear();
            bm.clear();
            k = 0;
          }
        }
      }
    }
  }
}

void Define_Dipoles::Get4Mom(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector mom) {
  Vec4D_Vector P;
  for(size_t i = 2; i < fl.size(); ++i)
  {
    if (fl[i].IntCharge()!=0) {
      m_test_dip[fl[i]] = mom[i];
      P.push_back(mom[i]);
      if (P.size() == 2) break;
    }
  }
  if (P.size() != 2) {
    PRINT_VAR(P.size());
    THROW(fatal_error, "Wrong size dipole");
  }
}


void Define_Dipoles::Dipole_II(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom) {
  CleanInParticles();
  Flavour_Vector ff;
  Vec4D_Vector mm, bm;
  for(size_t i = 0; i < 2; i++) {
    ff.push_back(fl[i]);
    mm.push_back(mom[i]);
    bm.push_back(m_bornmomenta[i]);
  }
  Dipole D(ff, mm, bm, dipoletype::initial,m_alpha);
  m_g=D.m_gamma;
  m_gp=D.m_gammap;
  D.SetFlavLab(0, 1);
  m_olddipoles.push_back(D);
  m_dipolesII.push_back(D);
}


void Define_Dipoles::Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom) {
  CleanOutParticles();
  if (fl.size() != mom.size()) {
    THROW(fatal_error, "Inconsistent flavour vector for Dipole_FF Momenta");
  }
  for(size_t i = 0; i < fl.size(); ++i)
  {
    if (fl[i].IsQED()) {
      m_chargedoutparticles.push_back(mom[i]);
      m_massOutC.push_back(mom[i].Mass());

    }
    else {
      m_neutraloutparticles.push_back(mom[i]);
      m_massOutN.push_back(mom[i].Mass());
    }
  }
}



void Define_Dipoles::Dipole_IF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  CleanInParticles();
  if (fl.size() != mom.size()) {
    THROW(fatal_error, "Inconsistent flavour vector for Dipole_IF Momenta");
  }
  Flavour_Vector ff;
  Vec4D_Vector mm, bm;
  //create IF dipoles
    for(size_t i = 0; i < 2; ++i)
    {
      for(size_t j = 2; j < fl.size(); ++j)
      {
        if(fl[i].IntCharge()==0) continue;
        if(fl[j].IntCharge()==0) continue;
        ff.clear();
        mm.clear();
        bm.clear();
        ff.push_back(fl[i]);
        ff.push_back(fl[j]);

        mm.push_back(mom[i]);
        mm.push_back(mom[j]);


        bm.push_back(born[i]);
        bm.push_back(born[j]);
        Dipole D(ff, mm, bm, dipoletype::ifi,m_alpha);
        D.SetResonance(false);
        m_dipolesIF.push_back(D);
      }
    }
}



double Define_Dipoles::CalculateRealSub(const Vec4D &k) {
  double sub(0);
  Vec4D eik{0.,0.,0.,0.};
  for (auto &D : m_dipolesII) {
    for(size_t i = 0; i < D.GetBornMomenta().size(); ++i)
    {
      Vec4D p = D.GetBornMomenta(i);
      eik += D.m_Q[i]*p/(p*k);
    }
  }
  for (auto &D : m_dipolesFF) {
    for(size_t i = 0; i < D.GetBornMomenta().size(); ++i)
    {
      Vec4D p = D.GetBornMomenta(i);
      eik += -D.m_Q[i]*p/(p*k);
    }
  }
  sub = -m_alpha / (4 * M_PI * M_PI)*eik*eik;
  return sub;
}

double Define_Dipoles::CalculateRealSubIF(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_dipolesIF){
    // if(k.E() < sqrt(m_s)/100) sub -= D.EikonalMassless(k, D.GetMomenta(0), D.GetMomenta(1));
    if(m_massless_sub) sub += D.EikonalMassless(k, D.GetMomenta(0), D.GetMomenta(1));
    else sub +=  D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
    // sub += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  }
  return sub;
}


double Define_Dipoles::CalculateVirtualSub() {
  double sub(0);
  if(m_tchannel>=2) return CalculateVirtualSubTchannel();
  if(m_dim_reg==1) return CalculateVirtualSubEps();
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_dipolesFF) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
  }
  return sub;
}

double Define_Dipoles::CalculateVirtualSubEps() {
  DivArrD sub(0);
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }
  for (auto &D : m_dipolesFF) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }

  for (auto &D : m_dipolesIF){
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) {
      msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
      // THROW(fatal_error, "YFS Subtraction fails");
    }
  }
  m_virtSub=sub;
  return sub.Finite();
}


double Define_Dipoles::CalculateRealVirtualSubEps(const Vec4D &k) {
  DivArrD sub(0);
  if(m_tchannel>=2) return CalculateVirtualSubTchannel();
  // if(m_tchannel==2) return CalculateRealSub(k);
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVR_full_eps(D, sqrt(m_s) / 2., 0);
  }
  for (auto &D : m_dipolesFF) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVR_full_eps(D, sqrt(m_s) / 2., 0);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVR_full_eps(D, sqrt(m_s) / 2., 0);
  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVR_full_eps(D, sqrt(m_s) / 2., 0);
  }
  m_virtSub=sub;
  return sub.Finite();
}



double Define_Dipoles::FormFactor(){
  double form = 0;

  for(auto &D: m_dipolesII){
    form+= D.ChargeNorm()*p_yfsFormFact->BVR_full(D.GetBornMomenta(0), D.GetBornMomenta(1), sqrt(m_s) / 2.);
  }
  if(!m_hidephotons){
      for(auto &D: m_dipolesFF){
        form += D.ChargeNorm()*p_yfsFormFact->BVR_full(D.GetBornMomenta(0), D.GetBornMomenta(1), sqrt(m_s) / 2.);
      }
    if(m_ifisub==1){
      for(auto &D: m_dipolesIF){
        form += D.ChargeNorm()*p_yfsFormFact->R1(D);
      }
    }
  }
  if(FixedOrder()==fixed_order::nlo){
    return 1.+form;
  }
  return exp(form); 
}


double Define_Dipoles::TFormFactor(){
  double form = 0;
  for(auto &D: m_dipolesII){
    form+= D.ChargeNorm()*p_yfsFormFact->R1(D);
  }
  // if(!m_hidephotons){
    for(auto &D: m_dipolesFF){
      form += D.ChargeNorm()*p_yfsFormFact->R1(D);
    // }
  }
  if(m_ifisub==1){
    for(auto &D: m_dipolesIF){
      form += D.ChargeNorm()*p_yfsFormFact->R1(D);
    }
  }
  return exp(form); 
}

double Define_Dipoles::CalculateVirtualSubTchannel(){
   // YFSij = 2.d0*B0ij - B0ii - B0jj
   //   .         + 4.d0 * mi2 * C0singular(mi2,phmass)
   //   .         + 4.d0 * mj2 * C0singular(mj2,phmass)
   //   .         + 8.d0*pi_pj * C0ij
  if(m_dim_reg) return CalculateVirtualSubTchannelEps();
  double sub(0);
  // Vec4D_Vector pvirt;
  // std::vector<double> z,th;
  // pvirt.push_back(m_dipolesII[0].GetNewMomenta(0));
  // pvirt.push_back(m_dipolesII[0].GetNewMomenta(1));
  // pvirt.push_back(m_dipolesFF[0].GetBornMomenta(0));
  // pvirt.push_back(m_dipolesFF[0].GetBornMomenta(1));
  // th.push_back(1);
  // th.push_back(1);
  // th.push_back(-1);
  // th.push_back(-1);
  // z.push_back(m_dipolesII[0].m_Qi);
  // z.push_back(m_dipolesII[0].m_Qj);
  // z.push_back(m_dipolesFF[0].m_Qi);
  // z.push_back(m_dipolesFF[0].m_Qj);
  // // double m1 = pvirt[0].Mass();
  // // double m2 = pvirt[1].Mass();
  // // double m3 = pvirt[2].Mass();
  // // double m4 = pvirt[3].Mass();
  // for(size_t i = 0; i < pvirt.size(); ++i)
  // {
  //   for(size_t j=i; j<pvirt.size(); ++j ){
  //     double etaij = z[i]*z[j]*th[i]*th[j];
  //     Complex YFSij = 0.;
  //     double mi = pvirt[i].Mass();
  //     double mj = pvirt[j].Mass();
  //     double s = (pvirt[i]-pvirt[j]).Abs2();

  //     Complex bii = B0(0,mi*mi,mi*mi);
  //     Complex bjj = B0(0,mj*mj,mj*mj);
  //     Complex bij = B0(0,mi*mi,mj*mj);
  //     Complex cij = C0(mi*mi,(th[i]*pvirt[i]+th[j]*pvirt[j]).Abs2(),mj*mj,
  //                      0.,mi*mi,mj*mj);
  //     Complex cii = C0(mi*mi, 0., mi*mi, 0.0, mi*mi, mi*mi);
  //     // YFSij = 8*pvirt[i]*pvirt[j]*cij;
  //     if(i==j){
  //       YFSij = bii-4.*mi*mi*0.5*log(m_photonMass*m_photonMass/mi/mi);
  //     }
  //     else{
  //       YFSij = 2.*pvirt[i]*pvirt[j]*cij+0.5*bij;
  //       }
  //     if(IsBad(YFSij)){
  //       msg_Error()<<"YFS Virtual Sub is NaN"<<endl
  //                  <<"bii = "<<bii<<endl
  //                  <<"bij = "<<bij<<endl
  //                  <<"cii = "<<cii<<endl
  //                  <<"cij = "<<cij<<endl;
  //     }
  //     sub+=m_alpi*etaij*YFSij.real();
  // //     // PRINT_VAR(etaij*YFSij);
  //   }
  // }
  // PRINT_VAR(count);
  for (auto &D : m_dipolesII){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  for (auto &D : m_dipolesFF){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  for (auto &D : m_dipolesIF){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  // clearcache();
  return sub;
}

double Define_Dipoles::CalculateVirtualSubTchannelEps() {
  DivArrD sub(0);
  DivArrD massph(0,-1.,0,0,0,0);
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    #ifdef USING__LOOPTOOLS
      //ii term
      Vec4D p1 = D.GetBornMomenta(0);
      Vec4D p2 = D.GetBornMomenta(1);
      double m1 = D.GetMass(0);
      double m2 = D.GetMass(1);
      double irloop = p_yfsFormFact->p_virt->IRscale();
      double epsloop = p_yfsFormFact->p_virt->Eps_Scheme_Factor({p1,p1});
      DivArrD c0 = (-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop));
      // sub += m_alpi/8.*B0(0, m1*m1, m1*m1).real()-8.*m_alpi*(m1*m1)*(-massph-log(4.*M_PI*sqr(irloop)/m1/m1/epsloop));
      // sub += m_alpi/8.*B0(0, m2*m2, m2*m2).real()-8.*m_alpi*(m2*m2)*(-massph-log(4.*M_PI*sqr(irloop)/m2/m2/epsloop));
    #endif
  }
  for (auto &D : m_dipolesFF) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    #ifdef USING__LOOPTOOLS
      //ii term
      Vec4D p1 = D.GetBornMomenta(0);
      Vec4D p2 = D.GetBornMomenta(1);
      double m1 = D.GetMass(0);
      double m2 = D.GetMass(1);
      double irloop = p_yfsFormFact->p_virt->IRscale();
      double epsloop = p_yfsFormFact->p_virt->Eps_Scheme_Factor({p1,p1});
      DivArrD c0 = (-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop));
      // sub += m_alpi/8.*B0(0, m1*m1, m1*m1).real()-8.*m_alpi*(m1*m1)*(-massph-log(4.*M_PI*sqr(irloop)/m1/m1/epsloop));
      // sub += m_alpi/8.*B0(0, m2*m2, m2*m2).real()-8.*m_alpi*(m2*m2)*(-massph-log(4.*M_PI*sqr(irloop)/m2/m2/epsloop));
    #endif
  }
  m_virtSub=sub;
  return sub.Finite();
}

double Define_Dipoles::CalculateRealVirtualSub(const Vec4D & k) {
  double sub(0);
  for (auto &D : m_dipolesII) {
    // sub += -D.Eikonal(k);
    sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D.GetNewMomenta(0), D.GetNewMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_dipolesFF) {
    sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D.GetOldMomenta(0), D.GetOldMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);

  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.m_QiQj*p_yfsFormFact->BVV_full(D.GetBornMomenta(0), D.GetBornMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);
  }
  return sub;
}


double Define_Dipoles::CalculateEEX(){
  double eex=0;
  for (auto &D: m_dipolesII){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  for (auto &D: m_dipolesFF){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  for (auto &D: m_dipolesIF){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  return eex;
}

double Define_Dipoles::CalculateEEXVirtual(){
  double vint{1.}, vfin{1};
  for (auto &D: m_dipolesII){
    vint*=1+D.VirtualEEX(m_betaorder);
  }
  for (auto &D: m_dipolesFF){
    vfin*=1+D.VirtualEEX(m_betaorder);
  }
  return vint*vfin;
}

double Define_Dipoles::EEXRealVirtual(const Vec4D &k){
  double eex = 0;
  for(auto &D: m_dipolesII){
    D.m_betaorder = 2;
    eex += D.Beta1(k)/D.Eikonal(k);
  }
  for(auto &D: m_dipolesFF){
    D.m_betaorder = 2;
    eex += D.Beta1(k)/D.Eikonal(k);
  }
  // for(auto &D: m_dipolesIF){
  //   D.m_betaorder = 2;
  //   eex += D.Beta1(k)/D.Eikonal(k);
  // }
  if(IsNan(eex)) return 0;
  return eex;
}

double Define_Dipoles::CalculateRealSubEEX(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_dipolesII) {
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  for (auto &D : m_dipolesFF) {
    // double norm=1;
    // if(D.m_Qi == D.m_Qj) norm=-1;
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  // for (auto &D : m_dipolesIF) {
  //   sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  // }

  return sub;
}


void Define_Dipoles::CleanInParticles() {
  m_chargedinparticles.clear();
  m_neutralinparticles.clear();
  m_massInC.clear();
  m_massInN.clear();
}

void Define_Dipoles::CleanOutParticles() {
  m_chargedoutparticles.clear();
  m_neutraloutparticles.clear();
  m_massOutC.clear();
  m_massOutN.clear();
}

void Define_Dipoles::CleanUp() {
  m_dipoles.clear();
}

double Define_Dipoles::CalculateFlux(const Vec4D &k){
  if(!HasISR()) return 1;
  double sq, sx;
  double flux = 1;
  dipoletype::code fluxtype;
  Vec4D Q,QX;
  if(m_noflux==1) return 1;
  if(HasISR()&&HasFSR()){
    // fluxtype = WhichResonant(k);
    // PRINT_VAR(fluxtype);
    fluxtype = dipoletype::final;
  }
  else if(!HasFSR()){
    fluxtype = dipoletype::initial;
  }
  else if(!HasISR()){
    fluxtype = dipoletype::final;
  }
  else{
    msg_Error()<<"Unknown dipole type in "<<METHOD<<std::endl;
  }
  if(fluxtype==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2(); 
      sx = (Q-k).Abs2();
      flux = (sx/sq);
      return flux;
    }

  }
  if(fluxtype==dipoletype::final){
    flux=0;
    for (auto &D : m_dipolesFF) {
      Q  = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux += (sq/sx);
      // flux = Propagator(sx)/Propagator(sq);
    }
    return flux/m_dipolesFF.size();
  }
  return flux;
}

double Define_Dipoles::CalculateFlux(const Vec4D &k, dipoletype::code &fluxtype){
  double sq, sx;
  double flux = 1;
  Vec4D Q,QX;
  if(m_noflux==1) return 1;
  // if(HasISR()&&HasFSR()){
  //   fluxtype = WhichResonant(k);
  //   // PRINT_VAR(fluxtype);
  //   // fluxtype = dipoletype::final;
  // }
  // else if(!HasFSR()){
  //   fluxtype = dipoletype::initial;
  // }
  // else if(!HasISR()){
  //   fluxtype = dipoletype::final;
  // }
  // else{
  //   msg_Error()<<"Unknown dipole type in "<<METHOD<<std::endl;
  // }
  if(fluxtype==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2(); 
      sx = (Q-k).Abs2();
      flux = (sx/sq);
      return flux;
    }

  }
  if(fluxtype==dipoletype::final){
    flux=0;
    for (auto &D : m_dipolesFF) {
      Q  = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux += (sq/sx);
      // flux = Propagator(sx)/Propagator(sq);
    }
    return flux/m_dipolesFF.size();
  }
  return flux;
}

double Define_Dipoles::CalculateFlux(const Vec4D &k, const Vec4D &kk){
  double sq, sx;
  double flux = 1;
  Vec4D Q,QX;
  dipoletype::code fluxtype1, fluxtype2;
  if(m_noflux==1) return 1;
  fluxtype1 = WhichResonant(k);
  fluxtype2 = WhichResonant(kk);
  // if(HasISR()&&HasFSR()){
  //   // fluxtype = WhichResonant(k);
  //   // fluxtype = dipoletype::final;
  // }
  // else if(HasISR()){
  //   fluxtype = dipoletype::initial;
  // }
  // else if(HasFSR()){
  //   fluxtype = dipoletype::final;
  // }
  // else{
  //   msg_Error()<<"Unknown dipole type in "<<METHOD<<std::endl;
  // }
  if(fluxtype1==dipoletype::initial && fluxtype2==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-k-kk).Abs2();
      flux = sx/sq;
      return flux;
    }
  }
  else if(fluxtype1==dipoletype::final && fluxtype2==dipoletype::final){
    for (auto &D : m_dipolesFF) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k+kk).Abs2();
      flux = sq/sx;
    }
    return flux;
  }
  else if(fluxtype1==dipoletype::initial && fluxtype2==dipoletype::final){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-k).Abs2();
      flux = sx/sq;
    }
    for (auto &D : m_dipolesFF) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+kk).Abs2();
      flux*= sq/sx;
   }
  }
  else if(fluxtype1==dipoletype::final && fluxtype2==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-kk).Abs2();
      flux = sx/sq;
    }
    for (auto &D : m_dipolesFF) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux*= sq/sx;
    }
  }
  else{
    msg_Error()<<"Unknown flux type in "<<METHOD<<" with \nfluxtype1 = "<<fluxtype1<<"\n and fluxtype2 = "<< fluxtype2 <<std::endl;
  }
  return flux;
}


double Define_Dipoles::Propagator(const double &s, int width){
  Flavour fl;
  Complex Prop = Complex(0.,0.);///Complex(s,0.0);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      fl = v->in[0];
      if(IsZero(fl.Mass())) continue;
      Prop += Complex(1.,0.)/Complex(s-sqr(fl.Mass()),fl.Width()*fl.Mass());
    }
  }
  return ((Prop*conj(Prop)).real());
}

void Define_Dipoles::IsResonant(YFS::Dipole &D) {
double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist;
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(D.m_QiQj==1 || !D.IsDecayAllowed()){
        D.SetResonance(false);
        continue;
        }   
      mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
      if(mdist<m_resonace_max) {
        D.SetResonance(true);
        return;
      }
      else D.SetResonance(false);
    }
    D.SetResonance(false);
  }
}

bool Define_Dipoles::CheckResonant(YFS::Dipole &D) {
  double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist;
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
      if(mdist<5) {
        return true;
      }
    }
  }
  return false;
}

bool Define_Dipoles::IsResonant() {
  for(auto &D: m_dipolesFF){
    double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
    double mdist;
    for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
      for (auto *v : it->second) {
        mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
        if(mdist<5) {
          return true;
        }
      }
    }
  }
  return false;
}

bool Define_Dipoles::CheckResonant(){
  bool isres = false;
  for(auto &D: m_dipolesII){
    if(CheckResonant(D)) isres=true;
  }
  for(auto &D: m_dipolesFF){
    if(CheckResonant(D)) isres=true;
  }
  return isres;
}

double Define_Dipoles::ResonantDist(YFS::Dipole &D, const Vec4D &k){
  double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)-k).Mass();
  double mass_i = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist(100000000);
  double mcheck(100000000);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width()) ) continue;
      mcheck = abs(mass_d- mass_i);
      if(mcheck < mdist) {
        mdist = mcheck;
        D.SetResonaceFlavour(v->in[0]);
      }
    }
  }
  return mdist;
}

double Define_Dipoles::ResonantDist(YFS::Dipole &D){
  double mass_i = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist(100000000);
  double mcheck(100000000);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width()) ) continue;
      mcheck = abs(mass_i - v->in[0].Mass());
      if(mcheck < mdist) mdist = mcheck;
    }
  }
  return mcheck;
}

dipoletype::code Define_Dipoles::WhichResonant(const Vec4D &k){
  if(!HasFSR()) return dipoletype::initial;
  if(!HasISR()) return dipoletype::final;
  double mdistisr(10000), mdistfsr(100000),mdistifi(100000);
  double mindis(10000);
  dipoletype::code min(dipoletype::initial);
  for(auto &D: m_dipolesII){
    mdistisr = ResonantDist(D,k);
    mindis = mdistisr;
    min = dipoletype::initial;
  }
  for(auto &D: m_dipolesFF){
    mdistfsr = ResonantDist(D,k);  
    if(mdistfsr < mdistisr){
      min = dipoletype::final;
    }
  }
  // for(auto &D: m_dipolesIF){
  //   mdistifi = ResonantDist(D,k);  
  //   if(mdistifi < mdistfsr) {
  //     min = dipoletype::ifi;
  //     // return dipoletype::final; 
  //   }
  // }
  // PRINT_VAR(min);
  return min;
}

void Define_Dipoles::generate_pairings(std::vector<std::vector<int>>& pairings, std::vector<int>& curr_pairing, std::vector<int>& available_nums) {
  if (available_nums.empty()) {
    pairings.push_back(curr_pairing);
    return;
  }
  int curr_num = available_nums[0];
  available_nums.erase(available_nums.begin());
  for(size_t i = 0; i < available_nums.size(); i++) {
    int next_num = available_nums[i];
    available_nums.erase(available_nums.begin() + i);
    curr_pairing.push_back(curr_num);
    curr_pairing.push_back(next_num);
    generate_pairings(pairings, curr_pairing, available_nums);
    curr_pairing.pop_back();
    curr_pairing.pop_back();
    available_nums.insert(available_nums.begin() + i, next_num);
  }
  available_nums.insert(available_nums.begin(), curr_num);
}

std::ostream& Define_Dipoles::operator<<(std::ostream &out) {
  out << "N_in = " << m_in << "\n m_out = " << m_out <<
      "Number of Charged incoming particles = " << m_chargedinparticles.size() << std::endl <<
      "Number of Charged outgoing particles = " << m_chargedoutparticles.size() << std::endl <<
      "Number of Neutral incoming particles = " << m_neutralinparticles.size() << std::endl <<
      "Number of Neutral outgoing particles = " << m_neutraloutparticles.size() << std::endl;
  return out;
}
