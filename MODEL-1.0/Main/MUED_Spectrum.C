#include "MUED_Spectrum.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

MUED_Spectrum::MUED_Spectrum(Data_Reader * _dataread,Model_Base * _model) :
  Spectrum_Generator_Base(_dataread,_model)  
{
  msg_Out()<<"In "<<METHOD<<" for "<<p_model->Name()<<endl;
  m_generations = (*p_model->GetScalarNumbers())["N_Generations"];
  m_invR        = (*p_model->GetScalarConstants())["1/Radius"];
  m_Lambda      = (*p_model->GetScalarConstants())["Lambda"];
  m_Lambda2     = sqr(m_Lambda);

  m_sin2thetaW  = (*p_model->GetScalarConstants())["sin2_thetaW"];
  m_cos2thetaW  = (*p_model->GetScalarConstants())["cos2_thetaW"];
  m_v           = (*p_model->GetScalarConstants())["vev"];

  m_Riemann3    = 1.20206;
}

MUED_Spectrum::~MUED_Spectrum() 
{ 
  if (!m_decays.empty()) 
    m_decays.erase(m_decays.begin(),m_decays.end());
}

void MUED_Spectrum::Run(string output) {
  msg_Out()<<"In "<<METHOD<<endl;
  Neutral_KK_Bosons();
  Neutral_KK_Scalars();
  LR_KK_Fermions();
}

void MUED_Spectrum::FillMasses() {}

void MUED_Spectrum::Neutral_KK_Bosons() {
  msg_Out()<<"In "<<METHOD<<endl;
  double delta_m2_gn(0.);
  double diag_m2_Bn(0.),delta_m2_Bn(0.),diag_m2_Wn(0.),delta_m2_Wn(0.),off_diag2(0.);
  double a_mn2(0.),as_mn2(0.),mn2(0.);
  double mAn(0.), mZn(0.), mGn(0.);

  Matrix<2> ZNeutral2, Eigenvectors;
  double Eigenvalues[2];
  int    kfno(5000000);

  cout<<METHOD<<" : "<<m_generations<<endl;
  for (int i=1;i<m_generations+1;i++) {
    mn2         = sqr(i*m_invR);
    a_mn2       = (*p_model->GetScalarFunction("alpha_QED"))(mn2);
    as_mn2      = (*p_model->GetScalarFunction("alpha_S"))(mn2);
    delta_m2_gn = 
      as_mn2/(4.*M_PI) * 
      (23./2.*mn2*log(m_Lambda2/mn2) -
       3./2.*m_Riemann3*sqr(m_invR/M_PI));

    diag_m2_Bn   = 
      M_PI*a_mn2/m_cos2thetaW*sqr(m_v);
    delta_m2_Bn = 
      a_mn2/(4.*M_PI*m_cos2thetaW)*
      (-1./6.*mn2*log(m_Lambda2/mn2) -
       39./2.*m_Riemann3*sqr(m_invR/M_PI));
    diag_m2_Wn   = 
      M_PI*a_mn2/m_sin2thetaW*sqr(m_v);
    delta_m2_Wn =
      a_mn2/(4.*M_PI*m_sin2thetaW)*
      (15./2.*mn2*log(m_Lambda2/mn2) -
       5./2.*m_Riemann3*sqr(m_invR/M_PI));
    off_diag2 =
      M_PI*a_mn2/sqrt(m_sin2thetaW*m_cos2thetaW)*sqr(m_v);

    ZNeutral2[0][0] = mn2+diag_m2_Bn+delta_m2_Bn;
    ZNeutral2[0][1] = off_diag2;
    ZNeutral2[1][0] = off_diag2;
    ZNeutral2[1][1] = mn2+diag_m2_Wn+delta_m2_Wn;

    ZNeutral2.DiagonalizeSort(Eigenvalues,Eigenvectors);

    mGn = sqrt(mn2+delta_m2_gn);
    mAn = sqrt(Eigenvalues[0]);
    mZn = sqrt(Eigenvalues[1]);

    cout<<METHOD<<" : "<<as_mn2<<endl;
    cout<<"("<<ZNeutral2[0][0]<<", "<<ZNeutral2[0][1]<<", "<<ZNeutral2[1][0]<<", "<<ZNeutral2[1][1]
	<<") --> G_n = "<<mGn<<", A_n = "<<mAn<<", Z_n = "<<mZn<<endl;


    cout<<" Mix : ("<<Eigenvectors[0][0]<<"  "<<Eigenvectors[0][1]<<endl
	<<"        "<<Eigenvectors[1][0]<<"  "<<Eigenvectors[1][1]<<")."<<endl;

    Flavour Gn((kf_code)(int(kfno+i*100000+21))); Gn.SetMass(mGn);
    Flavour An((kf_code)(int(kfno+i*100000+22))); An.SetMass(mAn);
    Flavour Zn((kf_code)(int(kfno+i*100000+24))); Zn.SetMass(mZn);

    if (i==1) {
      std::cout<<"Check "<<Flavour(kf_KK_gluon_1).Mass()
	       <<"/"<<Flavour(kf_KK_B1_1).Mass()
	       <<"/"<<Flavour(kf_KK_W3_1).Mass()<<endl;
    }
    else if (i==2) {
      std::cout<<"Check "<<Flavour(kf_KK_gluon_2).Mass()
	       <<"/"<<Flavour(kf_KK_B1_2).Mass()
	       <<"/"<<Flavour(kf_KK_W3_2).Mass()<<endl;
    }
  }
}

void MUED_Spectrum::Charged_KK_Bosons() {
}

void MUED_Spectrum::Neutral_KK_Scalars() {
  msg_Out()<<"In "<<METHOD<<endl;
  double delta_m2_Hn(0.);
  double a_mn2(0.),mn(0.),mn2(0.);
  double mh2(sqr(p_model->ScalarConstant("MH")));
  double mz2(sqr(p_model->ScalarConstant("MZ")));
  double mw2(sqr(p_model->ScalarConstant("MW")));
  double lambdaH(mh2/sqr(m_v));
  Flavour flav;

  cout<<METHOD<<" : "<<m_generations<<endl;
  for (int i=1;i<m_generations+1;i++) {
    mn          = i*m_invR;
    mn2         = mn*mn;
    a_mn2       = (*p_model->GetScalarFunction("alpha_QED"))(mn2);
    delta_m2_Hn = 
      p_model->ScalarConstant("M2bar_H") +
      mn2/(4.*M_PI)*log(m_Lambda2/mn2)*
      (3./2.*a_mn2/m_sin2thetaW + 3./4.*a_mn2/m_cos2thetaW - lambdaH/(4.*M_PI));
    flav = Flavour((kf_code)(5000000+i*100000+35));
    flav.SetMass(sqrt(mn2+delta_m2_Hn+mh2));
    flav = Flavour((kf_code)(5000000+i*100000+36));
    flav.SetMass(sqrt(mn2+delta_m2_Hn+mz2));
    flav = Flavour((kf_code)(5000000+i*100000+37));
    flav.SetMass(sqrt(mn2+delta_m2_Hn+mw2));

    if (i==1) {
      std::cout<<"Check "<<Flavour(kf_KK_H_1).Mass()
	       <<"/"<<Flavour(kf_KK_A_1).Mass()
	       <<"/"<<Flavour(kf_KK_Hplus_1).Mass()<<endl;
    }
    else if (i==2) {
      std::cout<<"Check "<<Flavour(kf_KK_H_2).Mass()
	       <<"/"<<Flavour(kf_KK_A_2).Mass()
	       <<"/"<<Flavour(kf_KK_Hplus_2).Mass()<<endl;
    }
  }
}

void MUED_Spectrum::Charged_KK_Scalars() {
}

void MUED_Spectrum::LR_KK_Fermions() {
  msg_Out()<<"In "<<METHOD<<endl;
  double delta_m_Qn(0.),delta_m_un(0.),delta_m_dn(0.),offdiag_un(0.),offdiag_dn(0.);
  double delta_m_Ln(0.),delta_m_en(0.),offdiag_en(0.);
  double a_mn2(0.),as_mn2(0.),mn(0.),mn2(0.);

  Matrix<2> ZULR, ZDLR, ZLLR, ZULR2, ZDLR2, ZLLR2;
  Matrix<2> EigenvectorsU, EigenvectorsD, EigenvectorsL;
  double EigenvaluesU[2], EigenvaluesD[2], EigenvaluesL[2];
  int    kfnoL(5000000),kfnoR(6000000);

  Flavour flav;

  cout<<METHOD<<" : "<<m_generations<<endl;
  for (int i=1;i<m_generations+1;i++) {
    mn          = i*m_invR;
    mn2         = mn*mn;
    a_mn2       = (*p_model->GetScalarFunction("alpha_QED"))(mn2);
    as_mn2      = (*p_model->GetScalarFunction("alpha_S"))(mn2);
    for (int j=1;j<4;j++) {
      if (j==3) {
	offdiag_un = p_model->ScalarConstant("Yukawa_t");
	offdiag_dn = p_model->ScalarConstant("Yukawa_b");
	offdiag_en = p_model->ScalarConstant("Yukawa_tau");
      }
      else {
	offdiag_un = 0.;
	offdiag_dn = 0.;
	offdiag_en = 0.;
      }
      delta_m_Qn = 
	mn/(4.*M_PI)*log(m_Lambda2/mn2) *
	(3.*as_mn2 + 27./16.*a_mn2/m_sin2thetaW + 1./16.*a_mn2/m_cos2thetaW);
      delta_m_un = 
	mn/(4.*M_PI)*log(m_Lambda2/mn2) *
	(3.*as_mn2 + a_mn2/m_cos2thetaW);
      delta_m_dn = 
	mn/(4.*M_PI)*log(m_Lambda2/mn2) *
	(3.*as_mn2 + 1./4.*a_mn2/m_cos2thetaW);
      delta_m_Ln = 
	mn/(4.*M_PI)*log(m_Lambda2/mn2) *
	(27./16.*a_mn2/m_sin2thetaW + 9./16.*a_mn2/m_cos2thetaW);
      delta_m_en = 
	mn/(4.*M_PI)*log(m_Lambda2/mn2) *
	(9./4.*a_mn2/m_cos2thetaW);

      ZULR[0][0] = mn+delta_m_Qn-mn*(3./4.*sqr(offdiag_un/(4.*M_PI*m_v))*log(m_Lambda2/mn2));
      ZULR[0][1] = offdiag_un;
      ZULR[1][0] = offdiag_un;
      ZULR[1][1] = -mn-delta_m_un+mn*(3./2.*sqr(offdiag_un/(4.*M_PI*m_v))*log(m_Lambda2/mn2));
     
      ZDLR[0][0] = mn+delta_m_Qn-mn*(3./4.*sqr(offdiag_dn/(4.*M_PI*m_v))*log(m_Lambda2/mn2));
      ZDLR[0][1] = offdiag_dn;
      ZDLR[1][0] = offdiag_dn;
      ZDLR[1][1] = -mn-delta_m_dn+mn*(3./2.*sqr(offdiag_dn/(4.*M_PI*m_v))*log(m_Lambda2/mn2));

      ZLLR[0][0] = mn+delta_m_Ln-mn*(3./4.*sqr(offdiag_en/(4.*M_PI*m_v))*log(m_Lambda2/mn2));
      ZLLR[0][1] = offdiag_en;
      ZLLR[1][0] = offdiag_en;
      ZLLR[1][1] = -mn-delta_m_en+mn*(3./2.*sqr(offdiag_en/(4.*M_PI*m_v))*log(m_Lambda2/mn2));
      
      for (int k=0;k<2;k++) {
	for (int l=0;l<2;l++) {
	  ZULR2[k][l]    = 0.;
	  ZDLR2[k][l]    = 0.;
	  ZLLR2[k][l]    = 0.;
	  for (int m=0;m<2;m++) {
	    ZULR2[k][l] += ZULR[k][m]*ZULR[m][l];
	    ZDLR2[k][l] += ZDLR[k][m]*ZDLR[m][l];
	    ZLLR2[k][l] += ZLLR[k][m]*ZLLR[m][l];
	  }
	}
      }  
      ZULR2.Diagonalize(EigenvaluesU,EigenvectorsU);
      ZDLR2.Diagonalize(EigenvaluesD,EigenvectorsD); 
      ZLLR2.Diagonalize(EigenvaluesL,EigenvectorsL);

      cout<<"("<<ZULR[0][0]<<", "<<ZULR[0][1]<<", "<<ZULR[1][0]<<", "<<ZULR[1][1]
	  <<") --> Q_n = "<<sqrt(EigenvaluesU[0])<<", u_n = "<<sqrt(EigenvaluesU[1])<<endl;
      cout<<"("<<ZDLR[0][0]<<", "<<ZDLR[0][1]<<", "<<ZDLR[1][0]<<", "<<ZDLR[1][1]
	  <<") --> Q_n = "<<sqrt(EigenvaluesD[0])<<", d_n = "<<sqrt(EigenvaluesD[1])<<endl;
      cout<<"("<<ZLLR[0][0]<<", "<<ZLLR[0][1]<<", "<<ZLLR[1][0]<<", "<<ZLLR[1][1]
	  <<") --> L_n = "<<sqrt(EigenvaluesL[0])<<", e_n = "<<sqrt(EigenvaluesL[1])<<endl;
      

      Flavour F1((kf_code)(int(kfnoL+i*100000+1+(j-1)*2))); F1.SetMass(sqrt(EigenvaluesD[0]));
      Flavour F2((kf_code)(int(kfnoR+i*100000+1+(j-1)*2))); F2.SetMass(sqrt(EigenvaluesD[1]));
      Flavour F3((kf_code)(int(kfnoL+i*100000+2+(j-1)*2))); F3.SetMass(sqrt(EigenvaluesU[0]));
      Flavour F4((kf_code)(int(kfnoR+i*100000+2+(j-1)*2))); F4.SetMass(sqrt(EigenvaluesU[1]));
      Flavour F5((kf_code)(int(kfnoL+i*100000+11+(j-1)*2))); F5.SetMass(sqrt(EigenvaluesL[0]));
      Flavour F6((kf_code)(int(kfnoR+i*100000+11+(j-1)*2))); F6.SetMass(sqrt(EigenvaluesL[1]));
      Flavour F7((kf_code)(int(kfnoL+i*100000+12+(j-1)*2))); F7.SetMass(sqrt(EigenvaluesL[0]));    
    }
    cout<<"Check for "<<i<<"th KK resonance"<<endl;
    if (i==1) {
      cout<<" d_{LR} :"<<Flavour(kf_KK_dL_1).Mass()
	  <<"/"<<Flavour(kf_KK_dR_1).Mass()
	  <<" s_{LR} :"<<Flavour(kf_KK_sL_1).Mass()
	  <<"/"<<Flavour(kf_KK_sR_1).Mass()
	  <<" b_{21} :"<<Flavour(kf_KK_b2_1).Mass()
	  <<"/"<<Flavour(kf_KK_b1_1).Mass()<<endl
	  <<" u_{LR} :"<<Flavour(kf_KK_uL_1).Mass()
	  <<"/"<<Flavour(kf_KK_uR_1).Mass()
	  <<" c_{LR} :"<<Flavour(kf_KK_cL_1).Mass()
	  <<"/"<<Flavour(kf_KK_cR_1).Mass()
	  <<" t_{21} :"<<Flavour(kf_KK_t2_1).Mass()
	  <<"/"<<Flavour(kf_KK_t1_1).Mass()<<endl
	  <<" nu_e :"<<Flavour(kf_KK_nueL_1).Mass()
	  <<" nu_mu :"<<Flavour(kf_KK_numuL_1).Mass()
	  <<" nu_tau :"<<Flavour(kf_KK_nutauL_1).Mass()<<endl
	  <<" e_{LR} :"<<Flavour(kf_KK_eL_1).Mass()
	  <<"/"<<Flavour(kf_KK_eR_1).Mass()
	  <<" mu_{LR} :"<<Flavour(kf_KK_muL_1).Mass()
	  <<"/"<<Flavour(kf_KK_muR_1).Mass()
	  <<" tau_{21} :"<<Flavour(kf_KK_tau2_1).Mass()
	  <<"/"<<Flavour(kf_KK_tau1_1).Mass()<<endl;
    }
    else if (i==2) {
      cout<<" d_{LR} :"<<Flavour(kf_KK_dL_2).Mass()
	  <<"/"<<Flavour(kf_KK_dR_2).Mass()
	  <<" s_{LR} :"<<Flavour(kf_KK_sL_2).Mass()
	  <<"/"<<Flavour(kf_KK_sR_2).Mass()
	  <<" b_{21} :"<<Flavour(kf_KK_b2_2).Mass()
	  <<"/"<<Flavour(kf_KK_b1_2).Mass()<<endl
	  <<" u_{LR} :"<<Flavour(kf_KK_uL_2).Mass()
	  <<"/"<<Flavour(kf_KK_uR_2).Mass()
	  <<" c_{LR} :"<<Flavour(kf_KK_cL_2).Mass()
	  <<"/"<<Flavour(kf_KK_cR_2).Mass()
	  <<" t_{21} :"<<Flavour(kf_KK_t2_2).Mass()
	  <<"/"<<Flavour(kf_KK_t1_2).Mass()<<endl
	  <<" nu_e :"<<Flavour(kf_KK_nueL_2).Mass()
	  <<" nu_mu :"<<Flavour(kf_KK_numuL_2).Mass()
	  <<" nu_tau :"<<Flavour(kf_KK_nutauL_2).Mass()<<endl
	  <<" e_{LR} :"<<Flavour(kf_KK_eL_2).Mass()
	  <<"/"<<Flavour(kf_KK_eR_2).Mass()
	  <<" mu_{LR} :"<<Flavour(kf_KK_muL_2).Mass()
	  <<"/"<<Flavour(kf_KK_muR_2).Mass()
	  <<" tau_{21} :"<<Flavour(kf_KK_tau2_2).Mass()
	  <<"/"<<Flavour(kf_KK_tau1_2).Mass()<<endl;
    }
  }
}



