#include "Virtual_Amp.H"
#include <cmath>
#include <iostream>
#include <string.h>
#include <ginac/ginac.h>
#include <cln/cln.h>

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "EXTAMP/External_ME_Interface.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "ATOOLS/Phys/NLO_Types.H"

#include "PHASIC++/Process/Process_Info.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"


namespace DiBosonsNLO{

Virtual_Amp::Virtual_Amp(const PHASIC::Process_Info& pi,
                         const ATOOLS::Flavour_Vector& flavs,
                         const bool off_shell) :
    PHASIC::Virtual_ME2_Base(pi, flavs),
    amp(5,IRschemeQt), F(5,IRschemeQt),
    m_off_shell(off_shell),
    m_current_1_type{0,0,0}, m_current_2_type{0,0,0},
    Z(23), W(24)
{
    DEBUG_FUNC(this);


    m_symfac =pi.m_fi.FSSymmetryFactor();
    m_symfac*=pi.m_ii.ISSymmetryFactor();

    m_flavs = flavs;

    MODEL::Coupling_Map cpls;
    MODEL::s_model->GetCouplings(cpls);
    SetCouplings(cpls);

    m_mode = 1;

    m_gluon_induced = flavs[0].IsGluon()&&flavs[1].IsGluon();

    m_beta0 =  11. - 10./3.;

    if(m_off_shell)
    {
      Get_Current_Types();
      Select_Processes();
    }else
    {
      m_process_type_on_shell = Process_Type();
    }
}

Virtual_Amp::~Virtual_Amp()
{
}


void Virtual_Amp::Get_Current_Types()
{


    // Fill m_current_type for one boson

    m_current_1_type[0] = (!m_flavs[2].IsNeutrino() && !m_flavs[4].IsNeutrino())?1:0;   // P
    m_current_1_type[1] = ((m_flavs[2].IsNeutrino() && m_flavs[4].IsNeutrino()) ||
        m_current_1_type[0]==1)?1:0;                                                    // Z
    m_current_1_type[2] = (m_flavs[2].IsNeutrino() != m_flavs[4].IsNeutrino())?1:0;     // W

    // Fill m_current_type for the other boson

    m_current_2_type[0] = (!m_flavs[3].IsNeutrino() && !m_flavs[5].IsNeutrino())?1:0;   // P
    m_current_2_type[1] = ((m_flavs[3].IsNeutrino() && m_flavs[5].IsNeutrino()) ||
        m_current_1_type[0]==1)?1:0;                                                    // Z
    m_current_2_type[2] = (m_flavs[3].IsNeutrino() != m_flavs[5].IsNeutrino())?1:0;     // W
}

void Virtual_Amp::Select_Processes()
{
    m_process_type[0] = m_current_1_type[0]*m_current_2_type[0];                        // PP
    if((m_current_1_type[1]*m_current_2_type[0]) ||
            (m_current_1_type[0]*m_current_2_type[1])) m_process_type[1] = 1;           // ZP || PZ
    m_process_type[2] = m_current_1_type[1]*m_current_2_type[1];                        // ZZ
    m_process_type[3] = m_current_1_type[2]*m_current_2_type[2];                        // WW
    for(int i=0;i<4;i++) std::cout<<m_process_type[i]<<std::endl;
}






void Virtual_Amp::Calc(const ATOOLS::Vec4D_Vector& p)
{
  THROW(fatal_error, "Invalid call");
}

double Virtual_Amp::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& p)
{
    /* They pull out a factor of 4\pi^\Epsilon/\Gamma(1-\Epsilon) in
     *        frot of the virtual amplitude. 1/\Gamma(1-\Epsilon) is alsways
     *               assumed, but we need to correct for the rest */
    return 4.*M_PI;
}



// THIS IS THE MAIN FUNCTION WHERE THE SQUARED ME IS CALCULATED.

void Virtual_Amp::Calc(const ATOOLS::Vec4D_Vector& P, const double& Born)
{

    // Here I've copied the momentum vector in order to use the same momentum labeling
    // used in the ggvvamp paper and thus ease the comparison


    alphaqed =  AlphaQED();
    alphas   =  AlphaQCD();
    double V=0.;

    if(m_off_shell) V = Calc_OffShell(P,Born);
    else V = Calc_OnShell(P,Born);

    m_res.Finite() = V;
    m_res.IR()     = 0.;
    m_res.IR2()    = 0.;


}


double Virtual_Amp::Calc_OffShell(const ATOOLS::Vec4D_Vector& P, const double& Born)
{
  m_P = P;

  m_P[3] = P[4];
  m_P[4] = P[3];

  double ma2 = (m_P[2]+m_P[3]).Abs2();
  double mb2 = (m_P[4]+m_P[5]).Abs2();
  double s   = (m_P[0]+m_P[1]).Abs2();
  double t   = (m_P[0]-m_P[2]-m_P[3]).Abs2();


  double V0(0.);
  double V(0.);
  double B_mless(0.);
  // m_gluon_induced is a check in case sherpa also look for 
  // di-loop quark initiated. They are not in ggvvamp so I set them 
  // to 0
  m_gluon_induced=1;
  if(m_gluon_induced)
  {
      amp.compute(s, t, ma2, mb2);
      get_form_factors();

      calc_hel_amp();


      for(int i=0;i<4;i++)
          for(int j=0;j<4;j++)
          {
              V0+=hel_squared(i,j);
              B_mless+=born_hel_squared(i,j);
          }

      V0 = Born*V0/B_mless;
      V =  V0 + Born*(m_beta0*log(m_mur2/s));
  }else V=0;
  return V;
}


double Virtual_Amp::Calc_OnShell(const ATOOLS::Vec4D_Vector& P, const double& Born)
{
  double M32  = sqr(m_flavs[2].Mass());
  double M42  = sqr(m_flavs[3].Mass());
  double s    = (P[0]+P[1]).Abs2();
  double t    = (P[0]-P[2]).Abs2();
  double u    = M32 + M42 - s - t;
  double b3, b4, b3s, b4s, bf;

  b3  = 1 - t/M32;
  b4  = 1 - t/M42;
  b3s = 1 - b3*M32/s;
  b4s = 1 - b4*M42/s;
  bf  = (s - M32 - M42)/2.;

  int reim;
  Form_Factors_Computation(P);
  #include"Vfin"
  #include"Born"

  double V0 =2.*Vfin*EWCouplings(m_process_type_on_shell);
  B_mless = B_mless*EWCouplings(m_process_type_on_shell);

  V0 = Born*V0/B_mless;

  double V = V0 + m_beta0*Born*log(m_mur2/s);

  return V;
}



void Virtual_Amp::get_form_factors()
{
    for(int i=0; i<m_nHel; i++)
        for(int j=0; j<(m_nCoeff); j++)
            for(int k=0; k<m_nLoop; k++)
                for(int m=0; m<m_nreim; m++)
                    m_E[i][j][k][m]=amp.E[i][j+1][k][m];
}



double Virtual_Amp::EWCouplings(int process_type)
{

  double I3[6] = {-0.5, 0.5, -0.5, 0.5, -0.5, 0.5};
  double q[6] = {-1., 2., -1., 2., -1., 2.};

  std::complex<double> N(0.,0.);


  //-----------------------
  //    process_type(s):  |
  //      -- 1    yy      |
  //      -- 2    Zy      |
  //      -- 3    ZZ      |
  //      -- 4    WW      |
  //-----------------------

  // Only 1 and 3 implemented so far in the onshell!!
  // offshell has them all


  if(process_type==1){
      for(size_t i=0; i<5; i++){N += q[i]*q[i];}
      N /= 9.;
  }


  if(process_type==3){
      for(size_t i=0; i<5; i++){N += q[i]*q[i]*m_sin2w*m_sin2w - 3.*q[i]*I3[i]*m_sin2w;}
      // 13.5 is the number for 6 flavours in the loops, while for 5 is 11.25
      N = (2.*N + 11.25)/(18.*m_sin2w*m_cos2w);
  }

  // Building the prefactor
  double fac = 1.;
  // EW structure
  fac *= std::norm(N);
  fac *= 16.*M_PI*M_PI*alphaqed*alphaqed;
  // QCD perturbative order definition
  fac *= alphas*alphas/(4.*M_PI*M_PI);
  // Gluon polarizations average/states
  fac /= 4.*64.;
  // Color structure
  fac *= 8.;
  // Identical final states
  fac /= 2.;

  return fac;

}



int Virtual_Amp::Process_Type(){
    // Photon Photon
    if(m_flavs[2].Kfcode()==22 && m_flavs[3].Kfcode()==22) {return 1;}

    // Photon Z
    if( (m_flavs[2].Kfcode()==23 && m_flavs[3].Kfcode()==22)
     || (m_flavs[2].Kfcode()==22 && m_flavs[3].Kfcode()==23) ) {return 2;}

    // Z Z
    if(m_flavs[2].Kfcode()==23 && m_flavs[3].Kfcode()==23) {return 3;}

    // W W
    if(m_flavs[2].Kfcode()==24 || m_flavs[3].Kfcode()==24) {return 4;}

    // Else False
    else {return 0;}
}





void Virtual_Amp::Form_Factors_Computation(const ATOOLS::Vec4D_Vector& p){

    double M32  = ATOOLS::sqr(m_flavs[2].Mass());
    double M42  = ATOOLS::sqr(m_flavs[3].Mass());

    double s    = (p[0]+p[1]).Abs2();
    double t    = (p[0]-p[2]).Abs2();
    double u    = M32 + M42 - s - t;
    F.compute(s, t, M32, M42);
  }



double Virtual_Amp::hel_squared(int i, int j)
{
    return 2.*Couplings()*Sym_Factor()*(M[0][i][j]*std::conj(M[1][i][j])).real();
}

double Virtual_Amp::born_hel_squared(int i, int j)
{
    return Couplings()*Sym_Factor()*(M[0][i][j]*std::conj(M[0][i][j])).real();
}

void Virtual_Amp::calc_hel_amp()
{
    char hs[]={'-','+'};
    char hs_conj[]={'+','-'};

    for(int i=0;i<2;i++){
      // --LL
      kin_factors_calc(hs, 0, 1, 2, 3, 4, 5);
      M[i][0][0] = Current_Factor(0)*calc_hel_amp(0,i); // LL--
      M[i][1][0] = Current_Factor(0)*calc_hel_amp(1,i); // LR--
      kin_factors_calc(hs_conj, 0, 1, 3, 2, 5, 4);
      M[i][2][0] = Current_Factor(0)*calc_hel_amp(0,i); // RR--
      M[i][3][0] = Current_Factor(0)*calc_hel_amp(1,i); // RL--

      // --LR
      kin_factors_calc(hs, 0, 1, 2, 3, 5, 4);
      M[i][0][1] = Current_Factor(1)*calc_hel_amp(0,i);
      M[i][1][1] = Current_Factor(1)*calc_hel_amp(1,i);
      kin_factors_calc(hs_conj, 0, 1, 3, 2, 4, 5);
      M[i][2][1] = Current_Factor(1)*calc_hel_amp(0,i);
      M[i][3][1] = Current_Factor(1)*calc_hel_amp(1,i);

      // --RL
      kin_factors_calc(hs, 0, 1, 3, 2, 4, 5);
      M[i][0][2] = Current_Factor(2)*calc_hel_amp(0,i);
      M[i][1][2] = Current_Factor(2)*calc_hel_amp(1,i);
      kin_factors_calc(hs_conj, 0, 1, 2, 3, 5, 4);
      M[i][2][2] = Current_Factor(2)*calc_hel_amp(0,i);
      M[i][3][2] = Current_Factor(2)*calc_hel_amp(1,i);

      // --RR
      kin_factors_calc(hs, 0, 1, 3, 2, 5, 4);
      M[i][0][3] = Current_Factor(3)*calc_hel_amp(0,i);
      M[i][1][3] = Current_Factor(3)*calc_hel_amp(1,i);
      kin_factors_calc(hs_conj, 0, 1, 2, 3, 4, 5);
      M[i][2][3] = Current_Factor(3)*calc_hel_amp(0,i);
      M[i][3][3] = Current_Factor(3)*calc_hel_amp(1,i);
    }

}


// Angular brackets are hs[0], square brackets are hs[1]
void Virtual_Amp::kin_factors_calc(const char hs[], int i1, int i2, int i5, int i6, int i7, int i8)
{
    std::complex<double> factor = hm(hs[1],i2,i5)*hm(hs[0],i5,i1) + hm(hs[1],i2,i6)*hm(hs[0],i6,i1);

    std::complex<double> C[2];
    C[0] = hm(hs[0],i1,i2)*( hm(hs[1],i1,i5)*hm(hs[0],i5,i2) + hm(hs[1],i1,i6)*hm(hs[0],i6,i2) )/hm(hs[1],i1,i2);
    C[1] = hm(hs[1],i2,i5)*hm(hs[0],i5,i1) + hm(hs[1],i2,i6)*hm(hs[0],i6,i1);

    for(int i = 0; i<2; i++)
    {
        kin_factors[i][0]=( C[i]*factor*hm(hs[0],i5,i7)*hm(hs[1],i6,i8) );
        kin_factors[i][1]=( C[i]*factor*hm(hs[0],i1,i5)*hm(hs[0],i1,i7)*hm(hs[1],i1,i6)*hm(hs[1],i1,i8) );
        kin_factors[i][2]=( C[i]*factor*hm(hs[0],i1,i5)*hm(hs[0],i2,i7)*hm(hs[1],i1,i6)*hm(hs[1],i2,i8) );
        kin_factors[i][3]=( C[i]*factor*hm(hs[0],i2,i5)*hm(hs[0],i1,i7)*hm(hs[1],i2,i6)*hm(hs[1],i1,i8) );
        kin_factors[i][4]=( C[i]*factor*hm(hs[0],i2,i5)*hm(hs[0],i2,i7)*hm(hs[1],i2,i6)*hm(hs[1],i2,i8) );
        kin_factors[i][5]=( C[i]*hm(hs[0],i1,i5)*hm(hs[0],i1,i7)*hm(hs[1],i1,i6)*hm(hs[1],i2,i8) );
        kin_factors[i][6]=( C[i]*hm(hs[0],i1,i5)*hm(hs[0],i1,i7)*hm(hs[1],i2,i6)*hm(hs[1],i1,i8) );
        kin_factors[i][7]=( C[i]*hm(hs[0],i1,i5)*hm(hs[0],i2,i7)*hm(hs[1],i2,i6)*hm(hs[1],i2,i8) );
        kin_factors[i][8]=( C[i]*hm(hs[0],i2,i5)*hm(hs[0],i1,i7)*hm(hs[1],i2,i6)*hm(hs[1],i2,i8) );
    }
}

// hm stays for helicity momentum
std::complex<double> Virtual_Amp::hm(char hel, int i, int j)
{
    ATOOLS::Vec4D pi = hel_to_mom(i);
    ATOOLS::Vec4D pj = hel_to_mom(j);

    double pi_plus  = pi[0] + pi[3];
    double pi_minus = pi[0] - pi[3];
    double pj_plus  = pj[0] + pj[3];
    double pj_minus = pj[0] - pj[3];
    std::complex<double> pi_T(pi[1],pi[2]);
    std::complex<double> pj_T(pj[1],pj[2]);

    std::complex<double> phase_i=1.;
    std::complex<double> phase_j=1.;

    if((std::abs(pi_plus)>pi[0]*1.e-12)&&(std::abs(pi[1])>pi[0]*1.e-12 || std::abs(pi[2])>pi[0]*1.e-12)) phase_i = pi_T/sqrt(pi[1]*pi[1] + pi[2]*pi[2]);
    if((std::abs(pj_plus)>pj[0]*1.e-12)&&(std::abs(pj[1])>pj[0]*1.e-12 || std::abs(pj[2])>pj[0]*1.e-12)) phase_j = pj_T/sqrt(pj[1]*pj[1] + pj[2]*pj[2]);

    std::complex<double> h_mom(0.,0.);
    h_mom = sqrt(pi_plus*pj_minus)*phase_j - sqrt(pj_plus*pi_minus)*phase_i;


    if(hel == '+') h_mom = (-1.)*std::conj(h_mom);

    return h_mom;
}





std::complex<double> Virtual_Amp::calc_hel_amp(const int hel, const int nLoop)
{
    std::complex<double> hel_amp(0.,0.);
    for(int i=0; i<m_nCoeff; i++){
        hel_amp += kin_factors[hel][i]*(m_E[hel][i][nLoop][0] + I*m_E[hel][i][nLoop][1]);
    }

    return hel_amp;
}




std::complex<double> Virtual_Amp::Current_Factor(int helicity)
{
    std::complex<double> NJ(0.,0.);
    std::complex<double> JZ1;
    std::complex<double> JZ2;
    std::complex<double> JW1;
    std::complex<double> JW2;
    std::complex<double> Jy1;
    std::complex<double> Jy2;

    // --LL
    if(helicity==0)
    {
        JZ1 = LJZ(m_P[2]+m_P[3],0);
        JW1 = LJW(m_P[2]+m_P[3],0);
        Jy1 = LJy(m_P[2]+m_P[3],0);
        JZ2 = LJZ(m_P[4]+m_P[5],1);
        JW2 = LJW(m_P[4]+m_P[5],1);
        Jy2 = LJy(m_P[4]+m_P[5],1);
    }

    // --LR
    if(helicity==1)
    {
        JZ1 = LJZ(m_P[2]+m_P[3],0);
        JW1 = LJW(m_P[2]+m_P[3],0);
        Jy1 = LJy(m_P[2]+m_P[3],0);
        JZ2 = RJZ(m_P[4]+m_P[5],1);
        JW2 = 0.;
        Jy2 = RJy(m_P[4]+m_P[5],1);
    }

    // --RL
    if(helicity==2)
    {
        JZ1 = RJZ(m_P[2]+m_P[3],0);
        JW1 = 0.;
        Jy1 = RJy(m_P[2]+m_P[3],0);
        JZ2 = LJZ(m_P[4]+m_P[5],1);
        JW2 = LJW(m_P[4]+m_P[5],1);
        Jy2 = LJy(m_P[4]+m_P[5],1);
    }

    // --RR
    if(helicity==3)
    {
        JZ1 = RJZ(m_P[2]+m_P[3],0);
        JW1 = 0.;
        Jy1 = RJy(m_P[2]+m_P[3],0);
        JZ2 = RJZ(m_P[4]+m_P[5],1);
        JW2 = 0.;
        Jy2 = RJy(m_P[4]+m_P[5],1);
    }


    double I[6] = {-0.5, 0.5, -0.5, 0.5, -0.5, 0.5};
    double q[6] = {-1., 2., -1., 2., -1., 2.};
    std::complex<double> N(0.,0.);

    if(m_process_type[0])
    {
        for(size_t i=0; i<5; i++){N += q[i]*q[i];}
        N /= 9.;
        NJ += N*Jy1*Jy2;
        N = 0.;
    }

    if(m_process_type[1])
    {
        for(size_t i=0; i<5; i++){N += q[i]*(q[i]*m_sin2w/3. - I[i]/2.)/3.;}
        N /= sqrt(m_sin2w*m_cos2w);
        NJ += N*(JZ1*Jy2 + Jy1*JZ2);
        N = 0.;
    }

    if(m_process_type[2])
    {
        for(size_t i=0; i<5; i++){N += 9.*I[i]*I[i] + 2.*q[i]*q[i]*m_sin2w*m_sin2w - 6.*q[i]*I[i]*m_sin2w;}
        N /= 18.*m_sin2w*m_cos2w;
        NJ += N*JZ1*JZ2;
        N = 0.;
    }

    if(m_process_type[3])
    {
        N += 2./(4.*m_sin2w);
        NJ += N*JW1*JW2;
        N = 0.;
    }
    return NJ;

}

std::complex<double> Virtual_Amp::LJy(const ATOOLS::Vec4D& Q, const int which_boson)
{
    // which_boson identify whether this current is attached to the first boson or the second.
    // It can take only two values, 0 and 1, respetively for the first and the second.

    std::complex<double> Current;
    if(which_boson==0) Current = -1./Denominator(Q,0,0);
    else if(which_boson==1) Current = -1./Denominator(Q,0,0);
    return Current;
}

std::complex<double> Virtual_Amp::RJy(const ATOOLS::Vec4D& Q, const int which_boson)
{
    // which_boson identify whether this current is attached to the first boson or the second.
    // It can take only two values, 0 and 1, respetively for the first and the second.

    std::complex<double> Current;
    if(which_boson==0) Current = -1./Denominator(Q,0,0);
    else if(which_boson==1) Current = -1./Denominator(Q,0,0);
    return Current;
}

std::complex<double> Virtual_Amp::LJZ(const ATOOLS::Vec4D& Q, const int which_boson)
{
    // which_boson identify whether this current is attached to the first boson or the second.
    // It can take only two values, 0 and 1, respetively for the first and the second.

    std::complex<double> Current(0.,0.);
    if(which_boson==0) Current = ((m_flavs[2].IsoWeak() - m_sin2w*m_flavs[2].Charge())/sqrt(m_sin2w*m_cos2w))/Denominator(Q,m_MZ,m_WZ);
    else if(which_boson==1) Current = ((m_flavs[3].IsoWeak() - m_sin2w*m_flavs[3].Charge())/sqrt(m_sin2w*m_cos2w))/Denominator(Q,m_MZ,m_WZ);
    return Current;
}

std::complex<double> Virtual_Amp::RJZ(const ATOOLS::Vec4D& Q, const int which_boson)
{
    // which_boson identify whether this current is attached to the first boson or the second.
    // It can take only two values, 0 and 1, respetively for the first and the second.

    std::complex<double> Current;
    if(which_boson==0) Current = ((- m_sin2w*m_flavs[2].Charge())/sqrt(m_sin2w*m_cos2w))/Denominator(Q,m_MZ,m_WZ);
    else if(which_boson==1) Current = ((- m_sin2w*m_flavs[3].Charge())/sqrt(m_sin2w*m_cos2w))/Denominator(Q,m_MZ,m_WZ);
    return Current;
}

std::complex<double> Virtual_Amp::LJW(const ATOOLS::Vec4D& Q, const int which_boson)
{
    // which_boson identify whether this current is attached to the first boson or the second.
    // It can take only two values, 0 and 1, respetively for the first and the second.

    std::complex<double> Current;
    if(which_boson==0) Current = ( 1./sqrt(2)/sqrt(m_sin2w) )/Denominator(Q,m_MW,m_WW);
    else if(which_boson==1) Current = ( 1./sqrt(2)/sqrt(m_sin2w) )/Denominator(Q,m_MW,m_WW);
    return Current;
}

// Right hand side W current is 0

double Virtual_Amp::Couplings()
{
    double coup=1.;
    // QED couplings
    coup *= 16.*M_PI*M_PI*alphaqed*alphaqed;
    // it is squared for one side is attached to the loop and the other to the final state fermions
    coup *= coup;
    // QCD perturbative order definition
    coup *= alphas*alphas/(4.*M_PI*M_PI);
    return coup;
}

double Virtual_Amp::Sym_Factor()
{
    double sym =1.;
    // gluon polarization and color average
    sym /= 4.*64.;
    // Color factor. Since the final state is colorless the initial can only be one kind of gluon
    // hence the factor is 8
    sym *= 8.;
    return sym;
}


std::complex<double> Virtual_Amp::Denominator(const ATOOLS::Vec4D& Q, double M, double W)
{
    return ((Q.Abs2() - ATOOLS::sqr(M)) + I*W*M);
}


  bool Virtual_Amp::IsMappableTo(const PHASIC::Process_Info& pi){
    return false;
  }
}








DECLARE_VIRTUALME2_GETTER(DiBosonsNLO::Virtual_Amp,
		       "Virtual_Amp")

PHASIC::Virtual_ME2_Base* ATOOLS::Getter<PHASIC::Virtual_ME2_Base,
			      PHASIC::Process_Info,
			      DiBosonsNLO::Virtual_Amp>::
operator()(const PHASIC::Process_Info &pi) const
{
  if (pi.m_loopgenerator!="ggvvamp") return NULL;
  /* Check NLO type (allow only QCD, not EW)  */
  if (pi.m_fi.m_nlotype!=ATOOLS::nlo_type::loop ||
      (pi.m_fi.m_nlocpl[0]!=1 && pi.m_fi.m_nlocpl[1]!=0)) return NULL;
  if(pi.m_maxcpl[0]!=3) return NULL;

  ATOOLS::Flavour_Vector flavs = pi.ExtractFlavours();
  bool fs_on_shell=false, fs_off_shell=false;
  if(flavs.size()==6) fs_off_shell=true;
  if(flavs.size()==4) fs_on_shell=true;
  std::cout<<fs_on_shell<<' '<<fs_off_shell<<std::endl;
  if(!(fs_on_shell^fs_off_shell)) return NULL;

  //if(!flavs[0].IsGluon() && !flavs[1].IsGluon()) return NULL;

  return new DiBosonsNLO::Virtual_Amp(pi, flavs, fs_off_shell);
}
