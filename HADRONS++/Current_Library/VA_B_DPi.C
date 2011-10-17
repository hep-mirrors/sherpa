#include "HADRONS++/Current_Library/VA_B_DPi.H"

using namespace HADRONS;
using namespace ATOOLS;


void VA_B_DPi::SetModelParameters( struct GeneralModel _md )
{
  m_Vxx = _md("Vxx",Tools::Vcb);
}


void VA_B_DPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D pB = moms[p_i[0]];
  Vec4D pD = moms[p_i[1]];
  Vec4D pPi= moms[p_i[2]];
  double mB = p_masses[0];
  double mD = p_masses[1];

  Vec4D vB = pB/mB;         //4-velocity of B meson
  Vec4D vD = pD/mD;         //4-velocity of D
  double w = vB*vD;         //four velocity transfer
 
  Complex dmb = Complex(0.0460,-0.5*0.00001); 
  double g = 0.5;        
  double alpha3 =  0.690; // See table I in G&R's paper
  double alpha1 = -1.430;
  double alpha2 = -0.140;
  double f0=0.093;        // The pion decay constant set to 93 MeV

  Complex dmt3 = Complex (0.563,-0.5*0.191);  
  Complex dmt1 = Complex(0.392,-0.5*1.040);
  Complex dmt2 = Complex(0.709,-0.5*0.405);
                   
  double betas=0.285;      
  double betap=0.280;      
  double betad=0.260;      
  double betasp=betas*betas+betap*betap;
  double betasd=betas*betas+betad*betad;

  double lambdabar=0.750;  

  double xi = exp(lambdabar*lambdabar*(1.0-w*w)/(4*betas*betas));
  double xi1= -1.0*sqrt(2.0/3.0)*(lambdabar*lambdabar*(w*w-1.0)/(4*betas*betas))*
              exp(lambdabar*lambdabar*(1.0-w*w)/(4*betas*betas));
  double rho1= sqrt(1.0/2.0)*(lambdabar/betas)*
               pow((2*betas*betap/(betasp)),2.5)*
               exp(lambdabar*lambdabar*(1.0-w*w)/(2*betasp));
  double rho2= sqrt(1.0/8.0)*(lambdabar*lambdabar/(betas*betas))*
               pow((2*betas*betad/(betasd)),3.5)*
               exp(lambdabar*lambdabar*(1.0-w*w)/(2*betasd));

  Complex h,a1,a2,a3;
  Complex hnr,a1nr,a2nr,a3nr;
  Complex hr,a1r,a2r,a3r;

// Non-resonance part 
  hnr = g*xi*(1.0/(pPi*vB+dmb))/(2*f0*mB*mD);
  // in the original paper: hnr = g*xi*(1.0/(pPi*vB+dmb-i*eps)-1.0/(pPi*vD-dmd+i*eps))/(2*f0*mB*mD)
  a1nr= -1.0*g*xi*(1+w)*(1.0/(pPi*vB+dmb))/(2*f0);
  // a1nr= -1.0*g*xi*(1+w)*(1.0/(pPi*vB+dmb-i*eps)-1.0/(pPi*vD-dmd+i*eps))/(2*f0)
  a2nr= g*xi*((pPi*(vB+vD))/(pPi*vB+dmb))/(2*f0*mB);
  // a2nr= g*xi*((pPi*(vB+vD))/(pPi*vB+dmb-i*eps))/(2*f0*mB)
  a3nr= Complex(0.0,0.0);
  // a3nr= -1.0*g*xi*((pPi*(vB+vD))/(pPi*vD-dmd+i*eps))/(2*f0*mD)

// Resonance part 
  hr = alpha2*rho2*(w-1)*(1.0/(pPi*vB+dmt2))/(6*f0*mB*mD) +
       alpha3*xi1*(1.0/(pPi*vB+dmt3))/(2*f0*mB*mD);
  // hr = alpha2*rho2*(w-1)*(1.0/(pPi*vB+dmt2)-1.0/(pPi*vD-dmt2))/(6*f0*mB*mD) +
  //      alpha3*xi1*(1.0/(pPi*vB+dmt3)-1.0/(pPi*vD-dmt3))/(2*f0*mB*mD)
  a1r= -1.0*alpha2*rho2*(w*w-1)*(1.0/(pPi*vB+dmt2))/(6*f0) -
       alpha3*xi1*(1+w)*(1.0/(pPi*vB+dmt3))/(2*f0);
  // a1r= -1.0*alpha2*rho2*(w*w-1)*(1.0/(pPi*vB+dmt2)-1.0/(pPi*vD+dmt2))/(6*f0) -
  //      alpha3*xi1*(1+w)*(1.0/(pPi*vB+dmt3)-1.0/(pPi*vD-dmt3))/(2*f0)
  a2r= alpha1*rho1*((pPi*vB)/(pPi*vB+dmt1))/(2*f0*mB) +
       alpha2*rho2*(0.5*pPi*(w*vD-vB)+pPi*(vD-w*vB))/
                  (3*f0*mB*(pPi*vB+dmt2)) +
       alpha3*xi1*((pPi*(vB+vD))/(pPi*vB+dmt3))/(2*f0*mB);
  // a2r= alpha1*rho1*((pPi*vB)/(pPi*vD-dmt1)+(pPi*vB)/(pPi*vB+dmt1))/(2*f0*mB) +
  //      alpha2*rho2*(1/(pPi*vB+dmt2)*(1/6*(w*pPi*vD-pPi*vB)+1/3*(pPi*vD-*w*pPi*vB)) +
  //      0.5*1/(pPi*vD-dmt2)*(pPi*vB-w*pPi*vD))/(f0*mB) +
  //      alpha3*xi1*((pPi*(vB+vD))/(pPi*vB+dmt3))/(2*f0*mB)
  a3r= -1.0*alpha1*rho1*((pPi*vB)/(pPi*vB+dmt1))/(2*f0*mD) -
       alpha2*rho2*((pPi*(vD-w*vB))/(pPi*vB+dmt2))/(2*f0*mD);
  // a3r= -1.0*alpha1*rho1*((pPi*vB)/(pPi*vD-dmt1)+(pPi*vB)/(pPi*vB+dmt1))/(2*f0*mD) -
  //     alpha2*rho2*(1/(pPi*vD-dmt2)*(1/6*(w*pPi*vB-pPi*vD)+1/3*(pPi*vB-*w*pPi*vD)) +
  //      0.5*1/(pPi*vB+dmt2)*(pPi*vD-w*pPi*vB))/(2*f0*mD) - 
  //      alpha3*xi1*((pPi*(vB+vD))/(pPi*vD-dmt3))/(2*f0*mB)


// Sum
  h=hnr+hr;
  a1=a1nr+a1r;
  a2=a2nr+a2r;
  a3=a3nr+a3r;

  Insert( m_Vxx * sqrt(mB*mD) * (Complex(0.0,1.0)*h*mB*mD*cross(vB,vD,pPi)+
                                 a1*pPi+a2*mB*vB+a3*mD*vD) , 0);
}

DEFINE_CURRENT_GETTER(VA_B_DPi,VA_B_DPi_Getter,"VA_B_DPi")

void VA_B_DPi_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"\\paragraph{Parametrization} \n"
    <<"  \\[ \\mathcal{J^{B\\to D \\pi}_\\mu} = V_{cb} \\sqrt{M_B M_D} \\left( \n"
    <<"      -i h M_B M_D \\epsilon_{\\mu\\nu\\rho\\sigma}v_B^\\nu v_D^\\rho p_\\pi^\\sigma \n"
    <<"      + A_1 p_{\\pi\\mu} + A_2 M_B v_{B\\mu} + A_3 M_D v_{D\\mu} \\right) \\]\n"
    <<"  \\begin{itemize} \n"
    <<"    \\item Particle order: 0 = $B$, 1 = $D$, 2 = $\\pi$ \n"
    <<"    \\item Example: $ B \\to D \\; \\pi \\; l \\; \\nu_l $ \n"
    <<"    \\item Reference: Goity, Roberts \\cite{Goity:1994xn}; EvtGen\\cite{Lange:2001uf} \n"
    <<"  \\end{itemize} \n\n"
    <<"\\paragraph{Available parameters} \n"
    <<"  No parameters. \n"<<std::endl;
}
