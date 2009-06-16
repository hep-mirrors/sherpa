#include "HADRONS++/Current_Library/VA_B_DstarPi.H"
#include "HELICITIES/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace HELICITIES;


void VA_B_DstarPi::SetModelParameters( struct GeneralModel _md )
{
  m_Vxx = _md("Vxx",Tools::Vcb);
}


void VA_B_DstarPi::Calc(const ATOOLS::Vec4D * moms)
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
  Complex dmd = Complex(0.1421,-0.5*0.00006);
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

  Complex h1nr,h2nr,h3nr,f1nr,f2nr;
  Complex f3nr,f4nr,f5nr,f6nr,knr,g1nr,g2nr,g3nr,g4nr,g5nr;
  Complex h1r,h2r,h3r,f1r,f2r,f3r,f4r,f5r,f6r,kr,g1r,g2r,g3r,g4r,g5r;
  Complex h1,h2,h3,f1,f2,f3,f4,f5,f6,k,g1,g2,g3,g4,g5;

// Non-resonance part
  h1nr = -g*xi*(pPi*vB)/(f0*mB*mD*(pPi*vB+dmb));
  // h1nr = -g*xi*(pPi*vB)/(f0*mB*mD*(pPi*vB+dmb-i*eps))
  h2nr = -g*xi/(f0*mB*(pPi*vB+dmb));
  // h2nr = -g*xi/(f0*mB*(pPi*vB+dmb-i*eps))
  h3nr = -(g*xi/(f0*mD))*(1.0/(pPi*vB+dmb)
			  -(1.0+w)/(pPi*vD));
  // h3nr = -(g*xi/(f0*mD))*(1.0/(pPi*vB+dmb-i*eps)
  //			  -(1.0+w)/(pPi*vD+i*eps)) 
  f1nr = -(g*xi/(2*f0*mB))*(1.0/(pPi*vB+dmb) -
         1.0/(pPi*vD+dmd));
  // f1nr = -(g*xi/(2*f0*mB))*(1.0/(pPi*vB+dmb-i*eps) -
  //        1.0/(pPi*vD+dmd+i*eps));
  f2nr = f1nr*mB/mD;
  f3nr = Complex(0.0,0.0);
  f4nr = Complex(0.0,0.0);
  f5nr = (g*xi/(2*f0*mB*mD))*(Complex(1.0,0.0)
                 +(pPi*vB)/(pPi*vB+dmb));
  // f5nr = (g*xi/(2*f0*mB*mD))*(Complex(1.0,0.0)
  //                +(pPi*vB)/(pPi*vB+dmb+i*eps));
  f6nr = (g*xi/(2*f0*mB))*(1.0/(pPi*vB+dmb)
			   -(1.0/(pPi*vD)));
  // f6nr = (g*xi/(2*f0*mB))*(1.0/(pPi*vB+dmb-i*eps)
  //			   -(1.0/(pPi*vD+i*eps)));
  knr = (g*xi/(2*f0))*((pPi*(vD-w*vB))/(pPi*vB+dmb) +
                 (pPi*(vB-w*vD))/(pPi*vD));
  // knr = (g*xi/(2*f0))*((pPi*(vD-w*vB))/(pPi*vB+dmb-i*eps) +
  //               (pPi*(vB-w*vD))/(pPi*vD+i*eps));
  g1nr = Complex(0.0,0.0);
  g2nr = Complex(0.0,0.0);
  g3nr = Complex(0.0,0.0);
  g4nr = (g*xi)/(f0*mD*(pPi*vD));
  // g4nr = (g*xi)/(f0*mD*(pPi*vD+i*eps))
  g5nr = Complex(0.0,0.0);

// Resonance part (D** removed by hand - alainb)
  h1r = -alpha1*rho1*(pPi*vB)/(f0*mB*mD*(pPi*vB+dmt1)) +
        alpha2*rho2*(pPi*(vB+2.0*w*vB-vD))
        /(3*f0*mB*mD*(pPi*vB+dmt2)) -
        alpha3*xi1*pPi*vB/(f0*mB*mD*pPi*vB+dmt3); 
  // h1r = -alpha1*rho1*(pPi*vB)/(f0*mB*mD*(pPi*vB+dmt1)-(pPi*vD)/(f0*mB*mD*(pPi*vD-dmt1)) +
  //      alpha2*rho2*((pPi*(vB+2.0*w*vB-vD))
  //      /(3*f0*mB*mD*(pPi*vB+dmt2)+ 3/2*(pPi*(w*vD-vB))
  //      /(3*f0*mB*mD*(pPi*vD-dmt2))) -
  //      alpha3*xi1*pPi*vB/(f0*mB*mD*pPi*vB+dmt3
  h2r = -alpha2*(1+w)*rho2/(3*f0*mB*(pPi*vB+dmt2)) -
        alpha3*xi1/(f0*mB*(pPi*vB+dmt3));
  h3r = alpha2*rho2*(1+w)/(3*f0*mD*(pPi*vB+dmt2)) -
        alpha3*xi1/(f0*mD*(pPi*vB+dmt3));
  //  h3r = alpha2*rho2/(3*f0*mD)*((1+w)/(pPi*vB+dmt2)-(w*w-1)/(2*(pPi*vD-dmt2))) -
  //      alpha3*xi1/(f0*mD)*(1/(pPi*vB+dmt3)-(1+w)/(pPi*vD-dmt3) )
  f1r = -alpha2*rho2*(w-1.0)/(6*f0*mB*(pPi*vB+dmt2)) -
        alpha3*xi1/(2*f0*mB*(pPi*vB+dmt3));
  // f1r = -alpha2*rho2*(w-1.0)/(6*f0*mB)*(1/(pPi*vB+dmt2)-1/(pPi*vD-dmt2)) -
  //      alpha3*xi1/(2*f0*mB)*(1/(pPi*vB+dmt3)-1/(pPi*vD-dmt3));

  f2r = f1r*mB/mD;
  f3r = Complex(0.0,0.0);
  f4r = Complex(0.0,0.0);
  f5r = alpha1*rho1*(pPi*vB)/(2*f0*mB*mD*(pPi*vB+dmt1)) +
        alpha2*rho2*(pPi*(vD-vB/3.0-2.0/3.0*w*vB))/
        (2*f0*mB*mD*(pPi*vB+dmt2)) +
        alpha3*xi1*(pPi*vB)/(2*f0*mB*mD*(pPi*vB+dmt3));
  // f5r = alpha1*rho1*(pPi*vB)*(2*f0*mB*mD**((pPi*vB)/(pPi*vB+dmt1)-(pPi*vD)/(pPi*vD-dmt1)) +
  //      alpha2*rho2/
  //      (2*f0*mB*mD)*(pPi*(vD-vB/3.0-2.0/3.0*w*vB)/(pPi*vB+dmt2)-pPi*(vB-vD/3.0-2.0/3.0*w*vD)/(pPi*vD-dmt2)) +
  //      alpha3*xi1*(pPi*vB)/(2*f0*mB*mD*(pPi*vB+dmt3))
  f6r = alpha2*rho2*(w-1.0)/(6*f0*mB*(pPi*vB+dmt2)) +
        alpha3*xi1/(2*f0*mB*(pPi*vB+dmt3));
  // f6r = alpha2*rho2/(6*f0*mB)*(w-1.0)*(1/(pPi*vB+dmt2)-1/(pPi*vD-dmt2)) +
  //       alpha3*xi1/(2*f0*mB)*(1/(pPi*vB+dmt3)-1/(pPi*vD-dmt3)) 
  kr = -alpha1*rho1*(w-1.0)*(pPi*vB)/(2*f0*(pPi*vB+dmt1)) -
       alpha2*rho2*(w-1.0)*(pPi*(vD-w*vB))
       /(3*f0*(pPi*vB+dmt2)) +
       alpha3*xi1*(pPi*(vD-w*vB))/(2*f0*(pPi*vB+dmt3));
  // kr = -alpha1*rho1*(w-1.0)/(2*f0)*((pPi*vB)/(pPi*vB+dmt1)+(pPi*vD)/(pPi*vD-dmt1)) -
  //     alpha2*rho2*(w-1.0)/(3*f0)*((pPi*(vD-w*vB))
  //     /(pPi*vB+dmt2)+(pPi*(vB-w*vD))
  //     /(pPi*vD-dmt2) ) +
  //     alpha3*xi1)/(2*f0)*((pPi*(vD-w*vB)/(pPi*vB+dmt3)-(pPi*(vB-w*vD)/(pPi*vD-dmt3))
  
  g1r = Complex(0.0,0.0);
  g2r = Complex(0.0,0.0);
  g3r = -g2r;
  g4r = 2.0*alpha2*rho2/(3*f0*mD*(pPi*vB+dmt2));
  // g4r = alpha2*rho2/(3*f0*mD)*(2/(pPi*vB+dmt2)-(1+0.5*w)*(1/(pPI+vD-dmt2))+alpha3*xi1/(f0*(pPi+vD-dmt3))
  g5r = Complex(0.0,0.0);

  //Sum
  h1 = h1nr + h1r;
  h2 = h2nr + h2r;
  h3 = h3nr + h3r;

  f1 = f1nr + f1r;
  f2 = f2nr + f2r;
  f3 = f3nr + f3r;
  f4 = f4nr + f4r;
  f5 = f5nr + f5r;
  f6 = f6nr + f6r;

  k = knr+kr;
  
  g1 = g1nr + g1r;
  g2 = g2nr + g2r;
  g3 = g3nr + g3r;
  g4 = g4nr + g4r;
  g5 = g5nr + g5r;
 
  Polarization_Vector pol(pD, sqr(mD));
  for( int h_had=0; h_had<3; h_had++) {
    Vec4C eps = pol[h_had];
    Complex i = Complex(0.0,1.0);
    Insert( m_Vxx*sqrt(mB*mD)*
            (-1.0*i/2.0*(h1*mB*mD*cross(conj(eps),vB,vD) +
                         h2*mB*cross(conj(eps),vB,pPi) +
                         h3*mD*cross(conj(eps),vD,pPi)) +
             f1*mB*vB*(conj(eps)*pPi) + f2*mD*vD*(conj(eps)*pPi) +
             f3*pPi*(conj(eps)*pPi) + f4*mB*mB*vB*(conj(eps)*vB) +
             f5*mB*mD*vD*(conj(eps)*vB) + f6*mB*pPi*(conj(eps)*vB) +
             k*conj(eps) -
             i/2.0*cross(vB,vD,pPi)*(g1*conj(eps)*pPi+g2*mB*conj(eps)*vB) -
             i/2.0*conj(eps)*cross(vB,vD,pPi)*(g3*mB*vB+g4*mD*vD+g5*pPi)) , h_had);
  }
}

DEFINE_CURRENT_GETTER(VA_B_DstarPi,VA_B_DstarPi_Getter,"VA_B_DstarPi")

void VA_B_DstarPi_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"\\paragraph{Parametrization} \n"
    <<"  \\[ \\mathcal{J^{B\\to D \\pi}_\\mu} = V_{cb} \\sqrt{M_B M_D} \\left( \n"
    <<"      \\frac{i}{2} \\epsilon_{\\mu\\nu\\rho\\sigma} \\varepsilon_{D\\nu} \\left[ \n"
    <<"      h_1 M_B M_D v_B^\\rho v_D^\\sigma + h_2 M_B v_B^\\rho p_\\pi^\\sigma +\n"
    <<"      h_3 M_D v_D^\\rho p_\\pi^\\sigma \\right] + \n"
    <<"      f_1 M_B v_{B\\mu} (\\varepsilon_D \\cdot p_\\pi) + \n"
    <<"      f_2 M_D v_{D\\mu} (\\varepsilon_D \\cdot p_\\pi) + \n"
    <<"      f_3 p_{\\pi\\mu} (\\varepsilon_D \\cdot p_\\pi) + \n"
    <<"      f_4 M_B^2 v_{B\\mu} (\\varepsilon_D \\cdot v_B) + \n"
    <<"      f_5 M_B M_D v_{D\\mu} (\\varepsilon_D \\cdot v_B) + \n"
    <<"      f_6 M_B p_{\\pi\\mu} (\\varepsilon_D \\cdot v_B) + \n"
    <<"      k g_{\\mu\\nu} + \n"
    <<"      \\frac{i}{2} \\epsilon_{\\mu\\delta\\rho\\sigma} v_B^\\delta v_D^\\rho p_\\pi^\\sigma \n "
    <<"        \\left(g_1 \\varepsilon_D \\cdot p_\\pi + g_2 M_B \\varepsilon_D \\cdot v_B\\right) \n"
    <<"       \\right) \\]\n"
    <<"  \\begin{itemize} \n"
    <<"    \\item Particle order: 0 = $B$, 1 = $D^*$, 2 = $\\pi$ \n"
    <<"    \\item Example: $ B \\rightarrow D^* \\; \\pi \\; l \\; \\nu_l $ \n"
    <<"    \\item Reference: Goity, Roberts \\cite{Goity:1994xn}; EvtGen\\cite{Lange:2001uf} \n"
    <<"  \\end{itemize} \n\n"
    <<"\\paragraph{Available parameters} \n"
    <<"  No parameters. \n"<<std::endl;
}
