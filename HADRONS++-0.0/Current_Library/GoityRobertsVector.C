#include "GoityRobertsVector.H"
#include "Tools.H"

using namespace HADRONS;
using namespace ATOOLS;


void GoityRobertsVector::SetModelParameters( struct GeneralModel _md )
{
  m_Vcb = _md("Vcb",0.04);
}


void GoityRobertsVector::Calc()
{
  Vec4D pB = p_moms[0];
  Vec4D pD = p_moms[1];
  Vec4D pPi= p_moms[2];
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
 
  for( int h_had=0; h_had<3; h_had++) {
    ComplexVec4D eps = Tools::ComplexBosonPolarizationVectorC(p_moms[1], h_had);
    Complex i = Complex(0.0,1.0);
    p_results[h_had] = m_Vcb*sqrt(mB*mD)*
                             (-1.0*i/2.0*(h1*mB*mD*cross(eps.Conjugate(),vB,vD) + 
                                     h2*mB*cross(eps.Conjugate(),vB,pPi) + 
                                     h3*mD*cross(eps.Conjugate(),vD,pPi)) + 
                              f1*mB*vB*(eps.Conjugate()*pPi) + f2*mD*vD*(eps.Conjugate()*pPi) + 
                              f3*pPi*(eps.Conjugate()*pPi) + f4*mB*mB*vB*(eps.Conjugate()*vB) +
                              f5*mB*mD*vD*(eps.Conjugate()*vB) + f6*mB*pPi*(eps.Conjugate()*vB) +
                              k*eps.Conjugate() - 
                              i/2.0*cross(vB,vD,pPi)*(g1*eps.Conjugate()*pPi+g2*mB*eps.Conjugate()*vB) -
                              i/2.0*eps.Conjugate()*cross(vB,vD,pPi)*(g3*mB*vB+g4*mD*vD+g5*pPi));
  }
}


DECLARE_GETTER(GoityRobertsVector_Getter, "GoityRobertsVector",
               Current_Base,Flavour_Info);

Current_Base* GoityRobertsVector_Getter::operator()(const Flavour_Info &parameters) const
{
  return new GoityRobertsVector(parameters.flavs, parameters.nout, parameters.indices, "GoityRobertsVector");
}

void GoityRobertsVector_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}


