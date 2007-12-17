#include "Interaction_Model_sLepton_sQuark.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sLepton_sQuark::Interaction_Model_sLepton_sQuark(MODEL::Model_Base * _model,
							   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),Ecms2)));
  g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(string("sin2_thetaW"))));
  sintW    = Kabbala(string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(string("sin2_thetaW"))));
  costW    = Kabbala(string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(string("sin2_thetaW"))));
  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  vev      = Kabbala(string("v_{EW}"),ScalarConstant(string("vev")));
  v1       = Kabbala(string("v_1"),ScalarConstant(string("vev")) *
		     sqrt(1./(1.+sqr(ScalarConstant(string("tan(beta)"))))));
  v2       = Kabbala(string("v_2"),ScalarConstant(string("vev")) *
		     ScalarConstant(string("tan(beta)")) *
		     sqrt(1./(1.+sqr(ScalarConstant(string("tan(beta)"))))));
  mu       = Kabbala(string("h"),ScalarConstant(string("mu")));
  conj_mu  = Kabbala(string("h"),ScalarConstant(string("mu")));
  K_zero   = Kabbala(string("zero"),0.);
  num_2    = Kabbala(string("2"),2.);    	
  num_3    = Kabbala(string("3"),3.);    	
  num_4    = Kabbala(string("4"),4.);    		
  num_6    = Kabbala(string("6"),6.);    		
  root2    = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  invroot2 = Kabbala(string("1/\\sqrt{2}"),sqrt(.5));
}

void Interaction_Model_sLepton_sQuark::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1,num_1,num_5;
  num_1    = Kabbala(string("1"),1.);    	
  num_5    = Kabbala(string("5"),5.);    		
  
  //snu - squark - squark - snu
  for (short int i=81;i<84;++i) {
    Flavour snu = Flavour(kf::code(i));
    if (snu.IsOn()) {
      //uptypes
      for (short int j=51;j<57;++j) {
	Flavour squark1 = Flavour(kf::code(j));
	if (squark1.IsOn()) {
	  for (short int k=51;k<57;++k) {
	    Flavour squark2 = Flavour(kf::code(k));
	    if (squark2.IsOn()) {
      
	      vertex[vanz].in[0] = snu;
	      vertex[vanz].in[1] = snu;
	      vertex[vanz].in[2] = squark1;
	      vertex[vanz].in[3] = squark2.Bar();
	  
	      vertex[vanz].nleg  = 4;  
	      
	      Kabbala help = K_zero;
	      if (j==k) help = num_1;
	  
	      kcpl0 = -M_I*g1*g1/(num_3*costW*costW)*
		(help + (num_3-num_2*num_4*sintW*sintW)/(num_4*sintW*sintW)*
		 (K_Z_U(0,j-51)*K_Z_U(0,k-51) +   
		  K_Z_U(1,j-51)*K_Z_U(1,k-51) +   
		  K_Z_U(2,j-51)*K_Z_U(2,k-51)));
	      
	      kcpl1 = kcpl0;
	      
	      vertex[vanz].cpl[0]  = kcpl0.Value(); 
	      vertex[vanz].cpl[1]  = kcpl1.Value();
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	      
	      vertex[vanz].ncf   = 1;
	      vertex[vanz].Color = new Color_Function(cf::D);     
	      vertex[vanz].Color->SetParticleArg(2,3);     
	      vertex[vanz].Color->SetStringArg('2','3');     
	      
	      vertex[vanz].nlf     = 1;
	      vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	      
	      vertex[vanz].on      = 1;
	      if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	
	    }
	  }
	}
      }//downtypes
      for (short int j=61;j<67;++j) {
	Flavour squark1 = Flavour(kf::code(j));
	if (squark1.IsOn()) {
	  for (short int k=61;k<67;++k) {
	    Flavour squark2 = Flavour(kf::code(k));
	    if (squark2.IsOn()) {
      
	      vertex[vanz].in[0] = snu;
	      vertex[vanz].in[1] = snu;
	      vertex[vanz].in[2] = squark1;
	      vertex[vanz].in[3] = squark2.Bar();
	  
	      vertex[vanz].nleg  = 4;  
	      
	      Kabbala help = K_zero;
	      if (j==k) help = num_1;
	  
	      kcpl0 = M_I*g1*g1/(num_6*costW*costW)*
		(help + (num_3-num_4*sintW*sintW)/(num_2*sintW*sintW)*
		 (K_Z_D(0,j-61)*K_Z_D(0,k-61) +   
		  K_Z_D(1,j-61)*K_Z_D(1,k-61) +   
		  K_Z_D(2,j-61)*K_Z_D(2,k-61)));
	      
	      kcpl1 = kcpl0;
	      
	      vertex[vanz].cpl[0]  = kcpl0.Value(); 
	      vertex[vanz].cpl[1]  = kcpl1.Value();
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	      
	      vertex[vanz].ncf   = 1;
	      vertex[vanz].Color = new Color_Function(cf::D);     
	      vertex[vanz].Color->SetParticleArg(2,3);     
	      vertex[vanz].Color->SetStringArg('2','3');     
	      
	      vertex[vanz].nlf     = 1;
	      vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	      
	      vertex[vanz].on      = 1;
	      if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	    }
	  }
	}
      }
      //snu - sLep - sDown - sUp
      for (int k=71;k<77;k++) {
	Flavour slep = Flavour(kf::code(k));
	if (slep.IsOn()) {
	  for (int l=61;l<67;l++) {
	    Flavour sdown = Flavour(kf::code(l));
	    if (sdown.IsOn()) {
	      for (int m=51;m<57;m++) {
		Flavour sup = Flavour(kf::code(m));
		if (sup.IsOn()) {
		  
		  vertex[vanz].in[0] = snu;
		  vertex[vanz].in[1] = slep;
		  vertex[vanz].in[2] = sup;
		  vertex[vanz].in[3] = sdown.Bar();
		  
		  vertex[vanz].nleg  = 4;  
		  
		  Kabbala factor = K_zero;
	      
		  for (int J=0;J<3;J++) {
		    for (int K=0;K<3;K++) {
		      for (int L=0;L<3;L++) {
			factor += K_Z_Nu(J,i-81)*K_Z_U(L,m-51)*K_CKM(L,K)*
			  (g2*g2/num_2*K_Z_D(K,l-61)*K_Z_L(J,k-71) + 
			   K_l(J)*K_d(K)*K_Z_D(K+3,l-61)*K_Z_L(J+3,k-71));  
		      }
		    }
		  }
		  
		  kcpl0 = -M_I*factor;
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0.Value(); 
		  vertex[vanz].cpl[1]  = kcpl1.Value();
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		  
		  vertex[vanz].ncf   = 1;
		  vertex[vanz].Color = new Color_Function(cf::D);     
		  vertex[vanz].Color->SetParticleArg(2,3);     
		  vertex[vanz].Color->SetStringArg('2','3');     
		  
		  vertex[vanz].nlf     = 1;
		  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		  
		  vertex[vanz].on      = 1;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  //squark - sLep - sLep - squark 
  for (int i=71;i<77;i++) {
    Flavour slepi = Flavour(kf::code(i));
    if (slepi.IsOn()) {
      for (int j=71;j<77;j++) {
	Flavour slepj = Flavour(kf::code(j));
	if (slepj.IsOn()) {
	  //uptypes
	  for (int k=51;k<57;k++) {
	    Flavour supk = Flavour(kf::code(k));
	    if (supk.IsOn()) {
	      for (int l=51;l<57;l++) {
		Flavour supl = Flavour(kf::code(l));
		if (supl.IsOn()) {
	
		  vertex[vanz].in[0] = supk.Bar();
		  vertex[vanz].in[1] = supl.Bar();
		  vertex[vanz].in[2] = slepj;
		  vertex[vanz].in[3] = slepi.Bar();
		  
		  vertex[vanz].nleg  = 4;  
		  
		  Kabbala d_ij = K_zero;
		  if (i==j) d_ij = num_1;
		  
		  Kabbala d_kl = K_zero;
		  if (k==l) d_kl = num_1;
		  
		  Kabbala ZLsum_ij = K_zero;
		  Kabbala ZUsum_kl = K_zero;
		  
		  for (int t=0;t<3;t++) {
		    ZLsum_ij += K_Z_L(t,i-71)*K_Z_L(t,j-71);
		    ZUsum_kl += K_Z_U(t,k-51)*K_Z_U(t,l-51);
		  }
		
		  kcpl0 = M_I*g1*g1/(num_6*costW*costW)*((num_3+num_2*num_6*sintW*sintW)/(num_2*sintW*sintW)*
							 ZLsum_ij*ZUsum_kl -
							 num_6*d_kl*ZLsum_ij - 
							 num_5*d_ij*ZUsum_kl + 
							 num_4*d_ij*d_kl);
		  
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0.Value(); 
		  vertex[vanz].cpl[1]  = kcpl1.Value();
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		  
		  vertex[vanz].ncf   = 1;
		  vertex[vanz].Color = new Color_Function(cf::D);     
		  vertex[vanz].Color->SetParticleArg(0,1);     
		  vertex[vanz].Color->SetStringArg('0','1');     
		  
		  vertex[vanz].nlf     = 1;
		  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		  
		  vertex[vanz].on      = 1;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
		}
	      }
	    }
	  }
	  //downtypes
	  for (int k=61;k<67;k++) {
	    Flavour sdownk = Flavour(kf::code(k));
	    if (sdownk.IsOn()) {
	      for (int l=61;l<67;l++) {
		Flavour sdownl = Flavour(kf::code(l));
		if (sdownl.IsOn()) {
	
		  vertex[vanz].in[0] = sdownk;
		  vertex[vanz].in[1] = sdownl;
		  vertex[vanz].in[2] = slepj;
		  vertex[vanz].in[3] = slepi.Bar();
		  
		  vertex[vanz].nleg  = 4;  
		  
		  Kabbala d_ij = K_zero;
		  if (i==j) d_ij = num_1;
		  
		  Kabbala d_kl = K_zero;
		  if (k==l) d_kl = num_1;
		  
		  Kabbala ZLsum_ij = K_zero;
		  Kabbala ZDsum_kl = K_zero;
		  
		  for (int t=0;t<3;t++) {
		    ZLsum_ij += K_Z_L(t,i-71)*K_Z_L(t,j-71);
		    ZDsum_kl += K_Z_D(t,k-61)*K_Z_D(t,l-61);
		  }
		
		  Kabbala addendum = K_zero;

		  for (int I=0;I<3;I++) {
		    for (int J=0;J<3;J++) {
		      addendum += K_l(I)*K_d(J)*
			(K_Z_L(I+3,i-71)*K_Z_L(I,j-71)*K_Z_D(J,k-61)*K_Z_D(J+3,l-61) + 
			 K_Z_L(I,i-71)*K_Z_L(I+3,j-71)*K_Z_D(J+3,k-61)*K_Z_D(J,l-61));
		    }
		  }

		  kcpl0 = M_I*(g1*g1/(num_6*costW*costW)*(-num_3/(num_2*sintW*sintW)*
							  ZLsum_ij*ZDsum_kl +
							  num_3*d_kl*ZLsum_ij + 
							  d_ij*ZDsum_kl -
							  num_2*d_ij*d_kl) - addendum);
		  
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0.Value(); 
		  vertex[vanz].cpl[1]  = kcpl1.Value();
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		  
		  vertex[vanz].ncf   = 1;
		  vertex[vanz].Color = new Color_Function(cf::D);     
		  vertex[vanz].Color->SetParticleArg(0,1);     
		  vertex[vanz].Color->SetStringArg('0','1');     
		  
		  vertex[vanz].nlf     = 1;
		  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		  
		  vertex[vanz].on      = 1;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


Kabbala Interaction_Model_sLepton_sQuark::K_l(short int i) 
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("l^")+string(hi),
		 -Flavour(kf::code(2*i+11)).Yuk()/v1.Value()*sqrt(2.));

}

Kabbala Interaction_Model_sLepton_sQuark::K_u(short int i)
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("u^")+string(hi),
		 Flavour(kf::code(2*i+2)).Yuk()/v2.Value()*sqrt(2.));
}

Kabbala Interaction_Model_sLepton_sQuark::K_d(short int i)
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("d^")+string(hi),
		 -Flavour(kf::code(2*i+1)).Yuk()/v1.Value()*sqrt(2.));
}

Kabbala Interaction_Model_sLepton_sQuark::K_Z_Nu(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_\\nu"),
		 ComplexMatrixElement(string("Z_nu"),i,j));
}  
 
Kabbala Interaction_Model_sLepton_sQuark::K_Z_L(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 ComplexMatrixElement(string("Z_l"),i,j));
}  

Kabbala Interaction_Model_sLepton_sQuark::K_Z_D(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_D"),
		 ComplexMatrixElement(std::string("Z_d"),i,j));
}  

Kabbala Interaction_Model_sLepton_sQuark::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 ComplexMatrixElement(std::string("Z_u"),i,j));
}  

Kabbala Interaction_Model_sLepton_sQuark::K_yuk(Flavour fl) {
  return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
}

Kabbala Interaction_Model_sLepton_sQuark::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

int Interaction_Model_sLepton_sQuark::gen_sLep(Flavour fl)
{
  int gen_sL;

  if (fl.Kfcode() == 71 || fl.Kfcode() == 74)
    gen_sL = 0;
  if (fl.Kfcode() == 72 || fl.Kfcode() == 75)
    gen_sL = 1;
  if (fl.Kfcode() == 73 || fl.Kfcode() == 76)
    gen_sL = 2;

  if (fl.Kfcode() == 81)
    gen_sL = 0;
  if (fl.Kfcode() == 82)
    gen_sL = 1;
  if (fl.Kfcode() == 83)
    gen_sL = 2;
  
  return gen_sL;
}

int Interaction_Model_sLepton_sQuark::gen_sUp(Flavour fl)
{
  int gen_sUp;

  if (fl.Kfcode() == 51 || fl.Kfcode() == 54)
    gen_sUp = 0;
  if (fl.Kfcode() == 52 || fl.Kfcode() == 55)
    gen_sUp = 1;
  if (fl.Kfcode() == 53 || fl.Kfcode() == 56)
    gen_sUp = 2;

  return gen_sUp;
}

int Interaction_Model_sLepton_sQuark::gen_sDown(Flavour fl)
{
  int gen_sDown;

  if (fl.Kfcode() == 61 || fl.Kfcode() == 64)
    gen_sDown = 0;
  if (fl.Kfcode() == 62 || fl.Kfcode() == 65)
    gen_sDown = 1;
  if (fl.Kfcode() == 63 || fl.Kfcode() == 66)
    gen_sDown = 2;

  return gen_sDown;
}

Kabbala Interaction_Model_sLepton_sQuark::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_sLepton_sQuark::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
