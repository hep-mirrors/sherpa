#include "MODEL/Interaction_Models/Interaction_Model_sLepton_sQuark.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sLepton_sQuark::Interaction_Model_sLepton_sQuark(MODEL::Model_Base * _model,
							   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),rpa->gen.CplScale())));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
    vev      = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
  }else{
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ComplexConstant(std::string("csin2_thetaW"))));
    vev      = Kabbala(string("v_{EW}"),ComplexConstant(std::string("cvev")));
  }

  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  v1       = Kabbala(string("v_1"),vev.Value() *
		     sqrt(1./(1.+sqr(ScalarConstant(string("tan(beta)"))))));
  v2       = Kabbala(string("v_2"),vev.Value() *
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
  int s_qu,s_lep;

  //snu - squark - squark - snu
  for (short int i=0;i<3;i++) {
    Flavour snu = Flavour((kf_code)(1000012+2*i));
    if (snu.IsOn()) {
      //uptypes
      for (short int j=1;j<7;++j) {
	if (j<4) s_qu=1000000 + 2*j;
	else     s_qu=2000000 + 2*j - 6;
	Flavour squark1 = Flavour((kf_code)(s_qu));
	if (squark1.IsOn()) {
	  for (short int k=1;k<7;++k) {
	    if (k<4) s_qu=1000000 + 2*k;
	    else     s_qu=2000000 + 2*k - 6;
	    Flavour squark2 = Flavour((kf_code)(s_qu));
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
		 (K_Z_U(0,j-1)*K_Z_U(0,k-1) +   
		  K_Z_U(1,j-1)*K_Z_U(1,k-1) +   
		  K_Z_U(2,j-1)*K_Z_U(2,k-1)));
	      
	      kcpl1 = kcpl0;
	      
	      vertex[vanz].cpl[0]  = kcpl0; 
	      vertex[vanz].cpl[1]  = kcpl1;
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      
	      
	      vertex[vanz].Color.push_back(Color_Function(cf::D));     
	      vertex[vanz].Color.back().SetParticleArg(2,3);     
	      vertex[vanz].Color.back().SetStringArg('2','3');     
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	      
	      vertex[vanz].on      = 1;
	      vertex[vanz].oew     = 2;
	      if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	
	    }
	  }
	}
      }//downtypes
      for (short int j=1;j<7;++j) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour squark1 = Flavour((kf_code)(s_qu));
	if (squark1.IsOn()) {
	  for (short int k=1;k<7;++k) {
	    if (k<4) s_qu=1000000 + 2*k - 1;
	    else     s_qu=2000000 + 2*k - 7;
	    Flavour squark2 = Flavour((kf_code)(s_qu));
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
		 (K_Z_D(0,j-1)*K_Z_D(0,k-1) +   
		  K_Z_D(1,j-1)*K_Z_D(1,k-1) +   
		  K_Z_D(2,j-1)*K_Z_D(2,k-1)));
	      
	      kcpl1 = kcpl0;
	      
	      vertex[vanz].cpl[0]  = kcpl0; 
	      vertex[vanz].cpl[1]  = kcpl1;
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      
	      
	      vertex[vanz].Color.push_back(Color_Function(cf::D));     
	      vertex[vanz].Color.back().SetParticleArg(2,3);     
	      vertex[vanz].Color.back().SetStringArg('2','3');     
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	      
	      vertex[vanz].on      = 1;
	      vertex[vanz].oew     = 2;
	      if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	    }
	  }
	}
      }
      //snu - sLep - sDown - sUp
      for (int k=1;k<7;k++) {
	if (k<4) s_lep = 1000010 + 2*k - 1;
	else     s_lep = 2000010 + 2*k - 7;
	Flavour slep = Flavour((kf_code)(s_lep));
	if (slep.IsOn()) {
	  for (int l=1;l<7;l++) {
	    if (l<4) s_qu=1000000 + 2*l - 1;
	    else     s_qu=2000000 + 2*l - 7;
	    Flavour sdown = Flavour((kf_code)(s_qu));
	    if (sdown.IsOn()) {
	      for (int m=1;m<7;m++) {
		if (m<4) s_qu=1000000 + 2*m;
		else     s_qu=2000000 + 2*m - 6;
		Flavour sup = Flavour((kf_code)(s_qu));
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
			factor += K_Z_Nu(J,i)*K_Z_U(L,m-1)*K_CKM(L,K)*
			  (g2*g2/num_2*K_Z_D(K,l-1)*K_Z_L(J,k-1) + 
			   K_l(J)*K_d(K)*K_Z_D(K+3,l-1)*K_Z_L(J+3,k-1));  
		      }
		    }
		  }
		  
		  kcpl0 = -M_I*factor;
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0; 
		  vertex[vanz].cpl[1]  = kcpl1;
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  
		  
		  vertex[vanz].Color.push_back(Color_Function(cf::D));     
		  vertex[vanz].Color.back().SetParticleArg(2,3);     
		  vertex[vanz].Color.back().SetStringArg('2','3');     
		  
		  
		  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		  
		  vertex[vanz].on      = 1;
		  vertex[vanz].oew     = 2;
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
  for (int i=1;i<7;i++) {
    if (i<4) s_lep = 1000010 + 2*i - 1;
    else     s_lep = 2000010 + 2*i - 7;
    Flavour slepi = Flavour((kf_code)(s_lep));
    if (slepi.IsOn()) {
      for (int j=1;j<7;j++) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour slepj = Flavour((kf_code)(s_lep));
	if (slepj.IsOn()) {
	  //uptypes
	  for (int k=1;k<7;k++) {
	    if (k<4) s_qu=1000000 + 2*k;
	    else     s_qu=2000000 + 2*k - 6;
	    Flavour supk = Flavour((kf_code)(s_qu));
	    if (supk.IsOn()) {
	      for (int l=1;l<7;l++) {
		if (l<4) s_qu=1000000 + 2*l;
		else     s_qu=2000000 + 2*l - 6;
		Flavour supl = Flavour((kf_code)(s_qu));
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
		    ZLsum_ij += K_Z_L(t,i-1)*K_Z_L(t,j-1);
		    ZUsum_kl += K_Z_U(t,k-1)*K_Z_U(t,l-1);
		  }
		
		  kcpl0 = M_I*g1*g1/(num_6*costW*costW)*((num_3+num_2*num_6*sintW*sintW)/(num_2*sintW*sintW)*
							 ZLsum_ij*ZUsum_kl -
							 num_6*d_kl*ZLsum_ij - 
							 num_5*d_ij*ZUsum_kl + 
							 num_4*d_ij*d_kl);
		  
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0; 
		  vertex[vanz].cpl[1]  = kcpl1;
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  
		  
		  vertex[vanz].Color.push_back(Color_Function(cf::D));     
		  vertex[vanz].Color.back().SetParticleArg(0,1);     
		  vertex[vanz].Color.back().SetStringArg('0','1');     
		  
		  
		  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		  
		  vertex[vanz].on      = 1;
		  vertex[vanz].oew     = 2;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
		}
	      }
	    }
	  }
	  //downtypes
	  for (int k=1;k<7;k++) {
	    if (k<4) s_qu=1000000 + 2*k - 1;
	    else     s_qu=2000000 + 2*k - 7;
	    Flavour sdownk = Flavour((kf_code)(s_qu));
	    if (sdownk.IsOn()) {
	      for (int l=1;l<7;l++) {
		if (l<4) s_qu=1000000 + 2*l - 1;
		else     s_qu=2000000 + 2*l - 7;
		Flavour sdownl = Flavour((kf_code)(s_qu));
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
		    ZLsum_ij += K_Z_L(t,i-1)*K_Z_L(t,j-1);
		    ZDsum_kl += K_Z_D(t,k-1)*K_Z_D(t,l-1);
		  }
		
		  Kabbala addendum = K_zero;

		  for (int I=0;I<3;I++) {
		    for (int J=0;J<3;J++) {
		      addendum += K_l(I)*K_d(J)*
			(K_Z_L(I+3,i-1)*K_Z_L(I,j-1)*K_Z_D(J,k-1)*K_Z_D(J+3,l-1) + 
			 K_Z_L(I,i-1)*K_Z_L(I+3,j-1)*K_Z_D(J+3,k-1)*K_Z_D(J,l-1));
		    }
		  }

		  kcpl0 = M_I*(g1*g1/(num_6*costW*costW)*(-num_3/(num_2*sintW*sintW)*
							  ZLsum_ij*ZDsum_kl +
							  num_3*d_kl*ZLsum_ij + 
							  d_ij*ZDsum_kl -
							  num_2*d_ij*d_kl) - addendum);
		  
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0; 
		  vertex[vanz].cpl[1]  = kcpl1;
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  
		  
		  vertex[vanz].Color.push_back(Color_Function(cf::D));     
		  vertex[vanz].Color.back().SetParticleArg(0,1);     
		  vertex[vanz].Color.back().SetStringArg('0','1');     
		  
		  
		  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		  
		  vertex[vanz].on      = 1;
		  vertex[vanz].oew     = 2;
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
		 -K_yuk(Flavour((kf_code)(2*i+11))).Value()/v1.Value()*sqrt(2.));

}

Kabbala Interaction_Model_sLepton_sQuark::K_u(short int i)
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("u^")+string(hi),
		 K_yuk(Flavour((kf_code)(2*i+2))).Value()/v2.Value()*sqrt(2.));
}

Kabbala Interaction_Model_sLepton_sQuark::K_d(short int i)
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("d^")+string(hi),
		 -K_yuk(Flavour((kf_code)(2*i+1))).Value()/v1.Value()*sqrt(2.));
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
  if(ScalarNumber(std::string("WidthScheme"))==0){
    return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
  }else{
    return Kabbala(string("M_{"+fl.TexName()+"}"),sqrt(sqr(fl.Yuk())-
                   Complex(0.,1.)*fl.Width()*fl.Yuk()));
  }
}

Kabbala Interaction_Model_sLepton_sQuark::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

int Interaction_Model_sLepton_sQuark::gen_sLep(Flavour fl)
{
  int gen_sL;

  if (fl.Kfcode() == 1000011 || fl.Kfcode() == 2000011)
    gen_sL = 0;
  if (fl.Kfcode() == 1000013 || fl.Kfcode() == 2000013)
    gen_sL = 1;
  if (fl.Kfcode() == 1000015 || fl.Kfcode() == 2000015)
    gen_sL = 2;

  if (fl.Kfcode() == 1000012)
    gen_sL = 0;
  if (fl.Kfcode() == 1000014)
    gen_sL = 1;
  if (fl.Kfcode() == 1000016)
    gen_sL = 2;
  
  return gen_sL;
}

int Interaction_Model_sLepton_sQuark::gen_sUp(Flavour fl)
{
  int gen_sUp;
  
  if (fl.Kfcode() == 1000002 || fl.Kfcode() == 2000002)
    gen_sUp = 0;
  if (fl.Kfcode() == 1000004 || fl.Kfcode() == 2000004)
    gen_sUp = 1;
  if (fl.Kfcode() == 1000006 || fl.Kfcode() == 2000006)
    gen_sUp = 2;
  
  return gen_sUp;
}

int Interaction_Model_sLepton_sQuark::gen_sDown(Flavour fl)
{
  int gen_sDown;

  if (fl.Kfcode() == 1000001 || fl.Kfcode() == 2000001)
    gen_sDown = 0;
  if (fl.Kfcode() == 1000003 || fl.Kfcode() == 2000003)
    gen_sDown = 1;
  if (fl.Kfcode() == 1000005 || fl.Kfcode() == 2000005)
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
