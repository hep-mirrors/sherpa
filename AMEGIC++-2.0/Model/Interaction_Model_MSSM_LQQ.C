#include "Interaction_Model_MSSM_LQQ.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_MSSM_LQQ::Interaction_Model_MSSM_LQQ(MODEL::Model_Base * _model,
							 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),Ecms2)));
  sintW    = Kabbala(string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(string("sin2_thetaW"))));
  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  K_zero   = Kabbala(string("zero"),0.);
  num_half = Kabbala(string("1/2"),0.5);    	
  num_1    = Kabbala(string("1"),1.);    	
  num_2    = Kabbala(string("2"),2.);    	
  num_3    = Kabbala(string("3"),3.);    	
  num_4    = Kabbala(string("4"),4.);    		
  num_6    = Kabbala(string("6"),6.);
  root2    = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  invroot2 = Kabbala(string("1/\\sqrt{2}"),sqrt(.5));

  p_mssm   = new Interaction_Model_MSSM(p_model,_cplscheme,_yukscheme); 
}

Interaction_Model_MSSM_LQQ::~Interaction_Model_MSSM_LQQ()
{
  delete p_momssm;
}

void Interaction_Model_MSSM_LQQ::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_momssm->c_FFV(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_momssm->c_VVV(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_momssm->c_VVVV(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_momssm->c_FFS(vertex,vanz);
  Kabbala kcpl0,kcpl1; 
  //d-quark - lepton - sup
  //u-quark - lepton - sdown

  for (short int u=2;u<7;u+=2) {
    Flavour flavu = Flavour(kf::code(u));    
    for (short int d=61;d<67;d++) {
      Flavour flavd = Flavour(kf::code(d));
    for (short int l=11;l<16;l+=2) {
      Flavour flavl = Flavour(kf::code(l));
	if (flavu.IsOn() && flavd.IsOn() && flavl.IsOn()) {
	  vertex[vanz].in[0] = flavu;
	  vertex[vanz].in[1] = flavl;
	  vertex[vanz].in[2] = flavd.Bar();
	  
	  kcpl0 = K_zero;
	  kcpl1 = -M_I*num_half*conj_K_LQQ(l,u,d)*conj_K_CKM((u-2)/2,gen_sDown(flav3));
			
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::D);     
	  vertex[vanz].Color->SetParticleArg(0,1);     
	  vertex[vanz].Color->SetStringArg('0','1');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    } 
  }
}

void Interaction_Model_MSSM_LQQ::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_momssm->c_VVS(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_momssm->c_SSS(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_momssm->c_SSV(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_momssm->c_SSVV(vertex,vanz);
}

void Interaction_Model_MSSM_LQQ::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_momssm->c_SSSS(vertex,vanz);
}



Kabbala Interaction_Model_MSSM_LQQ::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_MSSM_LQQ::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 
Kabbala Interaction_Model_MSSM_LQQ::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_MSSM_LQQ::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 
Kabbala Interaction_Model_MSSM_LQQ::K_LQQ(short int l,short int u,short int d)       
{   
  char hl[2],hu[2],hd[2];
  sprintf(hl,"%i",l);
  sprintf(hu,"%i",u);
  sprintf(hd,"%i",d);
  Complex entry;
  ComplexMatrixElement(string("LQQ_")+string(hl),u,d);
  return Kabbala(string("V^{LQQ}_{")+string(hl)+string(hu)+string(hd)+string("}"),entry);
} 
  
Kabbala Interaction_Model_MSSM_LQQ::conj_K_LQQ(short int l,short int u,short int d)       
{   
  char hl[2],hu[2],hd[2];
  sprintf(hl,"%i",l);
  sprintf(hu,"%i",u);
  sprintf(hd,"%i",d);
  Complex entry;
  ComplexMatrixElement(string("LQQ_")+string(hl),u,d);
  return Kabbala(string("V^{LQQ}_{")+string(hl)+string(hu)+string(hd)+string("}"),conj(entry));
} 
 

Kabbala Interaction_Model_MSSM_LQQ::K_Z_D(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_D"),
		 ComplexMatrixElement(std::string("Z_d"),i,j));
}  

Kabbala Interaction_Model_MSSM_LQQ::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 ComplexMatrixElement(std::string("Z_u"),i,j));
}  


//we use transposed convention !!! 

int Interaction_Model_MSSM_LQQ::gen_sUp(Flavour fl)
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

int Interaction_Model_MSSM_LQQ::gen_sDown(Flavour fl)
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

