#include "Interaction_Model_EW.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_EW::Interaction_Model_EW(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
  g2    = Kabbala(string("g_1/\\sin\\theta_W"), 
		  g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  vev   = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
}

void Interaction_Model_EW::c_FFV(Single_Vertex * vertex,int & vanz)
{
  Flavour flphoton(kf::photon);
  Flavour flZ(kf::Z);
  Flavour flW(kf::W);
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav1               = Flavour(kf::code(i));
    Kabbala charge1             = Kabbala(string("Q_{")+flav1.TexName()+string("}"),flav1.Charge());
    Kabbala isoweak1            = Kabbala(string("T_{")+flav1.TexName()+string("}"),flav1.IsoWeak());

    Kabbala kcpl0,kcpl1;    
    if (flav1.IsOn()) {
      for (short int j=i;j<17;j++) {
	if (j==7) j=11;
	Flavour flav2           = Flavour(kf::code(j));
	Kabbala charge2         = Kabbala(string("Q_{")+flav2.TexName()+string("}"),flav2.Charge());
	Kabbala isoweak2        = Kabbala(string("T_{")+ flav2.TexName()+string("}"),flav2.IsoWeak());	
	
	if (flav2.IsOn()) {
	  if (flav1==flav2) {
	    //photon
	    if (flphoton.IsOn()) {
	      kcpl0             = -g1*M_I*charge1;
	      kcpl1             = kcpl0;
	      if (!ATOOLS::IsZero(kcpl0.Value())) {
		vertex[vanz].in[0]   = flav1;
		vertex[vanz].in[1]   = Flavour(kf::photon);
		vertex[vanz].in[2]   = flav2;
		vertex[vanz].cpl[0]  = kcpl0.Value();
		vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
		vertex[vanz].cpl[2]  = 0.;
		vertex[vanz].cpl[3]  = 0.;
		vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
		vertex[vanz].ncf     = 1;
		if (flav1.Strong()) {
		  vertex[vanz].Color   = new Color_Function(cf::D);     
		  vertex[vanz].Color->SetParticleArg(0,2);     
		  vertex[vanz].Color->SetStringArg('0','2');     
		}
		else 
		  vertex[vanz].Color = new Color_Function(cf::None);

		vertex[vanz].nlf     = 1;
		vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
		vertex[vanz].Lorentz->SetParticleArg(1);     

		vertex[vanz].on      = 1;
		vanz++;
	      }
	    }
	    //Z
	    if (flZ.IsOn()) {
	      
	      kcpl0             = M_I/costW*charge1*sintW*sintW*g2;
	      kcpl1             = -M_I/costW*(isoweak1-charge1*sintW*sintW)*g2;
	     
	      vertex[vanz].in[0]     = flav1;
	      vertex[vanz].in[1]     = Flavour(kf::Z);
	      vertex[vanz].in[2]     = flav2;
	      vertex[vanz].cpl[0]    = kcpl0.Value();
	      vertex[vanz].cpl[1]    = kcpl1.Value();
	      vertex[vanz].cpl[2]    = 0.;
	      vertex[vanz].cpl[3]    = 0.;
	      vertex[vanz].Str       = (kcpl0*PR+kcpl1*PL).String();
	
	      vertex[vanz].ncf       = 1;
	      if (flav1.Strong()) {
		vertex[vanz].Color     = new Color_Function(cf::D);     
		vertex[vanz].Color->SetParticleArg(0,2);     
		vertex[vanz].Color->SetStringArg('0','2');     
	      }
	      else 
		  vertex[vanz].Color = new Color_Function(cf::None);

	      vertex[vanz].nlf     = 1;
	      vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
	      vertex[vanz].Lorentz->SetParticleArg(1);     

	      vertex[vanz].on     = 1;
	      vanz++;
	    }
	  }
	  //W
	  if (flW.IsOn()) {
	    short int hit = 1;
	    Kabbala kcpl0,kcpl1;
	    kcpl0 = Kabbala(string("zero"),0.);
	    kcpl1 = Kabbala(string("1"),0.);

	    if (!((flav1.IsDowntype() && flav2.IsUptype()) ||
                  (flav2.IsDowntype() && flav1.IsUptype()))) hit = 0;
	    if ((flav1.IsLepton() && !flav2.IsLepton()) ||
		(flav1.IsQuark() && !flav2.IsQuark()) ) hit = 0;
	    if (hit==1) {
	      if (flav1.IsDowntype() && i>10 && j==i+1) 
		kcpl1 = -M_I/root2*g2;
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		kcpl1 = -M_I/root2*g2*K_CKM((i-1)/2,j/2-1);
		else 	    
		kcpl1 = -M_I/root2*g2*K_CKM(i/2-1,(j-1)/2);		
	      }
	      if (!ATOOLS::IsZero(kcpl1.Value())) {
		vertex[vanz].in[1] = Flavour(kf::W);
		if (flav1.IsDowntype()) {
		  vertex[vanz].in[0] = flav1;
		  vertex[vanz].in[2] = flav2;
		}
		else {
		  vertex[vanz].in[0] = flav2;
		  vertex[vanz].in[2] = flav1;
		}
		
		vertex[vanz].cpl[0]  = kcpl0.Value();
		vertex[vanz].cpl[1]  = kcpl1.Value();
		vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
		vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
		vertex[vanz].ncf   = 1;
		if (flav1.Strong()) {
		  vertex[vanz].Color = new Color_Function(cf::D);     
		  vertex[vanz].Color->SetParticleArg(0,2);     
		  vertex[vanz].Color->SetStringArg('0','2');     
		}
		else 
		  vertex[vanz].Color = new Color_Function(cf::None);

		vertex[vanz].nlf     = 1;
		vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
		vertex[vanz].Lorentz->SetParticleArg(1);     

		vertex[vanz].on      = 1;
		vanz++;
	      }
	    }
	  }
	}
      }
    }
  }
}

void Interaction_Model_EW::c_VVV(Single_Vertex* vertex,int& vanz)
{
  Flavour flW(kf::W);
  Flavour flZ(kf::Z);
  Flavour flP(kf::photon);
  Kabbala kcpl0,kcpl1,kcpl0_1,kcpl1_1,charge;

  charge = Kabbala(string("Q_{")+flW.TexName()+string("}"),flW.Charge());

  if (flW.IsOn() && flP.IsOn()) {

  // photon WW
  vertex[vanz].in[0] = flW;
  vertex[vanz].in[1] = Flavour(kf::photon);
  vertex[vanz].in[2] = flW;

  kcpl0 = M_I*g1*charge;
  kcpl1 = kcpl0;

  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

  vertex[vanz].ncf   = 1;
  vertex[vanz].Color = new Color_Function(cf::None);     

  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::Gauge3);
  vertex[vanz].Lorentz->SetParticleArg(0,1,2);     

  vertex[vanz].on      = 1;
  vanz++;
}
  if (flZ.IsOn()) {
 
  // ZWW
  vertex[vanz].in[0] = flW;
  vertex[vanz].in[1] = Flavour(kf::Z);
  vertex[vanz].in[2] = flW;

  kcpl0 = M_I*g2*charge*costW;
  kcpl1 = kcpl0;
  
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

  vertex[vanz].ncf   = 1;
  vertex[vanz].Color = new Color_Function(cf::None);     

  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::Gauge3);
  vertex[vanz].Lorentz->SetParticleArg(0,1,2);     

  vertex[vanz].on      = 1;
  vanz++;
  }
}


void Interaction_Model_EW::c_FFS(Single_Vertex* vertex,int& vanz)
{
  Flavour flh(kf::h);
  Kabbala kcpl0,kcpl1,M_h;
  if (!flh.IsOn()) return;

  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav = Flavour(kf::code(i));
    if (flav.IsOn() && flav.IsFermion() && (flav.Yuk() > 0.)) {
      
      M_h = Kabbala(string("M_{")+flav.TexName()+string("}(m_h^2)"),
		    ScalarFunction(std::string("m")+std::string(flav.Name()),sqr(flh.PSMass())));

      kcpl0 = -M_I*M_h/vev;
      kcpl1 = kcpl0;
      
      if (!ATOOLS::IsZero(kcpl0.Value())) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flh;
	vertex[vanz].in[2] = flav;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	vertex[vanz].ncf   = 1;
	if (flav.Strong()) {
	  vertex[vanz].Color = new Color_Function(cf::D);     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	}
	else 
	  vertex[vanz].Color = new Color_Function(cf::None);
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);

	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }
}


void Interaction_Model_EW::c_VVS(Single_Vertex* vertex,int& vanz)
{
  Flavour flh(kf::h);
  Kabbala kcpl0,kcpl1;  
  Kabbala num_2 = Kabbala(string("2"),2.);  
 
 if (!flh.IsOn()) return;
  
  Flavour flav(kf::W);
  // W h W
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*g2*flav.Yuk();
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     

    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vanz++;
  }

  flav = Flavour(kf::Z);
  // Z h Z
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*g2*flav.Yuk()/costW;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     

    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);  
    vertex[vanz].Lorentz->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vanz++;
  }
}


void Interaction_Model_EW::c_SSS(Single_Vertex* vertex,int& vanz)
{
  Flavour flh = Flavour(kf::h);
  Kabbala kcpl0,kcpl1,yuk;  
  Kabbala num_3 = Kabbala(string("3"),3.);  

 if (flh.IsOn()) {  
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flh;

    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*yuk*yuk*(num_3/vev);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     

    vertex[vanz].on      = 1;
    vanz++;
  }
}

void Interaction_Model_EW::c_SSSS(Single_Vertex* vertex,int& vanz)
{
  Flavour flh = Flavour(kf::h);
  Kabbala kcpl0,kcpl1,yuk;  
  Kabbala num_3 = Kabbala(string("3"),3.);  

  if (flh.IsOn()) {  
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flh;
    vertex[vanz].in[3] = flh;

    vertex[vanz].nleg  = 4;  
    
    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*yuk*yuk*(num_3/(vev*vev));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     

    vertex[vanz].on      = 1;
    vanz++;
  }
}

void Interaction_Model_EW::c_VVVV(Single_Vertex* vertex,int& vanz)
{

  Flavour flavW(kf::W);
  Flavour flavZ(kf::Z);
  Flavour flavP(kf::photon);
  Kabbala kcpl0,kcpl1;
  
  // Ph - W - W - Ph  
  if (flavW.IsOn() && flavP.IsOn()) {
  vertex[vanz].in[0] = flavP;
  vertex[vanz].in[1] = flavW;
  vertex[vanz].in[2] = flavW.Bar();
  vertex[vanz].in[3] = flavP;
  
  vertex[vanz].nleg     = 4;

  kcpl0 = -M_I*g1*g1;
  kcpl1 = kcpl0;
    
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
  vertex[vanz].ncf   = 1;
  vertex[vanz].Color = new Color_Function(cf::None);     
    
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::Gauge4);
  vertex[vanz].Lorentz->SetParticleArg(0,3,1,2);     

  vertex[vanz].on      = 1;
  vanz++;
  }

  // Ph - W - W - Z  
  if (flavW.IsOn() && flavP.IsOn() && flavZ.IsOn()) {
  vertex[vanz].in[0] = flavP;
  vertex[vanz].in[1] = flavW;
  vertex[vanz].in[2] = flavW.Bar();
  vertex[vanz].in[3] = flavZ;

  vertex[vanz].nleg     = 4;  

  kcpl0 = -M_I*g1*g1*costW/sintW;
  kcpl1 = kcpl0;
  
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
  vertex[vanz].ncf   = 1;
  vertex[vanz].Color = new Color_Function(cf::None);     
    
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::Gauge4);     
  vertex[vanz].Lorentz->SetParticleArg(0,3,1,2);     

  vertex[vanz].on      = 1;
  vanz++;
  }
  // Z - W - W - Z  
  if (flavW.IsOn() && flavZ.IsOn()) {
  vertex[vanz].in[0] = flavZ;
  vertex[vanz].in[1] = flavW;
  vertex[vanz].in[2] = flavW.Bar();
  vertex[vanz].in[3] = flavZ;
  
  vertex[vanz].nleg     = 4;

  kcpl0 = -M_I*g1*g1*costW*costW/(sintW*sintW);
  kcpl1 = kcpl0;
  
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
  vertex[vanz].ncf   = 1;
  vertex[vanz].Color = new Color_Function(cf::None);     
    
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::Gauge4);     
  vertex[vanz].Lorentz->SetParticleArg(0,3,1,2);     

  vertex[vanz].on      = 1;
  vanz++;
  }
  
 // W - W - W - W  
  if (flavW.IsOn()) {
  vertex[vanz].in[0] = flavW.Bar();
  vertex[vanz].in[1] = flavW;
  vertex[vanz].in[2] = flavW.Bar();
  vertex[vanz].in[3] = flavW.Bar();
  
  vertex[vanz].nleg     = 4;

  kcpl0 = M_I*g1*g1/(sintW*sintW);
  kcpl1 = kcpl0;
  
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
  vertex[vanz].ncf   = 1;
  vertex[vanz].Color = new Color_Function(cf::None);     
    
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::Gauge4); 
  vertex[vanz].Lorentz->SetParticleArg(0,1,2,3);     

  vertex[vanz].on      = 1;
  vanz++;
  }
}

void Interaction_Model_EW::c_SSVV(Single_Vertex* vertex,int& vanz)
{
  Kabbala num_2 = Kabbala(string("2"),2.);  

  Flavour flavW(kf::W);
  Flavour flavZ(kf::Z);
  Flavour flavh(kf::h);
  Kabbala kcpl0,kcpl1;
  
  // h - Z - Z - h  
  if (flavZ.IsOn() && flavh.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flavh;
    vertex[vanz].in[2] = flavh;
    vertex[vanz].in[3] = flavZ;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = M_I*g2*g2/(costW*costW*num_2);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
    vertex[vanz].Lorentz->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vanz++;
  }

  // h - W - W - h  
  if (flavW.IsOn() && flavh.IsOn()) {
    vertex[vanz].in[0] = flavW;
    vertex[vanz].in[1] = flavh;
    vertex[vanz].in[2] = flavh;
    vertex[vanz].in[3] = flavW;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = M_I*g2*g2/num_2;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);     
    vertex[vanz].Lorentz->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vanz++;
  }
}


Kabbala Interaction_Model_EW::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_EW::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 


