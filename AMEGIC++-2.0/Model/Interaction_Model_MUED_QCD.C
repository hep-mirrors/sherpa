#include "Interaction_Model_MUED_QCD.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_MUED_QCD::Interaction_Model_MUED_QCD(MODEL::Model_Base * _model,
						       std::string _cplscheme,
						       std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g3        = Kabbala(string("g_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),Ecms2)));
  PL        = Kabbala(string("P_L"),1.);
  PR        = Kabbala(string("P_R"),1.);
  M_I       = Kabbala(string("i"),Complex(0.,1.)); 
  sqrt2     = Kabbala(string("\\sqrt{2}"),Complex(sqrt(2.),0.)); 
  onehalf   = Kabbala(string("\\frac{1}{2}"),Complex(1./2.,0.)); 
  threehalf = Kabbala(string("\\frac{3}{2}"),Complex(3./2.,0.)); 

  m_generations.clear();
  for (int i=1;i<=ScalarNumber(std::string("Number_Of_KK_Excitations"));i++) {
    m_generations.insert(i);
  }
}

void Interaction_Model_MUED_QCD::c_FFV(std::vector<Single_Vertex>& vertex,int & vanz)
{
  p_moqcd->c_FFV(vertex,vanz);

  Kabbala kcpl0 = -g3*M_I;
  Kabbala kcpl1 = kcpl0;
  Kabbala mixing;

  std::set<int>::iterator git1;
  Flavour flav1,flav2;

  if (Flavour(kf::gluon).IsOn()) {
    for (git1=m_generations.begin();git1!=m_generations.end();git1++) {
      for (short int i=1;i<=6;i++) {
	mixing = Cos((*git1),i)*Cos((*git1),i)+Sin((*git1),i)*Sin((*git1),i);

	// N0N+n0n
	flav1  = Flavour(kf::code((50+(*git1))*100000+i));
	if (flav1.Strong() && flav1.IsOn()) { 
	  vertex[vanz].in[0]         = flav1;
	  vertex[vanz].in[1]         = Flavour(kf::gluon);
	  vertex[vanz].in[2]         = flav1;

	  vertex[vanz].cpl[0]        = (kcpl0*mixing).Value();
	  vertex[vanz].cpl[1]        = (kcpl1*mixing).Value();
	  vertex[vanz].cpl[2]        = 0.;
	  vertex[vanz].cpl[3]        = 0.;
	  vertex[vanz].Str           = (kcpl0*mixing*PR+kcpl1*mixing*PL).String();      

	  vertex[vanz].ncf           = 1;
	  vertex[vanz].Color         = new Color_Function(cf::T,1,2,0,'1','2','0');     

	  vertex[vanz].nlf           = 1;
	  vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gamma);     
	  vertex[vanz].Lorentz->SetParticleArg(1);     
                  
	  vertex[vanz].on            = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
	flav1  = Flavour(kf::code((60+(*git1))*100000+i));
	if (flav1.Strong() && flav1.IsOn()) { 
	  vertex[vanz].in[0]         = flav1;
	  vertex[vanz].in[1]         = Flavour(kf::gluon);
	  vertex[vanz].in[2]         = flav1;

	  vertex[vanz].cpl[0]        = (kcpl0*mixing).Value();
	  vertex[vanz].cpl[1]        = (kcpl1*mixing).Value();
	  vertex[vanz].cpl[2]        = 0.;
	  vertex[vanz].cpl[3]        = 0.;
	  vertex[vanz].Str           = (kcpl0*mixing*PR+kcpl1*mixing*PL).String();      

	  vertex[vanz].ncf           = 1;
	  vertex[vanz].Color         = new Color_Function(cf::T,1,2,0,'1','2','0');     

	  vertex[vanz].nlf           = 1;
	  vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gamma);     
	  vertex[vanz].Lorentz->SetParticleArg(1);     
                  
	  vertex[vanz].on            = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}

	// N0n-N0n
	mixing = Cos((*git1),i)*Sin((*git1),i)+Sin((*git1),i)*Cos((*git1),i);
	flav1  = Flavour(kf::code((50+(*git1))*100000+i));
	flav2  = Flavour(kf::code((60+(*git1))*100000+i));
	if (flav1.Strong() && flav1.IsOn() &&
	    flav2.Strong() && flav2.IsOn()) { 
	  vertex[vanz].in[0]         = flav1;
	  vertex[vanz].in[1]         = Flavour(kf::gluon);
	  vertex[vanz].in[2]         = flav2;

	  vertex[vanz].cpl[0]        = (kcpl0*mixing).Value();
	  vertex[vanz].cpl[1]        = (kcpl1*mixing).Value();
	  vertex[vanz].cpl[2]        = 0.;
	  vertex[vanz].cpl[3]        = 0.;
	  vertex[vanz].Str           = (kcpl0*mixing*PR+kcpl1*mixing*PL).String();      

	  vertex[vanz].ncf           = 1;
	  vertex[vanz].Color         = new Color_Function(cf::T,1,2,0,'1','2','0');     

	  vertex[vanz].nlf           = 1;
	  vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gamma);     
	  vertex[vanz].Lorentz->SetParticleArg(1);     
                  
	  vertex[vanz].on            = 1;
	  vertex.push_back(Single_Vertex());vanz++;

	  vertex[vanz].in[0]         = flav2;
	  vertex[vanz].in[1]         = Flavour(kf::gluon);
	  vertex[vanz].in[2]         = flav1;

	  vertex[vanz].cpl[0]        = -(kcpl0*mixing).Value();
	  vertex[vanz].cpl[1]        = -(kcpl1*mixing).Value();
	  vertex[vanz].cpl[2]        = 0.;
	  vertex[vanz].cpl[3]        = 0.;
	  vertex[vanz].Str           = (-kcpl0*mixing*PR-kcpl1*mixing*PL).String();      

	  vertex[vanz].ncf           = 1;
	  vertex[vanz].Color         = new Color_Function(cf::T,1,2,0,'1','2','0');     

	  vertex[vanz].nlf           = 1;
	  vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gamma);     
	  vertex[vanz].Lorentz->SetParticleArg(1);     
                  
	  vertex[vanz].on            = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
      }
    }
  }

    
}

void Interaction_Model_MUED_QCD::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);

  Kabbala kcpl0 = -g3;
  Kabbala kcpl1 = kcpl0; 
  
  std::set<int>::iterator git1,git2;
  if (Flavour(kf::gluon).IsOn()) {
    for (git1=m_generations.begin();git1!=m_generations.end();git1++) {
      Flavour flav = Flavour(kf::code((50+(*git1))*100000+21));
      if (flav.Strong() && flav.IsOn()) { 
	vertex[vanz].in[0]         = flav;
	vertex[vanz].in[1]         = Flavour(kf::gluon);
	vertex[vanz].in[2]         = flav;

	vertex[vanz].cpl[0]        = kcpl0.Value();
	vertex[vanz].cpl[1]        = kcpl1.Value();
	vertex[vanz].cpl[2]        = 0.;
	vertex[vanz].cpl[3]        = 0.;
	vertex[vanz].Str           = (kcpl0*PR+kcpl1*PL).String();

	vertex[vanz].ncf           = 1;
	vertex[vanz].Color         = new Color_Function(cf::F);     
	vertex[vanz].Color->SetParticleArg(0,2,1);     
	vertex[vanz].Color->SetStringArg('0','2','1');     

	vertex[vanz].nlf           = 1;
	vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gauge3);     
	vertex[vanz].Lorentz->SetParticleArg(0,1,2);     

	vertex[vanz].on            = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }

  for (git1=m_generations.begin();git1!=m_generations.end();git1++) {
    for (git2=git1;git2!=m_generations.end();git2++) {
      int gen3 = (*git1)+(*git2);
      if (gen3>ScalarNumber(std::string("Number_Of_KK_Excitations"))) continue;
      Flavour flav1 = Flavour(kf::code((50+(*git1))*100000+21));
      Flavour flav2 = Flavour(kf::code((50+(*git2))*100000+21));
      Flavour flav3 = Flavour(kf::code((50+gen3)*100000+21));
      if (flav1.Strong() && flav1.IsOn() && 
	  flav2.Strong() && flav2.IsOn() &&
	  flav3.Strong() && flav3.IsOn()) { 
	vertex[vanz].in[0]         = flav1;
	vertex[vanz].in[1]         = flav2;
	vertex[vanz].in[2]         = flav3;

	vertex[vanz].cpl[0]        = (kcpl0/sqrt2).Value();
	vertex[vanz].cpl[1]        = (kcpl1/sqrt2).Value();
	vertex[vanz].cpl[2]        = 0.;
	vertex[vanz].cpl[3]        = 0.;
	vertex[vanz].Str           = (kcpl0/sqrt2*PR+kcpl1/sqrt2*PL).String();

	vertex[vanz].ncf           = 1;
	vertex[vanz].Color         = new Color_Function(cf::F);     
	vertex[vanz].Color->SetParticleArg(0,2,1);     
	vertex[vanz].Color->SetStringArg('0','2','1');     

	vertex[vanz].nlf           = 1;
	vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gauge3);     
	vertex[vanz].Lorentz->SetParticleArg(0,1,2);     

	vertex[vanz].on            = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
}

void Interaction_Model_MUED_QCD::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);

  Kabbala kcpl0 = -M_I*g3*g3; 
  Kabbala kcpl1 = kcpl0; 
  
  // 2 gluons: 00nn
  std::set<int>::iterator git1,git2,git3;
  if (Flavour(kf::gluon).IsOn()) { 
    for (git1=m_generations.begin();git1!=m_generations.end();git1++) {
      Flavour flav = Flavour(kf::code((50+(*git1))*100000+21));
      if (flav.Strong() && flav.IsOn()) { 
	vertex[vanz].in[0]      = Flavour(kf::gluon);
	vertex[vanz].in[1]      = Flavour(kf::gluon);
	vertex[vanz].in[2]      = flav;
	vertex[vanz].in[3]      = flav;
    
	vertex[vanz].nleg       = 4;
	vertex[vanz].cpl[0]     = kcpl0.Value();
	vertex[vanz].cpl[1]     = kcpl1.Value();
	vertex[vanz].cpl[2]     = 0.;
	vertex[vanz].cpl[3]     = 0.;
	vertex[vanz].Str        = (kcpl0*PR+kcpl1*PL).String();
    
	vertex[vanz].ncf        = 3;
	vertex[vanz].nlf        = 3;
    
	vertex[vanz].Color      = new Color_Function[3];
	vertex[vanz].Lorentz    = new Lorentz_Function[3]; 
    
	vertex[vanz].Color[0]   = Color_Function(cf::F,0,2,4,'0','2','4',
						 new Color_Function(cf::F,1,3,4,'1','3','4'));
	vertex[vanz].Lorentz[0] = Lorentz_Function(lf::Gluon4);
	vertex[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     
    
	vertex[vanz].Color[1]   = Color_Function(cf::F,0,3,4,'0','3','4',
						 new Color_Function(cf::F,1,2,4,'1','2','4'));
	vertex[vanz].Lorentz[1] = Lorentz_Function(lf::Gluon4);
	vertex[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     
    
	vertex[vanz].Color[2]   = Color_Function(cf::F,0,1,4,'0','1','4',
						 new Color_Function(cf::F,3,2,4,'3','2','4')); 
	vertex[vanz].Lorentz[2] = Lorentz_Function(lf::Gluon4);     
	vertex[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
    
	vertex[vanz].on         = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    } 
    // 1 gluon: 0nm(n+m)
    for (git1=m_generations.begin();git1!=m_generations.end();git1++) {
      for (git2=git1;git2!=m_generations.end();git2++) {
	int gen3 = (*git1)+(*git2);
	if (gen3>ScalarNumber(std::string("Number_Of_KK_Excitations"))) continue;
	Flavour flav1 = Flavour(kf::code((50+(*git1))*100000+21));
	Flavour flav2 = Flavour(kf::code((50+(*git2))*100000+21));
	Flavour flav3 = Flavour(kf::code((50+gen3)*100000+21));
	if (flav1.Strong() && flav1.IsOn() && 
	    flav2.Strong() && flav2.IsOn() &&
	    flav3.Strong() && flav3.IsOn()) { 

	  vertex[vanz].in[0]      = Flavour(kf::gluon);
	  vertex[vanz].in[1]      = flav1;
	  vertex[vanz].in[2]      = flav2;
	  vertex[vanz].in[3]      = flav3;
	  
	  vertex[vanz].nleg       = 4;
	  vertex[vanz].cpl[0]     = (kcpl0/sqrt2).Value();
	  vertex[vanz].cpl[1]     = (kcpl1/sqrt2).Value();
	  vertex[vanz].cpl[2]     = 0.;
	  vertex[vanz].cpl[3]     = 0.;
	  vertex[vanz].Str        = (kcpl0/sqrt2*PR+kcpl1/sqrt2*PL).String();
	  
	  vertex[vanz].ncf        = 3;
	  vertex[vanz].nlf        = 3;
	  
	  vertex[vanz].Color      = new Color_Function[3];
	  vertex[vanz].Lorentz    = new Lorentz_Function[3]; 
	  
	  vertex[vanz].Color[0]   = Color_Function(cf::F,0,2,4,'0','2','4',
						   new Color_Function(cf::F,1,3,4,'1','3','4'));
	  vertex[vanz].Lorentz[0] = Lorentz_Function(lf::Gluon4);
	  vertex[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     
	  
	  vertex[vanz].Color[1]   = Color_Function(cf::F,0,3,4,'0','3','4',
						   new Color_Function(cf::F,1,2,4,'1','2','4'));
	  vertex[vanz].Lorentz[1] = Lorentz_Function(lf::Gluon4);
	  vertex[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     
	  
	  vertex[vanz].Color[2]   = Color_Function(cf::F,0,1,4,'0','1','4',
						   new Color_Function(cf::F,3,2,4,'3','2','4')); 
	  vertex[vanz].Lorentz[2] = Lorentz_Function(lf::Gluon4);     
	  vertex[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
	  
	  vertex[vanz].on         = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      } 
    }
  }
  for (git1=m_generations.begin();git1!=m_generations.end();git1++) {
    Flavour flav1 = Flavour(kf::code((50+(*git1))*100000+21));
    // 4KKs: nnnn
    if (flav1.Strong() && flav1.IsOn()) {
      vertex[vanz].in[0]      = flav1;
      vertex[vanz].in[1]      = flav1;
      vertex[vanz].in[2]      = flav1;
      vertex[vanz].in[3]      = flav1;
      
      vertex[vanz].nleg       = 4;
      vertex[vanz].cpl[0]     = (threehalf*kcpl0).Value();
      vertex[vanz].cpl[1]     = (threehalf*kcpl1).Value();
      vertex[vanz].cpl[2]     = 0.;
      vertex[vanz].cpl[3]     = 0.;
      vertex[vanz].Str        = (threehalf*kcpl0*PR+threehalf*kcpl1*PL).String();
      
      vertex[vanz].ncf        = 3;
      vertex[vanz].nlf        = 3;
      
      vertex[vanz].Color      = new Color_Function[3];
      vertex[vanz].Lorentz    = new Lorentz_Function[3]; 
      
      vertex[vanz].Color[0]   = Color_Function(cf::F,0,2,4,'0','2','4',
					       new Color_Function(cf::F,1,3,4,'1','3','4'));
      vertex[vanz].Lorentz[0] = Lorentz_Function(lf::Gluon4);
      vertex[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     
      
      vertex[vanz].Color[1]   = Color_Function(cf::F,0,3,4,'0','3','4',
					       new Color_Function(cf::F,1,2,4,'1','2','4'));
      vertex[vanz].Lorentz[1] = Lorentz_Function(lf::Gluon4);
      vertex[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     
      
      vertex[vanz].Color[2]   = Color_Function(cf::F,0,1,4,'0','1','4',
					       new Color_Function(cf::F,3,2,4,'3','2','4')); 
      vertex[vanz].Lorentz[2] = Lorentz_Function(lf::Gluon4);     
      vertex[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
      
      vertex[vanz].on         = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
    for (git2=git1;git2!=m_generations.end();git2++) {
      Flavour flav2 = Flavour(kf::code((50+(*git2))*100000+21));
      if (flav1.Strong() && flav1.IsOn() && 
	  flav2.Strong() && flav2.IsOn()) {
	// 4KKs: nnmm
	if (flav1!=flav2) {
	  vertex[vanz].in[0]      = flav1;
	  vertex[vanz].in[1]      = flav1;
	  vertex[vanz].in[2]      = flav2;
	  vertex[vanz].in[3]      = flav2;
	  vertex[vanz].nleg       = 4;
	  vertex[vanz].cpl[0]     = (onehalf*kcpl0).Value();
	  vertex[vanz].cpl[1]     = (onehalf*kcpl1).Value();
	  vertex[vanz].cpl[2]     = 0.;
	  vertex[vanz].cpl[3]     = 0.;
	  vertex[vanz].Str        = (threehalf*kcpl0*PR+threehalf*kcpl1*PL).String();
	    
	  vertex[vanz].ncf        = 3;
	  vertex[vanz].nlf        = 3;
	    
	  vertex[vanz].Color      = new Color_Function[3];
	  vertex[vanz].Lorentz    = new Lorentz_Function[3]; 
	    
	  vertex[vanz].Color[0]   = Color_Function(cf::F,0,2,4,'0','2','4',
						   new Color_Function(cf::F,1,3,4,'1','3','4'));
	  vertex[vanz].Lorentz[0] = Lorentz_Function(lf::Gluon4);
	  vertex[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     
	    
	  vertex[vanz].Color[1]   = Color_Function(cf::F,0,3,4,'0','3','4',
						   new Color_Function(cf::F,1,2,4,'1','2','4'));
	  vertex[vanz].Lorentz[1] = Lorentz_Function(lf::Gluon4);
	  vertex[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     
	    
	  vertex[vanz].Color[2]   = Color_Function(cf::F,0,1,4,'0','1','4',
						   new Color_Function(cf::F,3,2,4,'3','2','4')); 
	  vertex[vanz].Lorentz[2] = Lorentz_Function(lf::Gluon4);     
	  vertex[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
	    
	  vertex[vanz].on         = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
      for (git3=git2;git3!=m_generations.end();git3++) {
	Flavour flav3 = Flavour(kf::code((50+(*git3))*100000+21));
	int gen4 = (*git1)+(*git2)+(*git3);
	if (gen4>ScalarNumber(std::string("Number_Of_KK_Excitations"))) continue;
	Flavour flav4 = Flavour(kf::code((50+gen4)*100000+21));
	// 4KKs: nmk(n+m+k)
	if (flav1.Strong() && flav1.IsOn() && 
	    flav2.Strong() && flav2.IsOn() &&
	    flav3.Strong() && flav3.IsOn() &&
	    flav4.Strong() && flav4.IsOn()) {
	  vertex[vanz].in[0]      = flav1;
	  vertex[vanz].in[1]      = flav2;
	  vertex[vanz].in[2]      = flav3;
	  vertex[vanz].in[3]      = flav4;
	  vertex[vanz].nleg       = 4;
	  vertex[vanz].cpl[0]     = (onehalf*kcpl0).Value();
	  vertex[vanz].cpl[1]     = (onehalf*kcpl1).Value();
	  vertex[vanz].cpl[2]     = 0.;
	  vertex[vanz].cpl[3]     = 0.;
	  vertex[vanz].Str        = (threehalf*kcpl0*PR+threehalf*kcpl1*PL).String();
	    
	  vertex[vanz].ncf        = 3;
	  vertex[vanz].nlf        = 3;
	    
	  vertex[vanz].Color      = new Color_Function[3];
	  vertex[vanz].Lorentz    = new Lorentz_Function[3]; 
	    
	  vertex[vanz].Color[0]   = Color_Function(cf::F,0,2,4,'0','2','4',
						   new Color_Function(cf::F,1,3,4,'1','3','4'));
	  vertex[vanz].Lorentz[0] = Lorentz_Function(lf::Gluon4);
	  vertex[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     
	    
	  vertex[vanz].Color[1]   = Color_Function(cf::F,0,3,4,'0','3','4',
						   new Color_Function(cf::F,1,2,4,'1','2','4'));
	  vertex[vanz].Lorentz[1] = Lorentz_Function(lf::Gluon4);
	  vertex[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     
	    
	  vertex[vanz].Color[2]   = Color_Function(cf::F,0,1,4,'0','1','4',
						   new Color_Function(cf::F,3,2,4,'3','2','4')); 
	  vertex[vanz].Lorentz[2] = Lorentz_Function(lf::Gluon4);     
	  vertex[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
	    
	  vertex[vanz].on         = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }
}

Kabbala Interaction_Model_MUED_QCD::Cos(const int KK,const int gen)
{   
  char sup[1], kfc[1];
  sprintf(sup,"%i",KK);
  sprintf(kfc,"%i",gen);
  Flavour flav = Flavour(kf::code(gen));
  return Kabbala(string("\\cos\\gamma^{("+string(sup)+")}_{"+flav.TexName()+"}"),
		 ScalarConstant(string("cos(gamma)[")+string(sup)+"]["+string(kfc)+string("]")));
}

Kabbala Interaction_Model_MUED_QCD::Sin(const int KK,const int gen)
{   
  char sup[1], kfc[1];
  sprintf(sup,"%i",KK);
  sprintf(kfc,"%i",gen);
  Flavour flav = Flavour(kf::code(gen));
  return Kabbala(string("\\sin\\gamma^{("+string(sup)+")}_{"+flav.TexName()+"}"),
		 ScalarConstant(string("sin(gamma)[")+string(sup)+"]["+string(kfc)+string("]")));
}
