#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <stdio.h>

#include "Vertex.H"
#include "Interaction_Model_Base.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Vector.H"
#include "String_Tree.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace ATOOLS;
using namespace ATOOLS;
using namespace std;


void Vertex::GenerateVertex()
{
  int vanzsave = m_nvertex; 
  int vanz4save = m_n4vertex; 
  
  Single_Vertex dummy;
  
  for (int i=0;i<vanz4save;++i) {
    int hit = 1;
    if (hit) {
      //required by Interaction_Model_ADD due to small couplings
      if (m_v4[i].cpl[0]==Complex(0.,0.) && m_v4[i].cpl[1]==Complex(0.,0.))
	m_v4[i].on = 0;
      else { 
	if(m_v4[i].nleg==4) {
	  for (short int k=1;k<5;k++) {
	    for (short int l=1;l<5;l++) {
	      if (l!=k) {
		for (short int m=1;m<5;m++) {
		  if (m!=l && m!=k) {
		    for (short int n=1;n<5;n++) {
		      if (n!=l && n!=k && n!=m) {
			if (k!=1) k = -k;
			if (l==1) l = -l;
			if (m==1) m = -m;
			if (n==1) n = -n;
			if (SetVertex(m_v4[i],dummy,k,l,m,n)) {
			  m_v4[m_n4vertex] = dummy;
			  m_n4vertex++;
			}
			if (SetVertex(m_v4[i],dummy,-k,-l,-m,-n)) {
			  m_v4[m_n4vertex] = dummy;
			  m_n4vertex++;
			}
			k = abs(k);l=abs(l);m=abs(m);n=abs(n);
		      }
		    }
		  }
		}
	      } 
	    }
	  }
	}
      }
    }
  }
  for (int i=0;i<vanzsave;++i) {
    int hit = 1;
    if (hit) {
      //required by Interaction_Model_ADD due to small couplings
      if (m_v[i].cpl[0]==Complex(0.,0.) && m_v[i].cpl[1]==Complex(0.,0.))
	m_v[i].on = 0;
      if (m_v[i].nleg==3) {  
	for (short int k=1;k<4;k++) {
	  for (short int l=1;l<4;l++) {
	    if (l!=k) {
	      for (short int m=1;m<4;m++) {
		if (m!=l && m!=k) {
		  if (k!=1) k = -k;
		  if (l==1) l = -l;
		  if (m==1) m = -m;
		  if (SetVertex(m_v[i],dummy,k,l,m)) {
		    m_v[m_nvertex] = dummy;
		    m_nvertex++;
		  }
		  if (SetVertex(m_v[i],dummy,-k,-l,-m)) {
		    m_v[m_nvertex] = dummy;
		    m_nvertex++;
		  }
		  k = abs(k);l=abs(l);m=abs(m);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


int Vertex::CheckExistence(Single_Vertex& probe)
{
  //4 leg vertices
  if (probe.nleg==4) {
    for (short int i=0;i<m_n4vertex;++i) {     
      // 0 -> 1 2 3
      if (m_v4[i].in[0]==probe.in[0] &&
	  m_v4[i].in[1]==probe.in[1] &&
	  m_v4[i].in[2]==probe.in[2] &&
	  m_v4[i].in[3]==probe.in[3]) return 0;
    }
  }
  //3 leg vertices
  if (probe.nleg==3) {
    // 0 -> 1 2
    for (short int i=0;i<m_nvertex;++i) {
      if (m_v[i].in[0]==probe.in[0] &&
	  m_v[i].in[1]==probe.in[1] &&
	  m_v[i].in[2]==probe.in[2]) return 0;
      // 0 -> 2 1
      if (m_v[i].in[0]==probe.in[0] &&
	  m_v[i].in[1]==probe.in[2] &&
	  m_v[i].in[2]==probe.in[1]) return 0;
    }
  }
  return 1;
}
  
int Vertex::FermionRule(Single_Vertex& probe)
{
  // fermionic: particles left, anti-particles right
  if (probe.in[1].IsFermion() && !probe.in[1].IsAnti() && !probe.in[1].Majorana()) return 0;
  if (probe.in[2].IsFermion() &&  probe.in[2].IsAnti() && !probe.in[2].Majorana()) return 0;

  return 1;
}

int Vertex::SetVertex(Single_Vertex& orig, Single_Vertex& probe, int i0, int i1, int i2, int i3)
{
  probe = orig;
  
  if (i0<0) probe.in[0] = orig.in[-i0-1].Bar();
       else probe.in[0] = orig.in[i0-1];
  if (i1<0) probe.in[1] = orig.in[-i1-1].Bar();
       else probe.in[1] = orig.in[i1-1];
  if (i2<0) probe.in[2] = orig.in[-i2-1].Bar();
       else probe.in[2] = orig.in[i2-1];
  if (orig.nleg==4) {
    if (i3<0) probe.in[3] = orig.in[-i3-1].Bar();
         else if (i3<99) probe.in[3] = orig.in[i3-1];
  }
  
  if (CheckExistence(probe)==0) return 0;
  if (probe.nleg==3) {
    if (FermionRule(probe)==0) return 0;}
  
  int hc = 0;

  int cnt = 0;
  for(int i=0;i<orig.nleg;i++) {
    if(orig.in[i]!=orig.in[i].Bar()) cnt++;
  }
  if(cnt>0){
    Flavour flavlist[cnt];
    int flaglist[cnt];
    cnt = 0;
    for(int i=0;i<orig.nleg;i++){
      if(orig.in[i]!=orig.in[i].Bar()) {
	flavlist[cnt] = orig.in[i];
	if (i==0) flavlist[cnt] = flavlist[cnt].Bar();
	flaglist[cnt] = 0;
	cnt++;
      }
    }
    for (int i=0;i<cnt;i++) {
      for (int j=i+1;j<cnt;j++) {
	if (flavlist[i]==flavlist[j].Bar()) flaglist[i]=flaglist[j]=1;
      }
    }
    for (int i=0;i<cnt;i++) {
      if (flaglist[i]==0) hc = 1;
    }
  }
  
  int probehc = 0;
  if (hc) {
    // probe = h.c. ???
    for (short int i=0;i<orig.nleg;i++) { 
      Flavour flav = orig.in[i];
      if (flav!=flav.Bar()) {
	if (i==0) flav = flav.Bar();
	for (short int j=0;j<orig.nleg;j++) {
	  Flavour flav2 = probe.in[j];
	  if (j==0) flav2 = flav2.Bar();
	  if (flav2!=flav2.Bar()) {
	    if (flav==flav2.Bar()) {
	      probehc = 1;
	      break;
	    }
	  }
	}
	if (probehc) break;
      }
    }
    if (probehc) {
      int conjugate = 1;

      for (short int i=0;i<orig.nleg;i++) {
	//pseudoscalar
	if (orig.in[i]==Flavour(kf::A0)) conjugate *= -1;
      }
      
      if (orig.Lorentz->type==lf::SSV ||
	  orig.Lorentz->type==lf::Gauge3) conjugate *= -1;
      
      if (conjugate==-1) {
	for (short int i=0;i<4;i++) probe.cpl[i] = -probe.cpl[i];
      }

      Conjugate(probe.Color);

       if (probe.Lorentz->String()==string("1")) {
	//exchange left and right
	Complex help = probe.cpl[0];
	probe.cpl[0] = probe.cpl[1];
	probe.cpl[1] = help;
      }
    }
  }
  //Color and Lorentz structure changes....
  int newIndex[4];
  
  for (short int i=0;i<4;i++) newIndex[i] = -1;
  
  for (short int i=0;i<orig.nleg;i++) { 
      Flavour flav = orig.in[i];
      if (i==0) flav = flav.Bar();
      for (short int j=0;j<orig.nleg;j++) {
	Flavour flav2 = probe.in[j];
	if (j==0) flav2 = flav2.Bar();
	if (flav==flav2 || 
	    (flav==flav2.Bar() && probehc)) {
	  int hit = 1;
	  for (short int k=0;k<i;k++) {
	    if (newIndex[k]==j) {
	      hit = 0;
	      break;
	    }
	  } 
	  if (hit) {
	    newIndex[i] = j;
	    break;
	  }
	}
      }
    }
    
  ColorExchange(probe.Color,newIndex[0],newIndex[1],newIndex[2]);
  LorentzExchange(probe.Lorentz,newIndex[0],newIndex[1],newIndex[2],newIndex[3]);
  
  return 1;
}
    
void Vertex::ColorExchange(Color_Function* colfunc,int new0,int new1,int new2)
{
  // T[0,2,1]
  for (short int i=0;i<3;i++) {
    if (colfunc->type==cf::D && i==2) break;
    switch (colfunc->partarg[i]) {
      case 0: colfunc->partarg[i] = new0;
         colfunc->strarg[i]  = new0+48;
	 break;
      case 1: colfunc->partarg[i] = new1;
         colfunc->strarg[i]  = new1+48;
         break;
      case 2: colfunc->partarg[i] = new2;
         colfunc->strarg[i]  = new2+48;
         break;
    }
  }
}

void Vertex::LorentzExchange(Lorentz_Function* lorfunc,int new0,int new1,int new2,int new3)
{
  for (short int i=0;i<lorfunc->NofIndex();i++) {
    switch (lorfunc->partarg[i]) {
      case 0: lorfunc->partarg[i] = new0;break;
      case 1: lorfunc->partarg[i] = new1;break;
      case 2: lorfunc->partarg[i] = new2;break;
      case 3: lorfunc->partarg[i] = new3;break;
    }
  }
}

void Vertex::Conjugate(Color_Function* colfunc)
{
  if (colfunc->type!=cf::T) return;

  int help = colfunc->partarg[1];
  colfunc->partarg[1] = colfunc->partarg[2];
  colfunc->partarg[2] = help;

  char help2 = colfunc->strarg[1];
  colfunc->strarg[1] = colfunc->strarg[2];
  colfunc->strarg[2] = help2;  
}

Vertex::Vertex(Interaction_Model_Base * _model)
{
  /* 
     use (roughly) notation and Vertices of J. Rosiek, PRD41 (1990) 3464
     pull out common factor -i of all Vertices
  */ 
  ATOOLS::msg.Debugging()<<"In Vertex::Vertex()."<<endl;

  m_nvertex  = 10000;
  m_n4vertex = 10000;
  m_v  = new Single_Vertex[m_nvertex];
  m_v4 = new Single_Vertex[m_n4vertex];
  int vanz  = 0;
  int vanz4 = 0;

  ATOOLS::msg.Tracking()<<"   Setting vertices..."<<endl;
  _model->c_FFV(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   FFV  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_FFS(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   FFS  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVV(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   VVV  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSV(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   SSV  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVS(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   VVS  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSS(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   SSS  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVVV(m_v4,vanz4);
  ATOOLS::msg.Debugging()<<"   VVVV : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSVV(m_v4,vanz4);
  ATOOLS::msg.Debugging()<<"   SSVV : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSSS(m_v4,vanz4);
  ATOOLS::msg.Debugging()<<"   SSSS : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_FFT(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   FFT  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVT(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   VVT  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SST(m_v,vanz);
  ATOOLS::msg.Debugging()<<"   SST  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVVT(m_v4,vanz4);
  ATOOLS::msg.Debugging()<<"   VVVT : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_FFVT(m_v4,vanz4);
  ATOOLS::msg.Debugging()<<"   FFVT : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSST(m_v4,vanz4);
  ATOOLS::msg.Debugging()<<"   SSST : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;

  m_nvertex  = vanz;
  m_n4vertex = vanz4;
  //Print();
  GenerateVertex();
  Print();
  vanz = m_nvertex;
  vanz4 = m_n4vertex;
  //Kill_Off(); to be written....
  //TexOutput();
  m_nvertex = 10000;
  m_n4vertex = 10000;
  
  //should be improved and made dynamical

  if (vanz>m_nvertex) {
      msg.Error()<<"Number of vertices to large"<<endl;
      abort();
  }
  if (vanz4>m_n4vertex) {
      msg.Error()<<"Number of 4 leg vertices to large"<<endl;
      abort();
  }
  
  ATOOLS::msg.Debugging()<<"... done with it ("<<vanz<<")."<<endl
                          <<"... done with the 4 legs ("<<vanz4<<")."<<endl;
  m_nvertex  = vanz;
  m_n4vertex = vanz4;
}

void Vertex::CheckEqual(Flavour** fl,short int& count)
{
  for (short int i=0;i<count-1;++i) {
    if (fl[i][0]==fl[count-1][0] &&
	fl[i][1]==fl[count-1][1] &&
	fl[i][2]==fl[count-1][2]) {
      count--;
      break;
    }
  }
}

void Vertex::Print()
{
  if (!rpa.gen.Debugging()) return;
  //3 legs
  for (int i=0;i<m_nvertex;i++) {
    ATOOLS::msg.Out()<<i+1<<". vertex for :"<<m_v[i].in[0]<<":"<<m_v[i].in[1]<<":"<<m_v[i].in[2];
    if (m_v[i].on) ATOOLS::msg.Out()<<"...On  ";
    else  ATOOLS::msg.Out()<<"...Off ";
    ATOOLS::msg.Out()<<m_v[i].cpl[0]<<";"<<m_v[i].cpl[1];
    ATOOLS::msg.Out()<<"; "<<m_v[i].Color->String();
    ATOOLS::msg.Out()<<"; "<<m_v[i].Lorentz->String()<<endl;
  }
  //4 legs
  for (int i=m_nvertex;i<(m_n4vertex+m_nvertex);i++) {
    if (m_v4[i-m_nvertex].ncf==1) {
    ATOOLS::msg.Out()<<i+1<<". 4 leg vertex for :"<<m_v4[i-m_nvertex].in[0]<<":"
	<<m_v4[i-m_nvertex].in[1]<<":"<<m_v4[i-m_nvertex].in[2]<<":"<<m_v4[i-m_nvertex].in[3];
    if (m_v4[i-m_nvertex].on) ATOOLS::msg.Out()<<"...On  ";
    else  ATOOLS::msg.Out()<<"...Off ";
    ATOOLS::msg.Out()<<m_v4[i-m_nvertex].cpl[0]<<";"<<m_v4[i-m_nvertex].cpl[1];
    ATOOLS::msg.Out()<<"; "<<m_v4[i-m_nvertex].Color->String();
    ATOOLS::msg.Out()<<"; "<<m_v4[i-m_nvertex].Lorentz->String()<<endl;
    }
    else {
      for (short int k=0;k<m_v4[i-m_nvertex].ncf;k++) {
      ATOOLS::msg.Out()<<i+1<<". 4 leg vertex for :"<<m_v4[i-m_nvertex].in[0]<<":"
	  <<m_v4[i-m_nvertex].in[1]<<":"<<m_v4[i-m_nvertex].in[2]<<":"<<m_v4[i-m_nvertex].in[3];
      if (m_v4[i-m_nvertex].on) ATOOLS::msg.Out()<<"...On  ";
      else  ATOOLS::msg.Out()<<"...Off ";
      ATOOLS::msg.Out()<<m_v4[i-m_nvertex].cpl[0]<<";"<<m_v4[i-m_nvertex].cpl[1];
      ATOOLS::msg.Out()<<"; "<<m_v4[i-m_nvertex].Color[k].String();
      if (m_v4[i-m_nvertex].Color[k].Next!=0) ATOOLS::msg.Out()<<" "<<m_v4[i-m_nvertex].Color[k].Next->String();
      ATOOLS::msg.Out()<<"; "<<m_v4[i-m_nvertex].Lorentz[k].String()<<endl;
      }
    }
  }
}
Vertex::~Vertex() {delete[] m_v;delete[] m_v4;}

void Vertex::TexOutput()
{
  mkdir("./tex",448);
  
  ofstream sf;

  String_Tree st;
  
  //Print();
  
  int fmfcount = 0;  
  
  for (int i=0;i<m_nvertex;i++) {
    st.Reset();
    sknot* shelp = st.String2Tree(m_v[i].Str);
    st.Delete(shelp,string("zero"));

    string newstr = st.Tree2Tex(shelp,0);
    
    if (i%90==0) {
      if (i!=0) {
	sf<<"\\end{fmffile}"<<endl;
	sf<<"\\end{document}"<<endl;
	
	sf.close();
      }

      char help[17];
      sprintf(help,"./tex/Vertex_%i.tex",fmfcount);

      sf.open(help);
      sf<<"\\documentclass[a4paper,10pt]{article}"<<endl;
      sf<<"\\newcommand{\\nnb}{\\nonumber}"<<endl;
      sf<<"\\newcommand{\\bea}{\\begin{eqnarray*}}"<<endl;
      sf<<"\\newcommand{\\eea}{\\end{eqnarray*}}"<<endl;
      sf<<"\\newcommand{\\ba}{\\begin{array}}"<<endl;
      sf<<"\\newcommand{\\ea}{\\end{array}}"<<endl;
      sf<<"\\newcommand{\\bt}{\\begin{tabular}}"<<endl;
      sf<<"\\newcommand{\\et}{\\end{tabular}}"<<endl;
      sf<<"\\newcommand{\\m}{-}"<<endl;
      sf<<"\\newcommand{\\p}{+}"<<endl;
      sf<<"\\newcommand{\\ti}{*}"<<endl;
      sf<<"\\usepackage{feynmf}"<<endl;
      sf<<"\\begin{document}"<<endl;
      sf<<"\\begin{center}"<<endl;
      sf<<"\\underline{Amplitudes for job}";
      sf<<"\\end{center}"<<endl;
      sf<<"\\begin{fmffile}{vertpics"<<fmfcount<<"}"<<endl;
      fmfcount++;
    }
    
    sf<<"\\bt{cc}"<<endl;
    sf<<"\\parbox{4cm}{";
    sf<<"\\begin{picture}(100,100)"<<endl;
    sf<<"{\\begin{fmfgraph*}(80,80)"<<endl;
    sf<<" \\fmfleft{l1}\\fmfright{r1,r2}"<<endl; 
    sf<<"\\fmfpen{thin}"<<endl;
    sf<<"\\fmf{phantom}{r1,v1,r2}"<<endl;
    sf<<"\\fmf{phantom}{l1,v1}"<<endl;
    sf<<"\\fmffreeze"<<endl;

    if (m_v[i].in[1].IsVector()) {
      if (m_v[i].in[0].IsFermion() && m_v[i].in[2].IsFermion()) {
	sf<<"\\fmf{plain}{l1,v1,r1}"<<endl; 
	if (m_v[i].in[1]==Flavour(kf::gluon))
	  sf<<"\\fmf{curly}{v1,r2}"<<endl;
	else 
	  sf<<"\\fmf{photon}{v1,r2}"<<endl;
	if (m_v[i].in[1].Charge() != 0)
	  sf<<"\\fmf{phantom_arrow}{v1,r2}"<<endl;
	sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
	sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[0].TexName()<<"$}{l1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[2].TexName()<<"$}{r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[1].TexName()<<",\\;r^\\mu$}{r2}"<<endl;
      }
      
     
      if (m_v[i].in[0].IsVector() && m_v[i].in[2].IsVector() ){
	if (m_v[i].in[1]==Flavour(kf::gluon)) {
	  sf<<"\\fmf{curly}{r1,v1,r2}"<<endl;
	  sf<<"\\fmf{curly}{v1,l1}"<<endl;
	}
	else { 
	  sf<<"\\fmf{photon}{r1,v1,r2}"<<endl;
	  sf<<"\\fmf{photon}{v1,l1}"<<endl;}
	sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
	sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;
	sf<<"\\fmf{phantom_arrow}{v1,r2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[0].TexName()<<",\\;k^\\lambda $}{r1}"<<endl;
	sf<<"\\fmflabel{$"<<(m_v[i].in[2].Bar()).TexName()<<",\\;p^\\nu $}{r2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[1].TexName()<<",\\;r^\\mu $}{l1}"<<endl;
      }   
      
      if (m_v[i].in[0].IsScalar() && m_v[i].in[2].IsScalar()) {
	if (m_v[i].in[1]==Flavour(kf::gluon))
	  sf<<"\\fmf{curly}{v1,r2}"<<endl;
	else 
	  sf<<"\\fmf{photon}{v1,r2}"<<endl;
	sf<<"\\fmf{dashes}{l1,v1,r1}"<<endl; 
	if (m_v[i].in[1].Charge() != 0)
	  sf<<"\\fmf{phantom_arrow}{v1,r2}"<<endl;
	if (m_v[i].in[0].Charge() !=0 || m_v[i].in[2].Charge() !=0){
	  sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
	  sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl; }
	sf<<"\\fmflabel{$"<<m_v[i].in[2].TexName()<<"$}{r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[0].TexName()<<"$}{l1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[1].TexName()<<",\\;r^\\mu $}{r2}"<<endl;
      }
    }
    
    if (m_v[i].in[1].IsScalar()) {
      if (m_v[i].in[0].IsVector() && m_v[i].in[2].IsVector()) {
	sf<<"\\fmf{dashes}{r2,v1}"<<endl; 
	if (m_v[i].in[1].Charge() != 0) {
	  sf<<"\\fmf{phantom_arrow}{r2,v1}"<<endl;
	  sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;}
	if (m_v[i].in[0].Charge() !=0 || m_v[i].in[2].Charge() !=0)
	  {sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
	  sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;}
	sf<<"\\fmf{photon}{l1,v1,r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[0].TexName()<<",\\;p^\\nu $}{l1}"<<endl;	  
	sf<<"\\fmflabel{$"<<m_v[i].in[1].TexName()<<"$}{r2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[2].TexName()<<",\\;r^\\mu $}{r1}"<<endl;
      } 
      
      if (m_v[i].in[0].IsFermion() && m_v[i].in[2].IsFermion()) {
	sf<<"\\fmf{plain}{l1,v1,r1}"<<endl;
	sf<<"\\fmf{dashes}{v1,r2}"<<endl;
	if (m_v[i].in[1].Charge() != 0)
	  sf<<"\\fmf{phantom_arrow}{v1,r2}"<<endl;
	
	sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
	sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[0].TexName()<<"$}{l1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[1].TexName()<<"$}{r2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[2].TexName()<<"$}{r1}"<<endl;
      }
      
      if (m_v[i].in[0].IsScalar() && m_v[i].in[2].IsScalar()) {
	sf<<"\\fmf{dashes}{r1,v1,r2}"<<endl; 
	sf<<"\\fmf{dashes}{v1,l1}"<<endl;
	if (m_v[i].in[0].Charge() != 0) {
	  sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl; 
	  sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;}
	sf<<"\\fmflabel{$"<<m_v[i].in[0].TexName()<<"$}{l1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[2].TexName()<<"$}{r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v[i].in[1].TexName()<<"$}{r2}"<<endl;
      }
    }
    
    sf<<"\\fmfdot{v1}"<<endl;
    sf<<"\\end{fmfgraph*}} "<<endl;
    sf<<"\\end{picture}} &"<<endl;  
    //sf<<"} &"<<endl;  
    sf<<"\\begin{minipage}[tl]{8cm}{"<<endl;
    sf<<"$\\displaystyle "<<newstr<<"$"<<endl;
    sf<<"}\\end{minipage} \\\\"<<endl;  
    sf<<"\\et\\\\[5mm]"<<endl;
  }
  sf<<"\\end{fmffile}"<<endl;
  sf<<"\\end{document}"<<endl;
  
}

void Vertex::AddVertex(Single_Vertex* addv){
  Single_Vertex * oldv=m_v;
  m_v = new Single_Vertex[m_nvertex+1];
  for (int i=0;i<m_nvertex;++i) {
    m_v[i].on=oldv[i].on;
    for (int j=0;j<4;++j) m_v[i].in[j]=oldv[i].in[j];
    for (int j=0;j<4;++j) m_v[i].cpl[j]=oldv[i].cpl[j];
  }
  m_v[m_nvertex].on=addv->on;
  for (int j=0;j<4;++j) m_v[m_nvertex].in[j]=addv->in[j];
  for (int j=0;j<4;++j) m_v[m_nvertex].cpl[j]=addv->cpl[j];
  ++m_nvertex;
  //  for (int i=0;i<m_nvertex;++i) printvertex(&m_v[i]);
}

int Vertex::FindVertex(Single_Vertex* v_tofind)
{
  int nr=-1;
  Flavour help;
  if (v_tofind->nleg==3) {
    for (int j=0;j<3;++j){
      // rotate flavours
      help           =v_tofind->in[0];
      v_tofind->in[0]=v_tofind->in[1];
      v_tofind->in[1]=v_tofind->in[2];
      v_tofind->in[2]=help;
      //    printvertex(v_tofind);
      for (int i=0;i<m_nvertex;++i) {
	//   printvertex(&m_v[i]);
	if (v_tofind->in[0].Kfcode()==m_v[i].in[0].Kfcode())
	  if (v_tofind->in[1].Kfcode()==m_v[i].in[1].Kfcode())
	    if (v_tofind->in[2].Kfcode()==m_v[i].in[2].Kfcode()) {
	    //    v_tofind=&m_v[i];
	      nr=i;
	      //break;
	      v_tofind=&m_v[nr];
	      return nr;
	    }
      }
    }
    ATOOLS::msg.Debugging()<<"Vertex not found!"<<endl;
  }
  else ATOOLS::msg.Debugging()<<"no routine to search for 4 legs"<<endl;
  return 0;
}

ostream& AMEGIC::operator<<(ostream& s, const Single_Vertex& sv)
{
  return s<<'('<<sv.in[0]<<','<<sv.in[1]<<','<<sv.in[2]<<','<<sv.in[3]
          <<") with cpl["<<sv.cpl[0]<<','<<sv.cpl[1]<<','<<sv.cpl[2]<<','<<sv.cpl[3]<<']'
          <<" is "<<((sv.on) ? "on" : "off");
}
 
ostream& AMEGIC::operator<<(ostream& s, const Vertex& v)
{
  s<<"---------------- Vertices --------------------------------"<<endl;

  int n=v.MaxNumber();
  s<<n<<" verticies found"<<endl;
  for (int i=0;i<n;++i)
    s<<(*v[i])<<endl;
  s<<"-----------------------------------------------------------"<<endl;

  return s;
}

void AMEGIC::Single_Vertex2MPI(const Single_Vertex * v , MPI_Single_Vertex & mpi_v) {
  
    
  if (!v) {
    for (int i=0;i<4;++i)
      mpi_v.m_fl[i] = 0;
    
    return; 
  }
  /*
  Lorentz_Function2MPI(v->Lorentz,mpi_v.m_lf);
  Color_Function2MPI(v->Color,mpi_v.m_cf);
  */

  for (int i=0;i<4;++i)
    mpi_v.m_fl[i] = int(v->in[i]);

  /*  
  for (int i=0;i<7;i+=2) {
    mpi_v.m_cpl[i]   = real(v->cpl[i/2]);
    mpi_v.m_cpl[i+1] = imag(v->cpl[i/2]);
  }
  */
}


Single_Vertex * AMEGIC::MPI2Single_Vertex(const MPI_Single_Vertex & mpi_v ) {

  Single_Vertex * v ;
  
  v = new Single_Vertex();
  
  for (int i=0;i<4;++i) {
    v->in[i] = Flavour((kf::code)(abs(mpi_v.m_fl[i])));
    if (mpi_v.m_fl[i]<0) v->in[i]=v->in[i].Bar();
  }

  /*
  v->Lorentz = AMEGIC::MPI2Lorentz_Function(mpi_v.m_lf);
  v->Color   = AMEGIC::MPI2Color_Function(mpi_v.m_cf);
  
  for (int i=0;i<7;i+=2) {
    (v->cpl[i/2]) = Complex(mpi_v.m_cpl[i],mpi_v.m_cpl[i+1]);
  }
  */
  
  return v;

}


std::ostream & AMEGIC::operator<<(std::ostream & s, const MPI_Single_Vertex & sv) {
  s<<sv.m_fl[0]<<","<<sv.m_fl[1]<<","<<sv.m_fl[2]<<","<<sv.m_fl[3];
  return s;
}




