#include "MODEL/Interaction_Models/Vertex.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

// Constructor and Destructor
Vertex::Vertex(Interaction_Model_Base * _model)
{
  /* 
     use (roughly) notation and Vertices of J. Rosiek, PRD41 (1990) 3464
     pull out common factor -i of all Vertices
  */ 

  // for backward compatibility there allways has to be a dummy vertex at
  // the end of the arrys that's overwritten, then.
  m_v.resize(1); m_v4.resize(1);

  int vanz  = 0;
  int vanz4 = 0;  


  msg_Debugging()<<"   Setting vertices..."<<endl;
  _model->c_FFV(m_v,vanz);
  msg_Debugging()<<"   FFV  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_FFS(m_v,vanz);
  msg_Debugging()<<"   FFS  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVV(m_v,vanz);
  msg_Debugging()<<"   VVV  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSV(m_v,vanz);
  msg_Debugging()<<"   SSV  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVS(m_v,vanz);
  msg_Debugging()<<"   VVS  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSS(m_v,vanz);
  msg_Debugging()<<"   SSS  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVVV(m_v4,vanz4);
  msg_Debugging()<<"   VVVV : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSVV(m_v4,vanz4);
  msg_Debugging()<<"   SSVV : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSSS(m_v4,vanz4);
  msg_Debugging()<<"   SSSS : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_FFT(m_v,vanz);
  msg_Debugging()<<"   FFT  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVT(m_v,vanz);
  msg_Debugging()<<"   VVT  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SST(m_v,vanz);
  msg_Debugging()<<"   SST  : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_VVVT(m_v4,vanz4);
  msg_Debugging()<<"   VVVT : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_FFVT(m_v4,vanz4);
  msg_Debugging()<<"   FFVT : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;
  _model->c_SSST(m_v4,vanz4);
  msg_Debugging()<<"   SSST : vanz, vanz4: "<<vanz<<", "<<vanz4<<endl;


  // delete the dummy vertices at the end of the vectors
  m_v.resize(vanz);  m_v4.resize(vanz4);

  m_nvertex  = m_v.size();
  m_n4vertex = m_v4.size();
  //Print();
  //TexOutput();
  GenerateVertex();
  Print();
  
  msg_Debugging()<<"... done with it ("<<m_nvertex+m_n4vertex<<")."<<endl;
  msg_Tracking()<<"Initialized interaction model of MODEL : "<<m_nvertex+m_n4vertex<<" vertices."<<std::endl;
}

Vertex::~Vertex() { }

class Order_KF_C {
public:
  bool operator()(const Flavour &a,const Flavour &b)
  { return a.Kfcode()<b.Kfcode(); }
};// end of class Order_KF_C

class Order_Anti_C {
public:
  bool operator()(const Flavour &a,const Flavour &b)
  { return a.IsAnti()<b.IsAnti(); }
};// end of class Order_Anti_C

class Sort_Order_Out_Flavour {
public:
  bool operator()(const Single_Vertex &a,const Single_Vertex &b)
  {
    Flavour_Vector fla(&a.in[1],&a.in[a.nleg]), flb(&b.in[1],&b.in[b.nleg]);
    std::sort(fla.begin(),fla.end(),Order_Anti_C());
    std::stable_sort(fla.begin(),fla.end(),Order_KF_C());
    std::sort(flb.begin(),flb.end(),Order_Anti_C());
    std::stable_sort(flb.begin(),flb.end(),Order_KF_C());
    for (Flavour_Vector::const_iterator ait(fla.begin()), bit(flb.begin());
	 ait!=fla.end();++ait,++bit) {
      if (ait->IsAnti()<bit->IsAnti()) return true;
      if (ait->IsAnti()>bit->IsAnti()) return false;
      if (ait->Kfcode()<bit->Kfcode()) return true;
      if (ait->Kfcode()>bit->Kfcode()) return false;
    }
    return false;
  }
};// end of class Sort_Order_Out_Flavour

typedef std::set<Single_Vertex,Sort_Order_Out_Flavour> OSV_Map;

// General Methods
void Vertex::GenerateVertex()
{
  int vanzsave = m_nvertex; 
  int vanz4save = m_n4vertex; 
  
  Single_Vertex dummy;
  
  for (int i=0;i<vanz4save;++i) {
    int hit = 1;
    if (hit) {
      //required by Interaction_Model_ADD due to small couplings
      m_v4[i].on = m_v4[i].CheckCoupling();
      if (m_v4[i].on) { 
	if(m_v4[i].nleg==4) {
	  OSV_Map svs;
	  int id[4]={1,2,3,4};
	  for (size_t k(0);k<4;++k) {
	    if (SetVertex(m_v4[i],dummy,id[0],id[1],id[2],id[3],1)
		&& svs.find(dummy)==svs.end()) { 
	      m_a.push_back(dummy); 
	      svs.insert(dummy); 
	    }
	    if (SetVertex(m_v4[i],dummy,-id[0],-id[3],-id[2],-id[1],1)
		&& svs.find(dummy)==svs.end()) { 
	      m_a.push_back(dummy); 
	      svs.insert(dummy); 
	    }
	    int l=-id[0];
	    for (size_t j(0);j<3;++j) id[j]=id[j+1];
	    id[0]=-id[0];
	    id[3]=l;
	  }
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
			  m_v4.push_back(dummy);
			  m_n4vertex++;
			}
			if (SetVertex(m_v4[i],dummy,-k,-l,-m,-n)) {
			  m_v4.push_back(dummy);
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
      m_v[i].on = m_v[i].CheckCoupling();
      if (m_v[i].nleg==3) {  
	OSV_Map svs;
	int id[3]={1,2,3};
	for (size_t k(0);k<3;++k) {
	  if (SetVertex(m_v[i],dummy,id[0],id[1],id[2],0,1)
	      && svs.find(dummy)==svs.end()) { 
	    m_a.push_back(dummy); 
	    svs.insert(dummy); 
	  }
	  if (SetVertex(m_v[i],dummy,-id[0],-id[2],-id[1],0,1)
	      && svs.find(dummy)==svs.end()) { 
	    m_a.push_back(dummy); 
	    svs.insert(dummy); 
	  }
	  int l=-id[0];
	  for (size_t j(0);j<2;++j) id[j]=id[j+1];
	  id[0]=-id[0];
	  id[2]=l;
	}
	for (short int k=1;k<4;k++) {
	  for (short int l=1;l<4;l++) {
	    if (l!=k) {
	      for (short int m=1;m<4;m++) {
		if (m!=l && m!=k) {
		  if (k!=1) k = -k;
		  if (l==1) l = -l;
		  if (m==1) m = -m;
		  if (SetVertex(m_v[i],dummy,k,l,m)) {
		    m_v.push_back(dummy);
		    m_nvertex++;
		  }
		  if (SetVertex(m_v[i],dummy,-k,-l,-m)) {
		    m_v.push_back(dummy);
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
  for (int i=0;i<m_nvertex;i++) {
    for (size_t j=0;j<m_v[i].cpl.size();j++) {
      if (m_v[i].Coupling(j)!=Complex(0.,0.)) {
	string hstr=m_v[i].cpl[j].String();
	for (size_t k=0;(k=hstr.find("\\"))!=string::npos;) hstr.erase(k,1); 
	if (m_cplmap.find(hstr)==m_cplmap.end()) {
	  m_cplmap[hstr]=m_v[i].Coupling(j);
	}
	else {
	  if (m_cplmap[hstr]!=m_v[i].Coupling(j)) {
	    msg_Error()<<"coupling ID not unique: "<<m_v[i].cpl[j].String()<<endl;
	    abort();
	  }
	}
      }
    }
  }
  for (int i=0;i<m_n4vertex;i++) {
    for (size_t j=0;j<m_v4[i].cpl.size();j++) {
      if (m_v4[i].Coupling(j)!=Complex(0.,0.)) {
	string hstr=m_v4[i].cpl[j].String();
	for (size_t k=0;(k=hstr.find("\\"))!=string::npos;) hstr.erase(k,1); 
	if (m_cplmap.find(hstr)==m_cplmap.end()) {
	  m_cplmap[hstr]=m_v4[i].Coupling(j);
	}
	else {
	  if (m_cplmap[hstr]!=m_v4[i].Coupling(j)) {
	    msg_Error()<<"coupling ID not unique: "<<m_v[i].cpl[j].String()
		       <<" vs. "<<hstr<<endl;
	    abort();
	  }
	}
      }
    }
  }
}


int Vertex::CheckExistence(Single_Vertex& probe)
{ // Checks if a vertex with the same flavours at the same legs already
  // exists; RETURNS TRUE IF NO SUCH VERTEX EXISTS. 

  // check either m_v or m_v4 depending on the probe's number of legs
  if (probe.nleg==4)
    { if (find(m_v4.begin(), m_v4.end(), probe)==m_v4.end()) { return 1; }
                                                        else { return 0; }}

  if (probe.nleg==3) 
    { if (find(m_v.begin(), m_v.end(), probe) ==m_v.end()) { return 1; }
                                                      else { return 0; }}
  // Default return: Vertex not found
  return 1;
}

int Vertex::FermionRule(Single_Vertex& probe)
{
  // fermionic: particles left, anti-particles right
  int hit = 1;
  if (probe.in[1].IsFermion() && !probe.in[1].IsAnti() && !probe.in[1].Majorana()) hit = 0;
  if (probe.in[2].IsFermion() &&  probe.in[2].IsAnti() && !probe.in[2].Majorana()) hit = 0;
  
  //FNV chargino interactions
  if (hit==0) 
    if (probe.in[0].IsIno() || probe.in[1].IsIno() || probe.in[2].IsIno()) hit=1;
  
  return hit;
}

int Vertex::SetVertex(Single_Vertex& orig, Single_Vertex& probe, int i0, int i1, int i2, int i3,int mode)
{
  probe = orig;

  if ((orig.dec>0 && orig.dec&4) &&
      !((i0==1 && i1==2 && i2==3) || 
	(i0==-2 && i1==3 && i2==-1) ||
	(i0==-3 && i1==-1 && i2==2))) return 0;
  if ((orig.dec>0 && orig.dec&2) && !(i0==1 && i1==2 && i2==3)) return 0;
  
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
  if (mode==1) {
    if (find(m_a.begin(), m_a.end(), probe)!=m_a.end()) return 0;
  } 
  else {
  if (CheckExistence(probe)==0) return 0;
  if (probe.nleg==3) {
    if (FermionRule(probe)==0) return 0;}
  }
  int hc = 0;

  int cnt = 0;
  for(int i=0;i<orig.nleg;i++) {
    if(orig.in[i]!=orig.in[i].Bar()) cnt++;
  }
  if(cnt>0){
    Flavour *flavlist= new Flavour[cnt];
    int *flaglist= new int[cnt];
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
    delete[] flavlist;
    delete[] flaglist;
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
	if (orig.in[i]==Flavour(kf_A0)) conjugate *= -1;
      }
      
      if (orig.Lorentz.front()->Type()=="SSV" ||
	  orig.Lorentz.front()->Type()=="Gauge3" || 
	  orig.Lorentz.front()->Type()=="AGauge3") conjugate *= -1;
      
      if (conjugate==-1) {
	for (short int i=0;i<4;i++) probe.cpl[i] = -probe.cpl[i];
      }

      probe.Color.front().Conjugate();

       if (probe.Lorentz.front()->String()=="1") {
	//exchange left and right
	Kabbala help = probe.cpl[0];
	probe.cpl[0] = probe.cpl[1];
	probe.cpl[1] = help;
       }
    }
  } 
  if (orig.nleg==3 && (orig.Lorentz.front()->Type()=="Gamma")) {
    //exchange left and right for 'barred' FFV vertices
    if ((i0==-1 && (!orig.in[0].SelfAnti() || (!orig.in[2].SelfAnti())))  || 
	(i0==-3 && (!orig.in[2].SelfAnti() || (!orig.in[0].SelfAnti())))) {
      Kabbala help = -probe.cpl[0];
      probe.cpl[0] = -probe.cpl[1];
      probe.cpl[1] = help;
    }
  } 
  if (probe.dec>0 && probe.in[2].IsDummy()) 
    for (short int i=0;i<4;i++) probe.cpl[i] = -probe.cpl[i];
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

  for (size_t i=0;i<probe.Color.size();i++) {
    ColorExchange(&probe.Color[i],newIndex[0],newIndex[1],newIndex[2],newIndex[3]);
  }
  for (size_t i=0;i<probe.Lorentz.size();i++)
  LorentzExchange(probe.Lorentz[i],newIndex[0],newIndex[1],newIndex[2],newIndex[3]);
  
  return 1;
}
    
void Vertex::ColorExchange(Color_Function* colfunc,int new0,int new1,int new2,int new3)
{
  //T[0,1,2] -> T[new0,new1,new2]
  int  partarg[3]  ={-1,-1,-1};
  char strarg[3]   ={-1,-1,-1};
  int  partargn[3] ={-1,-1,-1};
  char strargn[3]  ={-1,-1,-1};
 
  for (short int i=0;i<3;i++) {
    if (colfunc->Type()==cf::D && i==2) break;
    switch (colfunc->ParticleArg(i)) {
    case 0: 
      partarg[i] = new0;
      strarg[i]  = new0+48;
      break;
    case 1: 
      partarg[i] = new1;
      strarg[i]  = new1+48;
      break;
    case 2: 
      partarg[i] = new2;
      strarg[i]  = new2+48;
      break;
    case 3: 
      partarg[i] = new3;
      strarg[i]  = new3+48;
      break;
    case 4: 
      partarg[i] = 4;
      strarg[i]  = '4';
      break;
    }
  }
  //T[0,1,4]T[3,4,2] -> T[new0,new1,4]T[new3,4,new2]
  if (colfunc->Next()) {
  for (short int i=0;i<3;i++) {
      if (colfunc->Next()->Type()==cf::D && i==2) break;
      switch (colfunc->Next()->ParticleArg(i)) {
      case 0: 
	partargn[i] = new0;
	strargn[i]  = new0+48;
	break;
      case 1: 
	partargn[i] = new1;
	strargn[i]  = new1+48;
	break;
      case 2: 
	partargn[i] = new2;
	strargn[i]  = new2+48;
	break;
      case 3: 
	partargn[i] = new3;
	strargn[i]  = new3+48;
	break;
      case 4: 
	partargn[i] = 4;
	strargn[i]  = '4';
	break;
      }
    }
  }
  colfunc->SetStringArg(strarg[0],strarg[1],strarg[2]);
  colfunc->SetParticleArg(partarg[0],partarg[1],partarg[2]);
  if (colfunc->Next()) {
    colfunc->Next()->SetStringArg(strargn[0],strargn[1],strargn[2]);
    colfunc->Next()->SetParticleArg(partargn[0],partargn[1],partargn[2]);
  }
}

void Vertex::LorentzExchange(Lorentz_Function* lorfunc,int new0,int new1,int new2,int new3)
{
  int partarg[4]={-1,-1,-1,-1};
  for (short int i=0;i<lorfunc->NofIndex();i++) {
    switch (lorfunc->ParticleArg(i)) {
      case 0: partarg[i] = new0;break;
      case 1: partarg[i] = new1;break;
      case 2: partarg[i] = new2;break;
      case 3: partarg[i] = new3;break;
    }
  }
  lorfunc->SetParticleArg(partarg[0],partarg[1],partarg[2],partarg[3]);
}

void Vertex::CheckEqual(Flavour** fl,short int& count)
{
  for (short int i=0;i<count-1;++i) {
    if (fl[i][0]==fl[count-1][0] &&
	fl[i][1]==fl[count-1][1] &&
	fl[i][2]==fl[count-1][2]) {
      --count;
      break;
    }
  }
}

int Vertex::FindVertex(Single_Vertex* v_tofind)
{
  // find a vertex in m_v - does not search in m_v4
  int nr=-1;
  Flavour help;
  if (v_tofind->nleg==3) {
    for (int j=0;j<3;++j){
      // rotate flavours
      help           =v_tofind->in[0];
      v_tofind->in[0]=v_tofind->in[1];
      v_tofind->in[1]=v_tofind->in[2];
      v_tofind->in[2]=help;

      // find vertex with the rotated flavours in m_v
      for (int i=0;i<m_nvertex;++i) {
	if (v_tofind->in[0].Kfcode()==m_v[i].in[0].Kfcode())
	  if (v_tofind->in[1].Kfcode()==m_v[i].in[1].Kfcode())
	    if (v_tofind->in[2].Kfcode()==m_v[i].in[2].Kfcode()) {
	      // vertex found in m_v
	      nr=i;
	      v_tofind=&m_v[nr];
	      return nr;
	    }
      }
    }
  }
  return 0;
}


// Output methods
void Vertex::Print()
{
  if (!msg_LevelIsDebugging()) return;
  //3 legs
  for (int i=0;i<m_nvertex;i++) {
    msg_Out()<<i+1<<". vertex for :"<<m_v[i].in[0]<<":"<<m_v[i].in[1]<<":"<<m_v[i].in[2];
    if (m_v[i].on) msg_Out()<<"...On  ";
    else  msg_Out()<<"...Off ";
    msg_Out()<<m_v[i].Coupling(0)<<";"<<m_v[i].Coupling(1);
    msg_Out()<<"; "<<m_v[i].Color.front().String();
    msg_Out()<<"; "<<m_v[i].Lorentz.front()->String()<<endl;
  }
  //4 legs
  for (int i=m_nvertex;i<(m_n4vertex+m_nvertex);i++) {
    if (m_v4[i-m_nvertex].Color.size()==1) {
      msg_Out()<<i+1<<". 4 leg vertex for :"<<m_v4[i-m_nvertex].in[0]<<":"
	       <<m_v4[i-m_nvertex].in[1]<<":"<<m_v4[i-m_nvertex].in[2]<<":"<<m_v4[i-m_nvertex].in[3];
      if (m_v4[i-m_nvertex].on) msg_Out()<<"...On  ";
      else  
	msg_Out()<<"...Off ";
      msg_Out()<<m_v4[i-m_nvertex].Coupling(0)<<";"<<m_v4[i-m_nvertex].Coupling(1);
      msg_Out()<<"; "<<m_v4[i-m_nvertex].Color.front().String();
      msg_Out()<<"; "<<m_v4[i-m_nvertex].Lorentz.front()->String()<<endl;
    }
    else {
      for (size_t k=0;k<m_v4[i-m_nvertex].Color.size();k++) {
	msg_Out()<<i+1<<". 4 leg vertex for :"<<m_v4[i-m_nvertex].in[0]<<":"
		 <<m_v4[i-m_nvertex].in[1]<<":"<<m_v4[i-m_nvertex].in[2]<<":"<<m_v4[i-m_nvertex].in[3];
	if (m_v4[i-m_nvertex].on) 
	  msg_Out()<<"...On  ";
	else  
	  msg_Out()<<"...Off ";
	msg_Out()<<m_v4[i-m_nvertex].Coupling(0)<<";"<<m_v4[i-m_nvertex].Coupling(1);
	msg_Out()<<"; "<<m_v4[i-m_nvertex].Color[k].String();
	if (m_v4[i-m_nvertex].Color[k].Next()!=0) 
	  msg_Out()<<" "<<m_v4[i-m_nvertex].Color[k].Next()->String();
	msg_Out()<<"; "<<m_v4[i-m_nvertex].Lorentz[k]->String()<<endl;
      }
    }
  }
  msg_Out()<<"Couplings:-------------------"<<endl;
  for (tscmap::iterator mit=m_cplmap.begin();mit!=m_cplmap.end();mit++)
    msg_Out()<<mit->first<<endl<<"         "<<mit->second<<endl; 
}

void Vertex::TexOutput()
{
  /*
  ATOOLS::MakeDir("./tex");
  
  remove("./tex/Vertex_*");  

  ofstream sf;
  
  String_Tree st;
  
  //Print();
  
  int fmfcount = 0;  
  for (int i=0;i<m_nvertex;i++) {
    if (m_v[i].Coupling(0)!=Complex(0.,0.) && m_v[i].Coupling(1)!=Complex(0.,0.)) {
    st.Reset();
    sknot* shelp = st.String2Tree(m_v[i].Str);
    //st.Delete(shelp,string("zero"));
    string newstr = st.Tree2Tex(shelp,0);
    
    if (i%200==0) {
      if (i!=0) {
	sf<<"\\end{fmffile}"<<endl;
	sf<<"\\end{document}"<<endl;
	
	sf.close();
      }
      
      char help[30];
      sprintf(help,"./tex/Vertex3_%i.tex",fmfcount);

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
	if (m_v[i].in[1]==Flavour(kf_gluon))
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
	if (m_v[i].in[1]==Flavour(kf_gluon)) {
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
	if (m_v[i].in[1]==Flavour(kf_gluon))
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
      //end of is vector 
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
      //end of isscalar
    }
    
    sf<<"\\fmfdot{v1}"<<endl;
    sf<<"\\end{fmfgraph*}} "<<endl;
    sf<<"\\end{picture}} &"<<endl; 
    sf<<"\\begin{minipage}[tl]{8cm}{"<<endl; 
    sf<<"$\\displaystyle "<<newstr<<"$"<<endl;
    sf<<"}\\end{minipage} \\\\"<<endl;  
    sf<<"\\et\\\\[5mm]"<<endl;
    }
  }
  sf<<"\\end{fmffile}"<<endl;
  sf<<"\\end{document}"<<endl;
  sf.close();  

  for (int i=0;i<m_n4vertex;i++) {
    if (m_v4[i].Coupling(0)!=Complex(0.,0.) && m_v4[i].Coupling(1)!=Complex(0.,0.)) {
    st.Reset();
    sknot* shelp = st.String2Tree(m_v4[i].Str);
    //st.Delete(shelp,string("zero"));
   
    string newstr = st.Tree2Tex(shelp,0);
    
    if (i%200==0) {
      if (i!=0) {
	sf<<"\\end{fmffile}"<<endl;
	sf<<"\\end{document}"<<endl;
	
	sf.close();
      }
      
      char help[30];
      sprintf(help,"./tex/Vertex4_%i.tex",fmfcount);
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
    sf<<" \\fmfleft{l1,l2}\\fmfright{r1,r2}"<<endl; 
    sf<<"\\fmfpen{thin}"<<endl;
    sf<<"\\fmf{phantom}{r1,v1,r2}"<<endl;
    sf<<"\\fmf{phantom}{l1,v1,l2}"<<endl;
    sf<<"\\fmffreeze"<<endl;
    
    if (m_v4[i].in[0].IsVector()) {
      if (m_v4[i].in[1].IsVector()) {
      if (m_v4[i].in[1]==Flavour(kf_gluon)) {
	sf<<"\\fmf{curly}{r1,v1,r2}"<<endl;
	sf<<"\\fmf{curly}{l1,v1,l2}"<<endl;
      }
      else { 
	sf<<"\\fmf{photon}{r1,v1,r2}"<<endl;
	sf<<"\\fmf{photon}{l1,v1,l2}"<<endl;}
      sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
      sf<<"\\fmf{phantom_arrow}{v1,l2}"<<endl;
      sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;
      sf<<"\\fmf{phantom_arrow}{v1,r2}"<<endl;
      sf<<"\\fmflabel{$"<<m_v4[i].in[0].TexName()<<",\\;k^\\lambda $}{l1}"<<endl;
      sf<<"\\fmflabel{$"<<m_v4[i].in[1].TexName()<<",\\;r^\\mu $}{l2}"<<endl;
      sf<<"\\fmflabel{$"<<m_v4[i].in[2].TexName()<<",\\;p^\\nu $}{r1}"<<endl;
      sf<<"\\fmflabel{$"<<m_v4[i].in[3].TexName()<<",\\;l^\\kappa $}{r2}"<<endl;
      }
      
      if (m_v4[i].in[1].IsScalar()) {
	sf<<"\\fmf{photon}{l1,v1,r1}"<<endl;
	sf<<"\\fmf{dashes}{l2,v1,r2}"<<endl;
	if (m_v4[i].in[0].Charge()!=0) {
	  sf<<"\\fmf{phantom_arrow}{l1,v1}"<<endl;
	  sf<<"\\fmf{phantom_arrow}{v1,r1}"<<endl;}
	if (m_v4[i].in[1].Charge()!=0) {
	  sf<<"\\fmf{phantom_arrow}{v1,l2}"<<endl;
	  sf<<"\\fmf{phantom_arrow}{v1,r2}"<<endl;}
	sf<<"\\fmflabel{$"<<m_v4[i].in[0].TexName()<<",\\;k^\\lambda $}{l1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[1].TexName()<<"$}{l2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[2].TexName()<<"$}{r2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[3].TexName()<<",\\;l^\\kappa $}{r1}"<<endl;
      }
      //end of Isvector 
    }  
    
    if (m_v4[i].in[0].IsScalar()) {
      if (m_v4[i].in[1].IsScalar() && m_v4[i].in[2].IsScalar() && m_v4[i].in[3].IsScalar()) {
	sf<<"\\fmf{dashes}{l1,v1,r1}"<<endl; 
	sf<<"\\fmf{dashes}{r2,v1,l2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[0].TexName()<<"$}{l1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[1].TexName()<<"$}{l2}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[2].TexName()<<"$}{r1}"<<endl;
	sf<<"\\fmflabel{$"<<m_v4[i].in[3].TexName()<<"$}{r2}"<<endl;
      }
      //end of Isscalar
    }
    
    sf<<"\\fmfdot{v1}"<<endl;
    sf<<"\\end{fmfgraph*}} "<<endl;
    sf<<"\\end{picture}} &"<<endl;  
    sf<<"\\begin{minipage}[tl]{8cm}{"<<endl;
    sf<<"$\\displaystyle "<<newstr<<"$"<<endl;
    sf<<"}\\end{minipage} \\\\"<<endl;  
    sf<<"\\et\\\\[5mm]"<<endl;
    }
  }
  sf<<"\\end{fmffile}"<<endl;
  sf<<"\\end{document}"<<endl;
  sf.close();
  */
}
 
ostream& MODEL::operator<<(ostream& s, const Vertex& v)
{
  s<<"---------------- Vertices --------------------------------"<<endl;

  int n=v.MaxNumber();
  s<<n<<" verticies found"<<endl;
  for (int i=0;i<n;++i)
    s<<(*v[i])<<endl;
  s<<"-----------------------------------------------------------"<<endl;

  return s;
}
