#include "AMEGIC++/Phasespace/Channel_Generator_Decays.H"
#include "AMEGIC++/Main/Topology.H"
#include "ATOOLS/Phys/Flavour.H"
#include "AMEGIC++/Main/Point.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <algorithm>
#include <stdio.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Channel_Generator_Decays::Channel_Generator_Decays(int _nin,int _nout,
                                     Point * _plist,int _aid) 
  : Channel_Generator_Base(_nin,_nout,_plist)
{
  m_idstr = string("");
  m_aid = _aid;
}

Channel_Generator_Decays::~Channel_Generator_Decays() { }

int Channel_Generator_Decays::MakeChannel(int& echflag,int n,string& path,string& pID)
{  
  if (m_idstr==string("")) CreateChannelID(echflag);
  int oldflag      = echflag;
  extrachannelflag = echflag;

  //add Channel
  char name[22];
  sprintf(name,"CD%i_%i",nout,n);

  if (echflag!=0) sprintf(name,"%s%c",name,'a'+extrachannelflag-1);
  
  string filename = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+path+string("/")+
                    string(name)+string(".C");
  
  int    rannum = 0;

  ofstream chf;
  chf.open((filename).c_str());

  chf<<"#include "<<'"'<<"PHASIC++/Channels/Single_Channel.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"ATOOLS/Org/Run_Parameter.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"PHASIC++/Channels/Channel_Elements.H"<<'"'<<endl;  
  chf<<"#include "<<'"'<<"PHASIC++/Channels/Vegas.H"<<'"'<<endl<<endl;  

  chf<<"using namespace PHASIC;"<<endl;  
  chf<<"using namespace ATOOLS;"<<endl<<endl;

  chf<<"namespace PHASIC {"<<endl
     <<"  class "<<name<<" : public Single_Channel {"<<endl;

  //actual Channel
  if (m_idc.size()>0) {
    chf <<"    Info_Key ";
    bool first=true;
    for (size_t i=0; i<m_idc.size();++i) {
      if (m_idc[i].find("M")==std::string::npos) {
	if (!first) chf<<",";
	chf <<"m_k"<<m_idc[i];
	first=false;
      }
    }
    chf <<";"<<endl;
  }

  chf	<<"    Vegas* p_vegas;"<<endl
	<<"    bool m_optfinished;"<<endl
	<<"  public:"<<endl
	<<"    "<<name<<"(int,int,Flavour*,Integration_Info * const);"<<endl
	<<"    ~"<<name<<"();"<<endl
	<<"    void   GenerateWeight(Vec4D *,Cut_Data *);"<<endl
	<<"    void   GeneratePoint(Vec4D *,Cut_Data *,double *);"<<endl
	<<"    void   AddPoint(double);"<<endl
	<<"    void   MPISync()                 { p_vegas->MPISync(); }"<<endl
	<<"    void   Optimize()                { p_vegas->Optimize(); } "<<endl
	<<"    void   EndOptimize()             { p_vegas->EndOptimize(); m_optfinished=1; } "<<endl
	<<"    bool   OptimizationFinished()    { return m_optfinished; } "<<endl
	<<"    void   WriteOut(std::string pId) { p_vegas->WriteOut(pId); } "<<endl
	<<"    void   ReadIn(std::string pId)   { p_vegas->ReadIn(pId); } "<<endl
	<<"    std::string ChID();"<<endl
	<<"  };"<<endl
	<<"}"<<endl<<endl;

  chf<<"extern "<<'"'<<"C"<<'"'<<" Single_Channel * Getter_"<<name
     <<"(int nin,int nout,Flavour* fl,Integration_Info * const info) {"<<endl
     <<"  return new "<<name<<"(nin,nout,fl,info);"<<endl
     <<"}"<<endl<<endl;

  //Momenta
  chf<<"void "<<name<<"::";
  chf<<"GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)"<<endl;
  chf<<"{"<<endl;
  //chf<<"std::cout<<\""<<name<<"\"<<std::endl;"<<endl;
  chf<<"  double *ran = p_vegas->GeneratePoint(_ran);"<<endl;
  chf<<"  for(int i=0;i<rannum;i++) rans[i]=ran[i];"<<endl;
  Flavour * flav    = new Flavour[nout];  
  int       maxnumb = 0;

  newchannel = 0;
  Step0(0,plist,rannum,chf,flav,maxnumb);
  ClearDeclarations(); 
  extrachannelflag = newchannel;
  chf<<"}"<<endl<<endl;

  int rannumber = rannum;
  rannum = 0;

  echflag          = extrachannelflag;
  extrachannelflag = oldflag;

  //Weight
  chf<<"void "<<name<<"::";
  chf<<"GenerateWeight(Vec4D* p,Cut_Data * cuts)"<<endl<<"{"<<endl;
  chf<<"  double wt = 1.;"<<endl;

  maxnumb = 0;

  Step0(1,plist,rannum,chf,flav,maxnumb);
  ClearDeclarations();
  chf<<"  double vw = p_vegas->GenerateWeight(rans);"<<endl;
  chf<<"  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,"<<nout<<"*3.-4.);"<<endl;
  chf<<endl<<"  weight = wt;"<<endl; 
  chf<<"}"<<endl<<endl;
  
  
  //Constructor
  chf	<<name<<"::"<<name<<"(int nin,int nout,Flavour* fl,Integration_Info * const info)"<<endl
	<<"       : Single_Channel(nin,nout,fl)"<<endl
	<<"{"<<endl
	<<"  name = std::string(\""<<name<<"\");"<<endl
	<<"  rannum = "<<rannumber<<";"<<endl
	<<"  rans  = new double[rannum];"<<endl
	<<"  m_optfinished=0;"<<endl;
  for (size_t i=0; i<m_idc.size();++i) {
    if (m_idc[i].find("M")==string::npos) {
      chf <<"  m_k"<<m_idc[i]<<".Assign(std::string(\""<<m_idc[i]<<"\"),2,0,info);"<<endl;
    }
  }
  chf	<<"  p_vegas = new Vegas(rannum,100,name);"<<endl;
  chf	<<"}"<<endl<<endl;

  //Destructor
  chf	<<name<<"::~"<<name<<"()"<<endl
	<<"{"<<endl;
  chf	<<"  delete p_vegas;"<<endl;
  chf	<<"}"<<endl<<endl;

  chf<<"void "<<name<<"::AddPoint(double Value)";
  chf<<endl<<"{"<<endl;  
  chf<<"  Single_Channel::AddPoint(Value);"<<endl;
  chf<<"  p_vegas->AddPoint(Value,rans);"<<endl;  
  chf<<"}"<<endl;

  chf<<"std::string "<<name<<"::ChID()";
  chf<<endl<<"{"<<endl;  
  chf<<"  return std::string(\""<<m_idstr<<"\");"<<endl;  
  chf<<"}"<<endl;

  chf.close();  

  delete[] flav;
  return rannumber;
}


void Channel_Generator_Decays::Step0(int flag,Point* p,int& rannum,ofstream& sf,
			      Flavour* flav,int& maxnumb) 
{
  if (nin != 1) return;

  string m = Order(LinkedMasses(p));
  string help("");

	if (flag<2) {
	  AddToVariables(flag,m,string("p[0]"),1,sf);
	  AddToVariables(flag,m,string("dabs(p")+m+string(".Abs2())"),0,sf);
    
	  flag += 10;
	  if (!StepS(flag,p,rannum,sf,flav,maxnumb)) {
	    msg_Error()<<"This seems to be a 1->1 process !!!"<<endl
		       <<"  "<<p->fl<<" -> { "<<p->left->fl<<" "<<p->right->fl<<" }. Abort the run."<<endl;
	    abort();
	  }
	}
}

bool Channel_Generator_Decays::StepS(int flag,Point* p,int& rannum,
			      ofstream& sf,Flavour* flav,int& maxnumb)
{ 
  if (p->left==0) return 0;
  
  Point* l     = p->left;
  Point* r     = p->right;
  string lm    = LinkedMasses(l);
  string rm    = LinkedMasses(r);
  string mummy = Order(lm+rm);
  string moml,momr;
  //Minima
  if (l->left==0) moml = string("p[") + GetMassIndex(lm) + string("]");
            else  moml = string("p") + Order(lm);
  if (r->left==0) momr = string("p[") + GetMassIndex(rm) + string("]");
             else momr = string("p") + Order(rm);

  Point** _plist = new Point*[2]; 
  _plist[0] = p->left;
  _plist[1] = p->right;

  // Count resonating props, massless s-channels cannot be resonant....
  for (short int i=0;i<2;i++) {
    if ((_plist[i]->left !=0) && (_plist[i]->fl.IsMassive())) {
      flav[maxnumb] = _plist[i]->fl;
      maxnumb++;
    }
  }

  GenerateMasses(flag,_plist,2,rannum,sf);

  delete[] _plist;


  bool first = 0;
  if (flag>9 || flag==-1) { first = 1; flag -= 10; }

  // Check for decay type.
  if ((!first) && (l->left==0) && (l->fl.IsVector()) && 
      (!(l->fl.IsMassive())) && (r->fl.IsFermion()) && m_aid) {
    switch (flag) {
    case -11: m_idc.push_back(string("AI_")+Order(lm)+string("_")+Order(rm)); break;
    case 0: 
      sf<<"  CE.Anisotropic2Momenta(p"<<Order(mummy)
	<<",s"<<Order(lm)<<",s"<<Order(rm)<<",1.,-1.,1."<<","<<moml<<","<<momr<<","
	<<"ran["<<rannum<<"],ran["<<rannum+1<<"]);"<<endl;
      break;
    default:
       string idh = string("AI_")+Order(lm)+string("_")+Order(rm);
       //sf<<"  std::cout<<\""<<idh<<"\";"<<endl;
        sf<<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT)"<<endl; 
        sf<<"    m_k"<<idh<<"<<CE.Anisotropic2Weight(1,-1.,1.,"<<moml<<","<<momr<<");"<<endl;
        sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;


	//             sf<<"  wt *= CE.Anisotropic2Weight(1,-1.,1.,"<<moml<<",";
	//             sf<<momr<<");"<<endl;
    }
  }
  else {
    if ((!first) && (r->left==0) && (r->fl.IsVector()) && 
	(!(r->fl.IsMassive())) && (l->fl.IsFermion())  && m_aid) {
      //anisotropic decay for left outgoing massless vectors
      switch (flag) {
      case -11: m_idc.push_back(string("AI_")+Order(rm)+string("_")+Order(lm)); break;
      case 0:
	sf<<"  CE.Anisotropic2Momenta(p"<<Order(mummy)
	  <<",s"<<Order(rm)<<",s"<<Order(lm)<<",1.,-1.,1."<<","<<momr<<","<<moml<<","
	  <<"ran["<<rannum<<"],"<<"ran["<<rannum+1<<"]);"<<endl;
	break;	
      default:      
 	string idh = string("AI_")+Order(rm)+string("_")+Order(lm);
	//sf<<"  std::cout<<\""<<idh<<"\";"<<endl;
  	sf<<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT)"<<endl; 
  	sf<<"    m_k"<<idh<<"<<CE.Anisotropic2Weight(1,-1.,1.,"<<momr<<","<<moml<<");"<<endl;
  	sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
	//	sf<<"  wt *= CE.Anisotropic2Weight(1.,-1.,1.,"<<momr<<","<<moml<<");"<<endl;
      }
    }
    else {
      if ((l->number) < (r->number)) {
	switch (flag) {
	case -11: m_idc.push_back(string("I_")+Order(lm)+string("_")+Order(rm)); break;
	case 0:
	  sf<<"  CE.Isotropic2Momenta(p"<<mummy<<",s"<<Order(lm)<<",s"<<Order(rm)
	    <<","<<moml<<","<<momr<<",ran["<<rannum<<"],"<<"ran["<<rannum+1<<"]);"<<endl;
	  break;	
	default:      
 	  string idh = string("I_")+Order(lm)+string("_")+Order(rm);
	  //sf<<"  std::cout<<\""<<idh<<"\";"<<endl;
  	  sf<<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT)"<<endl; 
  	  sf<<"    m_k"<<idh<<"<<CE.Isotropic2Weight("<<moml<<","<<momr<<",m_k"<<idh<<"[0],m_k"<<idh<<"[1]);"<<endl;
  	  sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
	  sf<<"  rans["<<rannum<<"]= m_k"<<idh<<"[0];"<<endl;
	  sf<<"  rans["<<rannum+1<<"]= m_k"<<idh<<"[1];"<<endl;
	  //	  	  sf<<"  wt *= CE.Isotropic2Weight("<<moml<<","<<momr<<");"<<endl;
	}
      }
      else {
	switch (flag) {
	case -11: m_idc.push_back(string("I_")+Order(rm)+string("_")+Order(lm)); break;
	case 0:
	  sf<<"  CE.Isotropic2Momenta(p"<<mummy<<",s"<<Order(rm)<<",s"<<Order(lm)
	    <<","<<momr<<","<<moml<<",ran["<<rannum<<"],"<<"ran["<<rannum+1<<"]);"<<endl;
	  break;	
	default:      
 	  string idh = string("I_")+Order(rm)+string("_")+Order(lm);
	  //sf<<"  std::cout<<\""<<idh<<"\";"<<endl;
  	  sf<<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT)"<<endl; 
  	  sf<<"    m_k"<<idh<<"<<CE.Isotropic2Weight("<<momr<<","<<moml<<",m_k"<<idh<<"[0],m_k"<<idh<<"[1]);"<<endl;
  	  sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
	  sf<<"  rans["<<rannum<<"]= m_k"<<idh<<"[0];"<<endl;
	  sf<<"  rans["<<rannum+1<<"]= m_k"<<idh<<"[1];"<<endl;
	  //	  	  sf<<"  wt *= CE.Isotropic2Weight("<<momr<<","<<moml<<");"<<endl;
	}
      }
    }
  }
  if (flag==0||flag==1) rannum += 2;
  
  StepS(flag,l,rannum,sf,flav,maxnumb);
  StepS(flag,r,rannum,sf,flav,maxnumb);

  return 1;
}



void Channel_Generator_Decays::GenerateMasses(int flag,Point** _plist,int pcount,
				       int& rannum,ofstream& sf)
{
  string * lm    = new string[pcount];
  string * momp  = new string[pcount];
  int    * sflag = new int[pcount];
  string mummy;
  string sum_s_i;
  string help;
  for (short int i=0;i<pcount;i++) {
    lm[i] = LinkedMasses(_plist[i]);
    mummy += lm[i];
    if (_plist[i]->left==0) {
      if (flag==0 || flag==10) AddToVariables(flag,lm[i],string("ms[")+GetMassIndex(lm[i])+string("]"),0,sf);
      //sf<<"  double s"<<lm[i]<<" = ms["<<lm[i]<<"];"<<endl;
      momp[i]  = string("p[") + GetMassIndex(lm[i]) + string("]");
      sflag[i] = 1;
      //sum_s_i  += string("-sqrt(s")+lm[i]+string(")");
      help    += lm[i];
    }
    else {
      CalcSmin(flag,"min",lm[i],sf,_plist[i]);
      momp[i]  = string("p") + Order(lm[i]);
      sflag[i] = 0;
    }
  }
  if (help.length()>0) {
    //CalcSmin(flag,"restmin",help,sf,0);
    //sum_s_i = string("-sqrt(s")+Order(help)+string("_restmin)");
    CalcSmin(flag,"min",help,sf,0);
    sum_s_i = string("-sqrt(s")+Order(help)+string("_min)");
  }
  int hit;
  double maxpole;
  double res;
  Flavour flav;
  string smax;
  for (;;) {
    hit = -1;
    maxpole = -1.;
    for (short int j=0;j<pcount;j++) {
      if (sflag[j]==0) {
	flav = _plist[j]->fl;
	res  = ATOOLS::sqr(flav.Width()*flav.Mass());
	if (!ATOOLS::IsZero(res) && Massive(flav)) {
	  if (1./res>maxpole) {
	    maxpole = 1./res;
	    hit = j;
	  }
	}
	else {
	  if (hit==-1) hit = j;
	}
      }
    }
    if (hit==-1) break;
    smax = string("sqr(sqrt(s")+Order(mummy)+string(")")+sum_s_i;
    
    for (short int j=0;j<pcount;j++) {
      if (sflag[j]==0 && j!=hit) {
	smax  += string("-sqrt(s")+Order(lm[j])+string("_min)");
      }
    }
    smax += string(")");
    
    AddToVariables(flag,lm[hit] +string("_max"),smax,0,sf);
    //sf<<"  double s"<<Order(lm[hit])<<"_max = "<<smax<<endl;

    int hi = 0;
    if (maxpole>0.) {
      hi = (_plist[hit]->fl).Kfcode();
      if (flag>=0) sf<<"  Flavour fl"<<lm[hit]<<" = "<<"Flavour((kf_code)("<<hi<<"));"<<endl;
    } 
    switch (flag) {
    case -11: case -1:
      if (maxpole>0.) {
	m_idc.push_back(string("MP")+ToString(hi)+string("_")+Order(lm[hit]));
      }
      else m_idc.push_back(string("MlP_")+Order(lm[hit]));
      break;
    case 0: case 10:
      sf<<"  Vec4D  "<<momp[hit]<<";"<<endl
	<<"  double s"<<Order(lm[hit])<<";"<<endl;
      if (maxpole>0.) {
	sf<<"  s"<<Order(lm[hit])
	  <<" = CE.MassivePropMomenta(fl"<<lm[hit]<<".Mass(),"<<"fl"<<lm[hit]<<".Width(),1,"
	  <<"s"<<Order(lm[hit])<<"_min,s"<<Order(lm[hit])<<"_max,ran["<<rannum<<"]);"<<endl;
      }
      else {
	sf<<"  s"<<Order(lm[hit])<<" = CE.MasslessPropMomenta(1.,s"<<Order(lm[hit])<<"_min,"
	  <<"s"<<Order(lm[hit])<<"_max,ran["<<rannum<<"]);"<<endl;
      }
      rannum++;
      break;
    default:
      string s(""); 
      for (int i=0;i<(int)lm[hit].length()-1;i++) s += string("p[")+GetMassIndex(lm[hit][i])+string("]+");
      s += string("p[")+GetMassIndex(lm[hit][lm[hit].length()-1])+string("]");
     
      AddToVariables(flag,lm[hit],s,1,sf);
      AddToVariables(flag,lm[hit],string("dabs(")+momp[hit]+string(".Abs2())"),0,sf);
      if (maxpole>0.) {
	sf<<"  wt *= CE.MassivePropWeight(fl"<<lm[hit]<<".Mass(),"<<"fl"<<lm[hit]<<".Width(),1,"
	  <<"s"<<Order(lm[hit])<<"_min,s"<<Order(lm[hit])<<"_max,"<<"s"<<Order(lm[hit])<<",rans["<<rannum<<"]);"<<endl;
      }
      else {
	sf<<"  wt *= CE.MasslessPropWeight(1.,s"<<Order(lm[hit])<<"_min,"
	  <<"s"<<Order(lm[hit])<<"_max,s"<<Order(lm[hit])<<",rans["<<rannum<<"]);"<<endl;
      }
      rannum++;
    }
    sum_s_i  += string("-sqrt(s")+Order(lm[hit])+string(")");
    sflag[hit] = 1;
  }
  delete[] lm;
  delete[] momp;
  delete[] sflag;
}


void Channel_Generator_Decays::CalcSmin(int flag,const char* min,string lm,ofstream& sf,Point* p)
{
  // sum_{i<j} sij + (2-n)*sum_i m_i^2

  AddToVariables(flag,Order(lm) + string("_") + string(min),
		 string("cuts->Getscut(std::string(\"") + Order(lm) + string("\"))"),0,sf);

}


string Channel_Generator_Decays::LinkedMasses(Point* p)
{
  if (p->left==0) {
    char help[4];
    sprintf(help,"%i",0);
    if (p->number<10) help[0]=p->number+48;
    else help[0]=p->number+55;
    return string(help);
  }
  return LinkedMasses(p->left)+LinkedMasses(p->right);
}


string Channel_Generator_Decays::IString(int i)
{
  MyStrStream sstr;
  sstr<<i;
  string istr;
  sstr>>istr;
  return istr;
}

string Channel_Generator_Decays::Order(string s)
{
  int beg = s.find("_");
  if (beg!=-1) {
    return Order(s.substr(0,beg)) + string("_") + Order(s.substr(beg+1));
  }
  if (s[0]>=85 || s[0]<='0') return s;

  for (size_t i=0;i<s.length();i++) 
    for (size_t j=i+1;j<s.length();j++) {
      if (s[i]>s[j]) {
	char help = s[i];
	s[i] = s[j];
	s[j] = help;
      } 
    }
  return s;
}

void  Channel_Generator_Decays::AddToVariables(int flag,const string& lhs,const string& rhs,const int& type,
					ofstream& sf)
{
  if (flag<0) return;
  string lhso = Order(lhs);
  std::string name;
  if (type ==0) name=string("s")+lhso;
  else name=string("p")+lhso;

  Decls::const_iterator cit=declarations.find(name);
  if (cit==declarations.end()) {
    // daoes not exist
    declarations[name]=rhs;

    if (type == 0) sf<<"  double s";
              else sf<<"  Vec4D  p";
    sf<<lhso<<" = "<<rhs<<";"<<endl;
  } 
  else {
    // already exists
    if (rhs != declarations[name]) {
      msg_Error()<<" ERROR in Channel_Generator_Decays::AddToVariables. Abort the run."<<endl;
      abort();
    }
  }
}

namespace AMEGIC {
  class Compare_String {
  public:
    int operator()(const string & a, const string & b) {
      return (a<b);
    }
  };
}

std::string Channel_Generator_Decays::CreateChannelID(int echflag)
{
  extrachannelflag = echflag;
  int    rannum = 1;
  int   maxnumb = 0;
  ofstream chf;
  Flavour * flav    = new Flavour[nout];  
  Step0(-11,plist,rannum,chf,flav,maxnumb);
  delete[] flav;

  string idstr;
  std::sort(m_idc.begin(),m_idc.end(),Compare_String());
  for (String_List::iterator it = m_idc.begin();it!=m_idc.end();++it) {
    idstr+=(*it);
    idstr+=string("$");
  }
  idstr = string("CG0$")+idstr;
  m_idstr = idstr;
  return idstr;
}

double Channel_Generator_Decays::PMassSum(Point* p,vector<int>* kfs)
{
  if (!p->left) return 0.;
  double m = 0.;
  if (p->fl.IsMassive()) {
    m = p->fl.Mass();
    if (kfs) kfs->push_back(p->fl.Kfcode());
  }
  return m + PMassSum(p->left,kfs) + PMassSum(p->right,kfs); 
}
