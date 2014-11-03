#include "AMEGIC++/Phasespace/Channel_Generator_UniV.H"
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

Channel_Generator_UniV::Channel_Generator_UniV(int _nin,int _nout,
                                     Point * _plist,int _aid) 
  : Channel_Generator_Base(_nin,_nout,_plist)
{
  IdentifyProps(plist);
  m_idstr = string("");
  GenerateTopos(plist);
  m_aid = _aid;
}

Channel_Generator_UniV::~Channel_Generator_UniV() 
{ 
  for (size_t i=0;i<m_pclist.size();++i) delete m_pclist[i];
}

void Channel_Generator_UniV::GenerateTopos(Point* p)
{
  if (nin != 2) return;
  Point* ph = p->left;
  if (ph->left==0) {
    ph = p->right;
    if (ph->left==0 && p->middle) ph = p->middle;
  }
  if (ph == 0) {
    msg_Error()<<METHOD<<"(): This seems to be a 2->1 process. "
	       <<p->fl<<" -> { "<<p->left->fl<<" "<<p->right->fl<<" }."
	       <<" Aborting."<<std::endl;
    abort();
  }
  switch (tcount) {
  case 0:
    m_topos.push_back(CopyPoints(p));
    break;
  default:
    m_topos.push_back(TransformTS(p));
  }
  MRPScan();
}

Point* Channel_Generator_UniV::CopyPoints(Point* p)
{
  Point* ph = NULL; 
  if (p!=NULL) { 
    ph = new Point(*p);
    m_pclist.push_back(ph);
    ph->middle = NULL;
    ph->left  = CopyPoints(p->left);
    ph->right = CopyPoints(p->right);
    ph->m = 1;
  }
  return ph;
}


Point* Channel_Generator_UniV::TransformTS(Point* p)
{
  Point **props = new Point*[tcount+1];
  Point **propt = new Point*[tcount+1];

  int counts = 0;

  SetProps(p,props,propt,counts);
  Point* ph = new Point(*p); 
  Point* ps = ph;
  m_pclist.push_back(ph);
  ph->m = 1;
  ph->right = CopyPoints(propt[tcount]);
  if (props[tcount]->number<99 && 
      (props[tcount]->fl.IsVector() || props[0]->number>99 || 
       !(props[0]->fl.Strong() || (props[0]->fl.IsLepton() && props[0]->fl.IntCharge()!=0)) )) {
    ph->left = new Point(*propt[tcount-1]);
    m_pclist.push_back(ph->left);
    ph = ph->left;
    for (int i=0;i<tcount-1;i++) {
      ph->m = 0;
      if (i%2==0) {
	ph->left = new Point(*propt[i/2]);
	m_pclist.push_back(ph->left);
	ph->right = CopyPoints(props[tcount-i/2]);
	ph = ph->left;
      }
      else {
	ph->right = new Point(*propt[tcount-(i+1)/2-1]);
	m_pclist.push_back(ph->right);
	ph->left = CopyPoints(props[(i-1)/2]);
	ph = ph->right;
      }
    }
    ph->m = 0;
    ph->left  = CopyPoints(props[(tcount-1)/2]);
    ph->right = CopyPoints(props[(tcount+1)/2]);
  }
  else {
    ph->left = new Point(*propt[0]);  
    m_pclist.push_back(ph->left);
    ph = ph->left;
    for (int i=0;i<tcount-1;i++) {
      ph->m = 0;
      if (i%2==0) {
	ph->right = new Point(*propt[tcount-i/2-1]);
	m_pclist.push_back(ph->right);
	ph->left = CopyPoints(props[i/2]);
	ph = ph->right;
      }
      else {
	ph->left = new Point(*propt[(i+1)/2]);
	m_pclist.push_back(ph->left);
	ph->right = CopyPoints(props[tcount-(i-1)/2]);
	ph = ph->left;
      }
    }
    ph->m = 0;
    ph->left  = CopyPoints(props[tcount/2]);
    ph->right = CopyPoints(props[tcount/2+1]);
  }
  
  delete[] props;
  delete[] propt;
  return ps;
}


void Channel_Generator_UniV::MRPScan()
{}

int Channel_Generator_UniV::MakeChannel(int& echflag,int n,string& path,string& pID)
{  
  if (m_idstr==string("")) CreateChannelID(echflag);

  //add Channel
  char name[22];
  sprintf(name,"C%i_%i",nout,n);

  string filename = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+path+string("/")+
                    string(name)+string(".C");

//   cout<<name<<endl<<endl;
  
  int    rannum = 0;

  ofstream chf;
  chf.open((filename).c_str());

  chf<<"#include "<<'"'<<"PHASIC++/Channels/Single_Channel.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"ATOOLS/Org/Run_Parameter.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"ATOOLS/Org/MyStrStream.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"PHASIC++/Channels/Channel_Elements.H"<<'"'<<endl;  
  chf<<"#include "<<'"'<<"PHASIC++/Channels/Vegas.H"<<'"'<<endl<<endl;  

  chf<<"using namespace PHASIC;"<<endl;  
  chf<<"using namespace ATOOLS;"<<endl<<endl;

  chf<<"namespace PHASIC {"<<endl
     <<"  class "<<name<<" : public Single_Channel {"<<endl;

  //actual Channel
  if (tcount>0) chf <<"    double m_amct,m_alpha,m_ctmax,m_ctmin;"<<endl;
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
	<<"  public:"<<endl
	<<"    "<<name<<"(int,int,Flavour*,Integration_Info * const);"<<endl
	<<"    ~"<<name<<"();"<<endl
	<<"    void   GenerateWeight(Vec4D *,Cut_Data *);"<<endl
	<<"    void   GeneratePoint(Vec4D *,Cut_Data *,double *);"<<endl
	<<"    void   AddPoint(double);"<<endl
	<<"    void   MPISync()                 { p_vegas->MPISync(); }"<<endl
	<<"    void   Optimize()                { p_vegas->Optimize(); } "<<endl
	<<"    void   EndOptimize()             { p_vegas->EndOptimize(); } "<<endl
	<<"    void   WriteOut(std::string pId) { p_vegas->WriteOut(pId); } "<<endl
	<<"    void   ReadIn(std::string pId)   { p_vegas->ReadIn(pId); } "<<endl
        <<"    void   ISRInfo(int &,double &,double &);"<<endl
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
  acount = 0;
  newchannel = 0;
  Step0(0,m_topos[echflag],rannum,chf);
  ClearDeclarations(); 
  chf<<"}"<<endl<<endl;
  
  int rannumber = rannum;
  rannum = 0;
  //Weight
  chf<<"void "<<name<<"::";
  chf<<"GenerateWeight(Vec4D* p,Cut_Data * cuts)"<<endl<<"{"<<endl;
  chf<<"  double wt = 1.;"<<endl;

  maxnumb = 0;
  acount = 0;

  Step0(1,m_topos[echflag],rannum,chf);
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
	<<"  rans  = new double[rannum];"<<endl;
  if (tcount>0) {
    chf	<<"  m_amct  = 1.0+ToType<double>(rpa->gen.Variable(\"AMEGIC_CHANNEL_EPSILON\"));"<<endl
	<<"  m_alpha = ToType<double>(rpa->gen.Variable(\"AMEGIC_CHANNEL_ALPHA\"));"<<endl
	<<"  m_ctmax = 1.;"<<endl
	<<"  m_ctmin = -1.;"<<endl;
  }
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

  chf<<"void "<<name<<"::ISRInfo(int & type,double & mass,double & width)";
  chf<<endl<<"{"<<endl;  
  Step0(2,m_topos[echflag],rannum,chf);
  chf<<"}"<<endl<<endl;

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



void Channel_Generator_UniV::Step0(int flag,Point* p,int& rannum,ofstream& sf) 
{
  if (nin != 2) return;
  Point* ph = p->left;
  if (ph->left==0) {
    ph = p->right;
    if (ph->left==0 && p->middle) ph = p->middle;
  }
  if (ph == 0) {
    msg_Error()<<METHOD<<"(): This seems to be a 2->1 process. "
	       <<p->fl<<" -> { "<<p->left->fl<<" "<<p->right->fl<<" }."
	       <<" Aborting."<<std::endl;
    abort();
  }

  string m = Order(LinkedMasses(ph));
  
  string help("");
  switch (flag) {
  case -11: 
    if (ph->m==0) {
      help+=string("ZS_");
    } 
    else {
      if (!IsZero(ph->fl.Mass())) {
	help+=string("ZR")+ToString(ph->fl.Kfcode())+string("_");
      }
      else help+=string("ZS_");
    }
    {
      int pt;
      double mass = PMassSum(ph,&pt);
      if (pt>1) mass*=pow(1.5,(double)pt-1.);
      help+=ToString((int)mass);
      m_idc.push_back(help);
    }
  case 0: case 1:
    sf<<"  Vec4D p"<<m<<"=p[0]+p[1];"<<endl;
    AddToVariables(flag,m+string("_max"),string("p")+m+string(".Abs2()"),0,sf);
 
    GenerateMassChain(flag,ph,ph,rannum,sf);
    break;
  case 2:
    if (ph->m!=0 && !IsZero(ph->fl.Mass())) {
      sf<<"  type  = 1;"<<endl
	<<"  mass  = Flavour((kf_code)("<<ph->fl.Kfcode()<<")).Mass();"<<endl
	<<"  width = Flavour((kf_code)("<<ph->fl.Kfcode()<<")).Width();"<<endl;
      return;
    }
    {
      int pt;
      double mass = PMassSum(ph,&pt);
      if (pt>1) mass*=pow(1.5,(double)pt-1.);
      sf<<"  type  = 2;"<<endl
	<<"  mass  = "<<ToString(mass)<<";"<<endl
	<<"  width = 0.;"<<endl;
      return;
    }
  }
  vector<string> pin0;
  vector<string> pin1;
  GenerateDecayChain(flag,ph,rannum,sf,pin0,pin1);
}

void Channel_Generator_UniV::GenerateDecayChain(int flag,Point* p,int& rannum,ofstream& sf,
					    vector<string> pin0, vector<string> pin1)
{
  if (p->left==0) return;
//    int sa = AntennaS(p);
//    if (sa>2) return QCDAntenna(flag,p,rannum,sf,sa);

  if (p->m==0) {
    int hi = p->fl.Kfcode();
    string tmstr;
    if (flag>=0 && !IsZero(p->fl.Mass())) {
      tmstr = string("tmass")+ToString(p->number);
      sf<<"  double "<<tmstr<<" = Flavour((kf_code)("<<hi<<")).Mass();"<<endl;
    }
    else tmstr = string("0.");
    string pin0sum(""),pin1sum("");
    if (pin0.size()>0) {
      for (size_t i=0;i<pin0.size();++i) pin0sum+=pin0[i];
 //      cout<<"GDC0: "<<pin0sum<<":";for (int i=0;i<pin0.size();i++)cout<<" "<<pin0[i];cout<<endl;
      pin0sum = Order(pin0sum);
      string help1(""),help2 = pin0[pin0.size()-1];
      for (size_t i=0;i<pin0.size()-1;++i) help1+=pin0[i];
      if (help1.length()>0) help1 = string("p0_") + help1;
      else help1 = string("p[0]");

       if (help2.length()>1) AddToVariables(flag,string("0_")+pin0sum,help1+string("-p")+help2,1,sf);
       else AddToVariables(flag,string("0_")+pin0sum,help1+string("-p[")+GetMassIndex(help2)+string("]"),1,sf);
    }
    if (pin1.size()>0) {
      for (size_t i=0;i<pin1.size();++i) pin1sum+=pin1[i];
 //      cout<<"GDC1: "<<pin1sum<<":";for (int i=0;i<pin1.size();i++)cout<<" "<<pin1[i];cout<<endl;
      pin1sum = Order(pin1sum);
      string help1(""),help2 = pin1[pin1.size()-1];
      for (size_t i=0;i<pin1.size()-1;++i) help1+=pin1[i];
      if (help1.length()>0) help1 = string("p1_") + help1;
      else help1 = string("p[1]");

      if (help2.length()>1) AddToVariables(flag,string("1_")+pin1sum,help1+string("-p")+help2,1,sf);
      else AddToVariables(flag,string("1_")+pin1sum,help1+string("-p[")+GetMassIndex(help2)+string("]"),1,sf);
    }
    pin0sum = string("0_") + pin0sum; 
    pin1sum = string("1_") + pin1sum; 
    string pout0sum = Order(LinkedMasses(p->left));
    string pout1sum = Order(LinkedMasses(p->right));
    string sctmax("m_ctmax");
    string sctmin("m_ctmin");
    if (flag>=0) {
      if (pin0.size()==0 && pout0sum.length()==1 && pin1.size()==0 && pout1sum.length()==1) {
	sf<<"  m_ctmax = Min(cuts->cosmax[0]["<<GetMassIndex(pout0sum)
	  <<"],cuts->cosmax[1]["<<GetMassIndex(pout1sum)<<"]);"<<endl;
	if (nout>2) sf<<"  m_ctmin = Max(cuts->cosmin[0]["<<GetMassIndex(pout0sum)
		      <<"],cuts->cosmin[1]["<<GetMassIndex(pout1sum)<<"]);"<<endl;
      }
      else if (pin0.size()==0 && pin1.size()==0 && pout0sum.length()==1) {
	sf<<"  m_ctmax = cuts->cosmax[0]["<<GetMassIndex(pout0sum)<<"];"<<endl;
	if (nout>2) sf<<"  m_ctmin = cuts->cosmin[0]["<<GetMassIndex(pout0sum)<<"];"<<endl;
      }
      else if (pin0.size()==0 && pin1.size()==0 && pout1sum.length()==1) {
	sf<<"  m_ctmax = cuts->cosmax[1]["<<GetMassIndex(pout1sum)<<"];"<<endl;
	if (nout>2) sf<<"  m_ctmin = cuts->cosmin[1]["<<GetMassIndex(pout1sum)<<"];"<<endl;
      }
      else {
	sctmax = string("1.");
	sctmin = string("-1.");
      }
    }
    string his("");
    switch (flag) {
    case -11:
      if(!IsZero(p->fl.Mass())) his=ToString(hi);
      m_idc.push_back(string("TC")+his+
		      string("_")+pin0sum+string("_")+pin1sum+
		      string("_")+pout0sum+string("_")+pout1sum);
      break;
    case 0:
      sf<<"  CE.TChannelMomenta(";
      if (pin0.size()==0) sf<<"p[0]"; else sf<<"p"<<pin0sum;
      if (pin1.size()==0) sf<<",p[1]"; else sf<<",p"<<pin1sum;
      if (pout0sum.length()==1) sf<<",p["<<GetMassIndex(pout0sum)<<"]"; else sf<<",p"<<pout0sum;
      if (pout1sum.length()==1) sf<<",p["<<GetMassIndex(pout1sum)<<"]"; else sf<<",p"<<pout1sum;
      sf<<",s"<<pout0sum<<",s"<<pout1sum;
      sf<<","<<tmstr<<",m_alpha,"<<sctmax<<","<<sctmin<<",m_amct,0,ran["<<rannum++<<"],ran[";
      sf<<rannum++<<"]);"<<endl;
      break;
    default:
      if(!IsZero(p->fl.Mass())) his=ToString(hi);
      string idh = string("TC")+his+
	string("_")+pin0sum+string("_")+pin1sum+string("_")+pout0sum+string("_")+pout1sum;
      sf<<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT)"<<endl; 
      sf<<"    m_k"<<idh<<"<<CE.TChannelWeight(";
      if (pin0.size()==0) sf<<"p[0]"; else sf<<"p"<<pin0sum;
      if (pin1.size()==0) sf<<",p[1]"; else sf<<",p"<<pin1sum;
      if (pout0sum.length()==1) sf<<",p["<<GetMassIndex(pout0sum)<<"]"; else sf<<",p"<<pout0sum;
      if (pout1sum.length()==1) sf<<",p["<<GetMassIndex(pout1sum)<<"]"; else sf<<",p"<<pout1sum;
      sf<<","<<tmstr<<",m_alpha,"<<sctmax<<","<<sctmin<<",m_amct,0,m_k"<<idh<<"[0],m_k"<<idh<<"[1]);"<<endl;
      sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
      sf<<"  rans["<<rannum++<<"]= m_k"<<idh<<"[0];"<<endl;
      sf<<"  rans["<<rannum++<<"]= m_k"<<idh<<"[1];"<<endl;
    }
    if (p->left->m==0)  pin1.push_back(Order(LinkedMasses(p->right)));
    if (p->right->m==0) pin0.push_back(Order(LinkedMasses(p->left)));
  }
  else {
    Point* l     = p->left;
    Point* r     = p->right;
    string lm    = LinkedMasses(l);
    string rm    = LinkedMasses(r);
    string mummy = Order(lm+rm);
    string moml,momr;
    //Minima
    if (l->left==0) moml = string("p[") + GetMassIndex(lm) + string("]");
    else moml = string("p") + Order(lm);
    if (r->left==0) momr = string("p[") + GetMassIndex(rm) + string("]");
    else momr = string("p") + Order(rm);

    bool first = p->prev->number==0;
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
      }
    }
    else {
      if ((!first) && (r->left==0) && (r->fl.IsVector()) && 
	  (!(r->fl.IsMassive())) && (l->fl.IsFermion()) && m_aid) {
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
  }
  
  GenerateDecayChain(flag,p->left,rannum,sf,pin0,pin1);
  GenerateDecayChain(flag,p->right,rannum,sf,pin0,pin1);
}


bool Channel_Generator_UniV::QCDAntenna(int flag,Point* p,int& rannum,ofstream& sf,int n)
{ 
  string lm    = LinkedMasses(p->left);
  string rm    = LinkedMasses(p->right);
  string mummy = Order(lm+rm);

  if (flag<0) {
    m_idc.push_back(string("AP_")+mummy); 
    return 1;
  }

  bool first = 0;
  if (flag>9 || flag==-1) { first = 1; flag -= 10; }

  switch(flag) {
  case 0:
    sf <<"  double s0"<<acount<<" = cuts->scut["<<GetMassIndex(mummy[0])
       <<"]["<<GetMassIndex(mummy[1])<<"];"<<endl;
    sf <<"  Vec4D ps"<<acount<<"["<<n<<"];"<<endl;
    sf <<"  CE.QCDAPMomenta(ps"<<acount<<",p"<<mummy<<","<<n<<",s0"<<acount<<");"<<endl;
    for (int i=0;i<n;i++) {
      sf<<"  p["<<GetMassIndex(mummy[i])<<"] = ps"<<acount<<"["<<i<<"];"<<endl;
    }
    break;
  default:
    string idh = string("AP_")+mummy;
    sf <<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT) {"<<endl; 
    sf <<"    double s0"<<acount<<" = cuts->scut["<<GetMassIndex(mummy[0])
       <<"]["<<GetMassIndex(mummy[1])<<"];"<<endl;
    sf <<"    Vec4D ps"<<acount<<"["<<n<<"];"<<endl;
    for (int i=0;i<n;i++) {
      sf<<"    ps"<<acount<<"["<<i<<"] = p["<<GetMassIndex(mummy[i])<<"];"<<endl;
    }
    sf <<"    m_k"<<idh<<"<<CE.QCDAPWeight(ps"<<acount<<","<<n<<",s0"<<acount<<");"<<endl;    
    sf <<"  }"<<endl;
    sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
  }
  acount++;
  return 1;
}





void Channel_Generator_UniV::GenerateMassChain(int flag,Point* p,Point* clmp,int& rannum,ofstream& sf)
{
  if (p->left==0) {
    string m = LinkedMasses(p);
    AddToVariables(flag,m,string("ms[")+GetMassIndex(m)+string("]"),0,sf);
    return;
  }
  string lm,rm;
  string mummy,prt;
  string clm = Order(LinkedMasses(clmp));
  lm = Order(LinkedMasses(p->left));
  rm = Order(LinkedMasses(p->right));

  mummy = Order(lm+rm);
  for (size_t i=0;i<clm.length();++i) {
    if (mummy.find(clm[i])==std::string::npos) prt+=clm[i];
  }

  Point* sclmp = clmp;
  //max
  if (prt.length()>1) {
    if (!CheckVariables(flag,prt+string("_min"),0)) CalcSmin(flag,"min",prt,sf,0);
    else if (flag>=0 && CheckVariables(flag,prt,0)) sf<<"  s"<< prt <<"_min = s"<< prt <<";"<<endl;
    AddToVariables(flag,mummy+string("_max"),string("sqr(sqrt(s") + clm +
		   string("_max)-sqrt(s") + prt + string("_min))"),0,sf);    
  }
  if (prt.length()==1) {
    AddToVariables(flag,mummy+string("_max"),string("sqr(sqrt(s") + clm +
		   string("_max)-sqrt(ms[") + GetMassIndex(prt) + string("]))"),0,sf);
  }
  
  if (rm.length()>1 && lm.length()>1) sclmp = p;
  else if (clmp->left==p && LinkedMasses(clmp->right).length()>1) sclmp = p;
  else if (clmp->right==p && LinkedMasses(clmp->left).length()>1) sclmp = p;

  //cout<<":GenerateMassChain: right: "<<LinkedMasses(p->right)<<endl;
  GenerateMassChain(flag,p->right,sclmp,rannum,sf);
  //cout<<":GenerateMassChain: left:  "<<LinkedMasses(p->left)<<endl;
  GenerateMassChain(flag,p->left,sclmp,rannum,sf);

  if (clmp==p) return; 

//   if (sclm!=mummy) {
//     if (prt.length()>1) {
//       if (!CheckVariables(flag,prt+string("_min"),0)) CalcSmin(flag,"min",prt,sf,0);
//       else if (flag>=0 && CheckVariables(flag,prt,0)) sf<<"  s"<< prt <<"_min = s"<< prt <<";"<<endl;
//       AddToVariables(flag,mummy+string("_max"),string("sqr(sqrt(s") + clm +
// 		     string("_max)-sqrt(s") + prt + string("_min))"),0,sf);
//     }
//     if (prt.length()==1) {
//       AddToVariables(flag,mummy+string("_max"),string("sqr(sqrt(s") + clm +
// 		     string("_max)-sqrt(ms[") + prt + string("]))"),0,sf);
//     }
//   }
  //min
  double dth=false;
  CalcSmin(flag,"min",mummy,sf,0);
  if (mummy.length()>2 && flag>=0) {
    sf <<"  s"<< mummy <<"_min = Max(s"<< mummy <<"_min,sqr(sqrt(s"<< lm <<")+sqrt(s"<< rm <<")));"<<endl;
    dth=true;
  }


  double maxpole = -1.;
  double res = ATOOLS::sqr(p->fl.Width()*p->fl.Mass());
  if (p->m>0 && !ATOOLS::IsZero(res) && Massive(p->fl)) maxpole = 1./res;

  int hi = 4;
  if (mummy.length()>2) hi = 2;
  if (maxpole>0.) {
    hi = (p->fl).Kfcode();
    if (flag>=0) sf<<"  Flavour fl"<<mummy<<" = "<<"Flavour((kf_code)("<<hi<<"));"<<endl;
  } 
  string mlexp(".5");
  if (dth) mlexp=string("1.");
  string thexp("1.5");
  dth=false;
  //mlexp=string("0.");
//   if (p->m==0) thexp = string("1.5");
  switch (flag) {
  case -11:
    if (maxpole>0.) {
      m_idc.push_back(string("MP")+ToString(hi)+string("_")+mummy);
    }
    else m_idc.push_back(string("MTH_")+Order(mummy));
    break;
  case 0:
    sf<<"  Vec4D  p"<<mummy<<";"<<endl;
    if (maxpole>0.) {
      sf<<"  double s"<< mummy
	<<" = CE.MassivePropMomenta(fl"<<mummy<<".Mass(),"<<"fl"<<mummy<<".Width(),1,"
	<<"s"<<mummy<<"_min,s"<<mummy<<"_max,ran["<<rannum<<"]);"<<endl;
    }
    else {
      if (!dth) {
	sf<<"  double s"<<mummy<<" = CE.MasslessPropMomenta("<<mlexp<<",s"<<mummy<<"_min,"
	  <<"s"<<mummy<<"_max,ran["<<rannum<<"]);"<<endl;
      } 
      else {
	sf<<"  double s"<<mummy<<" = CE.ThresholdMomenta("<<thexp<<","
	  <<hi<<".*sqrt(s"<<mummy<<"_min),s"<<mummy<<"_min,"
	  <<"s"<<mummy<<"_max,ran["<<rannum<<"]);"<<endl;
      }
    }
    AddToVariables(flag,mummy,string(""),0,sf);
    rannum++;
    break;
  default:
    string s; 
    if (mummy.length()>0) { 
      for (size_t i=0;i<mummy.length()-1;++i) s += string("p[")+GetMassIndex(mummy[i])+string("]+");
      s += string("p[")+GetMassIndex(mummy[mummy.length()-1])+string("]");
    }
    AddToVariables(flag,mummy,s,1,sf);
    AddToVariables(flag,mummy,string("dabs(p")+mummy+string(".Abs2())"),0,sf);
    if (maxpole>0.) {
      sf<<"  wt *= CE.MassivePropWeight(fl"<<mummy<<".Mass(),"<<"fl"<<mummy<<".Width(),1,"
	<<"s"<<mummy<<"_min,s"<<mummy<<"_max,"<<"s"<<mummy<<",rans["<<rannum<<"]);"<<endl;
    }
    else {
      if (!dth) {
	sf<<"  wt *= CE.MasslessPropWeight("<<mlexp<<",s"<<mummy<<"_min,"
	  <<"s"<<mummy<<"_max,s"<<mummy<<",rans["<<rannum<<"]);"<<endl;
      }
      else {
	sf<<"  wt *= CE.ThresholdWeight("<<thexp<<","<<hi<<".*sqrt(s"<<mummy<<"_min),s"<<mummy<<"_min,"
	  <<"s"<<mummy<<"_max,s"<<mummy<<",rans["<<rannum<<"]);"<<endl;
      }
    }
    rannum++;
  }
}


void Channel_Generator_UniV::SetProps(Point* p,Point** props,Point** propt, int& count)
{
  if (p->left==0) return;

  if (p->right->t) {
    props[count] = p->left;
    propt[count] = p->right;
  }
  else {
    if (p->left->t) {
      props[count] = p->right;
      propt[count] = p->left;
    }
    else {
      if (p->right->b == -1 && p->right->number<99) {
	props[count] = p->left;
	propt[count] = p->right;
      }
      else {
	props[count] = p->right;
	propt[count] = p->left;
      }
      return;
    }
  }
  
  count++; 
  SetProps(propt[count-1],props,propt,count);
}

void Channel_Generator_UniV::CalcTSmin(int flag,vector<string>& p,ofstream& sf)
{
  string help;
  for (size_t i=0;i<p.size();++i) {
    if (p[i].length()==1) help += p[i];
  }

  string psum;
  for (size_t i=0;i<p.size();++i) psum += p[i];


  string s;

  if (help.length()>0) {
    if (help.length()==psum.length()) {
      CalcSmin(flag,"min",help,sf,0);
      return;
    }
    else CalcSmin(flag,"min",help,sf,0);    
    s = string("sqr(sqrt(s")+Order(help)+string("_min)");
  }
  else s = string("sqr(");


  for (size_t i=0;i<p.size();++i) {
    if (p[i].length()>1) s += string("+sqrt(s") + Order(p[i]) + string(")");
  }
  s += string(")");

  if (!CheckVariables(flag,Order(psum)+string("_min"),0)) {
    AddToVariables(flag,Order(psum) + string("_min"),s,0,sf);
  }
  else sf<<"  s"<<Order(psum)<<"_min = Max(s"<<Order(psum)<<"_min,"<<s<<");"<<endl;
}



void Channel_Generator_UniV::CalcSmin(int flag,const char* min,string lm,ofstream& sf,Point* p)
{
  // sum_{i<j} sij + (2-n)*sum_i m_i^2

  if (lm.length()>1) {
    AddToVariables(flag,Order(lm) + string("_") + string(min),
		   string("cuts->Getscut(std::string(\"") + Order(lm) + string("\"))"),0,sf);
  }
  else {
    AddToVariables(flag,Order(lm) + string("_") + string(min),
		   string("ms[") + GetMassIndex(lm) + string("]"),0,sf);
  }
  /*  string s("");

  for (short int i=0;i<lm.length();i++) {
    for (short int j=i+1;j<lm.length();j++) 
      s += string("cuts->scut[")+lm[i]+string("][")+lm[j]+string("]+");
  }
  
  if (lm.length()==1) s += string("ms[")+lm[0]+string("]");
  else {
    if (lm.length()==2) s +=string("0");
    else {      
      s += string("(2-")+IString(lm.length())+string(")*(ms[")+lm[0]+string("]");
      for (short int i=1;i<lm.length();i++) s += string("+ms[")+lm[i]+string("]");
      s += string(")");
    }
  }

  if (lm.length()>2 && p!=0) {
    AddToVariables(flag,Order(lm) + string("_") + string(min) + string("1"),s,0,sf);
    //sf<<"  double s"<<Order(lm)<<"_"<<min<<"1 = "<<s<<";"<<endl;
    
    //Additional Testcut    
    string s2;
    s2 = string("sqr(");
    CalcSmin2(p,s2);
    s2 += string(")");
    AddToVariables(flag,Order(lm) + string("_") + string(min) + string("2"),s2,0,sf);
    AddToVariables(flag,Order(lm) + string("_") + string(min),
		   string("Max(s")+Order(lm)+string("_")+string(min)+string("1,s")
		                  +Order(lm)+string("_")+string(min)+string("2)"),0,sf);
  }
  else {
    AddToVariables(flag,lm + string("_") + string(min),s,0,sf);
    }*/
}

string Channel_Generator_UniV::GetFlMass(Point* p)
{
  if (p->left==0) return string("");
  if (p->fl.Mass()>PMassSum(p->left,0)+PMassSum(p->right,0)) {
    return string("Flavour((kf_code)(")+ToString((p->fl).Kfcode())+string(")).Mass()");
  } 
  string h1 = GetFlMass(p->left);
  string h2 = GetFlMass(p->right);
  if (h1.length()==0) return h2;
  if (h2.length()==0) return h1;
  return h1+string("+")+h2;
}

string Channel_Generator_UniV::LinkedMasses(Point* p)
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


void Channel_Generator_UniV::IdentifyProps(Point* _plist)
{
  InitT(_plist);
  Point* endp;
  tcount = 0;
  _plist->prev = 0;
  BackLinks(_plist,endp);
  Point* p = endp;
  if (p->prev != _plist) {
    for (;;) {
      p = p->prev;
      p->t = 1;
      tcount++;
      if (p->prev == _plist) break;
    }
  }
}

void Channel_Generator_UniV::BackLinks(Point* p,Point* &endp)
{
  if ((p->left==0) && (p->right==0)) {
    if (p->b == -1) endp = p;
    return;
  }  
  p->t = 0;
  p->left->prev = p;
  p->right->prev = p;
  BackLinks(p->left,endp);
  BackLinks(p->right,endp);
}

void Channel_Generator_UniV::InitT(Point* p)
{
  p->t = 0;
  if (p->left==0) return;
  InitT(p->left);
  InitT(p->right);
}

string Channel_Generator_UniV::IString(int i)
{
  MyStrStream sstr;
  sstr<<i;
  string istr;
  sstr>>istr;
  return istr;
}

string Channel_Generator_UniV::Order(string s)
{
  int beg = s.find("_");
  if (beg!=-1) {
    return Order(s.substr(0,beg)) + string("_") + Order(s.substr(beg+1));
  }
  if (s[0]>85 || s[0]<='0') return s;

  for (size_t i=0;i<s.length();++i) 
    for (size_t j=i+1;j<s.length();++j) {
      if (s[i]>s[j]) {
	char help = s[i];
	s[i] = s[j];
	s[j] = help;
      } 
    }
  return s;
}

void  Channel_Generator_UniV::AddToVariables(int flag,const string& lhs,const string& rhs,const int& type,
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
    if (rhs!=string("")) {
      declarations[name]=rhs;
      
      if (type == 0) sf<<"  double s";
      else           sf<<"  Vec4D  p";
      sf<<lhso<<" = "<<rhs<<";"<<endl;
    }
    else declarations[name]=string("dummy");
  } 
  else {
    // already exists
    if (rhs != declarations[name]) {
      msg_Error()<<" ERROR in Channel_Generator_UniV::AddToVariables ()"<<endl;
      abort();
    }
  }
}

bool  Channel_Generator_UniV::CheckVariables(int flag,const string& lhs,const int& type)
{
  if (flag<0) return true;
  string lhso = Order(lhs);
  std::string name;
  if (type ==0) name=string("s")+lhso;
  else name=string("p")+lhso;

  Decls::const_iterator cit=declarations.find(name);
  if (cit==declarations.end()) return false;
  return true; 
}

int Channel_Generator_UniV::AntennaS(Point* p)
{
  if (!p->fl.Strong() || p->fl.IsMassive() || p->m!=1) return 0;
  if (p->left==0) return 1;
  int ls = AntennaS(p->left);
  int rs = AntennaS(p->right);
  if (ls==0 || rs==0) return 0;
  return ls+rs;
}

namespace AMEGIC {
  class Compare_String {
  public:
    int operator()(const string & a, const string & b) {
      return (a<b);
    }
  };
}

std::string Channel_Generator_UniV::CreateChannelID(int echflag)
{
  extrachannelflag = echflag;
  int    rannum = 1;
  ofstream chf;
  Step0(-11,m_topos[echflag],rannum,chf);

  string idstr;
  std::sort(m_idc.begin(),m_idc.end(),Compare_String());
  for (String_List::iterator it = m_idc.begin();it!=m_idc.end();++it) {
    if ((*it).find("I")!=string::npos||(*it).find("TC")!=string::npos){
      idstr+=(*it);
      idstr+=string("$");
    }
  }
  m_mapstr = idstr;

  idstr=string("");
  std::sort(m_idc.begin(),m_idc.end(),Compare_String());
  for (String_List::iterator it = m_idc.begin();it!=m_idc.end();++it) {
    idstr+=(*it);
    idstr+=string("$");
  }
  idstr = string("CGND$")+idstr;
  m_idstr = idstr;
  return idstr;
}

double Channel_Generator_UniV::PMassSum(Point* p,int *pt)
{
  int tptl,tptr;
  *pt=0;
  if (!p->left) return 0.;
  double m = 0.;
  if (p->m>0 && p->fl.IsMassive()) {
    m = p->fl.Mass();
  }
  double mc = PMassSum(p->left,&tptl) + PMassSum(p->right,&tptr);
  if (m>mc) *pt=1; 
  else if (tptl+tptr>0) *pt=Max(tptl,tptr)+1;
  return Max(m,mc);  
}
