#include "Channel_Generator3_NP.H"
#include "Topology.H"
#include "Flavour.H"
#include "Point.H"
#include "Message.H"
#include "MathTools.H"
#include "MyStrStream.H"

#include <algorithm>
#include <stdio.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Channel_Generator3_NP::Channel_Generator3_NP(int _nin,int _nout,Flavour * _fl,
                                     Point * _plist,int _aid) 
  : Channel_Generator_Base(_nin,_nout,_fl,_plist)
{
  IdentifyProps(plist);
  m_idstr = string("");
  GenerateTopos(plist);
  m_aid = _aid;
}

Channel_Generator3_NP::~Channel_Generator3_NP() 
{ 
  for (int i=0;i<m_pclist.size();i++) delete m_pclist[i];
}

void Channel_Generator3_NP::GenerateTopos(Point* p)
{
  if (nin != 2) return;
  Point* ph = p->left;
  if (ph->left==0) {
    ph = p->right;
    if (ph->left==0 && p->middle) ph = p->middle;
  }
  if (ph == 0) {
    msg.Error()<<"This seems to be a 2->1 process !!!"<<endl
	       <<"  "<<p->fl<<" -> { "<<p->left->fl<<" "<<p->right->fl<<" }"<<endl;
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

Point* Channel_Generator3_NP::CopyPoints(Point* p)
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


Point* Channel_Generator3_NP::TransformTS(Point* p)
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

Point* Channel_Generator3_NP::GetMirrorTopo(Point* p)
{
  Point* ph = NULL; 
  if (p!=NULL) { 
    ph = new Point(*p);
    m_pclist.push_back(ph);
    ph->middle = NULL;
    if (ph->m!=2) {
      ph->left  = GetMirrorTopo(p->left);
      ph->right = GetMirrorTopo(p->right);
    }
    else {
      ph->left  = GetMirrorTopo(p->right);
      ph->right = GetMirrorTopo(p->left);
    }
  }
  return ph;
}

int Channel_Generator3_NP::MarkNP(Point* p)
{
  if (p->left==0) return 0;
  if (p->m>0 && p->fl.IsMassive() && p->fl.Width()>0.) {
    if (p->fl.Mass()<PMassSum(p,0)) {
      p->m=2;
      if (p->left->left==0 || p->right->left==0) return 1;
      return 2;
    }
  }
  return MarkNP(p->left) + MarkNP(p->right);
}

void Channel_Generator3_NP::MRPScan()
{
  Point* ph = m_topos[0]->left;
  if (ph->left==0) {
    ph = m_topos[0]->right;
    if (ph->left==0 && m_topos[0]->middle) ph = m_topos[0]->middle;
  }
  switch (MarkNP(ph->left)) {
  case 2:
    m_topos.push_back(GetMirrorTopo(m_topos[0]));
  case 1:
    return;
  }
  switch (MarkNP(ph->right)) {
  case 2:
    m_topos.push_back(GetMirrorTopo(m_topos[0]));
  case 1:
    return;
  }
  m_topos.clear();
}

int Channel_Generator3_NP::MakeChannel(int& echflag,int n,string& path,string& pID)
{  
  if (m_idstr==string("")) CreateChannelID(echflag);

  //add Channel
  char name[22];
  sprintf(name,"C%i",n);

  string filename = string("Process/")+path+string("/")+
                    string(name)+string(".C");
  
  ifstream from;
  from.open((filename).c_str());
  
  int    rannum = 1;

  ofstream chf;
  chf.open((filename).c_str());

  chf<<"#include "<<'"'<<"Single_Channel.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"Multi_Channel.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"Phase_Space_Handler.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"Run_Parameter.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"Channel_Elements.H"<<'"'<<endl<<endl;  

  chf<<"using namespace PHASIC;"<<endl;  
  chf<<"using namespace ATOOLS;"<<endl<<endl;

  chf<<"namespace PHASIC {"<<endl
     <<"  class "<<name<<" : public Single_Channel {"<<endl;

  //actual Channel
  if (tcount>0) chf <<"    double m_amct,m_alpha,m_ctmax,m_ctmin;"<<endl;
  if (m_idc.size()>0) {
    chf <<"    Info_Key ";
    bool first=true;
    for (int i=0; i<m_idc.size();i++) {
      if (m_idc[i].find("M")==-1) {
	if (!first) chf<<",";
	chf <<"m_k"<<m_idc[i];
	first=false;
      }
    }
    chf <<";"<<endl;
  }

  chf	<<"    Info_Key m_spkey,m_sppkey,m_y1key,m_y2key;"<<endl
	<<"    Multi_Channel *p_isrchannels,*p_beamchannels;"<<endl
	<<"    int m_beammode,m_isrmode;"<<endl
	<<"  public:"<<endl
	<<"    "<<name<<"(int,int,Flavour*,Integration_Info * const,Phase_Space_Handler*);"<<endl
    //<<"       : Single_Channel(nin,nout,fl)"<<endl
    //<<"    { name = std::string(\""<<name<<"\"); };"<<endl
	<<"    void   GenerateWeight(Vec4D *,Cut_Data *);"<<endl
	<<"    void   GeneratePoint(Vec4D *,Cut_Data *,double *);"<<endl
    //	<<"    int    CountResonances(Flavour*&);"<<endl
    //<<"    void   ISRInfo(int &,double &,double &);"<<endl
	<<"    std::string ChID();"<<endl
	<<"  };"<<endl
	<<"}"<<endl<<endl;

  chf<<"extern "<<'"'<<"C"<<'"'<<" Single_Channel * Getter_"<<name
     <<"(int nin,int nout,Flavour* fl,Integration_Info * const info,Phase_Space_Handler* psh) {"<<endl
     <<"  return new "<<name<<"(nin,nout,fl,info,psh);"<<endl
     <<"}"<<endl<<endl;


  //Momenta
  chf<<"void "<<name<<"::";
  chf<<"GeneratePoint(Vec4D * p,Cut_Data * cuts,double * ran)"<<endl;
  chf<<"{"<<endl;
  //chf<<"std::cout<<\""<<name<<"\"<<std::endl;"<<endl;
  Flavour * flav    = new Flavour[nout];  
  int       maxnumb = 0;
  
  acount = 0;
  newchannel = 0;
  Step0(0,m_topos[echflag],rannum,chf);
  ClearDeclarations(); 
  chf<<"}"<<endl<<endl;
  
  int rannumber = rannum;

  //Weight
  chf<<"void "<<name<<"::";
  chf<<"GenerateWeight(Vec4D* p,Cut_Data * cuts)"<<endl<<"{"<<endl;
  chf<<"  double wt = 1.;"<<endl;

  maxnumb = 0;
  acount = 0;

  Step0(1,m_topos[echflag],rannum,chf);
  ClearDeclarations();
  chf<<"  if (wt!=0.) wt = 1./wt/pow(2.*M_PI,"<<nout<<"*3.-4.);"<<endl;
  chf<<endl<<"  weight = wt;"<<endl; 
  chf<<"}"<<endl<<endl;
    
  //Constructor
  chf	<<name<<"::"<<name<<"(int nin,int nout,Flavour* fl,Integration_Info * const info,Phase_Space_Handler* psh)"<<endl
	<<"       : Single_Channel(nin,nout,fl)"<<endl
	<<"{"<<endl
	<<"  name = std::string(\""<<name<<"\");"<<endl
	<<"  rannum = "<<rannumber<<";"<<endl
	<<"  rans  = new double[rannum];"<<endl;
  if (tcount>0) {
    chf	<<"  m_amct  = 1.;"<<endl
	<<"  m_alpha = .9;"<<endl
	<<"  m_ctmax = 1.;"<<endl
	<<"  m_ctmin = -1.;"<<endl;
  }
  for (int i=0; i<m_idc.size();i++) {
    if (m_idc[i].find("M")==string::npos) {
      chf <<"  m_k"<<m_idc[i]<<".Assign(std::string(\""<<m_idc[i]<<"\"),0,0,info);"<<endl;
    }
  }
  chf   <<"  m_sppkey.Assign(\"s' isr\",5,0,info);"<<endl
	<<"  m_spkey.Assign(\"s' beam\",5,0,info);"<<endl
	<<"  m_y2key.Assign(\"y isr\",3,0,info);"<<endl
	<<"  m_y1key.Assign(\"y beam\",3,0,info);"<<endl
	<<"  p_beamchannels = psh->AddIncommingBeamChannel(std::string(\"";
  for (int i=0; i<m_idc.size();i++) {
    if (m_idc[i].find("Z")!=string::npos) {
      chf <<m_idc[i]<<"\"));"<<endl; 
      break;
    }
  }
  chf	<<"  p_isrchannels = psh->AddIncommingISRChannel(std::string(\"";
  for (int i=0; i<m_idc.size();i++) {
    if (m_idc[i].find("Z")!=string::npos) {
      chf <<m_idc[i]<<"\"));"<<endl; 
      break;
    }
  }
  chf	<<"  m_beammode = psh->GetBeamMode();"<<endl
	<<"  m_isrmode  = psh->GetISRMode();"<<endl;
  chf	<<"}"<<endl<<endl;


  chf<<"std::string "<<name<<"::ChID()";
  chf<<endl<<"{"<<endl;  
  chf<<"  return std::string(\""<<m_idstr<<"\");"<<endl;  
  chf<<"}"<<endl;

  chf.close();  

  delete[] flav;
  return rannumber;
}



void Channel_Generator3_NP::Step0(int flag,Point* p,int& rannum,ofstream& sf) 
{
  if (nin != 2) return;
  Point* ph = p->left;
  if (ph->left==0) {
    ph = p->right;
    if (ph->left==0 && p->middle) ph = p->middle;
  }
  if (ph == 0) {
    msg.Error()<<"This seems to be a 2->1 process !!!"<<endl
	       <<"  "<<p->fl<<" -> { "<<p->left->fl<<" "<<p->right->fl<<" }"<<endl;
    abort();
  }

  string m = Order(LinkedMasses(ph));
  
  string help("");
  switch (flag) {
  case -11: 
    char hs[8];
    if (ph->m==0) {
      help+=string("ZS_");
    } 
    else {
      if (!IsZero(ph->fl.Mass())) {
	sprintf(hs,"%i",ph->fl.Kfcode());
	help+=string("ZR")+string(hs)+string("_");
      }
      else help+=string("ZS_");
    }
    sprintf(hs,"%i",(int)PMassSum(ph,0));
    help+=string(hs);
    m_idc.push_back(help);
  case 0: case 1:
    GenerateIS(flag,m,rannum,sf);
    AddToVariables(flag,m+string("_max"),string("m_sppkey[3]"),0,sf);
 
    GenerateMassChain(flag,ph,ph,rannum,sf);
//     help = Order(LinkedMasses(ph->right));
//     if (help.length()>1) {
//       string prt("");
//       for (int i=0;i<m.length();i++) if (help.find(m[i])==-1) prt+=m[i];
//       if (prt.length()>1) {
// 	CalcSmin(flag,"min",prt,sf,0);
// 	AddToVariables(flag,help+string("_max"),string("sqr(sqrt(s") + m +
// 		       string("_max)-sqrt(s") + prt + string("_min))"),0,sf);
//       }
//       else AddToVariables(flag,help+string("_max"),string("sqr(sqrt(s") + m +
// 			  string("_max)-sqrt(ms[") + prt + string("]))"),0,sf);
//     }
//     GenerateMassChain(flag,ph->right,ph->right,rannum,sf);
	  
//     help = Order(LinkedMasses(ph->left));
//     if (help.length()>1) {
//       string prt("");
//       for (int i=0;i<m.length();i++) if (help.find(m[i])==-1) prt+=m[i];
//       if (prt.length()>1) {
// 	if (flag>=0) sf <<"  s"<< prt <<"_min = s"<< prt <<";"<<endl;
// 	AddToVariables(flag,help+string("_max"),string("sqr(sqrt(s") + m +
// 		       string("_max)-sqrt(s") + prt + string("_min))"),0,sf);
//       }
//       else AddToVariables(flag,help+string("_max"),string("sqr(sqrt(s") + m +
// 			  string("_max)-sqrt(ms[") + prt + string("]))"),0,sf);
//     }
//     GenerateMassChain(flag,ph->left,ph->left,rannum,sf);
	  
//     if (tcount==0) {
//       string lm = string("s")+Order(LinkedMasses(ph->left));
//       string rm = string("s")+Order(LinkedMasses(ph->right));
//       if (lm.length()<3) lm = string("ms[")+ LinkedMasses(ph->left)+string("]");
//       if (rm.length()<3) rm = string("ms[")+ LinkedMasses(ph->right)+string("]");
//       AddToVariables(flag,m+string("_min"),
// 		     string("Max(m_sppkey[0],sqr(sqrt(")+lm+
// 		     string(")+sqrt(")+rm+string(")))"),0,sf);
//       GenerateIS(flag,m,rannum,sf);
//     }
    break;
  }
  vector<string> pin0;
  vector<string> pin1;
  GenerateDecayChain(flag,ph,rannum,sf,pin0,pin1);
}

void Channel_Generator3_NP::GenerateDecayChain(int flag,Point* p,int& rannum,ofstream& sf,
					    vector<string> pin0, vector<string> pin1)
{
  if (p->left==0) return;
//    int sa = AntennaS(p);
//    if (sa>2) return QCDAntenna(flag,p,rannum,sf,sa);

  if (p->m==0) {
    int hi = p->fl.Kfcode();
    string tmstr;
    if (flag>=0 && !IsZero(p->fl.Mass())) {
      char hs[3];
      sprintf(hs,"%i",p->number);
      tmstr = string("tmass")+string(hs);
      sf<<"  double "<<tmstr<<" = Flavour(kf::code("<<hi<<")).Mass();"<<endl;
    }
    else tmstr = string("0.");
    string pin0sum(""),pin1sum("");
    if (pin0.size()>0) {
      for (int i=0;i<pin0.size();i++) pin0sum+=pin0[i];
      pin0sum = Order(pin0sum);
      if (pin0sum.length()>1) AddToVariables(flag,string("0_")+pin0sum,string("p[0]-p")+pin0sum,1,sf);
      else AddToVariables(flag,string("0_")+pin0sum,string("p[0]-p[")+pin0sum+string("]"),1,sf);
    }
    if (pin1.size()>0) {
      for (int i=0;i<pin1.size();i++) pin1sum+=pin1[i];
      pin1sum = Order(pin1sum);
      if (pin1sum.length()>1) AddToVariables(flag,string("1_")+pin1sum,string("p[1]-p")+pin1sum,1,sf);
      else AddToVariables(flag,string("1_")+pin1sum,string("p[1]-p[")+pin1sum+string("]"),1,sf);
    }
    pin0sum = string("0_") + pin0sum; 
    pin1sum = string("1_") + pin1sum; 
    string pout0sum = Order(LinkedMasses(p->left));
    string pout1sum = Order(LinkedMasses(p->right));
    string sctmax("m_ctmax");
    if (flag>=0) {
      if (pin0.size()==0 && pout0sum.length()==1 && pin1.size()==0 && pout1sum.length()==1) {
	sf<<"  m_ctmax = Min(cuts->cosmax[0]["<<pout0sum<<"],cuts->cosmax[1]["<<pout1sum<<"]);"<<endl;
      }
      else if (pin0.size()==0 && pin1.size()==0 && pout0sum.length()==1) {
	sf<<"  m_ctmax = cuts->cosmax[0]["<<pout0sum<<"];"<<endl;
      }
      else if (pin0.size()==0 && pin1.size()==0 && pout1sum.length()==1) {
	sf<<"  m_ctmax = cuts->cosmax[1]["<<pout1sum<<"];"<<endl;
      }
      else sctmax = string("1.");
    }
    char hs[4];
    switch (flag) {
    case -11:
      if(!IsZero(p->fl.Mass())) sprintf(hs,"%i",hi);
      else sprintf(hs,"");
      m_idc.push_back(string("TC")+string(hs)+
		      string("_")+pin0sum+string("_")+pin1sum+
		      string("_")+pout0sum+string("_")+pout1sum);
      break;
    case 0:
      sf<<"  CE.TChannelMomenta(";
      if (pin0.size()==0) sf<<"p[0]"; else sf<<"p"<<pin0sum;
      if (pin1.size()==0) sf<<",p[1]"; else sf<<",p"<<pin1sum;
      if (pout0sum.length()==1) sf<<",p["<<pout0sum<<"]"; else sf<<",p"<<pout0sum;
      if (pout1sum.length()==1) sf<<",p["<<pout1sum<<"]"; else sf<<",p"<<pout1sum;
      sf<<",s"<<pout0sum<<",s"<<pout1sum;
      sf<<","<<tmstr<<",m_alpha,"<<sctmax<<",m_ctmin,m_amct,0,ran["<<rannum++<<"],ran[";
      sf<<rannum++<<"]);"<<endl;
      break;
    default:
      if(!IsZero(p->fl.Mass())) sprintf(hs,"%i",hi);
      else sprintf(hs,"");
      string idh = string("TC")+string(hs)+
	string("_")+pin0sum+string("_")+pin1sum+string("_")+pout0sum+string("_")+pout1sum;
      sf<<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT)"<<endl; 
      sf<<"    m_k"<<idh<<"<<CE.TChannelWeight(";
      if (pin0.size()==0) sf<<"p[0]"; else sf<<"p"<<pin0sum;
      if (pin1.size()==0) sf<<",p[1]"; else sf<<",p"<<pin1sum;
      if (pout0sum.length()==1) sf<<",p["<<pout0sum<<"]"; else sf<<",p"<<pout0sum;
      if (pout1sum.length()==1) sf<<",p["<<pout1sum<<"]"; else sf<<",p"<<pout1sum;
      sf<<","<<tmstr<<",m_alpha,"<<sctmax<<",m_ctmin,m_amct,0);"<<endl;
      sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
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
    if (l->left==0) moml = string("p[") + lm + string("]");
    else moml = string("p") + Order(lm);
    if (r->left==0) momr = string("p[") + rm + string("]");
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
	    sf<<"    m_k"<<idh<<"<<CE.Isotropic2Weight("<<moml<<","<<momr<<");"<<endl;
	    sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
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
  	  sf<<"    m_k"<<idh<<"<<CE.Isotropic2Weight("<<momr<<","<<moml<<");"<<endl;
  	  sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
	  //	  	  sf<<"  wt *= CE.Isotropic2Weight("<<momr<<","<<moml<<");"<<endl;
	}
      }
      }
    }
    if (flag==0) rannum += 2;
  }
  
  GenerateDecayChain(flag,p->left,rannum,sf,pin0,pin1);
  GenerateDecayChain(flag,p->right,rannum,sf,pin0,pin1);
}


bool Channel_Generator3_NP::QCDAntenna(int flag,Point* p,int& rannum,ofstream& sf,int n)
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
    sf <<"  double s0"<<acount<<" = cuts->scut["<<mummy[0]<<"]["<<mummy[1]<<"];"<<endl;
    sf <<"  Vec4D ps"<<acount<<"["<<n<<"];"<<endl;
    sf <<"  CE.QCDAPMomenta(ps"<<acount<<",p"<<mummy<<","<<n<<",s0"<<acount<<");"<<endl;
    for (int i=0;i<n;i++) {
      sf<<"  p["<<mummy[i]<<"] = ps"<<acount<<"["<<i<<"];"<<endl;
    }
    break;
  default:
    string idh = string("AP_")+mummy;
    sf <<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT) {"<<endl; 
    sf <<"    double s0"<<acount<<" = cuts->scut["<<mummy[0]<<"]["<<mummy[1]<<"];"<<endl;
    sf <<"    Vec4D ps"<<acount<<"["<<n<<"];"<<endl;
    for (int i=0;i<n;i++) {
      sf<<"    ps"<<acount<<"["<<i<<"] = p["<<mummy[i]<<"];"<<endl;
    }
    sf <<"    m_k"<<idh<<"<<CE.QCDAPWeight(ps"<<acount<<","<<n<<",s0"<<acount<<");"<<endl;    
    sf <<"  }"<<endl;
    sf<<"  wt *= m_k"<<idh<<".Weight();"<<endl<<endl;
  }
  acount++;
  return 1;
}





void Channel_Generator3_NP::GenerateMassChain(int flag,Point* p,Point* clmp,int& rannum,ofstream& sf)
{
  if (p->left==0) {
    string m = LinkedMasses(p);
    AddToVariables(flag,m,string("ms[")+m+string("]"),0,sf);
    return;
  }
  string lm,rm;
  string mummy,prt;
  string clm = Order(LinkedMasses(clmp));
  lm = Order(LinkedMasses(p->left));
  rm = Order(LinkedMasses(p->right));

  mummy = Order(lm+rm);
  for (int i=0;i<clm.length();i++) {
    if (mummy.find(clm[i])==-1) prt+=clm[i];
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
		   string("_max)-sqrt(ms[") + prt + string("]))"),0,sf);
  }
  
  if (rm.length()>1 && lm.length()>1) sclmp = p;
  else if (clmp->left==p && LinkedMasses(clmp->right).length()>1) sclmp = p;
  else if (clmp->right==p && LinkedMasses(clmp->left).length()>1) sclmp = p;

  if (p->m!=2) {
    //cout<<":GenerateMassChain: right: "<<LinkedMasses(p->right)<<endl;
    GenerateMassChain(flag,p->right,sclmp,rannum,sf);
    //cout<<":GenerateMassChain: left:  "<<LinkedMasses(p->left)<<endl;
    GenerateMassChain(flag,p->left,sclmp,rannum,sf);
  }

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
  CalcSmin(flag,"min",mummy,sf,0);
  if (p->m!=2) {
    if (mummy.length()>2 && flag>=0) {
      sf <<"  s"<< mummy <<"_min = Max(s"<< mummy <<"_min,sqr(sqrt(s"<< lm <<")+sqrt(s"<< rm <<")));"<<endl;
    }
  }

  double maxpole = -1.;
  double res = ATOOLS::sqr(p->fl.Width()*p->fl.Mass());
  if (p->m>0 && !ATOOLS::IsZero(res) && Massive(p->fl)) maxpole = 1./res;

  int hi = 4;
  if (mummy.length()>2) hi = 2;
  if (maxpole>0.) {
    hi = (p->fl).Kfcode();
    if (flag>=0) sf<<"  Flavour fl"<<mummy<<" = "<<"Flavour(kf::code("<<hi<<"));"<<endl;
  } 
  string thexp("1.5");
  if (p->m==0) thexp = string("1.5");
  switch (flag) {
  case -11:
    if (maxpole>0.) {
      char hs[4];
      sprintf(hs,"%i",hi);
      m_idc.push_back(string("MP")+string(hs)+string("_")+mummy);
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
      sf<<"  double s"<<mummy<<" = CE.ThresholdMomenta("<<thexp<<","
	<<hi<<".*sqrt(s"<<mummy<<"_min),s"<<mummy<<"_min,"
	<<"s"<<mummy<<"_max,ran["<<rannum<<"]);"<<endl;
    }
    AddToVariables(flag,mummy,string(""),0,sf);
    rannum++;
    break;
  default:
    string s(""); 
    for (short int i=0;i<mummy.length()-1;i++) s += string("p[")+mummy[i]+string("]+");
    s += string("p[")+mummy[mummy.length()-1]+string("]");
    
    AddToVariables(flag,mummy,s,1,sf);
    AddToVariables(flag,mummy,string("dabs(p")+mummy+string(".Abs2())"),0,sf);
    if (maxpole>0.) {
      sf<<"  wt *= CE.MassivePropWeight(fl"<<mummy<<".Mass(),"<<"fl"<<mummy<<".Width(),1,"
	<<"s"<<mummy<<"_min,s"<<mummy<<"_max,"<<"s"<<mummy<<");"<<endl;
    }
    else {
      sf<<"  wt *= CE.ThresholdWeight("<<thexp<<","<<hi<<".*sqrt(s"<<mummy<<"_min),s"<<mummy<<"_min,"
	<<"s"<<mummy<<"_max,s"<<mummy<<");"<<endl;
    }
  }

  if (p->m==2) {
    if (rm.length()>1) {
      if (lm.length()>1) {
	CalcSmin(flag,"min",lm,sf,0);
	AddToVariables(flag,rm + string("_max"),string("sqr(sqrt(s") + mummy +
		       string(")-sqrt(s") + lm + string("_min))"),0,sf);    
      }
      else {
	AddToVariables(flag,rm + string("_max"),string("sqr(sqrt(s") + mummy +
		       string(")-sqrt(ms[") + lm + string("]))"),0,sf);
      }
    }
    GenerateMassFwd(flag,p->right,rannum,sf);
    if (lm.length()>1) {
      AddToVariables(flag,lm + string("_max"),string("sqr(sqrt(s") + mummy +
		     string(")-sqrt(s") + rm + string("))"),0,sf);    
    }
    GenerateMassFwd(flag,p->left,rannum,sf);
  }
}

void Channel_Generator3_NP::GenerateMassFwd(int flag,Point* p,int& rannum,ofstream& sf)
{
  if (p->left==0) {
    string m = LinkedMasses(p);
    AddToVariables(flag,m,string("ms[")+m+string("]"),0,sf);
    return;
  }
  string lm,rm;
  string mummy;
  lm = Order(LinkedMasses(p->left));
  rm = Order(LinkedMasses(p->right));
  mummy = Order(lm+rm);
  CalcSmin(flag,"min",mummy,sf,0);

  double maxpole = -1.;
  double res = ATOOLS::sqr(p->fl.Width()*p->fl.Mass());
  if (p->m>0 && !ATOOLS::IsZero(res) && Massive(p->fl)) maxpole = 1./res;

  int hi = 4;
  if (mummy.length()>2) hi = 2;
  if (maxpole>0.) {
    hi = (p->fl).Kfcode();
    if (flag>=0) sf<<"  Flavour fl"<<mummy<<" = "<<"Flavour(kf::code("<<hi<<"));"<<endl;
  } 
  switch (flag) {
  case -11: 
    if (maxpole>0.) {
      char hs[4];
      sprintf(hs,"%i",hi);
      m_idc.push_back(string("MP")+string(hs)+string("_")+mummy);
    }
    else m_idc.push_back(string("MlP_")+mummy);
    break;
  case 0: 
    sf<<"  Vec4D  p"<<mummy<<";"<<endl;
    if (maxpole>0.) {
      sf<<"  double s"<< mummy
	  <<" = CE.MassivePropMomenta(fl"<<mummy<<".Mass(),"<<"fl"<<mummy<<".Width(),1,"
	<<"s"<<mummy<<"_min,s"<<mummy<<"_max,ran["<<rannum<<"]);"<<endl;
    }
    else {
      sf<<"  double s"<<mummy<<" = CE.MasslessPropMomenta(1.,s"<<mummy<<"_min,"
	<<"s"<<mummy<<"_max,ran["<<rannum<<"]);"<<endl;
    }
    rannum++;
    break;
  default:
    string s(""); 
    for (short int i=0;i<mummy.length()-1;i++) s += string("p[")+mummy[i]+string("]+");
    s += string("p[")+mummy[mummy.length()-1]+string("]");
    
    AddToVariables(flag,mummy,s,1,sf);
    AddToVariables(flag,mummy,string("dabs(p")+mummy+string(".Abs2())"),0,sf);
    if (maxpole>0.) {
	sf<<"  wt *= CE.MassivePropWeight(fl"<<mummy<<".Mass(),"<<"fl"<<mummy<<".Width(),1,"
	  <<"s"<<mummy<<"_min,s"<<mummy<<"_max,"<<"s"<<mummy<<");"<<endl;
    }
    else {
      sf<<"  wt *= CE.MasslessPropWeight(1.,s"<<mummy<<"_min,"
	<<"s"<<mummy<<"_max,s"<<mummy<<");"<<endl;
    }
  }

  if (rm.length()>1) {
    if (lm.length()>1) {
      CalcSmin(flag,"min",lm,sf,0);
      AddToVariables(flag,rm + string("_max"),string("sqr(sqrt(s") + mummy +
		     string(")-sqrt(s") + lm + string("_min))"),0,sf);    
    }
    else {
      AddToVariables(flag,rm + string("_max"),string("sqr(sqrt(s") + mummy +
		     string(")-sqrt(ms[") + lm + string("]))"),0,sf);
    }
  }
  GenerateMassFwd(flag,p->right,rannum,sf);
  if (lm.length()>1) {
    AddToVariables(flag,lm + string("_max"),string("sqr(sqrt(s") + mummy +
		   string(")-sqrt(s") + rm + string("))"),0,sf);    
  }
  GenerateMassFwd(flag,p->left,rannum,sf);

}

void Channel_Generator3_NP::GenerateIS(int flag,string& m,int& rannum,ofstream& sf)
{
  switch (flag) {
  case 0:  
    sf <<"  if (p_beamchannels) {"<<endl
       <<"    p_beamchannels->GeneratePoint(m_spkey,m_y1key,m_beammode);"<<endl
       <<"    m_sppkey[2] = m_spkey[3];"<<endl
       <<"    m_sppkey[1] = m_sppkey[2]*m_sppkey[4];"<<endl
       <<"  }"<<endl
       <<"  if (p_isrchannels) p_isrchannels->GeneratePoint(m_sppkey,m_y2key,m_isrmode);"<<endl
       <<"  else m_sppkey[3] = m_spkey[3];"<<endl
       <<"  if (p_beamchannels || p_isrchannels) {"<<endl
       <<"    cuts->Update(m_y1key[2]+m_y2key[2]);"<<endl
       <<"    double E1=0.5*sqrt(m_sppkey[3])*(1.+(sqr(ms[0])-sqr(ms[1]))/m_sppkey[3]);"<<endl
       <<"    p[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-sqr(ms[0])));"<<endl
       <<"    p[1]=Vec4D(sqrt(m_sppkey[3])-E1,(-1.)*Vec3D(p[0]));"<<endl
       <<"  }"<<endl
       <<"  Vec4D p"<<m<<"=p[0]+p[1];"<<endl;
    break;
  case 1: 
    string idh=m_idc.back();
    sf <<"  if (m_k"<<idh<<".Weight()==ATOOLS::UNDEFINED_WEIGHT) {"<<endl 
       <<"    double iswt = 1.;"<<endl
       <<"    if (p_beamchannels) {"<<endl
       <<"      p_beamchannels->GenerateWeight(m_beammode);"<<endl
       <<"      iswt /= p_beamchannels->Weight();"<<endl
       <<"    }"<<endl
       <<"    if (p_isrchannels) {"<<endl
       <<"      p_isrchannels->GenerateWeight(m_isrmode);"<<endl
       <<"      iswt /= p_isrchannels->Weight();"<<endl
       <<"    }"<<endl
       <<"    m_k"<<idh<<"<<iswt;"<<endl
       <<"  }"<<endl
       <<"  wt *= m_k"<<idh<<".Weight();"<<endl;
  }
}

void Channel_Generator3_NP::SetProps(Point* p,Point** props,Point** propt, int& count)
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

void Channel_Generator3_NP::CalcTSmin(int flag,vector<string>& p,ofstream& sf)
{
  string help;
  for (short int i=0;i<p.size();i++) {
    if (p[i].length()==1) help += p[i];
  }

  string psum;
  for (short int i=0;i<p.size();i++) psum += p[i];


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


  for (short int i=0;i<p.size();i++) {
    if (p[i].length()>1) s += string("+sqrt(s") + Order(p[i]) + string(")");
  }
  s += string(")");

  if (!CheckVariables(flag,Order(psum)+string("_min"),0)) {
    AddToVariables(flag,Order(psum) + string("_min"),s,0,sf);
  }
  else sf<<"  s"<<Order(psum)<<"_min = Max(s"<<Order(psum)<<"_min,"<<s<<");"<<endl;
}



void Channel_Generator3_NP::CalcSmin(int flag,char* min,string lm,ofstream& sf,Point* p)
{
  // sum_{i<j} sij + (2-n)*sum_i m_i^2

  if (lm.length()>1) {
    AddToVariables(flag,Order(lm) + string("_") + string(min),
		   string("cuts->Getscut(std::string(\"") + Order(lm) + string("\"))"),0,sf);
  }
  else {
    AddToVariables(flag,Order(lm) + string("_") + string(min),
		   string("ms[") + Order(lm) + string("]"),0,sf);
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

string Channel_Generator3_NP::EstimatePeak(int flag,Point* p,ofstream& sf)
{
  string help,h1;
  if (p->m==0) {
    help = EstimatePeak(flag,p->left,sf);
    h1   = EstimatePeak(flag,p->right,sf);
    if (help.length()>0 && h1.length()>0) help+=string("+");
    help += h1;
  }
  else {
    string h1 = Order(LinkedMasses(p));
    if (h1.length()==1) {}//help = string("sqrt(ms[")+h1+string("])");
    else {
      //if (!CheckVariables(flag,h1+string("_min"),0)) CalcSmin(flag,"min",h1,sf,0);
      //help = string("4.*sqrt(s")+h1+string("_min)");
      if (PMassSum(p,0)>5.) help= GetFlMass(p);
    }
  }
  return help;
}

string Channel_Generator3_NP::GetFlMass(Point* p)
{
  if (p->left==0) return string("");
  if (p->fl.Mass()>PMassSum(p->left,0)+PMassSum(p->right,0)) {
    char hs[5];
    sprintf(hs,"%i",(p->fl).Kfcode());
    return string("Flavour(kf::code(")+string(hs)+string(")).Mass()");
  } 
  string h1 = GetFlMass(p->left);
  string h2 = GetFlMass(p->right);
  if (h1.length()==0) return h2;
  if (h2.length()==0) return h1;
  return h1+string("+")+h2;
}

string Channel_Generator3_NP::LinkedMasses(Point* p)
{
  if (p->left==0) {
    char help[4];
    sprintf(help,"%i",p->number);
    return string(help);
  }
  return LinkedMasses(p->left)+LinkedMasses(p->right);
}


void Channel_Generator3_NP::IdentifyProps(Point* _plist)
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

void Channel_Generator3_NP::BackLinks(Point* p,Point* &endp)
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

void Channel_Generator3_NP::InitT(Point* p)
{
  p->t = 0;
  if (p->left==0) return;
  InitT(p->left);
  InitT(p->right);
}

string Channel_Generator3_NP::IString(int i)
{
  MyStrStream sstr;
  sstr<<i;
  string istr;
  sstr>>istr;
  return istr;
}

string Channel_Generator3_NP::Order(string s)
{
  int beg = s.find("_");
  if (beg!=-1) {
    return Order(s.substr(0,beg)) + string("_") + Order(s.substr(beg+1));
  }
  if (s[0]>='9' || s[0]<='0') return s;

  for (short i=0;i<s.length();i++) 
    for (short j=i+1;j<s.length();j++) {
      if (s[i]>s[j]) {
	char help = s[i];
	s[i] = s[j];
	s[j] = help;
      } 
    }
  return s;
}

void  Channel_Generator3_NP::AddToVariables(int flag,const string& lhs,const string& rhs,const int& type,
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
      msg.Error()<<" ERROR in Channel_Generator3_NP::AddToVariables ()"<<endl;
      abort();
    }
  }
}

bool  Channel_Generator3_NP::CheckVariables(int flag,const string& lhs,const int& type)
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

int Channel_Generator3_NP::AntennaS(Point* p)
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

std::string Channel_Generator3_NP::CreateChannelID(int echflag)
{
  m_idc.clear();
  extrachannelflag = echflag;
  int    rannum = 1;
  int   maxnumb = 0;
  ofstream chf;
  Step0(-11,m_topos[echflag],rannum,chf);

  string idstr;
  std::sort(m_idc.begin(),m_idc.end(),Compare_String());
  for (String_List::iterator it = m_idc.begin();it!=m_idc.end();++it) {
    idstr+=(*it);
    idstr+=string("$");
  }
  if (echflag==0) idstr = string("CGPR$")+idstr;
  if (echflag==1) idstr = string("CGPL$")+idstr;
  m_idstr = idstr;
  return idstr;
}

double Channel_Generator3_NP::PMassSum(Point* p,vector<int>* kfs)
{
  if (!p->left) return 0.;
  double m = 0.;
  if (p->m>0 && p->fl.IsMassive()) {
    m = p->fl.Mass();
    if (kfs) kfs->push_back(p->fl.Kfcode());
  }
  double mc = PMassSum(p->left,kfs) + PMassSum(p->right,kfs);
  return Max(m,mc);  
}
