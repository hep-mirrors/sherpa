#include "Channel_Generator.H"
#include "Topology.H"
#include "Message.H"
#include "MathTools.H"
#include "MyStrStream.H"

#include <stdio.h>


using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Channel_Generator::Channel_Generator(int _nin,int _nout,Flavour * _fl,
                                     Point * _plist) 
  : nin(_nin), nout(_nout)
{
  Topology top;
  plist  = new Point[2*(nout+1)];
  int ll = 0;
  top.Copy(_plist,plist,ll);
  IdentifyProps(plist);
}

Channel_Generator::~Channel_Generator() { delete[] plist; }

int Channel_Generator::MakeChannel(int& echflag,int n,string& path,string& pID)
{  
  int oldflag      = echflag;
  extrachannelflag = echflag;

  //add Channel
  char name[22];
  sprintf(name,"C%i",n);

  if (echflag!=0) sprintf(name,"%s%c",name,'a'+extrachannelflag-1);
  
  string filename = string("Process/")+path+string("/")+
                    string(name)+string(".C");
  
  ifstream from;
  from.open((filename).c_str());
  
  int    rannum = 1;

  ofstream chf;
  chf.open((filename).c_str());

  chf<<"#include "<<'"'<<"P.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"Run_Parameter.H"<<'"'<<endl;
  chf<<"#include "<<'"'<<"Channel_Elements.H"<<'"'<<endl<<endl;  
  chf<<"using namespace PHASIC;"<<endl;  
  chf<<"using namespace ATOOLS;"<<endl<<endl;


  //Momenta
  chf<<"void "<<pID.c_str()<<"::"<<name;
  chf<<"Momenta(Vec4D * p,Cut_Data * cuts,double * ran)"<<endl;
  chf<<"{"<<endl;
  
  Flavour * flav    = new Flavour[nout];  
  int       maxnumb = 0;

  newchannel = 0;
  Step0(0,plist,rannum,chf,flav,maxnumb);
  ClearDeclarations(); 
  extrachannelflag = newchannel;
  chf<<"}"<<endl<<endl;

  int rannumber = rannum;

  echflag          = extrachannelflag;
  extrachannelflag = oldflag;

  //Weight
  chf<<"double "<<pID.c_str()<<"::"<<name;
  chf<<"Weight(Vec4D* p,Cut_Data * cuts)"<<endl<<"{"<<endl;
  chf<<"  double wt = 1.;"<<endl;

  maxnumb = 0;

  Step0(1,plist,rannum,chf,flav,maxnumb);
  ClearDeclarations();
  chf<<"  if (wt!=0.) wt = 1./wt/pow(2.*M_PI,"<<nout<<"*3.-4.);"<<endl;
  chf<<endl<<"  return wt;"<<endl; 
  chf<<"}"<<endl<<endl;
  
  
  chf<<"int "<<pID.c_str()<<"::"<<name<<"Resonances(Flavour*& res_fl)"<<endl<<"{"<<endl;  
  chf<<"  res_fl = new Flavour["<<maxnumb<<"];"<<endl;
  for (short int i=0;i<maxnumb;i++) {    
    int hi = flav[i].Kfcode();
    chf<<"  res_fl["<<i<<"] = Flavour(kf::code("<<hi<<"));"<<endl;
  }
  chf<<"  return "<<maxnumb<<";"<<endl;
  chf<<"}"<<endl<<endl;
  
  chf<<"void "<<pID.c_str()<<"::"<<name<<"ISRtype(int & type,double & mass,double & width)";
  chf<<endl<<"{"<<endl;  
  Step0(2,plist,rannum,chf,flav,maxnumb);
  chf<<"}"<<endl;

  chf.close();
  delete[] flav;
  return rannumber;
}


void Channel_Generator::Step0(int flag,Point* p,int& rannum,ofstream& sf,
			      Flavour* flav,int& maxnumb) 
{
  string m;
  if (nin == 2) {
    switch (tcount) {
      case 0 : {
	if (flag<2) {
	  m = LinkedMasses(p->left);
	  if (m.length()<2) m = LinkedMasses(p->right);
	  AddToVariables(flag,m,string("p[0] + p[1]"),1,sf);
	  AddToVariables(flag,m,string("dabs(p")+Order(m)+string(".Abs2())"),0,sf);
    
	  flag += 10;
	  if (!StepS(flag,p->left,rannum,sf,flav,maxnumb)) {
	    if (!StepS(flag,p->right,rannum,sf,flav,maxnumb)) {
	      if (p->middle == 0) {
		msg.Error()<<"This seems to be a 2->1 process !!!"<<endl
			   <<"  "<<p->fl<<" -> { "<<p->left->fl<<" "<<p->right->fl<<" }"<<endl;
		abort();
	      }
	    }
	  }
	  break;
	}
	if (flag==2) {
	  sf<<"  type  = 1;"<<endl
	    <<"  mass  = Flavour(kf::code("<<p->left->fl.Kfcode()<<")).Mass();"<<endl
	    <<"  width = Flavour(kf::code("<<p->left->fl.Kfcode()<<")).Width();"<<endl;
	  return;
	}
      }
      default : {
	if (flag<=2) {
	  StepNT(flag,tcount,p,rannum,sf,flav,maxnumb);
	  break;
	}
      }
    }
    return;
  }
}

bool Channel_Generator::StepS(int flag,Point* p,int& rannum,
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
  if (l->left==0) moml = string("p[") + lm + string("]");
            else  moml = string("p") + Order(lm);
  if (r->left==0) momr = string("p[") + rm + string("]");
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
  if (flag>9) { first = 1; flag -= 10; }

  // Check for decay type.
  if ((!first) && (l->left==0) && (l->fl.IsVector()) && 
      (!(l->fl.IsMassive())) && (r->fl.IsFermion()) ) {
    if (flag==0) {
      sf<<"  CE.Anisotropic2Momenta(p"<<Order(mummy)
	<<",s"<<Order(lm)<<",s"<<Order(rm)<<",1.,-1.,1."<<","<<moml<<","<<momr<<","
	<<"ran["<<rannum<<"],ran["<<rannum+1<<"]);"<<endl;
    }
    else {
      sf<<"  wt *= CE.Anisotropic2Weight(1,-1.,1.,"<<moml<<",";
      sf<<momr<<");"<<endl;
    }
  }
  else {
    if ((!first) && (r->left==0) && (r->fl.IsVector()) && 
	(!(r->fl.IsMassive())) && (l->fl.IsFermion()) ) {
      //anisotropic decay for left outgoing massless vectors
      if (flag==0) {
	sf<<"  CE.Anisotropic2Momenta(p"<<Order(mummy)
	  <<",s"<<Order(rm)<<",s"<<Order(lm)<<",1.,-1.,1."<<","<<momr<<","<<moml<<","
	  <<"ran["<<rannum<<"],"<<"ran["<<rannum+1<<"]);"<<endl;
      }
      else {
	sf<<"  wt *= CE.Anisotropic2Weight(1.,-1.,1.,"<<momr<<","<<moml<<");"<<endl;
      }
    }
    else {
      if ((l->number) < (r->number)) {
	if (flag==0) {
	  sf<<"  CE.Isotropic2Momenta(p"<<mummy<<",s"<<Order(lm)<<",s"<<Order(rm)
	    <<","<<moml<<","<<momr<<",ran["<<rannum<<"],"<<"ran["<<rannum+1<<"]);"<<endl;
	}
	else {
	  sf<<"  wt *= CE.Isotropic2Weight("<<moml<<","<<momr<<");"<<endl;
	}
      }
      else {
	if (flag==0) {
	  sf<<"  CE.Isotropic2Momenta(p"<<mummy<<",s"<<Order(rm)<<",s"<<Order(lm)
	    <<","<<momr<<","<<moml<<",ran["<<rannum<<"],"<<"ran["<<rannum+1<<"]);"<<endl;
	}
	else {
	  sf<<"  wt *= CE.Isotropic2Weight("<<momr<<","<<moml<<");"<<endl;
	}
      }
    }
  }
  if (flag==0) rannum += 2;
  
  StepS(flag,l,rannum,sf,flav,maxnumb);
  StepS(flag,r,rannum,sf,flav,maxnumb);

  return 1;
}


void Channel_Generator::StepNT(int flag,int tcount,Point* p,int& rannum,ofstream& sf,
			       Flavour* flav,int& maxnumb)
{
  if (p->left==0) return;

  Point **props = new Point*[tcount+1];
  Point **propt = new Point*[tcount+1];

  int counts = 0;

  SetProps(p,props,propt,counts);

  if (flag==2) {
    sf<<"  type  = 2;"<<endl
      <<"  mass  = ";
    for (int i=0;i<tcount;i++) {
      sf<<"Flavour(kf::code("<<props[i]->fl.Kfcode()<<")).Mass() + ";
    } 
    sf<<"Flavour(kf::code("<<props[tcount]->fl.Kfcode()<<")).Mass();"<<endl;
    sf<<"  width = 0.;"<<endl;
    return;
  }

  string* s    = new string[tcount+1];
  string  m;  

  for (short int i=0;i<tcount+1;i++) {
    s[i] = LinkedMasses(props[i]);
    m += s[i];
    flav[maxnumb] = props[i]->fl;maxnumb++;
  }

  AddToVariables(flag,m,string("p[0] + p[1]"),1,sf);
  AddToVariables(flag,m,string("dabs(p")+Order(m)+string(".Abs2())"),0,sf);

  sf<<"  double amct  = 1.;"<<endl;
  sf<<"  double alpha = 0.5;"<<endl;
  sf<<"  double ctmax = 0.;"<<endl;
  sf<<"  double ctmin = 2.;"<<endl;

  GenerateMasses(flag,props,tcount+1,rannum,sf);

  vector<string> pin0;pin0.push_back(string("0"));
  vector<string> pin1;pin1.push_back(string("1"));

  vector<string> pout0;
  vector<string> pout1;

    
  int count = 0;

  if (props[tcount]->number<99 && props[tcount]->fl.IsVector()) {
    for (short int i=0;i<tcount;i++) pout0.push_back(s[i]);
    pout1.push_back(s[tcount]);
    count = tcount - 1;
  }
  else {    
    pout0.push_back(s[0]);
    for (short int i=1;i<tcount+1;i++) pout1.push_back(s[i]);
  }
  
  AddToVariables(flag,string("0_"),string("p[0]"),1,sf);
  AddToVariables(flag,string("1_"),string("p[1]"),1,sf);

  SingleTStep(flag,s,propt,tcount,rannum,sf,count,pin0,pin1,pout0,pout1);

  for (short int i=0;i<tcount+1;i++) StepS(flag,props[i],rannum,sf,flav,maxnumb);

  delete[] props;
  delete[] propt;
  delete[] s;
}

void Channel_Generator::SingleTStep(int flag,string* s,Point** propt,int tcount,
				    int& rannum,ofstream& sf,int count,
				    vector<string> pin0, vector<string> pin1,
				    vector<string> pout0,vector<string> pout1)
{ 
  int hi = propt[count]->fl.Kfcode(); 
  sf<<"  double tmass"<<count<<" = Flavour(kf::code("<<hi<<")).Mass();"<<endl;

  string pout0sum;
  for (short int i=0;i<pout0.size();i++) pout0sum += pout0[i];
  string pout1sum;
  for (short int i=0;i<pout1.size();i++) pout1sum += pout1[i];
  string pin0sum = string("0_");
  for (short int i=1;i<pin0.size();i++) pin0sum += pin0[i]; 
  string pin1sum = string("1_");
  for (short int i=1;i<pin1.size();i++) pin1sum += pin1[i];

  //mins
  if (pout0.size()>1) CalcTSmin(flag,pout0,sf); 
  if (pout1.size()>1) CalcTSmin(flag,pout1,sf); 
  //maxs
  if (pout0.size()>1) { 
    string s = string("sqr(sqrt(dabs((p")+Order(pin0sum)+string("+p")+
               Order(pin1sum)+string(").Abs2()))-");
    if (pout1.size()==1) {
       if (pout1[0].length()==1) s += string("sqrt(ms[")+pout1[0]+string("])");
                           else s += string("sqrt(s")+Order(pout1[0])+string(")");    
    }
    else s += string("sqrt(s")+Order(pout1sum)+string("_min)");  
    s += ")";

    AddToVariables(flag,Order(pout0sum) + string("_max"),s,0,sf);
    //sf<<"  double s"<<Order(pout0sum)<<"_max = "<<s<<";"<<endl;
  }
  //first dice
  if (pout0.size()>1) {
    if (flag==0) {
      sf<<"  Vec4D  p"<<Order(pout0sum)<<";"<<endl;
      sf<<"  double s"<<Order(pout0sum);

      if (pout1.size()==1 && pin1.size()==1 && pout1[0].length()==1 && extrachannelflag==0) {
	//check for extra LL-Channel
        MyStrStream sstr;
	sstr<<pout1[0];
	int i;
	sstr>>i;
      }


      if (pout1.size()==1 && pin1.size()==1 && pout1[0].length()==1 && extrachannelflag==1) {
	sf<<" = CE.LLPropMomenta(0.99,1.001*sqr(rpa.gen.Ecms()),s"<<Order(pout0sum)<<"_min,"
	  <<"s"<<Order(pout0sum)<<"_max,ran["<<rannum<<"]);"<<endl;
	extrachannelflag = 0;
      }
      else {
	sf<<" = CE.MasslessPropMomenta(0.5,s"<<Order(pout0sum)<<"_min,"
	  <<"s"<<Order(pout0sum)<<"_max,ran["<<rannum<<"]);"<<endl;
      }
      rannum++;
    }
    else {
      string s;
      for (short int i=0;i<pout0sum.length()-1;i++)
	s += string("p[")+pout0sum[i]+string("]+");
      s += string("p[")+pout0sum[pout0sum.length()-1]+string("]");
    
      AddToVariables(flag,pout0sum,s,1,sf);
      AddToVariables(flag,pout0sum,string("dabs(p")+Order(pout0sum)+string(".Abs2())"),0,sf);

      //sf<<"  Vec4D p"<<Order(pout0sum)<<" = "<<s<<";"<<endl
      //<<"  double s"<<Order(pout0sum)<<" = dabs(p"<<Order(pout0sum)<<".Abs2());"<<endl;    

      if (pout1.size()==1 && pin1.size()==1 && pout1[0].length()==1 && extrachannelflag==1) {
	sf<<"  wt *= CE.LLPropWeight(0.99,1.001*sqr(rpa.gen.Ecms()),s"<<Order(pout0sum)<<"_min,"
	  <<"s"<<Order(pout0sum)<<"_max,s"<<Order(pout0sum)<<");"<<endl;
      }
      else {
	sf<<"  wt *= CE.MasslessPropWeight(0.5,s"<<Order(pout0sum)<<"_min,"
	  <<"s"<<Order(pout0sum)<<"_max,s"<<Order(pout0sum)<<");"<<endl;
      }
    }
  } 
  //second dice
  if (pout1.size()>1) {
    string s = string("sqr(sqrt(dabs((p")+Order(pin0sum)+string("+p")+
               Order(pin1sum)+string(").Abs2()))-");
    if (pout0.size()==1) {
      if (pout0[0].size()==1) s += string("sqrt(ms[")+pout0[0]+string("])");
                         else s += string("sqrt(s")+Order(pout0[0])+string(")");    
    }
    else s += string("sqrt(s")+Order(pout1sum)+string(")");  
    s += string(")");
    
    AddToVariables(flag,pout1sum + string("_max"),s,0,sf);
    //sf<<"  double s"<<Order(pout1sum)<<"_max = "<<s<<";"<<endl;

    if (flag==0) {
      sf<<"  Vec4D  p"<<Order(pout1sum)<<";"<<endl
	<<"  double s"<<Order(pout1sum)
	<<" = CE.MasslessPropMomenta(0.5,s"<<Order(pout1sum)<<"_min,"
	<<"s"<<Order(pout1sum)<<"_max,ran["<<rannum<<"]);"<<endl;
      rannum++;
    }
    else {
      string s;
      for (short int i=0;i<pout1sum.length()-1;i++)
	s += string("p[")+pout1sum[i]+string("]+");
      s += string("p[")+pout1sum[pout1sum.length()-1]+string("]");
      
      AddToVariables(flag,pout1sum,s,1,sf);
      AddToVariables(flag,pout1sum,string("dabs(p")+Order(pout1sum)+string(".Abs2())"),0,sf);

      sf<<"  wt *= CE.MasslessPropWeight(0.5,s"<<Order(pout1sum)<<"_min,"
	<<"s"<<Order(pout1sum)<<"_max,s"<<Order(pout1sum)<<");"<<endl;
    }
  }


  if (flag==0) {  
    sf<<"  CE.TChannelMomenta(p"<<Order(pin0sum)<<",p"<<Order(pin1sum);
    if (pout0.size()==1 && pout0[0].length()==1) sf<<",p["<<pout0[0]<<"]";
                                            else sf<<",p"<<Order(pout0sum);
    if (pout1.size()==1 && pout1[0].length()==1) sf<<",p["<<pout1[0]<<"]";
                                            else sf<<",p"<<Order(pout1sum);
    sf<<",s"<<Order(pout0sum)<<",s"<<Order(pout1sum);
    sf<<",tmass"<<count<<",alpha,ctmax,ctmin,amct,0,ran["<<rannum++<<"],ran[";
    sf<<rannum++<<"]);"<<endl;
  }
  else {
    sf<<"  wt *= CE.TChannelWeight(p"<<Order(pin0sum)<<",p"<<Order(pin1sum);
    if (pout0.size()==1 && pout0[0].length()==1) sf<<",p["<<pout0[0]<<"]";
                                              else sf<<",p"<<Order(pout0sum);
    if (pout1.size()==1 && pout1[0].length()==1) sf<<",p["<<pout1[0]<<"]";
                                              else sf<<",p"<<Order(pout1sum);
    sf<<",tmass"<<count<<",alpha,ctmax,ctmin,amct,0);"<<endl;
  }

  //dice new

  if (pout0.size()>1 || pout1.size()>1) {
    if (pout0.size()>1) {
      //pout1  -> pin1
      if (pout1.size()==1) {
	pin1.push_back(pout1[0]);
	string pin1sum = string("1_");
	for (short int i=1;i<pin1.size();i++) pin1sum += pin1[i];
	string s = string("p[1]");
	for (short int i=1;i<pin1.size();i++) {
	  if (pin1[i].length()==1) s += string("-p[")+pin1[i]+string("]");
	                      else s += string("-p")+Order(pin1[i]);
	}
	
	AddToVariables(flag,pin1sum,s,1,sf);
	//sf<<"  Vec4D p"<<Order(pin1sum)<<" = "<<s<<";"<<endl ;
	
	/*
	  for (short int i=1;i<pin1.size();i++) {
	  if (pin1[i].length()==1) sf<<"-p["<<pin1[i]<<"]";
	  else sf<<"-p"<<pin1[i];
	  }
	  sf<<";"<<endl;
	*/
      }
      else {
	//later
      }

      pout1.clear();
      for (short int i=1;i<pout0.size();i++) pout1.push_back(pout0[i]);
      string help = pout0[0];
      pout0.clear();
      pout0.push_back(help);    
      count = pin0.size()-1;
    }
    else {
      if (pout1.size()>1) {
	//pout0  -> pin0
	if (pout0.size()==1) {
	  pin0.push_back(pout0[0]);
	  string pin0sum = string("0_");
	  for (short int i=1;i<pin0.size();i++) pin0sum += pin0[i];

	  //sf<<"  Vec4D p"<<Order(pin0sum)<<" = p[0]";

	  string s = string("p[0]");
	  for (short int i=1;i<pin0.size();i++) {
	    if (pin0[i].length()==1) s+=string("-p[")+pin0[i]+string("]");
	                        else s+=string("-p")+Order(pin0[i]);
	  }
	  //sf<<";"<<endl;
	  
	  AddToVariables(flag,pin0sum,s,1,sf);
	  
	}
	else {
	  //later
	}

	pout0.clear();
	for (short int i=0;i<pout1.size()-1;i++) pout0.push_back(pout1[i]);
	string help = pout1[pout1.size()-1];
	pout1.clear();
	pout1.push_back(help); 
	count = tcount-(pin1.size()-1)-1;
      }
    }

    SingleTStep(flag,s,propt,tcount,rannum,sf,count,pin0,pin1,pout0,pout1);
  }
}

void Channel_Generator::GenerateMasses(int flag,Point** _plist,int pcount,
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
      AddToVariables(flag,lm[i],string("ms[")+lm[i]+string("]"),0,sf);
      //sf<<"  double s"<<lm[i]<<" = ms["<<lm[i]<<"];"<<endl;
      momp[i]  = string("p[") + lm[i] + string("]");
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
    CalcSmin(flag,"restmin",help,sf,0);
    sum_s_i = string("-sqrt(s")+Order(help)+string("_restmin)");
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
    smax += string(");");
    
    AddToVariables(flag,lm[hit] +string("_max"),smax,0,sf);
    //sf<<"  double s"<<Order(lm[hit])<<"_max = "<<smax<<endl;

    if (maxpole>0.) {
      int hi = (_plist[hit]->fl).Kfcode();
      sf<<"  Flavour fl"<<lm[hit]<<" = "<<"Flavour(kf::code("<<hi<<"));"<<endl;
    }
    if ((flag==0) || (flag==10)) {
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
    }
    //    if ((flag==1) || (flag==11)) {
    else {
      string s(""); 
      for (short int i=0;i<lm[hit].length()-1;i++) s += string("p[")+lm[hit][i]+string("]+");
      s += string("p[")+lm[hit][lm[hit].length()-1]+string("]");
     
      AddToVariables(flag,lm[hit],s,1,sf);
      AddToVariables(flag,lm[hit],string("dabs(")+momp[hit]+string(".Abs2())"),0,sf);
      if (maxpole>0.) {
	sf<<"  wt *= CE.MassivePropWeight(fl"<<lm[hit]<<".Mass(),"<<"fl"<<lm[hit]<<".Width(),1,"
	  <<"s"<<Order(lm[hit])<<"_min,s"<<Order(lm[hit])<<"_max,"<<"s"<<Order(lm[hit])<<");"<<endl;
      }
      else {
	sf<<"  wt *= CE.MasslessPropWeight(1.,s"<<Order(lm[hit])<<"_min,"
	  <<"s"<<Order(lm[hit])<<"_max,s"<<Order(lm[hit])<<");"<<endl;
      }
    }
    sum_s_i  += string("-sqrt(s")+Order(lm[hit])+string(")");
    sflag[hit] = 1;
  }
  delete[] lm;
  delete[] momp;
  delete[] sflag;
}

void Channel_Generator::SetProps(Point* p,Point** props,Point** propt, int& count)
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
      if (p->right->b == -1 && p->right->number<99) props[count] = p->left;
                                               else props[count] = p->right;
      return;
    }
  }
  
  count++; 
  SetProps(propt[count-1],props,propt,count);
}

void Channel_Generator::CalcTSmin(int flag,vector<string>& p,ofstream& sf)
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
    else CalcSmin(flag,"tmpmin",help,sf,0);    
    s = string("sqr(sqrt(s")+Order(help)+string("_tmpmin)");
  }
  else s = string("sqr(");


  for (short int i=0;i<p.size();i++) {
    if (p[i].length()>1) s += string("+sqrt(s") + Order(p[i]) + string(")");
  }
  s += string(")");

  sf<<"  double s"<<Order(psum)<<"_min = "<<s<<";"<<endl;
}


void Channel_Generator::CalcSmin2(Point* p,string& s)
{
  if (p->left==0) return;

  if (p->left->left==0 && p->right->left==0) {    
    s += string("+sqrt(cuts->scut[")+IString(p->left->number)+
         string("][")+IString(p->right->number)+string("])");
    return;
  }

  if (p->left->left==0 && p->right->left!=0) {
    s += string("+sqrt(ms[")+IString(p->left->number)+string("])");
    CalcSmin2(p->right,s);
    return;
  }

  if (p->left->left!=0 && p->right->left==0) {
    s += string("+sqrt(ms[")+IString(p->right->number)+string("])");
    CalcSmin2(p->left,s);
    return;
  }

  CalcSmin2(p->left,s);
  CalcSmin2(p->right,s);
}

void Channel_Generator::CalcSmin(int flag,char* min,string lm,ofstream& sf,Point* p)
{
  // sum_{i<j} sij + (2-n)*sum_i m_i^2

  string s("");

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
  }
}


string Channel_Generator::LinkedMasses(Point* p)
{
  if (p->left==0) {
    char help[4];
    sprintf(help,"%i",p->number);
    return string(help);
  }
  return LinkedMasses(p->left)+LinkedMasses(p->right);
}


void Channel_Generator::IdentifyProps(Point* _plist)
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

void Channel_Generator::BackLinks(Point* p,Point* &endp)
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

void Channel_Generator::InitT(Point* p)
{
  p->t = 0;
  if (p->left==0) return;
  InitT(p->left);
  InitT(p->right);
}

string Channel_Generator::IString(int i)
{
  MyStrStream sstr;
  sstr<<i;
  string istr;
  sstr>>istr;
  return istr;
}

string Channel_Generator::Order(string s)
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

void  Channel_Generator::AddToVariables(int flag,const string& lhs,const string& rhs,const int& type,
					ofstream& sf)
{
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
      msg.Error()<<" ERROR in Channel_Generator::AddToVariables ()"<<endl;
      abort();
    }
  }
}

