#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/IO_Handler.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

#define SQRT_05 0.70710678118654757

std::ostream& AMEGIC::operator<<(std::ostream& os, const Momfunc& mf) {
  os<<mf.type<<";"<<mf.argnum;
  for (int i=0;i<mf.argnum;i++) os<<","<<mf.arg[i];
  os<<","<<mf.angle<<","<<0<<","<<0<<";";
  return os;
}

std::istream& AMEGIC::operator>>(std::istream& is, Momfunc& mf) {
  std::string in;
  is>>in;
  if (in.length()==0) THROW(critical_error,"String to momfunc translation failed.");
  size_t pos=in.find(";");
  mf.type = mt::momtype(ToType<int>(in.substr(0,pos)));
  if (pos!=std::string::npos) in=in.substr(pos+1);
  else in="";
  pos=in.find(",");
  mf.argnum = ToType<int>(in.substr(0,pos));
  if (pos!=std::string::npos) in=in.substr(pos+1);
  else in="";
  if (mf.arg) delete[] mf.arg;
  mf.arg=new int[mf.argnum];
  for (int i=0;i<mf.argnum;i++) {
    pos=in.find(",");
    mf.arg[i] = ToType<int>(in.substr(0,pos));	
    if (pos!=std::string::npos) in=in.substr(pos+1);
    else in="";
  }
  pos=in.find(",");
  mf.angle = ToType<double>(in.substr(0,pos));
  if (pos!=std::string::npos) in=in.substr(pos+1);
  else in="";
  pos=in.find(",");
  mf.mass = ToType<double>(in.substr(0,pos));
  if (pos!=std::string::npos) in=in.substr(pos+1);
  else in="";
  pos=in.find(";");
  mf.cplxmass2 = ToType<Complex>(in.substr(0,pos+1));
  return is;
}

Basic_Sfuncs::Basic_Sfuncs(int _nmom,int _nvec, Flavour* flav,int* _b) 
  : fl(flav), nmom(_nmom), nvec(_nvec), b(_b)
{
  momcount = InitializeMomlist();
  _eta=_mu=0;
  _S0=_S1=0;
  m_precalc=0;
  Setk0(10);
  p_epol=NULL;
}

Basic_Sfuncs::Basic_Sfuncs(int _nmom,int _nvec, Flavour* flav,int* _b,string name,string name2) 
  : fl(flav), nmom(_nmom), nvec(_nvec), b(_b)
{
  _eta=_mu=0;
  _S0=_S1=0;
  m_precalc=1;
  Setk0(10);
  name+="/Sfunc.dat";

  IO_Handler ioh;
  ioh.SetFileNameRO(name);

  ioh.GetFstream()>>momcount;
  for (int i=0;i<momcount;i++) {
    Momfunc mf;
    ioh.GetFstream()>>mf;
    Momlist.push_back(mf);    
  }
  Initialize();
  for (int i=0;i<momcount;i++) delete[] calc_st[i];
  delete[] calc_st;
  calc_st = ioh.MatrixInput<int>("",momcount,momcount);
  UpdateMasses(name2);
}

void Basic_Sfuncs::UpdateMasses(string name)
{
  int cnt=0;
  ifstream is;
  is.open(name.c_str());
  string str;
  for (;is;) {
    getline(is,str); 
    if (str.find(string("fl"))==0) {
      int a=str.find("[");
      int b=str.find("]");
      int idx = ToType<int>(str.substr(a+1,a-b-1));
      if (idx<momcount) {
	cnt++;
	a = str.find("=");
	kf_code kfc = ToType<int>(str.substr(a+1));
	Momlist[idx].kfc = kfc;
	Flavour flav(kfc);
	Momlist[idx].mass = flav.Mass();
	if (idx>=nmom) Momlist[idx].cplxmass2 = Complex(sqr(flav.Mass()),-flav.Width()*flav.Mass());
      }
      else {
	THROW(critical_error,"Inconsistent flavour entry in *.map");
      }
    }
  }
  if (cnt!=momcount) THROW(critical_error,"Missing flavour in *.map");
  is.close();  
}

void Basic_Sfuncs::Output(string name)
{
  name+="/Sfunc.dat";
  IO_Handler ioh;
  ioh.SetFileName(name);
  ioh.Output("",momcount);
  for (int i=0;i<momcount;i++) ioh.GetFstream()<<Momlist[i]<<endl;
  ioh.MatrixOutput<int>("",calc_st,momcount,momcount);
}

void Basic_Sfuncs::WriteMomFlavs(ofstream& os)
{
  for (int i=0;i<momcount;i++) {
    os<<"fl["<<i<<"]="<<Momlist[i].kfc<<endl;
  }
}

Basic_Sfuncs::~Basic_Sfuncs() 
{
  if (_eta) delete[] _eta;
  if (_mu) delete[] _mu;
  if (_S0) {
    for (int i=0;i<momcount;i++) {
      delete[] _S0[i];
      delete[] _S1[i];
      delete[] calc_st[i];
    }
    delete[] _S0;
    delete[] _S1;
    delete[] calc_st;
  }
}

void Basic_Sfuncs::Initialize()
{
  //S--Functions
  _eta = new Complex[momcount];
  _mu  = new Complex[momcount];

  _S0  = new Complex*[momcount];
  _S1  = new Complex*[momcount];
  calc_st = new int*[momcount];
  for (int i=0;i<momcount;i++) {
    _S0[i] = new Complex[momcount];
    _S1[i] = new Complex[momcount];   
    calc_st[i] = new int[momcount];
  }

}


int Basic_Sfuncs::InitializeMomlist()
{
  for (int i=0;i<nvec;i++) {
    Momfunc Mom;
    Mom.argnum = 1;
    Mom.arg    = new int[Mom.argnum];
    Mom.arg[0] = i;
    Mom.type=mt::mom;
    Mom.mass   = fl[i].Mass();
    Mom.kfc    = fl[i].Kfcode();
    Momlist.push_back(Mom);
  }
  return nvec; 
}

int Basic_Sfuncs::GetMomNumber(Pfunc* p)
{
  for(size_t k=0;k<Momlist.size();k++) {
    if (Momlist[k].argnum==p->argnum) {
      int hit = 0;
      for (int i=1;i<Momlist[k].argnum;i++) { 
	hit = 0;
	for (int j=1;j<p->argnum;j++) {
	  if (Momlist[k].arg[i]==p->arg[j]) {
	    hit = 1;
	    break;
	  }
	}
	if (hit==0) break;
      }
      if (hit==1) return Momlist[k].arg[0];
    }
  }
  return -1;
}

int Basic_Sfuncs::BuildMomlist(Pfunc_List& pl) 
{
  for (Pfunc_Iterator pit=pl.begin();pit!=pl.end();++pit) {
    Pfunc* p = *pit;
    int n = GetMomNumber(p);
    if (n==-1) {
      Momfunc* Mom;
      Mom = new Momfunc;

      Mom->argnum = p->argnum;
      Mom->arg    = new int[Mom->argnum];
      Mom->arg[0] = momcount;
      Mom->type   = mt::prop;
      Mom->kfc    = p->fl.Kfcode();
      p->momnum   = momcount;
      for (int i=1;i<p->argnum;i++) 
	Mom->arg[i] = p->arg[i];

      if(p->argnum==(nmom-1) && b[1]==-1){
        int hit=1;
        for(int i=1;i<p->argnum;i++) if(p->arg[i]<2)hit=0;
	if(hit==1) Mom->type=mt::cmprop;
      }
      Momlist.push_back(*Mom);
      momcount++;
      n = Momlist.size()-1;
      
      while(momcount>93&&momcount<100){
       Mom->arg[0] = momcount;
       Mom->type=mt::p_none;
       Momlist.push_back(*Mom);
       momcount++;
      }
      delete Mom;
      
    }
    else p->momnum = n;
    if (p->haspol) BuildPolarisations(n,p->fl); 	
  }
  return momcount;
}

void Basic_Sfuncs::PropPolarisation(int pindex,Pfunc_List& pl,vector<int>& iargs)
{
  int momindex = -1;
  Flavour momfl;
  for (Pfunc_Iterator pit=pl.begin();pit!=pl.end();++pit) {
    Pfunc* p = *pit;
    if (p->arg[0]==pindex) {
      momindex = p->momnum;
      momfl    = p->fl; 
      break;
    }
  }

  if(!momfl.IsScalar()){
    for(size_t k=nmom;k<Momlist.size();k++) {
      if (Momlist[k].arg[1]==momindex) {
	switch(Momlist[k].type)
	  {
	  case mt::p_s:
	    if(ATOOLS::IsZero(momfl.Mass()-Momlist[k].mass))iargs.push_back(Momlist[k].type);
	    break;
	  case mt::p_si:
	    if(ATOOLS::IsZero(momfl.Mass()))iargs.push_back(Momlist[k].type);
	    break;
	  default:
	    iargs.push_back(Momlist[k].type);
	  }
      }
    }
  }
  else iargs.push_back(0);  
}


int Basic_Sfuncs::GetPolNumber(int momindex, int sign,double mass,int check)
{
  for(size_t k=nmom;k<Momlist.size();k++) {
    if (Momlist[k].type==sign) 
      if (Momlist[k].arg[1]==momindex && (sign!=mt::p_s || Momlist[k].mass==mass)) return k;
    if (Momlist[k].type==mt::p_spec && Momlist[k].arg[1]==momindex && Momlist[k].arg[0]==sign) return k;
  }
  if (check==0){
    msg_Error()<<"******Get_Pol_Number: Not Found! "<<momindex<<" "<<sign<<" Mass:"<<mass<<endl;
    abort();
  }
  return -1;
}

void Basic_Sfuncs::PrintMomlist()
{
  return;

  msg_Out()<<"Momlist: "<<endl;
  for(size_t k=0;k<Momlist.size();k++) {
    msg_Out()<<Momlist[k].arg[0]<<" --> ";
    for (int i=1;i<Momlist[k].argnum;i++) msg_Out()<<Momlist[k].arg[i]<<",";
    msg_Out()<<" type = "<<Momlist[k].type<<endl;
  }
}

int Basic_Sfuncs::BuildTensorPolarisations(int momindex) 
{
  //Add polarisation vectors to construct tensors for external spin 2 particles
  if (momindex>nvec) {
    msg_Error()<<"*****BuildTensorPolarisations: Not an external momentum!"<<endl;
    return 0;
  }
  Momfunc* Mom;
  Mom = new Momfunc;

  //Polarisation -1
  Mom->argnum = 2;
  Mom->arg    = new int[Mom->argnum];
  Mom->arg[0] = momcount;
  Mom->arg[1] = momindex;
  Mom->type   = mt::p_m;  
  Mom->mass   = Momlist[momindex].mass;
  Mom->kfc    = Momlist[momindex].kfc;
  momcount++;
  Momlist.push_back(*Mom);

  //Polarisation +1
  Mom->arg[0] = momcount;
  Mom->type=mt::p_p;
  momcount++;
  Momlist.push_back(*Mom);

  //Polarisation longitudinal
  Mom->arg[0] = momcount;
  Mom->type=mt::p_l;
  momcount++;
  Momlist.push_back(*Mom);
  return momcount;
}


int Basic_Sfuncs::BuildPolarisations(int momindex,char type,double angle) 
{
  //Add polarisation vectors for external particles
  if (momindex>nvec) {
    msg_Error()<<"*****BuildPolarisations: Not an external momentum!"<<endl;
    return 0;
  }
  Momfunc mom_func;

  if (type=='e'){
    mom_func.type = mt::p_spec;
    mom_func.kfc  = Momlist[momindex].kfc;
    mom_func.argnum = 2;
    mom_func.arg    = new int[mom_func.argnum];
    mom_func.arg[1]=momindex;
    mom_func.arg[0]=90;
    Momlist.push_back(mom_func);
    mom_func.arg[0]=91;
    Momlist.push_back(mom_func);
    momcount+=2;
    return momcount;
  }
  //Polarisation -1
  mom_func.argnum = 3;
  mom_func.arg    = new int[mom_func.argnum];
  mom_func.arg[0] = momcount;
  mom_func.arg[1] = momindex;
  mom_func.arg[2] = 0;
  switch(type){
  case 'l':
    mom_func.type=mt::p_l0;
    mom_func.angle=angle/180*M_PI;
    break;
  case '+':
    mom_func.arg[2] = +1;
    mom_func.type=mt::p_m;
    break;
  case '-':
    mom_func.arg[2] = -1;
    mom_func.type=mt::p_m;
    break;
  default:
    mom_func.type=mt::p_m;
  }
  mom_func.mass=Momlist[momindex].mass;
  mom_func.kfc =Momlist[momindex].kfc;
  momcount++;
  Momlist.push_back(mom_func);

  //Polarisation +1
  mom_func.arg[0] = momcount;
  switch(type){
  case 'l':
    mom_func.type=mt::p_l1;
    mom_func.angle=angle/180.*M_PI;
    break;
  default:mom_func.type=mt::p_p;
  }
  mom_func.mass=Momlist[momindex].mass;
  mom_func.kfc =Momlist[momindex].kfc;
  momcount++;
  Momlist.push_back(mom_func);
  if(ATOOLS::IsZero(Momlist[momindex].mass))  return momcount;

  //Polarisation longitudinal
  mom_func.arg[0] = momcount;
  mom_func.type=mt::p_l;
  mom_func.mass=Momlist[momindex].mass;
  mom_func.kfc =Momlist[momindex].kfc;
  momcount++;
  Momlist.push_back(mom_func);
  return momcount;
}

int Basic_Sfuncs::BuildPolarisations(int momindex, Flavour fl) 
  //Add polarisation vectors for cutted propagators
{
  if (momindex<nvec) {
    msg_Error()<<"*****BuildPolarisations: Not an internal momentum!"<<endl;
    return 0;
  }
  double Mass = fl.Mass();
  Complex Mass2= Complex(sqr(Mass),0.);
  if(!ATOOLS::IsZero(fl.Width()))
      Mass2-=Complex(0.,fl.Width()*Mass);
  Momfunc* Mom = new Momfunc;
  Mom->argnum = 2;
  Mom->arg    = new int[Mom->argnum];
  Mom->arg[1] = momindex;
  Mom->mass  = Mass;
  Mom->cplxmass2 = Mass2;
  Mom->kfc = fl.Kfcode(); 

  if (GetPolNumber(momindex,mt::p_lh,0,1)==-1) {
    //Polarisation -1
    Mom->arg[0] = momcount;
    Mom->type=mt::p_lh;
    momcount++;
    Momlist.push_back(*Mom);

    //Polarisation +1
    Mom->arg[0] = momcount;
    Mom->type=mt::p_lp;
    momcount++;
    Momlist.push_back(*Mom);
    if(momindex<nvec&&ATOOLS::IsZero(Mass))  return momcount;

    //Polarisation longitudinal
    Mom->arg[0] = momcount;
    Mom->type=mt::p_l;
    momcount++;
    Momlist.push_back(*Mom);
    }
  
  if(ATOOLS::IsZero(Mass)){ 
    if (GetPolNumber(momindex,mt::p_si,0,1)==-1 && rpa->gen.CutScheme()!=1 ) {
      Mom->arg[0] = momcount;
      Mom->type=mt::p_si;
      momcount++;
      Momlist.push_back(*Mom);
    }
    return momcount;
  }
  if (GetPolNumber(momindex,mt::p_s,Mass,1)==-1) {
    //Polarisation scalar
    Mom->arg[0] = momcount;
    Mom->type=mt::p_s;
    momcount++;
    Momlist.push_back(*Mom);
  }
 
  return momcount;
}

double Basic_Sfuncs::Norm(int i,int j)
{
  return 1./(sqrt(2.)*abs(S0(i,j)));
}

void Basic_Sfuncs::CalcMomlist()
{
  double ps,pt;
  Vec4D mom,vh1,vh2;
  Complex help;
  for(size_t j=0;j<Momlist.size();j++){
    Momlist[j].mom_img=Vec4D(0.,0.,0.,0.);
    switch(Momlist[j].type){
    case mt::p_none : break;
    case mt::mom :    Momlist[j].mom = p[j];break;
    case mt::prop : 
      Momlist[j].mom = Vec4D(0.,0.,0.,0.);
      for (int i=1;i<Momlist[j].argnum;i++) {
	Momlist[j].mom += b[Momlist[j].arg[i]]*p[Momlist[j].arg[i]];
      }
      break;
    case mt::cmprop :
      Momlist[j].mom = p[0]+p[1]; 
      if (ATOOLS::IsZero(Momlist[j].mom[3]/Momlist[j].mom[0])) Momlist[j].mom[3]=0.;
      break;
    case mt::p_m : {
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      double mom1=C1(mom), mom2=C2(mom), mom3=C3(mom);
      pt=sqrt(sqr(mom1)+sqr(mom2));
      if(!ATOOLS::IsZero(pt)){
	Momlist[j].mom = SQRT_05/ps*(mom3/pt*(mom1*K1()+mom2*K2())-pt*K3());
	Momlist[j].mom_img = -SQRT_05/pt*(mom1*K2()-mom2*K1());
      }
      else {
	Momlist[j].mom = SQRT_05*K1();
	if(mom3>0)Momlist[j].mom_img = -SQRT_05*K2();
	else Momlist[j].mom_img = SQRT_05*K2();
       }
      Momlist[j+1].mom = Momlist[j].mom;
      Momlist[j+1].mom_img = (-1.)*Momlist[j].mom_img;
      j++;
      break;
    }
    case mt::p_l0 : {
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      double mom1=C1(mom), mom2=C2(mom), mom3=C3(mom);
      pt=sqrt(sqr(mom1)+sqr(mom2));
      if(!ATOOLS::IsZero(pt)){
	vh1 = 1.0/ps*(mom3/pt*(mom1*K1()+mom2*K2())-pt*K3());
	vh2 = 1.0/pt*(mom1*K2()-mom2*K1());
      }
      else {
	if(mom3>0) vh1 = K1();
	else vh1 = -K1();
	vh2 = K2();
       }
      Momlist[j].mom = ::cos(Momlist[j].angle)*vh1+::sin(Momlist[j].angle)*vh2;
      Momlist[j+1].mom = ::sin(Momlist[j].angle)*vh1-::cos(Momlist[j].angle)*vh2;
      j++;
      break;
    }
    case mt::p_lh : {
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      double mom1=C1(mom), mom2=C2(mom), mom3=C3(mom);
      pt=sqrt(sqr(mom1)+sqr(mom2));
      if(!ATOOLS::IsZero(pt)){
	vh1 = 1.0/ps*(mom3/pt*(mom1*K1()+mom2*K2())-pt*K3());
	vh2 = 1.0/pt*(mom1*K2()-mom2*K1());
	if(mom1+mom2<0)vh2=-1.*vh2;
      }
      else {
	vh1 = SQRT_05*(K1()-K2());
	vh2 = SQRT_05*(K1()+K2());
      }
      Momlist[j].mom = vh1;
      Momlist[j+1].mom = vh2;      
      j++;
      break;
    }
    case mt::p_l :     
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      if(!ATOOLS::IsZero(ps)){
	Momlist[j].mom=1./sqrt(abs(sqr(mom[0])-sqr(ps)))*Vec4D(ps,mom[0]*mom[1]/ps,
							       mom[0]*mom[2]/ps,
							       mom[0]*mom[3]/ps);
      }
      else Momlist[j].mom=K3();
      if (mom.Abs2()<0.) {
	Momlist[j].mom_img = Momlist[j].mom;
	Momlist[j].mom = Vec4D(0.,0.,0.,0.);
      }
      break;
    case mt::p_s :
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=mom.Abs2();
 
      if(ATOOLS::IsZero(Momlist[j].mass))Momlist[j].mom=Vec4D(0.,0.,0.,0.);
      else {
	  help=sqrt((Complex(1.,0.)-Momlist[j].cplxmass2/ps)/Momlist[j].cplxmass2);
	  Momlist[j].mom=real(help)*mom;
	  Momlist[j].mom_img=imag(help)*mom;
      }
      break;
    case mt::p_si :
      mom                = Momlist[Momlist[j].arg[1]].mom;
      help               = csqrt(-1./mom.Abs2());
      Momlist[j].mom     = real(help)*mom;
      Momlist[j].mom_img = imag(help)*mom;
      break;
    case mt::p_spec:
      if (p_epol) Momlist[j].mom = (*p_epol)[Momlist[j].arg[0]-90];
      else msg_Error()<<"Error in Basic_Sfuncs::CalcMomlist(): no external polarization array!!!"<<endl;
      break;
    default : break;
    }
  }
}

void Basic_Sfuncs::InitGaugeTest(double theta)
{
  double ps,pt;
  double s,c,s0,c0,sf,cf;
  Momfunc* m;
  Vec4D mom;
  for(size_t j=0;j<Momlist.size();j++){
    if(Momlist[j].type==mt::p_m&&ATOOLS::IsZero(Momlist[j].mass)){
      if (Momlist[j].arg[2]==+1 || Momlist[j].arg[2]==-1) {
	// do nothing
      } 
      else {
	mom=Momlist[Momlist[j].arg[1]].mom;
	ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
	pt=sqrt(sqr(mom[1])+sqr(mom[2]));
	s=sqrt(.5*(1-mom[3]/ps));
	c=sqrt(.5*(1+mom[3]/ps));
	s0=::sin(theta*0.5);
	c0=::cos(theta*0.5);
	if(!ATOOLS::IsZero(pt)){
	  sf=mom[2]/pt;
	  cf=mom[1]/pt;
	}
	else {sf=0.;cf=1.;}

	Momlist[j].mom=1./sqrt(1.-mom[1]/ps*::sin(theta)-mom[3]/ps*::cos(theta))*
	  Vec4D(c0*c+s0*s*cf,s0*c+s*c0*cf,s*c0*sf,c0*c-s*s0*cf);
	Momlist[j].mom_img=1./sqrt(1.-mom[1]/ps*::sin(theta)-mom[3]/ps*::cos(theta))*
	  Vec4D(s0*s*sf,s*c0*sf,s0*c-s*c0*cf,-s*s0*sf);

	Momlist[j+1].mom = Momlist[j].mom;
	Momlist[j+1].mom_img = (-1.)*Momlist[j].mom_img;
      }
    }    
    if(ATOOLS::IsZero(Momlist[j].mass))if(Momlist[j].type==mt::p_m||Momlist[j].type==mt::p_p){
      m=&Momlist[j];
      switch(k0_n){
      case 10:
	if(!IsComplex(j))
	  _eta[j] = csqrt(2.*(m->mom[0]+m->mom[R3()]));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]+m->mom[R3()]),
				    2.*(m->mom_img[0]+m->mom_img[R3()])));
        break;
      case 1:
	if(!IsComplex(j)) 
	  _eta[j] = csqrt(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*SQRT_05));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*SQRT_05),
				    2.*(m->mom_img[0]-(m->mom_img[2]+m->mom_img[3])*SQRT_05)));
	break;
      case 2:
	if(!IsComplex(j)) 
	  _eta[j] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*SQRT_05));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*SQRT_05),
				    2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[2])*SQRT_05)));
	break;
      default:
	if(!IsComplex(j)) 
	  _eta[j] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*SQRT_05));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*SQRT_05),
				    2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[3])*SQRT_05)));
      }
    }
  }
}

int Basic_Sfuncs::CalcEtaMu(Vec4D* _p)
{
  //  PROFILE_HERE;
  //  PROFILE_LOCAL("int Basic_Sfuncs::setS(Vec4D* _p)");
  // _eta's and _mu's precalc
  
  p = _p;
  CalcMomlist();

  Momfunc* m;
  
  int etachk=1;
  double zchk;

  for(size_t i=0;i<Momlist.size();i++) {
    m=&Momlist[i];

    switch(k0_n){
    case 10:
      if(!IsComplex(i))
 	_eta[i] = csqrt(2. * (m->mom[0]+m->mom[R3()]));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]+m->mom[R3()]),
				  2.*(m->mom_img[0]+m->mom_img[R3()])));
      break;
    case 1:
      if(!IsComplex(i)) 
	_eta[i] = csqrt(2. * (m->mom[0] - (m->mom[2]+m->mom[3])*SQRT_05 ));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*SQRT_05),
				  2.*(m->mom_img[0]-(m->mom_img[2]+m->mom_img[3])*SQRT_05)));
      break;
    case 2:
      if(!IsComplex(i)) 
	_eta[i] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*SQRT_05));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*SQRT_05),
				  2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[2])*SQRT_05)));
      break;
    default:
      if(!IsComplex(i)) 
	_eta[i] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*SQRT_05));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*SQRT_05),
				  2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[3])*SQRT_05)));
    }
    if(ATOOLS::IsZero(_eta[i]))etachk=0;
    if ((int)i<nmom) {
      if ((i==0)) _mu[0]  = sqrt(dabs(Momlist[i].mom.Abs2()))/_eta[0];
      else _mu[i]  = Momlist[i].mass/_eta[i];
      //if (b[i]==1 && 
      if (fl[i].IsAnti() && i!=0) 
	_mu[i] = - _mu[i];
      if (b[i]==-1 && !(fl[i].IsAnti()) && i==0) 
	_mu[i] = - _mu[i];
      
      if (fl[i].MassSign()==-1) _mu[i] = - _mu[i];
    }
    else {
	switch(m->type){
	    case mt::p_p:
	    case mt::p_m: _mu[i] = Complex(0.,0.);
		break;
	    case mt::p_l: _mu[i] = Complex(0.,1.)/_eta[i];
	        break;
	    case mt::p_si:
		_mu[i]  = (csqrt((m->mom).Abs2())+csqrt(-(m->mom_img).Abs2()))/_eta[i];
		break;
	    case mt::p_s:
		_mu[i] = sqrt(Complex((m->mom).Abs2()-(m->mom_img).Abs2(),
				      2 * m->mom * m->mom_img))/_eta[i];
		break;
	    default: 
	      zchk = (m->mom).Abs2();
	      if (IsZero(zchk)) _mu[i] = Complex(0.,0.);
	      else _mu[i] = csqrt(zchk)/_eta[i];
	}
    }
  }
  if (!m_precalc){
    for(int i=0;i<momcount;i++)
      for(int j=0;j<momcount;j++)
	calc_st[i][j]=0;
  }
  else PrecalcS();
  return etachk;
}

void Basic_Sfuncs::PrecalcS()
{
  for(int i=1;i<momcount;i++)
    for(int j=0;j<i;j++)
      if (calc_st[i][j]) {
	Momfunc* m = &Momlist[i];
	Momfunc* m1 = &Momlist[j];
	Complex A= _eta[j]/_eta[i];

      switch (k0_n) {
      case 10:
	_S0[i][j] = Complex(m->mom[R1()],-m->mom[R2()])*A;
	_S1[i][j] = Complex(-m->mom[R1()],-m->mom[R2()])*A; 
	if (IsComplex(i)){
	  _S0[i][j] += Complex(m->mom_img[R2()],m->mom_img[R1()])*A;
	  _S1[i][j] += Complex(m->mom_img[R2()],-m->mom_img[R1()])*A;
	}      
	_S0[i][j] -= Complex(m1->mom[R1()],-m1->mom[R2()])/A;
	_S1[i][j] -= Complex(-m1->mom[R1()],-m1->mom[R2()])/A;
	if (IsComplex(j)){
	  _S0[i][j] -= Complex(m1->mom_img[R2()],m1->mom_img[R1()])/A;
	  _S1[i][j] -= Complex(m1->mom_img[R2()],-m1->mom_img[R1()])/A;
	}
	break;
      case 1:
	_S0[i][j] = Complex(m->mom[1],(m->mom[2]-m->mom[3])*SQRT_05)*A; 
	_S1[i][j] = Complex(-m->mom[1],(m->mom[2]-m->mom[3])*SQRT_05)*A;
	if (IsComplex(i)){
	  _S0[i][j] += Complex(-(m->mom_img[2]-m->mom_img[3])*SQRT_05,m->mom_img[1])*A;
	  _S1[i][j] += Complex(-(m->mom_img[2]-m->mom_img[3])*SQRT_05,-m->mom_img[1])*A;
	}      
	_S0[i][j] -= Complex(m1->mom[1],(m1->mom[2]-m1->mom[3])*SQRT_05)/A;
	_S1[i][j] -= Complex(-m1->mom[1],(m1->mom[2]-m1->mom[3])*SQRT_05)/A;
	if (IsComplex(j)){
	  _S0[i][j] -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*SQRT_05,m1->mom_img[1])/A;
	  _S1[i][j] -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*SQRT_05,-m1->mom_img[1])/A;
	}
	break;
      case 2:
	_S0[i][j] = Complex(m->mom[3],(m->mom[1]-m->mom[2])*SQRT_05)*A; 
	_S1[i][j] = Complex(-m->mom[3],(m->mom[1]-m->mom[2])*SQRT_05)*A;
	if (IsComplex(i)){
	  _S0[i][j] += Complex(-(m->mom_img[1]-m->mom_img[2])*SQRT_05,m->mom_img[3])*A;
	  _S1[i][j] += Complex(-(m->mom_img[1]-m->mom_img[2])*SQRT_05,-m->mom_img[3])*A;
	}
	_S0[i][j] -= Complex(m1->mom[3],(m1->mom[1]-m1->mom[2])*SQRT_05)/A;
	_S1[i][j] -= Complex(-m1->mom[3],(m1->mom[1]-m1->mom[2])*SQRT_05)/A;
	if (IsComplex(j)){
	  _S0[i][j] -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*SQRT_05,m1->mom_img[3])/A;
	  _S1[i][j] -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*SQRT_05,-m1->mom_img[3])/A;
	}
	break;
      default:
	_S0[i][j] = Complex(m->mom[2],(m->mom[3]-m->mom[1])*SQRT_05)*A; 
	_S1[i][j] = Complex(-m->mom[2],(m->mom[3]-m->mom[1])*SQRT_05)*A;
	if (IsComplex(i)){
	  _S0[i][j] += Complex(-(m->mom_img[3]-m->mom_img[1])*SQRT_05,m->mom_img[2])*A;
	  _S1[i][j] += Complex(-(m->mom_img[3]-m->mom_img[1])*SQRT_05,-m->mom_img[2])*A;
	}
	_S0[i][j] -= Complex(m1->mom[2],(m1->mom[3]-m1->mom[1])*SQRT_05)/A;
	_S1[i][j] -= Complex(-m1->mom[2],(m1->mom[3]-m1->mom[1])*SQRT_05)/A;
	if (IsComplex(j)){
	  _S0[i][j] -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*SQRT_05,m1->mom_img[2])/A;
	  _S1[i][j] -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*SQRT_05,-m1->mom_img[2])/A;
	}
      }
	_S0[j][i] = - _S0[i][j];
	_S1[j][i] = - _S1[i][j];	
#ifdef DEBUG__BS
	msg_Debugging()<<" - _S0["<<i<<"]["<<j<<"] = "<<-_S0[i][j]<<" (["<<i<<","<<j<<"])\n";
	msg_Debugging()<<" - _S1["<<j<<"]["<<i<<"] = "<<-_S1[j][i]<<" (<"<<j<<","<<i<<">)\n";
#endif
      } 
}

void Basic_Sfuncs::CalcS(int i, int j)
{
  //  PROFILE_HERE;
  if (i!=j) {
    Momfunc* m = &Momlist[i];
    Momfunc* m1 = &Momlist[j];
    Complex A= _eta[j]/_eta[i];

    switch(k0_n){
    case 10:
      _S0[i][j] = Complex(m->mom[R1()],-m->mom[R2()])*A; 
      _S1[i][j] = Complex(-m->mom[R1()],-m->mom[R2()])*A;
      if (IsComplex(i)){
	_S0[i][j] += Complex(m->mom_img[R2()],m->mom_img[R1()])*A;
	_S1[i][j] += Complex(m->mom_img[R2()],-m->mom_img[R1()])*A;
      }      
      _S0[i][j] -= Complex(m1->mom[R1()],-m1->mom[R2()])/A;
      _S1[i][j] -= Complex(-m1->mom[R1()],-m1->mom[R2()])/A;
      if (IsComplex(j)){
	_S0[i][j] -= Complex(m1->mom_img[R2()],m1->mom_img[R1()])/A;
	_S1[i][j] -= Complex(m1->mom_img[R2()],-m1->mom_img[R1()])/A;
      }
      break;
    case 1:
      _S0[i][j] = Complex(m->mom[1],(m->mom[2]-m->mom[3])*SQRT_05)*A; 
      _S1[i][j] = Complex(-m->mom[1],(m->mom[2]-m->mom[3])*SQRT_05)*A;
      if (IsComplex(i)){
	_S0[i][j] += Complex(-(m->mom_img[2]-m->mom_img[3])*SQRT_05,m->mom_img[1])*A;
	_S1[i][j] += Complex(-(m->mom_img[2]-m->mom_img[3])*SQRT_05,-m->mom_img[1])*A;
      }      
      _S0[i][j] -= Complex(m1->mom[1],(m1->mom[2]-m1->mom[3])*SQRT_05)/A;
      _S1[i][j] -= Complex(-m1->mom[1],(m1->mom[2]-m1->mom[3])*SQRT_05)/A;
      if (IsComplex(j)){
	_S0[i][j] -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*SQRT_05,m1->mom_img[1])/A;
	_S1[i][j] -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*SQRT_05,-m1->mom_img[1])/A;
      }
      break;
    case 2:
      _S0[i][j] = Complex(m->mom[3],(m->mom[1]-m->mom[2])*SQRT_05)*A; 
      _S1[i][j] = Complex(-m->mom[3],(m->mom[1]-m->mom[2])*SQRT_05)*A;
      if (IsComplex(i)){
	_S0[i][j] += Complex(-(m->mom_img[1]-m->mom_img[2])*SQRT_05,m->mom_img[3])*A;
	_S1[i][j] += Complex(-(m->mom_img[1]-m->mom_img[2])*SQRT_05,-m->mom_img[3])*A;
      }
      _S0[i][j] -= Complex(m1->mom[3],(m1->mom[1]-m1->mom[2])*SQRT_05)/A;
      _S1[i][j] -= Complex(-m1->mom[3],(m1->mom[1]-m1->mom[2])*SQRT_05)/A;
      if (IsComplex(j)){
	_S0[i][j] -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*SQRT_05,m1->mom_img[3])/A;
	_S1[i][j] -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*SQRT_05,-m1->mom_img[3])/A;
      }
      break;
    default:
      _S0[i][j] = Complex(m->mom[2],(m->mom[3]-m->mom[1])*SQRT_05)*A; 
      _S1[i][j] = Complex(-m->mom[2],(m->mom[3]-m->mom[1])*SQRT_05)*A;
      if (IsComplex(i)){
	_S0[i][j] += Complex(-(m->mom_img[3]-m->mom_img[1])*SQRT_05,m->mom_img[2])*A;
	_S1[i][j] += Complex(-(m->mom_img[3]-m->mom_img[1])*SQRT_05,-m->mom_img[2])*A;
      }
      _S0[i][j] -= Complex(m1->mom[2],(m1->mom[3]-m1->mom[1])*SQRT_05)/A;
      _S1[i][j] -= Complex(-m1->mom[2],(m1->mom[3]-m1->mom[1])*SQRT_05)/A;
      if (IsComplex(j)){
	_S0[i][j] -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*SQRT_05,m1->mom_img[2])/A;
	_S1[i][j] -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*SQRT_05,-m1->mom_img[2])/A;
      }
    }
    if (IsZero(_S0[i][j])) _S0[i][j] = Complex(0.,0.);
    if (IsZero(_S1[i][j])) _S1[i][j] = Complex(0.,0.);
    _S0[j][i] = - _S0[i][j];
    _S1[j][i] = - _S1[i][j];
  }
  else {	
    _S0[i][j] = Complex(0.,0.);
    _S0[j][i] = Complex(0.,0.);
    _S1[i][j] = Complex(0.,0.);
    _S1[j][i] = Complex(0.,0.);
  }
  calc_st[i][j]=calc_st[j][i]=1;
}

Complex Basic_Sfuncs::CalcS(ATOOLS::Vec4D& m, ATOOLS::Vec4D& m1)
{
  Complex A,S;

  switch(k0_n){
  case 10:
    A = csqrt((m1[0]+m1[R3()])/(m[0]+m[R3()]));
    S = Complex(m[R1()],-m[R2()])*A; 
    S -= Complex(m1[R1()],-m1[R2()])/A;
    break;
  case 1:
    A = csqrt((m1[0]-(m1[2]+m1[3])*SQRT_05)/(m[0]-(m[2]+m[3])*SQRT_05));
    S = Complex(m[1],(m[2]-m[3])*SQRT_05)*A; 
    S -= Complex(m1[1],(m1[2]-m1[3])*SQRT_05)/A;
    break;
  case 2:
    A = csqrt((m1[0]-(m1[1]+m1[2])*SQRT_05)/(m[0]-(m[1]+m[2])*SQRT_05));
    S = Complex(m[3],(m[1]-m[2])*SQRT_05)*A; 
    S -= Complex(m1[3],(m1[1]-m1[2])*SQRT_05)/A;
    break;
  default:
    A = csqrt((m1[0]-(m1[1]+m1[3])*SQRT_05)/(m[0]-(m[1]+m[3])*SQRT_05));
    S = Complex(m[2],(m[3]-m[1])*SQRT_05)*A; 
    S -= Complex(m1[2],(m1[3]-m1[1])*SQRT_05)/A;
  }
  return S;
}

 
std::pair<Complex, Complex> Basic_Sfuncs::GetS(ATOOLS::Vec4D v, int j) 
{
  // This calculation is directly taken from CalcS; only minor difference is that
  // v is non-complex by definition so its complexity doesn't have to be checked.

  std::pair<Complex, Complex> S;
  
  Complex etav(csqrt(2. * Getk0() * v));
  if (ATOOLS::IsZero(etav)) {
    msg_Error()<<"An error occured in Basic_Sfunc::GetS(). The variable 'etav' was zero."
	       <<endl<<"This will cause a division by zero"<<endl;
  }

    Momfunc* m1 = &Momlist[j];
    Complex A= _eta[j]/etav;

    switch(k0_n){
    case 10:
      S.first   = Complex(v[R1()],-v[R2()])*A; 
      S.second  = Complex(-v[R1()],-v[R2()])*A;
      S.first  -= Complex(m1->mom[R1()],-m1->mom[R2()])/A;
      S.second -= Complex(-m1->mom[R1()],-m1->mom[R2()])/A;
      if (IsComplex(j)){
 	S.first  -= Complex(m1->mom_img[R2()],m1->mom_img[R1()])/A;
 	S.second -= Complex(m1->mom_img[R2()],-m1->mom_img[R1()])/A;
      }
      break;
    case 1:
      S.first   = Complex(v[1],(v[2]-v[3])*SQRT_05)*A; 
      S.second  = Complex(-v[1],(v[2]-v[3])*SQRT_05)*A;
      S.first  -= Complex(m1->mom[1],(m1->mom[2]-m1->mom[3])*SQRT_05)/A;
      S.second -= Complex(-m1->mom[1],(m1->mom[2]-m1->mom[3])*SQRT_05)/A;
      if (IsComplex(j)){
	S.first  -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*SQRT_05,m1->mom_img[1])/A;
	S.second -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*SQRT_05,-m1->mom_img[1])/A;
      }
      break;
    case 2:
      S.first  = Complex(v[3],(v[1]-v[2])*SQRT_05)*A; 
      S.second = Complex(-v[3],(v[1]-v[2])*SQRT_05)*A;
      S.first  -= Complex(m1->mom[3],(m1->mom[1]-m1->mom[2])*SQRT_05)/A;
      S.second -= Complex(-m1->mom[3],(m1->mom[1]-m1->mom[2])*SQRT_05)/A;
      if (IsComplex(j)){
	S.first  -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*SQRT_05,m1->mom_img[3])/A;
	S.second -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*SQRT_05,-m1->mom_img[3])/A;
      }
      break;
    default:
      S.first  = Complex(v[2],(v[3]-v[1])*SQRT_05)*A; 
      S.second = Complex(-v[2],(v[3]-v[1])*SQRT_05)*A;
      S.first  -= Complex(m1->mom[2],(m1->mom[3]-m1->mom[1])*SQRT_05)/A;
      S.second -= Complex(-m1->mom[2],(m1->mom[3]-m1->mom[1])*SQRT_05)/A;
      if (IsComplex(j)){
	S.first  -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*SQRT_05,m1->mom_img[2])/A;
	S.second -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*SQRT_05,-m1->mom_img[2])/A;
      }
    }
    if (IsZero(S.first))  S.first  = Complex(0.,0.);
    if (IsZero(S.second)) S.second = Complex(0.,0.);
    return S;
}

void Basic_Sfuncs::Setk0(int i)
{
  k0_n=i;
  //i=0: k0=Vec4D(1.,sqrt(.5),0.,sqrt(.5));
  //     k1=Vec4D(0.,0.,1.,0.);

  //i=1: k0=Vec4D(1.,0.,sqrt(.5),sqrt(.5));
  //     k1=Vec4D(0.,1.,0.,0.);

  //i=2: k0=Vec4D(1.,sqrt(.5),sqrt(.5),0.);
  //     k1=Vec4D(0.,0.,0.,1.);
  m_k1=SQRT_05*(Vec4D(0.0,Vec3D(-Getk0()))+Getk1());
  m_k2=Vec4D(0.0,cross(Vec3D(Getk1()),Vec3D(Getk0())));
  m_k3=SQRT_05*(Vec4D(0.0,Vec3D(-Getk0()))-Getk1());
  // if (k0_n<10) {
    m_k1=Vec4D(0.0,1.0,0.0,0.0);
    m_k2=Vec4D(0.0,0.0,1.0,0.0);
    m_k3=Vec4D(0.0,0.0,0.0,1.0);
  // }
  DEBUG_FUNC(k0_n);
  msg_Debugging()<<"k_0 = "<<Getk0()<<", k_1 = "<<Getk1()<<"\n";
  msg_Debugging()<<"k_x       = "<<K1()<<" "<<K1().Abs2()<<"\n";
  msg_Debugging()<<"k_y x k_z = "<<Vec4D(0.0,cross(Vec3D(K2()),Vec3D(K3())))<<"\n";
  msg_Debugging()<<"k_y       = "<<K2()<<" "<<K2().Abs2()<<"\n";
  msg_Debugging()<<"k_z x k_x = "<<Vec4D(0.0,cross(Vec3D(K3()),Vec3D(K1())))<<"\n";
  msg_Debugging()<<"k_z       = "<<K3()<<" "<<K3().Abs2()<<"\n";
  msg_Debugging()<<"k_x x k_y = "<<Vec4D(0.0,cross(Vec3D(K1()),Vec3D(K2())))<<"\n";
}

void Basic_Sfuncs::StartPrecalc()
{ 
  m_precalc=1; 
  int cnt=0;
  for(int i=0;i<momcount;i++)
    for(int j=0;j<momcount;j++)
      if (calc_st[i][j]){
	if (_S0[i][j]==Complex(0.,0.) && _S1[i][j]==Complex(0.,0.)) calc_st[i][j]=0;
	else if(i>j) ++cnt;
      }
}

bool Basic_Sfuncs::IsMomSum(int x,int y,int z)
{
  x = iabs(x); y = iabs(y); z = iabs(z); 
  if (Momlist[x].type==mt::p_s || Momlist[x].type==mt::p_si) x = Momlist[x].arg[1];
   else 
    if (Momlist[x].type==mt::p_l)
      if (Momlist[Momlist[x].arg[1]].mom[1]==0. &&
	  Momlist[Momlist[x].arg[1]].mom[2]==0. &&
	  Momlist[Momlist[x].arg[1]].mom[3]==0. ) x = Momlist[x].arg[1];
  if (Momlist[x].type!=mt::prop && 
      Momlist[x].type!=mt::cmprop) return false;
  if (Momlist[y].type!=mt::prop && 
      Momlist[y].type!=mt::cmprop && 
      Momlist[y].type!=mt::mom) return false;
  if (Momlist[z].type!=mt::prop && 
      Momlist[z].type!=mt::cmprop && 
      Momlist[z].type!=mt::mom) return false;

  Vec4D sum;
  if (Momlist[y].type==mt::mom) sum = b[y]*Momlist[y].mom;
  else sum = Momlist[y].mom;
  if (Momlist[z].type==mt::mom) sum += b[z]*Momlist[z].mom;
  else sum += Momlist[z].mom;
  return (sum==Momlist[x].mom);
}

ATOOLS::Vec4D Basic_Sfuncs::Getk0() {
  switch(k0_n) {
  case 10: return Spinor<double>::GetK0();
     case 1: return ATOOLS::Vec4D(1., 0, SQRT_05, SQRT_05);
     case 2: return ATOOLS::Vec4D(1., SQRT_05, SQRT_05, 0);
     default: return ATOOLS::Vec4D(1., SQRT_05, 0, SQRT_05);
  }
}

ATOOLS::Vec4D Basic_Sfuncs::Getk1() {
  switch(k0_n) {
  case 10: return Spinor<double>::GetK1();
     case 1: return ATOOLS::Vec4D(0., 1., 0., 0.);
     case 2: return ATOOLS::Vec4D(0., 0., 0., 1.);
     default: return ATOOLS::Vec4D(0., 0., 1., 0.);
  }
}
