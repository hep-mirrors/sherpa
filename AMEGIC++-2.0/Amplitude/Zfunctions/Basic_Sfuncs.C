#include "Basic_Sfuncs.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MathTools.H"
#include "prof.hh"

using namespace AMEGIC;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;

using namespace std;

Basic_Sfuncs::Basic_Sfuncs(int _nmom,int _nvec, Flavour* flav,int* _b) 
  : fl(flav), Nmom(_nmom), nvec(_nvec), b(_b)
{
  momcount = Initialize_Momlist();
  _eta=_mu=0;
  _S0=_S1=0;
  k0_n=0;
}

Basic_Sfuncs::~Basic_Sfuncs() 
{
  if (_eta) delete[] _eta;
  if (_mu) delete[] _mu;
  short int i;
  if (_S0) {
    for (i=0;i<momcount;i++) {
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
  for (short int i=0;i<momcount;i++) {
    _S0[i] = new Complex[momcount];
    _S1[i] = new Complex[momcount];   
    calc_st[i] = new int[momcount];
  }

}


int Basic_Sfuncs::Initialize_Momlist()
{
  for (short int i=0;i<nvec;i++) {
    Momfunc Mom;
    Mom.on = 1;
    Mom.argnum = 1;
    Mom.arg    = new int[Mom.argnum];
    Mom.arg[0] = i;
    Mom.type=mt::mom;
    Mom.mass =fl[i].Mass();
    Momlist.push_back(Mom);
  }
  return nvec; 
}

Vec4D Basic_Sfuncs::Mom_dir(int i)
{
  double sign=1.;-Sign(i);
  //return Vec4D(Momlist[i].mom[0],sign*Momlist[i].mom[1],sign*Momlist[i].mom[2],sign*Momlist[i].mom[3]);
  return sign*Momlist[i].mom;
}

int Basic_Sfuncs::Get_Mom_Number(Pfunc* p)
{
  for(short int k=0;k<Momlist.size();k++) {
    if (Momlist[k].argnum==p->argnum) {
      int hit = 0;
      for (short int i=1;i<Momlist[k].argnum;i++) { 
	hit = 0;
	for (short int j=1;j<p->argnum;j++) {
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

int Basic_Sfuncs::Build_Momlist(list<Pfunc*>& pl) 
{
  for (list<Pfunc*>::iterator pit=pl.begin();pit!=pl.end();++pit) {
    Pfunc* p = *pit;
    int n = Get_Mom_Number(p);
    if (n==-1) {
      Momfunc* Mom;
      Mom = new Momfunc;
      Mom->on = 1;

      Mom->argnum = p->argnum;
      Mom->arg    = new int[Mom->argnum];
      Mom->arg[0] = momcount;
      Mom->type=mt::prop;
      p->momnum   = momcount;
      for (short int i=1;i<p->argnum;i++) 
	Mom->arg[i] = p->arg[i];

      if(p->argnum==(Nmom-1) && b[1]==-1){
        int hit=1;
        for(short int i=1;i<p->argnum;i++) if(p->arg[i]<2)hit=0;
	if(hit==1) Mom->type=mt::cmprop;
      }
	
      Momlist.push_back(*Mom);
      //      cout<<"******Build_Momlist: "<<momcount<<endl;
      momcount++;
      n = Momlist.size()-1;
      
      while(momcount>93&&momcount<100){
       Mom->arg[0] = momcount;
       Mom->type=mt::p_none;
       Momlist.push_back(*Mom);
       momcount++;
       //  cout<<"******Build_Momlist(dummy): "<<momcount<<endl;
      }
      delete Mom;
      
    }
    else p->momnum = n;
    if (p->haspol) Build_Polarisations(n,p->fl); 	
  }
  return momcount;
}

void Basic_Sfuncs::PropPolarisation(int pindex,list<Pfunc*>& pl,vector<int>& iargs)
{
  int momindex = -1;
  Flavour momfl;
  for (list<Pfunc*>::iterator pit=pl.begin();pit!=pl.end();++pit) {
    Pfunc* p = *pit;
    if (p->arg[0]==pindex) {
      momindex = p->momnum;
      momfl    = p->fl; 
      break;
    }
  }

  if(!momfl.IsScalar()){
    for(short int k=Nmom;k<Momlist.size();k++) {
      if (Momlist[k].arg[1]==momindex) {
	switch(Momlist[k].type)
	  {
	  case mt::p_s:
	    if(AMATOOLS::IsZero(momfl.Mass()-Momlist[k].mass))iargs.push_back(Momlist[k].type);
	    break;
	  case mt::p_si:
	    if(AMATOOLS::IsZero(momfl.Mass()))iargs.push_back(Momlist[k].type);
	    break;
	  default:
	    iargs.push_back(Momlist[k].type);
	    //cout<<"PropPolarisations"<<iargs[iargs.size()-1]<<endl;
	  }
      }
    }
  }
  else iargs.push_back(0);  
}


int Basic_Sfuncs::Get_Pol_Number(int momindex, int sign,double mass,int check)
{
  for(short int k=Nmom;k<Momlist.size();k++) 
    if (Momlist[k].type==sign) 
      if (Momlist[k].arg[1]==momindex && (sign!=mt::p_s || Momlist[k].mass==mass)) return k;

  if (check==0){
    AORGTOOLS::msg.Error()<<"******Get_Pol_Number: Not Found! "<<momindex<<" "<<sign<<" Mass:"<<mass<<endl;
    abort();
  }
  return -1;
}

void Basic_Sfuncs::Print_Momlist()
{
  return;

  AORGTOOLS::msg.Out()<<"Momlist: "<<endl;
  for(short int k=0;k<Momlist.size();k++) {
    AORGTOOLS::msg.Out()<<Momlist[k].arg[0]<<" --> ";
    for (short int i=1;i<Momlist[k].argnum;i++) AORGTOOLS::msg.Out()<<Momlist[k].arg[i]<<",";
    AORGTOOLS::msg.Out()<<"on = "<<Momlist[k].on<<", type = "<<Momlist[k].type<<endl;
  }
}

int Basic_Sfuncs::Build_TensorPolarisations(int momindex) 
{
  //Add polarisation vectors to construct tensors for external spin 2 particles
  if (momindex>nvec) {
    msg.Error()<<"*****Build_TensorPolarisations: Not an external momentum!"<<endl;
    return 0;
  }
  Momfunc* Mom;
  Mom = new Momfunc;
  Mom->on = 1;

  //Polarisation -1
  Mom->argnum = 2;
  Mom->arg    = new int[Mom->argnum];
  Mom->arg[0] = momcount;
  Mom->arg[1] = momindex;
  Mom->type=mt::p_m;  
  Mom->mass=Momlist[momindex].mass;
  msg.Debugging()<<"*****Build_TensorPolarisations: -("<<Mom->type<<") "<<momcount<<" M:"<<Mom->mass<<endl;
  momcount++;
  Momlist.push_back(*Mom);

  //Polarisation +1
  Mom->arg[0] = momcount;
  Mom->type=mt::p_p;
  msg.Debugging()<<"*****Build_TensorPolarisations: + ("<<Mom->type<<") "<<momcount<<endl;
  momcount++;
  Momlist.push_back(*Mom);

  //Polarisation longitudinal
  Mom->arg[0] = momcount;
  Mom->type=mt::p_l;
  msg.Debugging()<<"*****Build_TensorPolarisations: l "<<momcount<<endl;
  momcount++;
  Momlist.push_back(*Mom);
  return momcount;
}


int Basic_Sfuncs::Build_Polarisations(int momindex,char type,double angle) 
{
  //Add polarisation vectors for external particles
  if (momindex>nvec) {
    msg.Error()<<"*****Build_Polarisations: Not an external momentum!"<<endl;
    return 0;
  }
  Momfunc* Mom;
  Mom = new Momfunc;
  Mom->on = 1;

  //Polarisation -1
  Mom->argnum = 2;
  Mom->arg    = new int[Mom->argnum];
  Mom->arg[0] = momcount;
  Mom->arg[1] = momindex;
  switch(type){
  case 'l':
    Mom->type=mt::p_l0;
    Mom->angle=angle/180*M_PI;
    break;
  default:Mom->type=mt::p_m;
  }
  Mom->mass=Momlist[momindex].mass;
  msg.Debugging()<<"*****Build_Polarisations: -("<<Mom->type<<") "<<momcount<<" M:"<<Mom->mass<<endl;
  momcount++;
  Momlist.push_back(*Mom);

  //Polarisation +1
  Mom->arg[0] = momcount;
  switch(type){
  case 'l':
    Mom->type=mt::p_l1;
    Mom->angle=angle/180.*M_PI;
    break;
  default:Mom->type=mt::p_p;
  }
  Mom->mass=Momlist[momindex].mass;
  msg.Debugging()<<"*****Build_Polarisations: + ("<<Mom->type<<") "<<momcount<<endl;
  momcount++;
  Momlist.push_back(*Mom);
  if(AMATOOLS::IsZero(Momlist[momindex].mass))  return momcount;

  //Polarisation longitudinal
  Mom->arg[0] = momcount;
  Mom->type=mt::p_l;
  Mom->mass=Momlist[momindex].mass;
  msg.Debugging()<<"*****Build_Polarisations: l "<<momcount<<endl;
  momcount++;
  Momlist.push_back(*Mom);
  return momcount;
}

int Basic_Sfuncs::Build_Polarisations(int momindex, Flavour fl) 
  //Add polarisation vectors for cutted propagators
{
  if (momindex<nvec) {
    msg.Error()<<"*****Build_Polarisations: Not an internal momentum!"<<endl;
    return 0;
  }
    //cout<<"Build Pols"<<endl;
double Mass = fl.Mass();
  Complex Mass2= Complex(sqr(Mass),0.);
  if(!AMATOOLS::IsZero(fl.Width()))
      Mass2-=Complex(0.,fl.Width()*Mass);
  Momfunc* Mom = new Momfunc;
  Mom->on = 1;
  Mom->argnum = 2;
  Mom->arg    = new int[Mom->argnum];
  Mom->arg[1] = momindex;
  Mom->mass  = Mass;
  Mom->cplxmass2 = Mass2;

  if (Get_Pol_Number(momindex,mt::p_lh,0,1)==-1) {
    //Polarisation -1
    Mom->arg[0] = momcount;
    Mom->type=mt::p_lh;
    msg.Debugging()<<"*****Build_Polarisations: -("<<Mom->type<<") "<<momcount<<endl;
    momcount++;
    Momlist.push_back(*Mom);

    //Polarisation +1
    Mom->arg[0] = momcount;
    Mom->type=mt::p_lp;
    msg.Debugging()<<"*****Build_Polarisations: +("<<Mom->type<<") "<<momcount<<endl;
    momcount++;
    Momlist.push_back(*Mom);
    if(momindex<nvec&&AMATOOLS::IsZero(Mass))  return momcount;

    //Polarisation longitudinal
    Mom->arg[0] = momcount;
    Mom->type=mt::p_l;
    msg.Debugging()<<"*****Build_Polarisations: l "<<momcount<<endl;
    momcount++;
    Momlist.push_back(*Mom);
    }
  
  if(AMATOOLS::IsZero(Mass)){ 
    if (Get_Pol_Number(momindex,mt::p_si,0,1)==-1 && rpa.me.CutScheme()!=1 ) {
      Mom->arg[0] = momcount;
      Mom->type=mt::p_si;
      msg.Debugging()<<"*****Build_Polarisations: zero mass s "<<momcount<<endl;
      momcount++;
      //cout<<"******Build_Polarisations: "<<momcount<<endl;
      Momlist.push_back(*Mom);
    }
    return momcount;
  }
  if (Get_Pol_Number(momindex,mt::p_s,Mass,1)==-1) {
    //Polarisation scalar
    Mom->arg[0] = momcount;
    Mom->type=mt::p_s;
    msg.Debugging()<<"*****Build_Polarisations: s "<<momcount<<" mass:"<<Mass<<endl;
    momcount++;
    Momlist.push_back(*Mom);
  }
 
  return momcount;
}

double Basic_Sfuncs::N(int i,int j)
{
  return 1./(sqrt(2.)*abs(S0(i,j)));
}

int Basic_Sfuncs::epsilon(int i,int j,int k,int l)
{
  short int i1,j1,k1,l1,dummy,perm;
  perm = 0;
  i1=i;j1=j;k1=k;l1=l;
  for (;;) {
     if (i1>j1) {perm++;dummy = j1;j1=i1;i1=dummy;}
     if (i1==j1) return 0;
     if (j1>k1) {perm++;dummy = k1;k1=j1;j1=dummy;}
     if (j1==k1) return 0;
     if (k1>l1) {perm++;dummy = l1;l1=k1;k1=dummy;}
     if (k1==l1) return 0;
     if ((i1<j1)&&(j1<k1)&&(k1<l1)) break;
  }
  if (perm%2==0) return 1;
  return -1; 
}

void Basic_Sfuncs::Calc_Momlist()
{
  double ps,pt;
  Vec4D mom,vh1,vh2;
  for(short int j=0;j<Momlist.size();j++){
    Momlist[j].mom_img=Vec4D(0.,0.,0.,0.);
    switch(Momlist[j].type){
    case mt::p_none : break;
    case mt::mom :   Momlist[j].mom = p[j];break;
    case mt::prop : 
      Momlist[j].mom = Vec4D(0.,0.,0.,0.);
      for (short int i=1;i<Momlist[j].argnum;i++) {
	Momlist[j].mom += b[Momlist[j].arg[i]]*p[Momlist[j].arg[i]];
      }
      break;
    case mt::cmprop :
      Momlist[j].mom = Vec4D(p[0][0]+p[1][0],0.,0.,0.);
      break;
    case mt::p_m :
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      pt=sqrt(sqr(mom[1])+sqr(mom[2]));
      if(!AMATOOLS::IsZero(pt)){
	Momlist[j].mom = sqrt(.5)*Vec4D(0.,mom[1]*mom[3]/ps/pt,
					mom[2]*mom[3]/ps/pt,-pt/ps);
	Momlist[j].mom_img = -sqrt(.5)*Vec4D(.0,-mom[2]/pt,mom[1]/pt,0.);
      }
      else {
	//if(mom[3]>0) Momlist[j].mom = sqrt(.5)*Vec4D(0.,1.,0.,0.);
	//else Momlist[j].mom = sqrt(.5)*Vec4D(0.,-1.,0.,0.);
	//Momlist[j].mom_img = sqrt(.5)*Vec4D(0.,0.,1.,0.);
	Momlist[j].mom = sqrt(.5)*Vec4D(0.,1.,0.,0.);
	if(mom[3]>0)Momlist[j].mom_img = -sqrt(.5)*Vec4D(0.,0.,1.,0.);
	else Momlist[j].mom_img = -sqrt(.5)*Vec4D(0.,0.,-1.,0.);
       }
      Momlist[j+1].mom = Momlist[j].mom;
      Momlist[j+1].mom_img = (-1.)*Momlist[j].mom_img;
      j++;
      break;

    case mt::p_l0 :
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      pt=sqrt(sqr(mom[1])+sqr(mom[2]));
      if(!AMATOOLS::IsZero(pt)){
	vh1 = Vec4D(0.,mom[1]*mom[3]/ps/pt,mom[2]*mom[3]/ps/pt,-pt/ps);
	vh2 = Vec4D(0.,-mom[2]/pt,mom[1]/pt,0.);
      }
      else {
	if(mom[3]>0) vh1 = Vec4D(0.,1.,0.,0.);
	else vh1 = Vec4D(0.,-1.,0.,0.);
	vh2 = Vec4D(0.,0.,1.,0.);
       }
      Momlist[j].mom = ::cos(Momlist[j].angle)*vh1+::sin(Momlist[j].angle)*vh2;
      Momlist[j+1].mom = ::sin(Momlist[j].angle)*vh1-::cos(Momlist[j].angle)*vh2;
      j++;
      break;
 
    case mt::p_lh :
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      pt=sqrt(sqr(mom[1])+sqr(mom[2]));
      if(!AMATOOLS::IsZero(pt)){
	vh1 = Vec4D(0.,mom[1]*mom[3]/ps/pt,mom[2]*mom[3]/ps/pt,-pt/ps);
	vh2 = Vec4D(0.,-mom[2]/pt,mom[1]/pt,0.);
	if(mom[1]+mom[2]<0)vh2=-1.*vh2;
     }
      else {
	vh1 = sqrt(.5)*Vec4D(0.,1.,-1.,0.);
	//if(mom[3]>0) vh1 = Vec4D(0.,1.,0.,0.);
	//else vh1 = Vec4D(0.,-1.,0.,0.);
	vh2 = sqrt(.5)*Vec4D(0.,1.,1.,0.);
      }
      //Momlist[j].mom = sqrt(.5)*(vh1+vh2);
      //Momlist[j+1].mom = sqrt(.5)*(vh1-vh2);      
      Momlist[j].mom = vh1;
      Momlist[j+1].mom = vh2;      
      j++;
      break;

    case mt::p_l :     
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      if(!AMATOOLS::IsZero(ps)){
	Momlist[j].mom=1./sqrt(abs(sqr(mom[0])-sqr(ps)))*Vec4D(ps,mom[0]*mom[1]/ps,
							       mom[0]*mom[2]/ps,
							       mom[0]*mom[3]/ps);
      }else
	Momlist[j].mom=Vec4D(0.,0.,0.,1.);
      if (mom.Abs2()<0.) {
	Momlist[j].mom_img = Momlist[j].mom;
	Momlist[j].mom = Vec4D(0.,0.,0.,0.);
      }
      break;
    case mt::p_s :
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=mom.Abs2();
 
      if(AMATOOLS::IsZero(Momlist[j].mass))Momlist[j].mom=Vec4D(0.,0.,0.,0.);
      else {
	  Complex help=sqrt((Complex(1.,0.)-Momlist[j].cplxmass2/ps)/Momlist[j].cplxmass2);
	  Momlist[j].mom=real(help)*mom;
	  Momlist[j].mom_img=imag(help)*mom;
	  //cout<<"Pol 's'=("<<j<<") "<<Momlist[j].mom<<","<<Momlist[j].mom_img<<endl;
	  //cout<<"complex Mass2:"<<Momlist[j].cplxmass2<<help<<" ps: "<<ps<<" Mom:"<<mom<<endl;
      }
      break;
    case mt::p_si :
      mom=Momlist[Momlist[j].arg[1]].mom;
      Complex help=csqrt(-1./mom.Abs2());
      //cout<<help<<endl;
      Momlist[j].mom =    real(help)*mom;
      Momlist[j].mom_img= imag(help)*mom;
    } 
    //cout<<"MOM("<<j<<"): "<<Momlist[j].mom<<Momlist[j].mom_img<<endl;
   
  }
}

void Basic_Sfuncs::ResetS_GT(double theta)
{
  double sq05=sqrt(.5);
  double ps,pt;
  double s,c,s0,c0,sf,cf;
  Momfunc* m;
  Vec4D mom;
  for(short int j=0;j<Momlist.size();j++){
    if(Momlist[j].type==mt::p_m&&AMATOOLS::IsZero(Momlist[j].mass)){
      mom=Momlist[Momlist[j].arg[1]].mom;
      ps=sqrt(sqr(mom[1])+sqr(mom[2])+sqr(mom[3]));
      pt=sqrt(sqr(mom[1])+sqr(mom[2]));
      s=sqrt(.5*(1-mom[3]/ps));
      c=sqrt(.5*(1+mom[3]/ps));
      s0=::sin(theta*0.5);
      c0=::cos(theta*0.5);
      if(!AMATOOLS::IsZero(pt)){
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
    if(AMATOOLS::IsZero(Momlist[j].mass))if(Momlist[j].type==mt::p_m||Momlist[j].type==mt::p_p){
      m=&Momlist[j];
      //cout<<"recalc "<<j<<endl;
      switch(k0_n){
      case 1:
	if(!iscplx(j)) 
	  _eta[j] = csqrt(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*sq05));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*sq05),
				    2.*(m->mom_img[0]-(m->mom_img[2]+m->mom_img[3])*sq05)));
	break;
      case 2:
	if(!iscplx(j)) 
	  _eta[j] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*sq05));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*sq05),
				    2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[2])*sq05)));
	break;
      default:
	if(!iscplx(j)) 
	  _eta[j] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*sq05));
	else _eta[j] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*sq05),
				    2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[3])*sq05)));
      }
    }
  }
}

int Basic_Sfuncs::setS(Vec4D* _p)
{
  //  PROFILE_HERE;
  PROFILE_LOCAL("int Basic_Sfuncs::setS(Vec4D* _p)");
  // _eta's and _mu's precalc
  
  double sq05=sqrt(.5);
  p = _p;
  Calc_Momlist();

  Momfunc* m;
  
  int etachk=1;

  for(short int i=0;i<Momlist.size();i++) {
    m=&Momlist[i];

    switch(k0_n){
    case 1:
      if(!iscplx(i)) 
	_eta[i] = csqrt(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*sq05));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]-(m->mom[2]+m->mom[3])*sq05),
				  2.*(m->mom_img[0]-(m->mom_img[2]+m->mom_img[3])*sq05)));
      break;
    case 2:
      if(!iscplx(i)) 
	_eta[i] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*sq05));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[2])*sq05),
				  2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[2])*sq05)));
      break;
    default:
      if(!iscplx(i)) 
	_eta[i] = csqrt(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*sq05));
      else _eta[i] = sqrt(Complex(2.*(m->mom[0]-(m->mom[1]+m->mom[3])*sq05),
				  2.*(m->mom_img[0]-(m->mom_img[1]+m->mom_img[3])*sq05)));
    }

    if(AMATOOLS::IsZero(_eta[i]))etachk=0;
    if (i<Nmom) {
      _mu[i]  =Momlist[i].mass/_eta[i];
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
	    case mt::p_l:
	    case mt::p_si:
		_mu[i]  = (csqrt((m->mom).Abs2())+csqrt(-(m->mom_img).Abs2()))/_eta[i];
		break;
	    case mt::p_s:
		_mu[i] = sqrt(Complex((m->mom).Abs2()-(m->mom_img).Abs2(),
				      2 * m->mom * m->mom_img))/_eta[i];
		break;
	    default: _mu[i]  = csqrt((m->mom).Abs2())/_eta[i];
	}
	/*_mu[i]  = csqrt((m->mom).Abs2())/_eta[i];
      if(m->type==mt::p_p||m->type==mt::p_m) _mu[i] = Complex(0.,0.);
      if(m->type==mt::p_l||m->type==mt::p_s) _mu[i]+= csqrt(-(m->mom_img).Abs2())/_eta[i];*/
    }
  }
  for(short int i=0;i<momcount;i++)
    for(short int j=0;j<momcount;j++)
      calc_st[i][j]=0;
  return etachk;
}


void Basic_Sfuncs::calcS(int i, int j)
{
  //  PROFILE_HERE;
  double sq05=sqrt(.5);
  Momfunc* m = &Momlist[i];
  Momfunc* m1 = &Momlist[j];
  if (!(AMATOOLS::IsZero((m->mom+(-1.)*m1->mom).Abs2()) &&
	AMATOOLS::IsZero(m->mom[0]-m1->mom[0]))) {
    Complex A= _eta[j]/_eta[i];

    switch(k0_n){
    case 1:
      _S0[i][j] = Complex(m->mom[1],(m->mom[2]-m->mom[3])*sq05)*A; 
      _S1[i][j] = Complex(-m->mom[1],(m->mom[2]-m->mom[3])*sq05)*A;
      if (iscplx(i)){
	_S0[i][j] += Complex(-(m->mom_img[2]-m->mom_img[3])*sq05,m->mom_img[1])*A;
	_S1[i][j] += Complex(-(m->mom_img[2]-m->mom_img[3])*sq05,-m->mom_img[1])*A;
      }      
      _S0[i][j] -= Complex(m1->mom[1],(m1->mom[2]-m1->mom[3])*sq05)/A;
      _S1[i][j] -= Complex(-m1->mom[1],(m1->mom[2]-m1->mom[3])*sq05)/A;
      if (iscplx(j)){
	_S0[i][j] -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*sq05,m1->mom_img[1])/A;
	_S1[i][j] -= Complex(-(m1->mom_img[2]-m1->mom_img[3])*sq05,-m1->mom_img[1])/A;
      }
      break;
    case 2:
      _S0[i][j] = Complex(m->mom[3],(m->mom[1]-m->mom[2])*sq05)*A; 
      _S1[i][j] = Complex(-m->mom[3],(m->mom[1]-m->mom[2])*sq05)*A;
      if (iscplx(i)){
	_S0[i][j] += Complex(-(m->mom_img[1]-m->mom_img[2])*sq05,m->mom_img[3])*A;
	_S1[i][j] += Complex(-(m->mom_img[1]-m->mom_img[2])*sq05,-m->mom_img[3])*A;
      }
      _S0[i][j] -= Complex(m1->mom[3],(m1->mom[1]-m1->mom[2])*sq05)/A;
      _S1[i][j] -= Complex(-m1->mom[3],(m1->mom[1]-m1->mom[2])*sq05)/A;
      if (iscplx(j)){
	_S0[i][j] -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*sq05,m1->mom_img[3])/A;
	_S1[i][j] -= Complex(-(m1->mom_img[1]-m1->mom_img[2])*sq05,-m1->mom_img[3])/A;
      }
      break;
    default:
      _S0[i][j] = Complex(m->mom[2],(m->mom[3]-m->mom[1])*sq05)*A; 
      _S1[i][j] = Complex(-m->mom[2],(m->mom[3]-m->mom[1])*sq05)*A;
      if (iscplx(i)){
	_S0[i][j] += Complex(-(m->mom_img[3]-m->mom_img[1])*sq05,m->mom_img[2])*A;
	_S1[i][j] += Complex(-(m->mom_img[3]-m->mom_img[1])*sq05,-m->mom_img[2])*A;
      }
      _S0[i][j] -= Complex(m1->mom[2],(m1->mom[3]-m1->mom[1])*sq05)/A;
      _S1[i][j] -= Complex(-m1->mom[2],(m1->mom[3]-m1->mom[1])*sq05)/A;
      if (iscplx(j)){
	_S0[i][j] -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*sq05,m1->mom_img[2])/A;
	_S1[i][j] -= Complex(-(m1->mom_img[3]-m1->mom_img[1])*sq05,-m1->mom_img[2])/A;
      }
    }

    _S0[j][i] = - _S0[i][j];
    _S1[j][i] = - _S1[i][j];
  }
  else {	
    //cout<<"Zero: "<<i<<";"<<j<<endl;
    _S0[i][j] = Complex(0.,0.);
    _S0[j][i] = Complex(0.,0.);
    _S1[i][j] = Complex(0.,0.);
    _S1[j][i] = Complex(0.,0.);
  }
  //cout<<"S0: "<<_S0[m->arg[0]][m1->arg[0]]<<endl;
  calc_st[i][j]=calc_st[j][i]=1;
}

void Basic_Sfuncs::set_k0(int i)
{
  k0_n=i;
  //i=0: k0=Vec4D(1.,sqrt(.5),0.,sqrt(.5));
  //     k1=Vec4D(0.,0.,1.,0.);

  //i=1: k0=Vec4D(1.,0.,sqrt(.5),sqrt(.5));
  //     k1=Vec4D(0.,1.,0.,0.);

  //i=2: k0=Vec4D(1.,sqrt(.5),sqrt(.5),0.);
  //     k1=Vec4D(0.,0.,0.,1.);
}


