#include "Flavour.H"
#include "Run_Parameter.H"
#include "Decay.H"
#include "Interaction_Model_Base.H"
#include "Message.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace std;

//const double min_BR = 0.05;
const double min_BR = 0;

int Decay_Handler::CheckInVertex(Flavour fl)
{
  Vertex* v = mo->GetVertex();
  for (int i=0;i<v->MaxNumber();i++) {
    if ((*v)[i]->on) {
      for (short int j=0;j<3;j++) {
	if ((*v)[i]->in[j].Kfcode()==fl.Kfcode()) return 1;
      }
    }
  }
  return 0;
}

int Decay_Handler::CheckVertex(Flavour flin,Flavour flout1,Flavour flout2)
{
  Vertex* v = mo->GetVertex();
  for (int i=0;i<v->MaxNumber();i++) {
    if ((*v)[i]->on) {
	if ((*v)[i]->in[0]==flin &&
	    (*v)[i]->in[1]==flout1 &&
	    (*v)[i]->in[2]==flout2) return 1;
    }
  }
  return 0;
}

void Decay_Handler::Print()
{
  DecayTable* dt = dtfirst;
  msg.Out()<<"Decay Table: "<<endl;
  while (dt) {
    msg.Out()<<dt->fl[0]<<"(Width="<<dt->width<<" GeV) decays into: "<<endl;
    DecayTable* dtm = dt->More;
    while (dtm) {
      msg.Out()<<"       "<<dtm->BR*100.<<"%    ---> "<<dtm->fl[1]<<";"<<dtm->fl[2]<<endl;
      dtm = dtm->Next;
    }
    dt = dt->Next;
  }
}

void Decay_Handler::Add(DecayTable* dt,Flavour flin, 
			Flavour flout1,Flavour flout2)
{
  int sw = 0;
  sw = CheckVertex(flin,flout1,flout2);
  if (sw==0) {
    sw = CheckVertex(flin,flout2,flout1);
    if (sw==1) {
      Flavour flh = flout2;
      flout2      = flout1;
      flout1      = flh;
    }
  }
  if (sw==1) {
    DecayTable* dtnew = new DecayTable;
    dtnew->fn    = 3;
    dtnew->fl    = new Flavour[3];
    dtnew->fl[0] = flin;
    dtnew->fl[1] = flout1;
    dtnew->fl[2] = flout2;
    dtnew->width = 0.;
    dtnew->BR    = 1.;
    dtnew->Next  = 0;
    dtnew->More  = 0;

    if (dt->More==0) dt->More = dtnew;
    else {
      DecayTable* dtit = dt->More;
      while(dtit) {
	if (dtit->Next==0) {
	  dtit->Next = dtnew;
	  break;
	}
	dtit = dtit->Next;
      }
    }
  }
}

void Decay_Handler::FindDecayProducts2()
{
  DecayTable* dt = dtfirst;
  while (dt) {
    Flavour infl = dt->fl[0];
    int charge = infl.IntCharge();
    int spin   = infl.IntSpin()%2;
    Flavour flmatch1,flmatch2;
    for (int i=1;i<MAX_PARTICLES;i++) {
      if (particles[i].kfc==kf::none) break;
      if (particles[i].on==1) {
	if (CheckInVertex(Flavour(particles[i].kfc))==1) {
	  Flavour flout1 = Flavour(particles[i].kfc);
	  for (int j=i;j<MAX_PARTICLES;j++) {
	    if (particles[j].kfc==kf::none) break;
	    if (particles[j].on==1) {
	      Flavour flout2 = Flavour(particles[j].kfc);
	      //mass check later
	      if (flout1!=infl && flout2!=infl &&
		  flout1!=infl.Bar() && flout2!=infl.Bar() &&
		  flout1.Mass()<infl.Mass() &&
		  flout2.Mass()<infl.Mass() ) {
		//spin check
		if (spin == (flout1.IntSpin()+flout2.IntSpin())%2) {
		  //charge check
		  if (charge==flout1.IntCharge()+flout2.IntCharge())
		    Add(dt,infl,flout1,flout2);
		  if (charge==-flout1.IntCharge()+flout2.IntCharge() &&
		      flout1!=flout1.Bar()) 
		    Add(dt,infl,flout1.Bar(),flout2);
		  if (charge==flout1.IntCharge()-flout2.IntCharge() &&
		      flout2!=flout2.Bar() &&
		    flout1!=flout2)
		    Add(dt,infl,flout1,flout2.Bar());
		  if (charge==-flout1.IntCharge()-flout2.IntCharge() &&
		      flout1!=flout1.Bar() &&
		      flout2!=flout2.Bar())
		    Add(dt,infl,flout1.Bar(),flout2.Bar());
		}
	      }
	    }
	  }
	}
      }
    }
    dt = dt->Next;
  }
}

void Decay_Handler::FindDecayProducts3()
{
  //mass check + further decay
  DecayTable* dt = dtfirst;
  while (dt) {
    DecayTable* dtm = dt->More;
    while (dtm) {
      double sumM = 0.;
      for (short int i=1;i<dtm->fn;i++) sumM += (dtm->fl[i]).Mass();
      if (sumM>(dtm->fl[0]).Mass()) {
	NextDecay(dtm);
      }
      dtm = dtm->Next;
    }
    dt = dt->Next;
  }
}

void Decay_Handler::NextDecay(DecayTable* Ndt)
{  
  //look for further decay
  for (short int i=1;i<Ndt->fn;i++) {
    // find in the list
    DecayTable* dt = dtfirst;
    while (dt) {
      if (dt->fl[0]==Ndt->fl[i] || dt->fl[0]==(Ndt->fl[i]).Bar()) {
	//Loop over further Decays
	DecayTable* dtm = dt->More;
	while (dtm) {
	  double sumM = 0.;
	  for (short int j=1;j<dtm->fn;j++) sumM += (dtm->fl[j]).Mass();
	  // no further decay
	  if (sumM<(dtm->fl[0]).Mass()) {
	    DecayTable* fdt = new DecayTable;
	    fdt->fn = Ndt->fn-1 + dtm->fn-1; 
	    fdt->fl = new Flavour[fdt->fn];		
	    //old flav part 1
	    for (short int j=0;j<i;j++) fdt->fl[j] = Ndt->fl[j];
	    //new flav
	    for (short int j=i;j<i+dtm->fn-1;j++) {
	      if (dt->fl[0]==Ndt->fl[i]) fdt->fl[j] = dtm->fl[j-i+1];
	      else fdt->fl[j] = (dtm->fl[j-i+1]).Bar();
	    }
	    //old flav part 2
	    for (short int j=i+dtm->fn-1;j<fdt->fn;j++) 
	      fdt->fl[j] = Ndt->fl[j-dtm->fn+1+1];
	    
	    fdt->width = 0.;
	    fdt->BR    = 1.;
	    fdt->Next  = 0;
	    fdt->More  = 0;
	    //adding
	    if (Ndt->More==0) Ndt->More = fdt;
	    else {
	      DecayTable* dtit = Ndt->More;
	      while (dtit) {
		if (dtit->Next==0) {
		  dtit->Next = fdt;
		  break;
		}
		dtit = dtit->Next;
	      }
	    }
	  }
	  dtm = dtm->Next;
	}
      }
      dt = dt->Next;
    }
  }    
}

void Decay_Handler::FindUnstable()
{
  int count=0;
  dtfirst = 0;

  for (int i=1;i<MAX_PARTICLES;i++) {
    if (particles[i].kfc==kf::none) break;
    if (particles[i].stbl==0 && particles[i].on==1) {
      if (CheckInVertex(Flavour(particles[i].kfc))==1) {
	DecayTable* dt = new DecayTable;
	dt->fn    = 1;
	dt->fl    = new Flavour;
	dt->fl[0] = Flavour(particles[i].kfc);
	dt->width = 0.;
	dt->BR    = 1.;
	dt->Next  = 0;
	dt->More  = 0;

	if (dtfirst==0) dtfirst = dt;
	else {
	  DecayTable* dtit = dtfirst;
	  DecayTable* prev = 0;
	  while (dtit) {
	    if (dtit->fl[0].Mass()>dt->fl[0].Mass()) break;
	    prev = dtit;
	    dtit = dtit->Next;
	  }
	  if (prev==0) {
	    dt->Next = dtfirst;
	    dtfirst = dt;
	  }
	  else {
	    prev->Next = dt;
	    dt->Next   = dtit;
	  }
	}

	count++;
      }
    }
  }
  FindDecayProducts2();
  FindDecayProducts3();
}


void Decay_Handler::RecCalc(Topology* top,DecayTable* dt,int sw_2)
{
  DecayTable* dtit = dt;
  while (dtit) {
    int sw = 1;
    if (dtit->More==0) {
      Single_Process* pro;
      //pro = new Single_Process(1,dtit->fn-1,dtit->fl,0,1);
      //pro->InitDecay(top);
      /*
	if (!pro->InitAmplitude(top,moms,results,links)) 
 	msg.Error()<<"Error in InitAmplitude !"<<endl;
      */      
      pro->SetUpIntegrator();
      if (!pro->CalculateTotalXSec()) msg.Debugging()<<"pro->CalculateTotalXSec left with zero "<<endl; 
	dtit->width = (pro->FSRIntegrator())->Result()/(pro->FSRIntegrator())->N();
      delete pro;      
    }
    else {
      if (sw_2==1) {
	//possible kicker
	int hit = 0;
	double sumM = 0.;
	for (short int i=1;i<dtit->fn;i++) sumM += (dtit->fl[i]).Mass();
	for (short int i=1;i<dtit->fn;i++) {
	  Flavour fl = dtit->fl[i];
	  double value = fl.Mass()*fl.Width()/
	    (sqr(fl.Mass())-sqr((dtit->fl[0]).Mass()-sumM+fl.Mass()));
	  if (value>min_BR) {
	    hit = 1;
	    break;
	  }
	}
	if (hit==0) sw = 0;
      }
      dtit->width = 0;
      if (sw==1) {
	msg.Out()<<"-----------------------------------------"<<endl;
	RecCalc(top,dtit->More,sw_2);
	DecayTable* dtm = dtit->More;
	while (dtm) {
	  dtit->width += dtm->width;
	  dtm = dtm->Next;
	}
	msg.Out()<<"----------------SUM----------------------"<<endl;
      }
    }
    msg.Out()<<dtit->fl[0]<<" ---> "<<dtit->fl[1];
    for (short int i=2;i<dtit->fn;i++) msg.Out()<<";"<<dtit->fl[i];
    if (sw==1) msg.Out()<<",  Width: "<<dtit->width<<endl;
          else msg.Out()<<",  Kicked ! "<<endl;
    dtit=dtit->Next;
  }
}

void Decay_Handler::RecBR(DecayTable* dt,double refwidth)
{
  DecayTable* dtit = dt;
  while (dtit) {
    if (dtit->More==0) dtit->BR = dtit->width/refwidth;
    else {
      RecBR(dtit->More,refwidth);
      //Sum over Ratios
      dtit->BR = 0;
      DecayTable* dtm = dtit->More;
      while (dtm) {
	dtit->BR += dtm->BR;
	dtm = dtm->Next;
      }
    }
    dtit = dtit->Next;
  }
}

void Decay_Handler::BranchingRatios()
{
  DecayTable* dtit = dtfirst;
  while (dtit) {
    RecBR(dtit->More,dtit->width);
    dtit = dtit->Next;
  }
}

void Decay_Handler::Calculate(Topology* top)
{
  if (msg.Tracking()) Print();
  DecayTable* dtit = dtfirst;
  while (dtit) {
    int sw_2 = 0;
    DecayTable* dtm = dtit->More;
    dtit->width = 0;
    while (dtm) {
      double sumM = 0.;
      for (short int i=1;i<dtm->fn;i++) sumM += (dtm->fl[i]).Mass();
      if (sumM<(dtm->fl[0]).Mass()) {
	sw_2 = 1;
	break;
      }
      dtit->width += dtm->width;
      dtm = dtm->Next;
    }
    RecCalc(top,dtit->More,sw_2);
    dtit->width = 0;
    dtm = dtit->More;
    while (dtm) {
      dtit->width += dtm->width;
      dtm = dtm->Next;
    }
    msg.Out()<<"-------------------------------------------------------------------"<<endl;
    msg.Out()<<dtit->fl[0]<<" ---> All,  Width: "<<dtit->width<<endl;
    //if (IsZero((dtit->fl[0]).width())) 
    (dtit->fl[0]).SetWidth(dtit->width);
    dtit=dtit->Next;
  }
  BranchingRatios();
  Print();
}
