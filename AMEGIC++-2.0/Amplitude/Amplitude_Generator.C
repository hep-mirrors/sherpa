//#include <iomanip>
#include "Amplitude_Generator.H"
#include "Amplitude_Manipulator.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MathTools.H"

//Do not use this here !!! 
//THis is only for test purposes !!!
//#define _USE_MPI_

#ifdef _USE_MPI_
#include <mpi++.h>
#include <algorithm>

int AMEGIC::Amplitude_Generator::NMAX=15000;
//int AMEGIC::Amplitude_Generator::NMAX=500;

MPI::Datatype   mpi_lf_type;
MPI::Datatype   mpi_cf_type;
MPI::Datatype   mpi_sv_type;
MPI::Datatype   mpi_point_type;
#endif

namespace AMEGIC {
  class Compare_Pre_Amplitudes {
  public:
    int operator()(const AMEGIC::Pre_Amplitude & a, const AMEGIC::Pre_Amplitude & b) {
      if (a.perm<b.perm) return 1;
      return 0;
    }
  };
}

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace std;

Amplitude_Generator::Amplitude_Generator(int _no,Flavour* _fl,int* _b,
					 Interaction_Model_Base * _model,Topology* _top,
					 int _nQCD,int _nEW,
					 Basic_Sfuncs* _BS,String_Handler* _shand) 
  : N(_no), fl(_fl), b(_b), p_model(_model), top(_top), nEW(_nEW), nQCD(_nQCD),
    BS(_BS), shand(_shand), s_buffer(0)
{
  single_top = top->Get(N-2);
  
  // 2 incoming
  prenum  = 1000;
  prea    = new Pre_Amplitude[prenum];
  short int i;
  for(i=0;i<prenum;i++) prea[i].p = new Point[single_top->depth];

  // fill hash table
  
  Vertex* v = p_model->GetVertex();

  for (int i=0;i<v->MaxNumber();++i) {
    if ((*v)[i]->on) {
      v_table[(*v)[i]->in[0]].push_back((*v)[i]);
    }
  }

#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();
  if (rank==0) {
    for (int i=0; i<size; ++i) {
      r_buffers.push_back(Prea_Buffer());
      msg.Tracking()<<" create Buffer for "<<i<<" size "<<NMAX<<endl;
      r_buffers[i].buff     = new MPI_Point[NMAX];
      r_buffers[i].counters = new int[2*NMAX]; // too big!!!
      r_buffers[i].nbuff = 0;
    }
  }
  else {
    r_buffers.push_back(Prea_Buffer());
    msg.Tracking()<<" create Buffer for Slave "<<rank<<" size "<<NMAX<<endl;
    r_buffers[0].buff     = new MPI_Point[NMAX];
    r_buffers[0].counters = new int[2*NMAX]; // too big!!!
    r_buffers[0].nbuff = 0;
  }

  CommitMPITypes();
#endif
  
}

Amplitude_Generator::~Amplitude_Generator() 
{
#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();
  if (rank==0) {
    for (int i=0; i<size; ++i) {
      delete [] r_buffers[i].buff;
    }
  }
#endif

  for(short int i=0;i<prenum;i++) delete[] prea[i].p;
  delete[] prea;
}

void Amplitude_Generator::Set_End(Point* p,int* &perm,int& pnum)
{
  p->b     = 0;
  p->fl    = Flavour(kf::none);
  if ((p->left==0) && (p->right==0)) {
    p->number = *perm;
    p->fl = fl[*perm];
    p->b  = b[*perm];
    if (p->fl.IsBoson()) {
      p->Lorentz->type = lf::Pol;
      p->Lorentz->SetParticleArg(0);
    }
    else {
      p->Lorentz->type = lf::None;
      p->Lorentz->SetParticleArg();
    }

    perm++;
    return;
  }
  p->number = pnum;
  pnum++;
  Set_End(p->left,perm,pnum);
  Set_End(p->right,perm,pnum);
}

void Amplitude_Generator::Next_P(Point* p,Point* &hit)
{
  if (hit) return;
  if (p==0) return;
  if ((p->left!=0) && (p->right!=0)) {
    if ((p->left->fl==Flavour(kf::none)) || (p->right->fl==Flavour(kf::none))) {
      hit = p;
      return;
    }
  }
  Next_P(p->left,hit); 
  Next_P(p->right,hit);
}


void Amplitude_Generator::Print_P(Point* p)
{
  if (!(AORGTOOLS::rpa.gen.Debugging())) return;
  if ((p->left==0) && (p->right==0)) {
    AORGTOOLS::msg.Out()<<"EndPoint : "<<p->fl<<"("<<p->b<<")"<<endl;
    return;
  }
  AORGTOOLS::msg.Out()<<"left : ";
  Print_P(p->left);
  AORGTOOLS::msg.Out()<<"right : ";
  Print_P(p->right);
  if(p->middle){
    AORGTOOLS::msg.Out()<<" middle : ";
    Print_P(p->middle);
  }
}

int Amplitude_Generator::Match_Vertex(Single_Vertex* v,Flavour* flav,Complex* cpl)
{
  if (flav[0] == v->in[0]) {
    short int hit = 1;
    if (flav[1] != Flavour(kf::none)) {if (flav[1] != v->in[1]) hit = 0;}
    else {flav[1] = v->in[1];}
    if (flav[2] != Flavour(kf::none)) {if (flav[2] != v->in[2]) hit = 0;}
    else {flav[2] = v->in[2];}
    if (hit==1) {
      for (short int j=0;j<4;j++) cpl[j] = v->cpl[j];
      return 1;
    }
  }
  return 0;
}

int Amplitude_Generator::Check_End(Point* p,Flavour infl) 
{
  if (p==0) return 1;
  if (p->left==0) return 1;
  if (((p->left->fl)!=Flavour(kf::none)) && ((p->right->fl)!=Flavour(kf::none))) { 
    //AORGTOOLS::msg.Out()<<"It's a Check!!!"<<endl;
    short int j,k;
    Flavour flav[3];
    Flavour s_flav[3];
    Complex cpl[4];
    for (j=0;j<4;j++) cpl[j] = Complex(0.,0.);

    //if (infl.Majorana()) AORGTOOLS::msg.Out()<<"Incoming majorana -> "<<p->left->fl<<";"<<p->right->fl<<endl;

    flav[0] = infl;
    flav[1] = p->left->fl;
    flav[2] = p->right->fl;

    if (infl.Majorana()) {
      if (p->left->fl.IsFermion()) {
	if (p->b*p->left->b==1) flav[1] = flav[1].Bar(); 
	if (p->right->b == -1)  flav[2] = flav[2].Bar();
      }
      else {
	if (p->b*p->right->b==1) flav[2] = flav[2].Bar(); 
	if (p->left->b == -1)    flav[1] = flav[1].Bar();
      }	  
      //AORGTOOLS::msg.Out()<<"After Change: "<<flav[1]<<";"<<p->b<<";"<<p->left->b<<endl;
      //AORGTOOLS::msg.Out()<<"After Change: "<<flav[2]<<";"<<p->b<<";"<<p->right->b<<endl;
    }
    else {
      if (infl.IsBoson()) {
	if (p->left->b == -1) flav[1] = flav[1].Bar();
	if (p->right->b == -1) flav[2] = flav[2].Bar();
      }
      else {
	if (infl.IsAnti()) {
	  if (p->b*p->left->b==1) flav[1] = flav[1].Bar(); 
	  if (p->right->b == -1)  flav[2] = flav[2].Bar();
	}
	else {
	  if (p->b*p->right->b==1) flav[2] = flav[2].Bar(); 
	  if (p->left->b == -1)    flav[1] = flav[1].Bar();
	}
      }
    }

    s_flav[0]=flav[0];
    s_flav[1]=flav[1];
    s_flav[2]=flav[2];

    // look for flav[0]
    Vertex_List & vl= v_table[flav[0]];

    for (j=0;j<vl.size();j++) {
      if (Match_Vertex(vl[j],flav,cpl)) {
	for (k=0;k<4;k++) p->cpl[k] = cpl[k];
	p->v = vl[j];
	*(p->Color)   = *(vl[j]->Color);
	*(p->Lorentz) = *(vl[j]->Lorentz);
	  
	return 1;
      }
    }
  }
  else return 1;
  return 0;
}

void Amplitude_Generator::Set_Props(Point* pl,int dep,Single_Amplitude* &first,int* perm,int topcount, int permcount)
{
  int ap = 0;
  int lanz = 1;

  Point* preah;
  preah = new Point[dep];

  Point* p;

  short int i,k;
  int help = 0;
  top->Copy(pl,prea[0].p,help);

  prea[0].on = 1;

  int sw1;
  Flavour flav[3];
  Flavour s_flav[3];
  Complex cpl[4];

  int first_try = 1;
  
  for (;;) {    
    sw1 = 1;
    if (prea[ap].on) {
      p = 0;
      Next_P(&prea[ap].p[0],p);
      if (p==0) {
	if (first_try==1) p = &prea[ap].p[0];
	else sw1 = 0;
      }
      first_try = 0;
    }
    else sw1 = 0;
    if (sw1) {
      flav[0] = p->fl;
      flav[1] = p->left->fl;
      flav[2] = p->right->fl;

      if (p->left->fl  == Flavour(kf::none)) p->left->b  = 0;
      if (p->right->fl == Flavour(kf::none)) p->right->b = 0;

      if (flav[0].Majorana()) {
	if (p->left->fl != Flavour(kf::none)) {
	  if (p->left->fl.IsFermion()) {
	    if (p->b*p->left->b == 1)  flav[1] = flav[1].Bar();
	  }
	}
	if (p->right->fl != Flavour(kf::none)) {
	  if (p->right->fl.IsFermion()) {
	    if (p->b*p->right->b == 1)  flav[2] = flav[2].Bar();
	  }
	}
      }
      else {
	if (flav[0].IsBoson()) {
	  if (p->left->b   == -1) flav[1] = flav[1].Bar();
	  if (p->right->b  == -1) flav[2] = flav[2].Bar();
	  if (p->left->fl  == Flavour(kf::none)) p->left->b  = -1;
	  if (p->right->fl == Flavour(kf::none)) p->right->b = -1;
	}
	else {
	  if (flav[0].IsAnti()) {
	    if (p->b*p->left->b == 1)  flav[1] = flav[1].Bar();
	    if (p->right->b     ==-1)  flav[2] = flav[2].Bar();
	    if (p->left->fl     == Flavour(kf::none)) p->left->b = p->b;
	  }
	  else {
	    if (p->b*p->right->b == 1) flav[2] = flav[2].Bar();
	    if (p->left->b       ==-1) flav[1] = flav[1].Bar();
	    if (p->right->fl     == Flavour(kf::none)) p->right->b = p->b;
	  }
	}
      }
      s_flav[0]=flav[0];
      s_flav[1]=flav[1];
      s_flav[2]=flav[2];

      // look for flav[0]
      Vertex_List & vl= v_table[flav[0]];

      for (i=0;i<vl.size();i++) {
	sw1 = 0;
	flav[0]=s_flav[0];
	flav[1]=s_flav[1];
	flav[2]=s_flav[2];

	if (Match_Vertex(vl[i],flav,cpl)) {
	  if (flav[0].Majorana()) {    
	    if (flav[1].IsFermion() && p->left->b==0)  p->left->b  = p->b;
	    if (flav[2].IsFermion() && p->right->b==0) p->right->b = p->b;
	  }
	  sw1 = Check_End(p->left,flav[1]);	    
	  if (sw1) {
	    sw1 = Check_End(p->right,flav[2]);
	  }
	}
	if (sw1) {
	  //match
	  int ll = 0;
	  top->Copy(prea[ap].p,preah,ll);
	  if (p->left->fl==Flavour(kf::none))  p->left->fl  = flav[1];
	  if (p->right->fl==Flavour(kf::none)) p->right->fl = flav[2];
	  p->v          = vl[i];
	  *(p->Color)   = *(vl[i]->Color);
	  *(p->Lorentz) = *(vl[i]->Lorentz);
	  
	  for (k=0;k<4;k++) p->cpl[k] = cpl[k];
	  ll = 0;top->Copy(prea[ap].p,prea[lanz].p,ll);
	  ll = 0;top->Copy(preah,prea[ap].p,ll);
	  prea[lanz].on = 1;
	  lanz++;
	}
      }
      prea[ap].on = 0;
    }
    if (ap==lanz-1) break;
    ap++;
  }

  for (i=0;i<lanz;i++) {
    if (prea[i].on) {
      int ll = 0;
      Point* ph;
      ph = new Point[dep];
      top->Copy(prea[i].p,ph,ll);
      prea_table.push_back(Pre_Amplitude(ph,topcount,permcount));
      //      counter_table.push_back(counters);
    }
  }

  delete[] preah;
}


void Amplitude_Generator::CreateSingleAmplitudes(Single_Amplitude * & first) {
#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();
#else
  int rank=0;
  //int size=1;
#endif


  int count=0;
  Single_Amplitude* n;
  int dep=single_top->depth;

  n = first;
  if (n) while (n->Next) n = n->Next;
  Single_Amplitude* gra;

  
  for (int i=0;i<prea_table.size();i++) {
    int sw1 = 1;
    if (AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::QCD ||
	AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::pure_QCD) {
      sw1 = 0;
      for (int j=0;j<dep;j++) {
	if (((prea_table[i].p[j].fl).IsBoson()) && 
	    (prea_table[i].p[j].fl!=Flavour(kf::gluon))) sw1++;
      }
      if (AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::pure_QCD) {
	if (sw1>0) sw1 = 0;
	else  sw1 = 1;
      }
      else {
	if (sw1>1) sw1 = 0;
	else sw1 = 1;
      }   
    }
    // test if 3-Vertex
    for (int j=0;j<dep;j++) {
      if (prea_table[i].p[j].left!=0) {
	if ((prea_table[i].p[j].v)->in[1]==(prea_table[i].p[j].v)->in[2]) {
	  if (prea_table[i].p[j].left->number>prea_table[i].p[j].right->number) {
	    sw1 = 0;
	    break;
	  }
	}
      }
    }
    //test if 4-vertex
    if (sw1) {
      for (int j=0;j<dep;j++) {
	if (prea_table[i].p[j].left!=0) {
	  if ((prea_table[i].p[j].v)->in[1]==(prea_table[i].p[j].v)->in[2]) {
	    if (prea_table[i].p[j].left->number>prea_table[i].p[j].right->number) {
	      sw1 = 0;
	      break;
	    }
	  }
	}
      }
    }

    // test if 5-vertex - Starform
    // smallest number must be on the upper side
    if (sw1) {
      int numbers[4];
      for (int j=0;j<dep;j++) {
	if (prea_table[i].p[j].left!=0) {
	  if ((prea_table[i].p[j].v)->in[1]==(prea_table[i].p[j].v)->in[2]) {
	    for (int k=0;k<4;k++) numbers[k] = -1;
	    //first gluon
	    if (prea_table[i].p[j].left->left!=0) {
	      if ((prea_table[i].p[j].left->v)->in[1]==(prea_table[i].p[j].left->v)->in[2]) {
		numbers[0] = prea_table[i].p[j].left->left->number;
		numbers[1] = prea_table[i].p[j].left->right->number;
	      }
	    }
	    //second
	    if (prea_table[i].p[j].right->left!=0) {
	      if ((prea_table[i].p[j].right->v)->in[1]==(prea_table[i].p[j].right->v)->in[2]) {
		numbers[2] = prea_table[i].p[j].right->left->number;
		numbers[3] = prea_table[i].p[j].right->right->number;
	      }
	    }
	    // check if match
	    int sw3 = 1;
	    for (int k=0;k<4;k++) {
	      if (numbers[k]== -1) {
		sw3 = 0;
		break;
	      }
	    }
	    if (sw3) {
	      if (numbers[0]>numbers[2] || numbers[0]>numbers[3]) sw1 = 0;
	    }
	  }
	}
      }
    }
    if (sw1) {
      for (int k=0;k<dep;++k) {
	if ((prea_table[i].p[k].number>99) && ((prea_table[i].p[k].fl).IsBoson())) 
	  prea_table[i].p[k].number += 100;
      }
      gra = new Single_Amplitude(prea_table[i].p,prea_table[i].top,prea_table[i].perm,b,dep,N,top,BS,fl,shand);

      count++;
      if (rank==0) {
	AORGTOOLS::msg.Out()<<"*";AORGTOOLS::msg.Out().flush();
	if (count%50==0)      AORGTOOLS::msg.Out()<<std::endl;
      }
      if (first) n->Next = gra;
      else first   = gra; 
      n = gra;
    }
  }


}


void  Amplitude_Generator::CollectPreAmplitudes() {
#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();

  int dep  = single_top->depth;

  if (size == 1) return;

  // if slave
  if (rank>0) {
    int nbuf  = prea_table.size()*dep;
    int ncbuf = 2*prea_table.size();

    msg.Tracking()<<"CPU "<<rank<<" creates MPI_Points "<<nbuf<<"("<<prea_table.size()<<")"<<endl;
    s_buffer=new MPI_Point[nbuf];
    s_counter_buffer = new int[ncbuf];

    // - fuelle buffer
    int ibuf  =0;
    int icbuf =0;
    for (int i=0; i<prea_table.size();++i) {      
      Point2MPI(prea_table[i].p,&s_buffer[ibuf]);
      ibuf+=dep;
      s_counter_buffer[2*i]   = prea_table[i].top;
      s_counter_buffer[2*i+1] = prea_table[i].perm;
    }

    // check if receive-buffer is big enough !!!
    if (nbuf<=NMAX) {
    
      // send local list
      msg.Tracking()<<"CPU"<<rank<<"sends "<<nbuf<<" points now"<<endl;
      MPI::COMM_WORLD.Send(s_buffer,nbuf,mpi_point_type,0,Mpi_Tag::complete_buffer);
      MPI::COMM_WORLD.Send(s_counter_buffer,ncbuf,MPI::INT,0,Mpi_Tag::complete_buffer);  
    }
    else {
      ibuf  =0;
      icbuf =0;
      int pnbuf=NMAX/dep;
      int pncbuf=pnbuf*2;
      pnbuf*=dep;
      for (;;) {
	msg.Tracking()<<"CPU"<<rank<<" sends only "<<pnbuf<<" points now"<<endl;
	MPI::COMM_WORLD.Send(&s_buffer[ibuf],pnbuf,mpi_point_type,0,Mpi_Tag::partial_buffer);      
	MPI::COMM_WORLD.Send(&s_counter_buffer[icbuf],pncbuf,MPI::INT,0,Mpi_Tag::partial_buffer);  
	ibuf+=pnbuf;
	icbuf+=pncbuf;
	if (ibuf+pnbuf>=nbuf) break;
      }
      pnbuf = nbuf  - ibuf;
      pncbuf= ncbuf - icbuf;
      msg.Tracking()<<"CPU"<<rank<<" sending final "<<pnbuf<<" points now"<<endl;      
      MPI::COMM_WORLD.Send(&s_buffer[ibuf],pnbuf,mpi_point_type,0,Mpi_Tag::complete_buffer);
      MPI::COMM_WORLD.Send(&s_counter_buffer[icbuf],pncbuf,MPI::INT,0,Mpi_Tag::complete_buffer);  
    }
  }
  else {
    // if master
    msg.Tracking()<<" Master creates MPI_Points "<<"("<<prea_table.size()<<")"<<endl;
    // receive all buffers!
    for (int i=1;i<size;++i) {
      r_buffers[i].tag=Mpi_Tag::partial_buffer;
    }
 
    for (;;) {
      bool rerun=0;
      for (int i=1;i<size;++i) {
	if (r_buffers[i].tag==Mpi_Tag::partial_buffer) {
	  MPI::Status    r_status;
	  MPI::COMM_WORLD.Recv(r_buffers[i].buff,NMAX,mpi_point_type, i,MPI::ANY_TAG , r_status);
	  MPI::COMM_WORLD.Recv(r_buffers[i].counters,2*NMAX,MPI::INT, i,MPI::ANY_TAG);
	  r_buffers[i].nbuff=r_status.Get_count(mpi_point_type);
	  r_buffers[i].tag  =(Mpi_Tag::code)r_status.Get_tag();
	  msg.Tracking()<<"recved buffer from CPU"<<i
	      <<" ("<<r_buffers[i].nbuff<<")"<<endl;
	  if (r_buffers[i].tag==Mpi_Tag::partial_buffer) rerun=1;
	}
      }

      // translate in own list:
      for (int i=1;i<size;++i) {
	if (r_buffers[i].tag==Mpi_Tag::complete_buffer||r_buffers[i].tag==Mpi_Tag::partial_buffer) {
	  int jc=0;
	  for (int j=0;j<r_buffers[i].nbuff;j+=dep, jc+=2) {
	    Point* ph = new Point[dep];
	    MPI2Point(&(r_buffers[i].buff[j]),ph);
	    ReplaceVertex(ph);
	    prea_table.push_back(Pre_Amplitude(ph,
		 r_buffers[i].counters[jc],r_buffers[i].counters[jc+1]));
	  }
	}
	if (r_buffers[i].tag==Mpi_Tag::complete_buffer) r_buffers[i].tag=Mpi_Tag::empty_buffer;
      }
      if (!rerun) break;
    }
    // distribute 
  }
#endif
}


void  Amplitude_Generator::DistributePreAmplitudes() {
#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();

  int dep  = single_top->depth;

  if (size==1) return;

  if (rank==0) {
    // master sends combined result to all cpus
    int nbuf  = prea_table.size()*dep;
    int ncbuf = 2*prea_table.size();

    msg.tracking()<<"Master creates MPI_Points for final transmission "<<nbuf<<"("<<prea_table.size()<<")"<<endl;
    s_buffer=new MPI_Point[nbuf];
    s_counter_buffer = new int[ncbuf];

    // - fill buffer
    int ibuf  =0;
    int icbuf =0;
    for (int i=0; i<prea_table.size();++i) {      
      Point2MPI(prea_table[i].p,&s_buffer[ibuf]);
      ibuf+=dep;
      s_counter_buffer[2*i]   = prea_table[i].top;
      s_counter_buffer[2*i+1] = prea_table[i].perm;
    }

    // check if receive-buffer is big enough !!!
    if (nbuf<=NMAX) {
    
      // send local list
      for (int cpu=1;cpu<size;++cpu) {
	msg.Tracking()<<"Master sends to"<<cpu<<"  "<<nbuf<<" points now"<<endl;
	MPI::COMM_WORLD.Send(s_buffer,nbuf,mpi_point_type,cpu,Mpi_Tag::complete_buffer);
	MPI::COMM_WORLD.Send(s_counter_buffer,ncbuf,MPI::INT,cpu,Mpi_Tag::complete_buffer);  
      }
    }
    else {
      ibuf  =0;
      icbuf =0;
      int pnbuf=NMAX/dep;
      int pncbuf=pnbuf*2;
      pnbuf*=dep;
      for (;;) {
	for (int cpu=1;cpu<size;++cpu) {
	  msg.Tracking()<<"Master sends to"<<cpu<<"  "<<pnbuf<<" points now"<<endl;
	  MPI::COMM_WORLD.Send(&s_buffer[ibuf],pnbuf,mpi_point_type,cpu,Mpi_Tag::partial_buffer);      
	  MPI::COMM_WORLD.Send(&s_counter_buffer[icbuf],pncbuf,MPI::INT,cpu,Mpi_Tag::partial_buffer);  
	}
	ibuf+=pnbuf;
	icbuf+=pncbuf;
	if (ibuf+pnbuf>=nbuf) break;
      }
      pnbuf = nbuf  - ibuf;
      pncbuf= ncbuf - icbuf;
      for (int cpu=1;cpu<size;++cpu) {
	msg.Tracking()<<"Master sends to"<<cpu<<"  "<<pnbuf<<" (final) points now"<<endl;
	MPI::COMM_WORLD.Send(&s_buffer[ibuf],pnbuf,mpi_point_type,cpu,Mpi_Tag::complete_buffer);
	MPI::COMM_WORLD.Send(&s_counter_buffer[icbuf],pncbuf,MPI::INT,cpu,Mpi_Tag::complete_buffer);  
      }
    }
  }
  else {
    // clear old Preas:
    for (Pre_Ampl_List::iterator pit=prea_table.begin();pit!=prea_table.end();++pit)
      delete [] (*pit).p;
    prea_table.clear();


    msg.Tracking()<<" Slave "<<rank<<" receives MPI_Points "<<endl;
    r_buffers[0].tag=Mpi_Tag::partial_buffer;
 
    for (;;) {
      if (r_buffers[0].tag==Mpi_Tag::partial_buffer) {
	MPI::Status    r_status;
	MPI::COMM_WORLD.Recv(r_buffers[0].buff,NMAX,mpi_point_type, 0,MPI::ANY_TAG , r_status);
	MPI::COMM_WORLD.Recv(r_buffers[0].counters,2*NMAX,MPI::INT, 0,MPI::ANY_TAG);
	r_buffers[0].nbuff=r_status.Get_count(mpi_point_type);
	r_buffers[0].tag  =(Mpi_Tag::code)r_status.Get_tag();
	msg.Tracking()<<" Slave "<<rank<<" received buffer from Master"
	    <<" ("<<r_buffers[0].nbuff<<")"<<endl;
      }

      // translate in own list:
      if (r_buffers[0].tag==Mpi_Tag::complete_buffer||r_buffers[0].tag==Mpi_Tag::partial_buffer) {
	int jc=0;
	for (int j=0;j<r_buffers[0].nbuff;j+=dep, jc+=2) {
	  Point* ph = new Point[dep];
	  MPI2Point(&(r_buffers[0].buff[j]),ph);
	  ReplaceVertex(ph);
	  prea_table.push_back(Pre_Amplitude(ph,
		       r_buffers[0].counters[jc],r_buffers[0].counters[jc+1]));
	}
      }
      if (r_buffers[0].tag==Mpi_Tag::complete_buffer) break;
    }
  }
#endif
}


void Amplitude_Generator::ReplaceVertex(Point * p) {
  if (!p) return;
  Single_Vertex * sv = p->v;
  if (!sv) return;

  int hit=0;
  Vertex_List & vl= v_table[sv->in[0]];
  for (int i=0;i<vl.size();++i) {
    if (vl[i]->in[1] == sv->in[1]  && vl[i]->in[2] == sv->in[2] && vl[i]->in[3] == sv->in[3]) {
      delete sv;
      p->v=vl[i];
      hit=1;
      break;
    }
  }
  if (!hit) {
    msg.Error()<<" vertex not found , something wrong with MPI "<<endl;
    abort();
  }
}

void Amplitude_Generator::Unite(Point* p,Point* pdel)
{
  int depth = single_top->depth;
  Point psave;
  for (short int i=0;i<depth;i++) {
    if (p[i].number==pdel[i].number) {
      if (p[i].fl!=pdel[i].fl) {
	// Double Counting !!!!!!!!!!!!!!
	int sw1 = 1;
	for (short int j=0;j<p[i].nextra;j++) {
	  if (p[i].extrafl[j]==pdel[i].fl) {
	    sw1 = 0;
	    break;
	  }
	}
	if (sw1==1) {
	  psave = p[i];
	
	  if (p[i].nextra>0) delete[] p[i].extrafl;
	  delete[] p[i].cpl;
	  
	  int nfl  = 1+pdel[i].nextra+p[i].nextra;
	  int ncpl = pdel[i].ncpl+p[i].ncpl;
	  
	  p[i].extrafl = new Flavour[nfl];
	  p[i].cpl     = new Complex[ncpl];
	  
	  //Flavour
	  int count = 0;
	  for (short int j=0;j<psave.nextra;j++)
	    p[i].extrafl[j] = psave.extrafl[j];   
	  count += psave.nextra;
	  
	  p[i].extrafl[count] = pdel[i].fl;
	  count++;
	  for (short int j=0;j<pdel[i].nextra;j++)
	    p[i].extrafl[count+j] = pdel[i].extrafl[j];			
	  p[i].nextra = nfl;
	  
	  //Couplings
	  count = 0;
	  for (short int j=0;j<psave.ncpl;j++)
	    p[i].cpl[j] = psave.cpl[j];   
	  count += psave.ncpl;
	  for (short int j=0;j<pdel[i].ncpl;j++)
	    p[i].cpl[count+j] = pdel[i].cpl[j];			
	  p[i].ncpl = ncpl;
	  
	  //previous couplings too
	  int hit = -1;
	  for (short int j=0;j<depth;j++) {
	    if (p[j].left==&p[i] || p[j].right==&p[i]) {
	      hit = j;
	      break;
	    }
	  }
	  if (hit!=-1) {
	    psave = p[hit];
	    int ncpl = pdel[hit].ncpl+p[hit].ncpl;	  
	    delete[] p[hit].cpl;	  
	    p[hit].cpl = new Complex[ncpl];
	    
	    //Couplings
	    count = 0;
	    for (short int j=0;j<psave.ncpl;j++)
	      p[hit].cpl[j] = psave.cpl[j];   
	    count += psave.ncpl;
	    for (short int j=0;j<pdel[hit].ncpl;j++)
	      p[hit].cpl[count+j] = pdel[hit].cpl[j];			
	    p[hit].ncpl = ncpl;	  
	  }
	  else AORGTOOLS::msg.Error()<<"Error in Amplitude_Generator"<<endl;
	}
      }
    }
  } 
}

int Amplitude_Generator::Single_Compare(Point* p1, Point* p2)
{
  //zero check
  if (p1==0) {
    if (p2==0) return 1;
    else return 0;
  }
  else {
    if (p2==0) return 0;
  }
  
  //Flavour equal....
  if (p1->fl!=p2->fl) return 0;

  //outgoing number equal
  if ((p1->left==0) && (p2->left==0)) {
    if (p1->number!=p2->number) return 0;
                           else return 1;
  }

  //Check extended Color_Functions
  if (p1->Color->String()==p2->Color->String()) {
    if (p1->Color->Next && p2->Color->Next) {
      Color_Function* ctmp1 = p1->Color;
      Color_Function* ctmp2 = p2->Color;
      while(ctmp1->Next) {
	if((ctmp1->Next)->String()!=(ctmp2->Next)->String()) return 0;
	ctmp1 = ctmp1->Next;
	ctmp2 = ctmp2->Next;
      }
    }
    else {
      if (p1->Color->Next!=0 || p2->Color->Next!=0) return 0;
    }
  }
  else return 0;

  // return 1 if equal and 0 if different
  
  if (Single_Compare(p1->middle,p2->middle)) {
    int sw1 = Single_Compare(p1->left,p2->left);
    if (sw1) sw1 = Single_Compare(p1->right,p2->right);
    if (sw1==0) {
      sw1 = Single_Compare(p1->left,p2->right);
      if (sw1) sw1 = Single_Compare(p1->right,p2->left);
    }
    return sw1;
  }

  if (Single_Compare(p1->middle,p2->left)) {
    int sw1 = Single_Compare(p1->left,p2->middle);
    if (sw1) sw1 = Single_Compare(p1->right,p2->right);
    if (sw1==0) {
      sw1 = Single_Compare(p1->left,p2->right);
      if (sw1) sw1 = Single_Compare(p1->right,p2->middle);
    }
    return sw1;
  }

  if (Single_Compare(p1->middle,p2->right)) {
    int sw1 = Single_Compare(p1->right,p2->middle);
    if (sw1) sw1 = Single_Compare(p1->left,p2->left);
    if (sw1==0) {
      sw1 = Single_Compare(p1->right,p2->left);
      if (sw1) sw1 = Single_Compare(p1->left,p2->middle);
    }
    return sw1;
  }

  return 0;
}

int Amplitude_Generator::Kill_Off(Single_Amplitude* &first)
{
  Single_Amplitude* last;
  last = first;
  Single_Amplitude* f1 = first;
  Single_Amplitude* f2;
  int count=0;
  while (f1) {
    if (f1->on==0) {
      ++count;
      if (f1==first) {
	first = f1->Next;
	f1 = f1->Next;
	delete last;
	last = first;
      }
      else {
	last->Next = f1->Next;
	f2 = f1;
	f1 = f1->Next;
	delete f2;
      }
    }
    else {
      last = f1;
      f1 = f1->Next;
    }
  }
  return count;
}

int Amplitude_Generator::FindQEDOrder(Point * p,int & countQED)
{
  if (!p) return countQED;
  
  int hit = 0;
  
  //Vector-Boson propagators
  if (p->number>99 && p->fl.IsVector() && !(p->fl.IsGluon())) {
    countQED += 2;
    hit       = 1; 
  }    
  //External Vector-Boson 
  if (p->number<99 && p->fl.IsVector() && !(p->fl.IsGluon())) {
    countQED += 1;
    hit       = 1;
  }
  
  //triple and quartic Vector-Boson interactions
  if (hit) {
    if (p->left   && p->left->fl.IsVector() && !(p->left->fl.IsGluon()))     countQED -= 1;
    if (p->right  && p->right->fl.IsVector() && !(p->right->fl.IsGluon()))   countQED -= 1;
    if (p->middle && p->middle->fl.IsVector() && !(p->middle->fl.IsGluon())) countQED -= 1;
  }
  
  FindQEDOrder(p->left,countQED);
  FindQEDOrder(p->right,countQED);
  if (p->middle) FindQEDOrder(p->middle,countQED);
  return countQED;
}

int Amplitude_Generator::FindQCDOrder(Point * p,int & countQCD)
{
  if (!p) return countQCD;
  
  int hit = 0;
  
  //Gluon propagators
  if (p->number>99 && p->fl.IsGluon()) {
    countQCD += 2;
    hit       = 1;
  }
  //External gluon 
  if (p->number<99 && p->fl.IsGluon()) {
    countQCD += 1;
    hit       = 1;
  }

   //triple and quartic Gluon vertices 
  if (hit) {
    if (p->left   && p->left->fl.IsGluon())   countQCD -= 1;
    if (p->right  && p->right->fl.IsGluon())  countQCD -= 1;
    if (p->middle && p->middle->fl.IsGluon()) countQCD -= 1;
  }
  
  FindQCDOrder(p->left,countQCD);
  FindQCDOrder(p->right,countQCD);
  if (p->middle) FindQCDOrder(p->middle,countQCD);
  return countQCD;
}

void Amplitude_Generator::KillHigherOrders(Single_Amplitude * & first)
{
  msg.Tracking()<<"In KillHigherOrders() nEW = "<<nEW<<", nQCD = "<<nQCD<<endl; 
    
  Single_Amplitude* last;
  last = first;
  Single_Amplitude* f1 = first;
  Single_Amplitude* f2;
  int count=0;
  while (f1) {
    int hitQED = 0;
    int hitQCD = 0;
    hitQED = FindQEDOrder(f1->GetPointlist(),hitQED);
    hitQCD = FindQCDOrder(f1->GetPointlist(),hitQCD);
    msg.Tracking()<<"hitQED / hitQCD "<<hitQED<<" "<<hitQCD<<endl;
    if (hitQED > nEW || hitQCD > nQCD) {
      msg.Tracking()<<"Diagram has to be kicked "<<endl;
      ++count;
      if (f1==first) {
	first = f1->Next;
	f1 = f1->Next;
	delete last;
	last = first;
      }
      else {
	last->Next = f1->Next;
	f2 = f1;
	f1 = f1->Next;
	delete f2;
      }
    }
    else {
      last = f1;
      f1 = f1->Next;
    }
  }
  msg.Tracking()<<"Kicked number of diagrams (Amplitude_Generator::KillHigherOrders()) "<<count<<endl;
}


int *  Amplitude_Generator::DivideComparisons(int cpu_size, int nampl, int rank) {
  // determine share (multi cpu only)
  double total_compare = double(nampl)*(double(nampl)-1.)/2.;
  double  rest_share    = total_compare;
  double * compare_share = new double[cpu_size];
  for (int cpu=0;cpu<cpu_size;++cpu) {
    // ideal shares
    compare_share[cpu]=floor(rest_share/double(cpu_size-cpu));
    rest_share-=compare_share[cpu];
  }

  rest_share    = total_compare;
  int * start_no= new int[cpu_size];
  for (int cpu=0;cpu<cpu_size;++cpu) {
    // ideal shares
    rest_share-=compare_share[cpu];
    double est_no = 0.5 + sqrt(0.25+2.*rest_share);
    int no = int(est_no+0.5);      // round
    start_no[cpu]=no;
    // real shares
    if (rank==0) msg.Tracking()<<" cpu"<<cpu<<" "<<start_no[cpu];
    no = 0;
    if (cpu>0) no=nampl-start_no[cpu-1];
    if (rank==0) msg.Tracking()<<" -> "<<no<<"  ";
    double cmp=total_compare;
    if (cpu>0) cmp=double(start_no[cpu-1])*double((start_no[cpu-1]-1.))*0.5;
    cmp-=double(start_no[cpu])*double(start_no[cpu]-1.)*0.5;
    if (rank==0) msg.Tracking()<<" ("<<compare_share[cpu]<<" -> "<<cmp<<") "<<endl;
    compare_share[cpu]=cmp;
  }
  // translate to reverse shares
  for (int cpu=cpu_size-1;cpu>=0;--cpu) {
    int no = 0;
    if (cpu>0) no=nampl-start_no[cpu-1];
    start_no[cpu]=no;
  }
  return start_no;
}

void Amplitude_Generator::Compare(Single_Amplitude* &first)
{
  int rank=0;
  int size =1;
#ifdef _USE_MPI_
  rank = MPI::COMM_WORLD.Get_Rank();
  size = MPI::COMM_WORLD.Get_size();
#endif

  Single_Amplitude* f1;
  Single_Amplitude* f2;
  msg.Debugging()<<" in Amplitude_Generator::Compare "<<endl;

  Single_Amplitude*   start_ampl=first;
  Single_Amplitude*   stop_ampl =0;
  
  // count amplitudes
  int nampl=0;
  f1 = first;
  while (f1) {
    ++nampl;
    f1=f1->Next;
  } 

  msg.Debugging()<<"CPU"<<rank<<" got "<<nampl<<endl;

  int * start_no=DivideComparisons(size,nampl,rank);

  // determine new start_amplitude and 
  f1 = first;
  int iampl=0;
  while (f1) {
    if (iampl==start_no[rank]) {
      start_ampl=f1;
      if (rank==size-1) break;
    }
    if (rank<size-1) {
      if (iampl==start_no[rank+1]) {
	stop_ampl=f1;
	break;
      }
    }
    ++iampl;
    f1=f1->Next;
  }

  delete [] start_no;

  // start comparison
  Point* p1;
  Point* p2;  

  int ncomps=0;
  int noffs=0;
  f1 = start_ampl;
  while (f1) { 
    p1 = f1->GetPointlist();
    f2 = f1->Next;
    while (f2) {
      p2 = f2->GetPointlist();
      ++ncomps;
      int sw1 = Single_Compare(p1,p2);
      if (sw1==1) {
	if (f2->on) ++noffs;
	f2->on = 0;
      }
      f2 = f2->Next;
    }
    f1 = f1->Next;
    if (f1==stop_ampl) break;
  }
  msg.Debugging()<<"CPU"<<rank<<" did "<<ncomps<<" compares "<<endl;
  msg.Debugging()<<"CPU"<<rank<<" switch off "<<noffs<<" amplitudes"<<endl;


#ifdef _USE_MPI_
  // collect off-switches
  if (rank>0) {
    int * sw_buffer = new int[nampl];
    f1 = first;
    for (int iamp=0;iamp<nampl;++iamp) {
      sw_buffer[iamp]=f1->on;
      f1 = f1->Next;
    }
    // send local off-switches
    MPI::COMM_WORLD.Send(sw_buffer,nampl,MPI::INT,0,Mpi_Tag::on_switches);
    // receive global off-switches
    MPI::COMM_WORLD.Recv(sw_buffer,nampl,MPI::INT,0,Mpi_Tag::on_switches);
    // set switches in amplitudes
    f1=first;
    for (int iamp=0;iamp<nampl;++iamp) {
      f1->on=sw_buffer[iamp];
      f1 = f1->Next;
    }
    delete [] sw_buffer;
  }
  else {
    int * sw_buffer = new int[nampl];
    int * all_sw_buffer = new int[nampl];
    f1=first;
    for (int iamp=0;iamp<nampl;++iamp) {
      all_sw_buffer[iamp]=f1->on;
      f1 = f1->Next;
    }
    for (int cpu=1;cpu<size;++cpu) {
      // receive slave off-switches
      MPI::COMM_WORLD.Recv(sw_buffer,nampl,MPI::INT,cpu,Mpi_Tag::on_switches);
      // combine slave off-switches with master 
      for (int iamp=0;iamp<nampl;++iamp) {
	all_sw_buffer[iamp]=all_sw_buffer[iamp] & sw_buffer[iamp];
      }
    }
    for (int cpu=1;cpu<size;++cpu) {
      // send global off-switches to all slaves
      MPI::COMM_WORLD.Send(all_sw_buffer,nampl,MPI::INT,cpu,Mpi_Tag::on_switches);
    }

    // set switches in amplitudes
    f1=first;
    for (int iamp=0;iamp<nampl;++iamp) {
      f1->on=all_sw_buffer[iamp];
      f1 = f1->Next;
    }

    delete [] sw_buffer;
    delete [] all_sw_buffer;
  }
#endif
  
  int killed=Kill_Off(first);
  msg.Debugging()<<"CPU"<<rank<<" killed "<<killed<<endl;

  f1 = first;
  nampl=0;
  while (f1) {
    ++nampl;
    f1=f1->Next;
  } 
  msg.Debugging()<<"CPU"<<rank<<" left with "<<nampl<<endl;
}

//==================================================================

Point* Amplitude_Generator::FindNext(Point* p)
{
  if (p==0) return 0;
  if (p->left->m==0) return p;
  if (p->right->m==0) return p;
  if ((p->middle!=0) && (p->middle->m==0)) return p;
  
  FindNext(p->left);
  FindNext(p->right);
  FindNext(p->middle);

  return 0;

}

int Amplitude_Generator::ShrinkProps(Point*& p,Point*& pnext, Point*& pcopy, Point*& beg_pcopy,
				     vector<Point*>& pcollist)
{
  if (p->left==0 || pnext->left==0) return 0;
  
  if (p->v->nleg==4 || pnext->v->nleg==4) return 0;
  
  if (pnext->m==1) return 0;

  if (pnext->fl.IsFermion()) return 0;
    
  int hit = 0;
  int in[4] = {0};
    
  Flavour flav[4];
  
  flav[0] = p->fl;
  if (p->number<99 && p->b==1) in[0] = 1;

  if (p->left->number==pnext->number) {
    flav[1] = pnext->left->fl;
    if (pnext->left->number<99 && pnext->left->b==-1)   in[1] = 1;
    flav[2] = pnext->right->fl;
    if (pnext->right->number<99 && pnext->right->b==-1) in[2] = 1;
    flav[3] = p->right->fl;
    if (p->right->number<99 && p->right->b==-1)         in[3] = 1;
  }
  if (p->right->number==pnext->number) {
    flav[1] = p->left->fl;
    if (p->left->number<99 && p->left->b==-1)           in[1] = 1;
    flav[2] = pnext->left->fl;
    if (pnext->left->number<99 && pnext->left->b==-1)   in[2] = 1;
    flav[3] = pnext->right->fl;
    if (pnext->right->number<99 && pnext->right->b==-1) in[3] = 1;
  }  

  //barflags 
  
  Vertex* v = p_model->GetVertex();
  
  for (short int i=0;i<v->MaxNumber4();i++) {
    if ((*v)(i)->on) {
      
      Single_Vertex test;
      
      test = *((*v))(i);
      
      if (in[0]) test.in[0] = test.in[0].Bar();
      if (in[1]) test.in[1] = test.in[1].Bar();
      if (in[2]) test.in[3] = test.in[3].Bar();
      if (in[3]) test.in[2] = test.in[2].Bar();

      if (flav[0]==test.in[0] && 
	  flav[1]==test.in[1] &&
	  flav[2]==test.in[3] &&
	  flav[3]==test.in[2]) 
	{
	  //AORGTOOLS::msg.Out()<<"4 leg Vertex found !!!"<<endl;
	  //AORGTOOLS::msg.Out()<<"x";
	  
	  hit = 1;
	  
	  pcopy->v      = (*v)(i);
	  for (short int k=0;k<4;k++) pcopy->cpl[k] = (*v)(i)->cpl[k];
	  //set pcopy legs
	  
	  if (p->left->number==pnext->number) {	    
	    pcopy->middle = pcopy->left->right;
	    pcopy->left   = pcopy->left->left;
	    pcopy->right  = pcopy->right;
	  }
	  if (p->right->number==pnext->number) {	    
	    pcopy->left   = pcopy->left;
	    pcopy->middle = pcopy->right->left;
	    pcopy->right  = pcopy->right->right;
	  }
	  
	  //setting the contraction flags
	  pnext->m  = 1;
	  
	  if ((*v)(i)->ncf==1) {
	    
	    *(pcopy->Color)   = *((*v)(i)->Color);
	    *(pcopy->Lorentz) = *((*v)(i)->Lorentz);
	    
	    break;
	  }
	  else {
	    
	    for (short int k=0;k<(*v)(i)->ncf;k++) {

	      *(pcopy->Color)   = ((*v)(i)->Color)[k];

	      if (((*v)(i)->Color[k]).Next!=0) {
		Color_Function* cforig = ((*v)(i)->Color[k]).Next;
		Color_Function* cfcopy = pcopy->Color;
		while (cforig) {
		  cfcopy->Next = new Color_Function(*cforig);
		  cfcopy = cfcopy->Next;
		  cforig = cforig->Next;
		}
	      }	      
	      
	      *(pcopy->Lorentz) = ((*v)(i)->Lorentz)[k];
	      
	      Color_Function* cfmemo = pcopy->Color;
	      
	      //set the contraction indices: 4 -> (prop->number) 
	      while (cfmemo) {
		for (short int i=0;i<3;i++) {
		  if ((cfmemo->type==cf::D || cfmemo->type==cf::G) && i==2) break;
		  if (cfmemo->partarg[i]==4) {
		    cfmemo->partarg[i] = pnext->number;
		  }
		  }
		cfmemo = cfmemo->Next;
	      }

	      //fill the vector pcollist
	      int ll = 0;
	      Point* ptmp = new Point[single_top->depth];
	      top->Copy(beg_pcopy,ptmp,ll);
	      pcollist.push_back(ptmp);
	    }
	    break;
	  }
	}
      else hit = 0;
    }
  }
  return hit;
}

int Amplitude_Generator::EvalPointlist(Point*& porig, Point*& pcopy,Point*& beg_pcopy,
				       vector<Point*>& pcollist)
{
  if(porig==0) return 0;
  
  if (ShrinkProps(porig,porig->right,pcopy,beg_pcopy,pcollist))    return 1;
  
  if (ShrinkProps(porig,porig->left,pcopy,beg_pcopy,pcollist))     return 1;
  
  if (EvalPointlist(porig->left,pcopy->left,beg_pcopy,pcollist))   return 1;
  if (EvalPointlist(porig->right,pcopy->right,beg_pcopy,pcollist)) return 1;
  if (EvalPointlist(porig->middle,pcopy->middle,beg_pcopy,pcollist)) return 1;

  return 0;
}

void Amplitude_Generator::CheckFor4Vertices(Single_Amplitude* &first)
{
  Single_Amplitude* f1 = first;
  Single_Amplitude* extraAmpl;
  Single_Amplitude* f2;
  Point* p;
  int dep       = single_top->depth;
  Point* pcopy  = new Point[dep];
  vector<Point*> pcollist;//vector with different color structures
  int counter   = 0;
  int amplcount = 0;
 
  AORGTOOLS::msg.Debugging()<<"=============== CheckFor4Vertices ======================="<<endl;  
  
  while (f1) { 
    amplcount++;
    
    p = f1->GetPointlist();
    
    for (int i=0; i<dep; i++) p[i].m = 0; 
    
    while (p) {
      //p initial Pointlist with m flags set, pcopy is the contracted Pointlist !!!
      int ll = 0;
      top->Copy(p,pcopy,ll);
      if (EvalPointlist(p,pcopy,pcopy,pcollist)) {
	if (pcollist.size()==0) {
	  counter++;
	  //AORGTOOLS::msg.Out()<<"Counter: "<<counter<<endl;
	  
	  extraAmpl = new Single_Amplitude(pcopy,f1->topnum,f1->permnum,b,dep,N,top,BS,fl,shand);
	  extraAmpl->Next = 0;
	  
	  f2 = first;
	  
	  while (f2) {
	    if (f2->Next==0) {
	      f2->Next = extraAmpl;
	      break;
	    }
	    f2 = f2->Next;
	  }
	}
	else { 
	  for (short int j=0;j<pcollist.size();j++) {
	    counter++;
	    //AORGTOOLS::msg.Out()<<"Counter: "<<counter<<endl;
	    
	    extraAmpl = new Single_Amplitude(pcollist[j],f1->topnum,f1->permnum,b,dep,N,top,BS,fl,shand);
	    extraAmpl->Next = 0;
	    
	    f2 = first;
	    while (f2) {
	      if (f2->Next==0) {
		f2->Next = extraAmpl;
		break;
	      }
	      f2 = f2->Next;
	    }
	  }
	  //erase all elements of pcollist
	  pcollist.clear();
	}
	p = FindNext(p);
      }
      else break;
    }
    f1 = f1->Next;
  }
  AORGTOOLS::msg.Debugging()<<"=========================================================="<<endl;
}

//=====================================================================

int Amplitude_Generator::Count4G(Point * p) {
  if (!p) return 0;
  int v4=0;
  v4+=Count4G(p->left);
  v4+=Count4G(p->right);
  if (p->middle) {
    v4+=Count4G(p->middle);
    Flavour gluon=Flavour(kf::gluon);
    if ((p->fl==gluon)&&(p->left->fl==gluon)&&(p->middle->fl==gluon)&&(p->right->fl==gluon))
      v4+=1;
  }
  return v4;
}

int Amplitude_Generator::CountRealAmplitudes(Single_Amplitude* first) 
{
  Single_Amplitude * f1=first;
  int counts[4]={0,0,0,0};
  while (f1) {
    Point * p =f1->GetPointlist();
    int v4 = Count4G(p);
    if (v4<4) ++(counts[v4]);
    else {
      cerr<<" to many four verticies in one amplitude "<<endl;
    }
    f1=f1->Next; 
  }
  msg.Out()<<endl;
  msg.Out()<<" in CountRealAmplitudes: "
      <<counts[0]<<" "<<counts[1]<<" "<<counts[2]<<" "<<counts[3]<<endl;
  int ra=counts[0];
  int three=3;
  for (int i=1;i<4;++i) {
    ra+=counts[i]/three;
    three*=3;
  } 
  msg.Out()<<" real ampls = "<<ra<<endl;
  return ra;
}


Single_Amplitude* Amplitude_Generator::Matching()
{
  AORGTOOLS::msg.Out()<<"Matching of topologies..."<<endl;
  int nloop = N-1;
  short int i,j;
  int* ii = new int[nloop];
  for (i=0;i<nloop;i++) ii[i] = 0;
  int* perm = new int[N];
  int sw1;
  int qsum,l1sum,l2sum,l3sum; 
  //int chsum,neusum;           
  Single_Amplitude* first_amp;
  int over = 0;
  perm[0]  = 0;
  first_amp = 0;

  int count = 0;

#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();

  int top_range =  single_top->number/size;

  int start_top = 0;
  int end_top = single_top->number;

  start_top=rank*top_range;
  if (rank<size-1)
    end_top=(rank+1)*top_range;


#else
  int rank = 0;
  int size = 1;
  int start_top = 0;
  int end_top = single_top->number;
#endif

  msg.Tracking()<<" "<<rank<<"/"<<size<<" CPUS  ("<<start_top<<"-"<<end_top<<")"<<endl;

  for(;;) {
    sw1 = 1;
    for (j=0;j<nloop;j++) perm[j+1] = ii[j];
    
    for(i=0;i<N;i++) {
      for (j=i+1;j<N;j++) 
	if (perm[i]==perm[j]) {sw1 = 0;break;}
    }

    if (sw1) count++;

    if (sw1) { 
      if (fl[0].IsQuark() && !fl[0].IsAnti()) {
	if (!(fl[perm[N-1]].IsQuark() && 
	      (b[perm[N-1]]==-1 && (fl[perm[N-1]].IsAnti() || fl[perm[N-1]].Majorana())) ||
	      (b[perm[N-1]]==1  && (!(fl[perm[N-1]].IsAnti()) || fl[perm[N-1]].Majorana())) ))  sw1 = 0;
      }
      if (fl[0].IsQuark() && fl[0].IsAnti()) {
	if (!(fl[perm[1]].IsQuark() && 
	      (b[perm[1]]==-1 && 
	       (!fl[perm[1]].IsAnti() || fl[perm[1]].Majorana())) ||
	      (b[perm[1]]==1 &&  
	       ((fl[perm[1]].IsAnti()) || fl[perm[1]].Majorana())) ))  sw1 = 0;
      }
      //Anti-Leptons too!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (fl[0].IsLepton()) {
	if (!(fl[perm[N-1]].IsLepton() //|| fl[perm[N-1]].IsSlepton() 
	      && 
	      (b[perm[N-1]]==-1 && (fl[perm[N-1]].IsAnti() || fl[perm[N-1]].Majorana())) ||
	      (b[perm[N-1]]==1  && (!(fl[perm[N-1]].IsAnti()) || fl[perm[N-1]].Majorana())) )) 
	  sw1 = 0;
      }
    }
    if (sw1) {
      qsum=l1sum=l2sum=l3sum=0;
      for (j=0;j<N-1;j++) {
	if (fl[perm[j]].IsQuark()) {
	  if (fl[perm[j]].IsAnti()) {
	    if (b[perm[j]]==1) qsum--;
	    else qsum++;
	  }
	  else {
	    if (b[perm[j]]==1) qsum++;
	    else qsum--;
	  }
	  if (!fl[0].IsAnti()) {
	    if ( b[perm[N-1]]==-1 && fl[perm[N-1]].IsAnti() && qsum>0 ) 
	      {sw1=0;break;} // s-channel
	    if ( b[perm[N-1]]==1 && !(fl[perm[N-1]].IsAnti()) && qsum>0 )
	      {sw1=0;break;} // t-channel,
	  }
	  else {
	    if ( b[perm[1]]==-1 && fl[perm[1]].IsAnti() && qsum<0 ) 
	      {sw1=0;break;} // s-channel
	    if ( b[perm[1]]==1 && !(fl[perm[1]].IsAnti()) && qsum<0 )
	      {sw1=0;break;} // t-channel,	    
	  }
	}
	if (fl[perm[j]].IsLepton()) {
	  if ( (fl[perm[j]]==Flavour(kf::e)) || 
	       (fl[perm[j]]==Flavour(kf::nue)) ) {
	    if (b[perm[j]]==1) l1sum++;
	    else l1sum--;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::e).Bar()) || 
	       (fl[perm[j]]==Flavour(kf::nue).Bar()) ) {
	    if (b[perm[j]]==1) l1sum--;
	    else l1sum++;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::mu)) || 
	       (fl[perm[j]]==Flavour(kf::numu)) ) {
	    if (b[perm[j]]==1) l2sum++;
	    else l2sum--;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::mu).Bar()) || 
	       (fl[perm[j]]==Flavour(kf::numu).Bar()) )  {
	    if (b[perm[j]]==1) l2sum--;
	    else l2sum++;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::tau)) || 
	       (fl[perm[j]]==Flavour(kf::nutau)) ) {
	    if (b[perm[j]]==1) l3sum++;
	    else l3sum--;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::tau).Bar()) || 
	       (fl[perm[j]]==Flavour(kf::nutau).Bar()) )  {
	    if (b[perm[j]]==1) l3sum--;
	    else l3sum++;
	  }
	}
	if ( (b[perm[N-1]]==-1) && (fl[perm[N-1]].IsAnti()) && 
	     ( (l1sum>0) || (l2sum>0) || (l3sum>0))  ) {sw1=0;break;} // s-channel, e+e-.
	if ( (b[perm[N-1]]==1) && !(fl[perm[N-1]].IsAnti()) &&
	     ( (l1sum>0) || (l2sum>0) || (l3sum>0)) )  {sw1=0;break;} // t-channel, e+e-.
      }
    }
    if (sw1) {
      for (j=start_top;j<end_top;j++) {
	perm++;
	int pnum = 100;
	Set_End(&single_top->p[j][0],perm,pnum);
	perm -= N;
	single_top->p[j][0].number = *perm;
	single_top->p[j][0].fl     = fl[*perm];
	single_top->p[j][0].b      = b[*perm];
	
	Set_Props(single_top->p[j],2*N-3,first_amp,perm,j,count);
	
	//Print_P(&single_top->p[j][0]);
      }
    }     
    for (j=nloop-1;j>=0;j--) {
      if ((ii[j]+1)<N) {
	ii[j]++;            
	break;
      }
      else {
	ii[j] = 0;
	if (j==0) over = 1;
      }
    }
    if (over) break;
  }
  delete[] ii;
  delete[] perm;
  
#ifdef _USE_MPI_
  // now collect all Pre_Amplitudes and make Single_Amplitudes
  CollectPreAmplitudes();

  if (rank==0) {
    msg.Debugging()<<" sorting "<<prea_table.size()<<" PreAmpl "<<endl;
    std::stable_sort(prea_table.begin(),prea_table.end(),Compare_Pre_Amplitudes());
    msg.Debugging()<<" sorting done"<<endl;
  }

  DistributePreAmplitudes();
#endif

  //  if (rank==0)
  CreateSingleAmplitudes(first_amp);
//   else {
//     int dummy;
//     MPI::COMM_WORLD.Recv(&dummy,1,MPI::INT, 0,MPI::ANY_TAG);
// }

  msg.Debugging()<<"Now CheckFor4Vertices()"<<endl;

  CheckFor4Vertices(first_amp);
  
  if (nEW != 99 || nQCD != 99) KillHigherOrders(first_amp);
 
  msg.Debugging()<<"Now Amplitude_Generator::Compare"<<endl;

  Compare(first_amp);

  // count physical number of 4 Verticies (nampl- n4/3*2 -n44/9*8)
  msg.Out()<<" real number of amplitudes "<<CountRealAmplitudes(first_amp)<<endl<<endl;
  
#ifdef _USE_MPI_
  if (rank==0) {
    for (int cpu=1;cpu<size;++cpu) {
      int dummy;
      MPI::COMM_WORLD.Recv(&dummy,1,MPI::INT, cpu,42);
    }
    for (int cpu=1;cpu<size;++cpu) {
      int dummy;
      MPI::COMM_WORLD.Recv(&dummy,1,MPI::INT, cpu,43);
    }
  }
  else {
    int dummy;
    MPI::COMM_WORLD.Send(&dummy,1,MPI::INT, 0,42);
    MPI::COMM_WORLD.Send(&dummy,1,MPI::INT, 0,43);
  }
#endif

  Amplitude_Manipulator(N,fl,b).FixSign(first_amp);
  
  return first_amp;
}

void Amplitude_Generator::CommitMPITypes() {
#ifdef _USE_MPI_
  MPI_Lorentz_Function   sample_lf;

  MPI::Datatype lf_types[2]  = {MPI::INT,MPI::INT};
  int           lf_blocks[2] = {1       ,4       };
  MPI::Aint     lf_disp[2];
  lf_disp[0] = MPI::Get_address(&sample_lf.m_type);
  lf_disp[1] = MPI::Get_address(&sample_lf.m_partarg);
  for (short int i=1;i>=0;i--) lf_disp[i] -= lf_disp[0];

  mpi_lf_type = MPI::Datatype::Create_struct(2,lf_blocks,lf_disp,lf_types);
  mpi_lf_type.Commit();

  // ----------------------------------------

  MPI_Color_Function   sample_cf;

  MPI::Datatype cf_types[3]  = {MPI::INT,MPI::INT,MPI::CHAR};
  int           cf_blocks[3] = {1       ,3       ,4        };
  MPI::Aint     cf_disp[3];
  cf_disp[0] = MPI::Get_address(&sample_cf.m_type);
  cf_disp[1] = MPI::Get_address(&sample_cf.m_partarg);
  cf_disp[2] = MPI::Get_address(&sample_cf.m_strarg);
  for (short int i=2;i>=0;i--) cf_disp[i] -= cf_disp[0];

  mpi_cf_type = MPI::Datatype::Create_struct(3,cf_blocks,cf_disp,cf_types);
  mpi_cf_type.Commit();


  // ----------------------------------------

  MPI::Datatype sv_types[1]  = {MPI::INT};
  int           sv_blocks[1] = {4       };
  MPI::Aint     sv_disp[1];
  sv_disp[0] =  0;

  mpi_sv_type = MPI::Datatype::Create_struct(1,sv_blocks,sv_disp,sv_types);
  mpi_sv_type.Commit();

  // ----------------------------------------


  MPI_Point   sample_point;

  MPI::Datatype point_types[6]  = {MPI::INT,mpi_lf_type,mpi_cf_type,mpi_sv_type,MPI::DOUBLE,MPI::INT };
  int           point_blocks[6] = {1       ,1          ,1          ,1          ,8          ,5        };
  MPI::Aint     point_disp[6];
  point_disp[0] = MPI::Get_address(&(sample_point.m_fl));
  point_disp[1] = MPI::Get_address(&(sample_point.m_lf));
  point_disp[2] = MPI::Get_address(&(sample_point.m_cf));
  point_disp[3] = MPI::Get_address(&(sample_point.m_v));
  point_disp[4] = MPI::Get_address(&(sample_point.m_cpl));
  point_disp[5] = MPI::Get_address(&(sample_point.m_left));


  for (short int i=5;i>=0;i--) point_disp[i] -= point_disp[0];
  
  mpi_point_type = MPI::Datatype::Create_struct(6,point_blocks,point_disp,point_types);
  mpi_point_type.Commit();
#endif

}



















