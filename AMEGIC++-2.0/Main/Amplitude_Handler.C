#include "Amplitude_Handler.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Color_Group.H"
#include <iostream>
#include <stdio.h>
#include "prof.hh"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Amplitude_Handler::Amplitude_Handler(int N,Flavour* fl,int* b,Polarisation* pol,
				     Topology* top,Basic_Sfuncs* BS,
				     String_Handler* _shand,
				     std::string pID) : shand(_shand)
{
  groupname = string("All Amplitudes");

  gen = new Amplitude_Generator(N,fl,b,top,BS,shand);
  Single_Amplitude* firstgraph = gen->Matching();
  delete gen;

  Single_Amplitude* n;
  n = firstgraph;

  ngraph = 0;

  while (n) {
    ++ngraph;
    n->Zprojecting(fl,ngraph);
    //n->FillCoupling(shand); 

    if (n->on) {
      pol->Replace_Numbers(N,fl,n);
      BS->Build_Momlist(*n->GetPlist());
    } 
    n = n->Next;
  }
  
  CheckEqual(firstgraph);

  msg.Debugging()<<ngraph<<" Graph(s) found"<<endl;  

  if (ngraph==0) {
    msg.Error()<<"No Graph found for ";
    for (short int i=0;i<N;i++) msg.Error()<<fl[i]<<";";
    msg.Error()<<endl;
    return;
  }
  //Colors
  
  //* FK *  CFCol_Matrix   = new CFColor(N,firstgraph);

  CFCol_Matrix   = new CFColor(N,firstgraph,pID);


  /*
  int* colgroup = new int[ngraph];
  for (short int i=0;i<ngraph;i++) colgroup[i] = -1;
  int colcount = 0;

  for (short int i=0;i<ngraph;i++) {
    if (colgroup[i]==-1) {
      colgroup[i] = colcount;colcount++;
      for (short int j=i+1;j<ngraph;j++) {
	int hit = 1;
	for (short int k=0;k<ngraph;k++) {
	  if (CFCol_Matrix->Mij(i,k)!=CFCol_Matrix->Mij(j,k)) {
	    hit = 0;
	    break;
	  }
	}
	if (hit) colgroup[j] = colgroup[i];
      }
    }
  }
  cout<<CFCol_Matrix->MatrixSize()<<" different color group(s) found!"<<endl;
  */

  // dafuer :  ncount
  //           map    == colgroup[i];
  //        !! sign   == sign(id[i])    
  //           colfactors == CFC

  /*
  colfactors = new Complex*[colcount];
  for (int i=0;i<colcount;i++) colfactors[i] = new Complex[colcount];

  for (short int i=0;i<ngraph;i++) {
    for (short int j=0;j<ngraph;j++) {
      colfactors[colgroup[i]][colgroup[j]] = CFCol_Matrix->Mij(i,j);
    }
  }
  */


  for (int i=0;i<CFCol_Matrix->MatrixSize();i++) graphs.push_back(new Color_Group());

  //On-Switches
  int* switch_graphs = new int[ngraph];
  // *FS*  for(short int i=0;i<ngraph;i++) onswitch[i] = 1;
  Kicker(switch_graphs,ngraph,pID);

  n = firstgraph;

  // fill color groups
  int ncount = 0;
  while (n) {
    //Kicker does not work properly right now!!!!
    //if (switch_graphs[ncount])
    pointlist.push_back(n->GetPointlist()); 
    graphs[CFCol_Matrix->CFMap(ncount)]->Add(n,CFCol_Matrix->CFSign(ncount));
    //    graphs[colgroup[ncount]]->Add(n);
    n = n->Next;
    ncount++;	   
  }
  
  delete[] switch_graphs;
  
  // create Superamplitudes:
  for (int i=0;i<graphs.size();i++) graphs[i]->BuildGlobalString(b,N,BS,fl,shand);
  
  int dummy = 0;
  SetNumber(dummy);  
  namplitude = dummy;

  if (rpa.gen.Tracking()) {
    PrintGraph();
    BS->Print_Momlist();
  }

  CheckEqualInGroup();
  
  //Probabilities
  sw_probabs = 0;

  probs = 0;
  
  probabs = new double[graphs.size()];
  Mi      = new Complex[graphs.size()];
}


Amplitude_Handler::Amplitude_Handler(int N,Flavour* fl,int* b,Polarisation* pol,
				     Topology* top,Basic_Sfuncs* BS,
				     String_Handler* _shand) : shand(_shand)
{
  groupname = string("All Amplitudes");

  gen = new Amplitude_Generator(N,fl,b,top,BS,shand);
  Single_Amplitude* firstgraph = gen->Matching();
  delete gen;

  Single_Amplitude* n;
  n = firstgraph;

  ngraph = 0;

  while (n) {
    ++ngraph;
    n->Zprojecting(fl,ngraph);
    //n->FillCoupling(shand); 

    if (n->on) {
      pol->Replace_Numbers(N,fl,n);
      BS->Build_Momlist(*n->GetPlist());
    } 
    n = n->Next;
  }
  
  CheckEqual(firstgraph);

  msg.Debugging()<<ngraph<<" Graph(s) found"<<endl;  

  if (ngraph==0) {
    msg.Error()<<"No Graph found for ";
    for (short int i=0;i<N;i++) msg.Error()<<fl[i]<<";";
    msg.Error()<<endl;
    return;
  }
  //Colors
  
  CFCol_Matrix   = new CFColor(N,firstgraph);
  // *FK*   CFCol_Matrix   = new CFColor(N,firstgraph,pID);


  /*
  int* colgroup = new int[ngraph];
  for (short int i=0;i<ngraph;i++) colgroup[i] = -1;
  int colcount = 0;

  for (short int i=0;i<ngraph;i++) {
    if (colgroup[i]==-1) {
      colgroup[i] = colcount;colcount++;
      for (short int j=i+1;j<ngraph;j++) {
	int hit = 1;
	for (short int k=0;k<ngraph;k++) {
	  if (CFCol_Matrix->Mij(i,k)!=CFCol_Matrix->Mij(j,k)) {
	    hit = 0;
	    break;
	  }
	}
	if (hit) colgroup[j] = colgroup[i];
      }
    }
  }
  cout<<CFCol_Matrix->MatrixSize()<<" different color group(s) found!"<<endl;
  */

  // dafuer :  ncount
  //           map    == colgroup[i];
  //        !! sign   == sign(id[i])    
  //           colfactors == CFC

  /*
  colfactors = new Complex*[colcount];
  for (int i=0;i<colcount;i++) colfactors[i] = new Complex[colcount];

  for (short int i=0;i<ngraph;i++) {
    for (short int j=0;j<ngraph;j++) {
      colfactors[colgroup[i]][colgroup[j]] = CFCol_Matrix->Mij(i,j);
    }
  }
  */


  for (int i=0;i<CFCol_Matrix->MatrixSize();i++) graphs.push_back(new Color_Group());

  //On-Switches
  int* switch_graphs = new int[ngraph];
  for(short int i=0;i<ngraph;i++) switch_graphs[i] = 1;
  // *FS*  Kicker(switch_graphs,ngraph,pID);

  n = firstgraph;

  // fill color groups
  int ncount = 0;
  while (n) {
    //Kicker does not work properly right now!!!!
    //if (switch_graphs[ncount])
    pointlist.push_back(n->GetPointlist()); 
    graphs[CFCol_Matrix->CFMap(ncount)]->Add(n,CFCol_Matrix->CFSign(ncount));
    //    graphs[colgroup[ncount]]->Add(n);
    n = n->Next;
    ncount++;	   
  }
  
  delete[] switch_graphs;
  //  delete[] colgroup;
  
  // create Superamplitudes:
  //for (int i=0;i<graphs.size();i++) graphs[i]->BuildGlobalString(b,N,BS,fl,shand);
  
  int dummy = 0;
  SetNumber(dummy);  
  namplitude = dummy;
  PrintGraph();

  BS->Print_Momlist();

  CheckEqualInGroup();
  
  //Probabilities
  sw_probabs = 0;

  probs = 0;
  
  probabs = new double[graphs.size()];
  Mi      = new Complex[graphs.size()];
}


Amplitude_Handler::~Amplitude_Handler() 
{
  if (ngraph>0) {
    delete CFCol_Matrix;
    delete[] probabs;
    delete[] Mi;

    //    for (int i=0;i<graphs.size();i++) delete[] colfactors[i];
    //    delete[] colfactors;
  }
  //Single_Amplitude's

  for (int i=0;i<graphs.size();i++) delete graphs[i];

  
}

int Amplitude_Handler::PropProject(Amplitude_Base* f,int zarg)
{
  if (zarg<100) return zarg;

  list<Pfunc*>* pl = f->GetPlist();
  for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
    if ((*pit)->arg[0]==iabs(zarg)) return (*pit)->momnum; 
  }  
  msg.Error()<<"Bug in Amplitude_Handler::PropProject()"<<endl;
  abort();
  return 0;
}


int Amplitude_Handler::CompareZfunc(Amplitude_Base* f1,Zfunc* z1,Amplitude_Base* f2,Zfunc* z2)
{
  if (z1->type!=z2->type) return 0;
  
  if (z1->pn!=z2->pn) return 0;
  
  //Arguments
  for (short int i=0;i<z1->narg;i++) {
    if (PropProject(f1,z1->arg[i])!=PropProject(f2,z2->arg[i])) return 0;
  }

  //couplings
  for (short int i=0;i<z1->ncoupl;i++) {
    if (z1->coupl[i]!=z2->coupl[i]) return 0;
  }

  //Propagators
  for (short int i=0;i<z1->pn;i++) {
    if (PropProject(f1,z1->psnew[i].numb)!=PropProject(f2,z2->psnew[i].numb)) return 0;
    //Flavour of props
    if (iabs(z1->psnew[i].numb)>99) {
      
      Flavour flav1;
      list<Pfunc*>* pl = f1->GetPlist();
      for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc* p = *pit;
	if (p->arg[0]==iabs(z1->psnew[i].numb)) {
	  flav1 = p->fl;
	  break;
	}
      }

      Flavour flav2;
      pl = f2->GetPlist();
      for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc* p = *pit;
	if (p->arg[0]==iabs(z2->psnew[i].numb)) {
	  flav2 = p->fl;
	  break;
	}
      }
      if (flav1!=flav2) return 0;
    }
  }
}

string IString(int i)
{
  std::strstream sstream;
  sstream<<i;
  string istr;
  sstream>>istr;
  return istr;
}

void Amplitude_Handler::CheckEqual(Single_Amplitude* firstgraph)
{
  Single_Amplitude* f1 = firstgraph;
  Single_Amplitude* f2;

  int count  = 0;
  int zcount = 0;

  int g1 = 0;

  int basiczcount = 0;


  while (f1) {
    list<Zfunc*>* zlist = f1->GetZlist();
    for (list<Zfunc*>::iterator zit=zlist->begin();zit!=zlist->end();++zit) (*zit)->Equal = *zit;
    f1 = f1->Next;
  }

  f1 = firstgraph;

  while (f1) {
    list<Zfunc*>* zlist = f1->GetZlist();
    int cz1 = 0;
    for (list<Zfunc*>::iterator zit=zlist->begin();zit!=zlist->end();++zit) {
      Zfunc* z1 = (*zit);
      zcount++;
      if (z1->Equal==z1) {
	//setting the string of it
	z1->str = string("Z")+IString(basiczcount);basiczcount++;
        f2 = f1->Next;
        int g2 = g1+1;
        while (f2) {
          list<Zfunc*>* zlist2 = f2->GetZlist();
          int cz2 = 0;
          for (list<Zfunc*>::iterator zit2=zlist2->begin();zit2!=zlist2->end();++zit2) {
            Zfunc* z2 = (*zit2);
            if (z2->Equal==z2) {
              if (CompareZfunc(f1,z1,f2,z2)) {
                z2->Equal = z1;
		z2->str   = z1->str;
                count++;
              }
            }
            cz2++;
          }
          f2 = f2->Next;g2++;
        }
      }
      cz1++;
    }
    f1 = f1->Next;g1++;
  }
  msg.Debugging()<<"Equal Zfuncs: "<<count<<"/"<<zcount<<endl;
}

void Amplitude_Handler::CheckEqualInGroup()
{
  //Renew all Zfuncs
  for (int g1=0;g1<namplitude;g1++) {
    Amplitude_Base* f1 = GetAmplitude(g1);    
    list<Zfunc*>* zlist = f1->GetZlist();
    for (list<Zfunc*>::iterator zit=zlist->begin();zit!=zlist->end();++zit) {
      Zfunc* zl1 = (*zit);
      for (int i=0;i<zl1->GetSize();i++) (*zl1)[i]->Equal = (*zl1)[i];
    }
  }

  int count  = 0;
  int zcount = 0;

  for (int g1=0;g1<namplitude;g1++) {
    Amplitude_Base* f1 = GetAmplitude(g1);    
    list<Zfunc*>* zlist = f1->GetZlist();
    int cz1 = 0;
    for (list<Zfunc*>::iterator zit=zlist->begin();zit!=zlist->end();++zit) {
      Zfunc* zl1 = (*zit);
      for (int i=0;i<zl1->GetSize();i++) {
	Zfunc* z1 = (*zl1)[i];
	zcount++;
	if (z1->Equal==z1) {
	  for (int g2=g1+1;g2<namplitude;g2++) {
	    Amplitude_Base* f2   = GetAmplitude(g2);    
	    list<Zfunc*>* zlist2 = f2->GetZlist();
	    int cz2 = 0;
	    for (list<Zfunc*>::iterator zit=zlist2->begin();zit!=zlist2->end();++zit) {
	      Zfunc* zl2 = (*zit);
	      for (int j=0;j<zl2->GetSize();j++) {
		Zfunc* z2 = (*zl2)[j];
		if (z2->Equal==z2) {
		  if (CompareZfunc(f1,z1,f2,z2)) {
		    z2->Equal = z1;
		    count++;
		  }
		}
		cz2++;
	      }
	    }
          }
        }
	cz1++;
      }
    }
  }
}


bool Amplitude_Handler::ExistFourVertex(Point* p)
{
  if (p==0)      return 0;
  if (p->middle) return 1;
  
  if (ExistFourVertex(p->left)) return 1;
  return ExistFourVertex(p->right);
}

void Amplitude_Handler::Kicker(int* Switch_Vector,int ngraph,std::string pID)
{
  char name[100];
  sprintf(name,"%s.kick",(string("Process/")+pID).c_str());
  //  sprintf(name,"%s/Kicker.dat",(string("Process/")+pID).c_str());
  fstream test;
  test.open(name,ios::in);

  if (test) {
    test.close();
    fstream from;
    from.open(name,ios::in);
    //read in
    int i,sw;
    for(;from;) {
      from>>i>>sw;
      Switch_Vector[i-1] =sw;
      if (sw==0) msg.Tracking()<<"Diagram "<<i<<" kicked!"<<endl;
      if (i==ngraph) break;
    }
    msg.Tracking()<<"File "<<name<<" read."<<endl;  
    return;
  }

  test.close();

  msg.Tracking()<<"File "<<name<<" not found."<<endl;  

  for(short int i=0;i<ngraph;i++) Switch_Vector[i] = 1;
  

  ofstream to;
  to.open(name);

  for(short int i=0;i<ngraph;i++) to<<i+1<<"     "<<Switch_Vector[i]<<endl;
  
  msg.Tracking()<<"File "<<name<<" saved."<<endl;  
}

Point* Amplitude_Handler::GetPointlist(int n)
{ return pointlist[n];}

void Amplitude_Handler::Reset_ProbAbs()
{
  short int i;
  for (i=0;i<graphs.size();i++) probabs[i] = 0.;
}

double Amplitude_Handler::Get_Probab(int i) {return probabs[i];}


Complex Amplitude_Handler::Zvalue(String_Handler * sh, int ihel)
{
  for (int i=0;i<graphs.size();i++) Mi[i] = graphs[i]->Zvalue(sh, ihel);

  Complex M(0.,0.);
  for (short int i=0;i<graphs.size();i++) {
    for (short int j=0;j<graphs.size();j++) {
      M += Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
    }
  }

  return M;
}

Complex Amplitude_Handler::Zvalue(int ihel,int* sign)
{
  for (int i=0;i<graphs.size();i++) Mi[i] = graphs[i]->Zvalue(ihel,sign);

  Complex M(0.,0.);
  for (short int i=0;i<graphs.size();i++) {
    for (short int j=0;j<graphs.size();j++) {
      M += Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
    }
  }

  return M;
}

	/*
	for (short int i=0;i<ngraph;i++) {
	  for (short int j=0;j<ngraph;j++) M += Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);
	  if (probs) {
	    if (sw_probabs==2) probabs[i] += abs(Mi[i]*conj(Mi[i])*CFCol_Matrix->Mij(i,i));
	    Complex c(0.,0.);
	    if (sw_probabs==3) {
	      for (short int j=0;j<ngraph;j++)
		c += Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);
	    }
	    probabs[i] += abs(c);
	  }
	}
	*/











