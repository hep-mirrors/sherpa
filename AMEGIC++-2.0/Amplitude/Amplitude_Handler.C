#include "Amplitude_Handler.H"
#include "Zfunc_Calc.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Color_Group.H"
#include <iostream>
#include <stdio.h>
#include "prof.hh"

//#include "MyTiming.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

void out_pfunc(Pfunc & pf) {
  pf.operator<<(cout);
}


Amplitude_Handler::Amplitude_Handler(int N,Flavour* fl,int* b,Polarisation* pol,
				     Interaction_Model_Base * model,Topology* top,
				     int _orderQCD,int _orderEW,
				     Basic_Sfuncs* BS,String_Handler* _shand,
				     std::string pID) : shand(_shand)
{
  groupname = string("All Amplitudes");

  gen = new Amplitude_Generator(N,fl,b,model,top,_orderQCD,_orderEW,BS,shand);
  Single_Amplitude* firstgraph = gen->Matching();
  delete gen;

  Single_Amplitude* n;
  n = firstgraph;
  ngraph = 0;

  while (n){ 
    ++ngraph;
    n->Zprojecting(fl,ngraph);
    //n->FillCoupling(shand); 

    if (n->on) {
      pol->Replace_Numbers(N,fl,n);

      OptimizeProps(N,n);

      BS->BuildMomlist(*n->GetPlist());
    } 
    n = n->Next;
  }

  PreCluster(firstgraph); 
  CheckEqual(firstgraph);

  msg.Debugging()<<ngraph<<" Graph(s) found"<<endl;  

  if (ngraph==0) {
    msg.Error()<<"No Graph found for ";
    for (short int i=0;i<N;i++) msg.Error()<<fl[i]<<";";
    msg.Error()<<endl;
    return;
  }
  //Colors
  
  CFCol_Matrix   = new CFColor(N,firstgraph,pID);

  for (int i=0;i<CFCol_Matrix->MatrixSize();i++) graphs.push_back(new Color_Group());

  //On-Switches
  //int* switch_graphs = new int[ngraph];
  // *FS*  for(short int i=0;i<ngraph;i++) onswitch[i] = 1;
  //Kicker(switch_graphs,ngraph,pID);

  n = firstgraph;

  // fill color groups
  int ncount = 0;
  while (n) {
    while(TOrder(n)>0)n=n->Next;
    //Kicker does not work properly right now!!!!
    //if (switch_graphs[ncount])
    pointlist.push_back(n->GetPointlist()); 
    graphs[CFCol_Matrix->CFMap(ncount)]->Add(n,CFCol_Matrix->CFSign(ncount));
    //    graphs[colgroup[ncount]]->Add(n);
    n = n->Next;
    ncount++;	   
  }
  
  //delete[] switch_graphs;
  
  ngraph = ncount;
   
  for (int i=0;i<graphs.size();i++) graphs[i]->BuildGlobalString(b,N,BS,fl,shand);

  int dummy = 0;
  SetNumber(dummy);  
  namplitude = dummy;

  if (rpa.gen.Tracking()) {
    PrintGraph();
    //BS->PrintMomlist();
  }

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

  Pfunc_List* pl = f->GetPlist();
  for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
    if ((*pit)->arg[0]==iabs(zarg)) return (*pit)->momnum; 
  }  
  msg.Error()<<"Bug in Amplitude_Handler::PropProject()"<<endl;
  abort();
  return 0;
}


int Amplitude_Handler::CompareZfunc(Amplitude_Base* f1,Zfunc* z1,Amplitude_Base* f2,Zfunc* z2)
{
  if (z1->GetSize()!=z2->GetSize()) return 0;

  if (z1->GetSize()>1){
    for(int i=0;i<z1->GetSize();i++)
      if ( CompareZfunc(f1,(*z1)[i],f2,(*z2)[i])==0 )return 0;
    return 1;
  }

  if (z1->m_type!=z2->m_type) return 0;
  
  if (z1->m_nprop!=z2->m_nprop) return 0;
  
  //Arguments
  for (short int i=0;i<z1->m_narg;i++) {
    if (PropProject(f1,z1->p_arguments[i])!=PropProject(f2,z2->p_arguments[i])) return 0;
  }

  //couplings
  for (short int i=0;i<z1->m_ncoupl;i++) {
    if (z1->p_couplings[i]!=z2->p_couplings[i]) return 0;
  }

  //Propagators
  for (short int i=0;i<z1->m_nprop;i++) {
    if (PropProject(f1,z1->p_propagators[i].numb)!=PropProject(f2,z2->p_propagators[i].numb)) return 0;
    //Flavour of props
    if (iabs(z1->p_propagators[i].numb)>99) {
      
      Flavour flav1;
      Pfunc_List* pl = f1->GetPlist();
      for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc* p = *pit;
	if (p->arg[0]==iabs(z1->p_propagators[i].numb)) {
	  flav1 = p->fl;
	  break;
	}
      }

      Flavour flav2;
      pl = f2->GetPlist();
      for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc* p = *pit;
	if (p->arg[0]==iabs(z2->p_propagators[i].numb)) {
	  flav2 = p->fl;
	  break;
	}
      }
      if (flav1!=flav2) return 0;
    }
  }
  return 1;
}

string IString(int i)
{
  std::strstream sstream;
  sstream<<i;
  string istr;
  sstream>>istr;
  return istr;
}

void Amplitude_Handler::OptimizeProps(int N,Single_Amplitude* f1)
{
  Zfunc_List* zlist = f1->GetZlist();
  Pfunc_List* pl = f1->GetPlist();
  for(Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit){
    if((*pit)->argnum > (N/2)+1){
      //cout<<"OptimizeProps: ";out_pfunc(*(*pit));
      int nargnum=N+2-(*pit)->argnum;
      int* arg= new int[nargnum];
      int cnt=1,hit;
      arg[0]=(*pit)->arg[0];
      for(int i=0;i<N;i++){
	hit=0;
	for(int j=1;j<(*pit)->argnum;j++)if(i==(*pit)->arg[j])hit=1;
	if(hit==0){
	  arg[cnt]=i;
	  cnt++;
	}
      }
      (*pit)->argnum=nargnum;
      delete[] (*pit)->arg;
      (*pit)->arg=new int[nargnum];
      for(int j=0;j<(*pit)->argnum;j++)(*pit)->arg[j]=arg[j];
      delete[] arg;
      //cout<<"after: ";out_pfunc(*(*pit));
      if((*pit)->fl.IsFermion()){
	f1->SetSign(-(f1->GetSign()));
	(*pit)->fl=(*pit)->fl.Bar();
      }
      for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit){
	for(int j=0;j<(*zit)->m_nprop;j++)if((*zit)->p_propagators[j].numb==(*pit)->arg[0]){
	  if((*zit)->p_propagators[j].direction==Direction::Incoming)
	       (*zit)->p_propagators[j].direction=Direction::Outgoing;
	  else (*zit)->p_propagators[j].direction=Direction::Incoming;
	}
      }
    }
  }
}

void Amplitude_Handler::PreCluster(Single_Amplitude* firstgraph)
{
  int cnt=0;
  Single_Amplitude* f1 = firstgraph;
  while (f1) {
    Zfunc_List* zlist = f1->GetZlist();
    vector<int> propselect;
    Pfunc_List* pl = f1->GetPlist();
    for(Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit)
      if((*pit)->fl.Kfcode()==kf::photon||(*pit)->fl.Kfcode()==kf::Z) 
	propselect.push_back((*pit)->arg[0]);
    for(Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit)
      if((*pit)->fl.IsScalar()) 
	propselect.push_back((*pit)->arg[0]);

    for(int i=0;i<propselect.size();i++){
      int ia=0;
      Zfunc *zh[2];
      for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit){
	for(int j=0;j<(*zit)->m_narg;j++){
	  if((*zit)->p_arguments[j]==propselect[i]){
	    zh[ia]=(*zit);
	    ia++;
	    break;
	  }
	}
	if(ia==2)break;
      }
      
      if(ia==2) if(zh[0]->GetSize()==1 && zh[1]->GetSize()==1){
	int n0=0;
	for(int j=0;j<zh[0]->m_narg;j++)if(zh[0]->p_arguments[j]>99)n0++;
	int n1=0;
	for(int j=0;j<zh[1]->m_narg;j++)if(zh[1]->p_arguments[j]>99)n1++;
	if (n0<=1 || n1<=1 || 
	   ( (zh[0]->m_type==zl::Y || zh[0]->m_type==zl::Z) && 
	     (zh[0]->m_type==zl::Y || zh[0]->m_type==zl::Z)    ) ) {

	  // cout<<"Precluster : Amplitude "<<cnt<<endl;
	  //zh[0]->Print();
	  //zh[1]->Print();
	  Zfunc *zh0,*zh1;
	  
	  zh0=zh[0];zh1=zh[1];   //unique order of zfunctions
	  if(zh[0]->m_narg>zh[1]->m_narg){zh0=zh[1];zh1=zh[0];}
	  if(zh[0]->m_narg==zh[1]->m_narg){
	    //cout<<zh[0]->m_narg<<" "<<zh[1]->m_narg<<endl;
	    for(int j=0;j<zh[0]->m_narg;j++){
	      if(PropProject(f1,zh[0]->p_arguments[j])<PropProject(f1,zh[1]->p_arguments[j]))break;
	      if(PropProject(f1,zh[0]->p_arguments[j])>PropProject(f1,zh[1]->p_arguments[j])){zh0=zh[1];zh1=zh[0];break;}
	    }
	  }
	  //cout<<"before erase"<<endl;
	  /*for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit)
	    cout<<(*zit)<<endl;
	    cout<<"erase: "<<zh0<<" "<<zh1<<endl;*/
	  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();){
	    //cout<<(*zit)<<endl;
	    if((*zit)==zh0 || (*zit)==zh1){ zit=zlist->erase(zit);}
	    else zit++;
	  }
	  Zfunc_Group *sf=new Zfunc_Group(*zh0,*zh1,propselect[i],pl);
	  zlist->push_back(sf);
	  sf->Print();
	  /*cout<<endl;
	  zlist=f1->GetZlist();
	  for (list<Zfunc*>::iterator zit=zlist->begin();zit!=zlist->end();++zit)
	  cout<<(*zit)<<endl;*/
	  
	}
      }
    }
    cnt++;
    f1 = f1->Next;
  }
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
    Zfunc_List* zlist = f1->GetZlist();
    for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) (*zit)->p_equal = *zit;
    f1 = f1->Next;
  }

  f1 = firstgraph;

  while (f1) {
    Zfunc_List* zlist = f1->GetZlist();
    int cz1 = 0;
    for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
      Zfunc* z1 = (*zit);
      zcount++;
      if (z1->p_equal==z1) {
	//setting the string of it
	z1->m_str = string("Z")+IString(basiczcount);basiczcount++;
        f2 = f1->Next;
        int g2 = g1+1;
        while (f2) {
          Zfunc_List* zlist2 = f2->GetZlist();
          int cz2 = 0;
          for (Zfunc_Iterator zit2=zlist2->begin();zit2!=zlist2->end();++zit2) {
            Zfunc* z2 = (*zit2);
            if (z2->p_equal==z2) {
              if (CompareZfunc(f1,z1,f2,z2)) {
                z2->p_equal = z1;
		z2->m_str   = z1->m_str;
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
  // === new ========================================
  //
  // * build list of zfuncs and owner graphs
  // (* while building check for existence)

  // === old ========================================
  //Renew all Zfuncs
  for (int g1=0;g1<namplitude;g1++) {
    Amplitude_Base* f1 = GetAmplitude(g1);    
    Zfunc_List* zlist = f1->GetZlist();
    for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
      Zfunc* zl1 = (*zit);
      zl1->p_equal = zl1;
      for (int i=0;i<zl1->GetSize();i++) (*zl1)[i]->p_equal = (*zl1)[i];
    }
  }

  int count  = 0;
  int zcount = 0;

  for (int g1=0;g1<namplitude;g1++) {
    Amplitude_Base* f1 = GetAmplitude(g1);    
    Zfunc_List* zlist = f1->GetZlist();
    int cz1 = 0;
    for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
      Zfunc* zl1 = (*zit);
      for (int i=0;i<zl1->GetSize();i++) {
	Zfunc* z1 = (*zl1)[i];
	zcount++;
	if (z1->p_equal==z1) {
	  for (int g2=g1+1;g2<namplitude;g2++) {
	    Amplitude_Base* f2   = GetAmplitude(g2);    
	    Zfunc_List* zlist2 = f2->GetZlist();
	    int cz2 = 0;
	    for (Zfunc_Iterator zit=zlist2->begin();zit!=zlist2->end();++zit) {
	      Zfunc* zl2 = (*zit);
	      for (int j=0;j<zl2->GetSize();j++) {
		Zfunc* z2 = (*zl2)[j];
		if (z2->p_equal==z2) {
		  if (CompareZfunc(f1,z1,f2,z2)) {
		    z2->p_equal = z1;
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
  msg.Tracking()<<"Equal Zfuncs: "<<count<<"/"<<zcount<<endl;
}


bool Amplitude_Handler::ExistFourVertex(Point* p)
{
  if (p==0)      return 0;
  if (p->middle) return 1;
  
  bool sw = 0;
  sw = ExistFourVertex(p->left);
  if (sw) return 1;
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
  for (int i=0;i<graphs.size();i++) {
    Mi[i] = graphs[i]->Zvalue(ihel,sign);
  }

  Complex M(0.,0.);
  for (short int i=0;i<graphs.size();i++) {
    for (short int j=0;j<graphs.size();j++) {
      M += Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
    }
  }
   
  return M;
}

int Amplitude_Handler::TOrder(Single_Amplitude* a)
{  
  if(rpa.me.Model()!=Model_Type::ADD) return 0;
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());
  int maxorder = dr.GetValue<int>("Max_KK-Props");

  int cnt=0;
  Pfunc_List* pl = a->GetPlist();
  for(Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit)
    if((*pit)->fl.IsKK())cnt++;
  if (cnt>maxorder) return 1;
  return 0;
} 









