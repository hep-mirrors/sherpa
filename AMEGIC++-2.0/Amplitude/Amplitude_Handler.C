#include "Amplitude_Handler.H"
#include "Amplitude_Output.H"
#include "Zfunc_Calc.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Amplitude_Generator.H"
#include "Amplitude_Manipulator.H"
#include "Color_Group.H"
#include <iostream>
#include <stdio.h>
#include "prof.hh"
#include "MyStrStream.H"
#include "Process_Info.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Amplitude_Handler::Amplitude_Handler(int N,Flavour* fl,int* b,Process_Info* pinfo,
				     Interaction_Model_Base * model,Topology* top,
				     int & _orderQCD,int & _orderEW,Basic_Sfuncs* BS,
				     String_Handler* _shand, bool print_graph,bool create_4V) 
  : shand(_shand),CFCol_Matrix(0),probabs(0),Mi(0), m_print_graph(print_graph)
{
  int ndecays=pinfo->Ndecays();
  int nm = pinfo->Nmax(0);
  int nin = 1;
  if (b[1]==-1) nin=2;
  int * b_dec = new int[nm];
  b_dec[0] = -1;
  for (int i=1;i<nm;i++) b_dec[i] = 1;

  Single_Amplitude** subgraphlist = new Single_Amplitude*[ndecays+1];
  Amplitude_Generator * gen; 

  Flavour *sfl;
  if (ndecays>0) {
    sfl = new Flavour[pinfo->Nmax(nin)];
    sfl[0] = fl[0];
    sfl[1] = fl[1];
    pinfo->GetFlavList(sfl+nin);
  }
  else sfl=fl;

  //core process
  gen = new Amplitude_Generator(nin+pinfo->Nout(),sfl,b,model,top,_orderQCD,_orderEW,BS,shand,create_4V);
  subgraphlist[0] = gen->Matching();
  gen->GetOrders(_orderEW,_orderQCD);
  delete gen;

  //decay processes
  for (int i=1;i<=ndecays;i++) {
    int j=i;
    Process_Info *pi=pinfo->GetDecay(j);
//     pi->Print();cout<<endl;
    sfl[0] = *(pi->p_fl);
    pi->GetFlavList(sfl+1);
    gen = new Amplitude_Generator(1+pi->Nout(),sfl,b_dec,model,top,99,99,BS,shand);
    subgraphlist[i] = gen->Matching();
    if (subgraphlist[i]==NULL) {
      ndecays = 0;
      subgraphlist[0] = NULL;
    }
    int ew,qcd;
    gen->GetOrders(ew,qcd);
    _orderEW  += ew;
    _orderQCD += qcd;
    delete gen;
  }

  if (msg.LevelIsTracking()) {
    msg.Out()<<"Amplitude_Handler::Amplitude_Handler:"<<endl;
    int f=1;
    for(int i=0;i<ndecays+1;i++) {
      int j=0;
      Single_Amplitude* nn = subgraphlist[i];
      while (nn){ 
	++j;
	nn = nn->Next;
      }
      msg.Out()<<"Process "<<i;
      if (i==0)msg.Out()<<" (core)";
      else msg.Out()<<" (decay)";
      msg.Out()<<" has "<<j<<" Amplitudes"<<endl;
      f*=j;
    }
    msg.Out()<<"Total: "<<f<<" Amplitudes"<<endl;
  }

  if (ndecays==0 || subgraphlist[0]==0) firstgraph = subgraphlist[0];
  else ConstructSignalAmplitudes(N,fl,b,pinfo,subgraphlist,BS);

  Single_Amplitude* n = firstgraph;
  ntotal = 0;
  while (n){ 
    ++ntotal;
    n = n->Next;
  }
  msg.Tracking()<<"Total number of Amplitudes "<<ntotal<<endl;
  ngraph = ntotal;
  delete [] b_dec;
}

void Amplitude_Handler::ConstructSignalAmplitudes(int N,Flavour* fl,int* b,
						  Process_Info* pinfo,Single_Amplitude** sglist,
						  Basic_Sfuncs* BS)
{
  int ndecays=pinfo->Ndecays();
  firstgraph = NULL;
  Single_Amplitude *n=NULL,*next;
  Single_Amplitude** nl = new Single_Amplitude*[ndecays+1];
  for (int i=0;i<ndecays+1;i++) nl[i] = sglist[i];
  int over = 0;
  for (;;) {
    next = new Single_Amplitude(b,N,pinfo,nl,BS,fl,shand);
    if (n) n->Next=next;
    n = next;
    if (!firstgraph) firstgraph = n;

    for (int i=ndecays;i>=0;i--) {
      nl[i] = nl[i]->Next;
      if (nl[i]) break;
      nl[i] = sglist[i];
      if (i==0) over = 1;
    }
    if (over) break;
  }
  
  Amplitude_Manipulator(N,fl,b,ndecays).FixSign(firstgraph);

  delete [] nl;
  for (int i=0;i<ndecays+1;i++) {
    n = sglist[i];
    while (n) {
      next = n->Next;
      delete n;
      n = next;
    }
  }
  delete [] sglist;
}

void Amplitude_Handler::CompleteAmplitudes(int N,Flavour* fl,int* b,Polarisation* pol,
					   Topology* top,Basic_Sfuncs* BS,std::string pID)
{
  bool gen_colors=true;
  // look for file
  std::string name = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+pID+".col";
  fstream test;
  test.open(name.c_str(),ios::in); 
  if (test) { 
    test.close();
    gen_colors=false;
  }
  Single_Amplitude* n = firstgraph;
  ngraph = 0;
  while (n) { 
    ++ngraph;
    n->Zprojecting(fl,ngraph,gen_colors);
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

  if (ngraph==0) {
    if (msg.LevelIsTracking()) {
      msg.Out()<<"No graph found for ";
      for (short int i=0;i<N;i++) msg.Out()<<fl[i]<<";";
      msg.Out()<<endl;
    }
    return;
  }


  //Colors
  CFCol_Matrix   = new CFColor(N,firstgraph,gen_colors,pID);



  for (int i=0;i<CFCol_Matrix->MatrixSize();i++) graphs.push_back(new Color_Group());

  //On-Switches
  //int* switch_graphs = new int[ngraph];
  // *FS*  for(short int i=0;i<ngraph;i++) onswitch[i] = 1;
  //Kicker(switch_graphs,ngraph,pID);

  n = firstgraph;

  // fill color groups
  int ncount = 0;
  int maxorder = 1;
  for(int i=0; i<N; i++) if (fl[i].IsKK()) maxorder--;
  if (maxorder<0) {
    msg.Error()<<"ERROR in Amplitude_Handler::CompleteAmplitudes :"<<std::endl
	       <<"   Multiple external KK-particles not supported. Abort the run."<<std::endl;
    abort();
  }
  while (n) {
    while(TOrder(n)>maxorder){
      ncount++;	   
      n=n->Next;
      if (!n) break;
    }
    if (!n) break; 
    pointlist.push_back(n->GetPointlist()); 
    graphs[CFCol_Matrix->CFMap(ncount)]->Add(n,CFCol_Matrix->CFSign(ncount));
    n = n->Next;
    ncount++;	   
  }
  
  //delete[] switch_graphs;
  ngraph=pointlist.size();
   
  for (size_t i=0;i<graphs.size();i++) graphs[i]->BuildGlobalString(b,N,BS,fl,shand);

  int dummy = 0;
  SetNumber(dummy);  
  namplitude = dummy;

  if (msg.LevelIsTracking()) {
    PrintGraph();
//     BS->PrintMomlist();
  }
  if (m_print_graph) {
    Amplitude_Output ao(pID,top);
    for (int i=0;i<namplitude;i++) {
      Amplitude_Base * am = GetAmplitude(i);
      if (am->Size()==1) {
        ao.WriteOut(am->GetPointlist());
      }
      else {
        ao.BeginSuperAmplitude();
	Amplitude_Group * ag=dynamic_cast<Amplitude_Group*>(am);
        for (int j=0;j<ag->Size();++j) ao.WriteOut((*ag)[j]->GetPointlist());
        ao.EndSuperAmplitude();
      }
    }
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
  if (CFCol_Matrix) delete CFCol_Matrix;
  if (probabs)      delete[] probabs;
  if (Mi)           delete[] Mi;
  if (ngraph>0) {
    Single_Amplitude * n; 
    while (firstgraph) {
      n = firstgraph->Next;
      delete firstgraph;
      firstgraph = n;
    }
  }
//   if (p_SCT!=NULL) delete p_SCT;
}

int Amplitude_Handler::PropProject(Amplitude_Base* f,int zarg)
{
  if (zarg<100) return zarg;

  Pfunc_List* pl = f->GetPlist();
  for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
    if ((*pit)->arg[0]==iabs(zarg)) return (*pit)->momnum; 
  }  
  msg.Error()<<"ERROR in Amplitude_Handler::PropProject() :"<<endl
	     <<"   Did not find a mom-number for propagator. Abort the run."<<std::endl;
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
  MyStrStream sstr;
  sstr<<i;
  string istr;
  sstr>>istr;
  return istr;
}

void Amplitude_Handler::OptimizeProps(int N,Single_Amplitude* f1)
{
  Zfunc_List* zlist = f1->GetZlist();
  Pfunc_List* pl = f1->GetPlist();
  for(Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit){
    if((*pit)->argnum > (N/2)+1){
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

    for(size_t i=0;i<propselect.size();i++){
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

	  Zfunc *zh0,*zh1;
	  
	  zh0=zh[0];zh1=zh[1];   //unique order of zfunctions
	  if(zh[0]->m_narg>zh[1]->m_narg){zh0=zh[1];zh1=zh[0];}
	  if(zh[0]->m_narg==zh[1]->m_narg){
	    for(int j=0;j<zh[0]->m_narg;j++){
	      if(PropProject(f1,zh[0]->p_arguments[j])<PropProject(f1,zh[1]->p_arguments[j]))break;
	      if(PropProject(f1,zh[0]->p_arguments[j])>PropProject(f1,zh[1]->p_arguments[j])){zh0=zh[1];zh1=zh[0];break;}
	    }
	  }
	  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();){
	    if((*zit)==zh0 || (*zit)==zh1){ zit=zlist->erase(zit);}
	    else zit++;
	  }
	  Zfunc_Group *sf=new Zfunc_Group(*zh0,*zh1,propselect[i],pl);
	  zlist->push_back(sf);
	  sf->Print();
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
  std::string name =rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+pID+".kick";
  //  sprintf(name,"%s/Kicker.dat",(rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pID).c_str());
  fstream test;
  test.open(name.c_str(),ios::in);

  if (test) {
    test.close();
    fstream from;
    from.open(name.c_str(),ios::in);
    //read in
    int i,sw;
    for(;from;) {
      from>>i>>sw;
      Switch_Vector[i-1] =sw;
      if (sw==0) msg_Tracking()<<"Amplitude_Handler::Kicker : Diagram "<<i<<" kicked!"<<endl;
      if (i==ngraph) break;
    }
    return;
  }

  test.close();

  for(short int i=0;i<ngraph;i++) Switch_Vector[i] = 1;
  

  ofstream to;
  to.open(name.c_str());

  for(short int i=0;i<ngraph;i++) to<<i+1<<"     "<<Switch_Vector[i]<<endl;
}

Point* Amplitude_Handler::GetPointlist(int n)
{ return pointlist[n];}

void Amplitude_Handler::Reset_ProbAbs()
{
  for (size_t i=0;i<graphs.size();i++) probabs[i] = 0.;
}

double Amplitude_Handler::Get_Probab(int i) {return probabs[i];}


Complex Amplitude_Handler::Zvalue(String_Handler * sh, int ihel)
{ // Called when no libraries are present (compiled)
  for (size_t i=0;i<graphs.size();i++){
    Mi[i] = graphs[i]->Zvalue(sh, ihel);
  }
  Complex M(0.,0.);
  for (size_t i=0;i<graphs.size();i++) {
    for (size_t j=0;j<graphs.size();j++) {
      M += Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
    }
  }
  return M;
}

Complex Amplitude_Handler::Zvalue(int ihel)
{ // Called for actual calculation of the CS
  for (size_t i=0;i<graphs.size();i++) {
    Mi[i] = graphs[i]->Zvalue(ihel);
  }

  Complex M(0.,0.);
  for (size_t i=0;i<graphs.size();i++) {
    for (size_t j=0;j<graphs.size();j++) {
      M+= Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
    }
  }
  return M;
}

double Amplitude_Handler::Zvalue(Helicity* hel)
{ 
  // 2D array for the amplitudes.
  typedef std::vector<Complex> CVec;
  std::vector<CVec> A;
  A.resize(graphs.size());

  /* For all graphs: Calculate all the helicity formalisms amplitudes and transform them to
     desired polarisation states, if nessecary. */
  for (size_t col=0; col<graphs.size(); ++col) {
    for (size_t ihel=0; ihel<hel->MaxHel(); ++ihel) A[col].push_back(graphs[col]->Zvalue(ihel));
    hel->SpinorTransformation(A[col]);
  }

  /* Calculate the scattering matrix M out of the amplitudes using the color matrix. Sum up
     the weighted Ms to obtain a pre-cross section sigma. */
  double sigma=0;
  for (size_t ihel=0; ihel<hel->MaxHel(); ++ihel) {
    if (hel->On(ihel)) {
      Complex M(0., 0.);
      for (size_t i=0;i<graphs.size();i++) {
	for (size_t j=0;j<graphs.size();j++) {
	  M+= A[i][ihel]*conj(A[j][ihel])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
	}
      }
      sigma += M.real() * hel->Multiplicity(ihel) * hel->PolarizationFactor(ihel);
    }
  }
  return sigma;
}


Complex Amplitude_Handler::Zvalue(int ihel,int* sign)
{ // This is called for the gauge test
  for (size_t i=0;i<graphs.size();i++) Mi[i] = graphs[i]->Zvalue(ihel,sign);

  Complex mcm,M(0.,0.);
  double max = 0.;
  for (size_t i=0;i<graphs.size();i++) {
    for (size_t j=0;j<graphs.size();j++) {
      mcm = Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
      M+=mcm;
      max = ATOOLS::Max(max,abs(mcm));
    }
  }
  if (abs(M)/max<(ATOOLS::Accu()*1.e-2)) return Complex(0.,0.); 
  return M;
}

void Amplitude_Handler::FillAmplitudes(Amplitude_Tensor *atensor,Helicity* hel,double sfactor)
{
  atensor->SetColorMatrix(CFCol_Matrix->GetCMatrix());
  vector<Complex> amps; amps.resize(graphs.size());
  
  for (size_t ihel=0; ihel<hel->MaxHel(); ++ihel) {
    for (size_t i=0;i<graphs.size();i++) amps[i]=graphs[i]->Zvalue(ihel)*sfactor;
    atensor->InsertAmplitude(amps,ihel);
  }
}

// ATOOLS::Spin_Correlation_Tensor* Amplitude_Handler::GetSpinCorrelations(Helicity* hel, 
// 									size_t nIn)
// { 
//   // If there is no pre-SCT constructed, then do this now
//   if (p_SCT==NULL) {
//     Spin_Correlation_Tensor SCTmethods;
  
//     std::vector<int> pList;
//     for (size_t i=nIn; i<hel->Nflavs(); ++i)
//       if ( SCTmethods.PossibleParticle( hel->GetFlav(i).Kfcode() ) )
// 	pList.push_back(i);

//     std::vector<int> AmplNrs;
//     for (size_t i=0; i<hel->MaxHel(); ++i) AmplNrs.push_back(i);

//     p_SCT = new AMEGIC_SCT(AmplNrs, AmplNrs, hel, &pList);
//   }

//   // Create an SCT from the pre-SCT
//   return p_SCT->CreateSCT(&graphs, CFCol_Matrix, hel);
// }

int Amplitude_Handler::TOrder(Single_Amplitude* a)
{  
  if(rpa.gen.Model()!=Model_Type::ADD) return 0;
  int cnt=0;
  Pfunc_List* pl = a->GetPlist();
  for(Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit)
    if((*pit)->fl.IsKK())cnt++;
  return cnt;
} 

int Amplitude_Handler::CompareAmplitudes(Amplitude_Handler* c_ampl, double & sf)
{
  if (GetTotalGraphNumber()!=c_ampl->GetTotalGraphNumber()) return 0;
  sf = 1.;

  Single_Amplitude * n = firstgraph;
  Single_Amplitude * n_cmp = c_ampl->GetFirstGraph();
  for (int i=0;i<GetTotalGraphNumber();i++) {
    double factor = 1.;
    if (!SingleCompare(n->GetPointlist(),n_cmp->GetPointlist(),factor)) return 0;
    if (i==0) sf = factor;
    else if(!ATOOLS::IsEqual(sf,factor)) return 0;
    n     = n->Next;
    n_cmp = n_cmp->Next;
  }
  return 1;
}

int Amplitude_Handler::SingleCompare(Point* p1,Point* p2, double & sf)
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
  if (p1->fl.Mass()!=p2->fl.Mass()) return 0;
  if (p1->fl.Spin()!=p2->fl.Spin()) return 0;

  //outgoing number equal
  if ((p1->left==0) && (p2->left==0)) {
    if (p1->number!=p2->number) return 0;
                           else return 1;
  }

  if ((p1->left==0) && (p2->left!=0)) return 0;
  if ((p1->left!=0) && (p2->left==0)) return 0;

  //Check extended Color_Functions
  if (p1->Color->Type()!=p2->Color->Type()) return 0;
  
  //Couplings equal
  //if (p1->ncpl!=p2->ncpl) return 0;
  Complex ratio = Complex(0.,0.);
  for (int i=0;i<2;i++) {
    if (ratio==Complex(0.,0.) && p2->cpl[i]!=Complex(0.,0.)) ratio = p1->cpl[i]/p2->cpl[i];
    if (!ATOOLS::IsEqual(p2->cpl[i]*ratio,p1->cpl[i])) return 0;
  } 
  sf *= abs(ratio);
  // return 1 if equal and 0 if different
  
  if (SingleCompare(p1->middle,p2->middle,sf)) {
    int sw1 = SingleCompare(p1->left,p2->left,sf);
    if (sw1) sw1 = SingleCompare(p1->right,p2->right,sf);
    return sw1;
  }
  return 0;
}

void Amplitude_Handler::FillPointlist()
{
  Single_Amplitude* n = firstgraph;  
  while (n) {
   pointlist.push_back(n->GetPointlist()); 
   n = n->Next;
  }
}







