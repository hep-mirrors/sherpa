#include "AMEGIC++/Amplitude/Amplitude_Handler.H"
#include "AMEGIC++/Amplitude/Amplitude_Output.H"
#include "AMEGIC++/Amplitude/Zfunctions/Zfunc_Calc.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "AMEGIC++/Amplitude/Amplitude_Generator.H"
#include "AMEGIC++/Amplitude/Amplitude_Manipulator.H"
#include "AMEGIC++/Amplitude/Color_Group.H"
#include <iostream>
#include <stdio.h>
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AMEGIC++/Main/Process_Tags.H"
#include "ATOOLS/Org/IO_Handler.H"
#include "METOOLS/Main/Spin_Structure.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

Amplitude_Handler::Amplitude_Handler(int N,Flavour* fl,int* b,Process_Tags* pinfo,
				     Amegic_Model * model,Topology* top,
				     std::vector<int> & _maxcpl,
				     std::vector<int> & _mincpl,int & _ntchan,
				     MODEL::Coupling_Map *const cpls,
				     Basic_Sfuncs* BS,String_Handler* _shand, 
				     std::string print_graph,bool create_4V,
				     bool cutvecprop)
  : m_cutvecprop(cutvecprop), shand(_shand), CFCol_Matrix(0), Mi(0),
    m_print_graph(print_graph)
{
  m_maxcpl=_maxcpl;
  m_mincpl=_mincpl;
  for (size_t i(0);i<m_mincpl.size();++i) m_mincpl[i]*=2;
  for (size_t i(0);i<m_maxcpl.size();++i) m_maxcpl[i]*=2;
  groupname = "Amplitude_Handler";
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
  }
  else sfl=fl;

  //core process
  gen = new Amplitude_Generator(nin+pinfo->Nout(),sfl,b,model,top,m_maxcpl,_ntchan,BS,shand,create_4V);
  subgraphlist[0] = gen->Matching();
  delete gen;

  //decay processes
  for (int i=1;i<=ndecays;i++) {
    int j=i;
    Process_Tags *pi=pinfo->GetDecay(j);
//     pi->Print();cout<<endl;
    sfl[0] = *(pi->p_fl);
    pi->GetFlavList(sfl+1);
    gen = new Amplitude_Generator(1+pi->Nout(),sfl,b_dec,model,top,std::vector<int>(2,99),-99,BS,shand);
    subgraphlist[i] = gen->Matching();
    if (subgraphlist[i]==NULL) {
      ndecays = 0;
      subgraphlist[0] = NULL;
    }
    delete gen;
  }

  if (ndecays>0) {
    sfl[0] = fl[0];
    sfl[1] = fl[1];
    pinfo->GetFlavList(sfl+nin);
  }

  //core process
  gen = new Amplitude_Generator(nin+pinfo->Nout(),sfl,b,model,top,m_maxcpl,_ntchan,BS,shand,create_4V);
  subgraphlist[0] = gen->Matching();
  delete gen;

  if (msg_LevelIsTracking()) {
    msg_Out()<<"Amplitude_Handler::Amplitude_Handler:"<<endl;
    int f=1;
    for(int i=0;i<ndecays+1;i++) {
      int j=0;
      Single_Amplitude* nn = subgraphlist[i];
      while (nn){ 
	++j;
	nn = nn->Next;
      }
      msg_Out()<<"Process "<<i;
      if (i==0)msg_Out()<<" (core)";
      else msg_Out()<<" (decay)";
      msg_Out()<<" has "<<j<<" Amplitudes"<<endl;
      f*=j;
    }
    msg_Out()<<"Total: "<<f<<" Amplitudes"<<endl;
  }

  if (ndecays==0 || subgraphlist[0]==0) firstgraph = subgraphlist[0];
  else ConstructSignalAmplitudes(N,fl,b,pinfo,subgraphlist,BS);
  
  Amplitude_Manipulator(N,fl,b,ndecays).FixSign(firstgraph);

  Single_Amplitude* n = firstgraph;
  Single_Amplitude* prev = firstgraph;
  ntotal = 0;
  while (n){ 
    if (TOrder(n)>1) {
      Single_Amplitude* next = n->Next;
      if (n==firstgraph) firstgraph = next;
      else prev->Next = next;
      delete n;
      n = next;
    }
    else {
      ++ntotal;
      prev = n;
      n->GetPointlist()->GeneratePropID();
      n->SetOrder();
      n = n->Next;
    }
  }
  msg_Tracking()<<"Total number of Amplitudes "<<ntotal<<endl;
  ngraph = ntotal;

  if (ngraph!=0) {
    p_aqcd=cpls->Get("Alpha_QCD");
    p_aqed=cpls->Get("Alpha_QED");
  }
  
  delete [] subgraphlist;
  delete [] b_dec;
}

void Amplitude_Handler::ConstructSignalAmplitudes(int N,Flavour* fl,int* b,
						  Process_Tags* pinfo,Single_Amplitude** sglist,
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
  
  delete [] nl;
  for (int i=0;i<ndecays+1;i++) {
    n = sglist[i];
    while (n) {
      next = n->Next;
      delete n;
      n = next;
    }
  }
}

void Amplitude_Handler::CompleteAmplitudes(int N,Flavour* fl,int* b,Polarisation* pol,
					   Topology* top,Basic_Sfuncs* BS,std::string pID,
					   char emit,char spect)
{
  Single_Amplitude* n = firstgraph;
  ngraph = 0;
  while (n) { 
    ++ngraph;
    n->Zprojecting(fl,ngraph,true,m_cutvecprop);
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
    if (msg_LevelIsTracking()) {
      msg_Out()<<"No graph found for ";
      for (short int i=0;i<N;i++) msg_Out()<<fl[i]<<";";
      msg_Out()<<endl;
    }
    return;
  }


  //Colors
  if (emit!=spect && emit!=127) {
    char cemit=emit,cspect=spect;
    if (fl[(int)emit].IsGluon() || IsGluino(fl[(int)emit])) cemit+='A';
    else cemit+='i';
    if (fl[(int)spect].IsGluon()|| IsGluino(fl[(int)spect])) cspect+='A';
    else cspect+='i';
    CFCol_Matrix   = new CFColor(N,firstgraph,fl,cemit,cspect,pID);
  }
  else {
    CFCol_Matrix   = new CFColor(N,firstgraph,fl,emit,spect,pID);
    if (emit==127) {
      for (int i=0;i<N-1;i++) if (fl[i].Strong()) {
	for (int j=i+1;j<N;j++) if (fl[j].Strong()) {
	  char cemit=i,cspect=j;
	  if (fl[i].IsGluon() || IsGluino(fl[i])) cemit+='A';
	  else cemit+='i';
	  if (fl[j].IsGluon() || IsGluino(fl[j])) cspect+='A';
	  else cspect+='i';
	  string sij=pID+string("_S")+ToString(i)+string("_")+ToString(j);
	  //msg_Out()<<METHOD<<" new CFColor("<<sij<<")."<<std::endl;
	  CFColor* mcfc = new CFColor(N,firstgraph,fl,cemit,cspect,sij);
	  CFCol_MMatrixMap[i*100+j] = mcfc;
	}
      }
    }
  }
  


  for (int i=0;i<CFCol_Matrix->MatrixSize();i++) graphs.push_back(new Color_Group());

  n = firstgraph;

  // fill color groups
  int ncount = 0;
  while (n) {
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

  if (msg_LevelIsTracking()) {
    PrintGraph();
    //BS->PrintMomlist();
  }
  if (m_print_graph!="") {
    Amplitude_Output ao(pID,top,m_print_graph);
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

//   probabs = new double[graphs.size()];
  Mi      = new Complex[graphs.size()];
}

std::vector<int> Amplitude_Handler::MaxCpl() const
{
  std::vector<int> maxcpl;
  for (size_t i=0;i<graphs.size();i++)
    for (size_t j=0;j<graphs.size();j++) {
      if (maxcpl.size()<m_cplmatrix[i][j].size())
	maxcpl.resize(m_cplmatrix[i][j].size(),0);
      for (size_t k=0;k<m_cplmatrix[i][j].size();++k)
	maxcpl[k]=Max(maxcpl[k],m_cplmatrix[i][j][k]);
    }
  return maxcpl;
}

std::vector<int> Amplitude_Handler::MinCpl() const
{
  std::vector<int> mincpl;
  for (size_t i=0;i<graphs.size();i++)
    for (size_t j=0;j<graphs.size();j++) {
      if (mincpl.size()<m_cplmatrix[i][j].size())
	mincpl.resize(m_cplmatrix[i][j].size(),99);
      for (size_t k=0;k<m_cplmatrix[i][j].size();++k)
	mincpl[k]=Min(mincpl[k],m_cplmatrix[i][j][k]);
    }
  return mincpl;
}

void Amplitude_Handler::StoreAmplitudeConfiguration(std::string path)
{
  std::string name = path+"/Cluster.dat";
  IO_Handler ioh;
  ioh.SetFileName(name);
  ioh.Output("",int(graphs.size()));
  My_Out_File cplfile(path+"/Couplings.dat");
  My_Out_File sqrcplfile(path+"/SquaredCouplings.dat");
  cplfile.Open();
  sqrcplfile.Open();
  m_on.resize(graphs.size());
  m_aon.resize(graphs.size(),0);
  m_cplmatrix.resize(graphs.size());
  for (size_t i=0;i<graphs.size();i++) {
    m_on[i].resize(graphs.size(),1);
    m_cplmatrix[i].resize(graphs.size(),graphs[i]->GetOrder());
    for (size_t j=0;j<graphs.size();j++) {
      if (m_cplmatrix[i][j].size()<graphs[j]->GetOrder().size())
	m_cplmatrix[i][j].resize(graphs[j]->GetOrder().size(),0);
      for (size_t k=0;k<graphs[j]->GetOrder().size();++k)
	m_cplmatrix[i][j][k]+=graphs[j]->GetOrder()[k];
      *sqrcplfile<<i<<" "<<j<<" "<<m_cplmatrix[i][j]<<"\n";
      for (size_t k=0;k<Min(m_cplmatrix[i][j].size(),m_maxcpl.size());++k)
	if (m_cplmatrix[i][j][k]>m_maxcpl[k]) m_on[i][j]=0;
      for (size_t k=0;k<Min(m_cplmatrix[i][j].size(),m_mincpl.size());++k)
	if (m_cplmatrix[i][j][k]<m_mincpl[k]) m_on[i][j]=0;
      if (m_on[i][j]) m_aon[i]=1;
    }
    *cplfile<<i<<" "<<graphs[i]->GetOrder()<<"\n";
    int size=graphs[i]->Size();
    int *nums= new int[size];
    for (int j=0;j<size;j++) nums[j]=(*graphs[i])[j]->GetNumber();
    ioh.ArrayOutput<int>("",nums,size);
    delete [] nums;
  }
}

void Amplitude_Handler::RestoreAmplitudes(std::string path)
{
  std::string name = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+path+"/Cluster.dat";
  My_In_File cplfile(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+path+"/Couplings.dat");
  My_In_File sqrcplfile(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+path+"/SquaredCouplings.dat");
  if (!cplfile.Open()) THROW(fatal_error,"Missing coupling data");
  if (!sqrcplfile.Open()) THROW(fatal_error,"Missing squared coupling data");
  IO_Handler ioh;
  ioh.SetFileNameRO(name);
  size_t cg = ioh.Input<int>("");
  if (cg!=graphs.size()) {
    msg_Error()<<"ERROR in Amplitude_Handler::RestoreAmplitudes() :"<<endl
	       <<"   Stored Cluster and Color information incompatible! Abort the run."<<std::endl;
    abort();
  }
  int cnt=0;
  Amplitude_Base* ab;
  m_on.resize(graphs.size());
  m_aon.resize(graphs.size(),0);
  m_cplmatrix.resize(graphs.size());
  for (size_t i=0;i<graphs.size();i++) {
    int *nums, ci, cj;
    std::string ords;
    m_on[i].resize(graphs.size(),1);
    m_cplmatrix[i].resize(graphs.size());
    for (size_t j=0;j<graphs.size();j++) {
      *sqrcplfile>>ci>>cj>>ords;
      if (ci!=i || cj!=j) THROW(fatal_error,"Invalid coupling data");
      Data_Reader read(",",";",")","(");
      read.SetString(ords);
      read.VectorFromString(m_cplmatrix[i][j],"");
      for (size_t k=0;k<Min(m_cplmatrix[i][j].size(),m_maxcpl.size());++k)
	if (m_cplmatrix[i][j][k]>m_maxcpl[k]) m_on[i][j]=0;
      for (size_t k=0;k<Min(m_cplmatrix[i][j].size(),m_mincpl.size());++k)
	if (m_cplmatrix[i][j][k]<m_mincpl[k]) m_on[i][j]=0;
      if (m_on[i][j]) m_aon[i]=1;
    }
    *cplfile>>ci>>ords;
    Data_Reader read(",",";",")","(");
    read.SetString(ords);
    std::vector<int> ord;
    read.VectorFromString(ord,"");
    if (ci!=(int)i) THROW(fatal_error,"Invalid coupling data");
    nums=ioh.ArrayInput<int>("");
    int size=ioh.Nx();
    for (int j=0;j<size;j++) {
      ab=new Single_Amplitude_Base(shand,nums[j]);
      ab->DefineOrder(ord);
      graphs[i]->Add(ab);
      
      m_ramplist.push_back(ab);
    }
    cnt+=size;
    delete [] nums;
  }
  namplitude = cnt;
}

void Amplitude_Handler::CompleteLibAmplitudes(int N,std::string pID,std::string lib,
					      char emit,char spect,Flavour* fl)
{
  std::string name = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+pID+".map";
  My_In_File from(name);
  from.Open();
  shand->Get_Generator()->ReadCouplings(*from);
  from.Close();

  Single_Amplitude* n = firstgraph;
  ngraph = 0;
  while (n) { 
    ++ngraph;
    n = n->Next;
  }

  //Colors
  //Colors
  if (emit!=spect && emit!=127) {
    char cemit=emit,cspect=spect;
    if (fl[(int)emit].IsGluon() || IsGluino(fl[(int)emit])) cemit+='A';
    else cemit+='i';
    if (fl[(int)spect].IsGluon() || IsGluino(fl[(int)spect])) cspect+='A';
    else cspect+='i';
    CFCol_Matrix   = new CFColor(N,firstgraph,fl,cemit,cspect,pID,true);
  }
  else {
    CFCol_Matrix   = new CFColor(N,firstgraph,fl,emit,spect,pID,true);
    if (emit==127) {
      for (int i=0;i<N-1;i++) if (fl[i].Strong()) {
	for (int j=i+1;j<N;j++) if (fl[j].Strong()) {
	  char cemit=i,cspect=j;
	  if (fl[i].IsGluon() || IsGluino(fl[i])) cemit+='A';
	  else cemit+='i';
	  if (fl[j].IsGluon() || IsGluino(fl[j])) cspect+='A';
	  else cspect+='i';
	  string sij=pID+string("_S")+ToString(i)+string("_")+ToString(j);
	  
	  CFColor* mcfc = new CFColor(N,firstgraph,fl,cemit,cspect,sij,true);
	  CFCol_MMatrixMap[i*100+j] = mcfc;
	}
      }
    }
  }
  for (int i=0;i<CFCol_Matrix->MatrixSize();i++) {
    //msg_Out()<<METHOD<<" push_back new Colour_Group["<<i<<"]."<<std::endl;
    graphs.push_back(new Color_Group());
  }
  n = firstgraph;

  // fill color groups
  int ncount = 0;

  while (n) {
    pointlist.push_back(n->GetPointlist()); 
    n = n->Next;
    ncount++;	   
  }
 
  RestoreAmplitudes(lib);
  
  ngraph=pointlist.size();

  Mi      = new Complex[graphs.size()];
}

Amplitude_Handler::~Amplitude_Handler() 
{
  for (size_t i=0;i<graphs.size();i++) delete graphs[i];
  graphs.clear();

  for (size_t i=0;i<m_ramplist.size();i++) delete m_ramplist[i];
  m_ramplist.clear();

  if (CFCol_Matrix) delete CFCol_Matrix;
  if (Mi)           delete[] Mi;
  if (ngraph>0) {
    Single_Amplitude * n; 
    while (firstgraph) {
      n = firstgraph->Next;
      delete firstgraph;
      firstgraph = n;
    }
  }

  for(CFC_iterator it=CFCol_MMatrixMap.begin();it!=CFCol_MMatrixMap.end();++it)
    delete it->second;
}

int Amplitude_Handler::PropProject(Amplitude_Base* f,int zarg)
{
  if (zarg<100) return zarg;

  Pfunc_List* pl = f->GetPlist();
  for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
    if ((*pit)->arg[0]==iabs(zarg)) return (*pit)->momnum; 
  }  
  msg_Error()<<"ERROR in Amplitude_Handler::PropProject() :"<<endl
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
      if((*pit)->fl.Kfcode()==kf_photon||(*pit)->fl.Kfcode()==kf_Z) 
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
	    (zh[0]->m_type=="Y" || zh[0]->m_type=="Z") ) {

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

Point* Amplitude_Handler::GetPointlist(int n)
{ return pointlist[n];}


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

Complex Amplitude_Handler::CommonColorFactor()
{
  if (graphs.empty()) return Complex(0.0,0.0);
  Complex C(CFCol_Matrix->Mij(0,0));
  for (size_t i=0;i<graphs.size();i++)
    for (size_t j=0;j<graphs.size();j++)
      if (C!=CFCol_Matrix->Mij(i,j)) return Complex(0.,0.);
  return C;
}

Complex Amplitude_Handler::Zvalue(int ihel)
{ 
  // Called for actual calculation of the CS
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  for (size_t i=0;i<graphs.size();i++) {
    if (m_aon.size() && !m_aon[i]) continue;
    double cplfac(1.0);
    const std::vector<int> &order(graphs[i]->GetOrder());
    if (p_aqcd && order.size()>0 && order[0]) {
#ifdef DEBUG__BG
      msg_Debugging()<<"  qcd: "<<sqrt(p_aqcd->Factor())<<" ^ "<<order[0]
		     <<" = "<<pow(p_aqcd->Factor(),order[0]/2.0)<<"\n";
#endif     
      cplfac *= pow(p_aqcd->Factor(),order[0]/2.0);
    }  
    if (p_aqed && order.size()>1 && order[1]) {
#ifdef DEBUG__BG
      msg_Debugging()<<"  qed: "<<sqrt(p_aqed->Factor())<<" ^ "<<order[1]
		     <<" = "<<pow(p_aqed->Factor(),order[1]/2.0)<<"\n";
#endif   
      cplfac *= pow(p_aqed->Factor(),order[1]/2.0); 
    }
#ifdef DEBUG__BG
    msg_Debugging()<<"  graph "<<i<<" -> "<<cplfac<<"\n";
#endif  
    Mi[i] = cplfac*(graphs[i]->Zvalue(ihel));
  }
  Complex M(0.,0.);
  for (size_t i=0;i<graphs.size();i++) {
    for (size_t j=0;j<graphs.size();j++) {
      if (m_on.empty() || m_on[i][j]) 
      M+= Mi[i]*conj(Mi[j])*CFCol_Matrix->Mij(i,j);  //colfactors[i][j];
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif  
  return M;
}

Complex Amplitude_Handler::Zvalue(int ihel,int ci,int cj)
{// Called for actual calculation of the CS
  int cid = 100*ci+cj;
  if (cj<ci) cid = 100*cj+ci;
  CFColor *col = CFCol_Matrix; 
  if (cid!=0) {
    CFC_iterator cit = CFCol_MMatrixMap.find(cid);
    if (cit==CFCol_MMatrixMap.end()) {
      msg_Error()<<"ERROR in Amplitude_Handler::Zvalue :"<<std::endl
		 <<"   Color matrix ("<<ci<<"/"<<cj<<") not found! Abort the run."<<std::endl;
      abort();
    }
    col = cit->second;
  }
  
  for (size_t i=0;i<graphs.size();i++) {
    double cplfac(1.0);
    const std::vector<int> &order(graphs[i]->GetOrder());
    if (p_aqcd && order.size()>0 && order[0]) {
      cplfac *= pow(p_aqcd->Factor(),order[0]/2.0);
    }
    if (p_aqed && order.size()>1 && order[1]) {
      cplfac *= pow(p_aqed->Factor(),order[1]/2.0);
    }
    Mi[i] = cplfac*(graphs[i]->Zvalue(ihel));
  }
  
  Complex M(0.,0.);
  for (size_t i=0;i<graphs.size();i++) {
    for (size_t j=0;j<graphs.size();j++) {
      if (m_on[i][j]) 
      M+= Mi[i]*conj(Mi[j])*col->Mij(i,j);  //colfactors[i][j];
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
    double cplfac(1.0);
    const std::vector<int> &order(graphs[col]->GetOrder());
    if (p_aqcd && order.size()>0 && order[0]) {
      cplfac *= pow(p_aqcd->Factor(),order[0]/2.0);
    }
    if (p_aqed && order.size()>1 && order[1]) {
      cplfac *= pow(p_aqed->Factor(),order[1]/2.0);
    }
    for (size_t ihel=0; ihel<hel->MaxHel(); ++ihel) A[col].push_back(cplfac*(graphs[col]->Zvalue(ihel)));
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
	  if (m_on[i][j]) 
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
  for (size_t i=0;i<graphs.size();i++) {
    double cplfac(1.0);
    const std::vector<int> &order(graphs[i]->GetOrder());
    if (p_aqcd && order.size()>0 && order[0]) {
      cplfac *= pow(p_aqcd->Factor(),order[0]/2.0);
    }
    if (p_aqed && order.size()>1 && order[1]) {
      cplfac *= pow(p_aqed->Factor(),order[1]/2.0);
    }
    Mi[i] = cplfac*(graphs[i]->Zvalue(ihel,sign));
  }
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

void Amplitude_Handler::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                       std::vector<std::vector<Complex> >& cols,
                                       Helicity* hel, double sfactor)
{
  cols.resize(graphs.size(),std::vector<Complex>(graphs.size()));
  for (size_t i=0; i<graphs.size(); ++i) {
    for (size_t j=0; j<graphs.size(); ++j) {
      cols[i][j]=CFCol_Matrix->Mij(i,j);
    }
  }
  ////////////////////////////////////////////////////// BEGIN HACK
  if (m_hm.empty()) {
    m_hm.resize(hel->MaxHel());
    METOOLS::Spin_Structure<int> amp(hel->GetFlavs(),0);
    for (size_t ihel=0; ihel<hel->MaxHel(); ++ihel) {
      std::vector<int> ch(amp(ihel));
      for (size_t j(0);j<ch.size();++j){
	if (hel->GetFlavs()[j].IsScalar()) continue;
	if (ch[j]<2) ch[j]=1-ch[j];
      }
      m_hm[ihel]=amp(ch);
    }
  }
  ////////////////////////////////////////////////////// END HACK
  for (size_t i=0;i<graphs.size();i++) {
    amps.push_back(METOOLS::Spin_Amplitudes(hel->GetFlavs(),Complex(0.0,0.0)));
    for (size_t ihel=0; ihel<hel->MaxHel(); ++ihel) {
      amps.back().Insert(graphs[i]->Zvalue(ihel)*sfactor, m_hm[ihel]);
    }
  }
  
}

int Amplitude_Handler::TOrder(Single_Amplitude* a)
{  
  if(MODEL::s_model->Name()!="ADD") return 0;
  return a->GetPointlist()->CountKK();
} 

int Amplitude_Handler::CompareAmplitudes(Amplitude_Handler* c_ampl, double & sf, map<string,Complex> & cplmap)
{
  m_flavourmap.clear();
  if (GetTotalGraphNumber()!=c_ampl->GetTotalGraphNumber()) return 0;
  sf = 1.;

  Single_Amplitude * n = firstgraph;
  Single_Amplitude * n_cmp = c_ampl->GetFirstGraph();
  for (int i=0;i<GetTotalGraphNumber();i++) {
    double factor = 1.;
    if (!SingleCompare(n->GetPointlist(),n_cmp->GetPointlist(),factor,cplmap)) {
      m_flavourmap.clear();
      return 0;
    }
    if (i==0) sf = factor;
    else if(sf!=factor) {
      m_flavourmap.clear();
      return 0;
    }
    n     = n->Next;
    n_cmp = n_cmp->Next;
  }
  return 1;
}

int Amplitude_Handler::SingleCompare(Point* p1,Point* p2, double & sf, map<string,Complex> & cplmap)
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
  if (p1->fl.Width()!=p2->fl.Width()) return 0;
  if (p1->fl.Spin()!=p2->fl.Spin()) return 0;

  //outgoing number equal
  if ((p1->left==0) && (p2->left==0)) {
    if (p1->number!=p2->number) return 0;
    else {
      string pid=p2->GetPropID();
      if (m_flavourmap.find(pid)==m_flavourmap.end()) 
	m_flavourmap[pid]=p1->fl;
      return 1;
    }
  }

  if ((p1->left==0) && (p2->left!=0)) return 0;
  if ((p1->left!=0) && (p2->left==0)) return 0;

  //Check extended Color_Functions
  if (p1->Color->Type()!=p2->Color->Type()) return 0;
  
  //Couplings equal
  if (p1->v->cpl.size()!=p2->v->cpl.size()) return 0;
  for (int i=0;i<p1->v->cpl.size();i++) {
    if (p1->v->Coupling(i)!=p2->v->Coupling(i)) return 0;
    if (p1->cpl[i]!=p2->cpl[i]) return 0;
  }
  // return 1 if equal and 0 if different

  {
    string pid=p2->GetPropID();
    if (m_flavourmap.find(pid)==m_flavourmap.end()) {
      m_flavourmap[pid]=p1->fl;
    }
    else if (m_flavourmap[pid]!=p1->fl){
      return 0;  
    }
  }
  
  if (SingleCompare(p1->middle,p2->middle,sf,cplmap)) {
    int sw1 = SingleCompare(p1->left,p2->left,sf,cplmap);
    if (sw1) sw1 = SingleCompare(p1->right,p2->right,sf,cplmap);
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

bool Amplitude_Handler::CheckEFMap()
{
  bool mf(1);
  Single_Amplitude* n = firstgraph;  
  while (n) {
    mf=CheckSingleEFM(n->GetPointlist());     
    if (!mf) return 0;
    n = n->Next;
  }
  return mf;
}

bool Amplitude_Handler::CheckSingleEFM(Point* p)
{
  if (p->left==0) return 1;
  if (p->fl.IsBoson() && p->fl.IntCharge()!=0) return 0;
  bool mf=CheckSingleEFM(p->left);
  if (mf) mf=CheckSingleEFM(p->right);
  if (p->middle && mf) mf=CheckSingleEFM(p->middle);
  return mf;
}




