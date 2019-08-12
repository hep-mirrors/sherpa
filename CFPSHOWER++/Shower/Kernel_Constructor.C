#include "CFPSHOWER++/Shower/Kernel_Constructor.H"
#include "CFPSHOWER++/Shower/Shower.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Shower/Cluster_Definitions.H"
#include "CFPSHOWER++/Tools/CFP_Parameters.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Kernel_Constructor::Kernel_Constructor(Shower * shower): p_shower(shower) {}

Kernel_Constructor::~Kernel_Constructor() { }

bool Kernel_Constructor::Init(MODEL::Model_Base * const model,
			      PDF::ISR_Handler * const isr,const size_t & order)
{
  // Shower parameters and switches are pulled from the map in the CFP_Parameter
  // - parton shower cutoffs for FS and IS showering
  // - order of the splitting function and details of terms included
  // - scale setting schemes
  // - eventually recoil schemes
  m_asfactor[0] = (*cfp_pars)("k_alpha(FS)");
  m_asfactor[1] = (*cfp_pars)("k_alpha(IS)");
  m_muR2factor  = (*cfp_pars)("k_muR");
  m_muF2factor  = (*cfp_pars)("k_muF");
  m_kinscheme   = (*cfp_pars)["kinematics"];
  m_kfactor     = (*cfp_pars)["kfactor"];
  m_cplscheme   = (*cfp_pars)["couplings"];
  m_MEcorrs     = (*cfp_pars)["ME_corrections"];
  p_as    = (MODEL::Running_AlphaS*)(model->GetScalarFunction("alpha_S"));
  for (int i=0;i<2;++i) p_pdf[i]=isr->PDF(i);
  MakePermutations(order);
  InitializeKernels(model,order);
}

bool Kernel_Constructor::InitializeKernels(MODEL::Model_Base * const model,
					   const size_t & order) {
  msg_Out()<<"***************************************************************\n"
   	   <<METHOD<<" starts collecting splitting kernels.\n"
   	   <<"** available kernels:\n";
  SF_Getter::PrintGetterInfo(msg->Out(),25);
  GP_Getter::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n";
  // Going through vertex tables and translating 3-particle vertices into splitting
  // kernels.  The kernels consist of a splitting function and a gauge part, where
  // the latter handles the colour configuration and all aspects realted to the
  // coupling constants at higher orders, while the former encapsulates AP splitting
  // functions, triple-collinear parts, etc..
  // The set of flavour triplets makes sure we do not initialize splitting kernels for
  // any parton configuration more than once.

  for (size_t i=0;i<5;i++) p_kernels[i] = new map<Flavour, Kernels *>;
 
  MakeKernelsFromVertices(model);
  if (order>2) {
    msg_Out()<<"*********************************************************\n"
	     <<"Now start collating 1->2 kernels into 1->3 kernels.\n"
	     <<"*********************************************************\n";
    for (size_t type=1;type<2;type++) MakeKernelsFromKernels(kernel_type::code(type));
  }
  PrintKernels();
  return true;
}

void Kernel_Constructor::MakeKernelsFromVertices(MODEL::Model_Base * const model) {
  set<vector<Flavour> > vetoed;
  const Vertex_Table *vtab(model->VertexTable());
  for (Vertex_Table::const_iterator vlit=vtab->begin();
       vlit!=vtab->end();++vlit) {
    if (!vlit->first.IsQuark() && !vlit->first.IsGluon()) continue;
    for (Vertex_List::const_iterator vit=vlit->second.begin();
	 vit!=vlit->second.end();++vit) {
      MODEL::Single_Vertex * vertex = (*vit);
      // At LO only 3-particle vertices are allowed and we keep track of such
      // configurations only in the kernel_flavours set.  
      if (vertex->NLegs()>3) continue;
      if (vetoed.find(vertex->in)!=vetoed.end()) return;
      vetoed.insert(vertex->in);
      for (size_t type=1;type<2;type++) MakeKernels(vertex->in,kernel_type::code(type));
    }
  }
}

void Kernel_Constructor::MakeKernels(Flavour_Vector & flavs,const kernel_type::code & type) {
  // The kernels are initialised with the information stored in the Kernel_Info
  // struct, which carries information about:
  // - the flavours, 
  // - their configuration for splitter/spectator in the kernel_type struct (FF, FI, IF, and II)
  // - pointers to alphaS and alpha(QED)
  // - the tagging sequence
  // - to be implemented: order information, masses of particles, etc..
  // Initialised kernels are organised in 4 maps (one for each of the splitter/spectator
  // configurations), connecting splitting flavours with all possible kernels.
  Flavour split = flavs[0];
  Flavour_Vector newflavs;
  for (size_t i=1;i<flavs.size();i++) newflavs.push_back(flavs[i]); 
  size_t order  = newflavs.size()-2;
  //if (newflavs.size()>2) {
  msg_Out()<<METHOD<<"["<<type<<", order = "<<(order+1)<<"] for "<<split.Bar()<<" -> ";
  for (size_t i=0;i<newflavs.size();i++) msg_Out()<<newflavs[i]<<" ";
  msg_Out()<<" with "<<m_permutations[order].size()<<" permutations.\n";
  //}
  for (list<vector<size_t> >::iterator lit=m_permutations[order].begin();
       lit!=m_permutations[order].end();lit++) {
    Kernel_Info info(split,newflavs,kernel_type::code(type),(*lit));
    info.SetAlphaS(p_as);
    info.SetKFactor(m_kfactor);
    info.SetCplScheme(m_cplscheme);
    info.SetAsFactor(((type==kernel_type::FF || type==kernel_type::FI) ?
		      m_asfactor[0] : m_asfactor[1]));
    info.SetMuR2Factor(m_muR2factor);
    msg_Out()<<"   * looking for "<<info;
    Kernel * kernel = Kernel_Getter::GetObject("Kernel",info);
    //if (newflavs.size()>2)
    if (kernel!=0) {
      if (p_kernels[int(info.Type())]->find(info.GetSplit())==
	  p_kernels[int(info.Type())]->end()) {
	msg_Out()<<"Init new kernel list for type = "<<info.Type()
		 <<" & flav = "<<info.GetSplit()<<": "<<info;
	(*p_kernels[int(info.Type())])[info.GetSplit()] = new Kernels();
      }
      (*p_kernels[int(info.Type())])[info.GetSplit()]->push_back(kernel);
      msg_Out()<<"Add kernel to ["<<info.Type()<<"]["<<info.GetSplit()<<"]: "<<info;
      for (size_t beam = 0;beam<2;beam++) kernel->SetPDF(beam,p_pdf[beam]);
      kernel->SetPDFMinValue((*cfp_pars)("PDF_min"));
      kernel->SetPDFXMin((*cfp_pars)("PDF_min_X"));
      kernel->SetMSel(p_shower->GetMassSelector());
    }
  }
}


void Kernel_Constructor::MakeKernelsFromKernels(const kernel_type::code & type) {
  msg_Out()<<METHOD<<"["<<type<<"]:\n";
  if (type==kernel_type::IF || type==kernel_type::II) {
    msg_Out()<<METHOD<<" for initial state kernels.  Will have to think about this.\n";
    exit(1);
  }
  // Iterate over all 1->2 kernels and see if you can attach another 1->2 kernel
  // to the outgoing particles/flavours - based on all kernels of a given type and
  // kernels for FF splittings
  map<Flavour,Kernels *> * allkernels = p_kernels[int(type)];
  map<Flavour,Kernels *> * ffkernels  = p_kernels[int(kernel_type::FF)];
  // initiate a look-up table of "vetoed" flavour combination, to make sure we have every
  // combination of four flavours only once
  set<multiset<Flavour> > vetoed;
  for (map<Flavour,Kernels *>::iterator allkit=allkernels->begin();
       allkit!=allkernels->end();allkit++) {
    const Flavour & splitter    = allkit->first;
    const Kernels * flavkernels = allkit->second;
    for (Kernel_Vector::const_iterator kit=flavkernels->begin();
	 kit!=flavkernels->end();kit++) {
      if ((*kit)->Tags(0)==1) continue;
      Flavour_Vector outs = (*kit)->GetFlavs();
      // go over the two outgoing flavours and produce trial final states by keeping one
      // and replacing the other with two outgoing flavours from the already initialised
      // 1->2 ffkernels.
      for (size_t j=0;j<2;j++) {
	if (ffkernels->find(outs[j])==ffkernels->end()) continue;
	Kernel_Vector * replkernels = (*ffkernels)[outs[j]];
	for (size_t k=0;k<replkernels->size();k++) {
	  Kernel * rkernel = (*replkernels)[k];
	  // produce a trial final state:
	  // insert the flavour that is not going to be replaced and the outgoing flavours of
	  // the already existing ffkernel with the other flavour as incoming flavour
	  multiset<Flavour> trial;
	  trial.insert(outs[2-j]);
	  for (size_t f=0;f<rkernel->GetFlavs().size();f++) trial.insert(rkernel->GetFlavs()[f]);
	  // check if we already have this combination
	  if (vetoed.find(trial)==vetoed.end()) {
	    // new combination - insert into table of vetoed combinations
	    vetoed.insert(trial);
	    Flavour_Vector flavs;
	    flavs.push_back(splitter);
	    for (multiset<Flavour>::iterator fit=trial.begin();fit!=trial.end();fit++)
	      flavs.push_back(*fit);
	    MakeKernels(flavs,type);
	  }
	}
      }
    }
  }
}

void Kernel_Constructor::MakePermutations(const size_t & order) {
  size_t now = 2, pos;
  m_permutations.resize(order-1);
  while (now<=order) {
    vector<size_t> perm;
    for (size_t i=0;i<now;i++) perm.push_back(i);
    pos = now;
    GeneratePermutation(perm,pos,now);
    now++;
  }
  //if (msg_LevelIsDebugging()) {
  for (size_t i=0;i<m_permutations.size();i++) {
    msg_Out()<<"####### Permutations for "<<(i+2)<<" partons ##############\n";
    size_t j=1;
    msg_Out()<<"        found "<<m_permutations[i].size()<<" entries:\n";
    for (list<vector<size_t> >::iterator lit=m_permutations[i].begin();
	 lit!=m_permutations[i].end();lit++) {
      msg_Out()<<" "<<(j++)<<"th permutation = {";
      for (size_t k=0;k<lit->size();k++) msg_Out()<<" "<<(*lit)[k];
      msg_Out()<<" }\n";
    }
  }
  //}
}

void Kernel_Constructor::
GeneratePermutation(std::vector<size_t> perm,const size_t & pos,const size_t & length) {
  if (pos==1) {
    vector<size_t> newperm = perm;
    m_permutations[length-2].push_back(newperm);
  }
  else {
    GeneratePermutation(perm,pos-1,length);
    for (size_t i=0;i<pos-1;i++) {
      if (i%2) swap(perm[i], perm[pos-1]);
      else     swap(perm[0], perm[pos-1]);
      GeneratePermutation(perm,pos-1,length);
    }
  }
}

void Kernel_Constructor::PrintKernels() {
  for (size_t i=1;i<5;i++) {
    if (p_kernels[i]->size()==0) continue;
    msg_Out()<<"--------------------------------------------------\n"
	     <<"--- Kernels for type = "<<i<<":\n";
    for (map<Flavour,Kernels *>::iterator ksit=p_kernels[i]->begin();
	 ksit!=p_kernels[i]->end();ksit++) {
      msg_Out()<<"--- flavour = "<<ksit->first<<"\n";
      for (Kernel_Vector::iterator kit=ksit->second->begin();
	   kit!=ksit->second->end();kit++) {
	msg_Out()<<(**kit);
      }
    }
    msg_Out()<<"--------------------------------------------------\n";
  }
}

