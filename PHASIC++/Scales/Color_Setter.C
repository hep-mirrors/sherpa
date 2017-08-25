#include "PHASIC++/Scales/Color_Setter.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flow.H"

#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;

size_t s_clmaxtrials(900);

Color_Setter::Color_Setter(const int mode): m_cmode(mode)
{
  m_pmap[nlo_type::lo] = new StringProcess_Map();
}

Color_Setter::~Color_Setter()
{
  for (Flav_ME_Map::const_iterator xsit(m_xsmap.begin());
       xsit!=m_xsmap.end();++xsit) delete xsit->second;
  for (size_t i(0);i<m_procs.size();++i) delete m_procs[i];
  delete m_pmap[nlo_type::lo];
}

bool Color_Setter::SetRandomColors(Cluster_Amplitude *const ampl)
{
  size_t trials(0), vc(1);
  std::vector<ColorID> oc(ampl->Legs().size());
  for (size_t i(0);i<ampl->Legs().size();++i) {
    oc[i]=ampl->Leg(i)->Col();
    if (ampl->Leg(i)->Flav().Strong() &&
	oc[i].m_i==0 && oc[i].m_j==0) vc=0;
  }
  if (vc==0) {
    // select new color configuration
    SP(Color_Integrator) colint(p_xs->Integrator()->ColorIntegrator());
    while (!colint->GeneratePoint());
    PHASIC::Int_Vector ni(colint->I()), nj(colint->J());
    for (size_t i(0);i<ampl->Legs().size();++i)
      oc[i]=ColorID(ni[i],nj[i]);
  }
  msg_Debugging()<<*ampl<<"\n";
  for (;trials<s_clmaxtrials;++trials) {
    bool sing(false);
    std::set<size_t> cs;
    PHASIC::Int_Vector ci(ampl->Legs().size()), cj(ampl->Legs().size());
    for (size_t i(0);i<ampl->Legs().size();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      int col(oc[i].m_i);
      if (col==0) continue;
      std::vector<size_t> js;
      for (size_t j(0);j<ampl->Legs().size();++j)
	if (i!=j && oc[j].m_j==col && cs.find(j)==cs.end()) 
	  js.push_back(j);
      if (js.empty()) {
	msg_Debugging()<<"color singlet "<<*cl<<"\n";
	sing=true;
	break;
      }
      size_t j(js[Min((size_t)(ran->Get()*js.size()),js.size()-1)]);
      cs.insert(j);
      Cluster_Leg *cp(ampl->Leg(j));
      size_t nc(Flow::Counter());
      cl->SetCol(ColorID(ci[i]=nc,cl->Col().m_j));
      cp->SetCol(ColorID(cp->Col().m_i,cj[j]=nc));
      msg_Debugging()<<"set color "<<nc<<"\n";
      msg_Debugging()<<"  "<<*cl<<"\n";
      msg_Debugging()<<"  "<<*cp<<"\n";
    }
    if (!sing) {
      for (size_t i(0);i<ampl->Legs().size();++i)
	ampl->Leg(i)->SetCol(ColorID(ci[i],cj[i]));
      double csum(p_xs->Differential(*ampl,1|2|4));
      msg_Debugging()<<"sc: csum = "<<csum<<"\n";
      if (csum!=0.0) {
	CI_Map &cmap(ampl->ColorMap());
	for (size_t i(0);i<ampl->Legs().size();++i)
	  if (oc[i].m_i!=0) cmap[ci[i]]=oc[i].m_i;
	break;
      }
    }
    if ((trials%9==0 && trials>0) || sing) {
      // select new color configuration
      SP(Color_Integrator) colint(p_xs->Integrator()->ColorIntegrator());
      while (!colint->GeneratePoint());
      PHASIC::Int_Vector ni(colint->I()), nj(colint->J());
      for (size_t i(0);i<ampl->Legs().size();++i)
	oc[i]=ColorID(ni[i],nj[i]);
    }
  }
  if (trials>=s_clmaxtrials) {
    msg_Error()<<METHOD<<"(): No solution."<<std::endl;
    return false;
  }
  return true;
}

bool Color_Setter::SetLargeNCColors(Cluster_Amplitude *const ampl)
{
  Process_Base *proc(ampl->Proc<Process_Base>());
  if (proc==NULL || proc->AllProcs()==NULL) return false;
  StringProcess_Map *pm((*proc->AllProcs())[nlo_type::lo]);
  Process_Base::SortFlavours(ampl);
  std::string name(Process_Base::GenerateName(ampl));
  StringProcess_Map::const_iterator pit(pm->find(name));
  p_xs=NULL;
  if (pit!=pm->end() && pit->second->
      Integrator()->ColorIntegrator()!=NULL) p_xs=pit->second;
  if (p_xs==NULL) {
    pm=m_pmap[nlo_type::lo];
    if ((pit=pm->find(name))!=pm->end()) p_xs=pit->second;
    else {
      MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process",true);
      My_In_File::OpenDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
      Process_Info pi;
      pi.m_megenerator="Comix";
      for (size_t i(0);i<ampl->NIn();++i) {
	Flavour fl(ampl->Leg(i)->Flav().Bar());
	if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
	pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
      }
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
	Flavour fl(ampl->Leg(i)->Flav());
	if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
	pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
      }
      PHASIC::Process_Base *proc=
	ampl->Proc<Process_Base>()->
	Generator()->Generators()->InitializeProcess(pi,false);
      if (proc==NULL) {
	My_In_File::CloseDB
	  (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
	(*pm)[name]=NULL;
	return false;
      }
      m_procs.push_back(proc);
      Selector_Key skey(NULL,NULL,true);
      proc->SetSelector(skey);
      proc->SetScale
	(Scale_Setter_Arguments
	 (MODEL::s_model,"VAR{"+ToString(sqr(rpa->gen.Ecms()))+"}","Alpha_QCD 1"));
      proc->SetKFactor(KFactor_Setter_Arguments("NO"));
      proc->Get<COMIX::Process_Base>()->Tests();
      proc->FillProcessMap(&m_pmap);
      My_In_File::CloseDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
      if ((pit=pm->find(name))==pm->end()) THROW(fatal_error,"Internal error");
      p_xs=pit->second;
    }
    if (p_xs==NULL) return false;
  }
  DEBUG_FUNC(p_xs->Name());
  msg_Debugging()<<*ampl<<"\n";
  SP(Color_Integrator) colint(p_xs->Integrator()->ColorIntegrator());
  PHASIC::Int_Vector ci(colint->I()), cj(colint->J());
  bool sol(false);
  switch (m_cmode) {
  case 1: {
    sol=SetRandomColors(ampl);
    break;
  } 
  default:
    THROW(fatal_error,"Invalid colour setting mode");
  }
  colint->SetI(ci);
  colint->SetJ(cj);
  return sol;
}

void Color_Setter::SetColors(ATOOLS::Cluster_Amplitude *ampl)
{
  bool cs(true);
  Vec4D_Vector moms(ampl->Legs().size());
  Flavour_Vector fl(ampl->Legs().size());
  for (int i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *l(ampl->Leg(i));
    moms[i]=i<ampl->NIn()?-l->Mom():l->Mom();
    fl[i]=i<ampl->NIn()?l->Flav().Bar():l->Flav();
  }
  Flav_ME_Map::const_iterator xit(m_xsmap.find(fl));
  if (xit==m_xsmap.end() && ampl->Legs().size()==4) {
    Process_Info pi;
    pi.m_megenerator="Internal";
    pi.m_maxcpl[0]=pi.m_mincpl[0]=ampl->OrderQCD();
    pi.m_maxcpl[1]=pi.m_mincpl[1]=ampl->OrderEW();
    for (size_t i(0);i<ampl->NIn();++i)
      pi.m_ii.m_ps.push_back(Subprocess_Info(fl[i]));
    for (size_t i(ampl->NIn());i<fl.size();++i)
      pi.m_fi.m_ps.push_back(Subprocess_Info(fl[i]));
    EXTRAXS::ME2_Base *me2=dynamic_cast<EXTRAXS::ME2_Base*>
      (PHASIC::Tree_ME2_Base::GetME2(pi));
    m_xsmap[fl]=me2;
    xit=m_xsmap.find(fl);
  }
  if (xit!=m_xsmap.end() && xit->second!=NULL) {
    bool test(xit->second->SetColours(moms));
    for (size_t i(0);i<fl.size();++i) {
      ColorID c(xit->second->Colours()[i][0],
		xit->second->Colours()[i][1]);
      if (i<ampl->NIn())c=ColorID(c.m_j,c.m_i);
      ampl->Leg(i)->SetCol(c);
    }
  }
  else {
    if (m_cmode==0 || !SetLargeNCColors(ampl)) {
      std::vector<int> tids, atids;
      for (size_t i(0);i<ampl->Legs().size();++i)
	if (ampl->Leg(i)->Flav().StrongCharge()>0) {
	  tids.push_back(i);
	  if (ampl->Leg(i)->Flav().StrongCharge()==8)
	    atids.push_back(i);
	  else ampl->Leg(i)->SetCol(ColorID(ampl->Leg(i)->Col().m_i,0));
	}
	else if (ampl->Leg(i)->Flav().StrongCharge()<0) {
	  ampl->Leg(i)->SetCol(ColorID(0,ampl->Leg(i)->Col().m_j));
	  atids.push_back(i);
	}
	else {
	  ampl->Leg(i)->SetCol(ColorID(0,0));
	}
      while (true) {
	std::random_shuffle(atids.begin(),atids.end(),*ran);
	size_t i(0);
	for (;i<atids.size();++i) if (atids[i]==tids[i]) break;
	if (i==atids.size()) break;
      }
      for (size_t i(0);i<tids.size();++i) {
	int cl(Flow::Counter());
	ampl->Leg(tids[i])->SetCol(ColorID(cl,ampl->Leg(tids[i])->Col().m_j));
	ampl->Leg(atids[i])->SetCol(ColorID(ampl->Leg(atids[i])->Col().m_i,cl));
      }
    }
  }
  for (Cluster_Amplitude *campl(ampl->Prev());campl;campl=campl->Prev()) {
    Cluster_Amplitude *next(campl->Next()); 
    Cluster_Leg *lij=NULL;
    for (size_t i(0);i<next->Legs().size();++i)
      if (next->Leg(i)->K()) {
	lij=next->Leg(i);
	break;
      }
    if (lij==NULL) THROW(fatal_error,"Invalid amplitude");
    Cluster_Leg *li=NULL, *lj=NULL;
    for (size_t i(0);i<campl->Legs().size();++i) {
      if (campl->Leg(i)->Id()&lij->Id()) {
	if (li==NULL) li=campl->Leg(i);
	else if (lj==NULL) lj=campl->Leg(i);
	else THROW(fatal_error,"Invalid splitting");
      }
      else {
	campl->Leg(i)->SetCol(next->IdLeg(campl->Leg(i)->Id())->Col());
      }
    }
    Cluster_Amplitude::SetColours(lij,li,lj);
  }
}
