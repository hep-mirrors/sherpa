#include "COMIX/Phasespace/PS_Generator.H"

#include "COMIX/Phasespace/PS_Current.H"
#include "COMIX/Phasespace/PS_Vertex.H"
#include "COMIX/Phasespace/PS_Channel.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "COMIX/Main/Process_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Smart_Pointer.H"
#include "ATOOLS/Org/Smart_Pointer.C"

#include <iomanip>

using namespace COMIX;
using namespace ATOOLS;
using namespace PHASIC;

const double s_pwmin(1.0e-6);

namespace ATOOLS { template class SP(PS_Generator); }

PS_Generator::PS_Generator(Process_Base *const xs):
  p_xs(xs), m_n(0), m_zmode(1), m_pmsinit(0),
  m_thmass(0.0), m_chmass(0.0)
{
  Data_Reader read(" ",";","!","=");
  read.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(m_itmin,"CDXS_ITMIN")) m_itmin=5000;
  else msg_Info()<<METHOD<<"(): Set iteration minimum "<<m_itmin<<".\n";
  if (!read.ReadFromFile(m_itmax,"CDXS_ITMAX")) m_itmax=50000;
  else msg_Info()<<METHOD<<"(): Set iteration maximum "<<m_itmax<<".\n";
  if (!read.ReadFromFile(m_ecmode,"CDXS_ECMODE")) m_ecmode=2;
  else msg_Info()<<METHOD<<"(): Set extra channel mode "<<m_ecmode<<".\n";
  if (!read.ReadFromFile(m_chmass,"CDXS_PS_CHTH")) m_chmass=0.01;
  else msg_Info()<<METHOD<<"(): Set channel mass threshold "<<m_chmass<<".\n";
  m_chmass*=rpa->gen.Ecms();
  p_xs->ConstructPSVertices(this);
  AddSC();
  AddSTCC();
}

PS_Generator::~PS_Generator()
{
  CleanUp();
}

void PS_Generator::CleanUp()
{
  for (size_t i(0);i<m_cur.size();++i) 
    for (size_t j(0);j<m_cur[i].size();++j) delete m_cur[i][j]; 
  m_n=0;
  m_cl=Int_Matrix();
  m_cur=Current_Matrix();
  m_ctt=Current_Vector();
  m_swidths=m_smasses=Double_Vector();
  m_cmap=CB_MMap();
  m_cbmap=CB_Map();
  m_thmass=0.0;
}

size_t PS_Generator::NChannels() const
{
  size_t nch(0);
  for (size_t i(0);i<m_cur.size();++i)
    nch+=m_cur[i].size();
  return nch;
}    

void PS_Generator::SetColors(const Int_Vector &rc,
			     const Int_Vector &ac)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  for (size_t i(0);i<m_n;++i) {
    m_cl[i][0]=rc[i];
    m_cl[i][1]=ac[i];
#ifdef DEBUG__BG
    msg_Debugging()<<"m_cl["<<i<<"][0] = "<<m_cl[i][0]
		   <<", m_cl["<<i<<"][1] = "<<m_cl[i][1]<<"\n";
#endif
  }
}

void PS_Generator::CalcJL()
{
  for (size_t i(0);i<m_cur[1].size();++i) 
    m_cur[1][i]->ConstructJ(Vec4D(),0,m_cl[i][0],m_cl[i][1],0);
  if (m_zmode>0) {
    for (size_t n(2);n<m_n;++n) {
      for (size_t i(0);i<m_cur[n].size();++i) 
	m_cur[n][i]->Evaluate();
    }
    for (size_t n(m_n-2);n>=2;--n)
      for (size_t i(0);i<m_cur[n].size();++i)
	m_cur[n][i]->ResetZero();
  }
}

bool PS_Generator::Evaluate()
{
  if (m_zmode>0) {
    PHASIC::Process_Base *cur(p_xs->Process());
    while ((*cur)[0]!=cur) {
      double sum(0.0);
      Double_Vector sums;
      std::vector<PHASIC::Process_Base*> pss;
      for (size_t i(0);i<cur->Size();++i)
	if (!(*cur)[i]->Get<PHASIC::Single_Process>()->Zero()) {
	  if ((*cur)[i]->Integrator()->SumSqr()==0) sum+=1.0;
	  else sum+=(*cur)[i]->Integrator()->SumSqr();
	  sums.push_back(sum);
	  pss.push_back((*cur)[i]);
	}
      double disc(sum*ran->Get());
      for (size_t i(0);i<pss.size();++i)
	if (disc<=sums[i]) {
	  cur=pss[i];
	  break;
	}
    }
    SP(Color_Integrator) ci(cur->Integrator()->ColorIntegrator());
    if (ci==NULL) 
      THROW(fatal_error,"No color integrator for "+cur->Name());
    SetColors(ci->I(),ci->J());
  }
  CalcJL();
  return true;
}

bool PS_Generator::AddCurrent
(Current *const ref,const Flavour &fl,
 const size_t &n,const int mode,
 const double &m,const double &w,Current *const scc)
{
  Current_Key ckey(fl,NULL,0);
  PS_Current *cur(new PS_Current(ckey));
  if (mode&1) {
    cur->SetMass(m);
    cur->SetWidth(w);
  }
  cur->SetId(ref->Id());
  cur->SetKey(m_cur[n].size());
  cur->SetDirection(ref->Direction());
  cur->SetCut(ref->Cut());
  cur->SetOnShell(ref->OnShell());
  // cur->SetOrder(ref->Order());
  cur->SetNTChannel(ref->NTChannel());
  cur->SetSCC(scc);
  m_cur[n].push_back(cur);
  if (scc==NULL) m_tccs[cur->CId()].push_back(cur);
  m_cmap.insert(CB_Pair(ref,cur));
  m_cbmap.insert(CB_Pair(cur,ref));
  bool found(false);
  for (size_t i(0);i<m_ctt.size();++i)
    if (m_ctt[i]->PSInfo()==cur->PSInfo()) {
      found=true;
      break;
    }
  if (!found) m_ctt.push_back(cur);
  return found;
}

int PS_Generator::DecayType(const Current *jc,
			    const Current *ja,const Current *jb) const
{
  if (jc->CId()==(1<<m_n)-1-3) return 0;
  if ((jc->CId()&3)==1 || (jc->CId()&3)==2 ||
      jc->Flav().Mass()!=0.0 ||
      ja->Flav().Mass()!=0.0 || jb->Flav().Mass()!=0.0) return 0;
  if (ja->Flav()==jc->Flav() &&
      jb->Flav()==jc->Flav()) return 2|4;
  if (ja->Flav()==jb->Flav().Bar()) return 0;// temporary (->1)
  if (ja->Flav()==jc->Flav()) return 2;
  if (jb->Flav()==jc->Flav()) return 4;
  return 0;
}

class CB_PSSort {
public:

  bool operator()(const Vertex_Key &ka,const Vertex_Key &kb) const
  {
    if (ka.m_j.front()<kb.m_j.front()) return true;
    if (ka.m_j.front()>kb.m_j.front()) return false;
    if (ka.m_j.back()<kb.m_j.back()) return true;
    if (ka.m_j.back()>kb.m_j.back()) return false;
    if (ka.p_c<kb.p_c) return true;
    if (ka.p_c>kb.p_c) return false;
    return false;
  }

};

void PS_Generator::AddSubChannel
(NLO_subevtlist *const subs,const Vertex_Key &vkey)
{
  if (subs==NULL || subs->empty()) return;
  NLO_subevt *csub(NULL);
  for (int l=0;l<2;++l) {
    if (vkey.m_j[l]->Id().size()>1) continue;
    size_t ii(vkey.m_j[l]->Id().front());
    size_t jj(vkey.m_j[1-l]->Id().front());
    if (vkey.m_j[1-l]->Id().size()>1) {
      size_t idj=(1<<m_n)-1-(1<<ii)-vkey.m_j[1-l]->CId();
      if (IdCount(idj)>1) {
	idj=(1<<m_n)-1-vkey.p_c->CId();
	if (IdCount(idj)>1) continue;
      }
      jj=ID(idj).front();
    }
    for (size_t k(0);k<subs->size()-1;++k)
      if ((ii==(*subs)[k]->m_i && jj==(*subs)[k]->m_j) ||
	  (jj==(*subs)[k]->m_i && ii==(*subs)[k]->m_j)) {
	csub=(*subs)[k];
	break;
      }
    if (csub) break;
  }
  if (csub==NULL) return;
  Vertex_Key *dummy(Vertex_Key::New(Current_Vector(),NULL,NULL));
  PS_Vertex *vtx(new PS_Vertex(*dummy));
  dummy->Delete();
  vtx->AddJ(vkey.m_j);
  vtx->SetJC(vkey.p_c);
  vtx->SetDip(csub);
}

bool PS_Generator::Construct(Amplitude *const ampl,NLO_subevtlist *const subs)
{
  m_pmsinit=0;
  Current_Matrix curs(ampl->Currents());
  msg_Debugging()<<METHOD<<"(): '"<<ampl<<"' {\n";
  {
    msg_Indent();
  if (m_n>0) {
    if (m_n!=curs.size()) 
      THROW(fatal_error,"Invalid number of external particles");
  }
  else {
    m_n=curs.size();
    m_cl.resize(m_n,Int_Vector(2));
    m_cur.resize(m_n);
  }
  Vertex_Key *dummy(Vertex_Key::New(Current_Vector(),NULL,NULL));
  for (size_t n(1);n<m_n;++n) {
    for (size_t j(0);j<curs[n].size();++j) {
      if (curs[n][j]->Sub() ||
	  curs[n][j]->Flav().IsDummy()) continue;
      bool found(false);
      for (size_t i(0);i<m_cur[n].size();++i)
	if (m_cur[n][i]->Id()==curs[n][j]->Id() &&
	    (n==1 ||
	     (m_cur[n][i]->Flav().Mass()==curs[n][j]->Flav().Mass() &&
	      m_cur[n][i]->Flav().Width()==curs[n][j]->Flav().Width() &&
	      // m_cur[n][i]->Order()==curs[n][j]->Order() &&
	      m_cur[n][i]->NTChannel()==curs[n][j]->NTChannel()))) {
	  Current *ref(m_cbmap[m_cur[n][i]]);
	  for (CB_MMap::const_iterator cit(m_cmap.lower_bound(ref));
	       cit!=m_cmap.upper_bound(ref);++cit) {
	    m_cmap.insert(CB_Pair(curs[n][j],cit->second));
	  }
	  found=true;
	  break;
	}
      if (found) {
	for (CB_MMap::const_iterator cit(m_cmap.lower_bound(curs[n][j]));
	     cit!=m_cmap.upper_bound(curs[n][j]);++cit) {
	  std::set<Vertex_Key,CB_PSSort> v3;
	  const Vertex_Vector &in(curs[n][j]->In());
	  const Vertex_Vector &rin(cit->second->In());
	  for (size_t i(0);i<in.size();++i) {
	    if (curs[n][j]->In()[i]->J().size()>2) continue;
	    Current *ja(curs[n][j]->In()[i]->J(0));
	    Current *jb(curs[n][j]->In()[i]->J(1));
	    if (ja->PSInfo()<jb->PSInfo()) std::swap<Current*>(ja,jb);
	    for (CB_MMap::const_iterator ait(m_cmap.lower_bound(ja));
		 ait!=m_cmap.upper_bound(ja);++ait)
	      for (CB_MMap::const_iterator bit(m_cmap.lower_bound(jb));
		   bit!=m_cmap.upper_bound(jb);++bit) {
		Current_Vector jj(2);
		jj[0]=ait->second;
		jj[1]=bit->second;
		Vertex_Key *vkey(Vertex_Key::New(jj,cit->second,NULL));
		int type(DecayType(cit->second,ait->second,bit->second)), mtype(0);
		bool vf(false);
		for (size_t k(0);k<rin.size();++k)
		  if (jj[0]==rin[k]->J(0) && jj[1]==rin[k]->J(1) &&
		      (type==((PS_Vertex*)rin[k])->Type() ||
		       (type&((PS_Vertex*)rin[k])->Type()))) {
		    mtype|=((PS_Vertex*)rin[k])->Type();
		    vf=true;
		  }
		if ((vf && type==mtype) || v3.find(*vkey)!=v3.end()) continue;
		v3.insert(*vkey);
		PS_Vertex *vtx(new PS_Vertex(*dummy));
		vtx->AddJ(vkey->m_j);
		vtx->SetJC(vkey->p_c);
		vtx->SetType(type);
		if (type==(2|4)) {
		  vtx->SetType(2);
		  vtx = new PS_Vertex(*dummy);
		  vtx->AddJ(vkey->m_j);
		  vtx->SetJC(vkey->p_c);
		  vtx->SetType(4);
		}
		AddSubChannel(subs,*vkey);
		vkey->Delete();
	      }
	  }
	}
	continue;
      }
      const Vertex_Vector &in(curs[n][j]->In());
      bool valid(in.empty());
      for (size_t i(0);i<in.size();++i)
	if (in[i]->J().size()<3) {
	  valid=true;
	  break;
	}
      if (!valid) continue;
      curs[n][j]->Print();
      if (curs[n][j]->Flav().Width()<s_pwmin &&
	  !curs[n][j]->Cut() && curs[n][j]->Flav().Mass()>0.0 &&
	  curs[n][j]->Flav().Mass()<m_chmass && n>1 && n<m_n-1) 
	AddCurrent(curs[n][j],curs[n][j]->Flav(),n,1);
      else AddCurrent(curs[n][j],curs[n][j]->Flav(),n);
      std::set<Vertex_Key,CB_PSSort> v3;
      for (size_t i(0);i<in.size();++i) {
	if (curs[n][j]->In()[i]->J().size()>2) continue;
	Current *ja(curs[n][j]->In()[i]->J(0));
	Current *jb(curs[n][j]->In()[i]->J(1));
	if (ja->PSInfo()<jb->PSInfo()) std::swap<Current*>(ja,jb);
	for (CB_MMap::const_iterator ait(m_cmap.lower_bound(ja));
	     ait!=m_cmap.upper_bound(ja);++ait)
	  for (CB_MMap::const_iterator bit(m_cmap.lower_bound(jb));
	       bit!=m_cmap.upper_bound(jb);++bit)
	    for (CB_MMap::const_iterator cit(m_cmap.lower_bound(curs[n][j]));
		 cit!=m_cmap.upper_bound(curs[n][j]);++cit) {
	      Current_Vector jj(2);
	      jj[0]=ait->second;
	      jj[1]=bit->second;
	      Vertex_Key *vkey(Vertex_Key::New(jj,cit->second,NULL));
	      if (v3.find(*vkey)!=v3.end()) continue;
	      v3.insert(*vkey);
	      PS_Vertex *vtx(new PS_Vertex(*dummy));
	      vtx->AddJ(vkey->m_j);
	      vtx->SetJC(vkey->p_c);
	      int type(DecayType(curs[n][j],ja,jb));
	      vtx->SetType(type);
	      if (type==(2|4)) {
		vtx->SetType(2);
		vtx = new PS_Vertex(*dummy);
		vtx->AddJ(vkey->m_j);
		vtx->SetJC(vkey->p_c);
		vtx->SetType(4);
	      }
	      AddSubChannel(subs,*vkey);
	      vkey->Delete();
	    }
      }
      for (CB_MMap::const_iterator cit(m_cmap.lower_bound(curs[n][j]));
	   cit!=m_cmap.upper_bound(curs[n][j]);++cit)
	cit->second->Print();
    }
  }
  dummy->Delete();
  }
  msg_Debugging()<<"}\n";
  for (size_t j(m_n-2);j>1;--j)
    for (Current_Vector::iterator cit(m_cur[j].begin());
	 cit!=m_cur[j].end();++cit)
      if ((*cit)->Dangling()) {
#ifdef DEBUG__BG
	msg_Debugging()<<"  delete current "<<**cit
		       <<", "<<(*cit)->Dangling()<<"\n";
#endif
	for (Current_Vector::iterator ctit(m_ctt.begin());
	     ctit!=m_ctt.end();++ctit)
	  if (*ctit==*cit) {
	    m_ctt.erase(ctit);
	    break;
	  }
	delete *cit;
	cit=--m_cur[j].erase(cit);
      }
  msg_Debugging()<<METHOD<<"(): Phase space statistics (n="
		 <<m_n<<") {\n  level currents vertices\n"<<std::right;
  double ismass(0.0);
  size_t csum(0), vsum(0), itmin(1);
  std::map<size_t,double> mmin;
  for (size_t i(1);i<m_n;++i) {
    csum+=m_cur[i].size();
    size_t cvsum(0);
    for (size_t j(0);j<m_cur[i].size();++j) {
      size_t cid(m_cur[i][j]->CId());
      double mass(i==1?m_cur[i][j]->Flav().Mass():m_cur[i][j]->Mass());
      if (mass>rpa->gen.Ecms()) mass=0.0;
      if ((cid&3)==1 || (cid&3)==2) {
	if (i==1) ismass+=mass;
	mass=0.0;
      }
      for (size_t k(0);k<m_cur[i][j]->In().size();++k)
	mass=Max(mass,mmin[m_cur[i][j]->In()[k]->J(0)->CId()]+
		 mmin[m_cur[i][j]->In()[k]->J(1)->CId()]);
      if (mmin.find(cid)==mmin.end()) mmin[cid]=mass;
      else mmin[cid]=Max(mmin[cid],mass);
      if (i==m_n-1) {
	m_thmass=Max(m_thmass,Max(mass,ismass));
      }
#ifdef DEBUG__BG
      msg_Debugging()<<"  set min mass "
		     <<m_cur[i][j]->PSInfo()<<" -> "<<mass<<"\n";
#endif
      if ((cid==3 || cid==(size_t)((1<<m_n)-1-3)) &&
	  m_cur[i][j]->Mass()>0.0) {
	bool found(false);
	for (size_t k(0);k<m_smasses.size();++k)
	  if (m_smasses[k]==m_cur[i][j]->Mass()&&
	      m_swidths[k]==m_cur[i][j]->Width()) {
	    found=true;
	    break;
	  }
	if (!found) {
	  m_smasses.push_back(m_cur[i][j]->Mass());
	  m_swidths.push_back(m_cur[i][j]->Width());
	}
      }
#ifdef DEBUG__BG
      msg_Debugging()<<"  set min mass "
		     <<m_cur[i][j]->PSInfo()<<" -> "<<mass<<"\n";
#endif
      cvsum+=m_cur[i][j]->NIn();
    }
    msg_Debugging()<<"  "<<std::setw(5)<<i<<" "<<std::setw(8)
		   <<m_cur[i].size()<<" "<<std::setw(8)<<cvsum<<"\n";
    vsum+=cvsum;
    if (i>1 && m_cur[i].size()>0) itmin*=cvsum/m_cur[i].size();
  }
  msg_Debugging()<<std::left<<"} -> "<<csum<<" currents, "
		 <<vsum<<" vertices\n";
  itmin=Min(Max(itmin,m_itmin),m_itmax);
  msg_Tracking()<<METHOD<<"(): Set iteration minimum "<<itmin<<".\n";
  p_xs->Process()->Integrator()->SetItMin(itmin);
  msg_Tracking()<<METHOD<<"(): Masses {\n"
		<<"  threshold  : "<<m_thmass<<"\n";
  for (size_t i(0);i<m_smasses.size();++i)
    msg_Tracking()<<"  resonance "<<i<<": "
		  <<m_smasses[i]<<" / "<<m_swidths[i]<<"\n";
  msg_Tracking()<<"}\n";
  return true;
}

void PS_Generator::AddSC()
{
#ifdef DEBUG__BG
  DEBUG_FUNC("");
#endif
  for (size_t n(2);n<m_n-1;++n) {
    size_t oldsize(m_cur[n].size());
    for (size_t j(0);j<oldsize;++j) {
      if (((PS_Current*)m_cur[n][j])->Dip()) continue;
      for (size_t i(0);i<m_cur[n][j]->NIn();++i) {
	NLO_subevt *dip(((PS_Vertex*)m_cur[n][j]->In()[i])->Dip());
	if (dip) {
 	  delete m_cur[n][j]->In()[i];
	  ((Vertex_Vector*)&m_cur[n][j]->In())->erase
	    (((Vertex_Vector*)&m_cur[n][j]->In())->begin()+i);
	  AddExtraCurrent(m_cur[n][j],n,m_cur[n][j]->Flav().Mass(),
			  m_cur[n][j]->Flav().Width());
	  ((PS_Current*)m_cur[n].back())->SetDip(dip);
	  break;
	}
      }
    }
  }
}

void PS_Generator::AddSTCC()
{
#ifdef DEBUG__BG
  DEBUG_FUNC("");
#endif
  for (size_t n(2);n<m_n;++n) {
    size_t oldsize(m_cur[n].size());
    for (size_t j(0);j<oldsize;++j) {
      if (((PS_Current*)m_cur[n][j])->SCC()) continue;
      size_t cid(m_cur[n][j]->CId());
      if ((cid&3)==0 || (cid&3)==3) continue;
      std::set<std::string> added;
      size_t pid(~3&cid);
      if (IdCount(pid)>1) {
	TCC_Map::const_iterator it(m_tccs.find(pid));
	if (it==m_tccs.end()) continue;
	for (size_t i(0);i<it->second.size();++i) {
	  if (added.empty()) 
	    ((PS_Current*)m_cur[n][j])->SetSCC(it->second[i]);
	  else if (added.find(it->second[i]->PSInfo())==added.end())
	    AddExtraCurrent(m_cur[n][j],n,m_cur[n][j]->Flav().Mass(),
			    m_cur[n][j]->Flav().Width(),it->second[i]);
	  added.insert(it->second[i]->PSInfo());
	}
      }
    }
  }
}

void PS_Generator::AddExtraCurrent
(Current *const cur,const size_t &n,
 const double &m,const double &w,Current *const scc)
{
  AddCurrent(cur,cur->Flav(),n,1,m,w,scc);
#ifdef DEBUG__BG
  msg_Debugging()<<"  Add "<<m_cur[n].back()->PSInfo()
		 <<(scc?" ("+scc->PSInfo()+") ":"")<<" {\n";
#endif
  Vertex_Key *dummy(Vertex_Key::New(Current_Vector(),NULL,NULL));
  const Vertex_Vector &in(cur->In());
  for (size_t i(0);i<in.size();++i) {
    PS_Vertex *vtx(new PS_Vertex(*dummy));
    vtx->AddJ(in[i]->J());
    vtx->SetJC(m_cur[n].back());
    vtx->SetDip(((PS_Vertex*)in[i])->Dip());
    vtx->SetType(((PS_Vertex*)in[i])->Type());
#ifdef DEBUG__BG
    msg_Debugging()<<"    "<<*vtx<<"\n";
#endif
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"  } <-> {\n";
#endif
  const Vertex_Vector &out(cur->Out());
  for (size_t i(0);i<out.size();++i) {
    Current_Vector j(out[i]->J());
    if (j[0]==cur) j[0]=m_cur[n].back();
    else j[1]=m_cur[n].back();
    PS_Vertex *vtx(new PS_Vertex(*dummy));
    vtx->AddJ(j);
    vtx->SetJC(out[i]->JC());
    vtx->SetDip(((PS_Vertex*)out[i])->Dip());
    vtx->SetType(((PS_Vertex*)out[i])->Type());
#ifdef DEBUG__BG
    msg_Debugging()<<"    "<<*vtx<<"\n";
#endif
  }
  dummy->Delete();
#ifdef DEBUG__BG
  msg_Debugging()<<"  }\n";
#endif
}

void PS_Generator::SetPrefMasses(Cut_Data *const cuts)
{
  if (m_pmsinit) return;
  m_pmsinit=1;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"() {\n";
#endif
  std::map<size_t,double> mmin;
  for (size_t n(2);n<m_n-2;++n) {
    size_t oldsize(m_cur[n].size());
    for (size_t j(0);j<oldsize;++j) {
      size_t cid(m_cur[n][j]->CId());
      size_t pid((cid&3)?(1<<m_n)-1-cid:cid);
      double psmin(sqrt(cuts->Getscut(PSId(pid))));
      double mass(Max(psmin,m_cur[n][j]->Mass()));
      if (m_cur[n][j]->OnShell()) {
	mmin[cid]=mass;
	continue;
      }
      if ((cid&3)==1 || (cid&3)==2) mass=0.0;
      for (size_t k(0);k<m_cur[n][j]->In().size();++k) {
	mass=Max(mass,mmin[m_cur[n][j]->In()[k]->J(0)->CId()]+
		 mmin[m_cur[n][j]->In()[k]->J(1)->CId()]);
      }
      if (mmin.find(cid)==mmin.end()) mmin[cid]=mass;
      else mmin[cid]=Max(mmin[cid],mass);
      if (m_cur[n][j]->Cut() || (cid&3)==1 || (cid&3)==2) continue;
      if ((m_ecmode&1) && mass>m_chmass &&
	  mass<rpa->gen.Ecms() && !IsEqual(mass,m_cur[n][j]->Mass()))
	AddExtraCurrent(m_cur[n][j],n,mass,0.0);
      if ((m_ecmode&2) && m_cur[n][j]->Mass()>m_chmass &&
	  m_cur[n][j]->Width()>s_pwmin && IdCount(pid)>2)
	AddExtraCurrent(m_cur[n][j],n,m_cur[n][j]->Mass(),0.0);
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}
