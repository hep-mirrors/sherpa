#include "COMIX/Cluster/Cluster_Algorithm.H"

#include "COMIX/Main/Single_Process.H"
#include "COMIX/Main/Single_Dipole_Term.H"
#include "COMIX/Amplitude/Amplitude.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "COMIX/Phasespace/PS_Channel.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

using namespace COMIX;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm(ATOOLS::Mass_Selector *const ms):
  p_ms(ms), p_ampl(NULL), p_clus(NULL),
  m_lfrac(0.0)
{
  Data_Reader read(" ",";","#","=");
  m_corecheck=read.GetValue<int>("COMIX_CLUSTER_CORE_CHECK",0);
  m_nocluster=read.GetValue<int>("COMIX_NO_CLUSTER",0);
}

Cluster_Algorithm::~Cluster_Algorithm()
{
}

ColorID Cluster_Algorithm::GetColor(Current *const j,
				    Current *const fcur) const
{
  for (size_t i(0);i<j->J().size();++i) {
    const CObject_Vector &cs(j->J()[i]);
    if (!cs.empty()) {
      ColorID col((*cs.front())(0),(*cs.front())(1));
      if (j==fcur) col=col.Conj();
      return col;
    }
  }
  return ColorID();
}

ColorID_Vector Cluster_Algorithm::Connected
(const ColorID &ci,const ColorID &cj) const
{
  ColorID_Vector cij;
  int sj(cj.m_i==cj.m_j);
  int mij(ci.m_i==cj.m_j), mji(ci.m_j==cj.m_i);
  if (ci.m_i && ci.m_j) {
    if (cj.m_i && cj.m_j) {
      if (sj && mij && mji) return cij;
      if (mij) cij.push_back(ColorID(cj.m_i,ci.m_j));
      if (mji) cij.push_back(ColorID(ci.m_i,cj.m_j));
      return cij;
    }
    return Connected(cj,ci);
  }
  else if (ci.m_i) {
    if (cj.m_i && cj.m_j) {
      if (sj) cij.push_back(ci);
      if (mij) cij.push_back(ColorID(cj.m_i,0));
      return cij;
    }
    if (cj.m_i) return cij;
    if (!mij) cij.push_back(ColorID(ci.m_i,cj.m_j));
    else {
      for (int i(Color_Calculator::CIMin());
	   i<=Color_Calculator::CIMax();++i)
	cij.push_back(ColorID(i,i));
    }
    return cij;
  }
  else if (ci.m_j) {
    if (cj.m_i && cj.m_j) {
      if (sj) cij.push_back(ci);
      if (mji) cij.push_back(ColorID(0,cj.m_j));
      return cij;
    }
    if (cj.m_j) return cij;
    if (!mji) cij.push_back(ColorID(cj.m_i,ci.m_j));
    else {
      for (int i(Color_Calculator::CIMin());
	   i<=Color_Calculator::CIMax();++i)
	cij.push_back(ColorID(i,i));
    }
    return cij;
  }
  if (cj.m_i && cj.m_j) {
    cij.push_back(cj);
    return cij;
  }
  return cij;
}

ColorID_Vector Cluster_Algorithm::Connected
(const ColorID_Vector &ci,const ColorID_Vector &cj,
 const ColorID_Vector &ck,const Flavour &mo) const
{
  ColorID_Vector c;
  for (size_t i(0);i<ci.size();++i)
    for (size_t j(0);j<cj.size();++j)
      for (size_t k(0);k<ck.size();++k) {
	int mij(ci[i].m_i==cj[j].m_j), mji(cj[j].m_i==ci[i].m_j);
	if (!mo.Strong() && mij && mji) {
	  if (!(ck[k].m_i&&ck[k].m_j)) c.push_back(ColorID(0,0));
	  continue;
	}
	if (ck[k].m_i==0 && ck[k].m_j==0) continue;
	ColorID_Vector ckj(Connected(ck[k],cj[j]));
	if (ckj.size()) {
	  ColorID_Vector cij(Connected(ci[i],cj[j]));
	  c.insert(c.end(),cij.begin(),cij.end());
	}
      }
  return c;
}

CParam Cluster_Algorithm::GetMeasure
(const size_t &idi,const size_t &idj,const size_t &idk,
 const ATOOLS::Flavour &mofl,Double_Map &kt2,const SizeT_Map &cid,
 int cut,const size_t &step)
{
  Double_Map::const_iterator iit(kt2.find(idi));
  if (iit!=kt2.end()) {
    std::map<size_t,std::map<size_t,std::map<Flavour,CParam> > >::const_iterator 
      jit(iit->second.find(idj));
    if (jit!=iit->second.end()) {
      std::map<size_t,std::map<Flavour,CParam> >::const_iterator 
	kit(jit->second.find(idk));
      if (kit!=jit->second.end()) {
	std::map<Flavour,CParam>::const_iterator fit(kit->second.find(mofl));
	if (fit!=kit->second.end()) return fit->second;
      }
    }
  }
  int i(cid.find(idi)->second), j(cid.find(idj)->second);
  int k(cid.find(idk)->second);
  if (p_ampl->Leg(i)->Id()!=idi || p_ampl->Leg(j)->Id()!=idj || 
      p_ampl->Leg(k)->Id()!=idk) THROW(fatal_error,"Internal error 1");
  bool ismo(idi&((1<<p_xs->NIn())-1));
  Flavour mmofl(p_xs->ReMap(ismo?mofl.Bar():mofl,0));
  if (ismo) mmofl=mmofl.Bar();
  if (p_ampl->Legs().size()>p_ampl->NIn()+m_nmin) {
    int nlo((m_wmode&4096) && p_ampl->Prev()==NULL);
    kt2[idi][idj][idk][mofl]=
      p_clus->KPerp2(*p_ampl,i,j,k,mmofl,p_ms,
		     (m_wmode&1024)||nlo?1:-1,
		     (p_xs->Parent()->Info().m_fi.
		      m_nloqcdtype!=PHASIC::nlo_type::lo?16:0)|
		     (p_xs!=p_proc && step==2?32:0)|
		     ((cut||!mmofl.Strong())?1:0)|(nlo?4:0));
  }
  else {
    p_ampl->SetProc(p_xs);
    kt2[idi][idj][idk][mofl]=
      (p_xs->IsMapped()?p_xs->MapProc():p_xs)
      ->ScaleSetter()->CoreScale(p_ampl);
  }
  msg_Debugging()<<"calc Q_{"<<ID(idi)<<p_ampl->Leg(i)->Flav()
		 <<","<<ID(idj)<<""<<p_ampl->Leg(j)->Flav()
		 <<"->"<<mmofl<<";"
		 <<ID(idk)<<"} -> "<<kt2[idi][idj][idk][mofl]<<"\n";
  msg_Debugging()<<"  p_{"<<ID(idi)<<"} = "<<p_ampl->Leg(i)->Mom()
		 <<" "<<p_ampl->Leg(i)->Col()
		 <<" "<<p_ampl->Leg(i)->Flav()<<"\n";
  msg_Debugging()<<"  p_{"<<ID(idj)<<"} = "<<p_ampl->Leg(j)->Mom()
		 <<" "<<p_ampl->Leg(j)->Col()
		 <<" "<<p_ampl->Leg(j)->Flav()<<"\n";
  msg_Debugging()<<"  p_{"<<ID(idk)<<"} = "<<p_ampl->Leg(k)->Mom()
		 <<" "<<p_ampl->Leg(k)->Col()
		 <<" "<<p_ampl->Leg(k)->Flav()<<"\n";
  return kt2[idi][idj][idk][mofl];
}

void Cluster_Algorithm::CalculateMeasures
(const size_t &step,const Vertex_Set &nocl,
 const Current_Vector &ccurs,Current *const fcur,
 ClusterInfo_Map &cinfo,Double_Map &kt2,const SizeT_Map &cid)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  msg_Debugging()<<"p_ampl = "<<*p_ampl<<"\n";
  ClusterInfo_Map ccinfo(cinfo);
  cinfo.clear();
  for (size_t nc(2);nc<=step;++nc) {
    const Current_Vector &curs(p_bg->Currents()[nc]);
    for (size_t i(0);i<curs.size();++i) {
      const Vertex_Vector &in(curs[i]->In()); 
      for (size_t j(0);j<in.size();++j) {
	if (in[j]->Zero() || in[j]->J().size()>2) continue;
	if (in[j]->JC()->Sub()) continue;
	if (in[j]->JC()->Flav().IsDummy()) continue;
	if (in[j]->Order()[1]>p_ampl->OrderEW() ||
	    in[j]->Order()[0]>p_ampl->OrderQCD()) continue;
	if (find(ccurs.begin(),ccurs.end(),in[j]->J(0))==ccurs.end()) continue;
	if (find(ccurs.begin(),ccurs.end(),in[j]->J(1))==ccurs.end()) continue;
	size_t idi(in[j]->J(0)->CId()), idj(in[j]->J(1)->CId());
	if (p_xs!=p_proc && step==2) {
	  NLO_subevt *sub(p_proc->Get<Single_Dipole_Term>()->Sub());
	  if (!(m_id[idi]==(1<<sub->m_i) && m_id[idj]==(1<<sub->m_j)) &&
	      !(m_id[idj]==(1<<sub->m_i) && m_id[idi]==(1<<sub->m_j))) continue;
	}
	msg_Debugging()<<ID(m_id[idi])<<"&"<<ID(m_id[idj])<<": "
		       <<in[j]->J(0)->Flav()<<","<<in[j]->J(1)->Flav()
		       <<" -> "<<in[j]->JC()->Flav()<<" "
		       <<in[j]->Order()<<" {\n";
	{
	msg_Indent();
	const ColorID_Vector &coli(m_cols.find(m_id[idi])->second);
	const ColorID_Vector &colj(m_cols.find(m_id[idj])->second);
	for (size_t k(0);k<p_ampl->Legs().size();++k) {
	  size_t idk(p_ampl->Leg(k)->Id());
	  if (idk==m_id[idi] || idk==m_id[idj]) continue;
	  if (p_xs!=p_proc && step==2) {
	    NLO_subevt *sub(p_proc->Get<Single_Dipole_Term>()->Sub());
	    if (idk!=(1<<sub->m_k)) continue;
	  }
	  if (nocl.find(Cluster_Info(in[j],idk))!=nocl.end()) continue;
	  const ColorID_Vector &colk(m_cols.find(idk)->second);
	  for (int l(0);l<=1;++l) {
	    ColorID_Vector colij=Connected
	      (l?coli:colj,l?colj:coli,colk,in[j]->JC()->Flav());
	    if (colij.size()) {
	      CParam ckt2(GetMeasure(m_id[l?idi:idj],m_id[l?idj:idi],idk,
				     in[j]->JC()->Flav(),kt2,cid,
				     in[j]->JC()->Cut(),step));
	      cinfo.insert(ClusterInfo_Pair
			   (Cluster_Key(l?idi:idj,l?idj:idi),
			    Cluster_Info(in[j],idk,ckt2,in[j]->Order(1),
					 in[j]->Order(0),
					 in[j]->JC()->Flav(),colij)));
	    }
	  }
	}
	}
	msg_Debugging()<<"}\n";
      }
    }
  }
  const Vertex_Vector &in(fcur->In()); 
  for (size_t j(0);j<in.size();++j) {
    if (in[j]->Zero() || in[j]->J().size()>2) continue;
    for (size_t i(1);i<ccurs.size();++i) {
      if (in[j]->J(0)==ccurs[i] || in[j]->J(1)==ccurs[i]) {
	if (ccurs[i]->CId()&2) continue;
	Current *mocur(in[j]->J(0)==ccurs[i]?in[j]->J(1):in[j]->J(0));
	if (mocur->Sub()) continue;
	Flavour mofl(mocur->Flav().Bar());
	if (mofl.IsDummy()) continue;
	if (in[j]->Order()[1]>p_ampl->OrderEW() ||
	    in[j]->Order()[0]>p_ampl->OrderQCD()) continue;
	size_t idi(fcur->CId()), idj(ccurs[i]->CId());
	if (p_xs!=p_proc && step==2) {
	  NLO_subevt *sub(p_proc->Get<Single_Dipole_Term>()->Sub());
	  if (!(m_id[idi]==(1<<sub->m_i) && m_id[idj]==(1<<sub->m_j)) &&
	      !(m_id[idj]==(1<<sub->m_i) && m_id[idi]==(1<<sub->m_j))) continue;
	}
	msg_Debugging()<<ID(m_id[idi])<<"&"<<ID(m_id[idj])<<": "
		       <<fcur->Flav()<<","<<ccurs[i]->Flav()<<" -> "
		       <<mocur->Flav()<<" "<<in[j]->Order()<<" {\n";
	{
	  msg_Indent();
	  ColorID_Vector coli(m_cols.find(m_id[idi])->second);
	  ColorID_Vector colj(m_cols.find(m_id[idj])->second);
	  for (size_t k(0);k<p_ampl->Legs().size();++k) {
	    size_t idk(p_ampl->Leg(k)->Id());
	    if (idk==m_id[idi] || idk==m_id[idj]) continue;
	    if (p_xs!=p_proc && step==2) {
	      NLO_subevt *sub(p_proc->Get<Single_Dipole_Term>()->Sub());
	      if (idk!=(1<<sub->m_k)) continue;
	    }
	    if (nocl.find(Cluster_Info(in[j],idk))!=nocl.end()) continue;
	    ColorID_Vector colk(m_cols.find(idk)->second);
	    for (int l(0);l<=1;++l) {
	      ColorID_Vector colij=
		Connected(l?coli:colj,l?colj:coli,colk,mofl);
	      if (colij.size()) {
		CParam ckt2(GetMeasure(m_id[l?idi:idj],m_id[l?idj:idi],
				       idk,mofl,kt2,cid,0,step));
		cinfo.insert(ClusterInfo_Pair
			     (Cluster_Key(l?idi:idj,l?idj:idi),
			      Cluster_Info(in[j],idk,ckt2,in[j]->Order(1),
					   in[j]->Order(0),mofl,colij)));
	      }
	    }
	  }
	}
	msg_Debugging()<<"}\n";
      }
    }
  }
  msg_Debugging()<<"}\n";
}

bool Cluster_Algorithm::CombineWinner
(const Cluster_Info &ci,Current_Vector &ccurs,
 Current *&fcur,ClusterInfo_Map &cinfo)
{
  Vertex *v(ci.p_v);
  if (v->JC()!=fcur) {
    Current *ja(v->J(0)), *jb(v->J(1));
    m_id[v->JC()->CId()]=m_id[ja->CId()]+m_id[jb->CId()];
    if (v->J(0)->Id().front()>v->J(1)->Id().front()) 
      std::swap<Current*>(ja,jb);
    int found(0);
    for (Current_Vector::iterator 
	   cit(ccurs.begin());cit!=ccurs.end();++cit)
      if (*cit==ja) {
	*cit=v->JC();
	found+=1;
	break;
      }
    for (Current_Vector::iterator 
	   cit(ccurs.begin());cit!=ccurs.end();++cit)
      if (*cit==jb) {
	ccurs.erase(cit);
	found+=2;
	break;
      }
    if (found!=3) THROW(fatal_error,"Invalid clustering");
    msg_Debugging()<<"combine "<<ID(m_id[v->J(0)->CId()])
		   <<"&"<<ID(m_id[v->J(1)->CId()])<<" -> "
		   <<ID(m_id[v->JC()->CId()])<<" <-> "<<ID(ci.m_k)<<"\n";
  }
  else {
    bool found(false);
    for (Current_Vector::iterator 
	   cit(ccurs.begin());cit!=ccurs.end();++cit) 
      if (*cit==v->JC()) {
	Current_Vector::iterator fit(ccurs.begin());
	for (;fit!=ccurs.end();++fit) {
	  if (*fit==v->J(0)) {
	    m_id[v->J(1)->CId()]=m_id[v->JC()->CId()]+m_id[v->J(0)->CId()];
	    fcur=*cit=v->J(1);
	    found=true;
	    break;
	  }
	  if (*fit==v->J(1)) {
	    m_id[v->J(0)->CId()]=m_id[v->JC()->CId()]+m_id[v->J(1)->CId()];
	    fcur=*cit=v->J(0);
	    found=true;
	    break;
	  }
	}
	ccurs.erase(fit);
	break;
      }
    if (!found) THROW(fatal_error,"Invalid clustering");
    msg_Debugging()<<"combine "<<ID(m_id[v->JC()->CId()])
		   <<" -> "<<ID(m_id[v->J(0)->CId()])<<"&"
		   <<ID(m_id[v->J(1)->CId()])<<" <-> "<<ID(ci.m_k)<<"\n";
  }
  return true;
}

bool Cluster_Algorithm::ClusterStep
(const size_t &step,Vertex_Set &nocl,
 Current_Vector &ccurs,Current *&fcur,
 ClusterInfo_Map &cinfo,Double_Map &kt2)
{
  msg_Debugging()<<METHOD<<"(): step = "<<step<<" {\n";
  msg_Indent();
  SizeT_Map cid;
  for (size_t i(0);i<p_ampl->Legs().size();++i) 
    cid[p_ampl->Leg(i)->Id()]=i;
  CalculateMeasures(step,nocl,ccurs,fcur,cinfo,kt2,cid);
  if (cinfo.empty()) {
    msg_Debugging()<<"rejected configuration\n";
    return false;
  }
  double wmin(std::numeric_limits<double>::max());
  double rwmin(sqrt(std::numeric_limits<double>::max())), sum(0.0);
  ClusterInfo_Map::const_iterator win(cinfo.end()), rwin(win);
  for (ClusterInfo_Map::const_iterator cit(cinfo.begin());
       cit!=cinfo.end();++cit) {
    if (cit->second.m_mofl.IsDummy()) continue;
    if (cit->second.m_kt2.m_mode<0) continue;
    if (m_wmode&1) {
      if (cit->second.m_kt2.m_op2>=0.0 &&
	  cit->second.m_kt2.m_op2<wmin) {
	win=cit;
	wmin=cit->second.m_kt2.m_op2;
      }
      else if (cit->second.m_kt2.m_kt2>=0.0 &&
	       cit->second.m_kt2.m_kt2<rwmin) {
	rwin=cit;
	rwmin=cit->second.m_kt2.m_kt2;
      }
    }
    else {
      if (cit->second.m_kt2.m_op2>=0.0) {
	sum+=1.0/cit->second.m_kt2.m_op2;
      }
      else if (cit->second.m_kt2.m_kt2>=0.0 &&
	       cit->second.m_kt2.m_kt2<rwmin) {
	rwin=cit;
	rwmin=cit->second.m_kt2.m_kt2;
      }
    }
  }
  if (!(m_wmode&1)) {
    double disc(sum*ran->Get()), psum(0.0);
    for (ClusterInfo_Map::const_iterator cit(cinfo.begin());
	 cit!=cinfo.end();++cit) {
      if (cit->second.m_mofl.IsDummy()) continue;
      if (cit->second.m_kt2.m_op2>=0.0 &&
	  (psum+=1.0/cit->second.m_kt2.m_op2)>=disc) {
	win=cit;
	break;
      }
    }
    if (sum>0.0 && win==cinfo.end()) THROW(fatal_error,"Internal error 2"); 
  }
  if (win==cinfo.end() && !(m_wmode&512) &&
      !((m_wmode&4096) && step==2)) win=rwin;
  if (win==cinfo.end()) return false;
  Cluster_Key wkey(win->first);
  Cluster_Info winfo(win->second);
  if (p_xs==p_proc || step>m_nmin)
    nocl[winfo]=win->second.m_kt2.m_kt2-p_ampl->KT2();
  if (!CombineWinner(winfo,ccurs,fcur,cinfo)) return false;
  if (p_ampl->Legs().size()==p_ampl->NIn()+m_nmin) {
    bool match(false);
    if (p_ampl->NIn()==1) {
      if (ccurs.size()!=2) THROW(fatal_error,"Internal error 3");
      match=true;
    }
    else {
    if (ccurs.size()!=3) THROW(fatal_error,"Internal error 4");
    const Vertex_Vector &in(fcur->In());
    for (size_t i(0);i<in.size();++i) {
      size_t ncm(0);
      for (size_t j(0);j<ccurs.size();++j) {
	if (ccurs[j]==in[i]->J(0) ||
	    ccurs[j]==in[i]->J(1)) ++ncm;
      }
      if (ncm==2) {
	match=true;
	break;
      }
    }
    }
    if (!match) {
      msg_Debugging()<<"Invalid core\n";
      return false;
    }
  }
  Vec4D_Vector p;
  if (p_ampl->Legs().size()>p_ampl->NIn()+m_nmin) {
    p=p_clus->Combine(*p_ampl,cid[m_id[wkey.first]],
		      cid[m_id[wkey.second]],cid[winfo.m_k],
		      winfo.m_mofl,p_ms,winfo.m_kt2.m_kin,
		      winfo.m_kt2.m_mode);
    if (p.empty()) {
      msg_Debugging()<<"kinematics failed\n";
      return false;
    }
    if ((-p[0][0]>rpa->gen.PBeam(0)[0] &&
	 !IsEqual(-p[0][0],rpa->gen.PBeam(0)[0],1.0e-6)) ||
	(-p[1][0]>rpa->gen.PBeam(1)[0] &&
	 !IsEqual(-p[1][0]>rpa->gen.PBeam(1)[0],1.0e-6))) {
      msg_Debugging()<<"kinematics failed\n";
      return false;
    }
  }
  else if (p_ampl->Legs().size()==p_ampl->NIn()+m_nmin) {
    p.push_back(p_ampl->Leg(0)->Mom());
    if (p_ampl->NIn()==1) {
      p.push_back(p_ampl->Leg(0)->Mom());
    }
    else {
    p.push_back(p_ampl->Leg(1)->Mom());
    p.push_back(p_ampl->Leg(2)->Mom()+p_ampl->Leg(3)->Mom());
    }
  }
  else {
    THROW(fatal_error,"Invalid amplitude");
  }
  Cluster_Amplitude *ampl(p_ampl);
  ampl->SetKT2(winfo.m_kt2.m_kt2);
  ampl->SetMu2(winfo.m_kt2.m_mu2);
  p_ampl=p_ampl->InitNext();
  p_ampl->SetMS(p_ms);
  p_ampl->SetNIn(ampl->NIn());
  p_ampl->SetMuQ2(ampl->MuQ2());
  p_ampl->SetMuR2(ampl->MuR2());
  p_ampl->SetMuF2(ampl->MuF2());
  p_ampl->SetKT2(winfo.m_kt2.m_kt2);
  p_ampl->SetMu2(winfo.m_kt2.m_mu2);
  p_ampl->SetJF(ampl->JF<Selector_Base>());
  p_ampl->SetOrderEW(ampl->OrderEW()-winfo.p_v->Order(1));
  p_ampl->SetOrderQCD(ampl->OrderQCD()-winfo.p_v->Order(0));
  p_ampl->SetKin(winfo.m_kt2.m_kin);
  p_ampl->SetProcs(ampl->Procs<void>());
  p_ampl->Decays()=ampl->Decays();
  for (size_t i(0);i<ccurs.size();++i) {
    size_t cid(m_id[ccurs[i]->CId()]);
    Flavour flav(p_xs->ReMap(ccurs[i]->Flav(),0));
    if (ccurs[i]==fcur) flav=flav.Bar();
    ColorID col;
    for (size_t j(0);j<ampl->Legs().size();++j) {
      const Cluster_Leg *cli(ampl->Leg(j));
      if (cli->Id()&cid) {
	if (cli->Id()==cid) {
	  col=cli->Col();
	  break;
	}
      }
    }
    p_ampl->CreateLeg(p[i],flav,col,cid);
    p_ampl->Legs().back()->SetStat(1);
    if (col.m_i<0) {
      if (ccurs[i]->Cut()) {
	p_ampl->Legs().back()->SetStat
	  (p_ampl->Legs().back()->Stat()|2);
	SetNMax(p_ampl->Prev(),cid,ccurs[i]->Cut());
      }
      p_ampl->Legs().back()->SetCol(winfo.m_cols.front());
      p_ampl->Legs().back()->SetK(winfo.m_k);
      if (winfo.m_kt2.m_mode)
	p_ampl->Legs().back()->SetStat
	  (p_ampl->Legs().back()->Stat()|4);
    }
  }
  if (p_ampl->Legs().size()>3)
    m_cols[m_id[wkey.first]|m_id[wkey.second]]=winfo.m_cols;
  else  {
    const ColorID_Vector *ci(NULL), *cj(NULL);
    for (size_t i(0);i<p_ampl->Legs().size();++i)
      if (p_ampl->Leg(i)->Id()!=
  	  (m_id[wkey.first]|m_id[wkey.second])) {
  	if (ci==NULL) ci=&m_cols[p_ampl->Leg(i)->Id()];
  	else if (cj==NULL) cj=&m_cols[p_ampl->Leg(i)->Id()];
  	else THROW(fatal_error,"Internal error 5");
      }
    if (ci==NULL || cj==NULL) THROW(fatal_error,"Internal error 6");
    for (size_t i(0);i<ci->size();++i)
      for (size_t j(0);j<cj->size();++j) {
  	if (!winfo.m_mofl.Strong() && 
  	    (*ci)[i].m_i==(*cj)[j].m_j &&
  	    (*ci)[i].m_j==(*cj)[j].m_i) {
  	  return true;
  	}
  	ColorID_Vector cij(Connected((*ci)[i],(*cj)[j]));
	for (size_t l(0);l<cij.size();++l) {
  	  for (size_t k(0);k<winfo.m_cols.size();++k)
  	    if (winfo.m_cols[k].m_i==cij[l].m_j &&
  		winfo.m_cols[k].m_j==cij[l].m_i) {
  	      return true;
  	    }
	}
      }
    p_ampl=p_ampl->Prev();
    p_ampl->DeleteNext();
    msg_Debugging()<<"core color check failed\n";
    return false;
  }
  msg_Debugging()<<"} step = "<<step<<"\n";
  return true;
}

void Cluster_Algorithm::PreCluster
(Single_Process *const xs,Single_Dipole_Term *const dip,
 const Vec4D_Vector &p)
{
  if (p_clus==NULL) return;
  DEBUG_FUNC("");
  if (xs==NULL) THROW(fatal_error,"Internal error 7");
  p_proc=xs;
  p_bg=(p_xs=xs)->GetAmplitude();
  if (p_bg==NULL) THROW(fatal_error,"Internal error 8");
  p_bg->Differential();
}

bool Cluster_Algorithm::Cluster
(Single_Process *const xs,Single_Dipole_Term *const dip,
 const Vec4D_Vector &ip,const size_t &mode)
{
  m_wmode=mode;
  m_cols.clear();
  Vec4D_Vector p(ip);
  if (xs) {
    p_proc=xs;
    p_bg=(p_xs=xs)->GetAmplitude();
  }
  if (dip) {
    p_proc=dip;
    p_xs=dip->Process()->Get<COMIX::Single_Process>();
    p_bg=p_xs->GetAmplitude();
    p=p_bg->Momenta();
    for (size_t i(0);i<p_proc->NIn();++i) p[i]=-p[i];
  }
  if (p_bg==NULL) THROW(fatal_error,"Internal error 9");
  Selector_Base *jf=p_xs->Selector()->GetSelector("Jetfinder");
  msg_Debugging()<<METHOD<<"(mode = "<<mode<<"): {\n";
  m_nmin=Min((size_t)2,p_proc->Info().m_fi.NMinExternal());
  msg_Indent();
  m_id.clear();
  p_bg->ResetZero();
  Current_Vector ccurs(p_bg->Currents()[1]);
  Current *fcur(p_bg->Currents().back().front());
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(p_ms);
  p_ampl->SetJF(jf);
  p_ampl->SetNIn(p_proc->NIn());
  p_ampl->SetOrderEW(p_bg->MaxCpl()[1]/2);
  p_ampl->SetOrderQCD(p_bg->MaxCpl()[0]/2);
  p_ampl->SetProcs(p_xs->AllProcs());
  p_ampl->Decays()=p_bg->DecayInfos();
  PHASIC::Process_Base *pb(p_proc->IsMapped()?p_proc->MapProc():p_proc);
  double muf2(pb->ScaleSetter()->Scale(stp::fac));
  double mur2(pb->ScaleSetter()->Scale(stp::ren));
  double muq2(pb->ScaleSetter()->Scale(stp::res));
  NLO_subevt *sub
    (p_proc==p_xs?NULL:p_proc->
     Get<COMIX::Single_Dipole_Term>()->Sub());
  for (size_t i(0);i<ccurs.size();++i) {
    size_t cid(m_id[(i==0?fcur:ccurs[i])->CId()]=1<<p_ampl->Legs().size());
    Flavour flav(p_xs->ReMap(ccurs[i]->Flav(),0));
    Vec4D mom(i<2?-p[i]:p[i]);
    m_cols[cid]=ColorID_Vector(1,GetColor(ccurs[i],fcur));
    p_ampl->CreateLeg(mom,flav,m_cols[cid].front(),cid);
  }
  ccurs[0]=fcur;
  p_ampl->SetMuQ2(muq2);
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  ClusterInfo_Map cinfo;
  msg_Debugging()<<"}\n";
  if (p_clus==NULL) return true;
  if (m_nocluster) {
    p_ampl->SetProc(p_xs);
    p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		   ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
  }
  else {
  KT2Info_Vector kt2ord
    (1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
  const DecayInfo_Vector &decids(p_bg->DecayInfos());
  for (size_t i(0);i<decids.size();++i)
    kt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
  if (!Cluster(2,Vertex_Set(),ccurs,fcur,cinfo,kt2ord,
  	       (m_wmode&512)?1:((m_wmode&16384)?1:0))) {
    if (!(m_wmode&512)) {
    KT2Info_Vector kt2ord
      (1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
    const DecayInfo_Vector &decids(p_bg->DecayInfos());
    for (size_t i(0);i<decids.size();++i)
      kt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
    msg_Debugging()<<"trying all unordered configurations\n";
    if (!Cluster(2,Vertex_Set(),ccurs,fcur,cinfo,kt2ord,0))
      msg_Debugging()<<"no valid combination -> classify as core\n";
      p_ampl->SetProc(p_xs);
      p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		     ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
    }
  }
  }
  size_t nmax(p_proc->Info().m_fi.NMaxExternal());
  SetNMax(p_ampl,(1<<ccurs.size())-1,nmax);
  msg_Debugging()<<"Final configuration:\n";
  msg_Debugging()<<*p_ampl<<"\n";
  while (p_ampl->Prev()) {
    p_ampl=p_ampl->Prev();
    msg_Debugging()<<*p_ampl<<"\n";
  }
  return true;
}

KT2Info_Vector Cluster_Algorithm::UpdateKT2
(const KT2Info_Vector &kt2ord,const Cluster_Amplitude *ampl,
 const int mode) const
{
  KT2Info_Vector nkt2ord(kt2ord);
  Cluster_Leg *split(ampl->Next()->Splitter());
  size_t sid(split->Id()), lmin(100), li(0);
  for (size_t i(0);i<nkt2ord.size();++i) {
    if ((nkt2ord[i].first&sid)==sid &&
	IdCount(nkt2ord[i].first)<lmin) {
      lmin=IdCount(nkt2ord[i].first);
      li=i;
    }
  }
  if ((split->Stat()!=3 &&
       split->Flav().Strong()) ||
      ampl->Legs().size()==ampl->NIn()+m_nmin) {
    nkt2ord[li].second=(mode?ampl->Next():ampl)->KT2();
    msg_Debugging()<<"set last k_T = "<<sqrt(nkt2ord[li].second)
		   <<" for "<<ID(nkt2ord[li].first)
		   <<" from "<<ID(sid)<<"\n";
  }
  return nkt2ord;
}

bool Cluster_Algorithm::Cluster
(const size_t &step,const Vertex_Set &onocl,const Current_Vector &ccurs,
 Current *const fcur,const ClusterInfo_Map &cinfo,KT2Info_Vector &kt2ord,
 const int complete)
{
  if (p_ampl->Legs().size()==p_ampl->NIn()+1) {
    if (m_nmin==1 && CheckCore(p_ampl)) {
      p_ampl->SetProc(p_xs);
      p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		     ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
      if (p_ampl->Prev()) kt2ord=UpdateKT2(kt2ord,p_ampl->Prev(),1);
      return true;
    }
    if (p_ampl->Prev()) {
      p_ampl=p_ampl->Prev();
      p_ampl->DeleteNext();
    }
    return true;
  }
  for (int order(1);order>=0;--order) {
  DEBUG_FUNC("step = "<<step<<", order = "<<order);
  size_t oldsize(0);
  Double_Map kt2;
  Vertex_Set nocl;
  Cluster_Amplitude *ampl(p_ampl);
  do {
    oldsize=nocl.size();
    Current_Vector nccurs(ccurs);
    Current *nfcur(fcur);
    ClusterInfo_Map ncinfo(cinfo);
    if (ClusterStep(step,nocl,nccurs,nfcur,ncinfo,kt2)) {
      KT2Info_Vector nkt2ord(UpdateKT2(kt2ord,ampl)), nnkt2ord(nkt2ord);
      if (!order || ((m_wmode&4096) && step==2)) {
	nkt2ord=KT2Info_Vector(1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
	const DecayInfo_Vector &decids(p_bg->DecayInfos());
	for (size_t i(0);i<decids.size();++i)
	  nkt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
	nnkt2ord=nkt2ord;
      }
      if (Cluster(step+1,nocl,nccurs,nfcur,ncinfo,nnkt2ord,complete)) {
  	if (ampl->Legs().size()==ampl->NIn()+m_nmin) {
	  kt2ord=nkt2ord;
	  return true;
	}
	bool ord(true);
	msg_Debugging()<<"check ordering:\n";
	for (size_t i(0);i<nkt2ord.size();++i) {
	  msg_Debugging()<<"  "<<ID(nkt2ord[i].first)<<": "
			 <<sqrt(nnkt2ord[i].second)<<" vs. "
			 <<sqrt(nkt2ord[i].second)<<"\n";
	  if (nnkt2ord[i].second<nkt2ord[i].second) {
	    msg_Debugging()<<"unordered configuration\n";
	    ord=false;
	    break;
	  }
	}
	if (ord || (m_wmode&16)) {
	  kt2ord=nkt2ord;
	  return true;
	}
	msg_Debugging()<<"reject ordering\n";
      }
      else if ((m_wmode&512) && (m_wmode&4096) && step==2) {
	if (!CheckCore(p_ampl)) continue;
	msg_Debugging()<<"no valid combination -> classify as RS core\n";
	p_ampl->SetProc(p_xs);
	p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		       ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
	return true;
      }
      p_ampl=ampl;
      p_ampl->DeleteNext();
    }
    else if (nocl.empty()) {
      if (!CheckCore(p_ampl)) continue;
      msg_Debugging()<<"no valid combination -> classify as core\n";
      p_ampl->SetProc(p_xs);
      p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		     ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
      if (p_ampl->Prev()) kt2ord=UpdateKT2(kt2ord,p_ampl->Prev(),1);
      return true;
    }
  } while (oldsize<nocl.size());
  if (complete==1) return false;
  if (m_wmode&512) {
    if (!CheckCore(p_ampl)) return false;
    msg_Debugging()<<"no valid combination -> classify as core\n";
    p_ampl->SetProc(p_xs);
    p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		   ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
    if (p_ampl->Prev()) kt2ord=UpdateKT2(kt2ord,p_ampl->Prev(),1);
    return true;
  }
  if (order) msg_Debugging()<<"trying unordered configurations\n";
  }
  return false;
}

bool Cluster_Algorithm::CheckCore(ATOOLS::Cluster_Amplitude *const ampl) const
{
  if (!m_corecheck) return true;
  PHASIC::Process_Base::SortFlavours(ampl);
  std::string name(PHASIC::Process_Base::GenerateName(ampl));
  StringProcess_Map *pm((*p_xs->AllProcs())[nlo_type::lo]);
  StringProcess_Map::const_iterator pit(pm->find(name));
  if (pit==pm->end()) msg_Debugging()<<"invalid core configuration\n";
  return pit!=pm->end();
}

void Cluster_Algorithm::SetNMax(Cluster_Amplitude *const ampl,
				const size_t &id,const size_t &nmax) const
{
  if (ampl==NULL) return;
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cli(ampl->Leg(i));
    if (cli->Id()&id) {
      cli->SetNMax(nmax);
      if (cli->Stat()!=3) 
	SetNMax(ampl->Prev(),cli->Id(),nmax);
    }
  }
}
