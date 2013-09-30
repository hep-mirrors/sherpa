#include "COMIX/Cluster/Cluster_Algorithm.H"

#include "COMIX/Cluster/Color_Setter.H"
#include "COMIX/Main/Single_Process.H"
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
  m_cs(this), p_ms(ms), p_ampl(NULL), p_clus(NULL),
  m_lfrac(0.0)
{
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

int Cluster_Algorithm::Connected
(const Vertex_Vector &dvs,const size_t &idk) const
{
  for (size_t i(0);i<dvs.size();++i) {
    SizeT_Map::const_iterator ait(m_id.find(dvs[i]->JA()->CId()));
    SizeT_Map::const_iterator bit(m_id.find(dvs[i]->JB()->CId()));
    SizeT_Map::const_iterator cit(m_id.find(dvs[i]->JC()->CId()));
    if ((ait!=m_id.end() && ait->second==idk) ||
	(bit!=m_id.end() && bit->second==idk) ||
	(cit!=m_id.end() && cit->second==idk)) {
      return dvs[i]->JC()->Zero()?0:1;
    }
  }
  return -1;
}

CParam Cluster_Algorithm::GetMeasure
(const size_t &idi,const size_t &idj,const size_t &idk,
 const ATOOLS::Flavour &mofl,Double_Map &kt2,const SizeT_Map &cid,int cut)
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
      p_ampl->Leg(k)->Id()!=idk) THROW(fatal_error,"Internal error");
  bool ismo(idi&((1<<p_xs->NIn())-1));
  Flavour mmofl(p_xs->ReMap(ismo?mofl.Bar():mofl,0));
  if (ismo) mmofl=mmofl.Bar();
  if (p_ampl->Legs().size()>p_ampl->NIn()+2) {
    kt2[idi][idj][idk][mofl]=
      p_clus->KPerp2(*p_ampl,i,j,k,mmofl,p_ms,
		     (m_wmode&1024)||(m_wmode&4096)?1:-1,
		     ((cut||!mmofl.Strong())?1:0)|
		     (p_xs->Parent()->Info().m_fi.m_nloqcdtype!=
		      PHASIC::nlo_type::lo?2:0));
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
		 <<" "<<p_ampl->Leg(i)->Col()<<"\n";
  msg_Debugging()<<"  p_{"<<ID(idj)<<"} = "<<p_ampl->Leg(j)->Mom()
		 <<" "<<p_ampl->Leg(j)->Col()<<"\n";
  msg_Debugging()<<"  p_{"<<ID(idk)<<"} = "<<p_ampl->Leg(k)->Mom()
		 <<" "<<p_ampl->Leg(k)->Col()<<"\n";
  return kt2[idi][idj][idk][mofl];
}

void Cluster_Algorithm::CalculateMeasures
(const size_t &step,const Vertex_Set &nocl,
 const Current_Vector &ccurs,Current *const fcur,
 ClusterInfo_Map &cinfo,Double_Map &kt2,const SizeT_Map &cid)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  msg_Debugging()<<*p_ampl<<"\n";
  ClusterInfo_Map ccinfo(cinfo);
  cinfo.clear();
  for (size_t nc(2);nc<=step;++nc) {
    const Current_Vector &curs(p_bg->Currents()[nc]);
    for (size_t i(0);i<curs.size();++i) {
      const Vertex_Vector &in(curs[i]->In()); 
      for (size_t j(0);j<in.size();++j) {
	if (in[j]->Zero()&&!m_nosol) continue;
	if (in[j]->JC()->Flav().IsDummy()) continue;
	if (find(ccurs.begin(),ccurs.end(),in[j]->JA())==ccurs.end()) continue;
	if (find(ccurs.begin(),ccurs.end(),in[j]->JB())==ccurs.end()) continue;
	size_t idi(in[j]->JA()->CId()), idj(in[j]->JB()->CId());
	msg_Debugging()<<ID(m_id[idi])<<"&"<<ID(m_id[idj])<<": "
		       <<in[j]->JA()->Flav()<<","<<in[j]->JB()->Flav()
		       <<" -> "<<in[j]->JC()->Flav()<<" ["
		       <<in[j]->OrderEW()<<","<<in[j]->OrderQCD()<<"] {\n";
	{
	msg_Indent();
	ColorID coli(p_ampl->Leg(cid.find(m_id[idi])->second)->Col());
	ColorID colj(p_ampl->Leg(cid.find(m_id[idj])->second)->Col());
	for (size_t k(0);k<p_ampl->Legs().size();++k) {
	  size_t idk(p_ampl->Leg(k)->Id());
	  if (idk==m_id[idi] || idk==m_id[idj]) continue;
	  if (nocl.find(Cluster_Info(in[j],idk))!=nocl.end()) continue;
	  ColorID colk(p_ampl->Leg(k)->Col());
	  int cc[2]={0,0};
	  for (int l(0);l<=1;++l) {
	    cc[l]=l?Connected(in[j]->JB()->Out(),idk):
	      Connected(in[j]->JA()->Out(),idk);
	    if (cc[l]<=0) cc[l]=Connected(in[j]->JC()->Out(),idk);
	    if (in[j]->JC()->Cut()) cc[l]=1;
	    if (p_ampl->Legs().size()==p_ampl->NIn()+2 || cc[l]>0) {
	      if (in[j]->OrderQCD()) {
		ColorID ccj(l?colj:coli);
		if (!((ccj.m_i && ccj.m_i==colk.m_j) ||
		      (ccj.m_j && ccj.m_j==colk.m_i))) continue;
	      }
	      CParam ckt2(GetMeasure(m_id[l?idi:idj],m_id[l?idj:idi],idk,
				     in[j]->JC()->Flav(),kt2,cid,
				     in[j]->JC()->Cut()));
	      cinfo.insert(ClusterInfo_Pair
			   (Cluster_Key(l?idi:idj,l?idj:idi),
			    Cluster_Info(in[j],idk,ckt2,in[j]->OrderEW(),
					 in[j]->OrderQCD(),in[j]->JC()->Flav())));
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
    if (in[j]->Zero()&&!m_nosol) continue;
    for (size_t i(1);i<ccurs.size();++i) {
      if (in[j]->JA()==ccurs[i] || in[j]->JB()==ccurs[i]) {
	if (ccurs[i]->CId()&2) continue;
	Current *mocur(in[j]->JA()==ccurs[i]?in[j]->JB():in[j]->JA());
	Flavour mofl(mocur->Flav().Bar());
	if (mofl.IsDummy()) continue;
	size_t idi(fcur->CId()), idj(ccurs[i]->CId());
	msg_Debugging()<<ID(m_id[idi])<<"&"<<ID(m_id[idj])<<": "
		       <<fcur->Flav()<<","<<ccurs[i]->Flav()<<" -> "
		       <<mocur->Flav()<<" ["<<in[j]->OrderEW()
		       <<","<<in[j]->OrderQCD()<<"] {\n";
	{
	  msg_Indent();
	  ColorID coli(p_ampl->Leg(cid.find(m_id[idi])->second)->Col());
	  ColorID colj(p_ampl->Leg(cid.find(m_id[idj])->second)->Col());
	  for (size_t k(0);k<p_ampl->Legs().size();++k) {
	    size_t idk(p_ampl->Leg(k)->Id());
	    if (idk==m_id[idi] || idk==m_id[idj]) continue;
	    if (nocl.find(Cluster_Info(in[j],idk))!=nocl.end()) continue;
	    ColorID colk(p_ampl->Leg(k)->Col());
	    int cc[2]={0,0};
	    for (int l(0);l<=1;++l) {
	      cc[l]=l?Connected(ccurs[i]->Out(),idk):
		Connected(fcur->In(),idk);
	      if (cc[l]<=0) cc[l]=Connected(mocur->In(),idk);
	      if (p_ampl->Legs().size()==p_ampl->NIn()+2 || cc[l]>0) {
		if (in[j]->OrderQCD()) {
		  ColorID ccj(l?colj:coli);
		  if (!((ccj.m_i && ccj.m_i==colk.m_j) ||
			(ccj.m_j && ccj.m_j==colk.m_i))) continue;
		}
		CParam ckt2(GetMeasure(m_id[l?idi:idj],m_id[l?idj:idi],
				       idk,mofl,kt2,cid,0));
		cinfo.insert(ClusterInfo_Pair
			     (Cluster_Key(l?idi:idj,l?idj:idi),
			      Cluster_Info(in[j],idk,ckt2,in[j]->OrderEW(),
					   in[j]->OrderQCD(),mofl)));
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
    Current *ja(v->JA()), *jb(v->JB());
    m_id[v->JC()->CId()]=m_id[ja->CId()]+m_id[jb->CId()];
    if (v->JA()->Id().front()>v->JB()->Id().front()) 
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
    msg_Debugging()<<"combine "<<ID(m_id[v->JA()->CId()])
		   <<"&"<<ID(m_id[v->JB()->CId()])<<" -> "
		   <<ID(m_id[v->JC()->CId()])<<" <-> "<<ID(ci.m_k)<<"\n";
  }
  else {
    bool found(false);
    for (Current_Vector::iterator 
	   cit(ccurs.begin());cit!=ccurs.end();++cit) 
      if (*cit==v->JC()) {
	Current_Vector::iterator fit(ccurs.begin());
	for (;fit!=ccurs.end();++fit) {
	  if (*fit==v->JA()) {
	    m_id[v->JB()->CId()]=m_id[v->JC()->CId()]+m_id[v->JA()->CId()];
	    fcur=*cit=v->JB();
	    found=true;
	    break;
	  }
	  if (*fit==v->JB()) {
	    m_id[v->JA()->CId()]=m_id[v->JC()->CId()]+m_id[v->JB()->CId()];
	    fcur=*cit=v->JA();
	    found=true;
	    break;
	  }
	}
	ccurs.erase(fit);
	break;
      }
    if (!found) THROW(fatal_error,"Invalid clustering");
    msg_Debugging()<<"combine "<<ID(m_id[v->JC()->CId()])
		   <<" -> "<<ID(m_id[v->JA()->CId()])<<"&"
		   <<ID(m_id[v->JB()->CId()])<<" <-> "<<ID(ci.m_k)<<"\n";
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
    if (sum>0.0 && win==cinfo.end()) THROW(fatal_error,"Internal error"); 
  }
  if (win==cinfo.end() && !(m_wmode&512)) win=rwin;
  if (win==cinfo.end()) return false;
  Cluster_Key wkey(win->first);
  Cluster_Info winfo(win->second);
  nocl[winfo]=win->second.m_kt2.m_kt2-p_ampl->KT2();
  if (!CombineWinner(winfo,ccurs,fcur,cinfo)) return false;
  if (p_ampl->Legs().size()==p_ampl->NIn()+2) {
    bool match(false);
    if (p_ampl->NIn()==1) {
      if (ccurs.size()!=2) THROW(fatal_error,"Internal error");
      match=true;
    }
    else {
    if (ccurs.size()!=3) THROW(fatal_error,"Internal error");
    const Vertex_Vector &in(fcur->In());
    for (size_t i(0);i<in.size();++i) {
      size_t ncm(0);
      for (size_t j(0);j<ccurs.size();++j) {
	if (ccurs[j]==in[i]->JA() ||
	    ccurs[j]==in[i]->JB()) ++ncm;
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
  if (p_ampl->Legs().size()>p_ampl->NIn()+2) {
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
  else if (p_ampl->Legs().size()==p_ampl->NIn()+2) {
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
  p_ampl->SetQ2(ampl->Q2());
  p_ampl->SetMuR2(ampl->MuR2());
  p_ampl->SetMuF2(ampl->MuF2());
  p_ampl->SetKT2(winfo.m_kt2.m_kt2);
  p_ampl->SetMu2(winfo.m_kt2.m_mu2);
  p_ampl->SetJF(ampl->JF<Selector_Base>());
  p_ampl->SetOrderEW(ampl->OrderEW()-winfo.p_v->OrderEW());
  p_ampl->SetOrderQCD(ampl->OrderQCD()-winfo.p_v->OrderQCD());
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
    if (IdCount(m_id[ccurs[i]->CId()])==1) {
      p_ampl->Legs().back()->SetStat(1);
    }
    else if (col.m_i<0 && winfo.m_kt2.m_mode) {
      size_t dmax(ccurs[i]->Cut()?ccurs[i]->Cut():IdCount(cid));
      p_ampl->Legs().back()->SetStat(3);
      SetNMax(p_ampl->Prev(),cid,dmax);
    }
    if (col.m_i<0) {
      p_ampl->Legs().back()->SetCol(GetColor(ccurs[i],fcur));
      p_ampl->Legs().back()->SetK(winfo.m_k);
    }
  }
  msg_Debugging()<<"} step = "<<step<<"\n";
  return true;
}

bool Cluster_Algorithm::Cluster
(Single_Process *const xs,const size_t &mode)
{
  m_wmode=mode;
  p_bg=(p_xs=xs)->GetAmplitude();
  if (p_bg==NULL) THROW(fatal_error,"Internal error");
  Selector_Base *jf=p_xs->Selector()
    ->GetSelector("Jetfinder");
  msg_Debugging()<<METHOD<<"(mode = "<<mode<<"): {\n";
  msg_Indent();
  m_id.clear();
  Current_Vector ccurs(p_bg->Currents()[1]);
  Current *fcur(ccurs[0]=p_bg->Currents().back().front());
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(p_ms);
  p_ampl->SetJF(jf);
  p_ampl->SetNIn(xs->NIn());
  p_ampl->SetOrderEW(p_bg->MaxOrderEW());
  p_ampl->SetOrderQCD(p_bg->MaxOrderQCD());
  p_ampl->SetProcs(p_xs->AllProcs());
  p_ampl->Decays()=p_bg->DecayInfos();
  PHASIC::Process_Base *pb(xs->Process()->IsMapped()?
			   xs->Process()->MapProc():xs);
  double muf2(pb->ScaleSetter()->Scale(stp::fac));
  double mur2(pb->ScaleSetter()->Scale(stp::ren));
  double Q2(pb->ScaleSetter()->Scale(stp::res));
  for (size_t i(0);i<ccurs.size();++i) {
    size_t cid(m_id[ccurs[i]->CId()]=1<<p_ampl->Legs().size());
    Flavour flav(p_xs->ReMap(ccurs[i]->Flav(),0));
    if (ccurs[i]==fcur) flav=flav.Bar();
    size_t idx(i);
    Vec4D mom(i<2?-xs->Process()->Integrator()->Momenta()[idx]:
	      xs->Process()->Integrator()->Momenta()[idx]);
    p_ampl->CreateLeg(mom,flav,ColorID(),cid);
    p_ampl->Legs().back()->SetStat(1);
  }
  p_ampl->SetQ2(Q2);
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  Cluster_Amplitude *eampl(p_ampl);
  m_nosol=!m_cs.SetColors(xs);
  ClusterInfo_Map cinfo;
  msg_Debugging()<<"}\n";
  KT2Info_Vector kt2ord
    (1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
  const DecayInfo_Vector &decids(p_bg->DecayInfos());
  for (size_t i(0);i<decids.size();++i)
    kt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
  if (!Cluster(2,Vertex_Set(),ccurs,fcur,cinfo,kt2ord,1)) {
    if (!(m_wmode&512)) {
    KT2Info_Vector kt2ord
      (1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
    const DecayInfo_Vector &decids(p_bg->DecayInfos());
    for (size_t i(0);i<decids.size();++i)
      kt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
    msg_Debugging()<<"trying all unordered configurations\n";
    if (!Cluster(2,Vertex_Set(),ccurs,fcur,cinfo,kt2ord,0))
      THROW(fatal_error,"Internal error");
    }
    else {
      msg_Debugging()<<"no valid combination -> classify as core\n";
      p_ampl->SetProc(p_xs);
      p_ampl->SetKT2((p_xs->IsMapped()?p_xs->MapProc():p_xs)
		     ->ScaleSetter()->CoreScale(p_ampl).m_mu2);
      if (p_ampl->Prev()) kt2ord=UpdateKT2(kt2ord,p_ampl->Prev(),1);
      return true;
    }
  }
  size_t nmax(xs->Process()->Info().m_fi.NMaxExternal());
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
      ampl->Legs().size()==ampl->NIn()+2) {
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
    p_ampl=p_ampl->Prev();
    p_ampl->DeleteNext();
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
      if (!order) {
	nkt2ord=KT2Info_Vector(1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
	const DecayInfo_Vector &decids(p_bg->DecayInfos());
	for (size_t i(0);i<decids.size();++i)
	  nkt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
	nnkt2ord=nkt2ord;
      }
      if (Cluster(step+1,nocl,nccurs,nfcur,ncinfo,nnkt2ord,complete)) {
  	if (ampl->Legs().size()==ampl->NIn()+2) {
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
      p_ampl=ampl;
      p_ampl->DeleteNext();
    }
    else if (nocl.empty()) {
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
