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
  m_cnt(0), m_rej(0), m_lfrac(0.0)
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

bool Cluster_Algorithm::EWConnected
(const ATOOLS::Flavour &c,const ATOOLS::Flavour &s) const
{
  if (c.IntCharge()==0) return s.IntCharge()!=0;
  return c.IntCharge()*s.IntCharge()<0;
}

int Cluster_Algorithm::
ColorConnected(const ColorID &i,const ColorID &j,const ColorID &k) const
{
  if (i.m_i && i.m_j && j.m_i && j.m_j) {
    if (k.m_j==i.m_i) return -1;
    if (k.m_i==i.m_j) return -1;
    if (k.m_j==j.m_i) return 1;
    if (k.m_i==j.m_j) return 1;
  }
  else {
    if (i.m_i>0 && k.m_j==i.m_i) return 2;
    if (i.m_j>0 && k.m_i==i.m_j) return 2;
    if (j.m_i>0 && k.m_j==j.m_i) return 2;
    if (j.m_j>0 && k.m_i==j.m_j) return 2;
  }
  if ((k.m_i>0)^(k.m_j>0)) {
    // coloured singlet
    if (i.m_i>0 && i.m_i==j.m_j && 
	i.m_j>0 && i.m_j==j.m_i) return 3;
    // colourless singlet
    if (i.m_i==0 && i.m_j==0 &&
	j.m_i==0 && j.m_j==0) return 3;
  }
  // all colourless
  if (k.m_i==0 && k.m_j==0 &&
      i.m_i==0 && i.m_j==0 &&
      j.m_i==0 && j.m_j==0) return 3;
  return 0;
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
  if (p_ampl->Legs().size()>4) {
    kt2[idi][idj][idk][mofl]=
      p_clus->KPerp2(*p_ampl,i,j,k,mmofl,p_ms,(m_wmode&1024)?1:-1,
		     (cut||!mmofl.Strong())?1:0);
  }
  else {
    p_ampl->SetProcs(p_xs);
    kt2[idi][idj][idk][mofl]=p_clus->CoreScale(p_ampl);
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
	  int cc(ColorConnected(coli,colj,colk));
	  if (cc==2 && (idi&3)<(idj&3)) cc=-1;
	  if (p_ampl->Legs().size()==4 ||
	      (in[j]->OrderQCD()==0?
	       EWConnected(in[j]->JC()->Flav(),p_ampl->Leg(k)->Flav()):cc)) {
	    if (cc==0) cc=1;
	    CParam ckt2(GetMeasure(m_id[cc>0?idi:idj],m_id[cc>0?idj:idi],idk,
				   in[j]->JC()->Flav(),kt2,cid,
				   in[j]->JC()->Cut()));
	    cinfo.insert(ClusterInfo_Pair
			 (Cluster_Key(cc>0?idi:idj,cc>0?idj:idi),
			  Cluster_Info(in[j],idk,ckt2,in[j]->OrderEW(),
				       in[j]->OrderQCD(),in[j]->JC()->Flav())));
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
	Flavour mofl((in[j]->JA()==ccurs[i]?
		      in[j]->JB():in[j]->JA())->Flav().Bar());
	if (mofl.IsDummy()) continue;
	size_t idi(fcur->CId()), idj(ccurs[i]->CId());
	msg_Debugging()<<ID(m_id[idi])<<"&"<<ID(m_id[idj])<<": "
		       <<fcur->Flav()<<","<<ccurs[i]->Flav()<<" -> "
		       <<(in[j]->JA()==ccurs[i]?in[j]->JB()->Flav():
			  in[j]->JA()->Flav())<<" ["
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
	    int cc(ColorConnected(coli,colj,colk));
	    if (p_ampl->Legs().size()==4 ||
		(in[j]->OrderQCD()==0?
		 EWConnected(mofl,p_ampl->Leg(k)->Flav()):cc)) {
	      if (cc==0) cc=1;
	      CParam ckt2(GetMeasure(m_id[cc>0?idi:idj],m_id[cc>0?idj:idi],
				     idk,mofl,kt2,cid,0));
	      cinfo.insert(ClusterInfo_Pair
			   (Cluster_Key(cc>0?idi:idj,cc>0?idj:idi),
			    Cluster_Info(in[j],idk,ckt2,in[j]->OrderEW(),
					 in[j]->OrderQCD(),mofl)));
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
  double rwmin(std::numeric_limits<double>::max()), sum(0.0);
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
  DecayInfo_Vector decays(p_ampl->Decays());
  const DecayInfo_Vector &decids(p_bg->DecayInfos());
  for (size_t j(0);j<decids.size();++j)
    if (decids[j]->m_id==m_id[wkey.first]+m_id[wkey.second]) {
      decays.push_back(decids[j]);
      break;
    }
  if (p_ampl->Legs().size()==4) {
    if (ccurs.size()!=3) THROW(fatal_error,"Internal error");
    bool match(false);
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
    if (decays.size()!=p_bg->DecayInfos().size()) {
      msg_Debugging()<<"Unclustered decay\n"<<*p_ampl<<"\n";
      match=false;
    }
    if (!match) {
      msg_Debugging()<<"Invalid core\n";
      return false;
    }
  }
  Vec4D_Vector p;
  if (p_ampl->Legs().size()>4) {
    p=p_clus->Combine(*p_ampl,cid[m_id[wkey.first]],
		      cid[m_id[wkey.second]],cid[winfo.m_k],
		      winfo.m_mofl,p_ms,winfo.m_kt2.m_kin,
		      winfo.m_kt2.m_mode);
    if (p.empty()) {
      msg_Debugging()<<"kinematics failed\n";
      return false;
    }
    if ((-p[m_swap][0]>rpa->gen.PBeam(0)[0] &&
	 !IsEqual(-p[m_swap][0],rpa->gen.PBeam(0)[0],1.0e-6)) ||
	(-p[1-m_swap][0]>rpa->gen.PBeam(1)[0] &&
	 !IsEqual(-p[1-m_swap][0]>rpa->gen.PBeam(1)[0],1.0e-6))) {
      msg_Debugging()<<"kinematics failed\n";
      return false;
    }
  }
  else if (p_ampl->Legs().size()==4) {
    p.push_back(p_ampl->Leg(0)->Mom());
    p.push_back(p_ampl->Leg(1)->Mom());
    p.push_back(p_ampl->Leg(2)->Mom()+p_ampl->Leg(3)->Mom());
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
  p_ampl->SetMuR2(ampl->MuR2());
  p_ampl->SetMuF2(ampl->MuF2());
  p_ampl->SetKT2(winfo.m_kt2.m_kt2);
  p_ampl->SetMu2(winfo.m_kt2.m_mu2);
  p_ampl->SetJF(ampl->JF<Selector_Base>());
  p_ampl->SetOrderEW(ampl->OrderEW()-winfo.p_v->OrderEW());
  p_ampl->SetOrderQCD(ampl->OrderQCD()-winfo.p_v->OrderQCD());
  p_ampl->SetKin(winfo.m_kt2.m_kin);
  p_ampl->Decays()=decays;
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
(Single_Process *const xs,const size_t &mode,const double &kt2)
{
  m_wmode=mode;
  p_bg=(p_xs=xs)->GetAmplitude();
  m_swap=xs->Process()->Integrator()->InSwaped();
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
  p_ampl->SetKT2(kt2);
  p_ampl->SetMu2(kt2);
  p_ampl->SetNIn(xs->NIn());
  p_ampl->SetOrderEW(p_bg->MaxOrderEW());
  p_ampl->SetOrderQCD(p_bg->MaxOrderQCD());
  PHASIC::Process_Base *pb(xs->Process()->IsMapped()?
			   xs->Process()->MapProc():xs);
  double muf2(pb->ScaleSetter()->Scale(stp::fac));
  double mur2(pb->ScaleSetter()->Scale(stp::ren));
  for (size_t i(0);i<ccurs.size();++i) {
    size_t cid(m_id[ccurs[i]->CId()]=1<<p_ampl->Legs().size());
    Flavour flav(p_xs->ReMap(ccurs[i]->Flav(),0));
    if (ccurs[i]==fcur) flav=flav.Bar();
    size_t idx(i<2?(m_swap?1-i:i):i);
    Vec4D mom(i<2?-xs->Process()->Integrator()->Momenta()[idx]:
	      xs->Process()->Integrator()->Momenta()[idx]);
    p_ampl->CreateLeg(mom,flav,ColorID(),cid);
    p_ampl->Legs().back()->SetStat(1);
  }
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  Cluster_Amplitude *eampl(p_ampl);
  if (m_swap) xs->Process()->SwapInOrder();
  m_nosol=!m_cs.SetColors(xs);
  if (m_swap) xs->Process()->SwapInOrder();
  ClusterInfo_Map cinfo;
  ++m_cnt;
  KT2Info_Vector kt2ord
    (1,KT2_Info((1<<p_ampl->Legs().size())-1,0.0));
  const DecayInfo_Vector &decids(p_bg->DecayInfos());
  for (size_t i(0);i<decids.size();++i)
    kt2ord.push_back(std::make_pair(decids[i]->m_id,0.0));
  if (!Cluster(2,Vertex_Set(),ccurs,fcur,cinfo,kt2ord)) {
    msg_Debugging()<<METHOD<<"(): No valid PS history.\n";
    ++m_rej;
    double frac(m_rej/(double)m_cnt);
    if (frac>1.25*m_lfrac && m_cnt>1000) {
      m_lfrac=frac;
      msg_Error()<<METHOD<<"(): No valid PS history in >"
		 <<(int(m_lfrac*1000)/10.0)<<"% of calls.\n";
    }
    p_ampl->Delete();
    p_ampl=NULL;
    return false;
  }
  msg_Debugging()<<"}\n";
  size_t nmax(xs->Process()->Info().m_fi.NMaxExternal());
  SetNMax(p_ampl,(1<<ccurs.size())-1,nmax);
  msg_Debugging()<<*p_ampl<<"\n";
  while (p_ampl->Prev()) {
    if (m_swap) std::swap<Cluster_Leg*>
      (p_ampl->Legs()[0],p_ampl->Legs()[1]);
    p_ampl=p_ampl->Prev();
    msg_Debugging()<<*p_ampl<<"\n";
  }
  if (m_swap) std::swap<Cluster_Leg*>
    (p_ampl->Legs()[0],p_ampl->Legs()[1]);
  return true;
}

bool Cluster_Algorithm::Cluster
(const size_t &step,const Vertex_Set &onocl,const Current_Vector &ccurs,
 Current *const fcur,const ClusterInfo_Map &cinfo,KT2Info_Vector &kt2ord)
{
  if (p_ampl->Legs().size()==3) {
    p_ampl=p_ampl->Prev();
    p_ampl->Decays()=p_ampl->Next()->Decays();
    p_ampl->DeleteNext();
    return true;
  }
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
      Cluster_Leg *split(ampl->Next()->Splitter());
      size_t sid(split->Id()), lmin(100), li(0);
      for (size_t i(0);i<kt2ord.size();++i) {
	if ((kt2ord[i].first&sid)==sid &&
	    IdCount(kt2ord[i].first)<lmin) {
	  lmin=IdCount(kt2ord[i].first);
	  li=i;
	}
      }
      if ((split->Stat()!=3 &&
	   split->Flav().Strong()) ||
	  p_ampl->Legs().size()==4)
	kt2ord[li].second=ampl->Next()->KT2();
      msg_Debugging()<<"set last k_T = "<<sqrt(ampl->Next()->KT2())
		     <<" "<<ID(kt2ord[li].first)<<"\n";
      KT2Info_Vector nkt2ord(kt2ord);
      if (Cluster(step+1,nocl,nccurs,nfcur,ncinfo,nkt2ord)) {
  	if (ampl->Legs().size()==4) return true;
	bool ord(true);
	for (size_t i(0);i<kt2ord.size();++i)
	  if (nkt2ord[i].second<kt2ord[i].second) {
	    msg_Debugging()<<"unordered configuration: "
			   <<sqrt(nkt2ord[i].second)<<" vs. "
			   <<sqrt(kt2ord[i].second)<<" "
			   <<ID(kt2ord[i].first)<<"\n";
	    ord=false;
	    break;
	  }
	if (ord || (m_wmode&16)) return true;
	msg_Debugging()<<"reject ordering\n";
      }
    }
    p_ampl=ampl;
    p_ampl->DeleteNext();
  } while (oldsize<nocl.size());
  msg_Debugging()<<"trying unordered configurations\n";
  if (ampl->Legs().size()==4) return false;
  if (nocl.empty()) return false;
  Vertex_Set nonocl;
  while (true) {
    double nmin(std::numeric_limits<double>::max()), pmin(nmin);
    Vertex_Set::iterator nwin(nocl.end()), pwin(nwin);
    for (Vertex_Set::iterator
	   vit(nocl.begin());vit!=nocl.end();++vit) {
      if (nonocl.find(vit->first)!=nonocl.end()) continue;
      if (vit->second<0.0) {
	if (-vit->second<nmin) {
	  nmin=-vit->second;
	  nwin=vit;
	}
      }
      else {
	if (vit->second<pmin) {
	  pmin=vit->second;
	  pwin=vit;
	}
      }
    }
    if (nwin==nocl.end()) nwin=pwin;
    if (nwin==nocl.end()) {
      p_ampl=ampl;
      return false;
    }
    nonocl[nwin->first]=nwin->second;
    nocl.erase(nwin);
    Current_Vector nccurs(ccurs);
    Current *nfcur(fcur);
    ClusterInfo_Map ncinfo(cinfo);
    if (ClusterStep(step,nocl,nccurs,nfcur,ncinfo,kt2))
      if (p_ampl->KT2()<sqrt(std::numeric_limits<double>::max()))
	if (Cluster(step+1,nocl,nccurs,nfcur,ncinfo,kt2ord)) return true;
    p_ampl=ampl;
    p_ampl->DeleteNext();
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
