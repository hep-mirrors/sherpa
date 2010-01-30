#include "COMIX/Cluster/Cluster_Algorithm.H"

#include "COMIX/Cluster/Color_Setter.H"
#include "COMIX/Main/Single_Process.H"
#include "COMIX/Amplitude/Amplitude.H"
#include "COMIX/Amplitude/Matrix_Element.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "COMIX/Phasespace/PS_Channel.H"
#include "ATOOLS/Phys/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm(ATOOLS::Mass_Selector *const ms):
  m_cs(this), p_ms(ms), p_ampl(NULL), p_clus(NULL),
  m_cnt(0), m_rej(0), m_lfrac(0.0)
{
}

Cluster_Algorithm::~Cluster_Algorithm()
{
  if (p_ampl) p_ampl->Delete();
}

template <class SType> ColorID
Cluster_Algorithm::GetPColor(Current_Base *const j,
				  Current_Base *const fcur) const
{
  ColorID col;
  if (j->Type()=='S') {
    const std::vector<CScalar<SType> > &cs(j->Current_Base::J<CScalar<SType> >());
    if (cs.empty()) return col;
    const CScalar<SType> &jc(cs.front());
    col=ColorID(jc(0),jc(1));
  }
  else if (j->Type()=='F') {
    const std::vector<CSpinor<SType> > &cs(j->Current_Base::J<CSpinor<SType> >());
    if (cs.empty()) return col;
    col=ColorID(0,0);
    const CSpinor<SType> &jc(cs.front());
    if (j->Flav().IsAnti()) col.m_j=jc();
    else col.m_i=jc();
  }
  else if (j->Type()=='V') {
    const std::vector<CVec4<SType> > &cs(j->Current_Base::J<CVec4<SType> >());
    if (cs.empty()) return col;
    const CVec4<SType> &jc(cs.front());
    col=ColorID(jc(0),jc(1));
  }
  else if (j->Type()=='T') {
    const std::vector<CAsT4<SType> > &cs(j->Current_Base::J<CAsT4<SType> >());
    if (cs.empty()) return col;
    const CAsT4<SType> &jc(cs.front());
    col=ColorID(jc(0),jc(1));
  }
  else {
    THROW(fatal_error,"Internal error");
  }
  if (j==fcur) col=col.Conj();
  return col;
}

ColorID Cluster_Algorithm::GetColor(Current_Base *const j,
					 Current_Base *const fcur) const
{
  if (p_bg->PMode()=='Q') return GetPColor<long double>(j,fcur);
  return GetPColor<double>(j,fcur);
}

bool Cluster_Algorithm::EWConnected
(const ATOOLS::Flavour &c,const ATOOLS::Flavour &s) const
{
  if (c.IntCharge()==0) return s.IntCharge()!=0;
  return c.IntCharge()*s.IntCharge()<0;
}

bool Cluster_Algorithm::
ColorConnected(const ColorID &i,const ColorID &j,const ColorID &k) const
{
  if (i.m_i>0 && k.m_j==i.m_i) return true;
  if (i.m_j>0 && k.m_i==i.m_j) return true;
  if (j.m_i>0 && k.m_j==j.m_i) return true;
  if (j.m_j>0 && k.m_i==j.m_j) return true;
  if ((k.m_i>0)^(k.m_j>0)) {
    // coloured singlet
    if (i.m_i>0 && i.m_i==j.m_j && 
	i.m_j>0 && i.m_j==j.m_i) return true;
    // colourless singlet
    if (i.m_i==0 && i.m_j==0 &&
	j.m_i==0 && j.m_j==0) return true;
  }
  // all colourless
  if (k.m_i==0 && k.m_j==0 &&
      i.m_i==0 && i.m_j==0 &&
      j.m_i==0 && j.m_j==0) return true;
  return false;
}

void Cluster_Algorithm::SwapID
(Cluster_Leg *const li,Cluster_Leg *const lj) const
{
  size_t idi(li->Id());
  li->SetId(lj->Id());
  lj->SetId(idi);
}

CParam Cluster_Algorithm::GetMeasure
(const size_t &idi,const size_t &idj,const size_t &idk,
 const ATOOLS::Flavour &mofl,Double_Map &kt2,const SizeT_Map &cid)
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
  if (m_swap) SwapID(p_ampl->Leg(0),p_ampl->Leg(1));
  bool ismo(idi&((1<<p_xs->NIn())-1));
  Flavour mmofl(p_xs->ReMap(ismo?mofl.Bar():mofl));
  if (ismo) mmofl=mmofl.Bar();
  if (p_ampl->Legs().size()>4) {
    kt2[idi][idj][idk][mofl]=
      p_clus->KPerp2(*p_ampl,i,j,k,mmofl,p_ms);
  }
  else {
    Cluster_Leg *li(p_ampl->Leg(i)), *lj(p_ampl->Leg(j));
    double ckt2(0.0);
    if ((mmofl.Resummed() || mmofl.Strong()) &&
	(li->Flav().Resummed() || li->Flav().Strong()) && 
	(lj->Flav().Resummed() || lj->Flav().Strong()) && 
	(p_ampl->Leg(k)->Flav().Resummed() || 
	 p_ampl->Leg(k)->Flav().Strong()))
      ckt2=Max(li->Mom().MPerp2(),lj->Mom().MPerp2());
    else ckt2=dabs((li->Mom()+lj->Mom()).Abs2());
    kt2[idi][idj][idk][mofl]=CParam(ckt2,ckt2);
  }
  if (m_swap) SwapID(p_ampl->Leg(0),p_ampl->Leg(1));
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
 const Current_Vector &ccurs,Current_Base *const fcur,
 ClusterInfo_Map &cinfo,const SizeT_Map &cid)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  Double_Map kt2;
  ClusterInfo_Map ccinfo(cinfo);
  cinfo.clear();
  for (size_t nc(2);nc<=step;++nc) {
    const Current_Vector &curs(p_bg->Currents()[nc]);
    for (size_t i(0);i<curs.size();++i) {
      const Vertex_Vector &in(curs[i]->In()); 
      for (size_t j(0);j<in.size();++j) {
	if (in[j]->Zero()&&!m_nosol) continue;
	if (nocl.find(in[j])!=nocl.end()) continue;
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
	  ColorID colk(p_ampl->Leg(k)->Col());
	  if (p_ampl->Legs().size()==4 ||
	      (in[j]->OrderQCD()==0?
	       EWConnected(in[j]->JC()->Flav(),p_ampl->Leg(k)->Flav()):
	       ColorConnected(coli,colj,colk))) {
	    CParam ckt2(GetMeasure(m_id[idi],m_id[idj],idk,
				   in[j]->JC()->Flav(),kt2,cid));
	    cinfo.insert(ClusterInfo_Pair
			 (Cluster_Key(idi,idj),
			  Cluster_Info(in[j],idk,ckt2,in[j]->OrderEW(),
				       in[j]->OrderQCD(),in[j]->JC()->Flav())));
	  }
	}
	}
	msg_Debugging()<<"}\n";
      }
    }
    const Vertex_Vector &in(fcur->In()); 
    for (size_t j(0);j<in.size();++j) {
      if (in[j]->Zero()&&!m_nosol) continue;
      if (nocl.find(in[j])!=nocl.end()) continue;
      for (size_t i(1);i<ccurs.size();++i) {
	if (in[j]->JA()==ccurs[i] || in[j]->JB()==ccurs[i]) {
	  if (ccurs[i]->CId()&2) continue;
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
	      ColorID colk(p_ampl->Leg(k)->Col());
	      Flavour mofl((in[j]->JA()==ccurs[i]?
			    in[j]->JB():in[j]->JA())->Flav().Bar());
	      if (p_ampl->Legs().size()==4 ||
		  (in[j]->OrderQCD()==0?
		   EWConnected(mofl,p_ampl->Leg(k)->Flav()):
		   ColorConnected(coli,colj,colk))) {
		CParam ckt2(GetMeasure(m_id[idi],m_id[idj],idk,mofl,kt2,cid));
		cinfo.insert(ClusterInfo_Pair
			     (Cluster_Key(idi,idj),
			      Cluster_Info(in[j],idk,ckt2,in[j]->OrderEW(),
					   in[j]->OrderQCD(),mofl)));
	      }
	    }
	  }
	  msg_Debugging()<<"}\n";
	}
      }
    }
  }
  msg_Debugging()<<"}\n";
}

bool Cluster_Algorithm::CombineWinner
(Vertex_Base *const v,Current_Vector &ccurs,
 Current_Base *&fcur,ClusterInfo_Map &cinfo)
{
  if (v->JC()!=fcur) {
    Current_Base *ja(v->JA()), *jb(v->JB());
    m_id[v->JC()->CId()]=m_id[ja->CId()]+m_id[jb->CId()];
    if (v->JA()->Id().front()>v->JB()->Id().front()) 
      std::swap<Current_Base*>(ja,jb);
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
		   <<ID(m_id[v->JC()->CId()])<<"\n";
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
		   <<" -> "<<ID(m_id[v->JA()->CId()])
		   <<"&"<<ID(m_id[v->JB()->CId()])<<"\n";
  }
  return true;
}

bool Cluster_Algorithm::ClusterStep
(const size_t &step,Vertex_Set &nocl,
 Current_Vector &ccurs,Current_Base *&fcur,
 ClusterInfo_Map &cinfo)
{
  msg_Debugging()<<METHOD<<"(): step = "<<step<<" {\n";
  msg_Indent();
  SizeT_Map cid;
  for (size_t i(0);i<p_ampl->Legs().size();++i) 
    cid[p_ampl->Leg(i)->Id()]=i;
  CalculateMeasures(step,nocl,ccurs,fcur,cinfo,cid);
  if (cinfo.empty()) {
    msg_Debugging()<<"rejected configuration\n";
    return false;
  }
  double rwmin(std::numeric_limits<double>::max()), sum(0.0);
  ClusterInfo_Map::const_iterator win(cinfo.end()), rwin(win);
  for (ClusterInfo_Map::const_iterator cit(cinfo.begin());
       cit!=cinfo.end();++cit)
    if (cit->second.m_kt2.m_op2>=0.0) {
      sum+=1.0/cit->second.m_kt2.m_op2;
    }
    else if (cit->second.m_kt2.m_kt2>=0.0 &&
	     cit->second.m_kt2.m_kt2<rwmin) {
      rwin=cit;
      rwmin=cit->second.m_kt2.m_kt2;
    }
  double disc(sum*ran.Get()), psum(0.0);
  for (ClusterInfo_Map::const_iterator cit(cinfo.begin());
       cit!=cinfo.end();++cit)
    if (cit->second.m_kt2.m_op2>=0.0 &&
	(psum+=1.0/cit->second.m_kt2.m_op2)>=disc) {
      win=cit;
      break;
    }
  if (sum>0.0 && win==cinfo.end()) THROW(fatal_error,"Internal error"); 
  if (win==cinfo.end()) win=rwin;
  if (win==cinfo.end()) THROW(fatal_error,"Invalid amplitude");
  Cluster_Key wkey(win->first);
  Cluster_Info winfo(win->second);
  nocl.insert(winfo.p_v);
  if (!CombineWinner(winfo.p_v,ccurs,fcur,cinfo)) return false;
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
    if (!match) {
      msg_Debugging()<<"Invalid core\n";
      return false;
    }
  }
  Vec4D_Vector p;
  if (p_ampl->Legs().size()>4) {
    if (m_swap) SwapID(p_ampl->Leg(0),p_ampl->Leg(1));
    p=p_clus->Combine(*p_ampl,cid[m_id[wkey.first]],
		      cid[m_id[wkey.second]],cid[winfo.m_k],
		      winfo.m_mofl,p_ms,winfo.m_kt2.m_kin);
    if (m_swap) SwapID(p_ampl->Leg(0),p_ampl->Leg(1));
    if (p.empty()) return false;
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
  ampl->SetKT2QCD(winfo.m_kt2.m_kt2);
  p_ampl=p_ampl->InitNext();
  p_ampl->SetMS(p_ms);
  p_ampl->SetNIn(ampl->NIn());
  p_ampl->SetMuR2(ampl->MuR2());
  p_ampl->SetMuF2(ampl->MuF2());
  p_ampl->SetX1(ampl->X1());
  p_ampl->SetX2(ampl->X2());
  size_t nid(m_id[wkey.first]+m_id[wkey.second]);
  if (nid&3) {
    if (nid&1) p_ampl->SetX1(ampl->X1()*winfo.m_kt2.m_x);
    else p_ampl->SetX2(ampl->X2()*winfo.m_kt2.m_x);
  }
  else if (winfo.m_k&3) {
    if (winfo.m_k&1) p_ampl->SetX1(ampl->X1()*winfo.m_kt2.m_x);
    else p_ampl->SetX2(ampl->X2()*winfo.m_kt2.m_x);
  }
  p_ampl->SetJF(ampl->JF<Selector_Base>());
  p_ampl->SetOrderEW(ampl->OrderEW()-winfo.p_v->OrderEW());
  p_ampl->SetOrderQCD(ampl->OrderQCD()-winfo.p_v->OrderQCD());
  p_ampl->SetKin(winfo.m_kt2.m_kin);
  for (size_t i(0);i<ccurs.size();++i) {
    size_t cid(m_id[ccurs[i]->CId()]);
    Flavour flav(p_xs->ReMap(ccurs[i]->Flav()));
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
    else if (col.m_i<0 && ccurs[i]->Cut()) {
      p_ampl->Legs().back()->SetStat(3);
      SetNMax(p_ampl->Prev(),cid,ccurs[i]->Cut());
      p_ampl->Legs().back()->SetDMax(ccurs[i]->Cut());
    }
    if (col.m_i<0) {
      p_ampl->Legs().back()->SetCol(GetColor(ccurs[i],fcur));
      p_ampl->Legs().back()->SetK(winfo.m_k);
    }
  }
  // set colour connections 
  msg_Debugging()<<"} step = "<<step<<"\n";
  return true;
}

bool Cluster_Algorithm::Cluster(Single_Process *const xs)
{
  p_bg=(p_xs=xs)->GetAmplitude();
  m_swap=xs->Process()->Integrator()->InSwaped();
  if (p_bg==NULL) THROW(fatal_error,"Internal error");
  Selector_Base *jf=p_xs->Selector()
    ->GetSelector("Jetfinder");
  bool trig(true);
  if (jf) {
    Vec4D_Vector moms(xs->Process()->Integrator()->Momenta());
    if (m_swap) {
      std::swap<Vec4D>(moms[0],moms[1]);
      for (size_t i(0);i<moms.size();++i)
	moms[i]=Vec4D(moms[i][0],-moms[i]);
    }
    trig=jf->Trigger(moms);
  }
  msg_Debugging()<<METHOD<<"(): trig = "<<trig<<" {\n";
  msg_Indent();
  m_id.clear();
  Current_Vector ccurs(p_bg->Currents()[1]);
  Current_Base *fcur(ccurs[0]=p_bg->Currents().back().front());
  if (p_ampl) p_ampl->Delete();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(p_ms);
  p_ampl->SetJF(jf);
  p_ampl->SetNIn(xs->NIn());
  p_ampl->SetOrderEW(p_bg->MaxOrderEW());
  p_ampl->SetOrderQCD(p_bg->MaxOrderQCD());
  double muf2(xs->Process()->ScaleSetter()->Scale(stp::fac));
  double mur2(xs->Process()->ScaleSetter()->Scale(stp::ren));
  for (size_t i(0);i<ccurs.size();++i) {
    size_t cid(m_id[ccurs[i]->CId()]=1<<p_ampl->Legs().size());
    ColorID col(GetColor(i==0?p_bg->Currents()[1][0]:ccurs[i],fcur));
    Flavour flav(p_xs->ReMap(ccurs[i]->Flav()));
    if (ccurs[i]==fcur) flav=flav.Bar();
    size_t idx(i<2?(m_swap?1-i:i):i);
    Vec4D mom(i<2?-xs->Process()->Integrator()->Momenta()[idx]:
	      xs->Process()->Integrator()->Momenta()[idx]);
    p_ampl->CreateLeg(mom,flav,col,cid);
    p_ampl->Legs().back()->SetStat(1);
  }
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  p_ampl->SetX1(xs->Process()->Integrator()->ISR()->X1());
  p_ampl->SetX2(xs->Process()->Integrator()->ISR()->X2());
  if (m_swap) xs->Process()->SwapInOrder();
  m_nosol=!m_cs.SetColors(xs);
  if (m_swap) xs->Process()->SwapInOrder();
  ClusterInfo_Map cinfo;
  ++m_cnt;
  if (!Cluster(2,Vertex_Set(),ccurs,fcur,cinfo)) {
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
  SetNMax(p_ampl,(1<<ccurs.size())-1,
	  trig?xs->Process()->Info().m_fi.NMaxExternal():
	  xs->Process()->Info().m_fi.NExternal());
  if (msg_LevelIsDebugging()) p_ampl->Print();
  while (p_ampl->Prev()) {
    if (m_swap) p_ampl->SwapInOrder();
    p_ampl=p_ampl->Prev();
    if (msg_LevelIsDebugging()) p_ampl->Print();
  }
  if (m_swap) p_ampl->SwapInOrder();
  return true;
}

bool Cluster_Algorithm::Cluster
(const size_t &step,const Vertex_Set &onocl,const Current_Vector &ccurs,
 Current_Base *const fcur,const ClusterInfo_Map &cinfo)
{
  if (p_ampl->Legs().size()==3) {
    p_ampl=p_ampl->Prev();
    p_ampl->DeleteNext();
    return true;
  }
  size_t oldsize(0);
  Vertex_Set nocl;
  Cluster_Amplitude *ampl(p_ampl);
  do {
    oldsize=nocl.size();
    Current_Vector nccurs(ccurs);
    Current_Base *nfcur(fcur);
    ClusterInfo_Map ncinfo(cinfo);
    if (ClusterStep(step,nocl,nccurs,nfcur,ncinfo))
      if (Cluster(step+1,nocl,nccurs,nfcur,ncinfo)) {
#ifdef METS__reject_unordered
 	if (ampl->Legs().size()==4 || 
 	    ampl->KT2QCD()<p_ampl->KT2QCD()) return true;
#else
	return true;
#endif
      }
    p_ampl=ampl;
  } while (oldsize<nocl.size());
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
