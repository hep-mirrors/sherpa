#include "AMEGIC++/Amplitude/FullAmplitude_External.H"

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "AMEGIC++/Main/Helicity.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

FullAmplitude_External::FullAmplitude_External
(const Process_Info &info,Model_Base *model,Coupling_Map *const cpls,
 Helicity *hel,const size_t &emit,const size_t &spect): 
  p_model(model->GetInteractionModel()),
  m_emit(emit), m_spect(spect)
{ 
  DEBUG_FUNC("");
  p_calc = Tree_ME2_Base::GetME2(info);
  if (p_calc==NULL) return;
  msg_Debugging()<<p_calc<<" -> "<<info<<"\n";
  p_calc->SetCouplings(*cpls);
  m_oqcd=p_calc->OrderQCD();
  m_oew=p_calc->OrderEW();
  m_nin=info.m_ii.NExternal();
  m_fls=info.ExtractFlavours();
  for (size_t i(0);i<m_nin;++i) m_fls[i]=m_fls[i].Bar();
  m_amap=p_calc->GetFlavourHelicityMap();
  BuildHelicityMap(hel);
  BuildColorMatrix();
}

FullAmplitude_External::~FullAmplitude_External() 
{ 
  delete p_calc;
}

double FullAmplitude_External::Calc(const ATOOLS::Vec4D *_p)
{
  DEBUG_FUNC("");
  Vec4D_Vector p(m_fls.size());
  for (size_t i(0);i<p.size();++i) p[i]=_p[i];
  double res(p_calc->Calc(p));
  m_amps=p_calc->GetAmplitudes();
#ifdef DEBUG__External
  for (size_t i(0);i<p.size();++i)
    msg_Debugging()<<"p["<<i<<"]="<<p[i]<<"\n";  
  for (size_t i(0);i<m_amap.size();++i) {
    msg_Debugging()<<"A[";
    for (size_t k(0);k<m_fls.size();++k)
      msg_Debugging()<<m_amap[i].m_perm[k]%100;
    msg_Debugging()<<"](";
    for (size_t k(0);k<m_fls.size();++k)
      msg_Debugging()<<(k==0?"":",")<<m_amap[i].m_hels[k];
    msg_Debugging()<<") = "<<m_amps[i]
		   <<" -> "<<std::abs(m_amps[i])<<"\n";
  }
#endif
  return res;
}

double FullAmplitude_External::MSquare(const size_t &h) const
{
  DEBUG_FUNC(h);
  Complex ampsqpp(0.,0.), ampsqmm(0.,0.);
  Complex ampsqpm(0.,0.), ampsqmp(0.,0.);
  for (size_t k=0;k<m_hmap[h].size();++k) {
    size_t i=m_hmap[h][k];
    for (size_t l=0;l<m_hmap[h].size();++l) {
      size_t j=m_hmap[h][l];
      Complex amp(m_amps[i]);
      amp*=std::conj(m_amps[j]);	
      amp*=m_colfs[m_emit][m_spect][i][j];
      msg_Debugging()<<m_colfs[m_emit][m_spect][i][j]
		     <<"*"<<m_amps[i]<<"*"
		     <<std::conj(m_amps[j])<<" = "<<amp
		     <<(m_pmap[h][k]>0?" +":" -")
		     <<(m_pmap[h][l]>0?'+':'-')<<"\n";
      if (m_pmap[h][k]<0) {
	if (m_pmap[h][l]<0) ampsqmm+=amp;
	else ampsqmp+=amp;
      }
      else {
	if (m_pmap[h][l]>0) ampsqpp+=amp;
	else ampsqpm+=amp;
      }
    }
  }
  Complex sum(ampsqpp);
  if (m_emit==m_spect) {
    msg_Debugging()<<"AA* = "<<sum<<"\n"; 
  }
  else {
  if (m_fls[m_emit].IsVector()) {
    Complex pp(m_eip*m_eip);
    sum=0.5*(1.0+m_A)*(ampsqmm+ampsqpp)+
      0.5*(1.0-m_A)*(ampsqpm/pp+ampsqmp*pp);
  }
  msg_Debugging()<<"AA*: ++ = "<<ampsqpp
		 <<", +- = "<<ampsqpm 
		 <<", -+ = "<<ampsqmp 
		 <<", -- = "<<ampsqmm<<" -> "<<sum.real()
		 <<" ( A = "<<m_A<<", eip = "<<m_eip<<" )\n"; 
  }
  return sum.real();
}

void FullAmplitude_External::SetSqMatrix
(const double &A,const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &eps1)
{
  m_A=A;
  m_eip=p_calc->GetHelicityPhase(pijt,eps1);
}

double FullAmplitude_External::MSquare
(const size_t &h,const size_t &ci,const size_t &cj)
{
  DEBUG_FUNC("h="<<h<<", "<<ci<<"<>"<<cj);
  Complex sum(0.,0.);
  for (size_t k=0;k<m_hmap[h].size();++k) {
    size_t i=m_hmap[h][k];
    for (size_t l=0;l<m_hmap[h].size();++l) {
      size_t j=m_hmap[h][l];
      if (m_pmap[h][k]!=m_pmap[h][l]) continue;
      Complex amp(m_amps[i]);
      amp*=std::conj(m_amps[j]);	
      amp*=m_colfs[ci][cj][i][j];
      msg_Debugging()<<m_colfs[ci][cj][i][j]
		     <<"*"<<m_amps[i]<<"*"
		     <<std::conj(m_amps[j])<<" = "<<amp
		     <<(m_pmap[h][k]>0?" +":" -")
		     <<(m_pmap[h][l]>0?'+':'-')<<"\n";
      sum+=amp;
    }
  }
  msg_Debugging()<<"AA* = "<<sum.real()<<"\n"; 
  return sum.real();
}

void FullAmplitude_External::BuildHelicityMap(Helicity *const hel)
{
  DEBUG_FUNC("");
  m_hmap.resize(hel->MaxHel());
  m_pmap.resize(hel->MaxHel());
  for (size_t k=0;k<m_amap.size();++k) {
    std::vector<int> chel(m_fls.size());
    for (size_t j=0;j<chel.size();++j)
      chel[m_amap[k].m_perm[j]%100]=m_amap[k].m_hels[j];
    for (size_t i=0;i<hel->MaxHel();++i) {
      int epol(hel->GetEPol(i)==90?1:-1);
      std::vector<int> rhel((*hel)[i],&(*hel)[i][chel.size()]);
      for (size_t l=0;l<rhel.size();++l)
	if (rhel[l]>=90) epol=rhel[l]=chel[l];
      if (chel==rhel) {
	msg_Debugging()<<"  map "<<k<<" -> "<<i<<" <=> "
		       <<m_amap[k]<<" -> "<<rhel<<", ep = "<<epol<<"\n";
	m_hmap[i].push_back(k);
	m_pmap[i].push_back(epol);
      }
    }
  }
}

void FullAmplitude_External::GetPermutation
(const std::vector<size_t> &ids,std::vector<size_t> &cid,
 Flavour_Vector &cfl,int &nsub,int &psub,int &swap) const
{
  for (size_t i(0);i<ids.size();++i) {
    int id=ids[i]%100;
    if (m_fls[id].StrongCharge()==0) continue;
    cid.push_back(ids[i]);
    cfl.push_back(m_fls[id]);
  }
  psub=nsub=0;
  int rid(-1);
  for (size_t i(0);i<cid.size();++i) {
    int id=cid[i]%100, aid=cid[i]/100;
    if (cfl[i].StrongCharge()==3) rid=aid;
    if (cfl[i].StrongCharge()==-3 &&
	i<cid.size()-1 && aid==rid) {
      if (m_fls[cid[i+1]%100]==m_fls[id].Bar()) ++psub;
      else ++nsub;
    }
    cid[i]=id;
  }
  swap=0;
  std::vector<size_t> fid(cid);
  for (int oswap(-1);oswap<swap;) {
    oswap=swap;
    for (size_t i(1);i<fid.size();++i)
      if (fid[i]<fid[i-1]) {
	std::swap<size_t>(fid[i],fid[i-1]);
	if (m_fls[fid[i]].IsFermion() &&
	    m_fls[fid[i]]==m_fls[fid[i-1]]) ++swap;
      }
  }
}

void FullAmplitude_External::BuildColorMatrix() 
{
  DEBUG_FUNC(m_amap.size());
  m_amps.resize(m_amap.size());
  m_colfs.resize(m_fls.size());
  for (size_t i(0);i<m_colfs.size();++i)
    if (m_fls[i].Strong()) m_colfs[i].resize(m_fls.size());
  BuildColorMatrix(m_emit,m_spect);
  if (m_emit==m_spect) {
    for (size_t i(0);i<m_colfs.size();++i) {
      if (!m_fls[i].Strong()) continue;
      for (size_t j(0);j<m_colfs[i].size();++j) {
	if (i==j || !m_fls[j].Strong()) continue;
	BuildColorMatrix(i,j);
      }
    }
  }
}

void FullAmplitude_External::BuildColorMatrix
(const size_t &ci,const size_t &cj) 
{
  m_colfs[ci][cj].
    resize(m_amap.size(),std::vector<Complex>
	   (m_amap.size(),Complex(0.0,0.0)));
  for (size_t a=0;a<m_amap.size();++a) {
    const std::vector<int> &cp(m_amap[a].m_perm);
    int nsuba, psuba, swapa;
    Flavour_Vector fla;
    std::vector<size_t> ida, perma(cp.begin(),cp.end());
    GetPermutation(perma,ida,fla,nsuba,psuba,swapa);
    msg_Debugging()<<"A: "<<perma<<ida<<fla<<" sub = -"<<nsuba<<"/+"<<psuba
		   <<" swap = "<<swapa<<" ( em = "<<ci<<" )\n";
    msg_Indent();
    for (size_t b=a;b<m_amap.size();++b) {
      const std::vector<int> &cp(m_amap[b].m_perm);
      int nsubb, psubb, swapb;
      Flavour_Vector flb;
      std::vector<size_t> idb, permb(cp.begin(),cp.end());
      GetPermutation(permb,idb,flb,nsubb,psubb,swapb);
      msg_Debugging()<<"B: "<<permb<<idb<<flb<<" sub = -"<<nsubb<<"/+"<<psubb
		     <<" swap = "<<swapb<<" ( sp = "<<cj<<" )\n";
      msg_Indent();
      Expression expression(100,100);
      expression.pop_back();
      expression.SetTR(p_calc->TR());
      size_t ad(expression.AIndex()), lf(expression.FIndex()), fl(lf);
      for (size_t i=0;i<fla.size();++i) {
	if (fla[i].StrongCharge()==3) {
	  if (i>0) lf=expression.FIndex();
	  expression.push_back(Delta::New(ida[i]+1,lf));
	  if (ci!=cj) {
	  if (ci==ida[i]) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(ad,lf,nlf));
	    lf=nlf;
	  }
	  }
	}
	else if (fla[i].StrongCharge()==-3) {
	  if (ci!=cj) {
	  if (ci==ida[i]) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(ad,lf,nlf));
	    expression.push_back(CNumber::New(Complex(-1.0)));
	    lf=nlf;
	  }
	  }
	  expression.push_back(Delta::New(lf,ida[i]+1));
	}
	else {
	  if (ci==cj) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(ida[i]+1,lf,nlf));
	    lf=nlf;
	  }
	  else {
	  if (ci!=ida[i]) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(ida[i]+1,lf,nlf));
	    lf=nlf;
	  }
	  else {
	    size_t la=expression.AIndex();
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(la,lf,nlf));
	    expression.push_back(Adjoint::New(ad,ida[i]+1,la));
	    expression.push_back(CNumber::New(Complex(0.0,-1.0)));
	    lf=nlf;
	  }
	  }
	}
      }
      if (fla.back().StrongCharge()==3 ||
	  fla.back().StrongCharge()==8)
	expression.push_back(Delta::New(lf,fl));
      fl=lf=expression.FIndex();
      for (size_t i=0;i<flb.size();++i) {
	if (flb[i].StrongCharge()==3) {
	  if (i>0) lf=expression.FIndex();
	  expression.push_back(Delta::New(lf,idb[i]+1));
	  if (ci!=cj) {
	  if (cj==idb[i]) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(ad,nlf,lf));
	    lf=nlf;
	  }
	  }
	}
	else if (flb[i].StrongCharge()==-3) {
	  if (ci!=cj) {
	  if (cj==idb[i]) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(ad,nlf,lf));
	    expression.push_back(CNumber::New(Complex(-1.0)));
	    lf=nlf;
	  }
	  }
	  expression.push_back(Delta::New(idb[i]+1,lf));
	}
	else {
	  if (ci==cj) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(idb[i]+1,nlf,lf));
	    lf=nlf;
	  }
	  else {
	  if (cj!=idb[i]) {
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(idb[i]+1,nlf,lf));
	    lf=nlf;
	  }
	  else {
	    size_t la=expression.AIndex();
	    size_t nlf=expression.FIndex();
	    expression.push_back(Fundamental::New(la,nlf,lf));
	    expression.push_back(Adjoint::New(ad,la,idb[i]+1));
	    expression.push_back(CNumber::New(Complex(0.0,-1.0)));
	    lf=nlf;
	  }
	  }
	}
      }
      if (fla.back().StrongCharge()==3 ||
	  flb.back().StrongCharge()==8)
	expression.push_back(Delta::New(fl,lf));
      if (msg_LevelIsDebugging()) expression.Print();
      expression.Evaluate();
      Complex col=expression.Result();
      if (nsuba || nsubb) col/=pow(-expression.NC(),nsuba+nsubb);
      if (psuba || psubb) col/=pow(expression.NC(),psuba+psubb);
      if ((swapa+swapb)%2) col=-col;
      double T2=1.0;
      if (ci!=cj) {
      if (m_fls[ci].StrongCharge()==8) {
	T2=expression.NC();
	T2/=CSC.CA;
      }
      else {
	T2=expression.TR()*(expression.NC()-1.0/expression.NC());
	T2/=CSC.CF;
      }
      }
      msg_Debugging()<<"A B^* = "<<col<<" <- sub = -"<<nsuba+nsubb
		     <<"/+"<<psuba+psubb<<", swap = "<<swapa+swapb
		     <<", T2 = "<<T2<<" ("<<a<<","<<b<<")\n";
      m_colfs[ci][cj][a][b]=col/T2;
      m_colfs[ci][cj][b][a]=conj(col)/T2;
    }
  }
}
