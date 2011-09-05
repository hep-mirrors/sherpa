#include "COMIX/Cluster/Color_Setter.H"

#include "COMIX/Cluster/Cluster_Algorithm.H"
#include "COMIX/Main/Single_Process.H"
#include "COMIX/Amplitude/Matrix_Element.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flow.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

size_t s_clmaxtrials(900);

Color_Setter::Color_Setter(Cluster_Algorithm *const ca): 
  p_ca(ca) 
{
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(m_cmode,"COMIX_CCMODE")) m_cmode=1;
}

bool Color_Setter::SetRandomColors
(const Int_Vector &ci,const Int_Vector &cj)
{
  size_t trials(0);
  Cluster_Amplitude *ampl(p_ca->GetAmplitude());
  std::vector<ColorID> oc(ampl->Legs().size());
  for (size_t i(0);i<ampl->Legs().size();++i)
    ampl->Leg(i)->SetCol(oc[i]=ColorID(ci[i],cj[i]));
  msg_Debugging()<<*ampl<<"\n";
  for (;trials<s_clmaxtrials;++trials) {
    bool sing(false);
    std::set<size_t> cs;
    Int_Vector ci(ampl->Legs().size()), cj(ampl->Legs().size());
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
    double csum(p_xs->GetME()->Differential
		(p_xs->Integrator()->Momenta(),ci,cj,true));
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
      Int_Vector ni(colint->I()), nj(colint->J());
      for (size_t i(0);i<ampl->Legs().size();++i)
	oc[i]=ColorID(ni[i],nj[i]);
    }
  }
  if (trials<s_clmaxtrials) {
    p_xs->GetAmplitude()->ResetZero();
  }
  else {
    msg_Error()<<METHOD<<"(): No solution."<<std::endl;
    return false;
  }
  return true;
}

bool Color_Setter::SetLargeNCColors
(const Int_Vector &ci,const Int_Vector &cj)
{
  Cluster_Amplitude *ampl(p_ca->GetAmplitude());
  std::vector<ColorID> oc(ampl->Legs().size());
  for (size_t i(0);i<ampl->Legs().size();++i)
    ampl->Leg(i)->SetCol(oc[i]=ColorID(ci[i],cj[i]));
  SP(Color_Integrator) colint(p_xs->Integrator()->ColorIntegrator());
  bool sotfcc(colint->OTFCC());
  colint->SetOTFCC(true);
  Double_Vector psum;
  std::map<size_t,size_t> cmap;
  size_t ntrials(0);
  while (true) {
    size_t cc(0);
    psum.clear();
    while (colint->NextOrder()) {
      ++cc;
      bool valid(true);
      Idx_Vector perm(colint->Orders().front());
      Int_Vector ci(perm.size(),0), cj(perm.size(),0);
      int lc(p_xs->Flavours()[perm.back()].StrongCharge());
      if (perm.back()<p_xs->NIn() && lc!=8) lc=-lc;
      for (size_t i(0);i<perm.size();++i) {
	size_t cur(perm[i]), next(i<perm.size()-1?perm[i+1]:perm[0]);
	int sc(p_xs->Flavours()[cur].StrongCharge());
	if (cur<p_xs->NIn() && sc!=8) sc=-sc;
	if (sc==3 || sc==8)
	  cj[next]=ci[cur]=Flow::Counter();
	if (lc!=0) 
	  if (((lc==8 || lc==3) && sc==3) ||
	      (lc==-3 && (sc==-3 || sc==8))) {
	    valid=false;
	    break;
	  }
	if (sc!=0) lc=sc;
      }
      if (!valid) continue;
      double part(p_xs->GetME()->Differential
		  (p_xs->Integrator()->Momenta(),ci,cj,true));
      part*=sqr(colint->Weights().front());
      if (psum.size()) part+=psum.back();
      if (part>0.0) {
	cmap[psum.size()]=cc;
	psum.push_back(part);
      }
      if (psum.size()) 
	msg_Debugging()<<"psum["<<psum.size()-1<<"]/["<<cc<<"]"
		       <<perm<<" = "<<psum.back()<<"\n";
    }
    if (psum.size()) break;
    // select new color configuration
    while (!colint->GeneratePoint());
    Int_Vector ni(colint->I()), nj(colint->J());
    for (size_t i(0);i<ampl->Legs().size();++i)
      oc[i]=ColorID(ni[i],nj[i]);
    if ((ntrials+=10)>=s_clmaxtrials) {
      msg_Error()<<METHOD<<"(): No nonzero amplitude. Randomize."<<std::endl;
      return false;
    }
  }
  msg_Debugging()<<"sum = "<<psum.back()<<"\n";
  size_t l(0), r(psum.size()-1), c((l+r)/2);
  double a(psum[c]), disc(ran->Get()*psum.back());
  while (r-l>1) {
    if (disc<a) r=c;
    else l=c;
    c=(l+r)/2;
    a=psum[c];
  }
  if (disc<psum[l]) r=l;
  size_t ck(cmap[r]), cc(0);
  msg_Debugging()<<"selected r = "<<r<<", ck = "<<ck<<": "<<psum[r]<<"\n";
  while (colint->NextOrder()) {
    ++cc;
    if (ck==cc) {
      Idx_Vector perm(colint->Orders().front());
      Int_Vector ci(perm.size(),0), cj(perm.size(),0);
      int lc(p_xs->Flavours()[perm.back()].StrongCharge());
      if (perm.back()<p_xs->NIn() && lc!=8) lc=-lc;
      for (size_t i(0);i<perm.size();++i) {
	size_t cur(perm[i]), next(i<perm.size()-1?perm[i+1]:perm[0]);
	int sc(p_xs->Flavours()[cur].StrongCharge());
	if (cur<p_xs->NIn() && sc!=8) sc=-sc;
	if (sc==3 || sc==8)
	  cj[next]=ci[cur]=Flow::Counter();
	if (lc!=0) 
	  if (((lc==8 || lc==3) && sc==3) ||
	      (lc==-3 && (sc==-3 || sc==8))) {
	    THROW(fatal_error,"Internal error");
	  }
	if (sc!=0) lc=sc;
      }
      for (size_t i(0);i<ampl->Legs().size();++i)
	ampl->Leg(i)->SetCol(ColorID(ci[i],cj[i]));
      if (p_xs->GetME()->Differential
	  (p_xs->Integrator()->Momenta(),ci,cj,true)==0.0)
	THROW(fatal_error,"Internal error");
      CI_Map &cmap(ampl->ColorMap());
      for (size_t i(0);i<ampl->Legs().size();++i)
	if (oc[i].m_i!=0) cmap[ci[i]]=oc[i].m_i;
      break;
    }
  }
  colint->SetOTFCC(sotfcc);
  p_xs->GetAmplitude()->ResetZero();
  return true;
}

bool Color_Setter::SetColors(Single_Process *const xs)
{
  p_xs=xs;
  SP(Color_Integrator) colint(p_xs->Integrator()->ColorIntegrator());
  Int_Vector ci(colint->I()), cj(colint->J());
  bool sol(false);
  switch (m_cmode) {
  case 1: {
    sol=SetRandomColors(ci,cj);
    if (!sol) sol=SetLargeNCColors(ci,cj);
    break;
  } 
  case 2: 
    sol=SetLargeNCColors(ci,cj);
    if (!sol) sol=SetRandomColors(ci,cj);
    break;
  default:
    THROW(fatal_error,"Invalid colour setting mode");
  }
  colint->SetPoint(ci,cj);
  return sol;
}

