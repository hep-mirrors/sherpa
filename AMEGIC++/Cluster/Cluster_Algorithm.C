#include "AMEGIC++/Cluster/Cluster_Algorithm.H"

#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "PHASIC++/Process/POWHEG_Process.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace EXTRAXS;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm(ATOOLS::Mass_Selector *const ms):
  p_ms(ms), p_ampl(NULL), p_clus(NULL), p_combi(NULL) {}

Cluster_Algorithm::~Cluster_Algorithm()
{
  for (Flav_ME_Map::const_iterator xsit(m_xsmap.begin());
       xsit!=m_xsmap.end();++xsit) delete xsit->second;
  if (p_combi) delete p_combi;
}

bool Cluster_Algorithm::Cluster
(Process_Base *const xs,const size_t mode,const double &kt2)
{
  p_proc=xs;
  p_ampl=NULL;
  int nampl=p_proc->NumberOfDiagrams();
  int nlegs=p_proc->NIn()+p_proc->NOut();
  Leg **legs(CreateLegs(nampl,nlegs));
  CreateTables(legs,nampl,mode,kt2);
  ++m_cnt;
  if (p_ct==NULL) {
    msg_Debugging()<<METHOD<<"(): No valid PS history.\n";
    ++m_rej;
    double frac(m_rej/(double)m_cnt);
    if (frac>1.25*m_lfrac && m_cnt>1000) {
      m_lfrac=frac;
      msg_Error()<<METHOD<<"(): No valid PS history in >"
		 <<(int(m_lfrac*1000)/10.0)<<"% of calls.\n";
    }
    p_combi=NULL;
    return false;
  }
  Vec4D_Vector moms(4);
  ATOOLS::Flavour_Vector fl(4);
  for (int i(0);i<4;++i) {
    moms[i]=p_ct->Momentum(i);
    fl[i]=p_ct->Flav(i);
  }
  p_xs=GetXS(fl);
  SetColours(p_xs,moms,&fl.front());
  Convert();
  return true;
}

bool Cluster_Algorithm::FillLegs(Leg * alegs, Point * root, int & l, int maxl) 
{
  if (l>= maxl) {
    msg_Error()<<" Error in FillLegs() !!! "<<std::endl;
    return 0;
  }
  if (l==0) {
    size_t id(1<<root->number);
    alegs[root->number]=Leg(root);
    alegs[root->number].SetExternal(1);
    alegs[root->number].SetID(id);    
    l++;
  }
  if (root->left) {
    if (root->middle) return 0; // four vertex 
    return FillLegs(alegs,root->left,l,maxl)*FillLegs(alegs,root->right,l,maxl);
  } 
  else {
    size_t id(1<<root->number);
    alegs[root->number]=Leg(root);
    alegs[root->number].SetExternal(1);
    alegs[root->number].SetID(id);    
    l++;
    return 1;
  }
}

Leg **Cluster_Algorithm::CreateLegs(int &nampl,const int nlegs)
{
  Leg **legs(NULL);
  if (p_combi) delete p_combi;
  p_combi = 0;
  legs = new Leg *[nampl];
  for (int k=0;k<nampl;) {
    legs[k] = new Leg[nlegs];
    int l   = 0;
    if (FillLegs(legs[k],p_proc->Diagram(k),l,nlegs)) ++k;
    else {
      delete [] legs[k];
      --nampl;
    }
  }
  for (int k=0;k<nampl;++k) {
    for (int i(0);i<nlegs;++i) {
      Flavour fl(p_proc->Flavours()[i]);
      if (i<2 && p_proc->Integrator()->InSwaped()) {
	fl=p_proc->Flavours()[1-i];
      }
      legs[k][i].SetMapFlavour(fl);
//       msg_Debugging()<<"set mapfl: "<<k<<", "<<i<<": "<<fl<<" "
// 		     <<p_proc->Integrator()->InSwaped()<<"\n";
    }
  }
  return legs;
}

void Cluster_Algorithm::CreateTables
(Leg ** legs,const int nampl,const size_t mode,const double &kt2) 
{
  p_ct = 0;
  // if no combination table exist, create it
  int nin(p_proc->NIn()), nout(p_proc->NOut()), nlegs(nin+nout);
  Vec4D * amoms = new Vec4D[nlegs];
  for (int i=0;i<nin+nout;++i)  
    amoms[i]     = p_proc->Integrator()->Momenta()[i];
  if (p_proc->Integrator()->InSwaped()) {
    Vec4D help=amoms[0];
    amoms[0]=amoms[1];
    amoms[1]=help;
  }
  if (!p_combi) {
    /*
      - copy moms to insert into Combine_Table (will be delete there)
      - create new Combine_Table with given momenta and given Jet-measure
      - initialise Combine_Table
      - determine best combination sheme
    */ 
    m_decids=p_proc->DecayInfos();
    p_combi = new Combine_Table(p_proc,p_ms,p_clus,amoms,0,&m_decids);
    p_combi->FillTable(legs,nlegs,nampl);   
    p_ct = p_combi->CalcJet(nlegs,NULL,mode,kt2); 
  }
  else {
    // use the existing combination table and determine best combination sheme
    p_ct = p_combi->CalcJet(nlegs,amoms,mode,kt2);
  }
  //  delete [] amoms;
}

int Cluster_Algorithm::SetDecayColours(const Vec4D_Vector& p, Flavour * fl,int col1,int col2)
{
  int ncol   = 0;
  int nquark = 0;
  int ngluon = 0;
  for (int i=0; i<3; ++i) {
    if (fl[i].Strong()) {
      ++ncol;
      if (fl[i].IsQuark() || fl[i].IsSquark()) ++nquark;
      if (fl[i].IsGluon() || fl[i].IsGluino()) ++ngluon;
    }
  }  
  m_colors[0][0] = col1; m_colors[0][1] = col2; 
  for (int i=1; i<3; ++i) m_colors[i][0] = m_colors[i][1] = 0;
  switch (ncol) {
  case 0: 
    // no colours at all.
    return 0;
  case 2:
    //3->31 8->81
    if (fl[0].Strong()) {
      for (short int i=1;i<3;i++) {
	if (fl[i].Strong()) {
	  for (short int j=0;j<2;j++) m_colors[i][j] = m_colors[0][j];
	  return 0;
	}
      }
    }
    // 1->33 1->88
    if (col1==0 && col2==0) {
      if ((fl[1].IsQuark()||fl[1].IsSquark()) && (fl[2].IsQuark()||fl[2].IsSquark())) {
	if (fl[1].IsAnti() && !(fl[2].IsAnti())) {
	  m_colors[1][1] = m_colors[2][0] = ATOOLS::Flow::Counter();
	  return 0;
	}
	if (fl[2].IsAnti() && !(fl[1].IsAnti())) {
	  m_colors[1][0] = m_colors[2][1] = ATOOLS::Flow::Counter();
	  return 0;
	}
      }
      if ((fl[1].IsGluon()||fl[1].IsGluino()) && (fl[2].IsGluon()||fl[2].IsGluino())) {
	m_colors[1][1] = m_colors[2][0] = ATOOLS::Flow::Counter();
	m_colors[1][0] = m_colors[2][1] = ATOOLS::Flow::Counter();
	return 0;
      }
    }
  case 3:
  default :
    msg_Error()<<"Error in Cluster_Algorithm::SetDecayColours:"<<std::endl
	       <<"   Cannot handle single color in 1 -> 2 process :"
               <<"   "<<fl[0]<<" "<<fl[1]<<" "<<fl[2]<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
}

int Cluster_Algorithm::SetColours(const Vec4D_Vector& p,Flavour * fl)
{
  // *** 2 -> 2 processes with unambiguous coulor structure
  // (a) no colors
  // (b) two (s)quarks
  // (c) two (s)quarks and one gluon/gluino
  // (d) two gluons (ADD-Model) 
  // (e) three gluons (ADD-Model)

  int ncol   = 0;
  int nquark = 0;
  int ngluon = 0;
  for (int i=0; i<4; ++i) {
    if (fl[i].Strong()) {
      ++ncol;
      if (fl[i].IsQuark() || fl[i].IsSquark()) ++nquark;
      if (fl[i].IsGluon() || fl[i].IsGluino()) ++ngluon;
    }
    m_colors[i][0]=m_colors[i][1]=0;
  }

  switch (ncol) {
  case 4:
    return Set4Colours(nquark,ngluon,p,fl);
  case 3:
    return Set3Colours(nquark,ngluon,p,fl);
  case 2:
    return Set2Colours(nquark,ngluon,p,fl);
  case 1:
    msg_Out()<<"Error in Cluster_Algorithm::SetColours() : called for 1 coloured object. \n"
	     <<"   Don't know how to handle this ! Abort the run."<<std::endl;
    for (int i=0; i<4; ++i) msg_Out()<<i<<" : "<<fl[i]<<std::endl;
    abort();
  case 0:
    return 0;
  }
  return 1;
}

int Cluster_Algorithm::Set4Colours(const int nquark,const int ngluon,
                                   const Vec4D_Vector& p,Flavour * fl)
{
  double scale;
  int prop(p_ct->IdentifyHardPropagator(scale));
  if (fl[0].IsGluon() || fl[0].IsGluino() || 
      fl[1].IsGluon() || fl[1].IsGluino() || 
      fl[2].IsGluon() || fl[2].IsGluino() || 
      fl[3].IsGluon() || fl[3].IsGluino() || prop<0) {

    msg_Out()<<METHOD<<"(): Cannot set colours for "<<std::endl;
    Combine_Table *ct(p_ct);
    while (ct->Up()!=NULL) ct=ct->Up();
    msg_Error()<<*ct<<std::endl;
    return false;
  }
  switch (prop) {
  case 1:
    if (!fl[0].IsAnti()) m_colors[0][0]=m_colors[1][1]=500;
    else m_colors[0][1]=m_colors[1][0]=500;
    if (!fl[2].IsAnti()) m_colors[2][0]=m_colors[3][1]=501;
    else m_colors[2][1]=m_colors[3][0]=501;
    break;
  case 2:
    if (!fl[0].IsAnti()) m_colors[0][0]=m_colors[2][0]=500;
    else m_colors[0][1]=m_colors[2][1]=500;
    if (!fl[1].IsAnti()) m_colors[1][0]=m_colors[3][0]=501;
    else m_colors[1][1]=m_colors[3][1]=501;
    break;
  case 3:
    if (!fl[0].IsAnti()) m_colors[0][0]=m_colors[3][0]=500;
    else m_colors[0][1]=m_colors[3][1]=500;
    if (!fl[1].IsAnti()) m_colors[1][0]=m_colors[2][0]=501;
    else m_colors[1][1]=m_colors[2][1]=501;
    break;
  }
  return true;
}

int Cluster_Algorithm::Set2Colours(const int nquark,const int ngluon,
                                   const Vec4D_Vector& p,Flavour * fl)
{
  if (nquark+ngluon>2) {
    msg_Error()<<"ERROR in Cluster_Algorithm::Set2Colours("<<nquark<<","<<ngluon<<")"<<std::endl
	       <<"   Wrong number of colours, abort."<<std::endl;
    abort();
  }
  int connected[2] = {-1,-1};
  int j(0);
  for (int i=0;i<4;i++) {
    if (!fl[i].Strong()) continue;
    if (fl[i].IsQuark() || fl[i].IsSquark()) {
      m_colors[i][0+int(fl[i].IsAnti())] = 500;
    }
    else if (fl[i].IsGluon()) {
      m_colors[i][j] = 500; m_colors[i][1-j] = 501;
    }
    connected[j++]=i;
  }    
  return 0;
}

int Cluster_Algorithm::Set3Colours(const int nquark,const int ngluon,
                                   const Vec4D_Vector& p,Flavour * fl)
{
  int connected[3] = {-1,-1,-1};
  int singlet      = -1;
  int j(0);
  if (ngluon==3) {
    for (int i=0;i<4;i++) {
      if (!fl[i].Strong()) { 
	singlet = i; 
	continue; 
      }
      if (fl[i].IsGluon() || fl[i].IsGluino()) {
	connected[j] = i;
	m_colors[i][0+(i>1)] = 500+j; 
	if (j==2) j=-1;
	m_colors[i][1-(i>1)] = 501+j;
	j++;
      }    
    }
  }
  else if (ngluon==1) {
    for (int i=0;i<4;i++) {
      if (!fl[i].Strong()) { 
	singlet = i; 
	continue; 
      }
      if (fl[i].IsQuark() || fl[i].IsSquark()) {
	connected[j] = i;
	m_colors[i][0+int(fl[i].IsAnti())] = 500+j;
	j++;
      }
    }
    bool tmode = (connected[0]<2) ^ (connected[1]<2);
    for (int i=0;i<4;i++) {
      if (fl[i].IsGluon() || fl[i].IsGluino()) {
	if (tmode) {
	  if (i<2) {
	    m_colors[i][1-int(fl[connected[1]].IsAnti())] = 500;
	    m_colors[i][0+int(fl[connected[0]].IsAnti())] = 501;
	  }
	  else {
	    m_colors[i][1-int(fl[connected[0]].IsAnti())] = 501;
	    m_colors[i][0+int(fl[connected[1]].IsAnti())] = 500;
	  }
	}
	else {
	  for (int j=0;j<2;j++) {
	    m_colors[i][j] += m_colors[connected[0]][j] + 
	      m_colors[connected[1]][j];
	  }
	}
      }
    }    
  }
  for (int i=0;i<4;i++) {
    msg_Debugging()<<METHOD<<" "<<i<<" "<<p_ct->GetLeg(i).Point()->fl<<" "
		   <<m_colors[i][0]<<" "<<m_colors[i][1]<<std::endl;
  }
  return 0;
}

ME2_Base *Cluster_Algorithm::GetXS(const ATOOLS::Flavour_Vector &fl)
{
  Flav_ME_Map::const_iterator xit(m_xsmap.find(fl));
  if (xit!=m_xsmap.end()) return xit->second;
  Process_Info pi;
  pi.m_oqcd=2;
  pi.m_oew=0;
  pi.m_ii.m_ps.push_back(Subprocess_Info(fl[0]));
  pi.m_ii.m_ps.push_back(Subprocess_Info(fl[1]));
  pi.m_fi.m_ps.push_back(Subprocess_Info(fl[2]));
  pi.m_fi.m_ps.push_back(Subprocess_Info(fl[3]));
  ME2_Base* me2=dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(pi));
  m_xsmap[fl]=me2;
  return me2;
}

int Cluster_Algorithm::SetColours(EXTRAXS::ME2_Base * xs, 
				  const Vec4D_Vector& p, ATOOLS::Flavour * fl)
{
  p_xs=xs;
  if (!p_xs) return SetColours(p,fl);
  bool test(p_xs->SetColours(p)), check(true);
  for (int i=0; i<4; ++i) {
    if (fl[i].IsQuark() || fl[i].IsSquark()) {
      if ( fl[i].IsAnti() && 
	   (p_xs->Colours()[i][0]!=0 || p_xs->Colours()[i][1]==0)) check=false;
      if (!fl[i].IsAnti() && 
	  (p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]!=0)) check=false;
    }
    if ((fl[i].IsGluon() || fl[i].IsGluino()) && 
	(p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]==0))   check=false;
    if (!check) {
      msg_Error()<<"Cluster_Algorithm::SetColours(..): \n"
		 <<"Colour check failed for the following combination:"
		 <<std::endl;
      for (int i=0; i<4; ++i) 
	msg_Error()<<"   "<<i<<" : "<<fl[i]<<" ("
		   <<p_xs->Colours()[i][0]<<","
		   <<p_xs->Colours()[i][1]<<")"<<std::endl;
      msg_Error()<<"Abort."<<std::endl;
      abort();
    }
  }
  for (int i=0;i<4;i++) {
    m_colors[i][0] = p_xs->Colours()[i][0];
    m_colors[i][1] = p_xs->Colours()[i][1];
  }
  return test;
}

void Cluster_Algorithm::Convert()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  Selector_Base *jf=p_proc->Selector()
    ->GetSelector("Jetfinder");
  Combine_Table *ct_tmp(p_ct);
  while (ct_tmp->Up()) ct_tmp=ct_tmp->Up();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(p_ms);
  p_ampl->SetJF(jf);
  p_ampl->SetNIn(p_proc->NIn());
  p_ampl->SetOrderEW(p_proc->OrderEW());
  p_ampl->SetOrderQCD(p_proc->OrderQCD());
  PHASIC::Process_Base *pb(p_proc->IsMapped()?
			   p_proc->MapProc():p_proc);
  double muf2(pb->ScaleSetter()->Scale(stp::fac));
  double mur2(pb->ScaleSetter()->Scale(stp::ren));
  for (int i(0);i<ct_tmp->NLegs();++i) {
    size_t id(ct_tmp->GetLeg(i).ID());
    Flavour flav(i<2?ct_tmp->Flav(i).Bar():ct_tmp->Flav(i));
    Vec4D mom(i<2?-ct_tmp->Momentum(i):ct_tmp->Momentum(i));
    p_ampl->CreateLeg(mom,flav,ColorID(),id);
    p_ampl->Legs().back()->SetStat(1);
  }
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  Cluster_Amplitude *eampl(p_ampl);
  while (ct_tmp->Down()) {
    int iwin, jwin, kwin;
    double mu2;
    double kt2qcd(ct_tmp->GetWinner(iwin,jwin,kwin,mu2));
    if (iwin>jwin) std::swap<int>(iwin,jwin);
    ct_tmp=ct_tmp->Down();
    const Leg &win(ct_tmp->GetLeg(iwin));
    Cluster_Amplitude *ampl(p_ampl);
    p_ampl=p_ampl->InitNext();
    p_ampl->SetMS(p_ms);
    p_ampl->SetJF(jf);
    p_ampl->SetNIn(ampl->NIn());
    ampl->SetKT2(kt2qcd);
    ampl->SetMu2(mu2);
    for (int i(0);i<ct_tmp->NLegs();++i) {
      size_t id(ampl->Leg(i<jwin?i:i+1)->Id());
      Flavour flav(i<2?ct_tmp->Flav(i).Bar():ct_tmp->Flav(i));
      Vec4D mom(i<2?-ct_tmp->Momentum(i):ct_tmp->Momentum(i));
      if (i==iwin) id+=ampl->Leg(jwin)->Id();
      p_ampl->CreateLeg(mom,flav,ColorID(),id);
      if (IdCount(id)==1) {
	p_ampl->Legs().back()->SetStat(1);
      }
      else if (i==iwin) {
	p_ampl->Legs().back()->SetK(ampl->Leg(kwin)->Id());
	ampl->SetIdNew(ct_tmp->Up()->GetLeg(jwin).ID());
	if (win.Point()->t>10) {
	  p_ampl->Legs().back()->SetStat(3);
	  SetNMax(p_ampl->Prev(),id,win.Point()->t-10);
	}
      }
    }
    p_ampl->SetMuR2(ampl->MuR2());
    p_ampl->SetMuF2(ampl->MuF2());
    p_ampl->Decays()=ct_tmp->Decays();
    p_ampl->SetOrderEW(ampl->OrderEW()-win.OrderQED());
    p_ampl->SetOrderQCD(ampl->OrderQCD()-win.OrderQCD());
    p_ampl->SetKin(win.Kin());
  }
  p_ampl->SetProcs(p_proc);
  PDF::CParam scale(p_clus->CoreScale(p_ampl));
  p_ampl->SetKT2(scale.m_kt2);
  p_ampl->SetMu2(scale.m_mu2);
  size_t nmax(p_proc->Info().m_fi.NMaxExternal());
  if (dynamic_cast<POWHEG_Process*>(p_proc->Parent())) {
    if (p_proc->Parent()->Info().m_fi.NExternal()<=nmax) ++nmax;
    else nmax=eampl->Legs().size()-eampl->NIn();
    eampl->SetNLO(1);
    while (eampl->Next()) {
      eampl=eampl->Next();
      eampl->SetNLO(1);
    }
  }
  p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
  SetNMax(p_ampl,(1<<(p_proc->NIn()+p_proc->NOut()))-1,nmax);
  for (size_t i(0);i<2;++i)
    p_ampl->Leg(i)->SetCol(ColorID(m_colors[i][1],m_colors[i][0]));
  for (size_t i(2);i<4;++i)
    p_ampl->Leg(i)->SetCol(ColorID(m_colors[i][0],m_colors[i][1]));
  while (p_ampl->Prev()) {
    Cluster_Amplitude *ampl(p_ampl->Prev());
    size_t jwin(std::numeric_limits<size_t>::max());
    for (size_t i(0);i<ampl->Legs().size();++i) {
      if (i==jwin) continue;
      Cluster_Leg *li(ampl->Leg(i));
      Cluster_Leg *lij(p_ampl->Leg(i<jwin?i:i-1));
      if (li->Id()==lij->Id()) {
	li->SetCol(lij->Col());
      }
      else {
	for (size_t k(i+1);k<ampl->Legs().size();++k) {
	  Cluster_Leg *lk(ampl->Leg(k));
	  if (lk->Id()&lij->Id()) {
	    Cluster_Amplitude::SetColours(lij,li,lk);
	    jwin=k;
	    break;
	  }
	}
      }
    }
    if (p_proc->Integrator()->InSwaped())
      std::swap<Cluster_Leg*>(p_ampl->Legs()[0],p_ampl->Legs()[1]);
    msg_Debugging()<<*p_ampl<<"\n";
    p_ampl=p_ampl->Prev();
  }
  if (p_proc->Integrator()->InSwaped())
    std::swap<Cluster_Leg*>(p_ampl->Legs()[0],p_ampl->Legs()[1]);
  msg_Debugging()<<*p_ampl<<"\n";
  msg_Debugging()<<"}\n";
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
