//bof
//Version: 2 ADICIC++-0.0/2004/08/10

//Implementation of Chain_Handler.H.



//#include "Random.H"
#include "Poincare.H"
#include "Chain_Handler.H"
#include "Dipole_Parameter.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



//ostream& ADICIC::operator<<(ostream& ost, const ADICIC::...&) {
//}



//=============================================================================



//So far there is no static Chain_Handler.
int Chain_Handler::s_count=0;
const int& Chain_Handler::InStore=Chain_Handler::s_count;

int Chain_Handler::s_param=Chain_Handler::AdjustParameters();



//=============================================================================



template
const bool Chain_Handler::FindDecDipole<Chain_Evolution_Strategy::Unknown>();
//template const bool
//Chain_Handler::FindDecDipole<Chain_Evolution_Strategy::Production>();
template
const bool Chain_Handler::FindDecDipole<Chain_Evolution_Strategy::Emission>();
template
const bool Chain_Handler::FindDecDipole<Chain_Evolution_Strategy::Mass>();



//=============================================================================



Chain_Handler::Chain_Handler()
  : f_below(false), p_cix(NULL),
    m_k2tcomp(0.0), p_cha(NULL),
    m_dh1(), m_dh2(),
    p_dhwait(&m_dh1), p_dhaciv(&m_dh2), p_dhtemp(NULL),
    i_fix(NULL), i_run(NULL) {

  ++s_count;

  PresetCompScale();

  fp_finddecdip[0]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Unknown>;
  fp_finddecdip[1]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Production>;
  fp_finddecdip[2]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Emission>;
  fp_finddecdip[3]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Mass>;

}





Chain_Handler::Chain_Handler(Chain& cha)
  :  f_below(false), p_cix(NULL),
     m_k2tcomp(0.0), p_cha(NULL),
     m_dh1(), m_dh2(),
     p_dhwait(&m_dh1), p_dhaciv(&m_dh2), p_dhtemp(NULL),
     i_fix(NULL), i_run(NULL) {

  ++s_count;

  PresetCompScale();

  fp_finddecdip[0]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Unknown>;
  fp_finddecdip[1]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Production>;
  fp_finddecdip[2]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Emission>;
  fp_finddecdip[3]=&ADICIC::Chain_Handler::
    FindDecDipole<Chain_Evolution_Strategy::Mass>;

  if(cha|*this) {
    if(cha.IsHandledBy(*this)); else {
      cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
      assert(cha.IsHandledBy(*this));
    }
  }
  else {
    cerr<<"\nMethod: ADICIC::Chain_Handler::Chain_Handler(ADICIC::Chain&): "
	<<"Warning: Attaching Chain failed!\n"<<endl;
  }

}





Chain_Handler::~Chain_Handler() {

  --s_count;

  if(p_cix) delete p_cix;

  if(!p_cha) return;
  if(p_cha->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
    assert(p_cha->IsHandledBy(*this));
  }

  *p_cha|0;

}



//=============================================================================



void Chain_Handler::ShowParameters() {    //Static.
  cout<<endl;
  cout<<"======================================"<<endl;
  cout<<"    Valid Chain_Handler parameters"<<endl;
  cout<<"--------------------------------------"<<endl;
  cout<<"Chain evolution strategy is set to "<<s_param<<".\n";
  cout<<"======================================"<<endl;
}





const int Chain_Handler::AdjustParameters() {    //Static.
  static bool firsttime=true;
  if(firsttime) {
    firsttime=false;
    Dipole_Parameter::ForceFirstInit();
#ifdef CHAIN_HANDLER_OUTPUT
    cout<<"ADICIC::Chain_Handler: initialize parameter for the 1st time.\n";
#endif
  }
  if(0 < Dipole_Parameter::ChainEvolutionStrategy() &&
     Dipole_Parameter::ChainEvolutionStrategy() < 4)
    s_param=Dipole_Parameter::ChainEvolutionStrategy();
  else
    s_param=0;
  return s_param;
}



//=============================================================================



void Chain_Handler::Reset() {

  f_below=false;
  if(p_cix) { delete p_cix; p_cix=NULL;}

  PresetCompScale();

  p_dhwait=&m_dh1;
  p_dhaciv=&m_dh2;
  p_dhtemp=NULL;

  i_fix=NULL;
  i_run=NULL;

  m_dh1.ResetStatus();
  m_dh2.ResetStatus();

  if(m_dh1.IsDocked() || m_dh2.IsDocked()) {
    assert(p_cha);
    for(list<Dipole*>::const_iterator it=p_cha->DipolePointerList().begin();
	it!=p_cha->DipolePointerList().end(); ++it)
      **it|0;
  }

}





const bool Chain_Handler::EvolveChainByOneStep() {

  if(!p_cha) return false;
  if(p_cha->IsEmpty()) return false;

  assert(!m_dh1.IsDocked() && !m_dh2.IsDocked());

  PresetCompScale();

  if(p_cha->Status() && p_cha->ChainType()!=Chain::incorrect &&
     p_cha->LastScale()>m_k2tcomp);
  else return false;

  f_below=false;
  if(p_cix) { delete p_cix; p_cix=NULL;}

  //No testing of global parameters.
  //assert( p_dip->InvMass() > Sudakov_Calculator::MinOfK2t() );

  if(this->FindDipole()) {
    this->SplitDipole();
    double tem=m_k2tcomp;
    p_cha->SetLastScale()=tem;
    return true;
  }

#ifdef CHAIN_HANDLER_OUTPUT
#endif

  return false;

}





const bool Chain_Handler::EvolveChain() {
  return false;
}



//=============================================================================



template<class _Strategy>
const bool Chain_Handler::FindDecDipole() {
  cerr<<"\nMethod: const bool ADICIC::Chain_Handler::FindDipole(): "
      <<"Warning: Method has not been specified!\n"<<endl;
  return false;
}





template<> const bool
Chain_Handler::FindDecDipole<Chain_Evolution_Strategy::Production>() {

  static const bool noaffirm=ProductionStrategyInfo();    //Gives warning.

  i_fix=NULL;
  i_run=p_cha->DipolePointerList().begin();

  for(; i_run!=p_cha->DipolePointerList().end(); ++i_run) {

    Dipole& dip=**i_run;
    assert(dip|*p_dhaciv);

    if(p_dhaciv->InduceGluonEmission()==noaffirm) {    //Removes the warning.
      dip|0; continue;
    }

    if(dip.EmitScale() >= p_cha->LastScale()) {
      bool result;
      do {
	dip.SetBootScale()=dip.EmitScale();
	result=p_dhaciv->InduceGluonEmission();
	if(result==false) break;
      }
      while(dip.EmitScale() >= p_cha->LastScale());
      if(result==false) {
	dip.SetBootScale()=dip.ProdScale();
	dip|0; continue;
      }
    }

    if(dip.EmitScale() < m_k2tcomp) {
      dip.SetBootScale()=dip.ProdScale();
      dip|0; continue;
    }

    m_k2tcomp=dip.EmitScale();
    p_dhtemp=p_dhwait;
    p_dhwait=p_dhaciv;
    p_dhaciv=p_dhtemp;
    if(p_dhtemp->IsDocked()) {
      Dipole& loser=**i_fix;
      loser.SetBootScale()=loser.ProdScale();
      loser|0;
    }
    i_fix=i_run;

  }

  return bool(i_fix!=NULL);

}



//=============================================================================



const bool Chain_Handler::SplitDipole() {

  Dipole& winner=**i_fix;

  Vec4D momm=winner.GetTopBranchPointer()->Momentum();
  momm+=winner.GetBotBranchPointer()->Momentum();
  p_cha->UpdateMomentum(-1.0,momm);

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<momm<<"\t"<<momm.Abs2()<<endl;
  cout<<p_cha->Momentum()<<"\t"
      <<p_cha->Momentum().Abs2()<<"\t"
      <<p_cha->InvMass()<<endl;
#endif

  assert(p_dhwait->FinishGluonEmission());

  bool below;
  Dipole* pnewdip=NULL;
  Dipole::Glubranch* pglu=NULL;

  p_dhwait->DecoupleNewDipole(pnewdip,below);
  assert(pnewdip);
  p_dhwait->DecoupleGlubranch(pglu);
  assert(pglu);

  winner|0;

  momm=pglu->Momentum();
  p_cha->GlubranchPointerList().push_back(pglu);

  if(below) {
    i_run=i_fix;
    ++i_run;
    i_run=p_cha->DipolePointerList().insert(i_run,pnewdip);
    momm+=winner.GetTopBranchPointer()->Momentum();
    momm+=pnewdip->GetBotBranchPointer()->Momentum();
  } else {
    i_run=p_cha->DipolePointerList().insert(i_fix,pnewdip);
    ++i_run;
    i_fix=i_run;
    momm+=pnewdip->GetTopBranchPointer()->Momentum();
    momm+=winner.GetBotBranchPointer()->Momentum();
  }

  p_cha->UpdateMomentum(1.0,momm);

  return true;

}



//=============================================================================



//=============================================================================





//eof
