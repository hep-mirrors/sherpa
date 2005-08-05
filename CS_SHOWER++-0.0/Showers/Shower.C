#include "Shower.H"
#include "Parton.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

Shower::Shower() : p_actual(NULL), p_all(NULL) {}

Shower::~Shower() 
{
  if (p_actual) { delete p_actual; p_actual=NULL;}
  if (p_all)    { delete p_all;    p_all=NULL;}
}

bool Shower::EvolveShower(All_Singlets * all)
{
  p_all  = all;
  ASiter asiter = p_all->begin();
  do {
    if (!EvolveSinglet((*asiter))) return false;
    asiter++;
  } while (asiter!=p_all->end());
  return true;
}

bool Shower::EvolveSinglet(Singlet * act)
{
  p_actual        = act;
  bool mustshower = true;
  PLiter   splitter;
  Parton * newpB, * newpC;
  Vec4D mom;
  double kt2win, last(8500.);
  int mustsplit;
  while (mustshower) {
    //cout<<"==========================================================="<<endl<<(*p_actual);
    kt2win = 0.;
    splitter = SelectSplitting(kt2win);
    //cout<<"-----------------------------------------------------------"
    //	<<"Winner is "<<(*splitter)<<" "<<(*splitter)->KtTest()<<" "<<kt2win<<" -> "<<(*splitter)->GetSpect()<<endl;
    if (splitter==p_actual->end()) mustshower = false;
    else {
      if (last<kt2win) {
	cout<<"Last < kt2win : Winning splitting with "<<kt2win<<" vs. "<<last<<endl
	    <<"Split the parton : "<<endl
	    <<"Split : "<<(**splitter)<<"Spect : "<<(*(*splitter)->GetSpect())<<endl
	    <<"z = "<<(*splitter)->ZTest()<<", y = "<<(*splitter)->YTest()<<", kt = "<<(*splitter)->KtTest()<<endl
	    <<"-------------------------------"<<endl;
      }
      //cout<<"Before MakeKinematics : "<<(*splitter)->KtTest()<<", "<<(*splitter)->ZTest()<<", "<<(*splitter)->YTest()<<endl;
      newpC     = m_kin.MakeKinematics((*splitter), m_flavC);
      mom       = (*splitter)->Momentum();
      newpB     = new Parton(m_flavB,mom,(*splitter)->GetType());
      //cout<<"Split the parton : "<<endl<<(*newpB)<<(*newpC)<<(*(*splitter)->GetSpect())<<endl
      //	  <<"==============================="<<endl;
      mustsplit = p_actual->SplitParton(splitter,newpB,newpC);
      newpB->SetStart(kt2win);
      newpC->SetStart(kt2win);
      last = kt2win;
      for (PLiter plit=p_actual->begin();plit!=p_actual->end();plit++) (*plit)->SetVeto(kt2win);
      if (mustsplit) {
	//cout<<"#########################################################"<<endl<<(*p_all)<<endl;
	splitter--;
	p_all->push_back(act->SplitList(splitter));
	//cout<<(*p_all)<<"#########################################################"<<endl;
      }
    }
    //cout<<"-----------------------------------------------------------"<<endl<<(*p_actual);
  }
  return true;
}

PLiter Shower::SelectSplitting(double & kt2win) {
  PLiter winner = p_actual->end(), specter;
  Parton * split, * spect;
  for (PLiter splitter = p_actual->begin(); splitter!=p_actual->end();splitter++) {
    split   = (*splitter);
    if (split->GetCPartner()==1) {
      specter = splitter;
      specter++;
      spect   = (*specter);
      if (m_sudakov.Dice(split,spect) && split->KtTest()>kt2win) {
	kt2win  = split->KtTest();
	winner  = splitter;
	m_flavB = m_sudakov.GetFlavourB();
	m_flavC = m_sudakov.GetFlavourC();
	split->SetSpect(spect);
      }
    } 
    else if (split->GetCPartner()==-1) {
      specter = splitter;
      specter--;
      spect   = (*specter);
      if (m_sudakov.Dice(split,spect) && split->KtTest()>kt2win) {
	kt2win  = split->KtTest();
	winner  = splitter;
	m_flavB = m_sudakov.GetFlavourB();
	m_flavC = m_sudakov.GetFlavourC();
	split->SetSpect(spect);
      }
    } 
    else if (split->GetCPartner()==0) {
      specter = splitter;
      specter++;
      spect   = (*specter);
      if (m_sudakov.Dice(split,spect) && split->KtTest()>kt2win) {
	kt2win  = split->KtTest();
	winner  = splitter;
	m_flavB = m_sudakov.GetFlavourB();
	m_flavC = m_sudakov.GetFlavourC();
	split->SetSpect(spect);
      }
      specter = splitter;
      specter--;
      spect   = (*specter);
      if (m_sudakov.Dice(split,spect) && split->KtTest()>kt2win) {
	kt2win  = split->KtTest();
	winner  = splitter;
	m_flavB = m_sudakov.GetFlavourB();
	m_flavC = m_sudakov.GetFlavourC();
	split->SetSpect(spect);
      }
    } 
  }
  //if (winner!=p_actual->end())
  //   cout<<"Winner is "<<(*winner)<<" "<<sqrt((*winner)->KtTest())<<" "<<sqrt(kt2win)<<" -> "<<spect<<endl;
  
  return winner;
}
