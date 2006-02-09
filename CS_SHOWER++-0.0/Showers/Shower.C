#include "Shower.H"
#include "Parton.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

Shower::Shower(PDF::ISR_Handler * isr) : 
  p_actual(NULL), p_all(NULL), m_sudakov(isr) 
{
}

Shower::~Shower() 
{
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
  Parton * newpB, * newpC, *spect;
  Vec4D mom;
  double kt2win, last(8500.);
  int mustsplit;
  while (mustshower) {
    cout<<"==========================================================="<<endl<<(*p_actual);
    kt2win = 0.;
    splitter = SelectSplitting(kt2win);
    
    //no shower anymore 
    if (splitter==p_actual->end()) return true; 
    else {
      if (last<kt2win) {
	cout<<"Last < kt2win : Winning splitting with "<<kt2win<<" vs. "<<last<<endl
	    <<"Split the parton : "<<endl
	    <<"Split : "<<(**splitter)<<"Spect : "<<(*(*splitter)->GetSpect())<<endl
	    <<"z = "<<(*splitter)->ZTest()<<", y = "<<(*splitter)->YTest()<<", kt = "<<(*splitter)->KtTest()<<endl
	    <<"-------------------------------"<<endl;
      }
      cout<<"Before MakeKinematics : "<<(*splitter)->KtTest()<<", "<<(*splitter)->ZTest()
	  <<", "<<(*splitter)->YTest()<<" with "
	  <<(*splitter)->GetFlavour()<<" -> "<<m_flavB<<" + "<<m_flavC
	  <<" / "<<(*splitter)->GetSpect()->GetFlavour()<<endl;
      
      Vec4D splitorig((*splitter)->Momentum()), spectorig((*splitter)->GetSpect()->Momentum());
      //the FF case 
      if ((*splitter)->GetType()==pst::FS && (*splitter)->GetSpect()->GetType()==pst::FS) {
	newpC     = m_kinFF.MakeKinematics((*splitter),m_flavC);
	mom       = (*splitter)->Momentum();
	newpB     = new Parton(m_flavB,mom,(*splitter)->GetType());
	spect     = (*splitter)->GetSpect();
	mustsplit = p_actual->SplitParton(splitter,newpB,newpC);
	newpB->SetStart(kt2win);
	newpC->SetStart(kt2win);
	last = kt2win;
      }
      else {
	//the FI case 
	if ((*splitter)->GetType()==pst::FS && (*splitter)->GetSpect()->GetType()==pst::IS) {
	  newpC     = m_kinFI.MakeKinematics((*splitter),m_flavC);
	  mom       = (*splitter)->Momentum();
	  newpB     = new Parton(m_flavB,mom,(*splitter)->GetType());
	  spect     = (*splitter)->GetSpect();
	  mustsplit = p_actual->SplitParton(splitter,newpB,newpC);
	  newpB->SetStart(kt2win);
	  newpC->SetStart(kt2win);
	  last = kt2win;
	}
	else {
	  //the IF case 
	  if ((*splitter)->GetType()==pst::IS && (*splitter)->GetSpect()->GetType()==pst::FS) {
	    newpC     = m_kinIF.MakeKinematics((*splitter),m_flavC);
	    mom       = (*splitter)->Momentum();
	    newpB     = new Parton(m_flavB,mom,(*splitter)->GetType());
	    spect     = (*splitter)->GetSpect();
	    mustsplit = p_actual->SplitParton(splitter,newpB,newpC);
	    newpB->SetStart(kt2win);
	    newpC->SetStart(kt2win);
	    last = kt2win;
	  }
	}
      }
      std::cout<<"Kinematics : "<<splitorig<<" + "<<spectorig<<" --> "<<endl
	       <<newpC->Momentum()<<newpB->Momentum()<<spect->Momentum()
	       <<"    ("<<splitorig+spectorig-newpC->Momentum()-newpB->Momentum()-spect->Momentum()<<")"<<std::endl;

      for (PLiter plit=p_actual->begin();plit!=p_actual->end();plit++) (*plit)->SetVeto(kt2win);
      if (mustsplit) {
	cout<<"#########################################################"<<endl<<(*p_all)<<endl;
	splitter--;
	p_all->push_back(act->SplitList(splitter));
	cout<<(*p_all)<<"#########################################################"<<endl;
      }
    }
    cout<<"-----------------------------------------------------------"<<endl<<(*p_actual);
  }
  return true;
}

PLiter Shower::SelectSplitting(double & kt2win) {
  PLiter winner = p_actual->end(), specter;
  Parton * split, * spect;
  std::cout<<"In SelectSplitting :"<<(*p_actual)<<std::endl;
  for (PLiter splitter = p_actual->begin(); splitter!=p_actual->end();splitter++) {
    split   = (*splitter);
    std::cout<<"Splitter : "<<(*split)<<" -> "<<split->GetCPartner()<<std::endl;
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
  if (winner!=p_actual->end())
    cout<<"Winner with "<<" "<<sqrt((*winner)->KtTest())<<" "<<sqrt(kt2win)<<" -> "
        <<(*winner)->GetSpect()->Xbj()<<(*winner)->Xbj()<<" from "<<endl
	<<"Winner : "<<(**winner)<<"Spect : "<<(*(*winner)->GetSpect())<<" "<<endl;
  return winner;
}
