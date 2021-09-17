#include "CFPSHOWER++/Tools/Splitting.H"
#include "CFPSHOWER++/Tools/Configuration.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;

Configuration::Configuration(const double & t, const double & t0) :
  p_ampl(NULL), m_t(t), m_t0(t0) {}

Configuration::Configuration(Cluster_Amplitude * const ampl,
			     map<Cluster_Leg*,Parton*> & lmap,
			     Cluster_Definitions * cluster) :
  p_ampl(ampl), m_t(p_ampl->KT2()), m_t0(-1.), m_pos(Vec4D(0.,0.,0.,0.))
{
  //msg_Out()<<"************************************************************\n"
  //	   <<"************************************************************\n"
  //	   <<"************************************************************\n"
  //	   <<METHOD<<":\n"<<(*ampl)<<"\n";
  if (!p_ampl) return;
  Fill(lmap);
  if (p_ampl->Prev()) m_tveto = p_ampl->Prev()->KT2();
}

void Configuration::Fill(map<Cluster_Leg*,Parton*> & lmap) {
  SanitizeSingletGluons();
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
    Cluster_Leg * leg = p_ampl->Leg(i);
    Parton * parton   = new Parton(leg->Flav(),leg->Mom());
    parton->SetColor(Color(leg->Col().m_i,leg->Col().m_j));
    if (i<p_ampl->NIn()) {
      parton->SetBeam(leg->Mom()[3]>0.?2:1);
      parton->SetXB();
    }
    push_back(parton);
    lmap[leg] = parton;
  }
  for (Parton_List::iterator pit1=begin();pit1!=end();pit1++) {
    if ((*pit1)->Beam()!=0) continue;
    if ((*pit1)->Flav().IsGluon() &&
	(*pit1)->GetColor()[0]==(*pit1)->GetColor()[1]) {
      msg_Out()<<"Gotcha: "
	       <<"["<<(*pit1)->GetColor()[0]<<" "<<(*pit1)->GetColor()[1]<<"], "
	       <<(*pit1)->Mom()<<"\n";
      exit(1);
    }
  }
  EstablishRelations();
}

void Configuration::SanitizeSingletGluons() {
  if (p_ampl->Legs().size()!=4) return;
  Cluster_Leg * leg2 = p_ampl->Leg(2), * leg3 = p_ampl->Leg(3);  
  if (leg2->Flav().IsGluon() && leg2->Col().Singlet()) {
    if (leg3->Flav().IsGluon() && leg3->Col().Singlet()) {
      if (leg2->Col().m_i!=leg3->Col().m_i) {
	leg2->SetCol(ATOOLS::ColorID(leg2->Col().m_i,leg3->Col().m_i));
      }
      else {
	leg2->SetCol(ATOOLS::ColorID(leg2->Col().m_i,leg2->Col().m_i+1));
      }
      leg3->SetCol(leg2->Col().Conj());
      //msg_Out()<<METHOD<<" fixed colour:\n"<<(*p_ampl)<<"\n";
    }
    else {
      msg_Error()<<METHOD<<" found a singlet gluon in 2->2 amplitude,\n"
		 <<"  but no partner gluon to compensate colours properly.\n"
		 <<"  Will continue and hope for the best.\n";
    }
  }
}

size_t Configuration::NPartons() {
  size_t npartons = 0;
  for (Parton_List::const_iterator pit=begin();pit!=end();pit++) {
    if ((*pit)->On() && (*pit)->Flav().Strong()) npartons++;
  }
  return npartons;
}

void Configuration::EstablishRelations() {
  Parton_List::iterator pit1=begin(), pit2;
  while (pit1!=end()) {
    pit2 = pit1; pit2++;
    while (pit2!=end()) {
      for (size_t i=0;i<2;i++) {
	if ((*pit1)->GetColor()[i]!=0 &&
	    (*pit1)->GetColor()[i]==(*pit2)->GetColor()[1-i]) {
	  (*pit1)->AddSpectator((*pit2));
	  (*pit2)->AddSpectator((*pit1));
	}
      }
      pit2++;
    }
    pit1++;
  }
}

void Configuration::EstablishHistories(map<Cluster_Leg*,Parton*> & lmap,
				       Cluster_Definitions * cluster) {
  // "inner" is the Cluster_Amplitude before one leg splitting results in p_ampl.
  // The name of the game is to reconstruct this splitting and translate it into
  // the language of the parton shower, to fix the starting conditions for the
  // configuration emerging from p_ampl.  If there is no "inner" Cluster_Amplitude,
  // p_ampl is the hardest one, and the starting conditions are given by the
  // scale mu_Q set externally.
  Cluster_Amplitude * inner = p_ampl->Next();
  if (inner) {
    // lij denotes the splitting leg, ic, jc, and kc are the id numbers for the
    // two offsprings and the spectator:  ij+k -> i+j+k
    int ic=-1,jc=-1,kc=-1;
    Cluster_Leg * splitter(NULL);
    Cluster_Leg * newleg(p_ampl->IdLeg(p_ampl->IdNew()));
    // Look for splitter in inner - the parton where a spectator K has been assigned.
    for (size_t i=0; i<inner->Legs().size()&&splitter==NULL; ++i) {
      if (inner->Leg(i)->K()) splitter = inner->Leg(i);
    }
    if (splitter==NULL) THROW(fatal_error,"Invalid PS tree");
    // Iterate through the legs of the "harder" amplitude ("inner") and the emerging
    // amplitude ("p_ampl").  inner should provide a leg splitting into one of the
    // p_ampl legs, allowing us to successively reconstruct ic, jc, and kc.
    for (size_t i=0;i<inner->Legs().size();++i) {
      Cluster_Leg * inleg(inner->Leg(i));
      Parton      * inpart = lmap[inleg];
      for (size_t j=0;j<p_ampl->Legs().size();++j) {
	Cluster_Leg * outleg(p_ampl->Leg(j));
	// If the out leg comes from the in leg, they have a common
	// bit in their Id - this is when we start reconstructing the splitting.
	// This is somewhat tricky, as there are two cases:
	// - inleg==splitter:
	//   then outleg is either the newleg or not, and, correspondingly, outleg is
	//   either j or or i.
	// - outleg = spectator (splitter->K()):
	//   then outleg is k.
	// In both cases we establish that inpart is the incoming parton to outpart,
	// while outpart is the incoming parton to inpart.
	if (inleg->Id() & outleg->Id()) {
	  // Is the out leg the spectator?
	  if (outleg->Id()==splitter->K()) kc=j;
	  Parton * outpart = lmap[outleg];
	  // Establish the relations between the splitter and out partons
	  if (inpart->Out(0)) {
	    if (inpart->Out(1)) THROW(fatal_error,"Invalid PS tree");
	    if (inleg==splitter) (outleg==newleg?jc:ic)=j;
	    inpart->SetOut(1,outpart);
	    outpart->SetIn(inpart);
	  }
	  else {
	    if (inleg==splitter) (outleg==newleg?jc:ic)=j;
	    inpart->SetOut(0,outpart);
	    outpart->SetIn(inpart);
	  }
	}
      }
    }
    if (ic<0 || jc<0 || kc<0) THROW(fatal_error,"Invalid PS tree");
    // Now we are in the position to reconstruct the actual splitting
    double ws, mu2;
    int flip(jc<ic), swap(jc<p_ampl->NIn() && flip);
    if (swap) std::swap<int>(ic,jc);

    kernel_type::code type = GetCode(ic<p_ampl->NIn(),kc<p_ampl->NIn());
    Splitting split  = cluster->KT2(*p_ampl,ic,jc,kc,splitter->Flav(),inner->Kin(),
				    type,1|(swap?2:0)|(inner->NLO()?16<<2:0),ws,mu2);
    split.SetSpectator(lmap[inner->IdLeg(splitter->K())]);
    split.SetSplitter(lmap[splitter]);
  }
}


Configuration::~Configuration() {
  while (!empty()) { if (front()) delete (front()); pop_front(); }
}

std::ostream &CFPSHOWER::operator<<(std::ostream & s,const Configuration & config) {
  s<<"Configuration with "<<config.size()<<" partons for evolution "
   <<"from "<<config.T()<<" to "<<config.T0()<<":\n";
  for (Parton_List::const_iterator pit=config.begin();pit!=config.end();pit++) s<<(**pit);
  return s;
}
