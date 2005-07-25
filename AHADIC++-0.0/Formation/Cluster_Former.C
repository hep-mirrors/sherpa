#include "Cluster_Former.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Cluster_Former::Cluster_Former() { }

Cluster_Former::~Cluster_Former() { }

void Cluster_Former::ConstructClusters(Part_List * plin, Cluster_List * clout)
{
  if (clout==NULL) {
    msg.Error()<<"ERROR in Cluster_Former::FormClusters : "<<std::endl
	       <<"   Funny Cluster_List, abort the run."<<std::endl;
    abort();
  }

  int col1, lead;
  Cluster * cluster=NULL;
  
  for (Part_Iterator pit1=plin->begin();pit1!=plin->end();++pit1) {
    col1 = (*pit1)->GetFlow(1);
    if (col1!=0) {
      if ((*pit1)->GetFlow(2)!=0) {
	msg.Error()<<"ERROR in Cluster_Former::FormClusters : "<<std::endl
		   <<"   Colour octet left in particle list."<<std::endl
		   <<"   "<<(*pit1)->Number()<<" : "<<(*pit1)->Flav()
		   <<" ("<<(*pit1)->GetFlow(1)<<","<<(*pit1)->GetFlow(2)
		   <<"), abort the run."<<std::endl;
	abort();
      }
      for (Part_Iterator pit2=plin->begin();pit2!=plin->end();++pit2) {
	if (int((*pit2)->GetFlow(2))==col1) {
	  if (int((*pit2)->GetFlow(1))!=0) {
	    msg.Error()<<"ERROR in Cluster_Former::FormClusters : "<<std::endl
		       <<"   Colour octet left in particle list."<<std::endl
		       <<"   "<<(*pit2)->Number()<<" : "<<(*pit2)->Flav()
		       <<" ("<<(*pit2)->GetFlow(1)<<","<<(*pit2)->GetFlow(2)
		       <<"), abort the run."<<std::endl;
	    abort();
	  }
	  cluster = new Cluster((*pit1)->Flav(),(*pit1)->Momentum(),(*pit2)->Flav(),(*pit2)->Momentum());
	  lead = 0;
	  if ((*pit1)->Info()=='L') lead+=1; 
	  if ((*pit2)->Info()=='L') lead+=2;
	  if (lead>0) cluster->SetLeads(ltp::code(lead));
	  clout->push_back(cluster);
	  (*pit1)->SetStatus(2);
	  (*pit2)->SetStatus(2);
	}
      }
    }
  }
}

