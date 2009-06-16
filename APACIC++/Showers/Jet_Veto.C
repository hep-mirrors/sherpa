#include "APACIC++/Showers/Jet_Veto.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <limits>

using namespace APACIC;
using namespace ATOOLS;

Jet_Veto::Jet_Veto():
  p_jf(NULL)
{
}

Jet_Veto::~Jet_Veto()
{
}

int Jet_Veto::TestISKinematics(Knot *const knot,Knot *const partner)
{
  if (p_jf==NULL) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t,jv} = "<<knot->qjv<<"\n";
  double pt23(p_jf->Qij2
	      (knot->part->Momentum(),
	       knot->left->part->Momentum(),
	       partner->part->Momentum(),
	       knot->part->Flav(),knot->left->part->Flav()));
  msg_Debugging()<<"kt = "<<sqrt(knot->pt2lcm)
 		 <<" pt_3 = "<<sqrt(pt23)<<" / "
		 <<knot->part->Momentum()<<" left = "
		 <<knot->left->part->Momentum()<<"\n";
  if (knot->part->Info()!='H') {
    if (pt23>sqr(knot->qjv)) {
      msg_Debugging()<<"jet veto\n";
      return 0;
    }
  }
  return 1;
}

int Jet_Veto::TestFSKinematics(Knot *const knot)
{
  if (p_jf==NULL) return 1;
  if (knot->left==NULL) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t,jv} = "<<knot->qjv<<"\n";
  Knot *d[2]={knot->left,knot->right};
  for (short unsigned int i(0);i<2;++i) {
    if (d[i]->left==NULL) continue;
    double pt2(p_jf->Qij2(d[i]->left->part->Momentum(),
			  d[i]->right->part->Momentum(),
			  d[1-i]->part->Momentum(),
			  d[i]->left->part->Flav(),
			  d[i]->right->part->Flav()));
    if (d[i]->left->part->Info()!='H' || 
	d[i]->right->part->Info()!='H') { 
      msg_Debugging()<<"  jv ("<<d[1-i]->kn_no<<","
		     <<d[1-i]->part->Flav()<<","<<d[1-i]->part->Info()
		     <<"),("<<d[i]->kn_no<<","<<d[i]->part->Flav()
		     <<","<<d[i]->part->Info()
		     <<")->("<<d[i]->left->kn_no
		     <<","<<d[i]->left->part->Flav()
		     <<","<<d[i]->left->part->Info()
		     <<"),("<<d[i]->right->kn_no
		     <<","<<d[i]->right->part->Flav()
		     <<","<<d[i]->right->part->Info()
		     <<"), jpt = "<<sqrt(pt2)<<"\n";
      if (d[i]->part->Flav().Strong()&&
	  d[i]->left->part->Flav().Strong()&&
	  d[i]->right->part->Flav().Strong()&&
	  pt2>=sqr(knot->qjv)) {
	msg_Debugging()<<"jet veto\n";
	return 0;
      }
    }
  }
  return 1;
}

