#include "Amplitude_Group.H"
#include "Super_Amplitude.H"
#include "Message.H"
#include "Zfunc.H"
#include "Basic_Sfuncs.H"

#include "MyTiming.H"
#include "prof.hh"

#include <algorithm>

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;


// compare to graph_families
class Compare_Graph_Families {
public:
  int operator()(const Graph_Family * a, const Graph_Family * b) {
    if (a->znumber<b->znumber) return 1;
    else if (a->znumber>b->znumber) return 0;
    else {
      if (a->permnumber<b->permnumber) return 1;
      else if (a->permnumber>b->permnumber) return 0;
      else {
	if (a->topnumber<b->topnumber) return 1;
	else return 0;
      }
    }
  }
};

Graph_Family *  Amplitude_Group::FindFamily(int zn, int tn, int pn) {
  for (Graph_Families::iterator git=family_table.begin(); git!=family_table.end();++git) {
    if (((*git)->znumber==zn) && ((*git)->topnumber==tn) && ((*git)->permnumber==pn)) {
      return (*git);
    }
  }
  return 0;
}


void Amplitude_Group::BuildGlobalString(int* _b,int _n,
					Basic_Sfuncs* _BS,
					Flavour* _fl,
					String_Handler* _shand)
{
  //  PROFILE_HERE;
  // fill map
  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {    
    Zfunc_List* zl = (*g)->GetZlist();
    for (Zfunc_Iterator zit=zl->begin();zit!=zl->end();++zit) {
      graph_table[(*zit)->m_str].push_back((*g));
    }
  }

  // fill family table
  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {    
    int zn = (*g)->GetZlist()->size();
    int tn = (*g)->topnum;
    int pn = (*g)->permnum;
    Graph_Family * gf = FindFamily(zn,tn,pn);
    if (!gf) {
      gf = new Graph_Family;
      gf->was_clustered = 0;
      gf->znumber       = zn;
      gf->topnumber     = tn;
      gf->permnumber    = pn;
      family_table.push_back(gf);
    }
    else {
      gf->banner += string("+");
    }
    gf->graphs.push_back(*g);

    Zfunc_List * zl = (*g)->GetZlist();
    for (Zfunc_Iterator zit=zl->begin();zit!=zl->end();++zit) {
      if (zit!=zl->begin()) gf->banner += string("*");
      gf->banner += (*zit)->m_str;
    }

  }

  // sort family table acording to znumber, permnumber, and topnumber,
  sort(family_table.begin(), family_table.end(),Compare_Graph_Families());

  // family output
  cout<<"families="<<family_table.size()<<endl;
  for (Graph_Families::iterator git=family_table.begin(); git!=family_table.end();++git) {
    cout<<" "<<(*git)->znumber<<" "<<(*git)->topnumber<<" "<<(*git)->permnumber<<" - "
	<<(*git)->graphs.size()<<endl;
    //    cout<<"  : "<<(*git)->banner<<endl;
  } 

  // count zn
  int zncount=0, last_zn=0;
  for (Graph_Families::iterator git=family_table.begin(); git!=family_table.end();++git) {
    if ((*git)->znumber!=last_zn) {
      last_zn=(*git)->znumber;
      ++zncount;
    }
  }
  cout<<" zncount="<<zncount<<endl;

  //using the Kabbala value of the Zfunc
  String_Tree st;
  string globalstr;  
  list<sknot*> addend_list;


  int combine_step=0;

  MyTiming watch;
  watch.Start();


  for (;;) {
    // basic output
    int nsuper=0;
    int nampl=0;
    //   nfamily   nsuper

    // cluster each family
    for (Graph_Families::iterator git=family_table.begin(); git!=family_table.end();++git) {

      if ((*git)->banner==string("")) {
	cerr<<" Empty String for family = "
	    <<(*git)->znumber<<" "
	    <<(*git)->topnumber<<" "
	    <<(*git)->permnumber<<endl;
      }
      else {
	if ( (*git)->graphs.size()>1 ) {
	  if ((*git)->was_clustered==0 ) {// has not been combined in last step, no need to try it again!!!
	    sknot* sh = st.String2Tree((*git)->banner);
	    st.DeleteMinus(sh);
	    st.Cluster(sh,0,1);
	    st.DeleteMinus(sh);
	    st.Delete(sh,string("Z[0]"));
	    (*git)->banner=st.Tree2String(sh,0);
	    list<sknot*> local_addend_list;
	    st.Addends(sh,local_addend_list);
	    (*git)->was_clustered=local_addend_list.size();
	  }
	  nsuper+=(*git)->was_clustered;
	  nampl+=(*git)->graphs.size();
	}
	else {
	  ++nsuper;
	  ++nampl;
	}
      }
    }

    watch.PrintTime();
    msg.Out()<<"  nsuper="<<nsuper<<"  ("<<nampl<<")  families="<<family_table.size()<<endl;

    // if only one family cluster finished
    if (family_table.size()==zncount) break;

    // combine families
    if (combine_step==0) {
      // combine all perms
      Graph_Families::iterator agit=family_table.begin();
      Graph_Families::iterator bgit=family_table.begin();      
      for (; agit!=family_table.end();) {
	if (agit!=bgit) {
	  if ( (*agit)->znumber==(*bgit)->znumber &&
	     (*agit)->permnumber==(*bgit)->permnumber ) {
	    //combine a to b
	    (*bgit)->was_clustered=0;
	    for (Amplitude_List::iterator g=(*agit)->graphs.begin();g!=(*agit)->graphs.end();++g) {
	      (*bgit)->graphs.push_back((*g));
	    }
	    (*bgit)->banner+=string("+")+(*agit)->banner;
	    agit=family_table.erase(agit);

	  }
	  else {
	    bgit=agit;
	    ++agit;
	  }
	}
	else 
	  ++agit;    
      }
    }
    else {
      // combine neighbors
      Graph_Families::iterator agit=family_table.begin();
      Graph_Families::iterator bgit=family_table.begin();      
      for (; agit!=family_table.end();) {
	if (agit!=bgit) {
	  if ( (*agit)->znumber==(*bgit)->znumber) {
	    //combine a to b
	    (*bgit)->was_clustered=0;
	    for (Amplitude_List::iterator g=(*agit)->graphs.begin();g!=(*agit)->graphs.end();++g) {
	      (*bgit)->graphs.push_back((*g));
	    }
	    (*bgit)->banner+=string("+")+(*agit)->banner;
	    agit=family_table.erase(agit);
	    bgit=agit;
	    if (agit!=family_table.end())
	      ++agit;
	  }
	  else {
	    bgit=agit;
	    ++agit;
	  }
	}
	else 
	  ++agit;    
      }

    }
    ++combine_step;

  }

  globalstr=family_table.front()->banner;
  if (family_table.size()>1) {
    for (int i=1;i<family_table.size();++i)
      globalstr+=string("+") + family_table[i]->banner;
  }

//   cout<<"Global Nachher: "<<globalstr<<endl;
  sknot* sh = st.String2Tree(globalstr);
  st.Addends(sh,addend_list);
  cout<<" nsuper="<<addend_list.size()<<endl;

   watch.Stop();
  watch.PrintTime();


  for (list<sknot*>::iterator it=addend_list.begin();it!=addend_list.end();++it) {
    string newaddend = st.Tree2String(*it,0); 
    st.Expand(*it);
    //msg.Tracking()<<"Expanded  : "<<st.Tree2String(*it,0)<<endl;
    
    list<sknot*> superlist;
    st.Addends(*it,superlist);
    
    if (superlist.size()>1) {
      Super_Amplitude* sg = new Super_Amplitude(_b,_n,_BS,_fl,_shand);

      for (list<sknot*>::iterator sit=superlist.begin();sit!=superlist.end();++sit) {
	//msg.Tracking()<<"   Single Amplitude: "<<st.Tree2String(*sit,0)<<endl;

	list<sknot*> zfunclist;
	st.Factors(*sit,zfunclist);
	
	sg->Add(GetSingleGraph(zfunclist));
      }  
      sg->Init(newaddend); 
      graphs.push_back(sg);
    }
  }

  //abort();
}

Amplitude_Base* Amplitude_Group::GetSingleGraph(list<sknot*>& zfunclist)
{
  String_Tree st;
  Amplitude_List & al=graph_table[st.Tree2String(zfunclist.front(),0)];

  for (vector<Amplitude_Base*>::iterator g=al.begin();g!=al.end();++g) {    
    int hit = 1;
    for (list<sknot*>::iterator it=zfunclist.begin();it!=zfunclist.end();++it) {      
      Zfunc_List* zl = (*g)->GetZlist();
      if (zl==0) {
	//SuperGraph....
	hit = 0;
	break;
      }
      int hit2 = 0;
      for (Zfunc_Iterator zit=zl->begin();zit!=zl->end();++zit) {
	if ((*zit)->m_str==st.Tree2String(*it,0)) {
	  hit2 = 1;
	  break;
	}
      }
      if (hit2==0) {
	hit = 0;
	break;
      }
    }
    if (hit) {
      Amplitude_Base* g2 = *g;
      al.erase(g);

      for (vector<Amplitude_Base*>::iterator git=graphs.begin();git!=graphs.end();++git) {          
	if (g2==(*git)) {
	  graphs.erase(git);
	  break;
	}
      }

      return g2;
    }
  }

  cerr<<"Error: No Amplitude found in Amplitude_Group::GetSingleGraph()!"<<endl;
  abort();
  
  return 0;
}
