#include "Amplitude_Group.H"
#include "Super_Amplitude.H"
#include "Message.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

void Amplitude_Group::BuildGlobalString(int* _b,int _n,
					Basic_Sfuncs* _BS,
					Flavour* _fl,
					String_Handler* _shand)
{
  string globalstr;

  //using the Kabbala value of the Zfunc
  
  int maxz = 0;
  int minz = 100;
  for (int f=0;f<graphs.size();f++) {
    int zsize = (graphs[f]->GetZlist())->size();
    if (zsize>maxz) maxz = zsize;
    if (zsize<minz) minz = zsize;
  }

  msg.Tracking()<<"MaxZ: "<<maxz<<endl;
  msg.Tracking()<<"MinZ: "<<minz<<endl;

  String_Tree st;
  
  list<sknot*> addend_list;

  for (int zn=minz;zn<=maxz;zn++) {
    string localstr = string("");
    msg.Tracking()<<"Znnumber: "<<zn<<endl;
    for (int f=0;f<graphs.size();f++) {
      list<Zfunc*>* zl = graphs[f]->GetZlist();
      if (zl==0) {
	msg.Error()<<"Error: Amplitude_Group::BuildGlobalString(). "<<endl;
	abort();
      }

      if (zl->size()==zn) {
	localstr += string("+");
	for (list<Zfunc*>::iterator zit=zl->begin();zit!=zl->end();++zit) {
	  if (zit!=zl->begin()) localstr += string("*");
	  localstr += (*zit)->str;
	}
      }
    }
    msg.Tracking()<<"Local String: "<<localstr<<endl;
    if (localstr==string("")) {
      msg.Tracking()<<" Empty String for zn = "<<zn<<endl;
    }
    else {
      sknot* sh = st.String2Tree(localstr);
      st.DeleteMinus(sh);
      st.Cluster(sh,0,1);
      st.DeleteMinus(sh);
      st.Delete(sh,string("Z[0]"));
      msg.Tracking()<<"Local Nachher: "<<st.Tree2String(sh,0)<<endl;
      st.Addends(sh,addend_list);
      globalstr += string("+")+st.Tree2String(sh,0);
    }
  }

  msg.Tracking()<<"Global Nachher: "<<globalstr<<endl;

  for (list<sknot*>::iterator it=addend_list.begin();it!=addend_list.end();++it) {
    msg.Tracking()<<"New Addend: "<<st.Tree2String(*it,0)<<endl;

    list<sknot*> expandlist;

    st.ExpandToDepth(*it,2,expandlist);
    if (expandlist.size()>1) {
      //msg.Tracking()<<"Expanded list"<<endl;
      *it = expandlist.front();
      for (list<sknot*>::iterator eit=++(expandlist.begin());eit!=expandlist.end();++eit) 
	it = addend_list.insert(it,*eit);
      msg.Tracking()<<"Now new Addend: "<<st.Tree2String(*it,0)<<endl;
    }

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

  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {    
    int hit = 1;
    for (list<sknot*>::iterator it=zfunclist.begin();it!=zfunclist.end();++it) {      
      list<Zfunc*>* zl = (*g)->GetZlist();
      if (zl==0) {
	//SuperGraph....
	hit = 0;
	break;
      }
      int hit2 = 0;
      for (list<Zfunc*>::iterator zit=zl->begin();zit!=zl->end();++zit) {
	if ((*zit)->str==st.Tree2String(*it,0)) {
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
      graphs.erase(g);
      return g2;
    }
  }

  cerr<<"Error: No Amplitude found in Amplitude_Group::GetSingleGraph()!"<<endl;
  abort();
  
  return 0;
}
