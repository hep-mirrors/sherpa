#include "AMEGIC++/Amplitude/CFColor.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/IO_Handler.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;


std::string CFColor::noname=string("noname");

CFColor::CFColor(int N,Single_Amplitude* first,ATOOLS::Flavour * fl,char emit,char spect,string pID,bool force)
{
  Single_Amplitude* m1 = first;

  mcount = 0;

  // on will be changed

  while (m1) {
    mcount++;
    m1->on = 1;
    m1 = m1->Next;  
  }
  
  if (pID!=noname) {
    CSC.Init();
    int ncol((int)CSC.Nc);
    if (ncol!=3) pID+="_NC"+ToString(ncol);
    std::string name=ATOOLS::rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+pID+".col";
    IO_Handler ioh;
    bool gc(ioh.SetFileNameRO(name)==0);
    if (gc&&force) {
      msg_Error()<<"Color matrix for process "<<pID<<" not found!"<<endl
		 <<" Rerun with option 'ME_LIBCHECK=1'."<<endl;
      THROW(critical_error,"Failed to load color matrix.");
    }
    if (!gc) { 
      int rmcount;

      rmcount = ioh.Input<int>("");
      ncount = ioh.Input<int>("");

      if (mcount==rmcount || (force && rmcount>0)) {
	mcount = rmcount;
	id  = ioh.ArrayInput<int>("",mcount);
	CFC = new CMatrix(ioh.MatrixInput<Complex>("",ncount,ncount),ncount);

	// generate map
	map = new int[mcount];
	int cc=0;
	for (int m=0; m<mcount; ++m) {
	  if (id[m]==mcount) {
	    map[m]=cc;
	    cc++;
	  }
	  else {
	    map[m]=map[iabs(id[m])];
	  }
	}
	return; 
      }
    }
  }
  Color_Function* cm1;
  Color_Function* cm2;

  //looking for zero color structure
  int sw1;  

  sw1 = 0;
  Single_Amplitude* m2;
  m1 = first;
  while (m1) {
    if (m1->Get_CFlist()==NULL) {
      sw1 = 1;
      break;
    }
    m1 = m1->Next;
  }
  if (sw1) {
    //no color structure
    
    // matrix
    CFC = new CMatrix(1);
    (*CFC)[0][0]=1.;

    ncount=1;

    // map table
    id = new int[mcount];
    map = new int[mcount];

    for (int i=0;i<mcount;++i) {
      id[i]=0;
      map[i]=0;
    }
    id[0]=mcount;

  }
  else {
    int prop;
  
    m1 = first;
    while (m1) {
      prop = 120;
      short int j = 0;
      for (;;) {
	int hit = 0;
	cm1 = m1->Get_CFlist();
	//looking for j
	while (cm1) {
	  for(short int i=0;i<3;i++) {
	    if ((cm1->Type()==cf::D || cm1->Type()==cf::G) && i==2) break;    
	    if (cm1->ParticleArg(i)==j) {
	      for (short int k=0;k<3;k++) {
		if ((cm1->Type()==cf::D || cm1->Type()==cf::G) && k==2) break;    
		if ((k!=i) && (cm1->ParticleArg(k)>99) && 
		    ((cm1->ParticleArg(k)<120) || (cm1->ParticleArg(k)>150)))
		  hit = cm1->ParticleArg(k);
	      }
	    }
	  } 
	  if (hit>0) break;
	  cm1 = cm1->Next();
	}
	if (hit>0) {
	  //replace hit with prop
	  cm1 = m1->Get_CFlist();
	  while (cm1) {
	    cm1->Replace(hit,prop);
	    cm1 = cm1->Next();
	  }
	  prop++;
	}
	else j++;
	if (j==N) break;
      }
      m1 = m1->Next;
    }
    
    //find all identical
    id = new int[mcount];
    for (short int i=0;i<mcount;i++) id[i] = mcount;
    m1 = first;
    int n1,n2,c1,c2;
    n1 = 0;
    
    ncount = mcount;
    
    while (m1) {
      if (m1->on) {
	cm1 = m1->Get_CFlist();
	c1 = 0;
	while (cm1) {
	  c1++;
	  cm1 = cm1->Next();
	}
	m2 = m1->Next;
	n2 = n1+1;
	while (m2) {
	  if (m2->on) { 
	    cm2 = m2->Get_CFlist();
	    c2 = 0;
	    while (cm2) {
	      c2++;
	      cm2 = cm2->Next();
	    }
	    if (c1==c2) {
	      cm1 = m1->Get_CFlist();
	      //cm1=cm2 ??
	      int hit = 1;
	      int hit2 = 0;
	      while (cm1) {
		cm2 = m2->Get_CFlist();
		hit2 = 0;
		while (cm2) {
		  hit2 = Compare(cm1,cm2);
		  if (hit2!=0) {
		    hit *= hit2;
		    break;
		  }
		  cm2 = cm2->Next();
		}
		if (hit2==0) {
		  hit = 0;
		  break;
		}
		cm1 = cm1->Next();
	      }
	      if (hit && hit2 && !(hit<0 && n1==0 )) {
		//equal
		ncount--;
		m2 ->on = 0;
		id[n2] = hit*n1;
	      }
	    }
	  }
	  n2++;
	  m2 = m2->Next;
	}
      }
      n1++;
      m1 = m1->Next;
    }

    map = new int[mcount];
    int cc=0;
    for (int m=0; m<mcount; ++m) {
      if (id[m]==mcount) {
	map[m]=cc;
	cc++;
      }
      else {
	map[m]=map[iabs(id[m])];
      }
    }

    // generate "reduced matrix"
    CFC = new CMatrix(ncount);
    Complex cffactor;
    
    m1 = first;
    c1 = 0;
    while (m1) { 
      if (m1->on) {
	m2 = m1;
	c2 = c1;
	while (m2) {
	  if (m2->on) {
	    st.Reset();
	    CharNum_Map indices;
	    indices.insert(std::make_pair('~',1)); 
	    indices.insert(std::make_pair(']',1)); 
	    indices.insert(std::make_pair('[',1)); 
	    indices.insert(std::make_pair('+',1)); 
	    indices.insert(std::make_pair('-',1)); 
	    indices.insert(std::make_pair('*',1)); 
	    char c = 'V';

	    std::string cs1=m1->CFColstring;
	    std::string cs2=m2->CFColstringC;

	    DEBUG_VAR(cs1<<" * "<<cs2);
	    for (size_t pos(cs1.find("D["));
		 pos!=string::npos;
		 pos=cs1.find("D[",pos+2)) {
	      char p1=cs1[pos+2];
	      if (cs1.find(p1)==pos+2 &&
		  cs1.find(p1,pos+3)==string::npos) {
		int anti=0;
		if (fl[int(p1-'i')].IsAnti() && int(p1-'i')>=2) anti=1;
		else if (!fl[int(p1-'i')].IsAnti() && int(p1-'i')<2) anti=1;
		if (anti) std::swap<char>(cs1[pos+2],cs1[pos+4]);
	      }
	      char p2=cs1[pos+4];
	      if (cs1.find(p1)==pos+4 &&
		  cs1.find(p1,pos+5)==string::npos) {
		int anti=0;
		if (fl[int(p2-'i')].IsAnti() && int(p2-'i')>=2) anti=1;
		else if (!fl[int(p2-'i')].IsAnti() && int(p2-'i')<2) anti=1;
		if (!anti) std::swap<char>(cs1[pos+4],cs1[pos+2]);
	      }
	    }
	    for (size_t pos(cs2.find("D["));
		 pos!=string::npos;
		 pos=cs2.find("D[",pos+2)) {
	      char p1=cs2[pos+2];
	      if (cs2.find(p1)==pos+2 &&
		  cs2.find(p1,pos+3)==string::npos) {
		int anti=0;
		if (fl[int(p1-'i')].IsAnti() && int(p1-'i')>=2) anti=1;
		else if (!fl[int(p1-'i')].IsAnti() && int(p1-'i')<2) anti=1;
		if (!anti) std::swap<char>(cs2[pos+2],cs2[pos+4]);
	      }
	      char p2=cs2[pos+4];
	      if (cs2.find(p1)==pos+4 &&
		  cs2.find(p1,pos+5)==string::npos) {
		int anti=0;
		if (fl[int(p2-'i')].IsAnti() && int(p2-'i')>=2) anti=1;
		else if (!fl[int(p2-'i')].IsAnti() && int(p2-'i')<2) anti=1;
		if (anti) std::swap<char>(cs2[pos+4],cs2[pos+2]);
	      }
	    }
	    DEBUG_VAR(cs1<<" * "<<cs2);
	    
	    sknot* s1 = st.String2Tree(cs1);
	    sknot* s2 = st.String2Tree(cs2);
	    
	    
            //begin -- subtraction term specifics
	    
	    if (emit!=spect) {
	      
	      list<sknot*>  fhelp_list;
	      st.Factors(s2,fhelp_list);
	      	      
	      /*
		char nc=(char)(4*N-9);
		char ce=emit+nc,cs=spect+nc,hc='A'+nc,cc;
		if (hc>=c||(ce<'a'&&ce>=c)||(cs<'a'&&cs>=c)) 
		cerr<<"Possible color index runout! Check color matrix!!!"<<endl;
		ce=DeliverIndex(indices,ce);
		cs=DeliverIndex(indices,cs);
		cc=DeliverIndex(indices,hc);
	      */
	      
	      char ce='$',cs='#',cc='%';
	      
	      int te=-1;
	      int ts=-1;
	      for (list<sknot*>::iterator itf=fhelp_list.begin();itf!=fhelp_list.end();++itf) {
		size_t p1=(*itf)->Str().find(emit,1);
		size_t p2=(*itf)->Str().find(spect,1);
		
		if (p1!=string::npos) { 
		  string shelp =(*itf)->Str();
		  shelp[p1]=ce;
		  (*itf)->SetString(shelp);
		  te = 0;
		  switch ((*itf)->Str()[0]) {
		  case 'T':
		    te=p1-2;
		    break;
		  case 'D': {
		    te=p1;
		    break;
		  }
		  }
		}
		if (p2!=string::npos) {
		  string shelp =(*itf)->Str();
		  shelp[p2]=cs;
		  (*itf)->SetString(shelp);
		  ts = 0;
		  switch ((*itf)->Str()[0]) {
		  case 'T':
		    ts=p2-2;
		    break;
		  case 'D': {
		    ts=p2;
		    break;
		  }
		  }
		}
	      }
	      fhelp_list.clear();
	      
	      string dpc("");
	      bool xcc=false;
	      
	      if (xcc) {if (ts+te==2) dpc=string("-");}
	      else if (ts+te==2||ts+te==6) dpc=string("-");
	      switch (te) {
	      case 0:
		if (xcc) dpc+=string("i*F[")+emit+string(",")+cc+string(",")+ce+string("]*");
		else dpc+=string("i*F[")+ce+string(",")+cc+string(",")+emit+string("]*");
		break;
	      case 2:   
		dpc+=string("T[")+cc+string(",")+emit+string(",")+ce+string("]*");
		break;
	      case 4:   
		dpc+=string("T[")+cc+string(",")+ce+string(",")+emit+string("]*");
		break;
	      default:
		Abort();
	      }
	      switch (ts) {
	      case 0:
		if (xcc) dpc+=string("i*F[")+spect+string(",")+cc+string(",")+cs+string("]");
		else dpc+=string("i*F[")+cs+string(",")+cc+string(",")+spect+string("]");
		break;
	      case 2:   
		dpc+=string("T[")+cc+string(",")+spect+string(",")+cs+string("]");
		break;
	      case 4:   
		dpc+=string("T[")+cc+string(",")+cs+string(",")+spect+string("]");
		break;
	      default:
		Abort();
	      }
	      sknot* sp = st.String2Tree(dpc);
	      sknot* s2p = s2;
	      s2 = st.newsk();
	      s2->op = '*';
	      s2->right = s2p;
	      s2->left  = sp;
	    }
	    
	    //end -- subtraction term specifics
	    
	    list<sknot*>  f1_list;
	    list<sknot*>  f2_list;
	    st.Factors(s1,f1_list);
	    st.Factors(s2,f2_list);
	    vector<string>  fstring_list;
	    string help;
	    int hit=0;
	    
	    //look for pure F products 
	    for (list<sknot*>::iterator itf=f1_list.begin();itf!=f1_list.end();++itf) {
	      help = (*itf)->Str();
	      if (help[0]=='F') fstring_list.push_back(help);
	      else hit = 1;
	    }
	    if (hit==0) {
	      for (list<sknot*>::iterator itf=f2_list.begin();itf!=f2_list.end();++itf) {
		help = (*itf)->Str();
		if (help[0]=='F') fstring_list.push_back(help);
		else hit = 1;
	      } 
	    }
	    if (hit) fstring_list.clear();
	    hit = 0;
		    
	    f1_list.clear();
	    f2_list.clear();
	    
	    string fchain;
	    Complex valuef;
	    
	    if (fstring_list.size()>0) {
	      fchain = MapFChain(fstring_list); 
	      //lookup string key in map 
	      TF_Iterator fit = f_table.find(fchain);
	      if (fit!=f_table.end()) {
		valuef = f_table[fchain];
		hit = 1;
	      }
	    }
	    //evaluate the color structure 
	    if (hit==0) {
	      //fill list of used indices
	      ExtractIndices(s1,indices);
	      ExtractIndices(s2,indices);
 	      ReplaceF(s1,indices,c);
 	      ReplaceF(s2,indices,c);

	      sknot* m = st.newsk();
	      m->op = '*';
	      m->right = s1;
	      m->left  = s2;
	      ReplaceF(m,indices,c);	    
	      st.Expand(m);
	      st.Linear(m);
	      
	      list<sknot*> addend_list;
	      st.Addends(m,addend_list);
	      
	      for (list<sknot*>::iterator it=addend_list.begin();it!=addend_list.end();++it) {
		string newaddend = st.Tree2String(*it,0); 
		ReplaceG(*it);
		ReplaceD(*it,*it);

		list<sknot*>    factor_list;
		st.Factors(*it,factor_list);
		
		vector<string>  string_list;
		string help;
		Complex total  = Complex(1.,0.);
		Complex factor = Complex(1.,0.);;
		
		int foundd = 0;
	      
		//extract string 
		for (list<sknot*>::iterator it2=factor_list.begin();it2!=factor_list.end();++it2) {
		  help = (*it2)->Str();
		  if (help[0]=='D') foundd=1;
		  if (help[0]=='T') string_list.push_back(help);
		  else {
		    factor = st.evalcolor(*it2);
		    total *= factor;
		    Kabbala* newone = new Kabbala(help,factor);
		    (*it2)->value = newone;
		    (*it2)->op = 0;
		  }
		}
		factor_list.clear(); 
		//can not deal with deltas right now
		if (foundd && string_list.size()>0) string_list.clear();
		
		Complex value;
		
		if (string_list.size()>0) {
		  newaddend = BuildTChain(string_list);
		  string_list.clear();
		  factor_list.clear();
		  
		  //lookup string key in map 
		  TF_Iterator tit = t_table.find(newaddend);
		  if (tit!=t_table.end()) value = total*t_table[newaddend];
		  else {
		    st.Sort(*it);
		    ReplaceT(*it);
		    st.Expand(*it);
		    st.Linear(*it);
		    ReplaceD(*it,*it);
		    value = st.evalcolor(*it);
		    t_table.insert(std::make_pair(newaddend,value/total));
		  }
		  Kabbala* newone = new Kabbala(newaddend,value);
		  (*it)->value = newone;
		  (*it)->op = 0;
		}
		else {
		  st.Sort(*it);
		  ReplaceT(*it);
		  st.Expand(*it);
		  st.Linear(*it);
		  ReplaceD(*it,*it);
		  value = st.evalcolor(*it);
		  Kabbala* newone = new Kabbala(newaddend,value);
		  (*it)->value = newone;
		  (*it)->op = 0;
		}
	      }
	      cffactor = st.Evaluate(m);
	      if (abs(cffactor)<rpa->gen.Accu()) cffactor = Complex(0.,0.);
	      if (fstring_list.size()>0) {
		f_table.insert(std::make_pair(fchain,cffactor));
	      }
	    }
	    else cffactor = valuef;

	    (*CFC)[map[c1]][map[c2]] = cffactor;
	    (*CFC)[map[c2]][map[c1]] = conj((*CFC)[map[c1]][map[c2]]);
	    
	    //clean up the string tree ...
	    st.CleanValues();
	    fstring_list.clear();	      
	  }
	  m2 = m2->Next;
	  c2++;
	}
      }
      m1 = m1->Next;
      c1++;
    }
  }
  
  if (pID!=noname && pID[0]!='N') Output(pID);

  // check if Matrix can be reduce even further!
  int idcc=0;
  
  int * idid = new int[ncount];
 
  for (int i=0; i<ncount; ++i) 
    idid[i]=-1;
  for (int i=0; i<ncount; ++i) {
    if (idid[i]==-1) { idid[i]=idcc; ++idcc; }
    for (int j=i+1; j<ncount; ++j) {
      int hit=1;
      Complex factor=(*CFC)[i][0]/(*CFC)[j][0];
      for (int k=0; k<ncount; ++k) {
	if ((*CFC)[i][k]!=factor*(*CFC)[j][k]) {
	  hit=0;
	  break;
	}
      }
      if (hit) {
	idid[j] =idid[i];
      }
    }
  }

  delete [] idid;

  //make m's on again
  m1 = first;
  while (m1) {   
    m1->on = 1;
    m1 = m1->Next;
  }
  
}

CFColor::~CFColor()
{
  if (CFC) delete CFC;
  if (id)  delete [] id;
  if (map) delete [] map;
}

string CFColor::BuildTChain(vector<string>  string_list) 
{
  string    key;
  Char_Map  translator;  
  char tmp,ca = 'A';
  vector<string> tmp_list;
  
  //generate the traces  
  for(;;) {
    if (tmp_list.size()==0) {
      tmp = string_list[0][2];
      tmp_list.push_back(string_list[0]);
      if (translator.insert(std::make_pair(tmp,ca)).second) ca++;
      key+=translator[tmp];
      vector<string>::iterator it = string_list.begin();
      string_list.erase(it);
    }
    for (size_t i=0;i<string_list.size();i++) {
      if (tmp_list[tmp_list.size()-1][6]==string_list[i][4]) { 
	tmp_list.push_back(string_list[i]);
	tmp = string_list[i][2];
	if (translator.insert(std::make_pair(tmp,ca)).second) ca++;
	key+=translator[tmp];
	vector<string>::iterator it = string_list.begin();
	string_list.erase(it+=i);
	i--;
      }
      if (tmp_list[0][4]==tmp_list[tmp_list.size()-1][6]) {
	key += ToString(tmp_list.size()); 
	tmp_list.clear();
	break;
      }
    }
    if (string_list.size()==0) break; 
  }
  return key;
}

string CFColor::MapFChain(vector<string> fstring_list) 
{
  string    key;
  Char_Map  translator;  
  char tmp,ca = 'A';
  vector<string> tmp_list;
  
  for (size_t i=0;i<fstring_list.size();i++) {
    for (int j=2;j<8;j+=2) {
      tmp = fstring_list[i][j];
      if (translator.insert(make_pair(tmp,ca)).second) ca++;
      key+=translator[tmp];
    }
  }
  return key;
}

void CFColor::Output(string & dirname) {
  std::string name;
  name=ATOOLS::rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+dirname+".col";
  IO_Handler ioh;
  ioh.SetFileName(name);
  ioh.Output("",mcount);          // no of ampls
  ioh.Output("",ncount);          // size of colormatrix

  ioh.ArrayOutput("",id,mcount);
  ioh.MatrixOutput("",CFC->GetMatrix(),ncount,ncount);
}

int CFColor::CompareArg(int a,int b, int c,Color_Function* cm1,Color_Function* cm2)
{
  if (cm1->ParticleArg(a)!= cm2->ParticleArg(0)) return 0;
  if (cm1->ParticleArg(b)!= cm2->ParticleArg(1)) return 0;
  if (cm1->Type()==cf::D || cm1->Type()==cf::G) return 1;
  if (cm1->ParticleArg(c)!=cm2->ParticleArg(2)) return 0;

  return 1;
}

int CFColor::Compare(Color_Function* cm1,Color_Function* cm2)
{
  if (cm1->Type()!=cm2->Type()) return 0;
  if (CompareArg(0,1,2,cm1,cm2)) return 1;
  if (cm1->Type()==cf::D) 
    if (CompareArg(1,0,2,cm1,cm2)) return 1;
  if (cm1->Type()==cf::F) {
    if (CompareArg(2,0,1,cm1,cm2)) return 1;
    if (CompareArg(1,2,0,cm1,cm2)) return 1;

    if (CompareArg(2,1,0,cm1,cm2)) return -1;
    if (CompareArg(0,2,1,cm1,cm2)) return -1;
    if (CompareArg(1,0,2,cm1,cm2)) return -1;
  }
  return 0;
}

void CFColor::ReplaceT(sknot* m)
{
  if (m==0) return;
  if (m->op=='*') {
    sknot* s1 = m->right;
    sknot* s2 = 0;
    if (m->left->op=='*') s2 = m->left->right;
    else {
      if (m->left->op==0) s2 = m->left;
    }
    if (s2!=0) {
      if (s1->Str().length()==8 && s2->Str().length()==8) { 
	if (s1->Str()[0]=='T' && s2->Str()[0]=='T' && s1->Str()[2]==s2->Str()[2]) {
	  const string* st1 = &s1->Str();
	  const string* st2 = &s2->Str();
	  
	  if ((*st1)[4]==(*st2)[6] || (*st1)[6]==(*st2)[4]) {
	    // == CF * D[]
	    char c1[2];
	    char c2[2];
	    c1[1] = 0;
	    c2[1] = 0;
	    if ((*st1)[4]==(*st2)[6]) {
	      c1[0] = (*st1)[6];
	      c2[0] = (*st2)[4];
	    }
	    else {
	      c1[0] = (*st1)[4];
	      c2[0] = (*st2)[6];
	    }
	    s1->SetString(string("CF"));
	    s2->SetString(string("D[") + string(c1) + string(",") +
			  string(c2) + string("]"));
	  }
	  else {
	    // == 1/2*(DD-1/3DD)
	    string s;
	    char c1[2];
	    char c2[2];
	    c1[1] = 0;
	    c2[1] = 0;

	    c1[0] = (*st1)[6];c2[0] = (*st2)[4];
	    s  = string("D[") + string(c1) + string(",") +
	      string(c2) + string("]*");
	    c1[0] = (*st1)[4];c2[0] = (*st2)[6];
	    s += string("D[") + string(c1) + string(",") +
	      string(c2) + string("]");
	    s2->left = st.String2Tree(s);

	    c1[0] = (*st1)[4];c2[0] = (*st1)[6];
	    s  = string("iNc*");
	    s += string("D[") + string(c1) + string(",") +
	      string(c2) + string("]*");
	    c1[0] = (*st2)[4];c2[0] = (*st2)[6];
	    s += string("D[") + string(c1) + string(",") +
	      string(c2) + string("]");
	    s2->right = st.String2Tree(s);
	    
	    s1->SetString(string("0.5"));
	    s2->SetString(string(""));
	    s2->op = '-';

	    s = string("D");
	  }
	}
      }
    }
  }
  ReplaceT(m->left);
  ReplaceT(m->right);
}

void CFColor::ReplaceD(sknot* m, sknot* start)
{
  if (m==0) return;
  if (m->op=='*') {
    sknot* s1 = m->right;
    sknot* s2 = 0;
    if (m->left->op=='*') {
      s2 = m->left->right;
    }
    else {
      if (m->left->op==0) {
	s2 = m->left;
      }
    }
    if (s2!=0) {
      if (s1->Str().length()==6) {
	if (s1->Str()[0]=='D') {
	  if (s1->Str()[2]==s1->Str()[4]) s1->SetString(string("Nc"));
	  else {
	    // kill D's
	    // replace s1->Str()[2] -> s1->Str()[4]
	    char from = s1->Str()[2];
	    char to   = s1->Str()[4];
	    SingleReplaceD(start,s1,from,to);
	    s1->SetString(string("1"));
	  }
	}
      }
      if (s2->Str().length()==6) {
	if (s2->Str()[0]=='D') {
	  if (s2->Str()[2]==s2->Str()[4]) s2->SetString(string("Nc"));
	  else {
	    // kill D's
	    // replace s2->Str()[2] -> s2->Str()[4]
	    char from = s2->Str()[2];
	    char to   = s2->Str()[4];
	    SingleReplaceD(start,s2,from,to);
	    s2->SetString(string("1"));
	  }
	}
      }
    }
  }
  if (m->op=='*') {
    ReplaceD(m->left,start);
    ReplaceD(m->right,start);
  }
  else {
    ReplaceD(m->left,m->left);
    ReplaceD(m->right,m->right);
  }
}


int CFColor::SingleReplaceD(sknot * m,sknot * orig,char from, char to)
{
  if (m==0) return 0;
  
  if (m->op==0) {
    if (m->Str().length()==6) {
      if (m->Str()[0]=='D') {
	string shelp = m->Str();
	if (from==shelp[2]) {
	  if (m!=orig) {
	    shelp[2]=to;
	    m->SetString(shelp);
	    return 1;
	  }
	}
	if (from==shelp[4]) {
	  string shelp = m->Str();
	  shelp[4]=to;
	  m->SetString(shelp);
	  return 1;
	}
      }
    }
    //new
    if (m->Str().length()==8) {
      if (m->Str()[0]=='T') {
	string shelp = m->Str();
	if (shelp[4]==from) {
	    shelp[4] = to;
	    m->SetString(shelp);
	    return 1;
	  }
	  else {
	    if (shelp[6]==from) {
	      shelp[6] = to;
	      m->SetString(shelp);
	      return 1;
	    }
	  }
      }
    }
  }
  if (SingleReplaceD(m->left,orig,from,to))  return 1;
  if (SingleReplaceD(m->right,orig,from,to)) return 1;
  return 0;
}


void CFColor::ReplaceG(sknot* m,sknot* m0)
{
  if (m==0) return;
  if (m->op=='*') {
    if (m0==0)m0=m;
    sknot* s1 = m->right;
    sknot* s2 = 0;
    if (m->left->op=='*') s2 = m->left->right;
    else {
      if (m->left->op==0) s2 = m->left;
    }
    if (s2!=0) {
      if (s1->Str().length()==6) {
	if (s1->Str()[0]=='G') {
	  if (s1->Str()[2]==s1->Str()[4]) s1->SetString(string("Ng"));
	  else {
	    // kill G's
	    // replace s1->Str()[2] -> s1->Str()[4]
	    char c = s1->Str()[2];
	    sknot* akt = m0;
	    do {
	      if (m0->left->op=='*' || m0->left->op==0) {
		// right...
		if(m0->right->Str().length()==8) {
		  string shelp = m0->right->Str();
		  for(short int k=1;k<4;k++) if(shelp[2*k]==c)shelp[2*k] = s1->Str()[4];
		  m0->right->SetString(shelp);
		}		
	      }
	      if (m0->left->op==0) {
		// left...
		if(m0->left->Str().length()==8) { 
		  string shelp = m0->left->Str();
		  for(short int k=1;k<4;k++) if(shelp[2*k]==c)shelp[2*k] = s1->Str()[4];
		  m0->left->SetString(shelp);
		}
	      }
	      m0 = m0->left;
	    }
	    while (m0->op=='*');
	    m0 = akt;
	  
	    akt = m;
	    do {
	      
	      if (m->left->op=='*' || m->left->op==0) {
		// right...
		if (m->right->Str().length()==6) {
		  if (m->right->Str()[0]=='G') {
		    string shelp = m->right->Str();
		    if (shelp[2]==c) shelp[2] = s1->Str()[4];
		    else {
		      if (shelp[4]==c) shelp[4] = s1->Str()[4];
		    }
		    m->right->SetString(shelp);
		  }
		}
	      }
	      if (m->left->op==0) {
		// left...
		if (m->left->Str().length()==6) {
		  if (m->left->Str()[0]=='G') {
		    string shelp = m->left->Str();
		    if (shelp[2]==c) shelp[2] = s1->Str()[4];	   
		    else {
		      if (shelp[4]==c) shelp[4] = s1->Str()[4];
		    }
		    m->left->SetString(shelp); 
		  }
		}
	      }
	      
	      m = m->left;
	    }
	    while (m->op=='*');
	    m = akt;
	    s1->SetString(string("1"));
	  }
	}
      }
      if (s2->Str().length()==6) {
	if (s2->Str()[0]=='G') {
	  if (s2->Str()[2]==s2->Str()[4]) s2->SetString(string("Ng"));
	  else {
	    // kill G's
	    // replace s2->Str()[2] -> s2->Str()[4]
	    char c = s2->Str()[2];
	    sknot* akt = m0;
	    do {
	      if (m0->left->op=='*' || m0->left->op==0) {
		// right...
		if(m0->right->Str().length()==8) {
		  string shelp = m0->right->Str();
		  for(short int k=1;k<4;k++) if(shelp[2*k]==c)shelp[2*k] = s2->Str()[4];
		  m0->right->SetString(shelp);
		}		
	      }
	      if (m0->left->op==0) {
		// left...
		if(m0->left->Str().length()==8) { 
		  string shelp = m0->left->Str();
		  for(short int k=1;k<4;k++) if(shelp[2*k]==c)shelp[2*k] = s2->Str()[4];
		  m0->left->SetString(shelp);
		}
	      }
	      m0 = m0->left;
	    }
	    while (m0->op=='*');
	    m0 = akt;
	  
	    akt = m;
	    do {
	      
	      if (m->left->op=='*' || m->left->op==0) {
		// right...
		if (m->right->Str().length()==6) {
		  if (m->right->Str()[0]=='G') {
		    string shelp = m->right->Str();
		    if (shelp[2]==c) shelp[2] = s2->Str()[4];
		    else {
		      if (shelp[4]==c) shelp[4] = s2->Str()[4];
		    }
		    m->right->SetString(shelp);
		  }
		}
	      }
	      if (m->left->op==0) {
		// left...
		if (m->left->Str().length()==6) {
		  if (m->left->Str()[0]=='G') {
		    string shelp = m->left->Str();
		    if (shelp[2]==c) shelp[2] = s2->Str()[4];	   
		    else {
		      if (shelp[4]==c) shelp[4] = s2->Str()[4];
		    }
		    m->left->SetString(shelp); 
		  }
		}
	      }
	      
	      m = m->left;
	    }
	    while (m->op=='*');
	    m = akt;
	    s2->SetString(string("1"));
	  }
	}
      }
    }
  }
  ReplaceG(m->left,m0);
  ReplaceG(m->right,m0);
}

void CFColor::ReplaceF(sknot* m,CharNum_Map& indices,char& c)
{
  int hit;
  do {
    do {
      st.Expand(m);st.Linear(m);
      hit = 0;
      SingleReplaceFT(m,hit,indices,c);
    }
    while (hit>0);
    //any single F ?
    hit = st.Tree2String(m,0).find("F[");
    if (hit==-1) break;
    hit = 0;
    SingleReplaceF(m,hit,indices,c);	   
  }
  while (hit>0);
}

void CFColor::SingleReplaceF(sknot* m,int& hit,CharNum_Map& indices, char& c)
{
  //replace fabc -> T's
  if (m==0) return;
  if (hit>0) return;
  if (m->op=='*') {
    sknot* s = 0;
    if (m->right->op==0) {
      if (m->right->Str()[0]=='F') s = m->right; 
    }
    if (m->left->op==0) {
      if (m->left->Str()[0]=='F')  s = m->left; 
    }
    if (s!=0) {
      hit = 1;
      char A[2],B[2],C[2];
      A[1] = B[1] = C[1] = 0;
      A[0] = s->Str()[2];
      B[0] = s->Str()[4];
      C[0] = s->Str()[6];
      char ii[2],jj[2],kk[2];
      ii[1] = jj[1] = kk[1] = 0;
      ii[0] = DeliverIndex(indices,c);
      c++;
      jj[0] = DeliverIndex(indices,c);
      c++;
      kk[0] = DeliverIndex(indices,c);
      c++;
      s->op = '*';
      string ss;
      ss = string("2*i");
      s->left = st.String2Tree(ss);
           
      ss  = string("T[") + string(A) + string(",") + string(ii) + 
	string(",") + string(jj) + string("]*");
      ss += string("T[") + string(C) + string(",") + string(jj) + 
	string(",") + string(kk) + string("]*");
      ss += string("T[") + string(B) + string(",") + string(kk) + 
	string(",") + string(ii) + string("]-");
      //-
      ss += string("T[") + string(A) + string(",") + string(ii) + 
	string(",") + string(jj) + string("]*");
      ss += string("T[") + string(B) + string(",") + string(jj) + 
	string(",") + string(kk) + string("]*");
      ss += string("T[") + string(C) + string(",") + string(kk) + 
	string(",") + string(ii) + string("]");
      s->right = st.String2Tree(ss);
    }
  }
  SingleReplaceF(m->left,hit,indices,c);
  SingleReplaceF(m->right,hit,indices,c);  
} 

void CFColor::ExtractIndices(sknot* m,CharNum_Map& indices)
{
  list<sknot*> addends;
  st.Addends(m,addends);
  for (list<sknot*>::iterator ita=addends.begin();ita!=addends.end();++ita) {
    //loop over factors and extract
    string argument;
    list<sknot*> factors;
    st.Factors(*ita,factors);
    for (list<sknot*>::iterator itf=factors.begin();itf!=factors.end();++itf) {
      //fill arguments in map
      argument = (*itf)->Str();
      if (argument[0]=='T' || argument[0]=='F') {
	indices.insert(std::make_pair(argument[2],1)); 
	indices.insert(std::make_pair(argument[4],1));
	indices.insert(std::make_pair(argument[6],1));
      }
      if (argument[0]=='D' || argument[0]=='G') {
	indices.insert(std::make_pair(argument[2],1));
	indices.insert(std::make_pair(argument[4],1));
      }
    }
  }
}

char CFColor::DeliverIndex(CharNum_Map& indices,char& c) 
{
  for(;;) {
    if (indices.insert(std::make_pair(c,1)).second) return c;
    else c++;
  }
}

void CFColor::SingleReplaceFT(sknot* m,int& hit,CharNum_Map& indices,char& c)
{
  // change fabc Ta -> T's
  if (m==0) return;
  if (hit>0) return;
  if (m->op=='*') {
    if (m->right->op==0) {
      sknot* s1=0;
      sknot* s2=0;
      if (m->right->Str()[0]=='T') {
	//search F
	s1 = m->right;
	char cT = m->right->Str()[2];
	sknot* akt = m;
	do {
	  if (m->left->op=='*' || m->left->op==0) {
	    if (m->right->Str()[0]=='F') {
	      for (short int i=2;i<7;i+=2)
		if (m->right->Str()[i]==cT) hit = i/2;
	      if (hit>0) s2 = m->right;
	    }
	  }
	  if (hit==0) {
	    if (m->left->op==0) {
	      if (m->left->Str()[0]=='F') {
		for (short int i=2;i<7;i+=2)
		  if (m->left->Str()[i]==cT) hit = i/2;
		if (hit>0) s2 = m->left;
	      }
	    }
	  }
	  m = m->left;
	}
	while (m->op=='*' && hit==0);
	m = akt;
      }
      else {
	if (m->right->Str()[0]=='F') {
	  s2 = m->right;
	  //search F
	  for (short int i=2;i<7;i+=2) {
	    char cT = m->right->Str()[i];
	    sknot* akt = m;
	    do {
	      if (m->left->op=='*' || m->left->op==0) {
		if (m->right->Str()[0]=='T') {
		  if (m->right->Str()[2]==cT) hit = i/2;
		  if (hit>0) s1 = m->right;
		}
	      }
	      if (hit==0) {
		if (m->left->op==0) {
		  if (m->left->Str()[0]=='T') {
		    if (m->left->Str()[2]==cT) hit = i/2;
		    if (hit>0) s1 = m->left;
		  }
		}
	      }
	      m = m->left;
	    }
	    while (m->op=='*' && hit==0);
	    m = akt;
	  }
	}
      }
      if (hit>0) {
	// s1 = T ; s2 = F
	char f1[2];
	char f2[2];
	f1[1] = 0;f2[1] = 0;
	switch (hit) {
	case 1:
	  f1[0] = s2->Str()[4];
	  f2[0] = s2->Str()[6];
	  break;
	case 2:
	  //extra sign
	  f1[0] = s2->Str()[6];
	  f2[0] = s2->Str()[2];
	  break;
	case 3:
	  f1[0] = s2->Str()[2];
	  f2[0] = s2->Str()[4];
	  break;
	default:cerr<<"Wrong hit: "<<hit<<endl;
	}
	char t1[2];
	char t2[2];
	char cc[2];
	cc[0] = DeliverIndex(indices,c);
	c++;
	cc[1] = 0;
	t1[1] = 0;t2[1] = 0;
	t1[0] = s1->Str()[4];t2[0] = s1->Str()[6];
	
	s1->SetString(string("i"));
	s2->op = '-';

	string s;
	s  = string("T[") + string(f1) + string(",") + string(t1) + 
	  string(",") + string(cc) + string("]*T[");
	s += string(f2) + string(",") + string(cc) + 
	  string(",") + string(t2) + string("]");
	s2->right = st.String2Tree(s);
	//-
	s  = string("T[") + string(f2) + string(",") + string(t1) + 
	  string(",") + string(cc) + string("]*T[");
	s += string(f1) + string(",") + string(cc) + 
	  string(",") + string(t2) + string("]");
	s2->left = st.String2Tree(s);
      }
    }
  }
  SingleReplaceFT(m->left,hit,indices,c);
  SingleReplaceFT(m->right,hit,indices,c);  
} 
