#include <stdio.h>
#include <iostream>
#include "CFColor.H"
#include "Run_Parameter.H"
#include "Message.H"

#include "IO_Handler.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace std;


std::string CFColor::noname=string("noname");

CFColor::CFColor(int N,Single_Amplitude* first,string& pID)
{
  Single_Amplitude* m1 = first;

  mcount = 0;

  // on will be changed

  while (m1) {
    mcount++;
    m1->on = 1;
    m1 = m1->Next;  
  }
  
  /*
  // too big !!!
  CFC = new Complex*[mcount];
  for (short int j=0;j<mcount;j++) CFC[j] = new Complex[mcount];
  */  

  //if already on
  

  if (pID!=noname) {
    // look for file
    char name[100];
    //    sprintf(name,"%s/CFColor.dat",pID.c_str());
    sprintf(name,"%s.col",(string("Process/")+pID).c_str());

    fstream test;
    test.open(name,ios::in); 

    if (test) { 
      test.close();

      IO_Handler ioh;
      ioh.SetFileNameRO(name);
      int model, rmcount;

      model  = ioh.Input<int>("");
      rmcount = ioh.Input<int>("");
      ncount = ioh.Input<int>("");


      if (mcount==rmcount) {

	id  = ioh.ArrayInput<int>("",mcount);
	CFC = ioh.MatrixInput<Complex>("",ncount,ncount);



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

	/*
	  fstream from;
	  from.open(name,ios::in);
	  //read in
	  from.precision(15);
	  short int i,j;
	  for(;from;) {      
	  from>>i>>j;
	  from>>CFC[i][j];
	  if ((i==mcount-1) && (j==mcount-1)) break;
	  }
	*/

	msg.Debugging()<<"File "<<name<<" read."<<endl;   
	return; 
      }
      else {
	msg.Tracking()<<"File "<<name<<" has to be recreated!"<<endl;   

      }
 
    }
    else {
      test.close(); 
    
      msg.Tracking()<<"File "<<name<<" not found."<<endl;   
    }
  }

  Color_Function* cm1;
  Color_Function* cm2;

  //looking for zero color structure
  int sw1;  

  sw1 = 0;
  Single_Amplitude* m2;
  m1 = first;
  int counter=0;
  while (m1) {
    if (m1->Get_CFlist()==NULL) {
      sw1 = 1;
      break;
    }
    
    counter++;
    m1 = m1->Next;
  }
  if (sw1) {
    //no color structure
    /*
      for(short int i=0;i<mcount;i++) {
      for (short int j=0;j<mcount;j++) CFC[i][j] = 1.;      
      }
    */

    // matrix
    CFC = new Complex*[1];
    CFC[0] = new Complex[1];
    CFC[0][0]=1.;

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
    //changing props
    msg.Debugging()<<"Finding diagrams with same color structure..."<<endl;
    
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
	  for (short int i=0;i<3;i++) {
	    if (cm1->type==cf::D && i==2) break;    
	    if (cm1->partarg[i]==j) {
	      for (short int k=0;k<3;k++) {
		if (cm1->type==cf::D && k==2) break;
		if ((k!=i) && (cm1->partarg[k]>99) && 
		    ((cm1->partarg[k]<120) || (cm1->partarg[k]>150)))
		  hit = cm1->partarg[k];
	      }
	    }
	  } 
	  if (hit>0) break;
	  cm1 = cm1->Next;
	}
	if (hit>0) {
	  //replace hit with prop
	  cm1 = m1->Get_CFlist();
	  while (cm1) {
	    for (short int i=0;i<3;i++) {
	      if (cm1->type==cf::D && i==2) break;    
	      if (cm1->partarg[i]==hit) cm1->partarg[i] = prop;
	    }
	    cm1 = cm1->Next;
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
	  cm1 = cm1->Next;
	}
	m2 = m1->Next;
	n2 = n1+1;
	while (m2) {
	  if (m2->on) { 
	    cm2 = m2->Get_CFlist();
	    c2 = 0;
	    while (cm2) {
	      c2++;
	      cm2 = cm2->Next;
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
		  cm2 = cm2->Next;
		}
		if (hit2==0) {
		  hit = 0;
		  break;
		}
		cm1 = cm1->Next;
	      }
	      if (hit && hit2) {
		//equal
		ncount--;
		m2 ->on = 0;
		id[n2] = hit*n1;
		if (hit<0 && n1==0) {
		  msg.Error()<<" ERROR in CFColor Constructor "<<endl;
		  abort();
		}
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
    msg.Tracking()<<endl;
    msg.Tracking()<<ncount<<" different color structures left"<<endl;

    map = new int[mcount];
    int cc=0;
    msg.Debugging()<<" id=";
    for (int m=0; m<mcount; ++m) {
      msg.Debugging()<<id[m]<<" ";
    }
    msg.Debugging()<<endl;
    msg.Debugging()<<" map=";
    for (int m=0; m<mcount; ++m) {
      if (id[m]==mcount) {
	map[m]=cc;
	cc++;
      }
      else {
	map[m]=map[iabs(id[m])];
      }
      
      msg.Debugging()<<map[m]<<" ";
    }
    msg.Debugging()<<endl;

    msg.Tracking()<<cc<<" different color structures left (cc)"<<endl;

    // generate "reduced matrix"
    CFC = new Complex*[ncount];
    for (short int j=0;j<ncount;j++) CFC[j] = new Complex[ncount];
    
    string Cstr;

    m1 = first;
    c1 = 0;
    while (m1) { 
      if (m1->on) {
	m2 = m1;
	c2 = c1;
	while (m2) {
	  if (m2->on) {
	    st.Reset();
	    char c = 'M';
	    sknot* s1 = st.String2Tree(m1->CFColstring);
	    sknot* s2 = st.String2Tree(m2->CFColstringC);
	    
	    ReplaceF(s1,c);
	    ReplaceF(s2,c);
	    sknot* m = st.newsk();
	    m->op = '*';
	    m->right = s1;
	    m->left  = s2;
	    
	    st.Expand(m);st.Linear(m);st.Sort(m);
	    ReplaceT(m);
	    st.Expand(m);
	    st.Linear(m);
	    ReplaceD(m);
	    //AORGTOOLS::msg.Out()<<"After D's"<<endl;
	    CFC[map[c1]][map[c2]] = st.eval(m);
	    CFC[map[c2]][map[c1]] = conj(CFC[map[c1]][map[c2]]);
	    msg.Debugging()<<"+";AORGTOOLS::msg.Out().flush();
	  }
	  m2 = m2->Next;
	  c2++;
	}
	msg.Debugging()<<endl;
      }
      m1 = m1->Next;
      c1++;
    }

    /*
    // too big !!
    int a,b;
    
    for (short int i=0;i<mcount;i++) {
    if (id[i]==mcount) a = i;
    else a = iabs(id[i]); 
    for (short int j=0;j<mcount;j++) {
    if (id[j]==mcount) b = j;
    else b = iabs(id[j]);
    int sign = 1;
    if (a!=i && id[i]<0) sign *= -1;
    if (b!=j && id[j]<0) sign *= -1;
    //AORGTOOLS::msg.Out()<<"Negative Sign....!"<<endl;
    CFC[i][j] = double(sign)*CFC[a][b];
    }
    }
    */

  }

  if (pID!=noname) Output(pID);

  // check if Matrix can be reduce even further!
  int idcc=0;
  int * idid = new int[ncount];
  for (int i=0; i<ncount; ++i) 
    idid[i]=-1;
  for (int i=0; i<ncount; ++i) {
    if (idid[i]==-1) { idid[i]=idcc; ++idcc; }
    for (int j=i+1; j<ncount; ++j) {
      int hit=1;
      Complex factor(0.,0.); //CFC[i][0]/CFC[j][0];
      for (int k=0; k<ncount; ++k) {
	if (factor==Complex(0.,0.) && CFC[j][k]!=Complex(0.,0.))
	  factor =CFC[i][k]/CFC[j][k];
	//	if (CFC[i][k]!=factor*CFC[j][k]) {
	Complex diff=(CFC[i][k]-factor*CFC[j][k]);
	if ( dabs(diff.real())<1.e-10 && dabs(diff.imag())<1.e-10 ) {
	  hit=0;
	  break;
	}
      }
      if (hit) {
	msg.Events()<<"Color Matrix could  be further simplified ("<<j<<" ->"<<i<<" with "<<factor<<endl;
	idid[j] =idid[i];
      }
    }
  }
  if (ncount!=idcc) msg.Events()<<"Color Matrix could  be further simplified to "<<idcc<<" structure(s)"<<endl;

  delete [] idid;


  //make m's on again
  m1 = first;
  while (m1) {   
    m1->on = 1;
    m1 = m1->Next;
  }
  
}

void CFColor::Output(string & dirname) {
  // check dir!!!!

  // write out
  //   "pocess name"
  // new_process name (ie. number of graphs, contributing helicities, max number of Z-functions)

  char name[100];
  //  sprintf(name,"%s/CFColor.dat",(string("Process/")+dirname).c_str());
  sprintf(name,"%s.col",(string("Process/")+dirname).c_str());
  msg.Tracking()<<" Writing Color Information to : "<<name<<endl;
  //write out
  IO_Handler ioh;
  ioh.SetFileName(name);

  ioh.Output("",int(rpa.me.Model()));  
  ioh.Output("",mcount);          // no of ampls
  ioh.Output("",ncount);          // size of colormatrix

  ioh.ArrayOutput("",id,mcount);
  ioh.MatrixOutput("",CFC,ncount,ncount);


  /*
    ofstream to;
    //  to.open(name,ios::out,436);
    to.open(name,ios::out);
    to.precision(15);
  
    for(short int i=0;i<mcount;i++) {
    for (short int j=0;j<mcount;j++) {
    to<<i<<"     "<<j<<"     "<<CFC[i][j]<<endl;
    }
    }
    msg.Debugging()<<"File "<<name<<" saved."<<endl;  
  */

}


int CFColor::CompareArg(int a,int b, int c,Color_Function* cm1,Color_Function* cm2)
{
  if (cm1->partarg[a]!=cm2->partarg[0]) return 0;
  if (cm1->partarg[b]!=cm2->partarg[1]) return 0;
  if (cm1->type==cf::D) return 1;
  if (cm1->partarg[c]!=cm2->partarg[2]) return 0;

  return 1;
}

int CFColor::Compare(Color_Function* cm1,Color_Function* cm2)
{
  if (cm1->type!=cm2->type) return 0;
  if (CompareArg(0,1,2,cm1,cm2)) return 1;
  if (cm1->type==cf::D) 
    if (CompareArg(1,0,2,cm1,cm2)) return 1;
  if (cm1->type==cf::F) {
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
	    s  = string("0.33*");
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

void CFColor::ReplaceD(sknot* m)
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
      if (s1->Str().length()==6) {
	if (s1->Str()[0]=='D') {
	  if (s1->Str()[2]==s1->Str()[4]) s1->SetString(string("3"));
	  else {
	    // kill D's
	    // replace s1->Str()[2] -> s1->Str()[4]
	    char c = s1->Str()[2];
	    sknot* akt = m;
	    do {
	      
	      if (m->left->op=='*' || m->left->op==0) {
		// right...
		if (m->right->Str().length()==6) {
		  if (m->right->Str()[0]=='D') {
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
		  if (m->left->Str()[0]=='D') {
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
	if (s2->Str()[0]=='D') {
	  if (s2->Str()[2]==s2->Str()[4]) s2->SetString(string("3"));
	}
      }
    }
  }

  ReplaceD(m->left);
  ReplaceD(m->right);
}

void CFColor::ReplaceF(sknot* m,char& c)
{
  int hit;
  do {
    do {
      st.Expand(m);st.Linear(m);
      hit = 0;
      Single_ReplaceFT(m,hit,c);
    }
    while (hit>0);
    //any single F ?
    hit = st.Tree2String(m,0).find("F");
    if (hit==-1) break;
    hit = 0;
    Single_ReplaceF(m,hit,c);	   
  }
  while (hit>0);
}

void CFColor::Single_ReplaceF(sknot* m,int& hit,char& c)
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
      if (c==']' || c=='[') c++;
      ii[0] = c;c++;
      if (c==']' || c=='[') c++;
      jj[0] = c;c++;
      if (c==']' || c=='[') c++;
      kk[0] = c;c++;

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
  Single_ReplaceF(m->left,hit,c);
  Single_ReplaceF(m->right,hit,c);  
} 

void CFColor::Single_ReplaceFT(sknot* m,int& hit,char& c)
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
	char c = m->right->Str()[2];
	sknot* akt = m;
	do {
	  if (m->left->op=='*' || m->left->op==0) {
	    if (m->right->Str()[0]=='F') {
	      for (short int i=2;i<7;i+=2)
		if (m->right->Str()[i]==c) hit = i/2;
	      if (hit>0) s2 = m->right;
	    }
	  }
	  if (hit==0) {
	    if (m->left->op==0) {
	      if (m->left->Str()[0]=='F') {
		for (short int i=2;i<7;i+=2)
		  if (m->left->Str()[i]==c) hit = i/2;
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
	    char c = m->right->Str()[i];
	    sknot* akt = m;
	    do {
	      if (m->left->op=='*' || m->left->op==0) {
		if (m->right->Str()[0]=='T') {
		  if (m->right->Str()[2]==c) hit = i/2;
		  if (hit>0) s1 = m->right;
		}
	      }
	      if (hit==0) {
		if (m->left->op==0) {
		  if (m->left->Str()[0]=='T') {
		    if (m->left->Str()[2]==c) hit = i/2;
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
	if (c==']' || c=='[') c++;
	cc[0] = c;cc[1] = 0;
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
	c++;
      }
    }
  }
  Single_ReplaceFT(m->left,hit,c);
  Single_ReplaceFT(m->right,hit,c);  
} 














