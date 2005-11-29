#include "MHV_Calculator.H"
#include "Pfunc.H"


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;



MHV_Calculator::MHV_Calculator(int part, Basic_Sfuncs* BS,int* signlist,int* plist,int p) :
n_part(part), p_BS(BS) ,print(p)
{
  m_dummyarg = new int[2*part];
  m_signlist = new int[part];
  m_plist = new int[part];
  for (int i=0;i<part;i++) {
    m_dummyarg[i]=i;
    m_dummyarg[part+i]=i;
    m_signlist[i]=signlist[i];
    m_plist[i]=plist[i];
  }
  m_arg = new int[part];

  m_qlist[0] = 0;
  for (int i=1;i<5;i++) m_qlist[i]=-1;
  Make_Qlist(m_dummyarg,m_plist,m_qlist,n_part);  

#ifndef Mom_H
  if (n_part>5) {
    Pfunc_List pl;
    int tt = (1<<(n_part-1))-1;
    for (int i=3;i<tt;i++) {
      int kk = 0;
      for (int j=0;j<n_part-1;j++) if (i&(1<<j)) kk++;
      if (kk>1) {
	int n=1;
	Pfunc* pf = new Pfunc(kk+1);
	for (int j=0;j<n_part-1;j++) if (i&(1<<j)) {
	  pf->arg[n]=j;
	  n++;
  	}
	pl.push_back(pf);
      }
    }
    p_BS->BuildMomlist(pl);
    for (Pfunc_Iterator it=pl.begin();it!=pl.end();it++) delete (*it);
  }
#endif
}

MHV_Calculator::~MHV_Calculator() 
{ 
  delete[] m_arg;
  delete[] m_signlist;
  delete[] m_plist;
}


Complex MHV_Calculator::Differential(int* perm)
{
  int* signlist   = new int[n_part];
  for (int i=0;i<n_part;i++) signlist[i]=m_signlist[perm[i]];
  int sum = 0;
  for (int i=0;i<n_part;i++) sum+=signlist[i];
  for (int i=0;i<n_part;i++) m_arg[i]=perm[i];
  int nh=(n_part-sum)/2;
  int ph=n_part - nh;


  if ((nh<2 || ph<2) && (n_part!=3 || nh==0 || ph ==0)) return 0.;

  if (Min(nh,ph)>3){
    std::cout<<"Error in MHV_Calculator::Differential: Amplitude not implemented!"<<std::endl;
    abort();
  }


  int qlist[] = {0,-1,-1,-1,-1};
  Make_Qlist(m_arg,m_plist,qlist,n_part); 

// pure gluonic amilitude 
  if (!qlist[0]) {
    if (nh==2) return Elementary_MHV_Amplitude(m_arg,signlist,n_part);
    if (ph==2) return Elementary_MHVbar_Amplitude(m_arg,signlist,n_part);
    
    if (nh<ph) return MHV_Amplitude(m_arg,signlist,nh);
    return MHVbar_Amplitude(m_arg,signlist,ph);
  }    
 
 // amplitudes with quarks
  if (qlist[0]==1 || qlist[0]==3) {
    std::cout<<"Error in MHV_Calculator::Differential: Odd number of quarks"<<std::endl;
    abort();
  }
  if (!Check_Qlist(perm,signlist,m_plist,qlist)) {
    std::cout<<"Error in MHV_Calculator::Differential: Wrong flavors or helicities for quarks"<<std::endl;
    abort();
  }

  if (qlist[0]==2) {
    if (nh==2) return Elementary_MHVQ2_Amplitude(m_arg,signlist,qlist,n_part);
    if (ph==2) return Elementary_MHVQ2bar_Amplitude(m_arg,signlist,qlist,n_part);
    
    if (nh<ph)  {
      if (signlist[qlist[1]]==-1) {
	for (int y=0;y<n_part;y++) {
	  m_arg[y]=perm[(qlist[1]+y)%n_part];
	  signlist[y]=m_signlist[m_arg[y]];
	}
	return MHVQ2_Amplitude(m_arg,signlist,nh);
      }
      else {
	for (int y=0;y<n_part;y++) {
	  m_arg[y]=perm[(qlist[2]+y)%n_part];
	  signlist[y]=m_signlist[m_arg[y]];
	}
	return -MHVQ2_Amplitude(m_arg,signlist,nh);
      }
    }
    if (signlist[qlist[1]]==1) {
      for (int y=0;y<n_part;y++) {
	m_arg[y]=perm[(qlist[1]+y)%n_part];
	signlist[y]=m_signlist[m_arg[y]];
      }
      return MHVQ2bar_Amplitude(m_arg,signlist,nh);
    }
    else {
      for (int y=0;y<n_part;y++) {
	m_arg[y]=perm[(qlist[2]+y)%n_part];
	signlist[y]=m_signlist[m_arg[y]];
      }
      return -MHVQ2bar_Amplitude(m_arg,signlist,nh);
    }
  }  

 
  if (qlist[0]==4) {
    //if (nh==2) return Elementary_MHVQ4_Amplitude(m_arg,signlist,qlist,n_part);
    //if (ph==2) return Elementary_MHVQ4bar_Amplitude(m_arg,signlist,qlist,n_part);
    
    //if (nh<ph) return MHVQ4_Amplitude(m_arg,signlist,nh);
    
    return MHVQ4bar_Amplitude(m_arg,signlist,ph);
  }     
  return 0;
}


int* MHV_Calculator::GetQlist() {
  return m_qlist;
}

int* MHV_Calculator::GetPlist() {
  return m_plist;
}



Complex MHV_Calculator::Elementary_MHV_Amplitude(int* perm,int* signlist,int part)
{
  int l;
  int v1=-1,v2=-1;
  for (l=0;l<part-1 && v1<0;l++) if (signlist[l]==-1) v1=perm[l];
  for (;l<part && v2<0;l++) if (signlist[l]==-1) v2=perm[l];
  if (v2<0) abort();
  if (print==2) cout<<"<"<<v1<<","<<v2<<">^4 "; //print
  Complex sm = p_BS->S0(v1,v2);
  sm = pow(sm,4);
  for (l=0;l<part-1;l++){ 
    if (print==2) cout<<"/<"<<perm[l]<<","<<perm[l+1]<<">"; //print
    sm/=p_BS->S0(perm[l],perm[l+1]);
  }
  if (print==2) cout<<"/<"<<perm[part-1]<<","<<perm[0]<<">"; //print
  sm/=p_BS->S0(perm[part-1],perm[0]);
 
  if (print==1) {                                           //print
    cout<<"A(";                                             //print
    for (int i=0;i<part-1;i++) {                            //print
      cout<<perm[i];                                        //print
      if (signlist[i]==1) cout<<"+,";                       //print
      else cout<<"-,";                                      //print
    }                                                       //print
    cout<<perm[part-1];                                     //print
    if (signlist[part-1]==1) cout<<"+)";                    //print
    else cout<<"-)";	                                    //print 
  }                                                         //print  

  delete[] signlist;
  return sm;
}


Complex MHV_Calculator::Elementary_MHVbar_Amplitude(int* perm,int* signlist,int part)
{
  int l;
  int v1=-1,v2=-1;
  for (l=0;l<part-1 && v1<0;l++) if (signlist[l]==1) v1=perm[l];
  for (;l<part && v2<0;l++) if (signlist[l]==1) v2=perm[l];
  if (v2<0) abort();
  if (print==2) cout<<"["<<v1<<","<<v2<<"]^4 "; //print
  Complex sm = p_BS->S1(v1,v2);
  sm = pow(sm,4);
  for (l=0;l<part-1;l++) {
    if (print==2) cout<<"/["<<perm[l]<<","<<perm[l+1]<<"]"; //print
    sm/=p_BS->S1(perm[l],perm[l+1]);
  }
  if (print==2) cout<<"/["<<perm[part-1]<<","<<perm[0]<<"]"; //print
  sm/=p_BS->S1(perm[part-1],perm[0]);
 
  if (print==1) {                                           //print
    cout<<"A(";                                             //print
    for (int i=0;i<part-1;i++) {                            //print
      cout<<perm[i];                                        //print
      if (signlist[i]==1) cout<<"+,";                       //print
      else cout<<"-,";                                      //print
    }                                                       //print
    cout<<perm[part-1];                                     //print
    if (signlist[part-1]==1) cout<<"+)";                    //print
    else cout<<"-)";	                                    //print
  }                                                         //print 

  delete[] signlist;
  return sm;
}

int MHV_Calculator::CountSign(int* signlist,int start,int end,int sg)
{
  int sum = 0;
  for (int i=start;i<end;i++) sum+=signlist[i];
  if (sg==-1) return (end-start-sum)/2;
  return (end-start+sum)/2;
}

Complex MHV_Calculator::MHV_Amplitude(int* perm,int* signlist,int vhel)
{
  if (vhel==2) return Elementary_MHV_Amplitude(perm,signlist,n_part);
  Complex amp(0.,0.);
  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2) {
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[i+k];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[(j+k)%n_part];
	switch (cs) {
	case 1:
	  sl1[j-i]=-1;
	  sl2[n_part-j+i]=1;
	  break;
	case 2:
	  sl1[j-i]=1;
	  sl2[n_part-j+i]=-1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	if (print && (i || (j-2))) cout<<" + ";                  //print
	Complex v = Elementary_MHV_Amplitude(&m_dummyarg[i],sl1,j-i+1);
	m_dummyarg[j] = perm[j];
	m_dummyarg[n_part+i] = pn;
	if (print) cout<<"*";                                  //print
	Complex vp = Elementary_MHV_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1);
	v*=vp;
	m_dummyarg[n_part+i] = perm[i];
	if (print) cout<<"*1/("<<pn<<")^2";                    //print
	v /= (p_BS->Momentum(pn)).Abs2();
	amp+=v;
      } 
    }    
  }
  delete[] signlist;
  return amp;
}


Complex MHV_Calculator::MHVbar_Amplitude(int* perm,int* signlist,int vhel)
{
  if (vhel==2) return Elementary_MHVbar_Amplitude(perm,signlist,n_part);  
  Complex amp(0.,0.);

  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2) {
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[i+k];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[(j+k)%n_part];
	switch (cs) {
	case 1:
	  sl1[j-i]=1;
	  sl2[n_part-j+i]=-1;
	  break;
	case 2:
	  sl1[j-i]=-1;
	  sl2[n_part-j+i]=1;
	} 	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	if (print && (i || (j-2))) cout<<" + ";                  //print
	Complex v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],sl1,j-i+1);
	m_dummyarg[j] = perm[j];
	m_dummyarg[n_part+i] = pn;
	if (print) cout<<"*";                                  //print
	v *= Elementary_MHVbar_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1);
	m_dummyarg[n_part+i] = perm[i];
	if (print) cout<<"*1/("<<pn<<")^2";                    //print
	v /= (p_BS->Momentum(pn)).Abs2();
	amp+=v;
      } 
    }    
  }
  delete[] signlist;
  return amp;
}


void MHV_Calculator::Make_Qlist(int* perm,int* plist,int* qlist,int part) 
{
  int nq=0; 
  for (int i=0;i<part;i++) {
    if ( !(plist[perm[i]]/6) && (plist[perm[i]]!=0)) {
       nq++;
       qlist[nq]=i;
    }
    else  {
      if ((plist[perm[i]]-21)*(plist[perm[i]]+21))  {
	std::cout<<"Error in MHV_Calculator::Make_Qlist: Amplitude not implemented!"<<std::endl;
	abort();  
	}
    } 
    if (nq==5) {
      std::cout<<"Error in MHV_Calculator::Make_Qlist: Too many quarks"<<std::endl;
      abort();
    }
  }
  qlist[0]=nq;
  return;
}

bool MHV_Calculator::Check_Qlist(int* perm,int* signlist,int* plist,int* qlist) 
{
  if (qlist[0]==2) {
    if ((signlist[qlist[1]]==-signlist[qlist[2]]) && (plist[perm[qlist[1]]]==-plist[perm[qlist[2]]])) return 1;
  }
  if (qlist[0]==4) {
    int i;
    for (i=2;(signlist[qlist[1]]+signlist[qlist[i]])*(plist[perm[qlist[1]]]+plist[perm[qlist[i]]]);i++) if (i>3) return 0;     
    int sh=signlist[qlist[1]], sf=plist[perm[qlist[1]]];
    for (i=2;i<5;i++) {
      sh+=signlist[qlist[i]];
      sf+=plist[perm[qlist[i]]];
    }
    if (sh==0 && sf==0 ) return 1;
    return 0;
  }
  return 0;
}


Complex MHV_Calculator::Elementary_MHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part)
{
  int l;
  int v1=-1;  
  for (l=0;l<part && v1<0;l++) if (signlist[l]==-1 && l!=qlist[1] && l!=qlist[2]) v1=perm[l];
  if (v1<0) abort();  
  if (print==2) cout<<"<"<<v1<<","<<perm[qlist[1]]<<">"; //print
  Complex sm = p_BS->S0(v1,perm[qlist[1]]);
  if (signlist[qlist[1]]==-1) {
    sm = pow(sm,3);
    if (print==2) cout<<"^3 "; //print
  }
  if (print==2) cout<<"<"<<v1<<","<<perm[qlist[2]]<<">"; //print
  Complex sm2 = p_BS->S0(v1,perm[qlist[2]]);
  if (signlist[qlist[2]]==-1) { 
    sm2 = -pow(sm2,3);
    if (print==2) cout<<"^3 (-1)"; //print
  }
  sm *= sm2;  
 
  for (l=0;l<part-1;l++){
    if (print==2) cout<<"/"<<"<"<<perm[l]<<","<<perm[l+1]<<">"; //print
    sm/=p_BS->S0(perm[l],perm[l+1]);
  }
  if (print==2) cout<<"/"<<"<"<<perm[part-1]<<","<<perm[0]<<">"; //print
  sm/=p_BS->S0(perm[part-1],perm[0]);

  if (print==1) {                                           //print
    cout<<"A(";                                             //print
    for (int i=0;i<part-1;i++) {                            //print
      if (i==qlist[1] || i==qlist[2])  cout<<"q";           //print
      cout<<perm[i];                                        //print
      if (signlist[i]==1) cout<<"+,";                       //print
      else cout<<"-,";                                      //print
    }                                                       //print
    if (part-1==qlist[1] || part-1==qlist[2])  cout<<"q";   //print
    cout<<perm[part-1];                                     //print
    if (signlist[part-1]==1) cout<<"+)";                    //print
    else cout<<"-)";	                                    //print
  }                                                         //print  

  delete[] signlist;    
  return sm;
}



Complex MHV_Calculator::Elementary_MHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part)
{ 
  int l;
  int v1=-1;  
  for (l=0;l<part && v1<0;l++) if (signlist[l]==1 && l!=qlist[1] && l!=qlist[2]) v1=perm[l];
  if (v1<0) abort();  
  if (print==2) cout<<"["<<v1<<","<<perm[qlist[1]]<<"]"; //print
  Complex sm = p_BS->S1(v1,perm[qlist[1]]);
  if (signlist[qlist[1]]==1) {
    sm = pow(sm,3);
    if (print==2) cout<<"^3 "; //print
  }
  if (print==2) cout<<"["<<v1<<","<<perm[qlist[2]]<<"]"; //print
  Complex sm2 = p_BS->S1(v1,perm[qlist[2]]);
  if (signlist[qlist[2]]==1) { 
    sm2 = -pow(sm2,3);
    if (print==2) cout<<"^3 (-1)"; //print
  }
  sm *= sm2;  
  
  for (l=0;l<part-1;l++){
    if (print==2) cout<<"/"<<"["<<perm[l]<<","<<perm[l+1]<<"]"; //print
    sm/=p_BS->S1(perm[l],perm[l+1]);
  }
  if (print==2) cout<<"/"<<"["<<perm[part-1]<<","<<perm[0]<<"]"; //print
  sm/=p_BS->S1(perm[part-1],perm[0]);

  if (print==1) {                                           //print
    cout<<"A(";                                             //print
    for (int i=0;i<part-1;i++) {                            //print
      if (i==qlist[1] || i==qlist[2])  cout<<"q";           //print
      cout<<perm[i];                                        //print
      if (signlist[i]==1) cout<<"+,";                       //print
      else cout<<"-,";                                      //print
    }                                                       //print
    if (part-1==qlist[1] || part-1==qlist[2])  cout<<"q";   //print
    cout<<perm[part-1];                                     //print
    if (signlist[part-1]==1) cout<<"+)";                    //print
    else cout<<"-)";	                                    //print
  }                                                         //print 

  delete[] signlist; 
  return sm;
}


Complex MHV_Calculator::Elementary_MHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part)
{
  int l;
  int v1=-1, v2=-1;
  for (l=0;l<3 && v1<0;l++) if (signlist[qlist[l+1]]==-1) v1=perm[qlist[l+1]];
  if (v1<0) abort();
  Complex sm(-1.0,0.0), sm2;
  sm=pow(sm,l+1);
  if (l==1 && signlist[qlist[4]]==-1) sm*= (-1);
  if (print==2 && sm==Complex(-1.0,0.0)) cout<<"(-1)"; //print
  for (int i=0;i<3;i++,l++) {
    if (signlist[qlist[l%4+1]]==-1) {
      if (print==2) cout<<"<"<<v1<<","<<perm[qlist[l%4+1]]<<">^3 "; //print
      sm *= p_BS->S0(v1,perm[qlist[l%4+1]]);
      sm = pow(sm,3) ;
      if (v2<0) {
	if (print==2) cout<<"(-1)"; //print
	sm*=(-1) ;
      }
    }
    else {
      if (v2<0)  v2=perm[qlist[l%4+1]];
      else {
	if (print==2) cout<<"<"<<v2<<","<<perm[qlist[l%4+1]]<<">"; //print
	sm2 = p_BS->S0(v2,perm[qlist[l%4+1]]);
      }
    }
  } 
  sm *= sm2;  
  for (l=0;l<part-1;l++) {
    if (print==2) cout<<"/<"<<perm[l]<<","<<perm[l+1]<<">"; //print
    sm/=p_BS->S0(perm[l],perm[l+1]);
  }
  if (print==2) cout<<"/<"<<perm[part-1]<<","<<perm[0]<<">"; //print
  sm/=p_BS->S0(perm[part-1],perm[0]);
  
  if (print==1) {                                                                                 //print
    cout<<"A(";                                                                                   //print
    for (int i=0;i<part-1;i++) {                                                                  //print
      if (i==qlist[1] || i==qlist[2] || i==qlist[3] || i==qlist[4])  cout<<"q";                   //print
      cout<<perm[i];                                                                              //print
      if (signlist[i]==1) cout<<"+,";                                                             //print
      else cout<<"-,";                                                                            //print
    }                                                                                             //print
    if (part-1==qlist[1] || part-1==qlist[2] || part-1==qlist[3] || part-1==qlist[4])  cout<<"q"; //print
    cout<<perm[part-1];                                                                           //print
    if (signlist[part-1]==1) cout<<"+)";                                                          //print
    else cout<<"-)";	                                                                          //print
  }                                                                                               //print 

  delete[] signlist; 
  return sm;
}



Complex MHV_Calculator::Elementary_MHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part)
{
  int l;
  int v1=-1, v2=-1;
  for (l=0;l<3 && v1<0;l++) if (signlist[qlist[l+1]]==1) v1=perm[qlist[l+1]];
  if (v1<0) abort();
  Complex sm(-1.0,0.0), sm2;
  sm=pow(sm,l+1);
  if (l==1 && signlist[qlist[4]]==1) sm*= (-1);
  if (print==2 && sm==Complex(-1.0,0.0)) cout<<"(-1)"; //print
  for (int i=0;i<3;i++,l++) {
    if (signlist[qlist[l%4+1]]==1) {
      if (print==2) cout<<"["<<v1<<","<<perm[qlist[l%4+1]]<<"]^3 "; //print
      sm *= p_BS->S1(v1,perm[qlist[l%4+1]]);
      sm = pow(sm,3) ;
      if (v2<0) {
	if (print==2) cout<<"(-1)"; //print
	sm*=(-1) ;
      }
    }
    else {
      if (v2<0)  v2=perm[qlist[l%4+1]];
      else {
	if (print==2) cout<<"["<<v2<<","<<perm[qlist[l%4+1]]<<"]"; //print
	sm2 = p_BS->S1(v2,perm[qlist[l%4+1]]);
      }
    }
  } 
  sm *= sm2;  
  for (l=0;l<part-1;l++){
    if (print==2) cout<<"/["<<perm[l]<<","<<perm[l+1]<<"]"; //print
    sm/=p_BS->S1(perm[l],perm[l+1]);
  }
  if (print==2) cout<<"/["<<perm[part-1]<<","<<perm[0]<<"]"; //print
  sm/=p_BS->S1(perm[part-1],perm[0]);
 
  if (print==1) {                                                                                 //print
    cout<<"A(";                                                                                   //print
    for (int i=0;i<part-1;i++) {                                                                  //print
      if (i==qlist[1] || i==qlist[2] || i==qlist[3] || i==qlist[4])  cout<<"q";                   //print
      cout<<perm[i];                                                                              //print
      if (signlist[i]==1) cout<<"+,";                                                             //print
      else cout<<"-,";                                                                            //print
    }                                                                                             //print
    if (part-1==qlist[1] || part-1==qlist[2] || part-1==qlist[3] || part-1==qlist[4])  cout<<"q"; //print
    cout<<perm[part-1];                                                                           //print
    if (signlist[part-1]==1) cout<<"+)";                                                          //print
    else cout<<"-)";	                                                                          //print
  }                                                                                               //print
                                                       
  delete[] signlist; 
  return sm;
}


Complex MHV_Calculator::MHVQ2_Amplitude(int* perm,int* signlist,int vhel)
{
  Complex amp(0.,0.);  
  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2) {
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[i+k];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[(j+k)%n_part];
	switch (cs) {
	case 1:
	  sl1[j-i]=-1;
	  sl2[n_part-j+i]=1;
	  break;
	case 2:
	  sl1[j-i]=1;
	  sl2[n_part-j+i]=-1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;

	int qlist1[]={0,-1,-1};
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	bool empty(false);
	Complex v;
// pure gluonic amilitude 
	if (!qlist1[0]) {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHV_Amplitude(&m_dummyarg[i],sl1,j-i+1);
	}

// amplitudes with quarks
	else if (qlist1[0]==2)   {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = -sl1[qlist1[1]];
	  v *= Elementary_MHVQ2_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && !(sl1[qlist1[1]]==-1 && cs==1) && !(sl1[qlist1[1]]==1 && cs==2)) {
	  if (print && (i || (j-2))) cout<<" + (-1)";                  //print
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else empty=true;

	m_dummyarg[j] = perm[j];
	if (!empty) {
	  m_dummyarg[n_part+i] = pn;
	  int qlist2[]={0,-1,-1};
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,n_part-j+i);
	  if (print) cout<<"*";                                  //print
	  Complex vp;
// pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHV_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1);

// amplitudes with quarks
	  if (qlist2[0]==2)   {
	    vp = -sl2[qlist2[1]];
	    vp *= Elementary_MHVQ2_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  }
	  if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=n_part-j+i;
	    vp = Elementary_MHVQ2_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  }
	  v*=vp;
	  m_dummyarg[n_part+i] = perm[i];
	  if (print) cout<<"*1/("<<pn<<")^2";                    //print
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
      } 
    }    
  }
  delete[] signlist;
  return amp;
}
 

Complex MHV_Calculator::MHVQ2bar_Amplitude(int* perm,int* signlist,int vhel)
{
  Complex amp(0.,0.);
  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2) {
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[i+k];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[(j+k)%n_part];
	switch (cs) {
	case 1:
	  sl1[j-i]=1;
	  sl2[n_part-j+i]=-1;
	  break;
	case 2:
	  sl1[j-i]=-1;
	  sl2[n_part-j+i]=1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;

	int qlist1[]={0,-1,-1};
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	bool empty(false);
	Complex v;
// pure gluonic amilitude 
	if (!qlist1[0]) {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],sl1,j-i+1);
	}

// amplitudes with quarks
	else if (qlist1[0]==2)   {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = sl1[qlist1[1]];
	  v *= Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && !(sl1[qlist1[1]]==1 && cs==1) && !(sl1[qlist1[1]]==-1 && cs==2)) {
	  if (print && (i || (j-2))) cout<<" + (-1)";                  //print
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else empty=true;

	m_dummyarg[j] = perm[j];
	if (!empty) {
	  m_dummyarg[n_part+i] = pn;
	  int qlist2[]={0,-1,-1};
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,n_part-j+i);
	  if (print) cout<<"*";                                  //print
	  Complex vp;
// pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHVbar_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1);

// amplitudes with quarks
	  if (qlist2[0]==2)   {
	    vp = sl2[qlist2[1]];
	    vp *= Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  }
	  if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=n_part-j+i;
	    vp = Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  }
	  v*=vp;
	  m_dummyarg[n_part+i] = perm[i];
	  if (print) cout<<"*1/("<<pn<<")^2";                    //print
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
      } 
    }
  }
  delete[] signlist;
  return amp;
}
 
Complex MHV_Calculator::MHVQ4_Amplitude(int* perm,int* signlist,int vhel)
{
  Complex amp(0.,0.);  
  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2) {
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[i+k];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[(j+k)%n_part];
	switch (cs) {
	case 1:
	  sl1[j-i]=-1;
	  sl2[n_part-j+i]=1;
	  break;
	case 2:
	  sl1[j-i]=1;
	  sl2[n_part-j+i]=-1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;

	int qlist1[]={0,-1,-1,-1,-1};
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	Complex v;
	int qsign=0;
	for (int l=1;l<qlist1[0]+1;l++) qsign+= sl1[qlist1[l]];
// pure gluonic amilitude 
	if (!qlist1[0]) {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHV_Amplitude(&m_dummyarg[i],sl1,j-i+1);
	}

// amplitudes with quarks
	else if (qlist1[0]==4)   {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHVQ4_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==3 && ((3-qsign)/2)==cs) {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  qlist1[0]=4;
	  qlist1[4]=j-i;
	  v = -Elementary_MHVQ4_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==2 && qsign==0)   {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHVQ2_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && !(sl1[qlist1[1]]==-1 && cs==1) && !(sl1[qlist1[1]]==1 && cs==2)) {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else {
	  m_dummyarg[j] = perm[j];
	  break;
	}

	m_dummyarg[j] = perm[j];
	m_dummyarg[n_part+i] = pn;


	int qlist2[]={0,-1,-1,-1,-1};
	Make_Qlist(&m_dummyarg[j],m_plist,qlist2,n_part-j+i);
	if (print) cout<<"*";                                  //print
// pure gluonic amilitude 
	if (!qlist2[0]) Complex vp = Elementary_MHV_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1);


// amplitudes with quarks
	Complex vp;
	if (qlist2[0]==4)   vp = Elementary_MHVQ4_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	if (qlist2[0]==3) {
	  qlist2[0]=4;
	  qlist2[4]=n_part-j+i;
	  vp = Elementary_MHVQ4_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	}
	if (qlist2[0]==2)   vp = Elementary_MHVQ2_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	if (qlist2[0]==1) {
	  qlist2[0]=2;
	  qlist2[2]=n_part-j+i;
	  vp = Elementary_MHVQ2_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	}
	v*=vp;
	m_dummyarg[n_part+i] = perm[i];
	if (print) cout<<"*1/("<<pn<<")^2";                    //print
	v /= (p_BS->Momentum(pn)).Abs2();
	amp+=v;
      } 
    }    
  }
  delete[] signlist;
  return amp;
}
 

Complex MHV_Calculator::MHVQ4bar_Amplitude(int* perm,int* signlist,int vhel)
{
  Complex amp(0.,0.);  
  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2) {
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[i+k];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[(j+k)%n_part];
	switch (cs) {
	case 1:
	  sl1[j-i]=1;
	  sl2[n_part-j+i]=-1;
	  break;
	case 2:
	  sl1[j-i]=-1;
	  sl2[n_part-j+i]=1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;

	int qlist1[]={0,-1,-1,-1,-1};
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	Complex v;
	bool empty(false);
	int qsign=0;
	for (int l=1;l<qlist1[0]+1;l++) qsign+= sl1[qlist1[l]];
// pure gluonic amilitude 
	if (!qlist1[0]) {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],sl1,j-i+1);
	}

// amplitudes with quarks
	else if (qlist1[0]==4)   {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHVQ4bar_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==3 && ((3+qsign)/2)==cs) {
	  if (print) {if  (i || (j-2)) cout<<" + (-1)"; else cout<<"(-1)";}                  //print
	  qlist1[0]=4;
	  qlist1[4]=j-i;
	  v = -Elementary_MHVQ4bar_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==2 && qsign==0)   {
	  if (print && (i || (j-2))) cout<<" + ";                  //print
	  v = Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && !(sl1[qlist1[1]]==1 && cs==1) && !(sl1[qlist1[1]]==-1 && cs==2)) {
	  if (print) {if  (i || (j-2)) cout<<" + (-1)"; else cout<<"(-1)";}                  //print
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],sl1,qlist1,j-i+1);
	}
	else empty=true;
	m_dummyarg[j] = perm[j];
	if (!empty) {
	  Complex vp;
	  m_dummyarg[n_part+i] = pn;
	  int qlist2[]={0,-1,-1,-1,-1};
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,n_part-j+i);
	  if (print) cout<<"*";                                  //print
// pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHVbar_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1);


// amplitudes with quarks
	  if (qlist2[0]==4)   vp = Elementary_MHVQ4bar_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  if (qlist2[0]==3) {
	    qlist2[0]=4;
	    qlist2[4]=n_part-j+i;
	    vp = Elementary_MHVQ4bar_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  }
	  if (qlist2[0]==2)   vp = Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=n_part-j+i;
	    vp = Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],sl2,qlist2,n_part-j+i+1);
	  }
	  v*=vp;
	  m_dummyarg[n_part+i] = perm[i];
	  if (print) cout<<"*1/("<<pn<<")^2";                    //print
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;cout<<"    "<<v<<endl;
	}
      } 
    } 
  }
  delete[] signlist;
  return amp;
}
