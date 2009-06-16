#include "AMEGIC++/Amplitude/Zfunctions/MHVCalculator.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;




// class MHVCalculator


// constructor

MHVCalculator::MHVCalculator(int part,int* plist) :
    n_part(part), m_dummyarg(0), m_ndummyarg(0), m_dummysl(0), m_ndummysl(0), m_plist(0), m_signlist(0), p_BS(0) 
{
    m_dummyarg = new int[2*part];
    m_dummysl = new int[2*part];
    m_plist = new int[(1<<(part-1))-1];
    m_signlist = new int[part];
    for (int i=0;i<part;i++) {
      m_dummyarg[part+i]=m_dummyarg[i]=i;
      m_plist[i]=plist[i];
    } 
    m_ndummyarg = new int[2*part]; 
    m_ndummysl = new int[2*part];
    for (int i=0;i<part;i++) 	m_ndummyarg[part+i]=m_ndummyarg[i]=i;
    Make_Qlist(m_dummyarg,m_plist,m_qlist,n_part); 
}



// destructor

MHVCalculator::~MHVCalculator() 
{ 
    if (m_ndummyarg) delete[] m_ndummyarg;
    if (m_dummyarg)  delete[] m_dummyarg;
    if (m_ndummysl)  delete[] m_ndummysl;
    if (m_dummysl)   delete[] m_dummysl;
    if (m_plist)     delete[] m_plist; 
    if (m_signlist)  delete[] m_signlist;
}



// public functions

Complex MHVCalculator::Differential(int* perm, int* sg)
{
  for (int i=0;i<n_part;i++) m_signlist[i]=sg[perm[i]]; 
  int sum = 0;   
  for (int i=0;i<n_part;i++) sum+=sg[i];
  int nh= (n_part-sum)/2;
  int ph= n_part-nh;
  
  if ((nh<2 || ph<2) && (n_part!=3 || nh==0 || ph ==0)) return 0.;
  
  if (Min(nh,ph)>4) {
    std::cout<<"Error in MHVCalculator::Differential: Amplitude not implemented!"<<std::endl;
    abort();
  }

 
// pure gluonic amilitude 

  if (!m_qlist[0]) {
    if (nh==2) return Elementary_MHV_Amplitude(perm,m_signlist,n_part);
    if (ph==2) return Elementary_MHVbar_Amplitude(perm,m_signlist,n_part);

    if (nh==3) return NMHV_Amplitude(perm,m_signlist,n_part,nh);
    if (ph==3) return NMHVbar_Amplitude(perm,m_signlist,n_part,ph);

    if (nh==4) return NNMHV_Amplitude(perm,m_signlist,n_part,nh);
    if (ph==4) return NNMHVbar_Amplitude(perm,m_signlist,n_part,ph);
  }    
  
  
 // amplitudes with quarks
 
  int qlist[9];
  Make_Qlist(perm,m_plist,qlist,n_part);
  
  if (m_qlist[0]==1 || m_qlist[0]==3) {
    std::cout<<"Error in MHVCalculator::Differential: Odd number of quarks"<<std::endl;
    abort();
  }
  if (!Check_Qlist(perm,m_signlist,qlist)) {
    return 0;
    //std::cout<<"Error in MHVCalculator::Differential: Wrong flavors or helicities for quarks"<<std::endl;
    //abort();
  }
    
  if (m_qlist[0]==2) {
    if (nh==2)  return Elementary_MHVQ2_Amplitude(perm,m_signlist,qlist,n_part);
    if (ph==2)  return Elementary_MHVQ2bar_Amplitude(perm,m_signlist,qlist,n_part);
     
    if (nh==3)  return NMHVQ2_Amplitude(perm,m_signlist,qlist,n_part,nh);
    if (ph==3)  return NMHVQ2bar_Amplitude(perm,m_signlist,qlist,n_part,ph);
      
    if (nh==4) 	return NNMHVQ2_Amplitude(perm,m_signlist,qlist,n_part,nh);
    if (ph==4) 	return NNMHVQ2bar_Amplitude(perm,m_signlist,qlist,n_part,ph);
  }  
   
    
  if (m_qlist[0]==4) {
    if (nh==2)  return Elementary_MHVQ4_Amplitude(perm,m_signlist,qlist,n_part);
    if (ph==2)  return Elementary_MHVQ4bar_Amplitude(perm,m_signlist,qlist,n_part);
       
    if (nh==3)  return NMHVQ4_Amplitude(perm,m_signlist,qlist,n_part,nh);
    if (ph==3)  return NMHVQ4bar_Amplitude(perm,m_signlist,qlist,n_part,ph);

    if (nh==4)  return NNMHVQ4_Amplitude(perm,m_signlist,qlist,n_part,nh);
    if (ph==4)  return NNMHVQ4bar_Amplitude(perm,m_signlist,qlist,n_part,ph);
  }     
  return 0;
}



// private functions


void MHVCalculator::Make_Qlist(int* perm,int* plist,int* qlist,int part) 
{ 
    int nq=0;
    int qpos[4];
    for (int i=0;i<part;i++) {
	if ( !(plist[perm[i]]/20)  && (plist[perm[i]]!=0)) {
	    nq++;
	    qlist[nq]=i;	    
	    qpos[nq-1]=plist[perm[i]];
	}
	else  {
	    if ((plist[perm[i]]-21)*(plist[perm[i]]+21)*(plist[perm[i]]-25)*(plist[perm[i]]+25))  {
		std::cout<<"Error in MHVCalculator::Make_Qlist: Amplitude not implemented!"<<std::endl;
		abort();  
	    }
	} 
	if (nq>4) {
	    std::cout<<"Error in MHVCalculator::Make_Qlist: Too many quarks"<<std::endl;
	    abort();
	}
    }
    qlist[0]=nq; 
    if (nq<3) for (int i=0;i<nq;i++) qlist[i+3]=qpos[i];
    else if (nq<5) for (int i=0;i<nq;i++) qlist[i+5]=qpos[i];
    return;
}



bool MHVCalculator::Check_Qlist(int* perm,int* signlist,int* qlist) 
{
  if (qlist[0]==2) {
    if ((signlist[qlist[1]]==-signlist[qlist[2]]) && (qlist[3]==-qlist[4])) return 1;
  }
  if (qlist[0]==4) {
    int i;
    for (i=2;i<5;i++) if (qlist[i+4]==qlist[5]) { 
      THROW(fatal_error,"MHVCalculator::Check_Qlist: Amplitude with 2 identical pairs of quarks is not implemented");
    }
    for (i=2;signlist[qlist[1]]+signlist[qlist[i]] || qlist[5]+qlist[i+4];i++) if (i>3) return 0;  
    int sh=signlist[qlist[1]], sf=qlist[5];
    for (i=2;i<5;i++) {
      sh+=signlist[qlist[i]];
      sf+=qlist[i+4];
    }
    if (sh==0 && sf==0 ) return 1;
    return 0;
  }
  return 0;
}




///////////////////////////////////////////////////////////////////////////////////////////

Complex MHVCalculator::Elementary_MHV_Amplitude(int* perm,int* signlist,int part)
{
  int l;
  int v1=-1,v2=-1;
  for (l=0;l<part-1 && v1<0;l++) if (signlist[l]==-1) v1=perm[l];
  for (;l<part && v2<0;l++) if (signlist[l]==-1) v2=perm[l]; 
  if (v2<0) abort();
  Complex sm = p_BS->S0(v1,v2);
  sm = pow(sm,4);
  for (l=0;l<part-1;l++) sm/=p_BS->S0(perm[l],perm[l+1]);
  sm/=p_BS->S0(perm[part-1],perm[0]);

  return sm;
}


Complex MHVCalculator::NMHV_Amplitude(int* perm,int* signlist,int part,int vhel)
{ 
  if (vhel==2) return Elementary_MHV_Amplitude(perm,signlist,part);

  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_dummyarg[i+part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_dummysl[i+part] = m_dummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2) {
	switch (cs) {
	case 1:
	  m_dummysl[j]=-1;
	  m_dummysl[part+i]=1;
	  break;
	case 2:
	  m_dummysl[j]=1;
	  m_dummysl[part+i]=-1;
	} 
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i]; 
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	Complex v = Elementary_MHV_Amplitude(&m_dummyarg[i],&m_dummysl[i],j-i+1); 
	m_dummyarg[j] = perm[j];
	m_dummysl[j] = signlist[j];

	m_dummyarg[part+i] = pn;
	Complex vp = Elementary_MHV_Amplitude(&m_dummyarg[j],&m_dummysl[j],part-j+i+1);
	v*=vp;
	m_dummyarg[part+i] = perm[i];	
	m_dummysl[part+i] = signlist[i];
	v /= (p_BS->Momentum(pn)).Abs2();
	amp+=v;
      } 
    }
  }

  return amp;
}



Complex MHVCalculator::NNMHV_Amplitude(int* perm,int* signlist,int part,int vhel)
{
  if (vhel<4) return NMHV_Amplitude(perm,signlist,part,vhel);

  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_ndummyarg[i+part] = m_ndummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_ndummysl[i+part] = m_ndummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2 || cs==3) {
	for (int l=-1;l<0 || (l<2 && cs==2);l+=2) {
	  int vh=0;
	  switch (cs) {
	  case 1:
	    m_ndummysl[j]=-1;
	    m_ndummysl[part+i]=1;
	    vh=2;
	    break;
	  case 2:
	    m_ndummysl[j]=l;
	    m_ndummysl[part+i]=-l;
	    vh=(5-l)/2;
	    break;
	  case 3:
	    m_ndummysl[j]=1;
	    m_ndummysl[part+i]=-1;
	    vh=3;
	  }	
	  Pfunc pf(j-i+1);
	  for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	  int pn = p_BS->GetMomNumber(&pf);
	  m_ndummyarg[j] = pn;
	  Complex v = NMHV_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],j-i+1,vh);
	  m_ndummyarg[j] = perm[j];
	  m_ndummysl[j] = signlist[j];

	  m_ndummyarg[part+i] = pn;
	  v *= NMHV_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],part-j+i+1,5-vh);
	  m_ndummyarg[part+i] = perm[i];
	  m_ndummysl[part+i] = signlist[i];
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
      } 
    }    
  }

  amp/=2;
  return amp;
}






Complex MHVCalculator::Elementary_MHVbar_Amplitude(int* perm,int* signlist,int part)
{
  int l;
  int v1=-1,v2=-1;
  for (l=0;l<part-1 && v1<0;l++) if (signlist[l]==1) v1=perm[l];
  for (;l<part && v2<0;l++) if (signlist[l]==1) v2=perm[l];
  if (v2<0) abort();
  Complex sm = p_BS->S1(v1,v2);
  sm = pow(sm,4);
  for (l=0;l<part-1;l++)  sm/=p_BS->S1(perm[l],perm[l+1]);
  sm/=p_BS->S1(perm[part-1],perm[0]);

  return sm;
}


Complex MHVCalculator::NMHVbar_Amplitude(int* perm,int* signlist,int part,int vhel)
{
  if (vhel==2) return Elementary_MHVbar_Amplitude(perm,signlist,part); 

  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_dummyarg[i+part] = m_dummyarg[i] = perm[i]; 
  for (int i=0;i<part;i++) m_dummysl[i+part] = m_dummysl[i] = signlist[i]; 
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2) {
	switch (cs) {
	case 1:
	  m_dummysl[j]=1;
	  m_dummysl[part+i]=-1;
	  break;
	case 2:
	  m_dummysl[j]=-1;
	  m_dummysl[part+i]=1;
	} 	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	Complex v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],&m_dummysl[i],j-i+1);
	m_dummyarg[j] = perm[j];
	m_dummysl[j] = signlist[j];

	m_dummyarg[part+i] = pn;
	v *= Elementary_MHVbar_Amplitude(&m_dummyarg[j],&m_dummysl[j],part-j+i+1);
	m_dummyarg[part+i] = perm[i];
	m_dummysl[part+i] = signlist[i];
	v /= (p_BS->Momentum(pn)).Abs2();
	amp+=v;
      } 
    }    
  }

  return amp;
}



Complex MHVCalculator::NNMHVbar_Amplitude(int* perm,int* signlist,int part,int vhel)
{
  if (vhel<4) return NMHVbar_Amplitude(perm,signlist,part,vhel);
  
  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_ndummyarg[i+part] = m_ndummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_ndummysl[i+part] = m_ndummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2 || cs==3) {
	for (int l=-1;l<0 || (l<2 && cs==2);l+=2) {
	  int vh=0;
	  switch (cs) {
	  case 1:
	    m_ndummysl[j]=1;
	    m_ndummysl[part+i]=-1;
	    vh=2;
	    break;
	  case 2:
	    m_ndummysl[j]=-l;
	    m_ndummysl[part+i]=l;
	    vh=(5-l)/2;
	    break;
	  case 3:
	    m_ndummysl[j]=-1;
	    m_ndummysl[part+i]=1;
	    vh=3;
	  }	
	  Pfunc pf(j-i+1);
	  for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	  int pn = p_BS->GetMomNumber(&pf);
	  m_ndummyarg[j] = pn;
	  Complex v = NMHVbar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],j-i+1,vh);
	  m_ndummyarg[j] = perm[j];
	  m_ndummysl[j] = signlist[j];

	  m_ndummyarg[part+i] = pn; 
	  v *= NMHVbar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],part-j+i+1,5-vh); 
	  m_ndummyarg[part+i] = perm[i];
	  m_ndummysl[part+i] = signlist[i];
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
      } 
    }    
  }
  
  amp/=2;
  return amp;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


Complex MHVCalculator::Elementary_MHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part)
{
    int l;
    int v1=-1;  
    for (l=0;l<part && v1<0;l++) if (signlist[l]==-1 && l!=qlist[1] && l!=qlist[2]) v1=perm[l]; 
    if (v1<0) abort();
    Complex sm = p_BS->S0(v1,perm[qlist[1]]);
    if (signlist[qlist[1]]==-1) sm = pow(sm,3);
    Complex sm2 = p_BS->S0(v1,perm[qlist[2]]);
    if (signlist[qlist[2]]==-1) sm2 = -pow(sm2,3);
    sm *= sm2;  
   
    for (l=0;l<part-1;l++)  sm/=p_BS->S0(perm[l],perm[l+1]);
    sm/=p_BS->S0(perm[part-1],perm[0]);
        
    return sm;
}




Complex MHVCalculator::NMHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel)
{
  if (vhel==2) return Elementary_MHVQ2_Amplitude(perm,signlist,qlist,part);

  Complex amp(0.,0.);  
  for (int i=0;i<part;i++) m_dummyarg[i+part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_dummysl[i+part] = m_dummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2) {
	switch (cs) {
	case 1:
	  m_dummysl[j]=-1;
	  m_dummysl[part+i]=1;
	  break;
	case 2:
	  m_dummysl[j]=1;
	  m_dummysl[part+i]=-1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	
	int qlist1[5];
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	bool empty(false);
	Complex v;
	// pure gluonic amilitude 
	if (!qlist1[0]) v = Elementary_MHV_Amplitude(&m_dummyarg[i],&m_dummysl[i],j-i+1);
	
	// amplitudes with quarks
	else if (qlist1[0]==2) {
	  v = -m_dummysl[i+qlist1[1]];
	  v *= Elementary_MHVQ2_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && m_dummysl[i+qlist1[1]]==-m_dummysl[j]) {
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else empty=true;	
	m_dummyarg[j] = perm[j];
	m_dummysl[j] = signlist[j];

	if (!empty) {
	  m_dummyarg[part+i] = pn;
	  
	  int qlist2[5]; 
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,part-j+i);
	  Complex vp;
	  // pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHV_Amplitude(&m_dummyarg[j],&m_dummysl[j],part-j+i+1);
	  
	  // amplitudes with quarks
	  else if (qlist2[0]==2)   {
	    vp = -m_dummysl[j+qlist2[1]];
	    vp *= Elementary_MHVQ2_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  else if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=part-j+i;
	    vp = Elementary_MHVQ2_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  v*=vp;
	  m_dummyarg[part+i] = perm[i];
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
	m_dummysl[part+i] = signlist[i];
      } 
    }    
  }
  if (signlist[qlist[1]]>0) amp= -amp;
  return amp;
}



Complex MHVCalculator::NNMHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel)
{ 
  if (vhel<4) return NMHVQ2_Amplitude(perm,signlist,qlist,part,vhel);

  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_ndummyarg[i+part] = m_ndummyarg[i] = perm[i]; 
  for (int i=0;i<part;i++) m_ndummysl[i+part] = m_ndummysl[i] = signlist[i]; 
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2 || cs==3) {
	for (int l=-1;l<0 || (l<2 && cs==2);l+=2) {
	  int vh=0; 
	  switch (cs) {
	  case 1:
	    m_ndummysl[j]=-1;
	    m_ndummysl[part+i]=1;
	    vh=2;
	    break;
	  case 2:
	    m_ndummysl[j]=l;
	    m_ndummysl[part+i]=-l;
	    vh=(5-l)/2;
	    break;
	  case 3:
	    m_ndummysl[j]=1;
	    m_ndummysl[part+i]=-1;
	    vh=3;
	  }			
	  Pfunc pf(j-i+1);	
	  for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	  int pn = p_BS->GetMomNumber(&pf);
	  m_ndummyarg[j] = pn;

	  int qlist1[5];
	  Make_Qlist(&m_ndummyarg[i],m_plist,qlist1,j-i);
	  bool empty(false);
	  Complex v;

	  // pure gluonic amilitude 
	  if (!qlist1[0])   v = NMHV_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],j-i+1,vh);
	  
	  // amplitudes with quarks
	  else if (qlist1[0]==2)   {
	    m_plist[pn]=21;
	    v = -m_ndummysl[i+qlist1[1]];
	    v *= NMHVQ2_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==1 && m_ndummysl[i+qlist1[1]]==-m_ndummysl[j]) {
	    m_plist[pn]=-qlist1[3];
	    qlist1[0]=2;
	    qlist1[2]=j-i;
	    v = -NMHVQ2_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else empty=true;	  
	  m_ndummyarg[j] = perm[j]; 
	  m_ndummysl[j] = signlist[j];

	  if (!empty) {
	    m_ndummyarg[part+i] = pn;

	    int qlist2[5];
	    Make_Qlist(&m_ndummyarg[j],m_plist,qlist2,part-j+i);
	    Complex vp;
	    // pure gluonic amilitude 
	    if (!qlist2[0]) vp = NMHV_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],part-j+i+1,5-vh);
	    
	    // amplitudes with quarks
	    else if (qlist2[0]==2)   {
	      m_plist[pn]=21;
	      vp = -m_ndummysl[j+qlist2[1]];
	      vp *= NMHVQ2_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    else if (qlist2[0]==1) {
	      m_plist[pn]=-qlist2[3];
	      qlist2[0]=2;
	      qlist2[2]=part-j+i;
	      vp = NMHVQ2_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    v*=vp;
	    m_ndummyarg[part+i] = perm[i];
	    v /= (p_BS->Momentum(pn)).Abs2();
	    amp+=v;
	  }
	  m_ndummysl[part+i] = signlist[i];
	}
      } 
    }    
  }

  amp/=2;
  if (signlist[qlist[1]]>0) amp= -amp;
  return amp;
}






Complex MHVCalculator::Elementary_MHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part)
{   
    int l;
    int v1=-1;  
    for (l=0;l<part && v1<0;l++) if (signlist[l]==1 && l!=qlist[1] && l!=qlist[2]) v1=perm[l];
    if (v1<0) abort();  
    Complex sm = p_BS->S1(v1,perm[qlist[1]]);
    if (signlist[qlist[1]]==1) sm = pow(sm,3); 
    Complex sm2 = p_BS->S1(v1,perm[qlist[2]]);
    if (signlist[qlist[2]]==1) sm2 = -pow(sm2,3);
    sm *= sm2;  
    
    for (l=0;l<part-1;l++)  sm/=p_BS->S1(perm[l],perm[l+1]);
    sm/=p_BS->S1(perm[part-1],perm[0]); 
    
    return sm;
}



Complex MHVCalculator::NMHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel)
{
  if (vhel==2) return Elementary_MHVQ2bar_Amplitude(perm,signlist,qlist,part);
  
  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_dummyarg[i+part] = m_dummyarg[i] = perm[i]; 
  for (int i=0;i<part;i++) m_dummysl[i+part] = m_dummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2) {
	switch (cs) {
	case 1:
	  m_dummysl[j]=1;
	  m_dummysl[part+i]=-1;
	  break;
	case 2:
	  m_dummysl[j]=-1;
	  m_dummysl[part+i]=1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	
	int qlist1[5];
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	bool empty(false);
	Complex v;
	// pure gluonic amilitude 
	if (!qlist1[0])  v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],&m_dummysl[i],j-i+1);
	
	// amplitudes with quarks
	else if (qlist1[0]==2)   {
	  v = m_dummysl[i+qlist1[1]];
	  v *= Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && m_dummysl[i+qlist1[1]]==-m_dummysl[j]) {
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else empty=true;	
	m_dummyarg[j] = perm[j];
	m_dummysl[j] = signlist[j];

	if (!empty) {
	  m_dummyarg[part+i] = pn;
	  
	  int qlist2[5];
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,part-j+i);
	  Complex vp;
	  // pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHVbar_Amplitude(&m_dummyarg[j],&m_dummysl[j],part-j+i+1);
	  
	  // amplitudes with quarks
	  else if (qlist2[0]==2)   {	
	    vp = m_dummysl[j+qlist2[1]];
	    vp *= Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  else if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=part-j+i;
	    vp = Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  v*=vp;
	  m_dummyarg[part+i] = perm[i];
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
	m_dummysl[part+i] = signlist[i];
      } 
    }
  }
  if (signlist[qlist[1]]<0) amp=-amp;
  return amp;
}



Complex MHVCalculator::NNMHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel)
{ 
  if (vhel<4) return NMHVQ2bar_Amplitude(perm,signlist,qlist,part,vhel);

  Complex amp(0.,0.);
  for (int i=0;i<part;i++) m_ndummyarg[i+part] = m_ndummyarg[i] = perm[i]; 
  for (int i=0;i<part;i++) m_ndummysl[i+part] = m_ndummysl[i] = signlist[i]; 
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2 || cs==3) {
	for (int l=-1;l<0 || (l<2 && cs==2);l+=2) {
	  int vh=0; 
	  switch (cs) {
	  case 1:
	    m_ndummysl[j]=1;
	    m_ndummysl[part+i]=-1;
	    vh=2;
	    break;
	  case 2:
	    m_ndummysl[j]=-l;
	    m_ndummysl[part+i]=l;
	    vh=(5-l)/2;
	    break;
	  case 3:
	    m_ndummysl[j]=-1;
	    m_ndummysl[part+i]=1;
	    vh=3;
	  }			
	  Pfunc pf(j-i+1);	
	  for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	  int pn = p_BS->GetMomNumber(&pf);
	  m_ndummyarg[j] = pn;

	  int qlist1[5];
	  Make_Qlist(&m_ndummyarg[i],m_plist,qlist1,j-i);
	  bool empty(false);
	  Complex v;

	  // pure gluonic amilitude 
	  if (!qlist1[0])   v = NMHVbar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],j-i+1,vh);
	  
	  // amplitudes with quarks
	  else if (qlist1[0]==2)   {
	    m_plist[pn]=21;
	    v = m_ndummysl[i+qlist1[1]];
	    v *= NMHVQ2bar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==1 && m_ndummysl[i+qlist1[1]]==-m_ndummysl[j]) {
	    m_plist[pn]=-qlist1[3];
	    qlist1[0]=2;
	    qlist1[2]=j-i;
	    v = -NMHVQ2bar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else empty=true;
	  m_ndummyarg[j] = perm[j];
	  m_ndummysl[j] = signlist[j];

	  if (!empty) {
	    m_ndummyarg[part+i] = pn;

	    int qlist2[5];
	    Make_Qlist(&m_ndummyarg[j],m_plist,qlist2,part-j+i);
	    Complex vp;
	    // pure gluonic amilitude 
	    if (!qlist2[0]) vp = NMHVbar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],part-j+i+1,5-vh);
	    
	    // amplitudes with quarks
	    else if (qlist2[0]==2)   {
	      m_plist[pn]=21;
	      vp = m_ndummysl[j+qlist2[1]];
	      vp *= NMHVQ2bar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    else if (qlist2[0]==1) {
	      m_plist[pn]=-qlist2[3];
	      qlist2[0]=2;
	      qlist2[2]=part-j+i;
	      vp = NMHVQ2bar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    v*=vp;
	    m_ndummyarg[part+i] = perm[i];
	    v /= (p_BS->Momentum(pn)).Abs2();
	    amp+=v;
	  }
	  m_ndummysl[part+i] = signlist[i];
	}
      } 
    }    
  }

  amp/=2;
  if (signlist[qlist[1]]<0) amp= -amp;
  return amp;
}






/////////////////////////////////////////////////////////////////////////////////////////////////////

Complex MHVCalculator::Elementary_MHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part)
{
    int v1(-1), v2(-1);
    for (int i=1;i<5;i++) {
	if (v1>-1 && signlist[qlist[i]]==-1) v2=i;
	if (v1<0 && signlist[qlist[i]]==-1)  v1=i;
    }
    
    Complex sm(p_BS->S0(perm[qlist[v1]],perm[qlist[v2]]));
    sm = pow(sm,2);
    if ((v2-v1)%2)  sm= -sm;
    if (qlist[5]>0) sm= -sm;
    
    if (qlist[6]==-qlist[5]) {  
	sm *= p_BS->S0(perm[qlist[1]],perm[qlist[4]]);
	sm *= p_BS->S0(perm[qlist[3]],perm[qlist[2]]);
    }
    else {
	sm *= p_BS->S0(perm[qlist[1]],perm[qlist[2]]);
	sm *= p_BS->S0(perm[qlist[3]],perm[qlist[4]]);
    }
    
    for (int l=0;l<part-1;l++) {
	sm/=p_BS->S0(perm[l],perm[l+1]);
    } 
    sm/=p_BS->S0(perm[part-1],perm[0]);
    
    return sm;
}



 
Complex MHVCalculator::NMHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part, int vhel)
{ 
  if (vhel==2) return Elementary_MHVQ4_Amplitude(perm,signlist,qlist,part);
  
  Complex amp(0.,0.);  
  for (int i=0;i<part;i++) m_dummyarg[i+part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_dummysl[i+part] = m_dummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2) {
	switch (cs) {
	case 1:
	  m_dummysl[j]=-1;
	  m_dummysl[part+i]=1;
	  break;
	case 2:
	  m_dummysl[j]=1;
	  m_dummysl[part+i]=-1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	
	int qlist1[9];
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	Complex v;
	bool empty(false);
	int qsign=0;
	for (int l=1;l<qlist1[0]+1;l++) qsign+= m_dummysl[i+qlist1[l]];
	// pure gluonic amilitude 
	if (!qlist1[0])  v = Elementary_MHV_Amplitude(&m_dummyarg[i],&m_dummysl[i],j-i+1);
	
	// amplitudes with quarks
	else if (qlist1[0]==4) v = Elementary_MHVQ4_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	else if (qlist1[0]==3 && qsign==-m_dummysl[j]) {
	  qlist1[0]=4;
	  qlist1[4]=j-i;
	  v = -Elementary_MHVQ4_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else if (qlist1[0]==2 && qsign==0 && qlist1[3]==-qlist1[4]) {
	  v = -Elementary_MHVQ2_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && qsign==-m_dummysl[j]) {
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else empty=true;
	m_dummyarg[j] = perm[j];
	m_dummysl[j] = signlist[j];

	if (!empty) {
	  Complex vp;
	  m_dummyarg[part+i] = pn;
	  
	  int qlist2[9];
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,part-j+i);
	  // pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHV_Amplitude(&m_dummyarg[j],&m_dummysl[j],part-j+i+1);
	  
	  
	  // amplitudes with quarks
	  else if (qlist2[0]==4) vp = Elementary_MHVQ4_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  else if (qlist2[0]==3) {
	    qlist2[0]=4;
	    qlist2[4]=part-j+i;
	    vp = -m_dummysl[part+i];
	    vp *= Elementary_MHVQ4_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	    if (m_plist[m_dummyarg[j+qlist2[3]]]<0) vp*= -1;
	  }
	  else if (qlist2[0]==2 && qlist1[3]==-qlist1[4])   {
	    vp = Elementary_MHVQ2_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  else if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=part-j+i; 
	    vp = -m_dummysl[part+i];
	    vp *= Elementary_MHVQ2_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	    if (m_plist[m_dummyarg[j+qlist2[1]]]<0) vp*= -1;
	  }
	  v*=vp;
	  m_dummyarg[part+i] = perm[i];
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
	m_dummysl[part+i] = signlist[i];
      } 
    } 
  }

  if (qlist[5]>0) amp=-amp;
  return amp;
}




Complex MHVCalculator::NNMHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part, int vhel)
{ 
  if (vhel<4) return NMHVQ4_Amplitude(perm,signlist,qlist,part,vhel);

  Complex amp(0.,0.);  
  for (int i=0;i<part;i++) m_ndummyarg[i+part] = m_ndummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_ndummysl[i+part] = m_ndummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,-1);
      if (cs==1 || cs==2 || cs==3) {
	for (int l=-1;l<0 || (l<2 && cs==2);l+=2) {
	  int vh=0; 
	  switch (cs) {
	  case 1:
	    m_ndummysl[j]=-1;
	    m_ndummysl[part+i]=1;
	    vh=2;
	    break;
	  case 2:
	    m_ndummysl[j]=l;
	    m_ndummysl[part+i]=-l;
	    vh=(5-l)/2;
	    break;
	  case 3:
	    m_ndummysl[j]=1;
	    m_ndummysl[part+i]=-1;
	    vh=3;
	  }		
	  Pfunc pf(j-i+1);
	  for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	  int pn = p_BS->GetMomNumber(&pf);
	  m_ndummyarg[j] = pn;

	  int qlist1[9];
	  Make_Qlist(&m_ndummyarg[i],m_plist,qlist1,j-i);
	  Complex v;
	  bool empty(false);
	  int qsign=0;
	  for (int l=1;l<qlist1[0]+1;l++) qsign+= m_ndummysl[i+qlist1[l]];
	  // pure gluonic amilitude 
	  if (!qlist1[0]) v = NMHV_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],j-i+1,vh);
	  
	  // amplitudes with quarks
	  else if (qlist1[0]==4)   {
	    m_plist[pn]=21;
	    v = NMHVQ4_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==3 && qsign==-m_ndummysl[j]) {
	    if (qlist1[5]==-qlist1[6]) m_plist[pn]=-qlist1[7];
	    else m_plist[pn]=-qlist1[5];
	    qlist1[0]=4;
	    qlist1[4]=j-i;
	    v = -NMHVQ4_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==2 && qsign==0 && qlist1[3]==-qlist1[4])   {
	    m_plist[pn]=21;
	    v = -NMHVQ2_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==1 && qsign==-m_ndummysl[j]) {
	    m_plist[pn]=-qlist1[3];
	    qlist1[0]=2;
	    qlist1[2]=j-i;
	    v = -NMHVQ2_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else empty=true;
	  m_ndummyarg[j] = perm[j];
	  m_ndummysl[j] = signlist[j];

	  if (!empty) {
	    Complex vp;
	    m_ndummyarg[part+i] = pn;	  
	    int qlist2[9];
	    Make_Qlist(&m_ndummyarg[j],m_plist,qlist2,part-j+i);
	    
	    // pure gluonic amilitude 
	    if (!qlist2[0]) vp = NMHV_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],part-j+i+1,5-vh);	  
	    
	    // amplitudes with quarks
	    else if (qlist2[0]==4)  {
	      m_plist[pn]=21;
	      vp = NMHVQ4_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    else if (qlist2[0]==3) {
	      if (qlist2[5]==-qlist2[6]) m_plist[pn]=-qlist2[7];
	      else m_plist[pn]=-qlist2[5];
	      qlist2[0]=4;
	      qlist2[4]=part-j+i;
	      vp = -m_ndummysl[part+i];
	      vp *= NMHVQ4_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	      if (m_plist[m_ndummyarg[j+qlist2[3]]]<0) vp*= -1;
	    }
	    else if (qlist2[0]==2 && qlist1[3]==-qlist1[4])   {
	      m_plist[pn]=21;
	      vp = NMHVQ2_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    else if (qlist2[0]==1) {
	      m_plist[pn]=-qlist2[3];
	      qlist2[0]=2;
	      qlist2[2]=part-j+i; 
	      vp = -m_ndummysl[part+i];
	      vp *= NMHVQ2_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	      if (m_plist[m_ndummyarg[j+qlist2[1]]]<0) vp*= -1;
	    }
	    v*=vp;
	    m_ndummyarg[part+i] = perm[i];
	    v /= (p_BS->Momentum(pn)).Abs2();
	    amp+=v;
	  } 
	  m_ndummysl[part+i] = signlist[i];
	}
      } 
    } 
  }
  
  amp/=2;
  if (qlist[5]>0) amp=-amp;
  return amp;
}




Complex MHVCalculator::Elementary_MHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part)
{
    int v1(-1), v2(-1);
    for (int i=1;i<5;i++) {
	if (v1>-1 && signlist[qlist[i]]==1) v2=i;
	if (v1<0 && signlist[qlist[i]]==1) v1=i;
    }
    
    Complex sm(p_BS->S1(perm[qlist[v1]],perm[qlist[v2]]));
    sm = pow(sm,2);
    if ((v2-v1)%2)  sm= -sm;
    if (qlist[5]>0) sm= -sm;
    
    if (qlist[6]==-qlist[5]) {  
	sm *= p_BS->S1(perm[qlist[1]],perm[qlist[4]]);
	sm *= p_BS->S1(perm[qlist[3]],perm[qlist[2]]);
    }
    else {
	sm *= p_BS->S1(perm[qlist[1]],perm[qlist[2]]);
	sm *= p_BS->S1(perm[qlist[3]],perm[qlist[4]]);
    }
    
    for (int l=0;l<part-1;l++) {
	sm/=p_BS->S1(perm[l],perm[l+1]);
    } 
    sm/=p_BS->S1(perm[part-1],perm[0]); 
    
    return sm;
}




Complex MHVCalculator::NMHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel)
{ 
  if (vhel==2) return Elementary_MHVQ4bar_Amplitude(perm,signlist,qlist,part);
  
  Complex amp(0.,0.);  
  for (int i=0;i<part;i++) m_dummyarg[i+part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_dummysl[i+part] = m_dummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2) {
	switch (cs) {
	case 1:
	  m_dummysl[j]=1;
	  m_dummysl[part+i]=-1;
	  break;
	case 2:
	  m_dummysl[j]=-1;
	  m_dummysl[part+i]=1;
	}	
	Pfunc pf(j-i+1);
	for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	int pn = p_BS->GetMomNumber(&pf);
	m_dummyarg[j] = pn;
	
	int qlist1[9];
	Make_Qlist(&m_dummyarg[i],m_plist,qlist1,j-i);
	Complex v;
	bool empty(false);
	int qsign=0;
	for (int l=1;l<qlist1[0]+1;l++) qsign+= m_dummysl[i+qlist1[l]]; 
	// pure gluonic amilitude 
	if (!qlist1[0]) v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],&m_dummysl[i],j-i+1);
	
	// amplitudes with quarks
	else if (qlist1[0]==4) v = Elementary_MHVQ4bar_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	else if (qlist1[0]==3 && qsign==-m_dummysl[j]) {
	  qlist1[0]=4;
	  qlist1[4]=j-i;
	  v = -Elementary_MHVQ4bar_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else if (qlist1[0]==2 && qsign==0 && qlist1[3]==-qlist1[4]) {
	  v = -Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);
	}
	else if (qlist1[0]==1 && qsign==-m_dummysl[j]) {
	  qlist1[0]=2;
	  qlist1[2]=j-i;
	  v = -Elementary_MHVQ2bar_Amplitude(&m_dummyarg[i],&m_dummysl[i],qlist1,j-i+1);	  
	}
	else empty=true;
	m_dummyarg[j] = perm[j];
	m_dummysl[j] = signlist[j];

	if (!empty) {
	  Complex vp;
	  m_dummyarg[part+i] = pn;
	  
	  int qlist2[9];
	  Make_Qlist(&m_dummyarg[j],m_plist,qlist2,part-j+i);
	  // pure gluonic amilitude 
	  if (!qlist2[0]) vp = Elementary_MHVbar_Amplitude(&m_dummyarg[j],&m_dummysl[j],part-j+i+1);	  
	  
	  // amplitudes with quarks
	  else if (qlist2[0]==4)  {
	    vp = Elementary_MHVQ4bar_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  else if (qlist2[0]==3) {
	    qlist2[0]=4;
	    qlist2[4]=part-j+i;
	    vp = m_dummysl[part+i];
	    vp*= Elementary_MHVQ4bar_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	    if (m_plist[m_dummyarg[j+qlist2[3]]]<0) vp*= -1;
	  }
	  else if (qlist2[0]==2 && qlist2[3]==-qlist2[4])   {
	    vp = Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	  }
	  else if (qlist2[0]==1) {
	    qlist2[0]=2;
	    qlist2[2]=part-j+i; 
	    vp = m_dummysl[part+i];
	    vp*= Elementary_MHVQ2bar_Amplitude(&m_dummyarg[j],&m_dummysl[j],qlist2,part-j+i+1);
	    if (m_plist[m_dummyarg[j+qlist2[1]]]<0) vp*= -1;
	  }
	  v*=vp;
	  m_dummyarg[part+i] = perm[i];
	  v /= (p_BS->Momentum(pn)).Abs2();
	  amp+=v;
	}
	m_dummysl[part+i] = signlist[i];
      } 
    } 
  }

  if (qlist[5]>0) amp=-amp;
  return amp;
}



Complex MHVCalculator::NNMHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part, int vhel)
{ 
  if (vhel<4) return NMHVQ4bar_Amplitude(perm,signlist,qlist,part,vhel);

  Complex amp(0.,0.);  
  for (int i=0;i<part;i++) m_ndummyarg[i+part] = m_ndummyarg[i] = perm[i];
  for (int i=0;i<part;i++) m_ndummysl[i+part] = m_ndummysl[i] = signlist[i];
  for (int i=0;i<part-2;i++) {
    for (int j=i+2;j<part && j-i+1<part;j++) {
      int cs=CountSign(signlist,i,j,1);
      if (cs==1 || cs==2 || cs==3) {
	for (int l=-1;l<0 || (l<2 && cs==2);l+=2) {
	  int vh=0; 
	  switch (cs) {
	  case 1:
	    m_ndummysl[j]=1;
	    m_ndummysl[part+i]=-1;
	    vh=2;
	    break;
	  case 2:
	    m_ndummysl[j]=-l;
	    m_ndummysl[part+i]=l;
	    vh=(5-l)/2;
	    break;
	  case 3:
	    m_ndummysl[j]=-1;
	    m_ndummysl[part+i]=1;
	    vh=3;
	  }		
	  Pfunc pf(j-i+1);
	  for (int k=0;k<j-i;k++) pf.arg[k+1]=perm[k+i];
	  int pn = p_BS->GetMomNumber(&pf);
	  m_ndummyarg[j] = pn;

	  int qlist1[9];
	  Make_Qlist(&m_ndummyarg[i],m_plist,qlist1,j-i);
	  Complex v;
	  bool empty(false);
	  int qsign=0;
	  for (int l=1;l<qlist1[0]+1;l++) qsign+= m_ndummysl[i+qlist1[l]];
	  // pure gluonic amilitude 
	  if (!qlist1[0]) v = NMHVbar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],j-i+1,vh);
	  
	  // amplitudes with quarks
	  else if (qlist1[0]==4)   {
	    m_plist[pn]=21;
	    v = NMHVQ4bar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==3 && qsign==-m_ndummysl[j]) {
	    if (qlist1[5]==-qlist1[6]) m_plist[pn]=-qlist1[7];
	    else m_plist[pn]=-qlist1[5];
	    qlist1[0]=4;
	    qlist1[4]=j-i;
	    v = -NMHVQ4bar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==2 && qsign==0 && qlist1[3]==-qlist1[4])   {
	    m_plist[pn]=21;
	    v = -NMHVQ2bar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else if (qlist1[0]==1 && qsign==-m_ndummysl[j]) {
	    m_plist[pn]=-qlist1[3];
	    qlist1[0]=2;
	    qlist1[2]=j-i;
	    v = -NMHVQ2bar_Amplitude(&m_ndummyarg[i],&m_ndummysl[i],qlist1,j-i+1,vh);
	  }
	  else empty=true;
	  m_ndummyarg[j] = perm[j];
	  m_ndummysl[j] = signlist[j];

	  if (!empty) {
	    Complex vp;
	    m_ndummyarg[part+i] = pn;	  
	    int qlist2[9];
	    Make_Qlist(&m_ndummyarg[j],m_plist,qlist2,part-j+i);
	    
	    // pure gluonic amilitude 
	    if (!qlist2[0]) vp = NMHVbar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],part-j+i+1,5-vh);	  
	    
	    // amplitudes with quarks
	    else if (qlist2[0]==4)  {
	      m_plist[pn]=21;
	      vp = NMHVQ4bar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    else if (qlist2[0]==3) {
	      if (qlist2[5]==-qlist2[6]) m_plist[pn]=-qlist2[7];
	      else m_plist[pn]=-qlist2[5];
	      qlist2[0]=4;
	      qlist2[4]=part-j+i;
	      vp = m_ndummysl[part+i];
	      vp *= NMHVQ4bar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	      if (m_plist[m_ndummyarg[j+qlist2[3]]]<0) vp*= -1;
	    }
	    else if (qlist2[0]==2 && qlist1[3]==-qlist1[4])   {
	      m_plist[pn]=21;
	      vp = NMHVQ2bar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	    }
	    else if (qlist2[0]==1) {
	      m_plist[pn]=-qlist2[3];
	      qlist2[0]=2;
	      qlist2[2]=part-j+i; 
	      vp = m_ndummysl[part+i];
	      vp *= NMHVQ2bar_Amplitude(&m_ndummyarg[j],&m_ndummysl[j],qlist2,part-j+i+1,5-vh);
	      if (m_plist[m_ndummyarg[j+qlist2[1]]]<0) vp*= -1;
	    }
	    v*=vp;
	    m_ndummyarg[part+i] = perm[i];
	    v /= (p_BS->Momentum(pn)).Abs2();
	    amp+=v;
	  }
	  m_ndummysl[part+i] = signlist[i];
	}
      } 
    } 
  }
  
  amp/=2;
  if (qlist[5]>0) amp=-amp;
  return amp;
}

