#include "MHV_Calculator.H"
#include "Pfunc.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

MHV_Calculator::MHV_Calculator(int part, Basic_Sfuncs* BS) :
  n_part(part), p_BS(BS)
{ 
  m_dummyarg = new int[2*part];
  for (int i=0;i<part;i++) {
    m_dummyarg[i]=i;
    m_dummyarg[part+i]=i;
  }
  m_arg = new int[part];
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
// 	std::cout<<pl.size()<<":";
// 	for (int k=1;k<pf->argnum;k++)std::cout<<" "<<pf->arg[k];
// 	std::cout<<std::endl;
	pl.push_back(pf);
      }
    }
    p_BS->BuildMomlist(pl);
    for (Pfunc_Iterator it=pl.begin();it!=pl.end();it++) delete (*it);
  }
//   std::cout<<"MHV_Calculator: constructor finalized"<<std::endl;
}

MHV_Calculator::~MHV_Calculator() 
{ 
  delete[] m_arg;
}

double MHV_Calculator::Differential(int* perm,int* signlist)
{
//   std::cout<<"in MHV_Calculator::Differential"<<std::endl;
  int sum = 0;
  for (int i=0;i<n_part;i++) sum+=signlist[i];
  for (int i=0;i<n_part-1;i++) m_arg[i]=perm[i];
  m_arg[n_part-1] = n_part-1;
  int nh=(n_part-sum)/2;
  int ph=n_part - nh;
//       std::cout<<"Signlist: "<<nh<<"/"<<ph<<std::endl;
//         for (int i=0;i<n_part;i++) std::cout<<" "<<signlist[i];
//          std::cout<<std::endl<<"Permutation:"<<std::endl;
//          for (int i=0;i<n_part-1;i++) std::cout<<" "<<perm[i];
//          std::cout<<std::endl;

  if (nh==2) return norm(Elementary_MHV_Amplitude(m_arg,signlist,n_part));
  if (ph==2) return norm(Elementary_MHVbar_Amplitude(m_arg,signlist,n_part));
  if (nh<2 || ph<2) return 0.;

  if (Min(nh,ph)>3){
    std::cout<<"Error in MHV_Calculator::Differential: Amplitude not implemented!"<<std::endl;
    abort();
  }
  if (nh<ph) return norm(MHV_Amplitude(m_arg,signlist,nh));
  return norm(MHVbar_Amplitude(m_arg,signlist,ph));
//    std::cout<<"norm: "<<tmv<<" "<<dc<<std::endl;
//   return Max(tmv,dc);
}



Complex MHV_Calculator::Elementary_MHV_Amplitude(int* perm,int* signlist,int part,int psign)
{
  int l;
  int v1=-1,v2=-1;
  if (psign) {
    for (l=0;l<part-1 && v1<0;l++) if (signlist[perm[l]]==-1) v1=perm[l];
    for (;l<part && v2<0;l++) if (signlist[perm[l]]==-1) v2=perm[l];
  }
  else {
    for (l=0;l<part-1 && v1<0;l++) if (signlist[l]==-1) v1=perm[l];
    for (;l<part && v2<0;l++) if (signlist[l]==-1) v2=perm[l];
  }
  if (v2<0) abort();

  int icnt=0;
  if (v1<part && p_BS->Sign(v1)<0) icnt+=4;
  if (v2<part && p_BS->Sign(v2)<0) icnt+=4;
  for (l=0;l<part;l++) if (perm[l]<part && p_BS->Sign(perm[l])<0) icnt-=2;
  Complex sm = p_BS->S0(v1,v2);
//   std::cout<<"MHV S("<<psign<<"): "<<v1<<" "<<v2<<" = "<<sm<<" "<<p_BS->Momentum(v1)<<p_BS->Momentum(v2)<<std::endl;
  sm = pow(sm,4);
  for (l=0;l<part-1;l++){ sm/=p_BS->S0(perm[l],perm[l+1]);
//    std::cout<<"MHV S("<<psign<<"): "<<perm[l]<<" "<<perm[l+1]<<" = "<<p_BS->S0(perm[l],perm[l+1])<<std::endl;
  }
  sm/=p_BS->S0(perm[part-1],perm[0]);
//   std::cout<<icnt<<std::endl;
//   sm*=pow(Complex(0.,1.),icnt);
//   if(psign) std::cout<<"MHV: "<<sm<<std::endl;
  return sm;
}


Complex MHV_Calculator::Elementary_MHVbar_Amplitude(int* perm,int* signlist,int part,int psign)
{
  int l;
  int v1=-1,v2=-1;
  if (psign) {
    for (l=0;l<part-1 && v1<0;l++) if (signlist[perm[l]]==1) v1=perm[l];
    for (;l<part && v2<0;l++) if (signlist[perm[l]]==1) v2=perm[l];
  }
  else {
    for (l=0;l<part-1 && v1<0;l++) if (signlist[l]==1) v1=perm[l];
    for (;l<part && v2<0;l++) if (signlist[l]==1) v2=perm[l];
  }
  if (v2<0) abort();

  int icnt=0;
  if (v1<part && p_BS->Sign(v1)<0) icnt+=4;
  if (v2<part && p_BS->Sign(v2)<0) icnt+=4;
  for (l=0;l<part;l++) if (perm[l]<part && p_BS->Sign(perm[l])<0) icnt-=2;
  Complex sm = p_BS->S1(v1,v2);
  sm = pow(sm,4);
  for (l=0;l<part-1;l++) sm/=p_BS->S1(perm[l],perm[l+1]);
  sm/=p_BS->S1(perm[part-1],perm[0]);
//   std::cout<<icnt<<std::endl;
//   sm*=pow(Complex(0.,1.),icnt);
//   if(psign) std::cout<<"MHV: "<<sm<<std::endl;
  return sm;
}

int MHV_Calculator::CountSign(int* perm,int* signlist,int part,int sg)
{
  int sum = 0;
  for (int i=0;i<part;i++) sum+=signlist[perm[i]];
  if (sg==-1) return (part-sum)/2;
  return (part+sum)/2;
}

Complex MHV_Calculator::MHV_Amplitude(int* perm,int* signlist,int vhel)
{
  if (vhel==2) return Elementary_MHV_Amplitude(perm,signlist,n_part);
  
  Complex amp(0.,0.);
//   std::cout<<"in MHV_Calculator::MHV_Amplitude"<<std::endl;
  
  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(&perm[i],signlist,j-i,-1);
      if (cs==1 || cs==2) {
 //          std::cout<<i<<"/"<<j<<":   "<<cs<<std::endl;
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[perm[i+k]];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[perm[(j+k)%n_part]];
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
 //     	std::cout<<"sl1:"; for (int k=0;k<j-i+1;k++) std::cout<<" "<<sl1[k]<<"/"<<m_dummyarg[k+i]; std::cout<<std::endl;
	Complex v = Elementary_MHV_Amplitude(&m_dummyarg[i],sl1,j-i+1,0);
//   	std::cout<<"first "<<v<<std::endl;
	m_dummyarg[j] = perm[j];
	m_dummyarg[n_part+i] = pn;
//      	std::cout<<"sl2:"; for (int k=0;k<n_part-j+i+1;k++) std::cout<<" "<<sl2[k]<<"/"<<m_dummyarg[k+j]; std::cout<<std::endl;
	Complex vp = Elementary_MHV_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1,0);
	v*=vp;
//   	std::cout<<"second "<<vp<<std::endl;
	m_dummyarg[n_part+i] = perm[i];
	delete[] sl1;
	delete[] sl2;
	v /= (p_BS->Momentum(pn)).Abs2();
//    	std::cout<<"prop "<<(p_BS->Momentum(pn)).Abs2()<<std::endl;
//    	std::cout<<"tot: "<<v<<std::endl;;
	amp+=v;
      } 
    }    
  }
//   std::cout<<"amp: "<<amp<<std::endl;
  return amp;
}


Complex MHV_Calculator::MHVbar_Amplitude(int* perm,int* signlist,int vhel)
{
  if (vhel==2) return Elementary_MHVbar_Amplitude(perm,signlist,n_part);
  
  Complex amp(0.,0.);
//     std::cout<<"in MHV_Calculator::MHVbar_Amplitude"<<std::endl;

  for (int i=0;i<n_part;i++) m_dummyarg[i+n_part] = m_dummyarg[i] = perm[i];
  for (int i=0;i<n_part-2;i++) {
    for (int j=i+2;j<n_part && j-i+1<n_part;j++) {
      int cs=CountSign(&perm[i],signlist,j-i,1);
      if (cs==1 || cs==2) {
//        std::cout<<i<<"/"<<j<<":   "<<cs<<std::endl;
	int* sl1   = new int[j-i+1];
	int* sl2   = new int[n_part-j+i+1];
	for (int k=0;k<j-i;k++) sl1[k] = signlist[perm[i+k]];
	for (int k=0;k<n_part-(j-i);k++) sl2[k] = signlist[perm[(j+k)%n_part]];
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
//  	std::cout<<"sl1:"; for (int k=0;k<j-i+1;k++) std::cout<<" "<<sl1[k]<<"/"<<m_dummyarg[k+i]; std::cout<<std::endl;
	Complex v = Elementary_MHVbar_Amplitude(&m_dummyarg[i],sl1,j-i+1,0);
	m_dummyarg[j] = perm[j];
	m_dummyarg[n_part+i] = pn;
//  	std::cout<<"sl2:"; for (int k=0;k<n_part-j+i+1;k++) std::cout<<" "<<sl2[k]<<"/"<<m_dummyarg[k+j]; std::cout<<std::endl;
	v *= Elementary_MHVbar_Amplitude(&m_dummyarg[j],sl2,n_part-j+i+1,0);
	m_dummyarg[n_part+i] = perm[i];
	delete[] sl1;
	delete[] sl2;
	v /= (p_BS->Momentum(pn)).Abs2();
	amp+=v;
      } 
    }    
  }

  return amp;
}
