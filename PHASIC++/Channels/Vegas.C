#include "PHASIC++/Channels/Vegas.H"
#include <iostream>
#include "ATOOLS/Math/MathTools.H"
#include <vector>
#include <stdlib.h>
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/My_MPI.H"


using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

int Vegas::s_onext=-1, Vegas::s_on=-1;

Vegas::Vegas(int dim,int ndx,const std::string & name,int opt)
{
  if (s_on<0) {
    Data_Reader dr(" ",";","!","=");
    dr.AddComment("#");
    dr.AddWordSeparator("\t");
    dr.SetInputPath(rpa->GetPath());
    dr.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
    s_on = dr.GetValue<int>("VEGAS_MODE",2);
  }
  m_on=s_on?1:0;
  if (s_onext>-1) m_on=s_onext;
  m_dim  = dim;
  m_nopt = 0;
  m_nevt = 0;
  m_mnevt = 0;
  m_snevt = 0;
  m_cevt = 0;
  m_mcevt = 0;
  m_name = name;
  m_mode=0;
  m_nd=(s_on&2)?10:ndx;
  m_sint=(s_on&2)?1:0;
  m_scnt=0;
  m_alpha = 1.;
  m_autooptimize = -1;
  m_cmode = m_omode = 1;
  p_x   = new double[m_dim];
  p_bm=p_cx=NULL;
  p_cb=NULL;
  if (m_on!=0) {
    p_xi = new double*[m_dim];
    p_bestxi = new double*[m_dim];
    p_d  = new double*[m_dim];
    p_di  = new double*[m_dim];
    p_hit= new int*[m_dim]; 
#ifdef USING__MPI
    p_md  = new double*[m_dim];
    p_mdi  = new double*[m_dim];
    p_mhit= new int*[m_dim]; 
#endif
    for(int i=0;i<m_dim;i++) {
      p_xi[i] = new double[m_nd];
      p_bestxi[i] = new double[m_nd];
      p_d[i]  = new double[m_nd];
      p_di[i]  = new double[m_nd];
      p_hit[i]= new int[m_nd]; 
#ifdef USING__MPI
      p_md[i]  = new double[m_nd];
      p_mdi[i]  = new double[m_nd];
      p_mhit[i]= new int[m_nd]; 
#endif
    }
    p_dt  = new double[m_dim];
    p_chi = new double[m_dim];
    p_bestchi = new double[m_dim];
    p_xin = new double[m_nd];
    p_r   = new double[m_nd];
    p_ia  = new int[m_dim];
    p_opt = new int[m_dim];
    for (int i=0;i<m_dim;i++) {
      p_xi[i][0]=1.0;
      p_opt[i]=1;
      p_bestchi[i]=.0;
      for (int j=0;j<m_nd;j++) {
	p_d[i][j]=0.;
	p_di[i][j]=0.;
	p_hit[i][j]=0;
#ifdef USING__MPI
	p_md[i][j]=0.;
	p_mdi[i][j]=0.;
	p_mhit[i][j]=0;
#endif
      }
    }
    for (int i=0;i<m_nd;i++) p_r[i]=1.0;
    p_xin[m_nd-1] = 1.;
    for (int i=0;i<m_dim;i++) Rebin(1./double(m_nd),p_xi[i]);
    m_nc = pow(double(m_nd),double(m_dim));
    for (int j=0;j<m_dim;j++) for (int i=0;i<m_nd;i++) p_bestxi[j][i]=p_xi[j][i];
  }
}

Vegas::~Vegas() 
{
  delete[] p_x;
  if (p_cx) {
    delete[] p_bm;
    delete[] p_cx;
    delete[] p_cb;
  }
  if (m_on==0) return;
  for(int i=0;i<m_dim;i++) {
    delete[] p_xi[i];
    delete[] p_bestxi[i];
    delete[] p_d[i];
    delete[] p_di[i];
    delete[] p_hit[i];
#ifdef USING__MPI
    delete[] p_md[i];
    delete[] p_mdi[i];
    delete[] p_mhit[i];
#endif
  }
  delete[] p_xi;
  delete[] p_bestxi;
  delete[] p_d;
  delete[] p_di;
  delete[] p_hit;
#ifdef USING__MPI
  delete[] p_md;
  delete[] p_mdi;
  delete[] p_mhit;
#endif
  delete[] p_dt;
  delete[] p_xin;
  delete[] p_chi;
  delete[] p_bestchi;
  delete[] p_r;
  delete[] p_ia;
  delete[] p_opt;
}

void Vegas::InitBinInfo()
{
  p_bm  = new double[m_dim];
  p_cx  = new double[m_dim];
  p_cb  = new int[m_dim];
}

void Vegas::MPISync()
{
  if (!m_on) return;
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    int cn=3*m_dim*m_nd+2;
    double *values = new double[cn];
    double *rvalues = new double[cn];
    for (int i=0;i<m_dim;i++) {
      for (int j=0;j<m_nd;j++) {
	values[i*m_nd+j]=p_md[i][j];
	values[(m_dim+i)*m_nd+j]=p_mdi[i][j];
	values[(2*m_dim+i)*m_nd+j]=p_mhit[i][j];
      }
    }
    values[cn-2]=m_mnevt;
    values[cn-1]=m_mcevt;
    mpi->MPIComm()->Allreduce(values,rvalues,cn,MPI::DOUBLE,MPI::SUM);
    for (int i=0;i<m_dim;i++) {
      for (int j=0;j<m_nd;j++) {
	p_md[i][j]=rvalues[i*m_nd+j];
	p_mdi[i][j]=rvalues[(m_dim+i)*m_nd+j];
	p_mhit[i][j]=rvalues[(2*m_dim+i)*m_nd+j];
      }
    }
    m_mnevt=rvalues[cn-2];
    m_mcevt=rvalues[cn-1];
    delete [] values;
    delete [] rvalues;
  }
  for (int i=0;i<m_dim;i++) {
    for (int j=0;j<m_nd;j++) {
      p_d[i][j]+=p_md[i][j];
      p_di[i][j]+=p_mdi[i][j];
      p_hit[i][j]+=p_mhit[i][j];
      p_md[i][j]=p_mdi[i][j]=0.;
      p_mhit[i][j]=0;
    }
  }
  m_nevt+=m_mnevt;
  m_cevt+=m_mcevt;
  m_mnevt=m_mcevt=0;
#endif
}

void Vegas::Rebin(double rc, double * xi)
{
  int i,k=-1;
  double dr=0.0,xn=0.0,xo=0.;
  
  for (i=0;i<m_nd-1;i++) {
    while (rc > dr) {
      dr += p_r[++k];
      xo=xn;
      xn=xi[k];
    }
    dr -= rc;
    p_xin[i]=xn-(xn-xo)*dr/p_r[k];
  }
  for (i=0;i<m_nd;i++) xi[i]=p_xin[i];
}

double* Vegas::GeneratePoint(const double * ran) 
{
  if (m_on==0) {
    for (int i=0;i<m_dim;i++) p_x[i]=ran[i];
    return p_x;
  }
  m_mode=1;
  double xx;
  int ia;
  for (int i=0;i<m_dim;i++) {
    xx = ran[i]*(double)m_nd;
    ia = (int)xx;
    if (ia>=m_nd) {
      msg_Out()<<" WARNING Vegas::GeneratePoint(const double* ran)"
	       <<" called with ran["<<i<<"]="<<ran[i]<<"\n";
      ia=m_nd-1;
    }
    if (ia==0) {
      p_x[i] = xx*p_xi[i][0];
      if (p_cx) {
	p_cx[i] = p_xi[i][0]/2.0;
	p_cb[i] = 0;
      }
    }
    else {
      p_x[i] = p_xi[i][ia-1]+(xx-ia)*(p_xi[i][ia]-p_xi[i][ia-1]);
      if (p_cx) {
	p_cx[i] = (p_xi[i][ia]+p_xi[i][ia-1])/2.0;
	p_cb[i] = ia;
      }
    }
  }
  return p_x;
}

double Vegas::GenerateWeight(const double* xy) const
{
  if (m_on==0) return 1.;
  double weight = m_nc;
  for (int i=0;i<m_dim;i++) {
    size_t l(0), r(m_nd-1), c((l+r)/2);
    double a(p_xi[i][c]);
    while (r-l>1) {
      if (xy[i]<a) r=c;
      else l=c;
      c=(l+r)/2;
      a=p_xi[i][c];
    }
    if (xy[i]<p_xi[i][l]) r=l;
    int k(r);
    p_ia[i] = k;
    if (k==0) {
      weight *= p_xi[i][k];
    }
    else {
      weight *= p_xi[i][k]-p_xi[i][k-1];
    }
  }
  return weight;
}

double Vegas::GenerateBinWeight(int* xy) const
{
  if (m_on==0) return 1.;
  double weight = m_nc;
  for (int i=0;i<m_dim;i++) {
    int k(xy[i]);
    p_ia[i] = k;
    if (k==0) {
      weight *= p_xi[i][k];
    }
    else {
      weight *= p_xi[i][k]-p_xi[i][k-1];
    }
  }
  return weight;
}

double *Vegas::GetBinsMean(int *xy) const
{
  for (int i=0;i<m_dim;i++) {
    int k(xy[i]);
    if (k==0) p_bm[i]=p_xi[i][k]/2.0;
    else p_bm[i]=(p_xi[i][k-1]+p_xi[i][k])/2.0;
  }
  return p_bm;
}

void Vegas::AddPoint(double value,double *xy)
{ if (m_on==0) return;
  if (m_mode==1 && m_cmode&1){
    for (int i=0;i<m_dim;i++) {
      if (dabs(p_x[i]-xy[i])>1.e-4) {
	msg_Tracking()<<"Mapping error in Vegas for "<<m_name<<endl;
	for (int j=0;j<m_dim;j++) msg_Tracking()<<j<<": "<<p_x[j]<<"<->"<<xy[j]<<" ("<<dabs(p_x[j]-xy[j])<<")"<<endl;
//    	abort();
	i=m_dim;
      }
    }
  }
  for (int i=0;i<m_dim;i++) {
    size_t l(0), r(m_nd-1), c((l+r)/2);
    double a(p_xi[i][c]);
    while (r-l>1) {
      if (xy[i]<a) r=c;
      else l=c;
      c=(l+r)/2;
      a=p_xi[i][c];
    }
    if (xy[i]<p_xi[i][l]) r=l;
    int k(r);
    p_ia[i] = k;
  }
  AddPoint(value);
}

void Vegas::AddPoint(double value)
{
#ifdef USING__MPI
  ++m_mnevt;
  if (value>0.) ++m_mcevt;
  double v2 = value*value;
  for (int i=0;i<m_dim;i++) {
    p_md[i][p_ia[i]]+=v2;
    p_mdi[i][p_ia[i]]+=v2*v2;
    p_mhit[i][p_ia[i]]++;
  }
  m_mode=0;
  if (MPI::COMM_WORLD.Get_size()>1) {
    if (m_autooptimize>0)
      THROW(fatal_error,"Autooptimize not possible in MPI mode");
  }
  else {
    if (m_autooptimize>0&&m_nevt%m_autooptimize==0) {
      int v=(m_nevt-m_snevt)/m_autooptimize;
      if (m_cevt*10*v>(unsigned long)m_autooptimize) { 
	if (m_nopt==0) { 
	  if(m_cevt*2>(unsigned long)m_nd) Optimize();
	}
	else if(m_cevt>m_nd*m_nopt) Optimize();
      }
    }
  }
#else
  ++m_nevt;
  if (value>0.) ++m_cevt;
  double v2 = value*value;
  for (int i=0;i<m_dim;i++) {
    p_d[i][p_ia[i]]+=v2;
    p_di[i][p_ia[i]]+=v2*v2;
    p_hit[i][p_ia[i]]++;
  }
  m_mode=0;
  if (m_autooptimize>0&&m_nevt%m_autooptimize==0) {
    int v=(m_nevt-m_snevt)/m_autooptimize;
    if (m_cevt*10*v>(unsigned long)m_autooptimize) { 
      if (m_nopt==0) { 
	if(m_cevt*2>(unsigned long)m_nd) Optimize();
      }
      else if(m_cevt>m_nd*m_nopt) Optimize();
    }
  }
#endif
}

void Vegas::AddBinPoint(double value,int *xy)
{ 
  if (m_on==0) return;
  for (int i=0;i<m_dim;i++) p_ia[i] = xy[i];
  AddPoint(value);
}

void Vegas::Reset()
{
  for (int i=0;i<m_dim;i++) {
    for (int j=0;j<m_nd;j++) {
      p_d[i][j]=p_di[i][j]=0.;
      p_hit[i][j]=0;
#ifdef USING__MPI
      p_md[i][j]=p_mdi[i][j]=0.;
      p_mhit[i][j]=0;
#endif
    }
  }
  m_snevt=m_nevt;
  m_cevt=0;
#ifdef USING__MPI
  m_mnevt=m_mcevt=0;
#endif
}

void Vegas::Optimize()
{ 
  if (m_on==0) return;
  if (m_nevt-m_snevt<(unsigned int)m_nd*20) return;
  if (m_omode&1)
    msg_Tracking()<<"Vegas optimize "<<m_name<<" "<<m_nopt<<" |"
		  <<m_nevt-m_snevt<<" ("
		  <<(double)m_cevt/(m_nevt-m_snevt)*100.<<" % )"<<endl;
  for (int j=0;j<m_dim;j++) if(p_opt[j]) {
    for (int i=0;i<m_nd;i++) {if(p_hit[j][i]) p_d[j][i]/=p_hit[j][i];
     if (p_hit[j][i]<2) p_hit[j][i]=2;
    }
    double av=0.,av2=0.;
    for (int i=0;i<m_nd;i++)av+=p_d[j][i];
    av/=m_nd;
    for (int i=0;i<m_nd;i++)av2+=sqr(p_d[j][i]);
    av2/=m_nd;
    double ts = sqrt((av2-sqr(av))/(m_nd-1));
    double chi=0.;
    double s2;
    double cx;
    for (int i=0;i<m_nd;i++){
      s2=(p_di[j][i]/p_hit[j][i] - sqr(p_d[j][i]))/(p_hit[j][i]-1);
      cx=sqr(p_d[j][i]-av)/s2;
      if ((p_d[j][i]<av && cx>1.e4) || !(cx>=0.)) cx=1.e4;
      chi+=cx;
    }
    chi=sqrt(chi/m_nd);
    if (m_nopt==0||chi<p_chi[j]) {
      p_chi[j]=chi;
      for (int i=0;i<m_nd;i++) p_bestxi[j][i]=p_xi[j][i];
    }
    if (chi<=1.1) {
      p_opt[j]=0;
      for (int i=0;i<m_nd;i++) p_bestxi[j][i]=p_xi[j][i];
    }
    if ((chi>2.*p_chi[j] && ts>p_bestchi[j] && m_nopt>=10) 
	|| (m_nopt>20&&chi>p_chi[j])) {
      p_opt[j]=0;
      for (int i=0;i<m_nd;i++) p_xi[j][i]=p_bestxi[j][i];
    }
    if (ts<p_bestchi[j]) p_bestchi[j]=ts;
    if (m_omode&1)
      msg_Tracking()<<"Chi"<<j<<" ="<<chi<<"   "<<p_chi[j]
		    <<"("<<p_opt[j]<<")"<<endl;
  }
  double xo,xn,rc;
  for (int j=0;j<m_dim;j++) if(p_opt[j]) {
    xo=p_d[j][0];
    xn=p_d[j][1];
     p_d[j][0]=(xo+xn)/2.0;
     p_dt[j]=p_d[j][0];
     for (int i=1;i<m_nd-1;i++) {
       rc=xo+xn;
       xo=xn;
       xn=p_d[j][i+1];
       p_d[j][i] = (rc+xn)/3.0;
       p_dt[j] += p_d[j][i];
      }
    p_d[j][m_nd-1]=(xo+xn)/2.0;
    p_dt[j] += p_d[j][m_nd-1];
  }
  for (int j=0;j<m_dim;j++) if(p_opt[j]) {
    rc=0.0;
    for (int i=0;i<m_nd;i++) {
     if (p_d[j][i] < p_dt[j]*1.e-10) p_d[j][i]=p_dt[j]*1.e-10;
      p_r[i]=pow((1.0-p_d[j][i]/p_dt[j])/
	       (log(p_dt[j])-log(p_d[j][i])),m_alpha);
      rc += p_r[i];
    }
    Rebin(rc/double(m_nd),p_xi[j]);
  }
  m_nopt++;
  if (m_sint>0 && ++m_scnt==m_sint) Refine();
  Reset();
}

void Vegas::Refine()
{
  if (m_omode&1)
    msg_Tracking()<<"Refine '"<<m_name<<"' "<<m_nd<<" -> "
		  <<2*m_nd<<" ( int = "<<m_sint<<" )\n";
  m_scnt=0;
  ++m_sint;
  m_nd*=2;
  m_nc=pow(double(m_nd),double(m_dim));
  delete [] p_xin; p_xin = new double[m_nd];
  delete [] p_r; p_r = new double[m_nd];
  for(int i=0;i<m_dim;++i) {
    p_r[i]=1.0;
    p_xin[m_nd-1] = 1.;
    std::vector<double> xi(&p_xi[i][0],&p_xi[i][m_nd/2]);
    delete [] p_xi[i]; p_xi[i] = new double[m_nd];
    delete [] p_bestxi[i]; p_bestxi[i] = new double[m_nd];
    delete [] p_d[i]; p_d[i] = new double[m_nd];
    delete [] p_di[i]; p_di[i] = new double[m_nd];
    delete [] p_hit[i]; p_hit[i] = new int[m_nd]; 
#ifdef USING__MPI
    delete [] p_md[i]; p_md[i] = new double[m_nd];
    delete [] p_mdi[i]; p_mdi[i] = new double[m_nd];
    delete [] p_mhit[i]; p_mhit[i] = new int[m_nd]; 
#endif
    double lx=0.0;
    for (int j=0;j<m_nd;++j) {
      p_xi[i][j]=(j%2==0)?(lx+xi[j/2])/2.0:lx=xi[j/2];
      p_bestxi[i][j]=p_xi[i][j];
    }
  }
}

void Vegas::EndOptimize()
{ 
  if (m_on==0||m_nopt==0) return;
  msg_Tracking()<<"Vegas EndOptimize: "<<m_name<<endl;
  for (int j=0;j<m_dim;j++) msg_Tracking()<<" "<<p_chi[j];
  msg_Tracking()<<endl;
  for (int j=0;j<m_dim;j++)
    for (int i=0;i<m_nd;i++) p_xi[j][i]=p_bestxi[j][i]; 
  m_autooptimize=-1;
}

std::vector<double> Vegas::GetMaxPos() const
{
  std::vector<double> maxs(m_dim,0.0);
  for (int i=0;i<m_dim;i++) {
    double min(1.0);
    for (int j=0;j<m_nd-1;j++) {
      double d(p_xi[i][j+1]-p_xi[i][j]);
      if (d<min) {
	min=d;
	maxs[i]=p_xi[i][j]+d/2.0;
      }
    }
  }
  return maxs;
}

std::vector<double> Vegas::GetMeanPos() const
{
  std::vector<double> means(m_dim,0.0);
  for (int i=0;i<m_dim;i++) {
    for (int j=0;j<m_nd-1;j++)
      means[i]+=(sqr(p_xi[i][j+1])-sqr(p_xi[i][j]))/2.0;
  }
  return means;
}

void Vegas::WriteHistos(const std::string & pid)
{
  bool mh=false;
  double ar=1./m_nd;
  double x=0.;
  if (mh)  {
    for (int i=0;i<m_dim;i++) {
      std::string fn=pid+std::string("_")+m_name+std::string("_Vegas_")
	+ToString(i)+std::string(".dat");
      My_Out_File ofile(fn);
      ofile.Open();
      *ofile<<x<<" "<<ar/p_xi[i][0]<<endl;
      for (int j=0;j<m_nd-1;j++) {
	*ofile<<x+p_xi[i][j]<<" "<<ar/(p_xi[i][j+1]-p_xi[i][j])<<endl;
      }
      *ofile<<x+1.<<" 0."<<endl;
      ofile.Close();
    }
  }
  else {
    std::string fn=pid+std::string("_")+m_name+std::string("_Vegas_")+std::string(".dat");
    My_Out_File ofile(fn);
    ofile.Open();
    for (int i=0;i<m_dim;i++) {
      *ofile<<x<<" "<<ar/p_xi[i][0]<<endl;
      for (int j=0;j<m_nd-1;j++) {
	*ofile<<x+p_xi[i][j]<<" "<<ar/(p_xi[i][j+1]-p_xi[i][j])<<endl;
      }
      *ofile<<x+1.<<" 0."<<endl;
      x+=1.;
    }
    ofile.Close();
  }  
}

void Vegas::WriteOut(const std::string & pid)
{
  if (msg_LevelIsTracking() && m_on) 
    WriteHistos(pid);
  std::string fn=pid+std::string("_")+m_name+std::string("_Vegas");
  My_Out_File ofile(fn);
  ofile.Open();

  *ofile<<m_name<<" "<<m_dim<<" "<<m_nd<<" "
       <<m_autooptimize<<" "<<m_nopt<<" "
       <<m_sint<<" "<<m_scnt<<std::endl;
  if (m_nopt>0) {
    ofile->precision(12);
    for (int i=0;i<m_dim;i++) {
      *ofile<<"(";
      for (int j=0;j<m_nd;j++) {
	if (j!=0) *ofile<<",";
	*ofile<<p_xi[i][j];
      }
      *ofile<<")"<<endl;
    }
    for (int i=0;i<m_dim;i++) {
      *ofile<<p_opt[i]<<" "<<p_chi[i]<<" (";
      for (int j=0;j<m_nd;j++) {
	if (j!=0) *ofile<<",";
	*ofile<<p_bestxi[i][j];
      }
      *ofile<<")"<<endl;
    }
  }
  ofile.Close();
}

void Vegas::ReadIn(const std::string & pid)
{
  std::string fn=pid+std::string("_")+m_name+std::string("_Vegas"), tn;
  My_In_File ifile(fn);
  if (!ifile.Open()) return;
  *ifile>>tn;
  if (tn!=m_name) THROW(fatal_error,"Corrupted input file");
  int nd;
  *ifile>>m_dim>>nd>>m_autooptimize>>m_nopt>>m_sint>>m_scnt;
  if (nd!=m_nd) {
    m_nd=nd;
    m_nc=pow(double(m_nd),double(m_dim));
    delete [] p_xin; p_xin = new double[m_nd];
    delete [] p_r; p_r = new double[m_nd];
    for(int i=0;i<m_dim;++i) {
      p_r[i]=1.0;
      p_xin[m_nd-1] = 1.;
      delete [] p_xi[i]; p_xi[i] = new double[m_nd];
      delete [] p_bestxi[i]; p_bestxi[i] = new double[m_nd];
      delete [] p_d[i]; p_d[i] = new double[m_nd];
      delete [] p_di[i]; p_di[i] = new double[m_nd];
      delete [] p_hit[i]; p_hit[i] = new int[m_nd]; 
#ifdef USING__MPI
      delete [] p_md[i]; p_md[i] = new double[m_nd];
      delete [] p_mdi[i]; p_mdi[i] = new double[m_nd];
      delete [] p_mhit[i]; p_mhit[i] = new int[m_nd]; 
#endif
    }
  }
  if (m_nopt==0||m_on==0) return;
  std::string buffer;
  getline(*ifile,buffer);
  for (int i=0;i<m_dim;++i) {
    getline(*ifile,buffer);
    size_t  a=buffer.find("(")+1;
    size_t  b=buffer.find(")");
    char * err;
    buffer=buffer.substr(a,b-a);
    for (int j=0;j<m_nd;++j) {
      size_t c=buffer.find(",");
      p_xi[i][j]=strtod(buffer.substr(0,c).c_str(),&err);
      buffer=buffer.substr(c+1);
    }
  }
  for (int i=0;i<m_dim;++i) {
    *ifile>>p_opt[i]>>p_chi[i];
    getline(*ifile,buffer);
    size_t  a=buffer.find("(")+1;
    size_t  b=buffer.find(")");
    char * err;
    buffer=buffer.substr(a,b-a);    
    for (int j=0;j<m_nd;++j) {
      size_t c=buffer.find(",");
      p_bestxi[i][j]=strtod(buffer.substr(0,c).c_str(),&err);
      buffer=buffer.substr(c+1);
    }
  }
}

bool Vegas::Finished()
{
  if (m_on==0) return true;
  for (int j=0;j<m_dim;j++) if(p_opt[j]) return false;  
  return true;
}

