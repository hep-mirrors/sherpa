#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <stdio.h>

using namespace ATOOLS;

template <class Type>
Type Get(const std::string & in) 
{
  if (in!="nan") {
    Type value;
    MyStrStream str;
    str<<in;
    str>>value;
    return value;
  }
  return (Type)0;
}

Histogram::Histogram(int _type,double _lower,double _upper,int _nbin) :
  m_type(_type), m_nbin(_nbin), m_lower(_lower), m_upper(_upper), 
  m_yvalues(0),m_y2values(0), m_psvalues(0), m_tmp(0), m_fills(0), m_psfills(0), 
  m_finished(false), m_initialized(false), m_fuzzyexp(-1)
{
  m_mcb = 0.;
  if (m_type>1000) {
    m_type-=1000;
    m_fuzzyexp = int(m_type/100);
    m_type-=100*m_fuzzyexp;
  }
  if (m_type>=100) {
    m_mcb = 1.;
    m_type-=100;
  }
  m_logarithmic = int(m_type/10);
  m_depth       = m_type-m_logarithmic*10+1;
  m_logbase = 1;
  switch(m_logarithmic) {
    case 1: 
      m_logbase = log(10.);
      m_upper   = log(m_upper)/m_logbase; m_lower = log(m_lower)/m_logbase;
      break;
    case 2: 
      m_upper   = log(m_upper); m_lower = log(m_lower);
      break;
    default: break;
  }
  m_binsize     = (m_upper-m_lower)/(double(m_nbin));

  if (m_binsize<=0.) {
    msg_Error()<<"Error in Histogram : Tried to initialize a histogram with  binsize <= 0 !"<<std::endl;
    m_active = 0;
    return;
  }

  m_nbin += 2;
  m_active = 1;
  m_yvalues   = new double[m_nbin];

  if (m_depth>1) {
    m_y2values   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_y2values[i]=0.;
    }
  }

  if (m_depth>2) {
    m_psvalues   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_psvalues[i]=0.;
    }
  }

  if (m_mcb!=0.) {
    m_tmp   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_tmp[i]=0.;
    }
  }

  for (int i=0;i<m_nbin;i++) {
    m_yvalues[i]=0.;
  }
}

Histogram::Histogram(const Histogram * histo)
: m_yvalues(0), m_y2values(0), m_psvalues(0), m_tmp(0) {
  m_lower   = histo->m_lower;
  m_upper   = histo->m_upper;
  m_logbase = histo->m_logbase;
  m_logarithmic = histo->m_logarithmic;
  m_mcb     = histo->m_mcb;
  m_nbin    = histo->m_nbin;
  m_depth   = histo->m_depth;
  m_type    = histo->m_type;
  m_fills   = histo->m_fills;
  m_psfills = histo->m_psfills;

  m_binsize = histo->m_binsize;
  m_active  = 1;
  m_finished = histo->m_finished;

  m_yvalues   = new double[m_nbin];
  for (int i=0;i<m_nbin;i++) {
    m_yvalues[i]  = histo->m_yvalues[i]; 
  }
  if (m_depth>1) {
    m_y2values   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_y2values[i]=histo->m_y2values[i];
    }
  }
  if (m_depth>2) {
    m_psvalues   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_psvalues[i]=histo->m_psvalues[i];
    }
  }
  if (m_mcb!=0.) {
    m_tmp   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_tmp[i]=0.;
    }
  }
}


Histogram::Histogram(const std::string & pID)
  :  m_yvalues(0), m_y2values(0), m_psvalues(0), m_tmp(0), m_fills(0), m_mcb(0.)  {
  m_finished=true;
  std::ifstream ifile(pID.c_str());

  std::string dummy;
  getline(ifile,dummy);
  
  Data_Reader dr(" ",";","!");
  dr.AddWordSeparator("\t");
  if (dummy!="") { 
    std::vector<std::string> conf;
    dr.SetString(dummy);
    dr.VectorFromString(conf);
    size_t k=0;

    if (k>=conf.size()) {
      msg_Error()<<"Error in Histogram : reading file :"<<pID<<std::endl;
      m_active = 0;
      return;
    }

    MyStrStream str;
    str<<dummy;
    str>>m_type>>m_nbin>>m_lower>>m_upper;
    k=4;

    m_logarithmic = int(m_type/10);
    m_depth       = m_type-m_logarithmic*10+1;

    m_logbase = 1;
    switch(m_logarithmic) {
    case 1: 
      m_logbase = log(10.);
      break;
    default: break;
    }
    m_binsize     = (m_upper-m_lower)/(double(m_nbin-2));
    
    if (m_binsize<=0.) {
      msg_Error()<<"Error in Histogram : "
		 <<"Tried to initialize a histogram with m_binsize <= 0 !"
		 <<std::endl;
      m_active = 0;
      return;
    }
    
    m_active = 1;
    m_yvalues   = new double[m_nbin];
    m_yvalues[0]  = Get<double>(conf[k++]);
    if (m_depth>1) {
      m_y2values   = new double[m_nbin];
      m_y2values[0] = sqr(Get<double>(conf[k++]));
    }    
    if (m_depth>2) {
      m_psvalues   = new double[m_nbin];
      m_psvalues[0] = Get<double>(conf[k++]);
    }    
    if (k>=conf.size()) {
      msg_Error()<<"Error in Histogram : reading file :"<<pID<<std::endl;
      m_active = 0;
      return;
    }

    m_yvalues[m_nbin-1]  = Get<double>(conf[k++]);
    if (m_depth>1) {
      m_y2values[m_nbin-1] = sqr(Get<double>(conf[k++]));
    }    
    if (m_depth>2) {
      m_psvalues[m_nbin-1] = Get<double>(conf[k++]);
    }    
    if (k>=conf.size()) {
      msg_Error()<<"Error in Histogram : reading file :"<<pID<<std::endl;
      m_active = 0;
      return;
    }
    m_fills = Get<double>(conf[k++]);
  }
  else {
    msg_Error()<<"Error in Histogram : reading file :"<<pID<<std::endl;
    m_active = 0;
    return;    
  }
  

  std::vector<std::string> data;
  MyStrStream str;
  for (int i=0;i<m_nbin-1;i++) {
    getline(ifile,dummy);
    data.clear();
    dr.SetString(dummy);
    dr.VectorFromString(data);

    //    ifile>>value;
    m_yvalues[i+1] = Get<double>(data[1]);
    if (m_depth>1) {
      m_y2values[i+1] = Get<double>(data[2]);
      m_y2values[i+1] = sqr(m_y2values[i+1]);
    }    
    if (m_depth>2) {
      m_psvalues[i+1] = Get<double>(data[3]);
    }    
  }
  ifile.close();
}


Histogram::~Histogram() {
    
  if (m_yvalues!=0) { 
    delete [] m_yvalues; m_yvalues = 0; 
  }
  if (m_y2values!=0) { 
    delete [] m_y2values; m_y2values = 0; 
  }
  if (m_psvalues!=0) { 
    delete [] m_psvalues; m_psvalues = 0; 
  }
  if (m_tmp!=0) { 
    delete [] m_tmp; m_tmp = 0; 
  }
}


void Histogram::Finalize() {
  if (!m_finished) {
    m_finished=true;
    if (m_fills==0.) return;
    for (int i=0;i<m_nbin;++i) {
      m_yvalues[i]/=m_fills*m_binsize;
      if (m_depth>1) {
 	m_y2values[i]/=m_fills*sqr(m_binsize);
	if (m_fills>1) m_y2values[i]=(m_y2values[i]-sqr(m_yvalues[i]))/(m_fills-1);
      }
    }
    if (m_depth>2) {
      double itg = Integral()/(m_psfills*m_binsize);
      for (int i=0;i<m_nbin;++i) {
	m_psvalues[i]*=itg;
      }
    }
  }
}

void Histogram::Restore() {
  if (m_finished) {
    for (int i=0;i<m_nbin;++i) {
      if (m_depth>1) {
	if (m_fills>1) m_y2values[i]=(m_fills-1)*m_y2values[i]+sqr(m_yvalues[i]);
	m_y2values[i]*=m_fills*sqr(m_binsize);
	if (m_depth>2) {
	  m_psvalues[i]*=m_psfills*m_binsize;
	}
      }
      m_yvalues[i]*=m_fills*m_binsize;
    }
    m_finished=false;
  }
}

double Histogram::Mean() const
{
  double sum=0., range=0.;
  for (int i=1;i<m_nbin-1;++i) {
    double width=(m_upper-m_lower)/m_nbin;
    if (m_logarithmic) 
      width=pow(m_logbase,m_lower+i*width)-pow(m_logbase,m_lower+(i-1)*width);  
    sum+=m_yvalues[i]*width;
    range+=width;
  }
  return sum/range;
}

void Histogram::Reset() {
  for (int i=0;i<m_nbin;i++) { 
    m_yvalues[i]=0.;
    if (m_depth>1) {
      m_y2values[i]=0.;
    }
    if (m_depth>2) {
      m_psvalues[i]=0.;
    }
  }
  m_fills=0;
  m_psfills=0;
}

void Histogram::Scale(double scale) {
  for (int i=0;i<m_nbin;i++) { 
    m_yvalues[i]*= scale;
    if (m_depth>1) m_y2values[i]*=sqr(scale); 
    if (m_depth>2) m_psvalues[i]*=scale; 
  }
}

void Histogram::Output() {
  if (!msg_LevelIsDebugging()) return;
  msg_Out()<<"----------------------------------------"<<std::endl
	   <<"    "<<m_yvalues[0]<<std::endl
	   <<"----------------------------------------"<<std::endl;
  double result = 0.;
  for (int i=0;i<m_nbin-2;i++) {
    msg_Out()<<m_lower+i*m_binsize<<"  ";
    msg_Out()<<m_yvalues[i+1]<<"  ";
    if (m_depth>1) msg_Out()<<sqrt(m_y2values[i+1]);
    result += m_yvalues[i+1];
    msg_Out()<<std::endl;
  }
  msg_Out()<<m_lower+(m_nbin-2)*m_binsize<<" == "<< m_upper<<std::endl
	   <<"----------------------------------------"<<std::endl
	   <<"    "<<m_yvalues[m_nbin-1]<<std::endl
	   <<"----------------------------------------"<<std::endl
	   <<"Inside the range : "<<result<<std::endl;
}


void Histogram::Output(const std::string name) 
{
  if (!m_active) return;
  std::ofstream ofile;
  ofile.open(name.c_str());

  if (m_fills>=0) {
    ofile<<m_type<<" "<<m_nbin<<" "<<m_lower<<" "<<m_upper<<" ";
    ofile<<m_yvalues[0]<<"  ";
    if (m_depth>1) ofile<<m_y2values[0]<<"  ";
    ofile<<m_yvalues[m_nbin-1]<<"  ";
    if (m_depth>1) ofile<<m_y2values[m_nbin-1]<<"  ";
    ofile<<m_fills<<"\n";
  }
  for (int i=0;i<m_nbin-1;i++) {
    ofile<<m_lower+i*m_binsize<<"  ";
    ofile<<m_yvalues[i+1]<<"  ";
    if (m_depth>1) ofile<<sqrt(m_y2values[i+1])<<"  ";
    if (m_depth>2) ofile<<m_psvalues[i+1]<<"  ";
    ofile<<"\n";
  }
  ofile.close();
}







void Histogram::Insert(double coordinate) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }

  m_fills++;

  if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;
  if (coordinate<m_lower) { m_yvalues[0 ]       += double(1); return; }
  if (coordinate>m_upper) { m_yvalues[m_nbin-1] += double(1); return; }
  for (int i=1;i<m_nbin-1;i++) {
    if ( (coordinate >= m_lower + (i-1)*m_binsize) &&
	 (coordinate <  m_lower + i*m_binsize) ) {
      m_yvalues[i] += double(1); 
      return; 
    }
  }
}

void Histogram::Insert(int i,double value,double ncount) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  m_fills+=ncount;
  if (value==0.) return;
  m_psfills++;

  if (i<0) { 
    m_yvalues[0] += value;
    if (m_depth>1) {
      if (value>m_y2values[0]) m_y2values[0] = value;
      if (m_depth>2) m_psvalues[0] += 1.;
    }
    return; 
  }

  if (i>=m_nbin) { 
    m_yvalues[m_nbin-1] += value; 
    if (m_depth>1) {
      if (value>m_y2values[m_nbin-1]) m_y2values[m_nbin-1] = value;
      if (m_depth>2) m_psvalues[m_nbin-1] += 1.;
    }
    return; 
  }

  m_yvalues[i] += value;
  if (m_depth>1) {
    m_y2values[i] += value*value;
    if (m_depth>2) m_psvalues[i] += 1.;
  }
}

void Histogram::InsertMCB(double coordinate,double value,double ncount) {
  if (!m_tmp) {
    m_tmp   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_tmp[i]=0.;
    }
  }
  m_mcb = ncount;

  if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

  int bin = int((coordinate-m_lower)/m_binsize+1.);
  if (bin<0) bin=0;
  if (bin>=m_nbin) bin=m_nbin-1;
  if (bin==0||bin==m_nbin-1) {
    m_tmp[bin] += value;
    return;
  }

  double x = (coordinate-m_lower)/m_binsize-double(bin)+0.5;
  if ((bin==1&&x<0.)||(bin==m_nbin-2&&x>0.)) {
    m_tmp[bin] += value;
    return;
  }
  double ff=1.;
  if (m_fuzzyexp==0) ff=0.5;
  if (m_fuzzyexp>0) ff=1.-0.5*pow(2.*dabs(x),m_fuzzyexp);
  if (m_fuzzyexp==9) ff=1.-0.5*sqrt(2.*dabs(x));

  m_tmp[bin] += ff*value; 
  if (x>0.) m_tmp[bin+1] += (1.-ff)*value;
  if (x<0.) m_tmp[bin-1] += (1.-ff)*value;
}

void Histogram::InsertMCBIM(double coordinate,double value) {
  if (!m_tmp) {
    m_tmp   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_tmp[i]=0.;
    }
  }
  m_mcb = 1.;
  if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

  int bin = int((coordinate-m_lower)/m_binsize+1.);
  if (bin<0) bin=0;
  if (bin>=m_nbin) bin=m_nbin-1;
  for (int i=bin+1;i<m_nbin;i++) m_tmp[i] += value;
}

void Histogram::FinishMCB()
{
  m_fills+=m_mcb;

  for (int i=0;i<m_nbin;i++) {
    m_yvalues[i] += m_tmp[i];
    if (m_depth>1) {
      m_y2values[i] += m_tmp[i]*m_tmp[i];
      if (m_depth>2) m_psvalues[i] += 1.;
    }
    m_tmp[i] = 0.;
  }
}

void Histogram::Insert(double coordinate,double value,double ncount) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  m_fills+=ncount;
  if (value==0.) return;
  m_psfills++;

  if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

  int bin = int((coordinate-m_lower)/m_binsize+1.);
  if (bin<0) bin=0;
  if (bin>=m_nbin) bin=m_nbin-1;
  if (bin==0||bin==m_nbin-1) {
    m_yvalues[bin] += value;
    if (m_depth>1) {
      if (value>m_y2values[bin]) m_y2values[bin] = value;
      if (m_depth>2) m_psvalues[bin] += 1.;
    }
    return;
  }

  m_yvalues[bin] += value; 
  if (m_depth>1) {
    m_y2values[bin] += sqr(value);
    if (m_depth>2) m_psvalues[bin] += 1.;
  }
 
  if (m_fuzzyexp<0) return;

  double x = (coordinate-m_lower)/m_binsize-double(bin)+0.5;
  if (bin==1&&x<0.) return;
  if (bin==m_nbin-2&&x>0.) return;
  double ff=1.;
  if (m_fuzzyexp==0) ff=0.5;
  if (m_fuzzyexp>0) ff=0.5*pow(2.*dabs(x),m_fuzzyexp);
  if (m_fuzzyexp==9) ff=0.5*sqrt(2.*dabs(x));

  m_yvalues[bin] -= ff*value; 
  if (m_depth>1) {
    m_y2values[bin] += sqr(ff*(value))-sqr(value);
    if (m_depth>2) m_psvalues[bin] -= ff;
  }

  if (x>0.) {
    m_yvalues[bin+1] += ff*value;
    if (m_depth>1) {
      m_y2values[bin+1] += sqr(ff*value);
      if (m_depth>2) m_psvalues[bin+1] += ff;
    }
  }
  if (x<0.) {
    m_yvalues[bin-1] += ff*value;
    if (m_depth>1) {
      m_y2values[bin-1] += sqr(ff*value);
      if (m_depth>2) m_psvalues[bin-1] += ff;
    }
  }
}

void Histogram::InsertRange(double start, double end, double value) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  if (m_logarithmic>0) {
    if (start>0)
      start = log(start)/m_logbase;
    else
      start = -30;
    if (end>0)
      end = log(end)/m_logbase;
    else 
      end = -30;
  }
  m_fills++;

  // underrun
  if (start<m_lower) { 
    m_yvalues[0] += value;
    if (m_depth>1) {
      if (value>m_y2values[0]) m_y2values[0] = value;
    }
  }

  // overflow
  if (start>m_upper) { 
    m_yvalues[m_nbin-1] += value; 
    if (m_depth>1) {
      if (value>m_y2values[m_nbin-1]) m_y2values[m_nbin-1] = value;
    }
  }

  double low,up;
  int hit=0;
  low = m_lower; up = m_lower+m_binsize;
  for (int i=1;i<m_nbin-1;i++) {
    if ((start < up) && (end >= low) ) {
      double fac=1;
      if ((start<=low)&&(up<=end )) {
	m_yvalues[i] += value;
      } 
      else if ((low<start)&&(up<=end)) {
	fac = (start-low)/m_binsize;
	m_yvalues[i] += value *fac;
      }
      else if ((start<=low)&&(end < up)) {
	fac = (up-end)/m_binsize;
	m_yvalues[i] += value *fac;
      }
      else if ((low<start)&&(end <up)) {
	fac = (end-start)/m_binsize;
	m_yvalues[i] += value *fac;
      }
    

      hit=1;
      if (m_depth>1) {
	if (value>m_y2values[i]) m_y2values[i] = value;
      }
    }
    low = up;
    up += m_binsize;
  }
}



double Histogram::Bin(int bin) { return m_yvalues[bin]; }

double Histogram::Bin(double coordinate) 
{ 
  if (!m_active) {
    msg_Error()<<"Error in Histogram : Tried to access a histogram wit binsize <= 0 ! Return 0.."<<std::endl;
    return -1.0;
  }
  else {
    if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

    if (coordinate<m_lower) return m_yvalues[0];
    if (coordinate>m_upper) return m_yvalues[m_nbin-1];
    for (int i=1;i<m_nbin+1;i++) {
      if ( (coordinate >= m_lower + (i-1)*m_binsize) &&
	   (coordinate <  m_lower + i*m_binsize) ) 
	return m_yvalues[i];
    }
  }
  return -1.0;
}

void Histogram::Extrapolate(double coordinate,double * res,int mode) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram : Tried to access a histogram with binsize <= 0 ! Return 0.."<<std::endl;
  }
  else {
    if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

    for (int i=1;i<m_nbin;i++) {
      if ( (coordinate >= m_lower + (i-1)*m_binsize) &&
	   (coordinate <  m_lower + i*m_binsize) ) {
	res[0] = m_yvalues[i-1] +
	  (m_yvalues[i]-m_yvalues[i-1])/m_binsize *
	  (coordinate - m_lower - (i-1) * m_binsize);
	double over = 0.; double under = 0.;
	switch (mode) {
	case 1  : 
	  under = (coordinate - m_lower - (i-1) * m_binsize)/m_binsize * m_yvalues[i];
	  over = (m_lower + (i) * m_binsize - coordinate)/m_binsize * m_yvalues[i];
	  for (int j=0;j<i;j++) {
	    under += m_yvalues[j];
	  }
	  for (int j=i;j<m_nbin-1;j++) {
	    over += m_yvalues[j+1];
	  }
	  res[0]  = over;

	  if (m_depth>1) {
	    res[1]=0.;
	    for (int j=i;j<m_nbin;j++) {
	      res[1] = Max(res[1],m_y2values[j]);
	    }
	  }


	  break;
	case -1 :
	  for (int j=i-1;j>0;j--) {
	    over  += m_yvalues[j];
	    under += m_yvalues[j-1];
	  }
	  over    += m_yvalues[0];
	  res[0]  += (over+under)/2;
	  break;
	default : break;
	}
	/*
	if (m_depth>0) {
	  for (int j=1;j<=m_depth;j++) 
	    res[j] = Max(m_bins[i][j],m_bins[i-1][j]);
	}
	*/

      }
    }
  }
}

double Histogram::Integral() const 
{
  double total=0.;
  for (int i=0;i<m_nbin;i++) { 
    total += m_yvalues[i]; 
  }
  return total*m_binsize;
}

double Histogram::Ymax() const
{
  double ymax=m_yvalues[1];
  for (int i=1;i<m_nbin-1;i++) { 
    if (ymax<m_yvalues[i]) ymax=m_yvalues[i];
  }
  return ymax;
}

double Histogram::Ymin() const
{
  double ymin=1.e+65;
  for (int i=1;i<m_nbin-1;i++) { 
    if (ymin>m_yvalues[i] && m_yvalues[i]!=0) ymin=m_yvalues[i];
  }
  return ymin;
}

double Histogram::LogCoeff() const
{
  double ymax=m_yvalues[1];
  double ymin=1.e+65;
  double meany,meany2,meanly,meanly2;
  meany=meany2=meanly=meanly2=0.;

  int nl = 0;
  for (int i=1;i<m_nbin-1;i++) { 
    if (ymax<m_yvalues[i]) ymax=m_yvalues[i];
    if (ymin>m_yvalues[i] && m_yvalues[i]!=0.) ymin=m_yvalues[i];
    double y=m_yvalues[i];
    if (y!=0.) {
      meany   += y;
      meany2  += y*y;
      meanly  += log(y);
      meanly2 += sqr(log(y));
      ++nl;
    }
  }
  double rl = 0.;
  double rn = 0.;
  if (ymax!=0. && ymin!=0. && nl!=0) {
    double lymax=log(ymax);
    double lymin=log(ymin);
    meanly  = meanly/nl;
    meanly2 = meanly2/nl;
    double sigl2=meanly2-sqr(meanly);
    double ly0=0.5*(lymax+lymin);
    if (sigl2!=0.) rl=sigl2/sqr(ly0-meanly);
  }
  if (nl!=0) {
    meany   = meany/nl;
    meany2  = meany2/nl;
    double sig2=meany2-sqr(meany);
    double y0=0.5*(ymax+ymin);
    if (sig2!=0) rn=sig2/sqr(y0-meany);
  }
  double r;
  if (rl==0. && rn==0.) r=1.;
  else if (rl==0.)      r=0.;
  else if (rn==0.)      r=20.;
  else r=rl/rn;
  return r;
}



Histogram & Histogram::operator+=(const Histogram & histo)
{
  if (histo.m_nbin!=m_nbin) {
    msg_Error()<<"Error in Histogram : can not add histograms with different number of bins"<<std::endl;
    return *this;
  }
  for (int i=0;i<m_nbin;i++) { 
    m_yvalues[i]+= histo.m_yvalues[i]; 
  }
  if (m_depth>1) {
    for (int i=0;i<m_nbin;i++) { 
      m_y2values[i]+= histo.m_y2values[i]; 
    }
  }
  if (m_depth>2) {
    for (int i=0;i<m_nbin;i++) { 
      m_psvalues[i]+= histo.m_psvalues[i]; 
    }
  }
  
  m_fills+=histo.m_fills;
  m_psfills+=histo.m_psfills;
  return *this;
}


void Histogram::Addopt(const Histogram & histo)
{
  
  if (m_depth<=1) {
    msg_Error()<<"Error in Histogram : can not Addopt histograms without statistical errors"<<std::endl;
    return;
  }
  if (histo.m_nbin!=m_nbin) {
    msg_Error()<<"Error in Histogram : can not add histograms with different number of bins"<<std::endl;
    return;
  }
  for (int i=0;i<m_nbin;i++) { 
    double w1=sqr(m_yvalues[i])/m_y2values[i];
    double w2=sqr(histo.m_yvalues[i])/histo.m_y2values[i];
    if (!(w1>0. && w2>0.)) w1=w2=1.;
    m_yvalues[i]=(m_yvalues[i]*w1+histo.m_yvalues[i]*w2)/(w1+w2); 
    m_y2values[i]=sqr(m_yvalues[i])/(w1+w2); 

    if (m_depth>2) {
      m_psvalues[i]+= histo.m_psvalues[i];       
    }
  }  
  m_fills+=histo.m_fills;
  m_psfills+=histo.m_psfills;
  return;
}
