#include "Histogram.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>

using namespace AMATOOLS;
using namespace AORGTOOLS;

Histogram::Histogram(int _type,double _lower,double _upper,int _nbin) :
  m_type(_type), m_nbin(_nbin), m_lower(_lower), m_upper(_upper), m_bins(0), m_fills(0)
{
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
    msg.Error()<<"Error in Histogram : Tried to initialize a histogram with  binsize <= 0 !"<<std::endl;
    m_active = 0;
    return;
  }

  m_nbin += 2;
  m_active = 1;
  m_bins   = new double * [m_nbin];

  for (int i=0;i<m_nbin;i++) {
    m_bins[i]   = new double[m_depth];
    for (int j=0;j<m_depth;j++) m_bins[i][j] = 0.;
  }
}

Histogram::Histogram(Histogram * histo) {
  m_lower   = histo->m_lower;
  m_upper   = histo->m_upper;
  m_logbase = histo->m_logbase;
  m_logarithmic = histo->m_logarithmic;
  m_nbin   = histo->m_nbin;
  m_depth   = histo->m_depth;
  m_type    = histo->m_type;
  m_fills   = histo->m_fills;

  m_binsize = histo->m_binsize;
  m_active  = 1;

  m_bins    = new double*[m_nbin];
  for (int i=0;i<m_nbin;i++) {
    m_bins[i]  = new double[m_depth]; 
    for (int j=0;j<m_depth;j++) m_bins[i][j] = histo->m_bins[i][j];
  }
}


Histogram::Histogram(std::string _pID) {
  char pID[100];
  sprintf(pID,"%s",_pID.c_str());

  std::ifstream ifile(pID);

  int    _type, _nbins;
  double _lower, _upper;
  ifile>>_type>>_nbins>>_lower>>_upper;

  m_type = _type; m_nbin = _nbins; m_lower = _lower; m_upper = _upper;

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
    msg.Error()<<"Error in Histogram : "
	       <<"Tried to initialize a histogram with m_binsize <= 0 !"
	       <<std::endl;
    m_active = 0;
    return;
  }

  //  m_nbin += 2;
  m_active = 1;
  m_bins   = new double * [m_nbin];

  for (int i=0;i<m_nbin;i++) m_bins[i]   = new double[m_depth];

  for (int j=0;j<m_depth;j++) ifile>>m_bins[0][j];
  for (int j=0;j<m_depth;j++) ifile>>m_bins[m_nbin-1][j];

  double value;
  for (int i=0;i<m_nbin-1;i++) {
    ifile>>value;
    for (int j=0;j<m_depth;j++) ifile>>m_bins[i+1][j];
  }
  ifile.close();
}


Histogram::~Histogram() {
    
  if (m_bins!=0) { 
    for (int i=0;i<m_nbin;i++) 
      delete [] m_bins[i];
    delete [] m_bins; m_bins = 0; 
  }
    
};


void Histogram::Finalize() {
  double total = 0;
  for (int i=0;i<m_nbin;i++) { 
    //    total += m_bins[i][0]; 
    total += m_bins[i][0]/=m_fills; 
  }
}

void Histogram::Reset() {
  for (int i=0;i<m_nbin;i++) { 
    m_bins[i][0]=0;
  }
  m_fills=0;
}

void Histogram::Scale(double scale) {
  double total = 0;
  for (int i=0;i<m_nbin;i++) { 
    m_bins[i][0]*= scale;
    total += m_bins[i][0]; 
  }
  m_fills = int(double(m_fills)/scale);
}

void Histogram::Output() {
  if (!rpa.gen.Tracking()) return;
  msg.Out()<<"----------------------------------------"<<std::endl
	   <<"    "<<m_bins[0][0]<<std::endl
	   <<"----------------------------------------"<<std::endl;
  double result = 0.;
  for (int i=0;i<m_nbin-2;i++) {
    msg.Out()<<m_lower+i*m_binsize<<"  ";
    for (int j=0;j<m_depth;j++) msg.Out()<<m_bins[i+1][j]<<"  ";
    result += m_bins[i+1][0];
    msg.Out()<<std::endl;
  }
  msg.Out()<<m_lower+(m_nbin-2)*m_binsize<<" == "<< m_upper<<std::endl
	   <<"----------------------------------------"<<std::endl
	   <<"    "<<m_bins[m_nbin-1][0]<<std::endl
	   <<"----------------------------------------"<<std::endl
	   <<"Inside the range : "<<result<<std::endl;
}


void Histogram::Output(std::string name) {
  std::ofstream ofile;
  ofile.open(name.c_str());

  ofile<<m_type<<" "<<m_nbin<<" "<<m_lower<<" "<<m_upper<<" ";
  for (int j=0;j<m_depth;j++) ofile<<m_bins[0][j]<<"  ";
  for (int j=0;j<m_depth;j++) ofile<<m_bins[m_nbin-1][j]<<"  ";
  ofile<<std::endl;
  for (int i=0;i<m_nbin-1;i++) {
    ofile<<m_lower+i*m_binsize<<"  ";
    for (int j=0;j<m_depth;j++) ofile<<m_bins[i+1][j]<<"  ";
    ofile<<std::endl;
  }
  ofile.close();
}







void Histogram::Insert(double coordinate) {
  if (!m_active) {
    msg.Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }

  m_fills++;

  if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;
  if (coordinate<m_lower) { m_bins[0][0]       += double(1); return; }
  if (coordinate>m_upper) { m_bins[m_nbin-1][0] += double(1); return; }
  for (int i=1;i<m_nbin-1;i++) {
    if ( (coordinate >= m_lower + (i-1)*m_binsize) &&
	 (coordinate <  m_lower + i*m_binsize) ) {
      m_bins[i][0] += double(1); 
      return; 
    }
  }
}

void Histogram::Insert(double coordinate,double value) {
  if (!m_active) {
    msg.Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

  m_fills++;

  if (coordinate<m_lower) { 
    m_bins[0][0] += value;
    if (m_depth>1) {
      if (value>m_bins[0][1]) m_bins[0][1] = value;
    }
    return; 
  }

  if (coordinate>m_upper) { 
    m_bins[m_nbin-1][0] += value; 
    if (m_depth>1) {
      if (value>m_bins[m_nbin-1][1]) m_bins[m_nbin-1][1] = value;
    }
    return; 
  }

  double low,up;
  low = m_lower; up = m_lower+m_binsize;
  for (int i=1;i<m_nbin-1;i++) {
    if ( (coordinate >= low) && (coordinate < up) ) {
      m_bins[i][0] += value;
      if (m_depth>1) {
	if (value>m_bins[i][1]) m_bins[i][1] = value;
      }
      return; 
    }
    low = up;
    up += m_binsize;
  }
}

void Histogram::InsertRange(double start, double end, double value) {
  if (!m_active) {
    msg.Error()<<"Error in Histogram : Tried to access a "
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
    m_bins[0][0] += value;
    if (m_depth>1) {
      if (value>m_bins[0][1]) m_bins[0][1] = value;
    }
  }

  // overflow
  if (start>m_upper) { 
    m_bins[m_nbin-1][0] += value; 
    if (m_depth>1) {
      if (value>m_bins[m_nbin-1][1]) m_bins[m_nbin-1][1] = value;
    }
  }

  double low,up;
  int hit=0;
  low = m_lower; up = m_lower+m_binsize;
  for (int i=1;i<m_nbin-1;i++) {
    if ((start < up) && (end >= low) ) {
      double fac=1;
      if ((start<=low)&&(up<=end )) {
	m_bins[i][0] += value;
      } 
      else if ((low<start)&&(up<=end)) {
	fac = (start-low)/m_binsize;
	m_bins[i][0] += value *fac;
      }
      else if ((start<=low)&&(end < up)) {
	fac = (up-end)/m_binsize;
	m_bins[i][0] += value *fac;
      }
      else if ((low<start)&&(end <up)) {
	fac = (end-start)/m_binsize;
	m_bins[i][0] += value *fac;
      }
    

      hit=1;
      if (m_depth>1) {
	if (value>m_bins[i][1]) m_bins[i][1] = value;
      }
    }
    low = up;
    up += m_binsize;
  }
}



double * Histogram::Bin(int bin) { return m_bins[bin]; }

double * Histogram::Bin(double coordinate) { 
  if (!m_active) {
    msg.Error()<<"Error in Histogram : Tried to access a histogram wit binsize <= 0 ! Return 0.."<<std::endl;
    return NULL;
  }
  else {
    if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

    if (coordinate<m_lower) return m_bins[0];
    if (coordinate>m_upper) return m_bins[m_nbin-1];
    for (int i=1;i<m_nbin+1;i++) {
      if ( (coordinate >= m_lower + (i-1)*m_binsize) &&
	   (coordinate <  m_lower + i*m_binsize) ) 
	return m_bins[i];
    }
  }
  return NULL;
}


void Histogram::Extrapolate(double coordinate,double * res,int mode) {
  if (!m_active) {
    msg.Error()<<"Error in Histogram : Tried to access a histogram with binsize <= 0 ! Return 0.."<<std::endl;
  }
  else {
    if (m_logarithmic>0) coordinate = log(coordinate)/m_logbase;

    for (int i=1;i<m_nbin;i++) {
      if ( (coordinate >= m_lower + (i-1)*m_binsize) &&
	   (coordinate <  m_lower + i*m_binsize) ) {
	res[0] = m_bins[i-1][0] +
	  (m_bins[i][0]-m_bins[i-1][0])/m_binsize *
	  (coordinate - m_lower - (i-1) * m_binsize);
	double over = 0.; double under = 0.;
	switch (mode) {
	case 1  : 
	  under = (coordinate - m_lower - (i-1) * m_binsize)/m_binsize * m_bins[i][0];
	  over = (m_lower + (i) * m_binsize - coordinate)/m_binsize * m_bins[i][0];
	  for (int j=0;j<i;j++) {
	    under += m_bins[j][0];
	  }
	  for (int j=i;j<m_nbin-1;j++) {
	    over += m_bins[j+1][0];
	  }
	  res[0]  = over;

	  if (m_depth) {
	    res[1]=0.;
	    for (int j=i;j<m_nbin;j++) {
	      res[1] = AMATOOLS::Max(res[1],m_bins[j][1]);
	    }
	  }


	  break;
	case -1 :
	  for (int j=i-1;j>0;j--) {
	    over  += m_bins[j][0];
	    under += m_bins[j-1][0];
	  }
	  over    += m_bins[0][0];
	  res[0]  += (over+under)/2;
	  break;
	default : break;
	}
	/*
	if (m_depth>0) {
	  for (int j=1;j<=m_depth;j++) 
	    res[j] = AMATOOLS::Max(m_bins[i][j],m_bins[i-1][j]);
	}
	*/

      }
    }
  }
}
