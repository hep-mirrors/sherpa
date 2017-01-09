#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Default_Reader.H"
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

Histogram_2D::Histogram_2D(int type,double xmin, double xmax, int nbinx,
                                    double ymin, double ymax, int nbiny) :
  m_type(type), m_nbinx(nbinx), m_nbiny(nbiny),
  m_lowerx(xmin), m_upperx(xmax), m_lowery(ymin), m_uppery(ymax),
  m_zvalues(0),m_z2values(0),
  m_psvalues(0), m_tmp(0), m_fills(0), m_psfills(0),
  m_finished(false), m_initialized(false), m_fuzzyexp(-1)
{
  m_mcb = 0.;
  if (m_type>10000) {
    m_type-=10000;
    m_fuzzyexp = int(m_type/1000);
    m_type-=1000*m_fuzzyexp;
  }
  if (m_type>=1000) {
    m_mcb = 1.;
    m_type-=1000;
  }
  m_logarithmicx = int(m_type/100);
  m_logarithmicy = int((m_type-m_logarithmicx*100)/10);
  m_depth        = m_type-m_logarithmicx*100-m_logarithmicy*10+1;
  m_logbasex = m_logbasey = 1;
  switch(m_logarithmicx) {
    case 1:
      m_logbasex = log(10.);
      m_upperx   = log(m_upperx)/m_logbasex;
      m_lowerx   = log(m_lowerx)/m_logbasex;
      break;
    case 2:
      m_upperx   = log(m_upperx); m_lowerx = log(m_lowerx);
      break;
    default: break;
  }
  switch(m_logarithmicy) {
    case 1:
      m_logbasey = log(10.);
      m_uppery   = log(m_uppery)/m_logbasey;
      m_lowery   = log(m_lowery)/m_logbasey;
      break;
    case 2:
      m_uppery   = log(m_uppery); m_lowery = log(m_lowery);
      break;
    default: break;
  }
  m_binsizex = (m_upperx-m_lowerx)/(double(m_nbinx));
  m_binsizey = (m_uppery-m_lowery)/(double(m_nbiny));

  if (m_binsizex<=0. || m_binsizey<=0.) {
    msg_Error()<<"Error in Histogram_2D : Tried to initialize a "
               <<"histogram with binsize <= 0 ! :"
               <<m_binsizex<<" , "<<m_binsizey<<std::endl;
    m_active = 0;
    return;
  }

  m_nbin = m_nbinx*m_nbiny+2;
  m_active = 1;
  m_zvalues = new double[m_nbin];
  for (int i=0;i<m_nbin;i++) {
    m_zvalues[i]=0.;
  }
  if (m_depth>1) {
    m_z2values   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_z2values[i]=0.;
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
}

Histogram_2D::Histogram_2D(const Histogram_2D * histo) :
    m_zvalues(0), m_z2values(0), m_psvalues(0), m_tmp(0) {
  m_lowerx  = histo->m_lowerx;
  m_upperx  = histo->m_upperx;
  m_lowery  = histo->m_lowery;
  m_uppery  = histo->m_uppery;
  m_logbasex = histo->m_logbasex;
  m_logbasey = histo->m_logbasey;
  m_logarithmicx = histo->m_logarithmicx;
  m_logarithmicy = histo->m_logarithmicy;
  m_mcb     = histo->m_mcb;
  m_nbin    = histo->m_nbin;
  m_nbinx   = histo->m_nbinx;
  m_nbiny   = histo->m_nbiny;
  m_depth   = histo->m_depth;
  m_type    = histo->m_type;
  m_fills   = histo->m_fills;
  m_psfills = histo->m_psfills;

  m_binsizex = histo->m_binsizex;
  m_binsizey = histo->m_binsizey;
  m_active  = 1;
  m_finished = histo->m_finished;

  m_zvalues   = new double[m_nbin];
  for (int i=0;i<m_nbin;i++) {
    m_zvalues[i]  = histo->m_zvalues[i];
  }
  if (m_depth>1) {
    m_z2values   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_z2values[i]=histo->m_z2values[i];
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


Histogram_2D::Histogram_2D(const std::string & pID) :
    m_zvalues(0), m_z2values(0), m_psvalues(0),
    m_tmp(0), m_fills(0), m_mcb(0.)  {
  m_finished=true;
  std::ifstream ifile(pID.c_str());

  std::string dummy;
  getline(ifile,dummy);

  Default_Reader reader;
  if (dummy!="") {
    std::vector<std::string> conf;
    reader.SetString(dummy);
    reader.VectorFromString(conf);
    size_t k=0;

    if (k>=conf.size()) {
      msg_Error()<<"Error in Histogram_2D : reading file :"<<pID<<std::endl;
      m_active = 0;
      return;
    }

    MyStrStream str;
    str<<dummy;
    str>>m_type>>m_nbin>>m_nbinx>>m_lowerx>>m_upperx
                       >>m_nbiny>>m_lowery>>m_uppery;
    k=8;

    if (m_type>10000) {
      m_type-=10000;
      m_fuzzyexp = int(m_type/1000);
      m_type-=1000*m_fuzzyexp;
    }
    if (m_type>=1000) {
      m_mcb = 1.;
      m_type-=1000;
    }
    m_logarithmicx = int(m_type/100);
    m_logarithmicy = int((m_type-m_logarithmicx*100)/10);
    m_depth        = m_type-m_logarithmicx*100-m_logarithmicy*10+1;

    m_logbasex = m_logbasey = 1;
    switch(m_logarithmicx) {
      case 1:
        m_logbasex = log(10.);
        m_upperx   = log(m_upperx)/m_logbasex;
        m_lowerx   = log(m_lowerx)/m_logbasex;
        break;
      case 2:
        m_upperx   = log(m_upperx); m_lowerx = log(m_lowerx);
        break;
      default: break;
    }
    switch(m_logarithmicy) {
      case 1:
        m_logbasey = log(10.);
        m_uppery   = log(m_uppery)/m_logbasey;
        m_lowery   = log(m_lowery)/m_logbasey;
        break;
      case 2:
        m_uppery   = log(m_uppery); m_lowery = log(m_lowery);
        break;
      default: break;
    }
    m_binsizex = (m_upperx-m_lowerx)/(double(m_nbinx));
    m_binsizey = (m_uppery-m_lowery)/(double(m_nbiny));

    if (m_binsizex<=0. || m_binsizey<=0.) {
      msg_Error()<<"Error in Histogram_2D : Tried to initialize a "
                 <<"histogram with  binsize <= 0 ! :"
                 <<m_binsizex<<" , "<<m_binsizey<<std::endl;
      m_active = 0;
      return;
    }

    m_active = 1;
    m_zvalues   = new double[m_nbin];
    m_zvalues[0]  = Get<double>(conf[k++]);
    if (m_depth>1) {
      m_z2values   = new double[m_nbin];
      m_z2values[0] = sqr(Get<double>(conf[k++]));
    }
    if (m_depth>2) {
      m_psvalues   = new double[m_nbin];
      m_psvalues[0] = Get<double>(conf[k++]);
    }
    if (k>=conf.size()) {
      msg_Error()<<"Error in Histogram_2D : reading file :"<<pID<<std::endl;
      m_active = 0;
      return;
    }

    m_zvalues[m_nbin-1]  = Get<double>(conf[k++]);
    if (m_depth>1) {
      m_z2values[m_nbin-1] = sqr(Get<double>(conf[k++]));
    }
    if (m_depth>2) {
      m_psvalues[m_nbin-1] = Get<double>(conf[k++]);
    }
    if (k>=conf.size()) {
      msg_Error()<<"Error in Histogram_2D : reading file :"<<pID<<std::endl;
      m_active = 0;
      return;
    }
    m_fills = Get<double>(conf[k++]);
  }
  else {
    msg_Error()<<"Error in Histogram_2D : reading file :"<<pID<<std::endl;
    m_active = 0;
    return;
  }


  std::vector<std::string> data;
  MyStrStream str;
  for (int i=0;i<m_nbin-1;i++) {
    getline(ifile,dummy);
    data.clear();
    reader.SetString(dummy);
    reader.VectorFromString(data);

    m_zvalues[i+1] = Get<double>(data[2]);
    if (m_depth>1) {
      m_z2values[i+1] = Get<double>(data[3]);
      m_z2values[i+1] = sqr(m_z2values[i+1]);
    }
    if (m_depth>2) {
      m_psvalues[i+1] = Get<double>(data[4]);
    }
  }
  ifile.close();
}


Histogram_2D::~Histogram_2D() {
  if (m_zvalues!=0) {
    delete [] m_zvalues; m_zvalues = 0;
  }
  if (m_z2values!=0) {
    delete [] m_z2values; m_z2values = 0;
  }
  if (m_psvalues!=0) {
    delete [] m_psvalues; m_psvalues = 0;
  }
  if (m_tmp!=0) {
    delete [] m_tmp; m_tmp = 0;
  }
}


void Histogram_2D::Finalize() {
  if (!m_finished) {
    m_finished=true;
    if (m_fills==0.) return;
    for (int i=0;i<m_nbin;++i) {
      m_zvalues[i]/=m_fills*m_binsizex*m_binsizey;
      if (m_depth>1) {
        m_z2values[i]/=m_fills*sqr(m_binsizex*m_binsizey);
        if (m_fills>1)
          m_z2values[i]=(m_z2values[i]-sqr(m_zvalues[i]))/(m_fills-1);
      }
    }
    if (m_depth>2) {
      double itg = Integral()/(m_psfills*m_binsizex*m_binsizey);
      for (int i=0;i<m_nbin;++i) {
        m_psvalues[i]*=itg;
      }
    }
  }
}

void Histogram_2D::Restore() {
  if (m_finished) {
    for (int i=0;i<m_nbin;++i) {
      if (m_depth>1) {
        if (m_fills>1)
          m_z2values[i]=(m_fills-1)*m_z2values[i]+sqr(m_zvalues[i]);
        m_z2values[i]*=m_fills*sqr(m_binsizex*m_binsizey);
        if (m_depth>2)
          m_psvalues[i]*=m_psfills*m_binsizex*m_binsizey;
      }
      m_zvalues[i]*=m_fills*m_binsizex*m_binsizey;
    }
    m_finished=false;
  }
}

double Histogram_2D::Mean() const
{
  double sum(0.), range(0.);
  int index(0);
  for (int i=0;i<m_nbinx;++i) {
    for (int j=0;j<m_nbiny;++j) {
      index++;
      double widthx((m_upperx-m_lowerx)/m_nbinx);
      double widthy((m_uppery-m_lowery)/m_nbiny);
      if (m_logarithmicx)
        widthx=pow(m_logbasex,m_lowerx+i*widthx)-pow(m_logbasex,m_lowerx+(i-1)*widthx);
      if (m_logarithmicy)
        widthy=pow(m_logbasey,m_lowery+j*widthy)-pow(m_logbasey,m_lowery+(j-1)*widthy);
      sum+=m_zvalues[index]*widthx*widthy;
      range+=widthx*widthy;
    }
  }
  return sum/range;
}

void Histogram_2D::Reset() {
  for (int i=0;i<m_nbin;i++) {
    m_zvalues[i]=0.;
    if (m_depth>1) {
      m_z2values[i]=0.;
    }
    if (m_depth>2) {
      m_psvalues[i]=0.;
    }
  }
  m_fills=0;
  m_psfills=0;
}

void Histogram_2D::Scale(double scale) {
  for (int i=0;i<m_nbin;i++) {
    m_zvalues[i]*= scale;
    if (m_depth>1) m_z2values[i]*=sqr(scale);
    if (m_depth>2) m_psvalues[i]*=scale;
  }
}

void Histogram_2D::Output() {
  if (!msg_LevelIsDebugging()) return;
  msg_Out()<<"----------------------------------------"<<std::endl
           <<"    "<<m_zvalues[0]<<std::endl
           <<"----------------------------------------"<<std::endl;
  double result(0.);
  int index(0);
  for (int i=0;i<m_nbinx;++i) {
    for (int j=0;j<m_nbiny;++j) {
      index++;
      msg_Out()<<m_lowerx+i*m_binsizex<<"  ";
      msg_Out()<<m_lowery+j*m_binsizey<<"  ";
      msg_Out()<<m_zvalues[index]<<"  ";
      if (m_depth>1) msg_Out()<<sqrt(m_z2values[index]);
      result += m_zvalues[index];
      msg_Out()<<std::endl;
    }
  }
  msg_Out()<<m_lowerx+m_nbinx*m_binsizex<<" == "<<m_upperx<<std::endl;
  msg_Out()<<m_lowery+m_nbiny*m_binsizey<<" == "<<m_uppery<<std::endl
           <<"----------------------------------------"<<std::endl
           <<"    "<<m_zvalues[m_nbin-1]<<std::endl
           <<"----------------------------------------"<<std::endl
           <<"Inside the range : "<<result<<std::endl;
}

void Histogram_2D::Output(const std::string name) {
  if (!m_active) return;
  std::ofstream ofile;
  ofile.open(name.c_str());

  if (m_fills>=0) {
    ofile<<m_type<<" "<<m_nbin<<" "
         <<m_nbinx<<" "<<m_lowerx<<" "<<m_upperx<<" "
         <<m_nbiny<<" "<<m_lowery<<" "<<m_uppery<<" ";
    ofile<<m_zvalues[0]<<"  ";
    if (m_depth>1) ofile<<m_z2values[0]<<"  ";
    ofile<<m_zvalues[m_nbin-1]<<"  ";
    if (m_depth>1) ofile<<m_z2values[m_nbin-1]<<"  ";
    ofile<<m_fills<<"\n";
  }
  int index(0);
  for (int i=0;i<m_nbinx;++i) {
    for (int j=0;j<m_nbiny;++j) {
      index++;
      ofile<<m_lowerx+i*m_binsizex<<"  ";
      ofile<<m_lowery+j*m_binsizey<<"  ";
      ofile<<m_zvalues[index]<<"  ";
      if (m_depth>1) ofile<<sqrt(m_z2values[index])<<"  ";
      if (m_depth>2) ofile<<m_psvalues[index]<<"  ";
      ofile<<"\n";
    }
  }
  ofile.close();
}

void Histogram_2D::Insert(double coordinatex, double coordinatey) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram_2D : Tried to access a "
               <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }

  m_fills++;

  if (m_logarithmicx>0) coordinatex = log(coordinatex)/m_logbasex;
  if (m_logarithmicy>0) coordinatey = log(coordinatey)/m_logbasey;
  if (coordinatex<m_lowerx) { m_zvalues[0] += double(1); return; }
  if (coordinatey<m_lowery) { m_zvalues[0] += double(1); return; }
  if (coordinatex>m_upperx) { m_zvalues[m_nbin-1] += double(1); return; }
  if (coordinatey>m_uppery) { m_zvalues[m_nbin-1] += double(1); return; }
  int index(0);
  for (int i=0;i<m_nbinx;++i) {
    for (int j=0;j<m_nbiny;++j) {
      index++;
      if ((coordinatex >= m_lowerx + i*m_binsizex) &&
          (coordinatex <  m_lowerx + (i+1)*m_binsizex) &&
          (coordinatey >= m_lowery + j*m_binsizey) &&
          (coordinatey <  m_lowery + (j+1)*m_binsizey) ) {
        m_zvalues[index] += double(1);
        return;
      }
    }
  }
}

void Histogram_2D::Insert(int ix, int iy,double value,double ncount) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram_2D : Tried to access a "
               <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  m_fills+=ncount;
  if (value==0.) return;
  m_psfills++;

  if (ix<0 || iy<0) {
    m_zvalues[0] += value;
    if (m_depth>1) {
      if (value>m_z2values[0]) m_z2values[0] = value;
      if (m_depth>2) m_psvalues[0] += 1.;
    }
    return;
  }

  if (ix>=m_nbinx || iy>=m_nbiny) {
    m_zvalues[m_nbin-1] += value;
    if (m_depth>1) {
      if (value>m_z2values[m_nbin-1]) m_z2values[m_nbin-1] = value;
      if (m_depth>2) m_psvalues[m_nbin-1] += 1.;
    }
    return;
  }

  int index(iy+ix*m_nbiny+1);
  m_zvalues[index] += value;
  if (m_depth>1) {
    m_z2values[index] += value*value;
    if (m_depth>2) m_psvalues[index] += 1.;
  }
}

void Histogram_2D::Insert(double coordinatex,double coordinatey,
                          double value,double ncount) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram_2D : Tried to access a "
               <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  if (m_logarithmicx>0) coordinatex = log(coordinatex)/m_logbasex;
  if (m_logarithmicy>0) coordinatey = log(coordinatey)/m_logbasey;
  int binx((int)((coordinatex-m_lowerx)/m_binsizex));
  int biny((int)((coordinatey-m_lowery)/m_binsizey));
  Insert(binx,biny,value,ncount);

  if (m_fuzzyexp<0) return;

  double x = (coordinatex-m_lowerx)/m_binsizex-double(binx)+0.5;
  double y = (coordinatey-m_lowery)/m_binsizey-double(biny)+0.5;
  if (binx==0&&x<0.) return;
  if (biny==0&&y<0.) return;
  if (binx==m_nbinx&&x>0.) return;
  if (biny==m_nbiny&&x>0.) return;
  double ff=1.;
  if (m_fuzzyexp==0) ff=0.5;
  if (m_fuzzyexp>0) ff=0.5*pow(2.*dabs(x),m_fuzzyexp);
  if (m_fuzzyexp==9) ff=0.5*sqrt(2.*dabs(x));

  int bin(biny+binx*m_nbiny+1);
  m_zvalues[bin] -= ff*value;
  if (m_depth>1) {
    m_z2values[bin] += sqr(ff*(value))-sqr(value);
    if (m_depth>2) m_psvalues[bin] -= ff;
  }

  if (x>0.) {
    m_zvalues[biny+(binx+1)*m_nbiny+1] += ff*value;
    if (m_depth>1) {
      m_z2values[biny+(binx+1)*m_nbiny+1] += sqr(ff*value);
      if (m_depth>2) m_psvalues[biny+(binx+1)*m_nbiny+1] += ff;
    }
  }
  if (y>0.) {
    m_zvalues[bin+1] += ff*value;
    if (m_depth>1) {
      m_z2values[bin+1] += sqr(ff*value);
      if (m_depth>2) m_psvalues[bin+1] += ff;
    }
  }
  if (x<0.) {
    m_zvalues[biny+(binx-1)*m_nbiny+1] += ff*value;
    if (m_depth>1) {
      m_z2values[biny+(binx-1)*m_nbiny+1] += sqr(ff*value);
      if (m_depth>2) m_psvalues[biny+(binx-1)*m_nbiny+1] += ff;
    }
  }
  if (y<0.) {
    m_zvalues[bin-1] += ff*value;
    if (m_depth>1) {
      m_z2values[bin-1] += sqr(ff*value);
      if (m_depth>2) m_psvalues[bin-1] += ff;
    }
  }
}

void Histogram_2D::InsertMCB(double coordinatex,double coordinatey,
                             double value,double ncount) {
  if (!m_tmp) {
    m_tmp   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_tmp[i]=0.;
    }
  }
  m_mcb = ncount;

  if (m_logarithmicx>0) coordinatex = log(coordinatex)/m_logbasex;
  if (m_logarithmicy>0) coordinatey = log(coordinatey)/m_logbasey;

  int binx(int((coordinatex-m_lowerx)/m_binsizex));
  int biny(int((coordinatey-m_lowery)/m_binsizey));
  int bin(biny+binx*m_nbiny+1);
  if (bin<0) bin=0;
  if (bin>=m_nbin) bin=m_nbin-1;
  if (bin==0||bin==m_nbin-1) {
    m_tmp[bin] += value;
    return;
  }

  double x = (coordinatex-m_lowerx)/m_binsizex-double(binx)+0.5;
  double y = (coordinatey-m_lowery)/m_binsizey-double(biny)+0.5;
  if ((bin==1&&x<0.)||(bin==m_nbin-2&&x>0.)) {
    m_tmp[bin] += value;
    return;
  }
  double ff=1.;
  if (m_fuzzyexp==0) ff=0.5;
  if (m_fuzzyexp>0) ff=1.-0.5*pow(2.*dabs(x),m_fuzzyexp);
  if (m_fuzzyexp==9) ff=1.-0.5*sqrt(2.*dabs(x));

  m_tmp[bin] += ff*value;
  if (x>0.) m_tmp[biny+(binx+1)*m_nbiny+1] += (1.-ff)*value;
  if (y>0.) m_tmp[bin+1] += (1.-ff)*value;
  if (x<0.) m_tmp[biny+(binx-1)*m_nbiny+1] += (1.-ff)*value;
  if (y<0.) m_tmp[bin-1] += (1.-ff)*value;
}

void Histogram_2D::InsertMCBIM(double coordinatex,double coordinatey,
                               double value) {
  if (!m_tmp) {
    m_tmp   = new double[m_nbin];
    for (int i=0;i<m_nbin;i++) {
      m_tmp[i]=0.;
    }
  }
  m_mcb = 1.;
  if (m_logarithmicx>0) coordinatex = log(coordinatex)/m_logbasex;
  if (m_logarithmicy>0) coordinatey = log(coordinatey)/m_logbasey;

  int binx(int((coordinatex-m_lowerx)/m_binsizex));
  int biny(int((coordinatey-m_lowery)/m_binsizey));
  int bin(biny+binx*m_nbiny+1);
  if (bin<0) bin=0;
  if (bin>=m_nbin) bin=m_nbin-1;
  for (int i=bin+1;i<m_nbin;i++) m_tmp[i] += value;
}

void Histogram_2D::FinishMCB()
{
  m_fills+=m_mcb;

  for (int i=0;i<m_nbin;i++) {
    m_zvalues[i] += m_tmp[i];
    if (m_depth>1) {
      m_z2values[i] += m_tmp[i]*m_tmp[i];
      if (m_depth>2) m_psvalues[i] += 1.;
    }
    m_tmp[i] = 0.;
  }
}

void Histogram_2D::InsertRange(double startx, double endx,
                               double starty, double endy, double value) {
  if (!m_active) {
    msg_Error()<<"Error in Histogram_2D : Tried to access a "
               <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  if (m_logarithmicx>0) {
    if (startx>0)
      startx = log(startx)/m_logbasex;
    else
      startx = -30;
    if (endx>0)
      endx = log(endx)/m_logbasex;
    else
      endx = -30;
  }
  if (m_logarithmicy>0) {
    if (starty>0)
      starty = log(starty)/m_logbasey;
    else
      starty = -30;
    if (endy>0)
      endy = log(endy)/m_logbasey;
    else
      endy = -30;
  }
  m_fills++;

  // underrun
  if (startx<m_lowerx) {
    startx = m_lowerx;
    m_zvalues[0] += value;
    if (m_depth>1) {
      if (value>m_z2values[0]) m_z2values[0] = value;
    }
  }
  if (starty<m_lowery) {
    starty = m_lowery;
    m_zvalues[0] += value;
    if (m_depth>1) {
      if (value>m_z2values[0]) m_z2values[0] = value;
    }
  }
  if (endx<=m_lowerx || endy<=m_lowery) return;

  // overflow
  if (endx>m_upperx) {
    endx=m_upperx;
    m_zvalues[m_nbin-1] += value;
    if (m_depth>1) {
      if (value>m_z2values[m_nbin-1]) m_z2values[m_nbin-1] = value;
    }
  }
  if (endy>m_uppery) {
    endy=m_uppery;
    m_zvalues[m_nbin-1] += value;
    if (m_depth>1) {
      if (value>m_z2values[m_nbin-1]) m_z2values[m_nbin-1] = value;
    }
  }
  if (startx>=m_upperx || starty>=m_uppery) return;

  double lowx,upx,lowy,upy;
  lowx = m_lowerx; upx = m_lowerx+m_binsizex;
  lowy = m_lowery; upy = m_lowery+m_binsizey;
  int index(0);
  for (int i=1;i<m_nbinx;++i) {
    for (int j=1;j<m_nbiny;++j) {
      index++;
      if ((startx < upx) && (endx >= lowx) &&
          (starty < upy) && (endy >= lowy)) {
        if ((startx<=lowx) && (upx<=endx) &&
            (starty<=lowy) && (upy<=endy)) {
          double newval((std::max(startx,lowx)-std::min(endx,upx))/m_binsizex
                        *(std::max(starty,lowy)-std::min(endy,upy))/m_binsizex
                        *value);
          m_zvalues[index] += newval;
          if (m_depth>1) m_z2values[index] += sqr(newval);
          if (m_depth>2) m_psvalues[index] += newval;
        }
      }
      lowy = upy;
      upy += m_binsizey;
    }
    lowx = upx;
    upx += m_binsizex;
    lowy = m_lowery;
    upy  = m_lowery+m_binsizey;
  }
}

double Histogram_2D::Bin(double coordinatex, double coordinatey)
{
  if (!m_active) {
    msg_Error()<<"Error in Histogram_2D : Tried to access a histogram wit binsize <= 0 ! Return 0.."<<std::endl;
    return -1.0;
  }
  else {
    if (m_logarithmicx>0) coordinatex = log(coordinatex)/m_logbasex;
    if (m_logarithmicy>0) coordinatey = log(coordinatey)/m_logbasey;

    if (coordinatex<m_lowerx) return m_zvalues[0];
    if (coordinatey<m_lowery) return m_zvalues[0];
    if (coordinatex>m_upperx) return m_zvalues[m_nbin-1];
    if (coordinatey>m_uppery) return m_zvalues[m_nbin-1];
    int index(0);
    for (int i=0;i<m_nbinx;++i) {
      for (int j=0;j<m_nbiny;++j) {
        index++;
        if ((coordinatex >= m_lowerx+i*m_binsizex) &&
            (coordinatey >= m_lowery+j*m_binsizey) &&
            (coordinatex <  m_lowerx + (i+1)*m_binsizex) &&
            (coordinatey <  m_lowery + (j+1)*m_binsizey))
          return m_zvalues[index];
      }
    }
  }
  return -1.0;
}

double Histogram_2D::Integral() const
{
  double total=0.;
  for (int i=1;i<m_nbin-1;i++) {
    total += m_zvalues[i];
  }
  return total*m_binsizex*m_binsizey;
}

double Histogram_2D::Integral(int xminbin, int xmaxbin,
                              int yminbin, int ymaxbin) const
{
  double total=0.;
  int index(0);
  for (int i=0;i<m_nbinx;++i) {
    for (int j=0;j<m_nbiny;++j) {
      index++;
      if (i>=xminbin && i<xmaxbin && j>=yminbin && j<ymaxbin)
        total += m_zvalues[i];
    }
  }
  return total*m_binsizex*m_binsizey;
}


double Histogram_2D::Zmax() const
{
  double ymax=m_zvalues[1];
  for (int i=1;i<m_nbin-1;i++) {
    if (ymax<m_zvalues[i]) ymax=m_zvalues[i];
  }
  return ymax;
}

double Histogram_2D::Zmin() const
{
  double ymin=1.e+65;
  for (int i=1;i<m_nbin-1;i++) {
    if (ymin>m_zvalues[i] && m_zvalues[i]!=0) ymin=m_zvalues[i];
  }
  return ymin;
}

double Histogram_2D::LogCoeff() const
{
  double zmax(m_zvalues[1]);
  double zmin(1.e+65);
  double meanz(0.),meanz2(0.),meanlz(0.),meanlz2(0.);

  int nl(0);
  for (int i=1;i<m_nbin-1;i++) {
    if (zmax<m_zvalues[i]) zmax=m_zvalues[i];
    if (zmin>m_zvalues[i] && m_zvalues[i]!=0.) zmin=m_zvalues[i];
    double z(m_zvalues[i]);
    if (z!=0.) {
      meanz   += z;
      meanz2  += z*z;
      meanlz  += log(z);
      meanlz2 += sqr(log(z));
      ++nl;
    }
  }
  double rl(0.),rn(0.);
  if (zmax!=0. && zmin!=0. && nl!=0) {
    double lzmax(log(zmax));
    double lzmin(log(zmin));
    meanlz  = meanlz/nl;
    meanlz2 = meanlz2/nl;
    double sigl2(meanlz2-sqr(meanlz));
    double lz0(0.5*(lzmax+lzmin));
    if (sigl2!=0.) rl=sigl2/sqr(lz0-meanlz);
  }
  if (nl!=0) {
    meanz   = meanz/nl;
    meanz2  = meanz2/nl;
    double sig2(meanz2-sqr(meanz));
    double z0(0.5*(zmax+zmin));
    if (sig2!=0) rn=sig2/sqr(z0-meanz);
  }
  double r;
  if (rl==0. && rn==0.) r=1.;
  else if (rl==0.)      r=0.;
  else if (rn==0.)      r=20.;
  else r=rl/rn;
  return r;
}



Histogram_2D & Histogram_2D::operator+=(const Histogram_2D & histo)
{
  if (histo.m_nbinx!=m_nbinx && histo.m_nbiny!=m_nbiny) {
    msg_Error()<<"Error in Histogram_2D : can not add histograms with "
               <<"different number of bins"<<std::endl;
    return *this;
  }
  for (int i=0;i<m_nbin;i++) {
    m_zvalues[i]+= histo.m_zvalues[i];
  }
  if (m_depth>1) {
    for (int i=0;i<m_nbin;i++) {
      m_z2values[i]+= histo.m_z2values[i];
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


void Histogram_2D::Addopt(const Histogram_2D & histo)
{
  if (m_depth<=1) {
    msg_Error()<<"Error in Histogram_2D : can not Addopt histograms "
               <<"without statistical errors"<<std::endl;
    return;
  }
  if (histo.m_nbinx!=m_nbinx && histo.m_nbiny!=m_nbiny) {
    msg_Error()<<"Error in Histogram_2D : can not add histograms "
               <<"with different number of bins"<<std::endl;
    return;
  }
  for (int i=0;i<m_nbin;i++) {
    double w1=sqr(m_zvalues[i])/m_z2values[i];
    double w2=sqr(histo.m_zvalues[i])/histo.m_z2values[i];
    if (!(w1>0. && w2>0.)) w1=w2=1.;
    m_zvalues[i]=(m_zvalues[i]*w1+histo.m_zvalues[i]*w2)/(w1+w2);
    m_z2values[i]=sqr(m_zvalues[i])/(w1+w2);

    if (m_depth>2) {
      m_psvalues[i]+= histo.m_psvalues[i];
    }
  }
  m_fills+=histo.m_fills;
  m_psfills+=histo.m_psfills;
  return;
}
