#include "Heavy_Tensor_Current.H"
#include "Tools.H"

using namespace HADRONS;
using namespace ATOOLS;

void Heavy_Tensor_Current::SetModelParameters( struct GeneralModel _md )
{
  m_Vxx = _md("Vxx",0.04);
  switch( int(_md("HT_FORM_FACTOR", 2)+0.5) ) {
  case 1:
    p_ff = new ISGW( _md, p_masses );
    msg.Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  case 2:
    p_ff = new ISGW2( _md, p_masses );
    msg.Tracking()<<"    Using ISGW2 form factor model for "<<m_name<<std::endl;
    break;
  default:
    msg.Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}

void Heavy_Tensor_Current::Calc()
{
  Vec4D pB = p_moms[0];
  Vec4D pT = p_moms[1];
  Vec4D q = pB - pT;
  double q2 = q*q;
  double mB = p_masses[0];
  double mT = p_masses[1];

  // Get formfactors
  p_ff->CalcFFs(q2);
  double h = p_ff->h();
  double k = p_ff->k();
  double bplus = p_ff->bplus();
  double bminus = p_ff->bminus();

  for( int h_had=0; h_had<5; h_had++) {
    CMatrix eps = Tools::ComplexSpin2BosonPolarizationVectorC(p_moms[1], h_had);
    Complex i = Complex(0.0,1.0);
    p_results[h_had] = m_Vxx*( h*i*cross((eps.Conjugate()*pB),pB+pT,pB-pT)
                               - k*(eps.Conjugate()*pB) - bplus*(eps.Conjugate()*pB*pB)*(pB+pT)
                               - bminus*(eps.Conjugate()*pB*pB)*(pB-pT)); 
  }
}



Heavy_Tensor_Current::ISGW::ISGW( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_msb = _md("msb",5.2);
  m_msd = _md("msd",0.33);
  m_bb2 = _md("bb2",0.41*0.41);
  m_msq = _md("msq",1.82);
  m_bx2 = _md("bx2",0.34*0.34);
}

void Heavy_Tensor_Current::ISGW::CalcFFs( double q2 )
{
  double mtb = m_msb + m_msd;
  double mtx = m_msq + m_msd;
 
  double mup = 1.0/(1.0/m_msq+1.0/m_msb);
  double mum = 1.0/(1.0/m_msq-1.0/m_msb);
  double bbx2= 0.5*(m_bb2+m_bx2);

  double tm=(m_mB-m_mT)*(m_mB-m_mT);
  if (q2>tm) q2 = 0.99*tm;
  double kap = 0.7*0.7;

  double f5 = sqrt(mtx/mtb)*pow(sqrt(m_bx2*m_bb2)/bbx2,5.0/2.0)*
       exp(-1.0*((m_msd*m_msd*(tm-q2)/(4.0*mtb*mtx*kap*bbx2))));

  m_h      = f5*(m_msd/(sqrt(8.0*m_bb2)*mtb))*((1.0/m_msq)-(m_msd*m_bb2/(2.0*mum*
               mtx*bbx2)));
  
  m_k      = f5*m_msd*sqrt(2.0/m_bb2);

  m_bplus  = (-1.0*f5*m_msd/(sqrt(8.0*m_bb2)*m_msb*mtx))*(1.0-(m_msd*m_msb*m_bx2/(
             2.0*mup*mtb*bbx2))+(m_msd*m_msb*m_bx2*(1.0-(m_msd*m_bx2/(2.0*mtb*bbx2)))/
            (4.0*mtb*mum*bbx2)));

  m_bminus = 0.0;
 
  m_calced = true;
}


Heavy_Tensor_Current::ISGW2::ISGW2( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_msb = _md("msb",5.2);
  m_msd = _md("msd",0.33);
  m_bb2 = _md("bb2",0.431*0.431);
  m_mbb = _md("msq",5.31);
  m_nf  = _md("bx2",4.0);
  m_msq = _md("msb",1.82);
  m_bx2 = _md("msb",0.33*0.33);
  m_mbx = _md("msb",(5.0*2.46+3.0*2.42)/8.0);
  m_nfp = _md("msb",3.0);
}

void Heavy_Tensor_Current::ISGW2::CalcFFs( double q2 )
{
  double mtb = m_msb + m_msd;
  double mtx = m_msq + m_msd;
  
  double mup = 1.0/(1.0/m_msq+1.0/m_msb);
  double mum = 1.0/(1.0/m_msq-1.0/m_msb);
  double bbx2= 0.5*(m_bb2+m_bx2);
  double tm  = (m_mB-m_mT)*(m_mB-m_mT);
  if (q2>tm) q2 = 0.99*tm;
  double wt  = 1.0+(tm-q2)/(2.0*m_mbb*m_mbx);
  
  double mqm = 0.1;
  double r2  = 3.0/(4.0*m_msb*m_msq)+3*m_msd*m_msd/(2*m_mbb*m_mbx*bbx2)+
    (16.0/(m_mbb*m_mbx*(33.0-2.0*m_nfp)))*
    log(Getas(mqm)/Getas(m_msq));

  double f5  = sqrt(mtx/mtb)*pow(sqrt(m_bx2*m_bb2)/bbx2,5.0/2.0) /
       (pow((1.0+r2*(tm-q2)/18.0),3.0));
  
  double f5h = f5*pow(( m_mbb / mtb ),-1.5)*pow((m_mbx/mtx),-0.5);
  double f5k = f5*pow(( m_mbb / mtb ),-0.5)*pow((m_mbx/mtx),0.5);
  double f5bppbm = f5*pow(( m_mbb / mtb ),-2.5)*pow((m_mbx/mtx),0.5);
  double f5bpmbm = f5*pow(( m_mbb / mtb ),-1.5)*pow((m_mbx/mtx),-0.5);
  
  m_h = f5h*(m_msd/(sqrt(8.0*m_bb2)*mtb))*((1.0/m_msq)-(m_msd*m_bb2/(2.0*mum*
        mtx*bbx2)));
  
  m_k = f5k*(m_msd/(sqrt(2.0*m_bb2)))*(1.0+wt);
  
  double bppbm = ((m_msd*m_msd*f5bppbm*m_bx2)/(sqrt(32.0*m_bb2)*m_msq*m_msb*mtb*bbx2))*
          (1.0-(m_msd*m_bx2/(2.0*mtb*bbx2)));

  double bpmbm = -1.0*(m_msd*f5bpmbm/(sqrt(2.0*m_bb2)*m_msb*mtx))*(1.0-
          ((m_msd*m_msb*m_bx2)/(2.0*mup*mtb*bbx2))+((m_msd*m_bx2*(1.0-
          ((m_msd*m_bx2)/(2.0*mtb*bbx2))))/(4.0*m_msq*bbx2)));

  m_bplus = (bppbm + bpmbm)/2.0;
  m_bminus = (bppbm - bpmbm)/2.0;
  m_calced = true;
}

double Heavy_Tensor_Current::ISGW2::Getas ( double mass )
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( mass > 0.6 ) {
    if ( mass < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*M_PI / ( 33.0 - 2.0*nflav) /
      log( mass*mass/lqcd2);
  }
  return temp;
  
}



DECLARE_GETTER(Heavy_Tensor_Current_Getter, "Heavy_Tensor_Current",
               Current_Base,Flavour_Info);

Current_Base* Heavy_Tensor_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new Heavy_Tensor_Current(parameters.flavs, parameters.nout, parameters.indices, "Heavy_Tensor_Current");
}

void Heavy_Tensor_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}

