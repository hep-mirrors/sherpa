//bof
//Version: 1 ADICIC++-0.0/2004/06/07

//Implementation of the template structures of Sudakov_Calculator.H.



//#include ""





//using;    //is already done in Sudakov_Calculator.C





//=============================================================================



template<Dipole::Type DT>
const bool Sudakov<DT>::GenerateEfracsFor(const Dipole& dip) {

  InitWith(dip);

  while(m_x2t>m_x2tmin) {

    if(GenerateX2t()==false) return false;

    GenerateRap();

    m_x3=sqrt(m_x2t);
    m_x1=1.0-m_x3*std::exp(m_rap);
    m_x3=1.0-m_x3*std::exp(-m_rap);

    m_p2t=m_s*m_x2t;

    GenerateCorr();

    if( ATOOLS::ran.Get() < m_corr ) return TestEfracs();

  }

  return false;

}



//=============================================================================



template<Dipole::Type DT>
const bool Sudakov<DT>::GenerateX2t() {

  double ran=ATOOLS::ran.Get();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"\t\t\trandom="<<ran<<endl;
#endif
  double coeff=std::log(ran)*Sudakov_Info<DT>::Colourfactor/AlphaSApprox();
  double A=sqr(std::log(m_x2t));
  if( coeff < A-sqr(std::log(m_x2tmin)) ) return false;
  m_x2t=std::exp(-sqrt(A-coeff));
  return true;

}



//=============================================================================





//eof
