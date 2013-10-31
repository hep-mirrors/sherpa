
struct Value {
  Complex f, e, e2;
  inline Value(const Complex &_f,
	       const Complex &_e=Complex(0.,0.),
	       const Complex &_e2=Complex(0.,0.)):
    f(_f), e(_e), e2(_e2) {}
  inline Value operator+(const Value &v) const { return Value(f+v.f,e+v.e,e2+v.e2); }
  inline Value operator-(const Value &v) const { return Value(f-v.f,e-v.e,e2-v.e2); }
  inline Value Conj() const { return Value(std::conj(f),std::conj(e),std::conj(e2)); }
};

Value operator*(const Complex &c,const Value &v) { return Value(c*v.f,c*v.e,c*v.e2); }

std::ostream &operator<<(std::ostream &s,const Value &v)
{ return s<<"{f="<<v.f<<"|e="<<v.e<<"|e2="<<v.e2<<"}"; }

Complex ggqqb_pmpm(int i1,int i2,int i3,int i4) { return -spa(i2,i3)*sqr(spa(i2,i4))/(spa(i3,i1)*spa(i1,i2)*spa(i4,i3)); }
Complex ggqqb_mppm(int i1,int i2,int i3,int i4) { return pow(spa(i1,i4),3)/(spa(i1,i2)*spa(i2,i4)*spa(i4,i3)); }
Complex ggqqb_mpmp(int i1,int i2,int i3,int i4) { return -spb(i2,i3)*sqr(spb(i2,i4))/(spb(i3,i1)*spb(i1,i2)*spb(i4,i3)); }
Complex ggqqb_pmmp(int i1,int i2,int i3,int i4) { return pow(spb(i1,i4),3)/(spb(i1,i2)*spb(i2,i4)*spb(i4,i3)); }

double htheta(const double &x) { return x>0.0?1.0:0.0; }

Complex lnrat(const double &x,const double &y)
{
  return Complex(log(dabs(x/y)),-M_PI*(htheta(-x)-htheta(-y)));
}

Value Fc1mp(int i1,int i2,int i3,int i4)
{
  double s12(sij(i1,i2)), s13(sij(i1,i3)), s14(sij(i1,i4));
  Complex L=lnrat(-s13,-s12), L12=lnrat(mu_sq,-s12);
  Complex e2=-C_F, e=C_F*(-3./2.-L12);
  Complex f=C_F*(-4.-.5*s12/s14*(sqr(1.-s12/s14*L)+L+sqr(s12/s14)*sqr(M_PI)));
  f+=C_F*(-.5*sqr(L12)-3./2.*L12);
  return Value(f,e,e2);
}

Value Fc1pm(int i1,int i2,int i3,int i4)
{
  double s12(sij(i1,i2)), s13(sij(i1,i3)), s14(sij(i1,i4));
  Complex L=lnrat(-s13,-s12), L12=lnrat(mu_sq,-s12);
  Complex e2=-C_F, e=C_F*(-3./2.-L12);
  Complex f=C_F*(-4.-.5*s12/s14*(sqr(L)+sqr(M_PI)));
  f+=C_F*(-.5*sqr(L12)-3./2.*L12);
  return Value(f,e,e2);
}

Complex qqbyy(int h1,int h3,int h4)
{
  if (h3==h4) return 0.0;
  if (h1<0) {
    if (h3<0) return ggqqb_mppm(3,4,2,1)+ggqqb_pmpm(4,3,2,1); 
    return ggqqb_pmpm(3,4,2,1)+ggqqb_mppm(4,3,2,1); 
  }
  else {
    if (h3>0) return std::conj(ggqqb_mppm(3,4,2,1)+ggqqb_pmpm(4,3,2,1));
    return std::conj(ggqqb_pmpm(3,4,2,1)+ggqqb_mppm(4,3,2,1));
  }
}

Value qqbyy1l(int h1,int h3,int h4)
{
  if (h3==h4) return Value(0.0);
  if (h1<0) {
    if (h3<0) return ggqqb_mppm(3,4,2,1)*Fc1mp(3,4,2,1)+ggqqb_pmpm(4,3,2,1)*Fc1pm(4,3,2,1); 
    return ggqqb_pmpm(3,4,2,1)*Fc1pm(3,4,2,1)+ggqqb_mppm(4,3,2,1)*Fc1mp(4,3,2,1); 
  }
  else {
    if (h3>0) return (ggqqb_mppm(3,4,2,1)*Fc1mp(3,4,2,1)+ggqqb_pmpm(4,3,2,1)*Fc1pm(4,3,2,1)).Conj();
    return (ggqqb_pmpm(3,4,2,1)*Fc1pm(3,4,2,1)+ggqqb_mppm(4,3,2,1)*Fc1mp(4,3,2,1)).Conj();
  }
}

Value qbqyy1l(int h2,int h3,int h4)
{
  if (h3==h4) return Value(0.0);
  if (h2<0) {
    if (h3<0) return ggqqb_mppm(3,4,1,2)*Fc1mp(3,4,1,2)+ggqqb_pmpm(4,3,1,2)*Fc1pm(4,3,1,2); 
    return ggqqb_pmpm(3,4,1,2)*Fc1pm(3,4,1,2)+ggqqb_mppm(4,3,1,2)*Fc1mp(4,3,1,2); 
  }
  else {
    if (h3>0) return (ggqqb_mppm(3,4,1,2)*Fc1mp(3,4,1,2)+ggqqb_pmpm(4,3,1,2)*Fc1pm(4,3,1,2)).Conj();
    return (ggqqb_pmpm(3,4,1,2)*Fc1pm(3,4,1,2)+ggqqb_mppm(4,3,1,2)*Fc1mp(4,3,1,2)).Conj();
  }
}
