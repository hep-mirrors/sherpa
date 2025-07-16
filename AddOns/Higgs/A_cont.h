struct Value {
  Complex f, e, e2;
  inline Value(const Complex &_f,
               const Complex &_e=Complex(0.,0.),
               const Complex &_e2=Complex(0.,0.));
  inline Value operator+(const Value &v);
  inline Value operator-(const Value &v);
  inline Value Conj();
};

Complex qqbyy(int h1, int h3, int h4);
Complex qbqyy(int h1, int h3, int h4);
Value qqbyy1l(int h1, int h3, int h4);
Value qbqyy1l(int h1, int h3, int h4);
