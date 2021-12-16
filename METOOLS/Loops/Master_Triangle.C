#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! C_0(p12,p22,p32;m12,m22,m32)
/*! all momenta incoming

        p2
       -----|\              p1 + p2 + p3 = 0
            | \ m3
            |  \            s12 = p3^2 = (-p1-p2)^2 = (p1+p2)^2
            |   \   p3
         m2 |    |-----
            |   /
            |  /            s12 > 0    s-channel
            | / m1          s12 < 0    t-channel
       -----|/
        p1
*/

// Implement complex logarithm with identifier for small imaginary part
// Equivalent to C++ implementation with std::imag = +0. or -0. determining which side 
// of the branch cut x is on. This way is more clear though!
Complex METOOLS::CLog(const double& x, const int& ieps) {
  // normal real log for x > 0, but return as Complex
  if (x > 0.) { 
    return Complex(log(x),0.); 
  }
  // determine imaginary part by ieps prescription for x < 0
  else
    return Complex(log(-x),0.) + Complex(0.,ieps*M_PI);
}

Complex METOOLS::CLog(const Complex& x, const int& ieps) {
  // if argument on cut, determine Im by ieps prescription
  if (std::imag(x) == 0. && std::real(x) < 0.) { 
    return Complex(log(-std::real(x)),0.) + Complex(0.,ieps*M_PI); 
  }
  // normal log if Re > 0 or if Im != 0 - complex log picks out correct imaginary part
  else 
    return log(x);
}

Complex METOOLS::Eta(const Complex& x, const int& ix, const Complex& y, const int& iy) {
  // Veltman eta-function eta(a,b) - ix, iy are infinitesimal imaginary parts of x and y
  if (std::imag(x) != 0. && std::imag(y) != 0.) {
    return 2.*M_PI*Complex(0.,1.)*(Theta(-std::imag(x))*Theta(-std::imag(y))*Theta(std::imag(x*y))
				   -Theta(std::imag(x))*Theta(std::imag(y))*Theta(-std::imag(x*y)));
  }
  else if (std::imag(x) != 0 && std::imag(y) == 0.) {
    double imy = iy*1.;
    return 2.*M_PI*Complex(0.,1.)*(Theta(-std::imag(x))*Theta(-imy)*Theta(std::imag(x*y))
				   -Theta(std::imag(x))*Theta(imy)*Theta(-std::imag(x*y)));
  }
  else if (std::imag(x) == 0 && std::imag(y) != 0.) {
    double imx = ix*1.;
    return 2.*M_PI*Complex(0.,1.)*(Theta(-std::imag(y))*Theta(-imx)*Theta(std::imag(x*y))
				   -Theta(std::imag(y))*Theta(imx)*Theta(-std::imag(x*y)));
  }
  else if (std::imag(x) == 0 && std::imag(y) == 0.) {
    double imy = iy*1.;
    double imx = ix*1.;
    double imxy = ix*std::real(y)+iy*std::real(x);
    return 2.*M_PI*Complex(0.,1.)*(Theta(-imy)*Theta(-imx)*Theta(imxy)
				   -Theta(imy)*Theta(imx)*Theta(-imxy));
  }    
}

double METOOLS::Theta(const double& x) {
  if (x > 0.) return 1.;
  else return 0.;
}

// determine signs of real and imaginary part as ints - useful for ieps-prescription
int METOOLS::ReSign(const Complex& z) {
  if (real(z) >= 0.) return 1;
  else return -1;
}

int METOOLS::ImSign(const Complex& z) {
  if (imag(z) >= 0.) return 1;
  else return -1;
}

int METOOLS::Prod2Sign(const Complex& z1, const int& ieps1, const Complex& z2, const int& ieps2) {
  // returns small imaginary part of (z1+ieps1)*(z2+ieps2)
  double epsprod = double(ieps1)*real(z2) + double(ieps2)*real(z1);
  if (epsprod > 0.) return 1;
  else return -1;
}

int METOOLS::Prod3Sign(const Complex& z1, const int& ieps1, const Complex& z2, const int& ieps2, const Complex& z3, const int& ieps3) {
  // returns small imaginary part of (z1+ieps1)*(z2+ieps2)*(z3+ieps3)
  double epsprod = double(ieps1)*(real(z2)+real(z3)) + double(ieps2)*(real(z1)+real(z3)) + double(ieps3)*(real(z1)+real(z2));
  if (epsprod > 0.) return 1;
  else return -1;
}

int METOOLS::DivSign(const Complex& z1, const int& ieps1, const Complex& z2, const int& ieps2) {
  // returns small imaginary part of (z1+ieps1)/(z2+ieps2)
  double epsprod = double(ieps1)/real(z2) - double(ieps2)*real(z1)/sqr(real(z2));
  if (epsprod > 0.) return 1;
  else return -1;
}

int METOOLS::ProdDivSign(const Complex& z1, const int& ieps1, const Complex& z2, const int& ieps2, const Complex& z3, const int& ieps3) {
  // returns small imaginary part of (z1+ieps1)*(z2+ieps2)/(z3+ieps3)
  double epsprod = (double(ieps1)*real(z2)+double(ieps2)*real(z1))/real(z3) - double(ieps3)*real(z1*z2)/sqr(real(z3));
  if (epsprod > 0.) return 1;
  else return -1;
}

int METOOLS::DivProdSign(const Complex& z1, const int& ieps1, const Complex& z2, const int& ieps2, const Complex& z3, const int& ieps3) {
  // returns small imaginary part of (z1+ieps1)/((z2+ieps2)*(z3+ieps3))
  double epsprod = double(ieps1)/real(z2*z3) - real(z1)*(double(ieps3)*real(z2)+double(ieps2)*real(z3))/sqr(real(z2*z3));
  if (epsprod > 0.) return 1;
  else return -1;
}

void METOOLS::rfun(Complex& rr, Complex& dd, const Complex& qq) {
  // Pass rr such that qq = rr + 1/rr and Im(rr) has same sign as Im(qq)
  // If Im(qq) == 0 then Im(rr) <= 0, if Im(rr) = 0 then |rr| > 1/|rr|
  // also passes dd = rr - 1/rr
  dd = sqrt(qq*qq-4.);
  rr = qq+dd;
  Complex r2 = qq-dd;
  // Makes sure |rr| > 1/|rr|
  if (abs(r2) > abs(rr)) {
    rr = r2;
    dd = -dd;
  }
  if (imag(qq) == 0.) {
    if (imag(rr) <= 0.) {
      rr /= 2.;
    }
    else {
      rr = 2./rr;
      dd = -dd;
    }
  }
  else {
    if (ImSign(qq) == ImSign(rr)) {
      rr = rr/2.;
    }
    else {
      rr = 2./rr;
      dd = -dd;
    }
  }
}

void METOOLS::solabc(Complex& x1, Complex& x2, Complex& dd, const Complex& aa, const Complex& bb,
		     const Complex& cc) {
  // Passes solutions x1, x2 to a x^2 + b x + c = 0 as well as dd = aa*(x1-x2)
  if (aa == 0.) {
    if (bb == 0.) {
      msg_Out() << "No solution in solabc" << std::endl;
      x1 = 0.;
      x2 = 0.;
      dd = 0.;
    }
    else {
      x1 = -cc/bb;
      x2 = x1;
      dd = bb;
    }
  }
  else if (cc == 0.) {
    dd = -bb;
    x1 = dd/aa;
    x2 = 0.;
  }
  else {
    dd = sqrt(bb*bb-4.*aa*cc); // Discriminant
    if (abs(-bb+dd) > abs(-bb-dd)) {
      x1 = (-bb+dd)/(2.*aa);
      x2 = 2.*cc/(-bb+dd);
    }
    else {
      x2 = (-bb-dd)/(2.*aa);
      x1 = (2.*cc)/(-bb-dd);
    }
  }
}

Complex METOOLS::Dilog(const Complex& x, const int& ieps) {
  // Calculates the complex DiLog
  // Special values, taken from 
  // Weisstein, Eric W. "Dilogarithm." 
  // From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Dilogarithm.html
  if (x == 0.) return 0.;
  else if (x == 1.) return sqr(M_PI)/6.;
  else if (x == -1.) return -sqr(M_PI)/12.;
  else if (x == 1./2.) return sqr(M_PI)/12.-0.5*sqr(log(2.));
  // Algorithm based on "R. E. Crandall, “Note on fast polylogarithm computation,” 2006."
  // and arXiv:1601.02649v3
  double precision = 1E-20;
  double x_0(1./2.);
  Complex next(0.), sum(0.);
  
  // B_{2m}/(2m+1)! in arXiv:1601.02649v3, eq. 5.12
  const double coeffA[100] = {1., 0.0277777777777777777778, // alpha, alpha^3
			      -0.000277777777777777777778, 4.72411186696900982615e-6, // alpha^5, alpha^7
			      -9.18577307466196355085e-8, 1.8978869988970999072e-9, // ...
			      -4.06476164514422552681e-11, 8.92169102045645255522e-13,
			      -1.99392958607210756872e-14, 4.51898002961991819165e-16,
			      -1.03565176121812470145e-17, 2.39521862102618674574e-19,
			      -5.58178587432500933628e-21, 1.30915075541832128581e-22,
			      -3.08741980242674029324e-24, 7.31597565270220342036e-26,
			      -1.74084565723400074099e-27, 4.15763564461389971962e-29,
			      -9.96214848828462210319e-31, 2.39403442489616530052e-32,
			      -5.76834735536739008429e-34, 1.39317947964700797783e-35,
			      -3.37212196548508947047e-37, 8.17820877756210262176e-39,
			      -1.98701083115238592556e-40, 4.83577851804055089629e-42,
			      -1.17869372487183843267e-43, 2.877096408117257145e-45,
			      -7.03205909815602801496e-47, 1.72086031450331462909e-48,
			      -4.21607239056044549168e-50, 1.03404064051330395739e-51,
			      -2.53866306259946531616e-53, 6.23855317692459088784e-55,
			      -1.53443980691346503917e-56, 3.777294635578550234e-58,
			      -9.30586212480468658839e-60, 2.29434368222418732072e-61,
			      -5.66068873941414784857e-63, 1.39756872198540085454e-64,
			      -3.45267347330633889778e-66, 8.53498373696321710664e-68,
			      -2.11106753916379301954e-69, 5.22446788920815644493e-71,
			      -1.29363445303169121491e-72, 3.20479645174984764253e-74,
			      -7.94326694999713716406e-76, 1.96969401238294242761e-77,
			      -4.88642119356886536649e-79, 1.2127399993287380323e-80,
			      -3.01107647674682414162e-82, 7.47904589777194783774e-84,
			      -1.85837941991417435344e-85, 4.6193425842611789485e-87,
			      -1.1486235467130135708e-88, 2.8570740556428067213e-90,
			      -7.10896369091622518708e-92, 1.76940464275317093111e-93,
			      -4.40533971491237081972e-95, 1.09713120606298044082e-96,
			      -2.73313083816123243674e-98, 6.81053053662527748214e-100,
			      -1.69752549739977101562e-101, 4.23216763434037815455e-103,
			      -1.05540011102109065058e-104, 2.63254505953694582769e-106,
			      -6.56803912905139478216e-108, 1.6390562839944812655e-109,
			      -4.09116816995417822688e-111, 1.02139414028717896929e-112,
			      -2.55052340302039905225e-114, 6.3701938923623999053e-116,
			      -1.59133256352353330242e-117, 3.97605039816258639081e-119,
			      -9.93626602127902345807e-121, 2.48354935273027793979e-122,
			      -6.20866996208690084797e-124, 1.55238189965646443514e-125,
			      -3.88213719560680557378e-127, 9.70987570883174305432e-129,
			      -2.42898695457098901872e-130, 6.07720263182787102327e-132,
			      -1.52071433809476079929e-133, 3.8058825087439879074e-135,
			      -9.52632528675517774988e-137, 2.38482361980845745438e-138,
			      -5.97099262734982093322e-140, 1.49518472872052581311e-141,
			      -3.74455227876193378591e-143, 9.37908338404674472152e-145,
			      -2.34949819836149205611e-146, 5.88630640090960269216e-148,
			      -1.47489970707240971116e-149, 3.696007761196264006e-151,
			      -9.26302721807284788094e-153, 2.32178307158320397184e-154,
			      -5.82020071424204647756e-156, 1.45915330388174209927e-157,
			      -3.65855485706352121659e-159, 9.17408974648312441359e-161};

  // B_{2m}/(2m(2m+1)!) in arXiv:1601.02649v3, eq. 5.9
  const double coeffB[100] = {-1., 0.013888888888888888889, // alpha, alpha^3
			      -0.000069444444444444444444, 7.8735197782816830436e-7, // alpha^5, alpha^7
			      -1.1482216343327454439e-8, 1.8978869988970999072e-10, // ...
			      -3.3873013709535212723e-12, 6.3726364431831803966e-14,
			      -1.2462059912950672305e-15, 2.5105444608999545509e-17,
			      -5.1782588060906235072e-19, 1.0887357368300848844e-20,
			      -2.3257441143020872235e-22, 5.0351952131473895608e-24,
			      -1.1026499294381215333e-25, 2.4386585509007344735e-27,
			      -5.4401426788562523156e-29, 1.2228340131217352117e-30,
			      -2.7672634689679505842e-32, 6.3000905918320139487e-34,
			      -1.4420868388418475211e-35, 3.3170939991595428044e-37,
			      -7.6639135579206578874e-39, 1.7778714733830657873e-40,
			      -4.1396058982341373449e-42, 9.6715570360811017926e-44,
			      -2.2667187016766123705e-45, 5.3279563113282539722e-47,
			      -1.2557248389564335741e-48, 2.9670005422470941881e-50,
			      -7.0267873176007424861e-52, 1.6678074846988773506e-53,
			      -3.9666610353116645565e-55, 9.4523532983705922543e-57,
			      -2.2565291278139191752e-58, 5.3961351936836431914e-60,
			      -1.2924808506673175817e-61, 3.1004644354380909739e-63,
			      -7.4482746571238787481e-65, 1.7917547717761549417e-66,
			      -4.3158418416329236222e-68, 1.0408516752394167203e-69,
			      -2.5131756418616583566e-71, 6.0749626618699493546e-73,
			      -1.4700391511723763806e-74, 3.5608849463887196028e-76,
			      -8.6339858152142795261e-78, 2.0954191621095132209e-79,
			      -5.0900220766342347568e-81, 1.2374897952334061554e-82,
			      -3.0110764767468241416e-84, 7.3323979389921057233e-86,
			      -1.7869032883790138014e-87, 4.3578703625105461778e-89,
			      -1.0635403210305681211e-90, 2.5973400505843697466e-92,
			      -6.3472890097466296313e-94, 1.5521093357483955536e-95,
			      -3.7977066507865265687e-97, 9.2977220852794952612e-99,
			      -2.2776090318010270306e-100, 5.5824020792010471165e-102,
			      -1.3689721753223959803e-103, 3.3588632018574429798e-105,
			      -8.2453133673522707077e-107, 2.0250346611822660213e-108,
			      -4.9757872189783293804e-110, 1.2231763313391651235e-111,
			      -3.0082118896721898727e-113, 7.4014068136752099224e-115,
			      -1.8218024307288564659e-116, 4.4860520368749295108e-118,
			      -1.1050920580024536822e-119, 2.7233221905223194458e-121,
			      -6.713693257620961796e-123, 1.6556995684868519599e-124,
			      -4.0846512908466452947e-126, 1.0080401945821197631e-127,
			      -2.4885494843633369063e-129, 6.1454909549567994015e-131,
			      -1.5181168466068681367e-132, 3.7513596492764635946e-134,
			      -9.272648403016834142e-136, 2.2927003064722818719e-137,
			      -5.6704317183066534225e-139, 1.402837423416739679e-140,
			      -3.4715073414824540309e-142, 8.5930156823018724891e-144,
			      -2.1275865220238260147e-145, 5.2691479685655869222e-147,
			      -1.3052767768674955867e-148, 3.2342342862140674133e-150,
			      -8.0157592775674440824e-152, 1.9871009468797118312e-153,
			      -4.927142137272791426e-155, 1.2219910903069494589e-156,
			      -3.0313545386677325404e-158, 7.5214087828955778313e-160,
			      -1.8666096209507761309e-161, 4.633378659839961825e-163};
  
  if (std::real(x) < x_0 && std::abs(x) <= 1.) {
    Complex alpha = -CLog(1.-x,-ieps);
    int i(1);
    Complex prod(alpha*alpha*alpha);
    sum += alpha*coeffA[0]; // alpha term
    next = prod*coeffA[1]; // alpha^3 term
    // only add next term if precision not reached yet
    while (std::abs(next/sum) > precision && (i < 100)) {
      ++i;
      sum += next;
      prod *= alpha*alpha; 
      next = prod*coeffA[i]; // alpha^(2i+1) term
    }
    // add in alpha^2 term
    return sum - 0.25*alpha*alpha;
  }
  else if (std::real(x) > x_0 && std::abs(x-1.) <= 1.) {
    Complex alpha = -CLog(x,ieps);
    int aieps = -ieps*ReSign(x); // required if Re(x) > 0. and Im(x) == 0.
    int i(1);
    Complex prod(alpha*alpha*alpha);
    sum += alpha*coeffB[0]; // alpha term
    next = prod*coeffB[1]; // alpha term
    // only add next term if precision not reached yet
    while (std::abs(next/sum) > precision && (i < 100)) {
      ++i;
      sum += next;
      prod *= alpha*alpha; 
      next = prod*coeffB[i]; // alpha^(2i+1) term
    }
    // add in additional terms
    return sqr(M_PI)/6. - 0.25*alpha*alpha + alpha*CLog(alpha,aieps) + sum;
  }
  else if (std::abs(x) > 1.) {
    // use inversion relation to bring into convergent region
    return - Dilog(1./x,-ieps) - 0.5*sqr(CLog(-x,-ieps)) - sqr(M_PI)/6.;
  }
  else {
    //msg_Out() << "No convergent region found for Dilog. x = " << x;
  }
}






METOOLS::DivArrC
METOOLS::Master_Triangle(const double&  p12, const double&  p22, const double&  p32,
                         const Complex& m12, const Complex& m22, const Complex& m32,
                         double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;

  /// Adopt FF/QCDLoop approach: move masses to the end using symmetry relations
  /// Note: Bardin/Passarino is out by a factor (-1)^n due to choice of metric

  // massless internal lines
  if (IsZero(m12) && IsZero(m22) && IsZero(m32)) {
    // All zero
    if (IsZero(p12) && IsZero(p22) && IsZero(p32)) {
      // C_0(0,0,0;0,0,0) = ?
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    // one momentum squared non-zero - different cases related by symmetry (eq. 4.5)
    // analytic continuation by pi2 -> pi2 + ieps
    else if ((IsZero(p12) && IsZero(p22)) ||
	     (IsZero(p12) && IsZero(p32)) ||
	     (IsZero(p22) && IsZero(p32))) {
      double p2 = p12 + p22 + p32; 
      //! C_0(0,0,s;0,0,0) = 1/s*(1/epsIR2 + 1/epsIR*ln(mu2/-s) + 1/2*ln2(mu2/-s))         s < 0
      return 1./p2*DivArrC(0.,CLog(-mu2/p2,1),1.,0.5*sqr(CLog(-mu2/p2,1)),0.,0.);
    }
    
    // two momenta squared non-zero (eq. 4.6, 4.7)
    else if (IsZero(p12)) {
      if (IsZero(p22-p32)) {
	// C_0(0,p2,p2,0,0,0) = -1/p2*(1/epsIR+log(mu2/-p2))
	return -1./p22*DivArrC(0.,1.,0.,CLog(-mu2/p22,1),0.,0.);
      }
      else 
	// C_0(0,p22,p32,0,0,0) = 1/(p22-p32)*(1/epsIR*log(p32/p22) + 0.5*(log^2(mu2/-p22)-log^2(mu2/-p32))) 
	return 1./(p22-p32)*DivArrC(0.,CLog(-mu2/p22,1)-CLog(-mu2/p32,1),
				    0.,0.5*(sqr(CLog(-mu2/p22,1))
					    -sqr(CLog(-mu2/p32,1))),0.,0.);
      
    }
    // C_0(p12,0,p32,0,0,0) = C_0(0,p32,p12,0,0,0)
    else if (IsZero(p22)) {
      return Master_Triangle(0.,p32,p12,0.,0.,0.,mu2);
    }
    // C_0(p12,p22,0,0,0,0) = C_0(0,p12,p22,0,0,0)
    else if (IsZero(p32)) {
      return Master_Triangle(0.,p12,p22,0.,0.,0.,mu2);
    }
  
    // all momenta squared non zero (Bardin,Passarino eq. 5.70)
    // C_0(p12,p22,p32,0,0,0) 
    else {
      Complex ap = 0.5/p32*(p32+p12-p22-csqrt(sqr(p32)+sqr(p22)+sqr(p12)-2.*p12*p22-2.*p12*p32-2.*p22*p32));
      Complex am = 0.5/p32*(p32+p12-p22+csqrt(sqr(p32)+sqr(p22)+sqr(p12)-2.*p12*p22-2.*p12*p32-2.*p22*p32));
      return 1./(p32*(ap-am))*DivArrC(0.,0.,0.,log(ap*am)*log((1.-ap)/(1.-am))+2.*Dilog(ap)-2.*Dilog(am),0.,0.);
    }
  }
  
  // two massless internal lines - rotate all such that non-zero m2 is at end
  if (IsZero(m12) && IsZero(m32)) {
    return Master_Triangle(p32,p12,p22,0.,0.,m22,mu2);
  }
  if (IsZero(m22) && IsZero(m32)) {
    return Master_Triangle(p22,p32,p12,0.,0.,m12,mu2);
  }
  if (IsZero(m12) && IsZero(m22)) {
    Complex m2 = m12 + m22 + m32;
    if (IsZero(p12) && IsZero(p22) && IsZero(p32)) {
      // C_0(0,0,0;0,0,m2) = 1/m2*(1/epsIR + log(mu2/m2)+1)  from QCDLoop general exp. (m2 -> m2-ieps)
      return 1./m2*DivArrC(0.,1.,0.,CLog(mu2/m2,1)+1.,0.,0.);
    }
    // one momentum squared non-zero (QCDLoop eq. 4.8, 4.9)
    // continuation by pi2 -> pi2 + ieps
    else if ((IsZero(p12) && IsZero(p22)) ||
	     (IsZero(p12) && IsZero(p32))) {
      double p2 = p22 + p32;
      if (IsEqual(p2,real(m2))) {
	// C_0(0,0,m2;0,0,m2)
	return -1./m2*DivArrC(0.,0.5*CLog(mu2/m2,-1),0.5,sqr(M_PI)/12.+0.25*sqr(CLog(mu2/m2,-1)),0.,0.);
      }
      // C_0(0,0,p32;0,0,m2)
      else {
	return -1./p2*DivArrC(0.,CLog((m2-p2)/m2,-1),0.,CLog(mu2/m2,1)*CLog((m2-p2)/m2,-1)-Dilog(p2/m2,1)-sqr(CLog((m2-p2)/m2,-1)),0.,0.);
      }
    }
    else if (IsZero(p22) && IsZero(p32)) { 
      // Bardin/Passarino eq. 5.59
      // C_0(p12,0,0;0,0,m2)
      return 1./p12*DivArrC(0.,0.,0.,Dilog(1.)-Dilog(1.+p12/m2,1),0.,0.); 
    }    
    else if (IsZero(p12)) {
      if (IsEqual(p22,p32)) {
	// QCDLoop eq. 4.12
	if (IsEqual(p22,m2)) {
	  // C_0(0,m2,m2;0,0,m2)
	  return 1./m2*DivArrC(0.,-0.5,0.,-0.5*CLog(mu2/m2,1)+1.,0.,0.);
	}
	// QCDLoop eq. 4.9
	// C_0(0,p22,p22;0,0,m2)
	else 
	  return 1./(m2-p22)*DivArrC(0.,1.,0.,CLog(mu2/m2,1)+(m2+p22)/p22*CLog(m2/(m2-p22),1),0.,0.);
      }
      else if (!IsEqual(p22,p32)) {
	// QCDLoop eq. 4.11
	if (IsEqual(p32,m2)) {
	  // C_0(0,p22,m2;0,0,m2)
	  return 1./(p22-m2)*DivArrC(0.,0.5*CLog(mu2/m2,1)+CLog(m2/(m2-p22),1),0.5,sqr(M_PI)/12.+0.25*sqr(CLog(mu2/m2,1))+0.5*sqr(CLog(m2/(m2-p22),1))-Dilog(-p22/(m2-p22),-1),0.,0.);
	}
	else if (IsEqual(p22,m2)) {
	  // C_0(0,m2,p32;0,0,m2) = C_0(0,p32,m2,0,0,m2)
	  return Master_Triangle(0.,p32,p22,0.,0.,m2,mu2);
	}
	// QCDLoop eq. 4.8
	// C_0(0,p22,p32;0,0,m2)
	else {
	  return 1./(p22-p32)*DivArrC(0.,CLog((m2-p32)/(m2-p22),ReSign(p22-p32)),0.,Dilog(p22/m2,1)+sqr(CLog((m2-p22)/m2,-1))-Dilog(p32/m2,1)-sqr(CLog((m2-p32)/m2,-1))+CLog(mu2/m2,1)*CLog((m2-p32)/(m2-p22),ReSign(p22-p32)),0.,0.);
	}
      }
    }
    else if ((IsZero(p22) && !IsZero(p12) && !IsZero(p32)) ||
	     (IsZero(p32) && !IsZero(p12) && !IsZero(p22))) {
      // C_0(p12,p22,0,0,0,m2) = ?
      return DivArrC(0.,0.,0.,0.,0.,0.); 
    }	 
    // QCDLoop online, finite triangle 1.
    else if (IsEqual(p22,m2) && IsEqual(p32,m2)) {
      // C_0(p12,m2,m2;0,0,m2)
      Complex b = sqrt(1. - 4.*m2/p12);
      Complex lp = 0.5*(1.+b);
      Complex lm = 0.5*(1.-b);
      return 1./(p12*b)*DivArrC(0.,0.,0.,2./3.*sqr(M_PI) + 2.*Dilog(-lm/lp,ReSign(p12-4.*m2)) + 0.5*sqr(CLog(-lm/lp,ReSign(p12-4.*m2))),0.,0.);
    }
    else {
      // Finite triangle with one non-zero mass. Following Denner:1991kt
      // C_0(p12,p22,p32,0,0,m32)
      Complex res(0.,0.);
      Complex alpha = sqrt(sqr(p12)+sqr(p32)+sqr(p22)-2.*p12*p32-2.*p12*p22-2.*p32*p22);
      Complex alphai[] = {sqrt(sqr(p22)+sqr(m32)-2.*p22*m32),
			  sqrt(sqr(p32)+sqr(m32)-2.*p32*m32),
			  p12};
      // small imaginary sign of alphai: if alphai is real, ieps is sign(alphai*pjk^2)
      //                                 if alphai is imaginary, ieps is insignificant - choose 1
      double pjk2[] = {p22,p32,p12};
      int sgnalphai[] = {0,0,0};
      for (int i(0); i < 3; ++i) {
	if (std::imag(alphai[i]) > 0.) sgnalphai[i] = 1;
	else {
	  if (std::real(alphai[i]*pjk2[i]) > 0.) sgnalphai[i] = 1;
	  else sgnalphai[i] = -1;
	}
      }
      Complex xip[] = {0.5/p22*(p22+m32+alphai[0]),
		       0.5/p32*(p32-m32+alphai[1]),
		       0.5/p12*(p12+alphai[2])};
      Complex xim[] = {0.5/p22*(p22+m32-alphai[0]),
		       0.5/p32*(p32-m32-alphai[1]),
		       0.5/p12*(p12-alphai[2])};
      Complex y0i[] = {0.5/(alpha*p22)*(p22*(p22-p32-p12-m32)
					-(p32-p12)*(-m32)+alpha*(p22+m32)),
		       0.5/(alpha*p32)*(p32*(p32-p12-p22-m32)
					-(p12-p22)*m32+alpha*(p32-m32)),
		       0.5/(alpha*p12)*(p12*(p12-p22-p32+2.*m32)
					+alpha*p12)};
      Complex yip[] = {y0i[0]-xip[0],
		       y0i[1]-xip[1],
		       y0i[2]-xip[2]};
      Complex yim[] = {y0i[0]-xim[0],
		       y0i[1]-xim[1],
		       y0i[2]-xim[2]};
      if (std::imag(alpha) == 0 && std::imag(m12) == 0 && std::imag(m22) == 0 && std::imag(m32) == 0){
	for (int i = 0; i < 3; ++i) {
	  res += Dilog((y0i[i]-1.)/yim[i],-ReSign((y0i[i]-1.))*sgnalphai[i])
	    - Dilog(y0i[i]/yim[i],-ReSign(y0i[i])*sgnalphai[i]) 
	    + Dilog((y0i[i]-1.)/yip[i],ReSign((y0i[i]-1.))*sgnalphai[i]) 
	    - Dilog(y0i[i]/yip[i],ReSign(y0i[i])*sgnalphai[i]);
	}
	res /= alpha;
	return res*DivArrC(0.,0.,0.,1.,0.,0.);  
      }
      else {
	for (int i = 0; i < 3; ++i) {
	  // Figure out theta function in final term
	  // First Im(yip*yim) part
	  int theta = 1;
	  Complex yipyim = yip[i]*yim[i];
	  if ((std::imag(yipyim) == 0. && std::real(alphai[i]) > 0.) ||
	      (std::imag(yipyim) > 0.)) theta*=0;
	  // Second pjk2 part
	  if (pjk2[i] > 0.) theta*=0;
	  res += Dilog((y0i[i]-1.)/yim[i],-ReSign((y0i[i]-1.))*sgnalphai[i]) 
	    - Dilog(y0i[i]/yim[i],-ReSign(y0i[i])*sgnalphai[i]) 
	    + Dilog((y0i[i]-1.)/yip[i],ReSign((y0i[i]-1.))*sgnalphai[i]) 
	    - Dilog(y0i[i]/yip[i],ReSign(y0i[i])*sgnalphai[i])
	    + Eta(1.-xim[i],sgnalphai[i],1./yim[i],-sgnalphai[i])*CLog((y0i[i]-1.)/yim[i],-ReSign((y0i[i]-1.))*sgnalphai[i])
	    - Eta(-xim[i],sgnalphai[i],1./yim[i],-sgnalphai[i])*CLog(y0i[i]/yim[i],-ReSign(y0i[i])*sgnalphai[i])
	    + Eta(1.-xip[i],-sgnalphai[i],1./yip[i],sgnalphai[i])*CLog((y0i[i]-1.)/yip[i],ReSign((y0i[i]-1.))*sgnalphai[i])
	    - Eta(-xip[i],-sgnalphai[i],1./yip[i],sgnalphai[i])*CLog(y0i[i]/yip[i],ReSign(y0i[i])*sgnalphai[i])
	    -(Eta(-xip[i],-sgnalphai[i],-xim[i],sgnalphai[i])
	      -Eta(yip[i],-sgnalphai[i],yim[i],sgnalphai[i])
	      -theta*2.*M_PI*Complex(0.,1.))*CLog(-(1.-y0i[i])/y0i[i]);
	}
	res /= alpha;
	return res*DivArrC(0.,0.,0.,1.,0.,0.);
      }
    }
  }    
  if (IsZero(m12) && IsZero(m32)) {
    // This shouldn't happen - masses should be rotated already. But keep for sanity.
    return Master_Triangle(p32,p12,p22,0.,0.,m22,mu2);
  }
  if (IsZero(m22) && IsZero(m32)) {
    // This shouldn't happen - masses should be rotated already. But keep for sanity.
    return Master_Triangle(p22,p32,p12,0.,0.,m12,mu2);
  }
  

  // one massless internal line
  if (IsZero(m22)) {
    return Master_Triangle(p22,p32,p12,0.,m32,m12,mu2);
  }
  if (IsZero(m32)) {
    return Master_Triangle(p32,p12,p22,0.,m12,m22,mu2);
  }
  if (IsZero(m12)) {
    if (IsEqual(p12,m22) && IsEqual(p32,m32)) {
      Complex sm22 = sqrt(m22);
      Complex sm32 = sqrt(m32);
      // QCDLoop 4.16
      // C_0(m22,s,m32,0,m22,m32)
      if (IsEqual(p22,sqr(sm22-sm32))) {
	return 0.5/(sm22*sm32)*DivArrC(0.,1.,0.,CLog(mu2/(sm22*sm32),1)-2.
				       -(sm22+sm32)/(sm32-sm22)*CLog(sm22/sm32,ReSign(sm22-sm32)),
				       0.,0.);
      }
      // QCDLoop 4.14
      else {
	Complex xs = -(1.-sqrt(1.-4.*sm22*sm32/(p22-sqr(sm22-sm32))))/(1.+sqrt(1.-4.*sm22*sm32/(p22-sqr(sm22-sm32))));
	return xs/(sm22*sm32*(1.-sqr(xs)))*DivArrC(0.,-CLog(xs,1),0.,
						   -0.5*sqr(CLog(xs,1))
						   +2.*CLog(xs,1)*CLog(1.-sqr(xs),-1)
						   +CLog(xs,1)*CLog(sm22*sm32/mu2,-1)-sqr(M_PI)/6.
						   +Dilog(sqr(xs),1)
						   +0.5*sqr(CLog(sm22/sm32,DivSign(sm22,-1,sm32,-1)))
						   +Dilog(1.-xs*sm22/sm32,-ProdDivSign(xs,1,sm22,-1,sm32,-1))
						   +Dilog(1.-xs*sm32/sm22,-ProdDivSign(xs,1,sm32,-1,sm22,-1)),
						   0.,0.);
      }
    }
    // else finite ...
    
    // Bardin/Passarino eq. 5.67
    else if (IsEqual(m22,m32) && IsZero(p12) && IsZero(p32)) {
      // C_0(0,p22,0,0,m2,m2)
      Complex m2 = 0.5*(m22+m32);
      Complex b = sqrt(1.-4.*m2/p22);
      return 1./p22*DivArrC(0.,0.,0.,sqr(CLog((1.+b)/(b-1.),-ReSign(p22-4.*m2))),0.,0.);
    }
    // Bardin/Passarino eq. 5.66
    else if (IsZero(p12) && IsZero(p32)) {
      // C_0(0,p22,0,0,m22,m32)
      Complex x1 = 0.5/p22*(p22+m32-m22+sqrt(sqr(p22)+sqr(m32)+sqr(m22)-2.*p22*m22-2.*p22*m32-2.*m22*m32));
      Complex x2 = 0.5/p22*(p22+m32-m22-sqrt(sqr(p22)+sqr(m32)+sqr(m22)-2.*p22*m22-2.*p22*m32-2.*m22*m32));
      return 1./p22*DivArrC(0.,0.,0.,log(x2/(x2-1.))*log((x1-1.)/x1),0.,0);
    }
    // Bardin/Passarino eq. 5.68
    else if (IsZero(p12) && IsZero(p22)) {
      // C_0(0,0,p32,0,m22,m32)
      return 1./p32*DivArrC(0.,0.,0.,Dilog(1.-m32/m22,ReSign(m22-m32))-Dilog(1.-(m32-p32)/m22,1),0.,0.);
    }
    // QCDLoop online, finite triangle 2.
    else if (IsEqual(m22,m32) && !IsEqual(p12,m22) && IsZero(p22) && IsEqual(p32,m32)) { 
      // C_0(0,m2,p32,m2,m2,0)
      Complex m2 = 0.5*(m22 + m32);
      return -1./(m2-p12)*DivArrC(0.,0.,0.,sqr(M_PI)/6. - Dilog(p12/m2,1),0.,0.);
    }
    else if (IsEqual(m22,m32) && IsEqual(p12,m22) && IsZero(p22) && !IsEqual(p32,m32)) { 
      // C_0(0,m2,p32,m2,m2,0)
      Complex m2 = 0.5*(m22 + m32);
      return -1./(m2-p32)*DivArrC(0.,0.,0.,sqr(M_PI)/6. - Dilog(p32/m2,1),0.,0.);
    }
    else if (IsEqual(m22,m32) && IsEqual(p12,p22) && IsEqual(p12,m22) && IsZero(p32)) {
      // C_0(m2,m2,0,0,m2,m2)
      Complex m2 = 0.5*(m22+m32);
      return 0.5/m2*DivArrC(0.,1.,0.,CLog(mu2/m2,1),0.,0.);
    }
    else if (IsEqual(m22,m32) && IsEqual(p12,p22) && IsZero(p32) && !IsEqual(m22,p22)) {
      // Own calculation.
      double p2 = 0.5*(p12+p22);
      Complex m2 = 0.5*(m22+m32);
      Complex lambda = 1. + m2/p2*m2/p2-2.*m2/p2;
      Complex L1 = CLog(1.-(1.+m2/p2)/lambda,1);
      Complex L2 = CLog(1.+(1.+m2/p2)/lambda,1);
      Complex L3 = CLog(1.-(1.-m2/p2)/lambda,1);
      Complex L4 = CLog(1.+(1.-m2/p2)/lambda,1);
      return 1./(p2*lambda)*DivArrC(0.,0.,0.,(L1-L2+L3-L4),0.,0.);
    }
    else {
      // C_0(p12,p22,p32,0,m22,m32) following One-LOop implementation by Van Hameren 2010,
      // based on Denner, Nerste, Scharf 1991 D0 with one massless propagator and one mass (m1) taken
      // to infinity
      Complex res(0.,0.);
      Complex sm2 = sqrt(m32);
      Complex sm3 = abs(sm2);
      Complex sm4 = sqrt(m22);
      Complex r23 = (m32 - p32)/(sm2*sm3);
      int epsr23 = DivProdSign(m32-p32,-1,sm2,-1,sm3,1);
      Complex k24 = (m32+m22-p22)/(sm2*sm4);
      int epsk24 = DivProdSign(m32+m22-p22,-1,sm2,-1,sm4,-1);
      Complex r34 = (m22-p12)/(sm3*sm4);
      int epsr34 = DivProdSign(m22-p12,-1,sm3,1,sm4,-1);
      // Set r24, d24
      Complex r24, d24;
      rfun(r24,d24,k24);
      Complex aa = r34/r24-r23;
      
      if (aa == 0.) {
	msg_Out() << "Threshold singularity, returning 0" << std::endl;
	return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      
      Complex bb = -d24/sm3 + r34/sm2 - r23/sm4;
      Complex cc = (sm4/sm2-r24)/(sm3*sm4);
      Complex x1, x2, dd;
      solabc(x1,x2,dd,aa,bb,cc);
      x1 = -x1;
      x2 = -x2;
      Complex y1 = x1/r24;
      int epsy1 = DivSign(x1,1,r24,-1);
      Complex y2 = x2/r24;
      int epsy2 = DivSign(x2,1,r24,-1);
      
      res += (Dilog(y1*sm2,Prod2Sign(y1,epsy1,sm2,-1))-Dilog(y2*sm2,Prod2Sign(y2,epsy2,sm2,-1)))/((y1*sm2-y2*sm2)*r24*sm2);
      
      if (x2 != 0.) {
	res += (CLog(y1/y2,DivSign(y1,epsy1,y2,epsy2))/(1.-y1/y2)
		*CLog(y1*y2/(sm2*sm2),ProdDivSign(y1,epsy1,y2,epsy2,sm2*sm2,-1))
		-CLog(x1/x2,DivSign(x1,1,x2,1))/(1.-x1/x2)
		*CLog(x1*x2/(sm4*sm4),ProdDivSign(x1,1,x2,1,sm4*sm4,-1)))/(2.*x2);
      }

      res += - sm4*(Dilog(x1*sm4,Prod2Sign(x1,1,sm4,-1))-Dilog(x2*sm4,Prod2Sign(x2,1,sm4,-1)))/(x1*sm4-x2*sm4);
      
      if (abs(real(r23))+abs(imag(r23)) != 0.) {
	Complex qss = r23*sm3/r24;
	int epsqss = ProdDivSign(r23,-1,sm3,-1,r24,-1);
	
	res += -r23*sm3/r24*(Dilog(x1*qss,Prod2Sign(x1,1,qss,epsqss))-Dilog(x2*qss,Prod2Sign(x2,1,qss,epsqss)))/(x1*qss-x2*qss);
      }

      if (abs(real(r34))+abs(imag(r34)) != 0.) {
	Complex qss = r34*sm3;
	int epsqss = Prod2Sign(r34,-1,sm3,-1);
	
	res += r34*sm3*(Dilog(x1*qss,Prod2Sign(x1,1,qss,epsqss))-Dilog(x2*qss,Prod2Sign(x2,1,qss,epsqss)))/(x1*qss-x2*qss);
      }

      res /= (aa*sm2*sm3*sm4);
      return res*DivArrC(0.,0.,0.,1.,0.,0.);
    }
  }
  if (IsZero(m22)) {
    // This shouldn't happen - masses should be rotated already. But keep for sanity
    return Master_Triangle(p22,p32,p12,0.,m32,m12,mu2);
  }
  if (IsZero(m32)) {
    // This shouldn't happen - masses should be rotated already. But keep for sanity
    return Master_Triangle(p32,p12,p22,0.,m12,m22,mu2);
  }

  if (IsEqual(m12,m22) && IsEqual(m12,m32)) {
    Complex m2 = 0.5*(m12+m22);
    // QCDLoop online, finite triangles 3.
    if ((IsZero(p12) && IsZero(p22)) ||
	(IsZero(p12) && IsZero(p32)) ||
	(IsZero(p22) && IsZero(p32))) {
      // C_0(0,0,p2,m2,m2,m2)
      double p2 = p12 + p22 + p32;
      Complex b = sqrt(1.-4.*m2/p2);
      Complex lp = 0.5*(1.+b);
      Complex lm = 0.5*(1.-b);
      return 0.5/p2*DivArrC(0.,0.,0.,sqr(CLog(-lm/lp,ReSign(p2-4.*m2))),0.,0.);
    }
    // QCDLoop online, finite triangles 4.
    else if (IsZero(p12) && !IsEqual(p22,p32)) {
      // C_0(0,p22,p32,m2,m2,m2)
      return 1./(p22-p32)*(p22*Master_Triangle(0.,0.,p22,m2,m2,m2,mu2) - p32*Master_Triangle(0.,0.,p32,m2,m2,m2,mu2));
    }
    else if (IsZero(p22) && !IsEqual(p12,p32)) {
      return 1./(p32-p12)*(p32*Master_Triangle(0.,0.,p32,m2,m2,m2,mu2) - p12*Master_Triangle(0.,0.,p12,m2,m2,m2,mu2));
    }
    else if (IsZero(p32) && !IsEqual(p12,p22)) {
      return 1./(p12-p22)*(p12*Master_Triangle(0.,0.,p12,m2,m2,m2,mu2) - p22*Master_Triangle(0.,0.,p22,m2,m2,m2,mu2));
    }
  }
  // Bardin/Passarino, eq. 5.64   
  else if (IsZero(p12) && IsZero(p22)) {
    // C_0(0,0,p32,m12,m22,m32)
    Complex x0 = 1.-(m12-m22)/p32;
    Complex x1 = 0.5/p32*(p32+m12-m32+sqrt(sqr(p32)+sqr(m12)+sqr(m32)-2.*p32*m12-2.*p32*m32-2.*m12*m32));
    Complex x2 = 0.5/p32*(p32+m12-m32-sqrt(sqr(p32)+sqr(m12)+sqr(m32)-2.*p32*m12-2.*p32*m32-2.*m12*m32));
    Complex x3 = m32/(m32-m22);
    return 1./p32*DivArrC(0.,0.,0.,
			  Dilog((x0-1.)/(x0-x1))-Dilog(x0/(x0-x1))
			  +Dilog((x0-1.)/(x0-x2))-Dilog(x0/(x0-x2))
			  -Dilog((x0-1.)/(x0-x3))+Dilog(x0/(x0-x3)),0.,0.);
  }
  

  // full finite result - Denner eq. 4.26/4.27, for real masses
  if ((IsZero(m12) && IsZero(m22) && !IsEqual(p22,m32)) ||
	   (IsZero(m12) && (!IsEqual(p12,m22) || !IsEqual(p32,m32))) ||
	   (!IsZero(p12) && !IsZero(p22) && !IsZero(p32) &&
	    !IsZero(m12) && !IsZero(m22) && !IsZero(m32))){
    // C_0(p12,p22,p32,m12,m22,m32)
    Complex res(0.,0.);
    Complex alpha = sqrt(sqr(p12)+sqr(p32)+sqr(p22)-2.*p12*p32-2.*p12*p22-2.*p32*p22);
    Complex alphai[] = {sqrt(sqr(p22)+sqr(m22)+sqr(m32)-2.*p22*m22-2.*p22*m32-2.*m22*m32),
    			sqrt(sqr(p32)+sqr(m32)+sqr(m12)-2.*p32*m32-2.*p32*m12-2.*m32*m12),
    			sqrt(sqr(p12)+sqr(m12)+sqr(m22)-2.*p12*m12-2.*p12*m22-2.*m12*m22)};
    // small imaginary sign of alphai: if alphai is real, ieps is sign(alphai*pjk^2)
    //                                 if alphai is imaginary, ieps is insignificant - choose 1
    double pjk2[] = {p22,p32,p12};
    int sgnalphai[] = {0,0,0};
    for (int i(0); i < 3; ++i) {
      if (std::imag(alphai[i]) > 0.) sgnalphai[i] = 1;
      else {
    	if (std::real(alphai[i]*pjk2[i]) > 0.) sgnalphai[i] = 1;
    	else sgnalphai[i] = -1;
      }
    }
    Complex xip[] = {0.5/p22*(p22-m22+m32+alphai[0]),
    		     0.5/p32*(p32-m32+m12+alphai[1]),
    		     0.5/p12*(p12-m12+m22+alphai[2])};
    Complex xim[] = {0.5/p22*(p22-m22+m32-alphai[0]),
    		     0.5/p32*(p32-m32+m12-alphai[1]),
    		     0.5/p12*(p12-m12+m22-alphai[2])};
    Complex y0i[] = {0.5/(alpha*p22)*(p22*(p22-p32-p12+2.*m12-m22-m32)
    				      -(p32-p12)*(m22-m32)+alpha*(p22-m22+m32)),
    		     0.5/(alpha*p32)*(p32*(p32-p12-p22+2.*m22-m32-m12)
    				      -(p12-p22)*(m32-m12)+alpha*(p32-m32+m12)),
    		     0.5/(alpha*p12)*(p12*(p12-p22-p32+2.*m32-m12-m22)
    				      -(p22-p32)*(m12-m22)+alpha*(p12-m12+m22))};
    // yip take imaginary sign of -alphai
    // yim take imaginary sign of +alphai
    Complex yip[] = {y0i[0]-xip[0],
		     y0i[1]-xip[1],
		     y0i[2]-xip[2]};
    Complex yim[] = {y0i[0]-xim[0],
		     y0i[1]-xim[1],
		     y0i[2]-xim[2]};
    if (std::imag(alpha) == 0 && std::imag(m12) == 0 && std::imag(m22) == 0 && std::imag(m32) == 0){
      for (int i = 0; i < 3; ++i) {
	res += Dilog((y0i[i]-1.)/yim[i],-ReSign((y0i[i]-1.))*sgnalphai[i]) 
	  - Dilog(y0i[i]/yim[i],-ReSign(y0i[i])*sgnalphai[i]) 
	  + Dilog((y0i[i]-1.)/yip[i],ReSign((y0i[i]-1.))*sgnalphai[i]) 
	  - Dilog(y0i[i]/yip[i],ReSign(y0i[i])*sgnalphai[i]);
      }
      res /= alpha;
      return res*DivArrC(0.,0.,0.,1.,0.,0.);
    }
    else {
      for (int i = 0; i < 3; ++i) {
	// Figure out theta function in final term
	// First Im(yip*yim) part
	int theta = 1;
	Complex yipyim = yip[i]*yim[i];
	if ((std::imag(yipyim) == 0. && std::real(alphai[i]) > 0.) ||
	    (std::imag(yipyim) > 0.)) theta*=0;
	// Second pjk2 part
	if (pjk2[i] > 0.) theta*=0;
	res += Dilog((y0i[i]-1.)/yim[i],-ReSign((y0i[i]-1.))*sgnalphai[i]) 
	  - Dilog(y0i[i]/yim[i],-ReSign(y0i[i])*sgnalphai[i]) 
	  + Dilog((y0i[i]-1.)/yip[i],ReSign((y0i[i]-1.))*sgnalphai[i])
	  - Dilog(y0i[i]/yip[i],ReSign(y0i[i])*sgnalphai[i])
	  + Eta(1.-xim[i],sgnalphai[i],1./yim[i],-sgnalphai[i])*CLog((y0i[i]-1.)/yim[i],-ReSign((y0i[i]-1.))*sgnalphai[i])
	  - Eta(-xim[i],sgnalphai[i],1./yim[i],-sgnalphai[i])*CLog(y0i[i]/yim[i],-ReSign(y0i[i])*sgnalphai[i])
	  + Eta(1.-xip[i],-sgnalphai[i],1./yip[i],sgnalphai[i])*CLog((y0i[i]-1.)/yip[i],ReSign((y0i[i]-1.))*sgnalphai[i])
	  - Eta(-xip[i],-sgnalphai[i],1./yip[i],sgnalphai[i])*CLog(y0i[i]/yip[i],ReSign(y0i[i])*sgnalphai[i])
	  -(Eta(-xip[i],-sgnalphai[i],-xim[i],sgnalphai[i])
	    -Eta(yip[i],-sgnalphai[i],yim[i],sgnalphai[i])
	    -theta*2.*M_PI*Complex(0.,1.))*CLog(-(1.-y0i[i])/y0i[i]);
      }
      res /= alpha;
      return res*DivArrC(0.,0.,0.,1.,0.,0.);
    }
  }
  msg_Out() << "No triangle result found!\n"
	    << "(p12, p22, p32): (" << p12 << ", " << p22 << ", " << p32 << ")\n"
	    << "(m12, m22, m32): (" << m12 << ", " << m22 << ", " << m32 << ")\n";
  return DivArrC(0.,0.,0.,0.,0.,0.);
}





