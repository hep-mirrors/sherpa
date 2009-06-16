#ifndef ATOOLS_Math_Lorentz_Ten4_H
#define ATOOLS_Math_Lorentz_Ten4_H

#include "ATOOLS/Math/MathTools.H"

namespace ATOOLS {
  template<typename Scalar> class Lorentz_Ten3;

  template<typename Scalar>
  class Lorentz_Ten4 {
    Scalar m_x[4][4][4][4];
  public:
    inline Lorentz_Ten4() {}
    inline Lorentz_Ten4(const Scalar& s) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              m_x[i][j][k][l] = s;
    }
    template<typename Scalar2>
    inline Lorentz_Ten4(const Lorentz_Ten4<Scalar2> & ten) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              m_x[i][j][k][l] = ten.at(i,j,k,l);
    }
    inline Lorentz_Ten4(
          const Scalar& x0000, const Scalar& x1000, const Scalar& x2000, const Scalar& x3000,
          const Scalar& x0100, const Scalar& x1100, const Scalar& x2100, const Scalar& x3100,
          const Scalar& x0200, const Scalar& x1200, const Scalar& x2200, const Scalar& x3200,
          const Scalar& x0300, const Scalar& x1300, const Scalar& x2300, const Scalar& x3300,

          const Scalar& x0010, const Scalar& x1010, const Scalar& x2010, const Scalar& x3010,
          const Scalar& x0110, const Scalar& x1110, const Scalar& x2110, const Scalar& x3110,
          const Scalar& x0210, const Scalar& x1210, const Scalar& x2210, const Scalar& x3210,
          const Scalar& x0310, const Scalar& x1310, const Scalar& x2310, const Scalar& x3310,

          const Scalar& x0020, const Scalar& x1020, const Scalar& x2020, const Scalar& x3020,
          const Scalar& x0120, const Scalar& x1120, const Scalar& x2120, const Scalar& x3120,
          const Scalar& x0220, const Scalar& x1220, const Scalar& x2220, const Scalar& x3220,
          const Scalar& x0320, const Scalar& x1320, const Scalar& x2320, const Scalar& x3320,

          const Scalar& x0030, const Scalar& x1030, const Scalar& x2030, const Scalar& x3030,
          const Scalar& x0130, const Scalar& x1130, const Scalar& x2130, const Scalar& x3130,
          const Scalar& x0230, const Scalar& x1230, const Scalar& x2230, const Scalar& x3230,
          const Scalar& x0330, const Scalar& x1330, const Scalar& x2330, const Scalar& x3330,

          const Scalar& x0001, const Scalar& x1001, const Scalar& x2001, const Scalar& x3001,
          const Scalar& x0101, const Scalar& x1101, const Scalar& x2101, const Scalar& x3101,
          const Scalar& x0201, const Scalar& x1201, const Scalar& x2201, const Scalar& x3201,
          const Scalar& x0301, const Scalar& x1301, const Scalar& x2301, const Scalar& x3301,

          const Scalar& x0011, const Scalar& x1011, const Scalar& x2011, const Scalar& x3011,
          const Scalar& x0111, const Scalar& x1111, const Scalar& x2111, const Scalar& x3111,
          const Scalar& x0211, const Scalar& x1211, const Scalar& x2211, const Scalar& x3211,
          const Scalar& x0311, const Scalar& x1311, const Scalar& x2311, const Scalar& x3311,

          const Scalar& x0021, const Scalar& x1021, const Scalar& x2021, const Scalar& x3021,
          const Scalar& x0121, const Scalar& x1121, const Scalar& x2121, const Scalar& x3121,
          const Scalar& x0221, const Scalar& x1221, const Scalar& x2221, const Scalar& x3221,
          const Scalar& x0321, const Scalar& x1321, const Scalar& x2321, const Scalar& x3321,

          const Scalar& x0031, const Scalar& x1031, const Scalar& x2031, const Scalar& x3031,
          const Scalar& x0131, const Scalar& x1131, const Scalar& x2131, const Scalar& x3131,
          const Scalar& x0231, const Scalar& x1231, const Scalar& x2231, const Scalar& x3231,
          const Scalar& x0331, const Scalar& x1331, const Scalar& x2331, const Scalar& x3331,

          const Scalar& x0002, const Scalar& x1002, const Scalar& x2002, const Scalar& x3002,
          const Scalar& x0102, const Scalar& x1102, const Scalar& x2102, const Scalar& x3102,
          const Scalar& x0202, const Scalar& x1202, const Scalar& x2202, const Scalar& x3202,
          const Scalar& x0302, const Scalar& x1302, const Scalar& x2302, const Scalar& x3302,

          const Scalar& x0012, const Scalar& x1012, const Scalar& x2012, const Scalar& x3012,
          const Scalar& x0112, const Scalar& x1112, const Scalar& x2112, const Scalar& x3112,
          const Scalar& x0212, const Scalar& x1212, const Scalar& x2212, const Scalar& x3212,
          const Scalar& x0312, const Scalar& x1312, const Scalar& x2312, const Scalar& x3312,

          const Scalar& x0022, const Scalar& x1022, const Scalar& x2022, const Scalar& x3022,
          const Scalar& x0122, const Scalar& x1122, const Scalar& x2122, const Scalar& x3122,
          const Scalar& x0222, const Scalar& x1222, const Scalar& x2222, const Scalar& x3222,
          const Scalar& x0322, const Scalar& x1322, const Scalar& x2322, const Scalar& x3322,

          const Scalar& x0032, const Scalar& x1032, const Scalar& x2032, const Scalar& x3032,
          const Scalar& x0132, const Scalar& x1132, const Scalar& x2132, const Scalar& x3132,
          const Scalar& x0232, const Scalar& x1232, const Scalar& x2232, const Scalar& x3232,
          const Scalar& x0332, const Scalar& x1332, const Scalar& x2332, const Scalar& x3332,

          const Scalar& x0003, const Scalar& x1003, const Scalar& x2003, const Scalar& x3003,
          const Scalar& x0103, const Scalar& x1103, const Scalar& x2103, const Scalar& x3103,
          const Scalar& x0203, const Scalar& x1203, const Scalar& x2203, const Scalar& x3203,
          const Scalar& x0303, const Scalar& x1303, const Scalar& x2303, const Scalar& x3303,

          const Scalar& x0013, const Scalar& x1013, const Scalar& x2013, const Scalar& x3013,
          const Scalar& x0113, const Scalar& x1113, const Scalar& x2113, const Scalar& x3113,
          const Scalar& x0213, const Scalar& x1213, const Scalar& x2213, const Scalar& x3213,
          const Scalar& x0313, const Scalar& x1313, const Scalar& x2313, const Scalar& x3313,

          const Scalar& x0023, const Scalar& x1023, const Scalar& x2023, const Scalar& x3023,
          const Scalar& x0123, const Scalar& x1123, const Scalar& x2123, const Scalar& x3123,
          const Scalar& x0223, const Scalar& x1223, const Scalar& x2223, const Scalar& x3223,
          const Scalar& x0323, const Scalar& x1323, const Scalar& x2323, const Scalar& x3323,

          const Scalar& x0033, const Scalar& x1033, const Scalar& x2033, const Scalar& x3033,
          const Scalar& x0133, const Scalar& x1133, const Scalar& x2133, const Scalar& x3133,
          const Scalar& x0233, const Scalar& x1233, const Scalar& x2233, const Scalar& x3233,
          const Scalar& x0333, const Scalar& x1333, const Scalar& x2333, const Scalar& x3333
            ) {
   m_x[0][0][0][0]=x0000; m_x[1][0][0][0]=x1000; m_x[2][0][0][0]=x2000; m_x[3][0][0][0]=x3000;
   m_x[0][1][0][0]=x0100; m_x[1][1][0][0]=x1100; m_x[2][1][0][0]=x2100; m_x[3][1][0][0]=x3100;
   m_x[0][2][0][0]=x0200; m_x[1][2][0][0]=x1200; m_x[2][2][0][0]=x2200; m_x[3][2][0][0]=x3200;
   m_x[0][3][0][0]=x0300; m_x[1][3][0][0]=x1300; m_x[2][3][0][0]=x2300; m_x[3][3][0][0]=x3300;

   m_x[0][0][1][0]=x0010; m_x[1][0][1][0]=x1010; m_x[2][0][1][0]=x2010; m_x[3][0][1][0]=x3010;
   m_x[0][1][1][0]=x0110; m_x[1][1][1][0]=x1110; m_x[2][1][1][0]=x2110; m_x[3][1][1][0]=x3110;
   m_x[0][2][1][0]=x0210; m_x[1][2][1][0]=x1210; m_x[2][2][1][0]=x2210; m_x[3][2][1][0]=x3210;
   m_x[0][3][1][0]=x0310; m_x[1][3][1][0]=x1310; m_x[2][3][1][0]=x2310; m_x[3][3][1][0]=x3310;

   m_x[0][0][2][0]=x0020; m_x[1][0][2][0]=x1020; m_x[2][0][2][0]=x2020; m_x[3][0][2][0]=x3020;
   m_x[0][1][2][0]=x0120; m_x[1][1][2][0]=x1120; m_x[2][1][2][0]=x2120; m_x[3][1][2][0]=x3120;
   m_x[0][2][2][0]=x0220; m_x[1][2][2][0]=x1220; m_x[2][2][2][0]=x2220; m_x[3][2][2][0]=x3220;
   m_x[0][3][2][0]=x0320; m_x[1][3][2][0]=x1320; m_x[2][3][2][0]=x2320; m_x[3][3][2][0]=x3320;

   m_x[0][0][3][0]=x0030; m_x[1][0][3][0]=x1030; m_x[2][0][3][0]=x2030; m_x[3][0][3][0]=x3030;
   m_x[0][1][3][0]=x0130; m_x[1][1][3][0]=x1130; m_x[2][1][3][0]=x2130; m_x[3][1][3][0]=x3130;
   m_x[0][2][3][0]=x0230; m_x[1][2][3][0]=x1230; m_x[2][2][3][0]=x2230; m_x[3][2][3][0]=x3230;
   m_x[0][3][3][0]=x0330; m_x[1][3][3][0]=x1330; m_x[2][3][3][0]=x2330; m_x[3][3][3][0]=x3330;

   m_x[0][0][0][1]=x0001; m_x[1][0][0][1]=x1001; m_x[2][0][0][1]=x2001; m_x[3][0][0][1]=x3001;
   m_x[0][1][0][1]=x0101; m_x[1][1][0][1]=x1101; m_x[2][1][0][1]=x2101; m_x[3][1][0][1]=x3101;
   m_x[0][2][0][1]=x0201; m_x[1][2][0][1]=x1201; m_x[2][2][0][1]=x2201; m_x[3][2][0][1]=x3201;
   m_x[0][3][0][1]=x0301; m_x[1][3][0][1]=x1301; m_x[2][3][0][1]=x2301; m_x[3][3][0][1]=x3301;

   m_x[0][0][1][1]=x0011; m_x[1][0][1][1]=x1011; m_x[2][0][1][1]=x2011; m_x[3][0][1][1]=x3011;
   m_x[0][1][1][1]=x0111; m_x[1][1][1][1]=x1111; m_x[2][1][1][1]=x2111; m_x[3][1][1][1]=x3111;
   m_x[0][2][1][1]=x0211; m_x[1][2][1][1]=x1211; m_x[2][2][1][1]=x2211; m_x[3][2][1][1]=x3211;
   m_x[0][3][1][1]=x0311; m_x[1][3][1][1]=x1311; m_x[2][3][1][1]=x2311; m_x[3][3][1][1]=x3311;

   m_x[0][0][2][1]=x0021; m_x[1][0][2][1]=x1021; m_x[2][0][2][1]=x2021; m_x[3][0][2][1]=x3021;
   m_x[0][1][2][1]=x0121; m_x[1][1][2][1]=x1121; m_x[2][1][2][1]=x2121; m_x[3][1][2][1]=x3121;
   m_x[0][2][2][1]=x0221; m_x[1][2][2][1]=x1221; m_x[2][2][2][1]=x2221; m_x[3][2][2][1]=x3221;
   m_x[0][3][2][1]=x0321; m_x[1][3][2][1]=x1321; m_x[2][3][2][1]=x2321; m_x[3][3][2][1]=x3321;

   m_x[0][0][3][1]=x0031; m_x[1][0][3][1]=x1031; m_x[2][0][3][1]=x2031; m_x[3][0][3][1]=x3031;
   m_x[0][1][3][1]=x0131; m_x[1][1][3][1]=x1131; m_x[2][1][3][1]=x2131; m_x[3][1][3][1]=x3131;
   m_x[0][2][3][1]=x0231; m_x[1][2][3][1]=x1231; m_x[2][2][3][1]=x2231; m_x[3][2][3][1]=x3231;
   m_x[0][3][3][1]=x0331; m_x[1][3][3][1]=x1331; m_x[2][3][3][1]=x2331; m_x[3][3][3][1]=x3331;

   m_x[0][0][0][2]=x0002; m_x[1][0][0][2]=x1002; m_x[2][0][0][2]=x2002; m_x[3][0][0][2]=x3002;
   m_x[0][1][0][2]=x0102; m_x[1][1][0][2]=x1102; m_x[2][1][0][2]=x2102; m_x[3][1][0][2]=x3102;
   m_x[0][2][0][2]=x0202; m_x[1][2][0][2]=x1202; m_x[2][2][0][2]=x2202; m_x[3][2][0][2]=x3202;
   m_x[0][3][0][2]=x0302; m_x[1][3][0][2]=x1302; m_x[2][3][0][2]=x2302; m_x[3][3][0][2]=x3302;

   m_x[0][0][1][2]=x0012; m_x[1][0][1][2]=x1012; m_x[2][0][1][2]=x2012; m_x[3][0][1][2]=x3012;
   m_x[0][1][1][2]=x0112; m_x[1][1][1][2]=x1112; m_x[2][1][1][2]=x2112; m_x[3][1][1][2]=x3112;
   m_x[0][2][1][2]=x0212; m_x[1][2][1][2]=x1212; m_x[2][2][1][2]=x2212; m_x[3][2][1][2]=x3212;
   m_x[0][3][1][2]=x0312; m_x[1][3][1][2]=x1312; m_x[2][3][1][2]=x2312; m_x[3][3][1][2]=x3312;

   m_x[0][0][2][2]=x0022; m_x[1][0][2][2]=x1022; m_x[2][0][2][2]=x2022; m_x[3][0][2][2]=x3022;
   m_x[0][1][2][2]=x0122; m_x[1][1][2][2]=x1122; m_x[2][1][2][2]=x2122; m_x[3][1][2][2]=x3122;
   m_x[0][2][2][2]=x0222; m_x[1][2][2][2]=x1222; m_x[2][2][2][2]=x2222; m_x[3][2][2][2]=x3222;
   m_x[0][3][2][2]=x0322; m_x[1][3][2][2]=x1322; m_x[2][3][2][2]=x2322; m_x[3][3][2][2]=x3322;

   m_x[0][0][3][2]=x0032; m_x[1][0][3][2]=x1032; m_x[2][0][3][2]=x2032; m_x[3][0][3][2]=x3032;
   m_x[0][1][3][2]=x0132; m_x[1][1][3][2]=x1132; m_x[2][1][3][2]=x2132; m_x[3][1][3][2]=x3132;
   m_x[0][2][3][2]=x0232; m_x[1][2][3][2]=x1232; m_x[2][2][3][2]=x2232; m_x[3][2][3][2]=x3232;
   m_x[0][3][3][2]=x0332; m_x[1][3][3][2]=x1332; m_x[2][3][3][2]=x2332; m_x[3][3][3][2]=x3332;

   m_x[0][0][0][3]=x0003; m_x[1][0][0][3]=x1003; m_x[2][0][0][3]=x2003; m_x[3][0][0][3]=x3003;
   m_x[0][1][0][3]=x0103; m_x[1][1][0][3]=x1103; m_x[2][1][0][3]=x2103; m_x[3][1][0][3]=x3103;
   m_x[0][2][0][3]=x0203; m_x[1][2][0][3]=x1203; m_x[2][2][0][3]=x2203; m_x[3][2][0][3]=x3203;
   m_x[0][3][0][3]=x0303; m_x[1][3][0][3]=x1303; m_x[2][3][0][3]=x2303; m_x[3][3][0][3]=x3303;

   m_x[0][0][1][3]=x0013; m_x[1][0][1][3]=x1013; m_x[2][0][1][3]=x2013; m_x[3][0][1][3]=x3013;
   m_x[0][1][1][3]=x0113; m_x[1][1][1][3]=x1113; m_x[2][1][1][3]=x2113; m_x[3][1][1][3]=x3113;
   m_x[0][2][1][3]=x0213; m_x[1][2][1][3]=x1213; m_x[2][2][1][3]=x2213; m_x[3][2][1][3]=x3213;
   m_x[0][3][1][3]=x0313; m_x[1][3][1][3]=x1313; m_x[2][3][1][3]=x2313; m_x[3][3][1][3]=x3313;

   m_x[0][0][2][3]=x0023; m_x[1][0][2][3]=x1023; m_x[2][0][2][3]=x2023; m_x[3][0][2][3]=x3023;
   m_x[0][1][2][3]=x0123; m_x[1][1][2][3]=x1123; m_x[2][1][2][3]=x2123; m_x[3][1][2][3]=x3123;
   m_x[0][2][2][3]=x0223; m_x[1][2][2][3]=x1223; m_x[2][2][2][3]=x2223; m_x[3][2][2][3]=x3223;
   m_x[0][3][2][3]=x0323; m_x[1][3][2][3]=x1323; m_x[2][3][2][3]=x2323; m_x[3][3][2][3]=x3323;

   m_x[0][0][3][3]=x0033; m_x[1][0][3][3]=x1033; m_x[2][0][3][3]=x2033; m_x[3][0][3][3]=x3033;
   m_x[0][1][3][3]=x0133; m_x[1][1][3][3]=x1133; m_x[2][1][3][3]=x2133; m_x[3][1][3][3]=x3133;
   m_x[0][2][3][3]=x0233; m_x[1][2][3][3]=x1233; m_x[2][2][3][3]=x2233; m_x[3][2][3][3]=x3233;
   m_x[0][3][3][3]=x0333; m_x[1][3][3][3]=x1333; m_x[2][3][3][3]=x2333; m_x[3][3][3][3]=x3333;
    }
    inline Lorentz_Ten4(const Scalar ten[4][4][4][4]) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              m_x[i][j][k][l] = ten[i][j][k][l];
    }

    // element extraction
    inline const Scalar at(unsigned short int i, unsigned short int j,
                           unsigned short int k, unsigned short int l) const {
      return m_x[i][j][k][l];
    }

//     typedef Scalar (*ScalarArray)[4][4][4];
//     inline ScalarArray operator[] (unsigned short int i) {
//       return m_x[i];
//     }

//     inline Lorentz_Ten3<Scalar> operator[] (unsigned short int i) {
//       return Lorentz_Ten3<Scalar>(***m_x[i]);
//     }

//     inline const Scalar *** operator[] (unsigned short int i) const {
//       return m_x[i];
//     }

    // sign flip
    inline Lorentz_Ten4<Scalar> operator-() const {
      Scalar x[4][4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              x[i][j][k][l] = -m_x.at(i,j,k,l);
      return Lorentz_Ten4<Scalar>(x);
    }

    // transpose pair of indizes
    inline Lorentz_Ten4<Scalar> Transpose(unsigned short int a, unsigned short int b) const {
      Scalar x[4][4][4][4];
      if      (((a==1) && (b==2)) || ((a==2) && (b==1)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
              for (unsigned short int l=0; l<4; ++l)
                x[i][j][k][l] = m_x[j][i][k][l];
      else if (((a==1) && (b==3)) || ((a==3) && (b==1)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
              for (unsigned short int l=0; l<4; ++l)
                x[i][j][k][l] = m_x[k][j][i][l];
      else if (((a==1) && (b==4)) || ((a==4) && (b==1)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
              for (unsigned short int l=0; l<4; ++l)
                x[i][j][k][l] = m_x[l][j][k][i];
      else if (((a==2) && (b==3)) || ((a==3) && (b==2)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
              for (unsigned short int l=0; l<4; ++l)
                x[i][j][k][l] = m_x[i][k][j][l];
      else if (((a==2) && (b==4)) || ((a==4) && (b==2)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
              for (unsigned short int l=0; l<4; ++l)
                x[i][j][k][l] = m_x[i][l][k][j];
      else if (((a==3) && (b==4)) || ((a==4) && (b==3)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
              for (unsigned short int l=0; l<4; ++l)
                x[i][j][k][l] = m_x[i][j][l][k];
      return Lorentz_Ten4<Scalar>(x);
    }

    // addition/subtraction operators
    inline Lorentz_Ten4<Scalar>& operator+=(const Lorentz_Ten4<Scalar>& t) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              m_x[i][j][k][l] += t.at(i,j,k,l);
      return *this;
    }

    inline Lorentz_Ten4<Scalar>& operator-=(const Lorentz_Ten4<Scalar>& t) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              m_x[i][j][k][l] -= t.at(i,j,k,l);
      return *this;
    }

    template<typename Scalar2> 
    inline Lorentz_Ten4<PROMOTE(Scalar,Scalar2)>
    operator+ (const Lorentz_Ten4<Scalar2>& ten) const {
      PROMOTE(Scalar,Scalar2) x[4][4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              x[i][j][k][l] = m_x[i][j][k][l]+ten.at(i,j,k,l);
      return Lorentz_Ten4<PROMOTE(Scalar,Scalar2)>(x);
    }

    template<typename Scalar2> 
    inline Lorentz_Ten4<PROMOTE(Scalar,Scalar2)>
    operator- (const Lorentz_Ten4<Scalar2>& ten) const {
      PROMOTE(Scalar,Scalar2) x[4][4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              x[i][j][k][l] = m_x[i][j][k][l]-ten.at(i,j,k,l);
      return Lorentz_Ten4<PROMOTE(Scalar,Scalar2)>(x);
    }

    // multiplication operators
    inline Lorentz_Ten4<Scalar>& operator*= (const Scalar& scal) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              m_x[i][j][k][l] *= scal;
      return *this;
    }

    template<typename Scalar2> 
    inline Lorentz_Ten4<PROMOTE(Scalar,Scalar2)> operator* (const Scalar2 scal) const {
      PROMOTE(Scalar,Scalar2) x[4][4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              x[i][j][k][l] = scal*m_x[i][j][k][l];
      return Lorentz_Ten4<PROMOTE(Scalar,Scalar2)>(x);
    }

    template<typename Scalar2> 
    inline Lorentz_Ten4<PROMOTE(Scalar,Scalar2)> operator/ (const Scalar2 scal) const {
      return (*this)*(1./scal);
    }

    inline bool Nan() const {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              if (ATOOLS::IsNan<Scalar>(m_x[i][j][k][l])) return true;
      return false;
    }

    inline bool IsZero() const {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
              if (!ATOOLS::IsZero<Scalar>(m_x[i][j][k][l])) return false;
      return true;
    }
  };

  template<typename Scalar,typename Scal2> 
  inline Lorentz_Ten4<PROMOTE(Scalar,Scal2)> 
  operator* (const Scal2 s, const Lorentz_Ten4<Scalar>& t) {
    return Lorentz_Ten4<PROMOTE(Scalar,Scal2)>(t*s);
  }

  template<typename Scalar>
  std::ostream& operator<<(std::ostream& s, const Lorentz_Ten4<Scalar>& ten) {
    return s<<"not implemented for tensors of fourth rank ... "
            <<"please project to second rank and display then ..."<<std::endl;
  }
}

/*!
 \file
 \brief   contains class Lorentz_Ten4
*/

/*!
 \class Lorentz_Ten4<Scalar>
 \brief implementation of a 4 dimensional Minkowski 2nd rank tensor and its algebraic
 operations

 This class can be used as Minkowski 2nd rank tensor with arbitrary scalar types.
 All necessary operations, e.g. addition, scalar multiplication, contractions
 with other 2nd rank tensor or vectors, etc. are available.
 If you want to use this 2nd rank tensor for a new scalar type, you might have to
 implement specific functions for this type in MathTools.H.
 "Fallback" types for operations with two 2nd rank tensors of different types, e. g.
 "complex and double become complex" have to be specified in MathTools.H as
 well, see the DECLARE_PROMOTE macro.
*/

/*!
 \fn Lorentz_Ten4::Lorentz_Ten4()
 \brief Standard Constructor
*/

/*!
 \fn Lorentz_Ten4::Lorentz_Ten4(Scalar x00, Scalar x10, Scalar x20, Scalar x30, Scalar x01, Scalar x11, Scalar x21, Scalar x31, Scalar x02, Scalar x12, Scalar x22, Scalar x32, Scalar x03, Scalar x13, Scalar x23, Scalar x33){
 \brief Special Constructor, templated in Scalar, taking 16 single components.
*/

/*!
 \fn   inline Scalar* Lorentz_Ten4::operator[] (int i)
 \brief returns array \f$x_{ij}\f$, \f$j\in[0,3]\f$. (May be manipulated.)
*/

/*!
 \fn   inline const Scalar* Lorentz_Ten4::operator[] (int i) const
 \brief returns array \f$x_{ij}\f$, \f$j\in[0,3]\f$.
*/

#endif
