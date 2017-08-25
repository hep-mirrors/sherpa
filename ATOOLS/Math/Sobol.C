// Frances Y. Kuo
//
// Email: <f.kuo@unsw.edu.au>
// School of Mathematics and Statistics
// University of New South Wales
// Sydney NSW 2052, Australia
// 
// Last updated: 21 October 2008
//
//   You may incorporate this source code into your own program 
//   provided that you
//   1) acknowledge the copyright owner in your program and publication
//   2) notify the copyright owner by email
//   3) offer feedback regarding your experience with different direction numbers
//
//
// -----------------------------------------------------------------------------
// Licence pertaining to sobol.cc and the accompanying sets of direction numbers
// -----------------------------------------------------------------------------
// Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
// 
//     * Neither the names of the copyright holders nor the names of the
//       University of New South Wales and the University of Waikato
//       and its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -----------------------------------------------------------------------------

// Adapted for Sherpa by Stefan Hoeche

#include "ATOOLS/Math/Sobol.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Exception.H"

#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

using namespace ATOOLS;

Sobol::Sobol(unsigned _D,std::string file)
{
  file = rpa->gen.Variable("SHERPA_SHARE_PATH")+"/Sobol/"+file;
#ifdef USING__GZIP
  file+=".gz";
  igzstream infile(file.c_str(),std::ios::in);
#else
  std::ifstream infile(file.c_str(),std::ios::in);
#endif
  if (!infile) THROW(fatal_error,"Direction file not found '"+file+"'.");
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  n = 0;
  D = _D;
  L = 32;
  d.resize(D,0);
  s.resize(D,0);
  a.resize(D,0);
  m.resize(D);
  for (unsigned j=1;j<=D-1;j++) {
    infile>>d[j]>>s[j]>>a[j];
    m[j].resize(s[j]+1);
    for (unsigned i=1;i<=s[j];i++) infile>>m[j][i];
  }
  V.resize(D);
  V[0].resize(L+1);
  for (unsigned i=1;i<=L;i++) V[0][i] = 1 << (32-i);
  for (unsigned j=1;j<D;j++) {
    V[j].resize(L+1); 
    if (L <= s[j]) {
      for (unsigned i=1;i<=L;i++) V[j][i] = m[j][i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s[j];i++) V[j][i] = m[j][i] << (32-i); 
      for (unsigned i=s[j]+1;i<=L;i++) {
	V[j][i] = V[j][i-s[j]] ^ (V[j][i-s[j]] >> s[j]); 
	for (unsigned k=1;k<=s[j]-1;k++) 
	  V[j][i] ^= (((a[j] >> (s[j]-1-k)) & 1) * V[j][i-k]); 
      }
    }
  }
  X.resize(D,0);
  for (unsigned i=0;i<10;++i) Point();
}

unsigned Sobol::C(const unsigned i)
{
  unsigned C=1;
  for (unsigned value=i;value&1;value>>=1) C++;
  return C;
}

std::vector<double> Sobol::Point()
{
  std::vector<double> p(D,0.0);
  if (n++==0) return p;
  for (unsigned j=0;j<D;j++) {
    X[j] = X[j] ^ V[j][C(n-2)];
    p[j] = (double)X[j]/pow(2.0,32);
  }
  return p;
}
