BEGIN { pi=0; i=0; j=0; pb=3.89379656e8; end1=0; warnings=0; errors=0; cf=0; } 
{
  if (pi==0) { 
    if ($1!="") ++pi;
    filename=$1".xsd.dat"; 
    htmlname=$1"/index.html";
  } 
  else if (pi==1) {
    if ($1!="") ++pi;
    genone=$0;
  } 
  else if (pi==2) {
    if ($1!="") ++pi;
    gentwo=$0;
  } 
  else if (pi==3) {
    if ($1!="") ++pi;
    printf "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD " \
      "HTML 4.01 Transitional//EN\">\n" > htmlname
    printf "<html>\n  <head>\n" > htmlname
    printf "    <title>xs comparison</title>\n" > htmlname
    printf "  <head>\n  <body>\n    <hr size=\"1\">\n" > htmlname
    printf "    <center><font color=\"#0000ff\" size=\"+2\"><b>\n      "$0 \
      "</b></font></center>\n" > htmlname
    printf "    <hr size=\"1\"><br>\n    <table width=\"100%\" border=\"1\"" \
      " bordercolor=\"#888888\">\n" > htmlname 
    printf "      <tr bgcolor=\"#bbbbbb\"><td><b>Process</b></td>\n" > htmlname
    printf "        <td><center><b>XS from "genone \
      " [pb]</b></center></td>\n" > htmlname
    printf "        <td><center><b>XS from "gentwo \
      " [pb]</b></center></td>\n" > htmlname
    printf "        <td><center><b>w<sub>2</sub>/w<sub>1</sub>-1 [%]" \
      "</b></center></td>\n" > htmlname
    printf "        <td><center><b>(w<sub>2</sub>-w<sub>1</sub>)/\n" \
      "          (&sigma;<sub>1</sub><sup>2</sup>+" \
      "&sigma;<sub>2</sub><sup>2</sup>)\n          <sup>1/2</sup>" \
      "</b></center></td></tr>\n" > htmlname
  } 
  else { 
    if ($1=="end") {
      end1=1; 
      cf=0;
    }
    else {
    if (end1==0) { 
      proc1[i]=$1; 
      if (proc1[i]=="#") {
	proc1[i]=$0;
      }
      else if (proc1[i]=="%%") {
        cf=$2;
	--i;
      }
      else {
        file1[i]=cf;
        xs1[i]=$2;  
        max1[i]=$3;  
        err1[i]=$4;  
        relerr1[i]=$4/$2;  
#        printf "1st file: process: "$1", xs = "xs1[i]*pb" pb, err = " \
#          err1[i]*pb" ( "relerr1[i]*100"% ), max = "max1[i]*pb" pb\n"; 
      }
      ++i; 
    } 
    else { 
      proc2[j]=$1; 
      if (proc2[j]=="#") {
        proc2[j]=$0;
      }
      else if (proc2[j]=="%%") {
        cf=$2;
	--j;
      }
      else {
	file2[j]=cf;
	xs2[j]=$2;  
	max2[j]=$3;  
	err2[j]=$4;  
	relerr2[j]=$4/$2;  
#        printf "2nd file: process: "$1", xs = "xs2[j]*pb" pb, err = " \
#          err2[j]*pb" ( "relerr2[j]*100"% ), max = "max2[j]*pb" pb\n"; 
      }
      ++j; 
    } 
  } 
  }
} 
END { 
  min=-3.; max=3.; bins=31; binwidth=(max-min)/bins; 
  for (k=0;k<=bins;++k) { 
    histox[k]=min+binwidth*k; 
    histoy[k]=0; 
  } 
  devsum = 0.; 
  devavg = 0.; 
  pss=0;
  for (ii=0;ii<i;++ii) { 
    match(proc1[ii],"#");
    if (RSTART==1) {
      printf "      <tr><td colspan=\"5\"><center>"substr(proc1[ii],2) \
        "</center></td></tr>\n" > htmlname
      continue;
    }
    for (jj=0;jj<j;++jj) { 
      if (file1[ii]!=file2[jj] || proc1[ii]!=proc2[jj]) continue; 
      ++pss;
      meanerr=sqrt(err1[ii]*err1[ii]+err2[jj]*err2[jj]); 
      devvar=(xs2[jj]-xs1[ii])/meanerr; 
      devsum += devvar*devvar; 
      devavg += devvar; 
      reldev=devvar; 
      if (reldev<0) reldev=-reldev; 
      for (k=1;k<bins;++k) { 
        if (devvar>=histox[k] && devvar<histox[k+1]) ++histoy[k]; 
      } 
      if (devvar<histox[0]) ++histoy[0]; 
      if (devvar>=histox[bins]) ++histoy[bins]; 
      printf "test process: \033[1m"proc1[ii] \
        "\033[0m, rel deviation = "; 
      if (reldev>1.0) { 
	if (reldev>2.0) { 
	  printf "      <tr bgcolor=\"#ffcccc\">" > htmlname
	  printf "\033[41m";
	}
	else {
	  printf "      <tr bgcolor=\"#ffffcc\">" > htmlname
	  printf "\033[31m";
	}
      }
      else {
	printf "      <tr bgcolor=\"#ccffcc\">" > htmlname
	printf "\033[34m";
      }
      printf "<td><b>"proc1[ii]"</b></td>\n" > htmlname
      printf "        <td><center>"xs1[ii]*pb" +- "err1[ii]*pb \
        " ( "relerr1[ii]*100"% )</center></td>\n" > htmlname
      printf "        <td><center>"xs2[jj]*pb" +- "err2[jj]*pb \
        " ( "relerr2[jj]*100"% )</center></td>\n" > htmlname
      printf "        <td><center>"(xs2[jj]/xs1[ii]-1)*100 \
        "</center></td>\n        " > htmlname
      printf reldev"\033[0m sigma vs. rel errors = \033[32m" \
	relerr1[ii]*100"%\033[0m, \033[32m"relerr2[jj]*100"%\033[0m\n"; 
      if (reldev>1.0) { 
        if (reldev>2.0) { 
	  printf "<td><center><font color=\"#aa0000\"><b>"devvar \
            "</b></font></center></td></tr>\n" > htmlname
	  ++errors; 
        } 
        else { 
	  printf "<td><center><font color=\"#dddd00\"><b>"devvar \
            "</b></font></center></td></tr>\n" > htmlname
	  ++warnings; 
        }
      }
      else {
	printf "<td><center><font color=\"#00aa00\"><b>"devvar	\
	  "</b></font></center></td></tr>\n" > htmlname
      }
      break; 
    }
  }
  printf "finished test with "errors \
    " errors and "warnings" warnings in "pss" processes\n"; 
  if (pss>1) printf "mean sigma is "devavg/pss \
    ", delta sigma is "sqrt((devsum-devavg*devavg/pss)/(pss-1))"\n"; 
  else printf "Only one process."; 
  printf "write deviation histo to "filename"\n"; 
  for (k=0;k<=bins;++k) { 
    printf histox[k]" "histoy[k]"\n" > filename; 
  } 
  printf "    </table><br>\n    <hr size=\"1\">\n" > htmlname
  printf "    Test yields "errors \
    " errors and "warnings" warnings in "pss" processes" > htmlname 
  if (pss>1) printf ",\n    &lt;&sigma;&gt; = "devavg/pss \
    ", &Delta;&sigma; = " sqrt((devsum-devavg*devavg/pss)/(pss-1)) > htmlname
  printf "\n    <hr size=\"1\">\n    " > htmlname
  printf "<img src=\"Dev_Stat_1.gif\"/>\n    " > htmlname
  printf "<hr size=\"1\">\n    <br><br>\n  </body>\n</html>\n" > htmlname
  printf "wrote data to "htmlname"\n"; 
}
