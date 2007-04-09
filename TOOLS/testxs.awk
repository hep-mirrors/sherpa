BEGIN { pi=0; i=0; j=0; pb=3.89379656e8; end1=0; warnings=0; errors=0 } 
{
  if (pi==0) { 
    if ($1!="") pi=1;
    filename=$1".xsd.dat"; 
    htmlname=$1".dat.html";
    genone=$2;
    gentwo=$3;
    printf "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD " \
      "HTML 4.01 Transitional//EN\">\n" > htmlname
    printf "<html>\n  <head>\n" > htmlname
    printf "    <title>xs comparison</title>\n" > htmlname
    printf "  <head>\n  <body>\n" > htmlname
    printf "  <center><h1><font color=\"#0000ff\">"$1 \
      " cross section comparison</color></h1></center>\n" > htmlname
    printf "  <table width=\"100%\" border=\"1\"" \
      " bordercolor=\"#888888\">\n" > htmlname 
    printf "    <tr bgcolor=\"#bbbbbb\"><td><b>Process</b></td>" > htmlname
    printf "<td><b>XS from "genone" [pb]</b></td>" > htmlname
    printf "<td><b>XS from "gentwo" [pb]</b></td>" > htmlname
    printf "<td><center><b>w<sub>2</sub>/w<sub>1</sub>-1 [%]" \
      "</b></center></td>" > htmlname
    printf "<td><center><b>(w<sub>2</sub>-w<sub>1</sub>)/" \
      "(&sigma;<sub>1</sub>+&sigma;<sub>2</sub>)" \
      "<sup>1/2</sup></b></center></td></tr>\n" > htmlname
  } 
  else { 
  if ($1=="end") end1=1; 
  else { 
    if (end1==0) { 
      proc1[i]=$1; 
      xs1[i]=$2;  
      max1[i]=$3;  
      err1[i]=$4;  
      relerr1[i]=$4/$2;  
#      printf "1st file: process: "$1", xs = "xs1[i]*pb" pb, err = " \
#        err1[i]*pb" ( "relerr1[i]*100"% ), max = "max1[i]*pb" pb\n"; 
      ++i; 
    } 
    else { 
      proc2[j]=$1; 
      xs2[j]=$2;  
      max2[j]=$3;  
      err2[j]=$4;  
      relerr2[j]=$4/$2;  
#      printf "2nd file: process: "$1", xs = "xs2[j]*pb" pb, err = " \
#        err2[j]*pb" ( "relerr2[j]*100"% ), max = "max2[j]*pb" pb\n"; 
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
  for (ii=0;ii<i;++ii) { 
    for (jj=0;jj<j;++jj) { 
      if (proc1[ii]!=proc2[jj]) continue; 
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
        "\033[0m, rel deviation = \033[34m"reldev \
        "\033[0m sigma vs. rel errors = \033[32m"relerr1[ii]*100 \
        "%\033[0m, \033[32m"relerr2[jj]*100"%\033[0m\n"; 
      if (reldev>1.0) { 
	if (reldev>2.0) { 
	  printf "    <tr bgcolor=\"#ffcccc\">" > htmlname
	}
	else {
	  printf "    <tr bgcolor=\"#ffffcc\">" > htmlname
	}
      }
      else {
	printf "    <tr bgcolor=\"#ccffcc\">" > htmlname
      }
      printf "<td><b>"proc1[ii]"</b></td>" > htmlname
      printf "<td>"xs1[ii]*pb" +- "err1[ii]*pb \
        " ( "relerr1[ii]*100"% )</td>" > htmlname
      printf "<td>"xs2[jj]*pb" +- "err2[jj]*pb \
        " ( "relerr2[jj]*100"% )</td>" > htmlname
      printf "<td><center>"(xs2[jj]/xs1[ii]-1)*100"</center></td>" > htmlname
      if (reldev>1.0) { 
        if (reldev>2.0) { 
          printf "==================================================\n"; 
          printf "\033[41mERROR\033[0m: rel deviation > 2 sigma in \033[1m" \
            proc1[ii]"\033[0m\n"; 
          printf "==================================================\n"; 
	  printf "<td><center><font color=\"#aa0000\"><b>"devvar \
            "</b></font></center></td></tr>\n" > htmlname
	  ++errors; 
        } 
        else { 
          printf "--------------------------------------------------\n"; 
          printf "\033[31mWARNING\033[0m: relative deviation " \
            "> 1 sigma in \033[1m"proc1[ii]"\033[0m\n"; 
          printf "--------------------------------------------------\n"; 
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
    " errors and "warnings" warnings in "i" processes\n"; 
  if (k>1) printf "Average was "devavg/k \
    ", mean sigma^2 was "(devsum - devavg*devavg/i)/(i-1)"\n"; 
  else printf "Only one process."; 
  printf "write deviation histo to "filename"\n"; 
  for (k=0;k<=bins;++k) { 
    printf histox[k]" "histoy[k]"\n" > filename; 
  } 
  printf "    </table>\n  </body>\n</html>\n" > htmlname
  printf "wrote data to "htmlname"\n"; 
}
