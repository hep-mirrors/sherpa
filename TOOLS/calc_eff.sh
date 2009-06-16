#!/bin/bash

print_help() 
{
  echo "calc_eff version 1.0" && echo && \
  echo "options: -e <epsilon>   set epsilon <epsilon>" && \
  echo "         -p <path>      use results in <path>" && \
  echo "         -n <name>      use process name <name>" && \
  echo "         -v             increase verbosity" && \
  echo "         -q             report efficiency only" && \
  echo "         -t             compute efficiency for triggered events" && \
  echo "         -h             display this help and exit" && \
  echo
}

rpath="Results";
eps="0.001";
pname="*";
qmode=0;
tmode=0;
verb=0;
while getopts :e:p:n:vtqh OPT
do
  case $OPT in
  e) eps=$OPTARG ;;
  p) rpath=$OPTARG ;; 
  n) pname=$OPTARG ;; 
  v) let ++verb ;; 
  q) qmode=1 ;; 
  t) tmode=1 ;;
  h) print_help && exit 0 ;;
  \?)
    shift `expr $OPTIND - 1`
    echo -n "calc_eff: error: unrecognized option "
    if [ $OPTARG != "-" ]; then echo "'-$OPTARG'. try '--help'"
    else echo "'$1'. try '--help'"
    fi
    print_help && exit 1
  esac
done

if ! test -d $(echo $rpath | cut -d' ' -f1); then
  echo "calc_eff: error: result directory '$rpath' not found";
  exit;
fi;

pb=3.89379656e8

declare -a sum mean;

ress=$(echo $rpath/$pname.xs_tot);

sum[0]=0; mean[0]=0; n=0;
for i in $ress; do
  cpn=$(echo $i | awk '{ n=split($1,a,"/"); split(a[n],b,".xs_tot"); print b[1]; }');
  effs=$(echo $rpath/WD_$cpn/*);
  if ! test -f $(echo $effs | cut -d' ' -f1); then
    echo "calc_eff: error: weight histos not found in '$rpath/WD_$cpn/'";
    exit;
  fi;
  for j in $effs; do
    cur=$(echo $j | awk '{ n=split($1,a,"/"); 
      if (match(a[n],"SG") || match(a[n],"j") || match(a[n],"Q") || 
          match(a[n],"fermion") || match(a[n],"neutrino") || match(a[n],"lepton")) exit 0;
      print a[n]; }');
    if test -z "$cur"; then continue; fi;
    let ++n;
    ent=$(grep $cur $i);
    if test -z "$ent"; then 
      echo "calc_eff: error: subprocess xs for '$cur' not found"; 
      exit 1;
    fi;
    xs=$(echo $ent | awk '{ print $2*'$pb'; }' );
    max=$(awk 'BEGIN{ eps='$eps'; n=0; }{
	if (NF!=2) {
	  bw2=($4-$3)/($2-2)/2;
	  ent=$7;
	  ns=$5;
	  next;
        }
        ++n; w[n]=$1; e[n]=$2;
	ns+=e[n];
      }END{ wsum=nsum=0;
        for (i=n;i>=1;--i) {
	  nsum+=e[i];
	  wsum+=exp(log(10)*(w[i]+bw2))*'$pb'*e[i];
	  cwmax=exp(log(10)*w[i])*'$pb';
	  if (wsum-nsum*cwmax>eps*'$xs'*ent) break;
        }
        if (i==n) --i;
	if ('$tmode'==0) ns=ent;
	max=exp(log(10)*w[i+1])*ns/ent;	
	print max*'$pb'; }' $j);
    eff=$(echo $xs $max | awk '{ print $1/$2; }' );
    sum[n]=$(echo ${sum[$n-1]} $xs | awk '{ print $1+$2; }');
    mean[n]=$(echo ${mean[$n-1]} $xs $eff | awk '{ print $1+$2*$3; }');
    if test $verb -gt 0; then
      echo "calc_eff: '$cur' -> <w>/w_{max}^{$eps} = $eff, "\
        "\sigma = $xs pb"; fi;
    if test $verb -gt 1; then
      echo "          => \sum\sigma*<w>/w_{max}^{$eps} = ${mean[n]} pb,"\
        "\sum\sigma = ${sum[n]} pb ->"\
        "$(echo ${mean[n]} ${sum[n]} | awk '{ print $1/$2; }')"; fi;
  done;
done;
effeps=$(echo ${mean[n]} ${sum[n]} | awk '{ print $1/$2; }');
if test $qmode -eq 1; then
  echo $effeps;
else
  echo "calc_eff: '$rpath'/'$pname' -> "\
    "\sum\sigma*<w>/w_{max}^{$eps} / \sum\sigma = $effeps";
fi;
