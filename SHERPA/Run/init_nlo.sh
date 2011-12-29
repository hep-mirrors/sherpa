#!/bin/bash

if test $# -ne 2; then 
  echo "usage: $0 <sherpa exe> <input file>";
  exit 1;
fi;
tn=$(mktemp); nt=$(echo $tn | sed -e's/\/tmp\///g');
tp=$(grep NLO_QCD $2 | sed -e's/.*NLO_QCD_Part[ \t]*\(\w*\).*/\1/g');
if test $tp = RS; then
  if ! grep -q "MEH_RSADD[ \t=]*0" $2; then
    echo "Input file must contain 'MEH_RSADD 0;'";
    rm $tn;
    exit 1;
  fi;
fi;
sed -e's/}(run)/  ONLY_MAPPING_FILE 1;\n}(run)/g' < $2 > $2.$tp;
sed -e'/NLO_QCD/ d' < $2.$tp > $2.B;
export SHERPA_CPP_PATH=$PWD/$nt;
mkdir $nt; cp -r Process/ $nt/;
$1 -f$2.B;
for i in $nt/Process/Comix/*[^\)].map; do
  if grep -q x $i; then
    sed -e's/ /__QCD('$tp') /g' $i > $i.tmp;
    mv $i.tmp $(echo $i | sed -e's/.map/__QCD('$tp').map/g');
  else
    sed -e'1 s/ /__QCD('$tp') /g' -e'1 s/$/__QCD('$tp')/g' $i > $i.tmp;
    if awk '{ if ($1!=$2) exit 1; exit 0; }' < $i; then rm $i.tmp;
    else mv $i.tmp $(echo $i | sed -e's/.map/__QCD('$tp').map/g'); fi;
  fi
done;
$1 -f$2.$tp;
cp -ur $nt/Process .;
rm -rf $nt;
rm $tn