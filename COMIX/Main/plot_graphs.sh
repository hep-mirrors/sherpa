#!/bin/bash
plot(){
texs=(*.tex)
test -f ${texs[0]} || return
echo "" | sed '1 i \
%.aux: %.tex \
\	latex $(addsuffix .tex,$(basename $<)) \
%_fg.log: %.tex %.aux \
\	if test -f $(addsuffix .mp,$(basename $@)); then \\\
\	  mpost $(addsuffix .mp,$(basename $@)); \\\
\	elif test -f $(addsuffix .mf,$(basename $@)); then \\\
\	  mf $(addsuffix .mf,$(basename $@)); fi \
%.ps: %.tex %.aux %_fg.log \
\	make $(addsuffix _fg.log,$(basename $<)); lc=3; \\\
\	while egrep -s 'Process' $(addsuffix .log,$(basename $<)) \\\
\	  && [ $$lc -gt 0 ] ; do latex $<; lc=`expr $$lc - 1`; done; \\\
\	dvips -o $(addsuffix .ps,$(basename $<)) \\\
\	  $(addsuffix .dvi,$(basename $<))' > Makefile.Graphs
echo '    <h2>Contents of directory '$(echo $1 | sed 's|./||1')'</h2>
    <table border="1">
      <tr><td>Process</td><td>Graphs</td></tr>' >> $2/index.html
for I in *.tex; do
  bn=`echo $I | cut -d'.' -f1`
  make -f Makefile.Graphs $bn.ps
  ps2pdf $bn.ps
  convert -trim $bn.ps $bn.png
  num=`ls -a ${bn}_fg.t* | wc -l`
  echo -n "      <tr><td>"$num" graphs for<br>"$bn"<br><a href="$1/$bn".ps>[ps]</a> <a href="$1/$bn".pdf>[pdf]</a></td><td>" >> $2/index.html
  for i in $(ls -rc $bn*.png); do
    echo -n "<img src=\""$1/$i"\">" >> $2/index.html
  done
  echo "</td></tr>" >> $2/index.html
done
echo '    </table>' >> $2/index.html
}
sub(){
cd $3
plot $1 $2
for i in *; do
  test -d $i && (sub $1/$i ../$2 $i)
done
}
help(){
echo "usage: "$(basename $0)" <path>" 
exit 1
}
test -z "$1" && help
test -d $1 || help
echo '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <body>' > $PWD/$1/index.html
(sub . . $1)
echo '  </body>
</html>' >> $PWD/$1/index.html
