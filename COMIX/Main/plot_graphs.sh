#!/bin/bash
if ! test -d $1; then
  echo "directory "$1" does not exist"
  exit 1
fi 
cd $1
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
echo '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <body>
    <h2>Contents of '$1'</h2>
    <table border="1">
      <tr><td>Process</td><td>Graphs</td></tr>' > index.html
for I in *.tex; do
  bn=`echo $I | cut -d'.' -f1`
  make -f Makefile.Graphs $bn.ps
  convert -trim $bn.ps $bn.png
  echo -n "      <tr><td>"$bn"</td><td>" >> index.html
  for i in $(ls -rc $bn*.png); do
    echo -n "<img src=\""$i"\">" >> index.html
  done
  echo "</td></tr>" >> index.html
done
echo '    </table>
  </body>
</html>' >> index.html
