#!/bin/bash
test -z "$1" && exit 1
zipinfo $1 | awk 'function mapped(file){
  command = ("unzip -p '$1' "file" | head -n 1");
  (command | getline result); close(command);
  split(result,a); return a[1]!=a[2]; }
BEGIN{np=0}
{ if (match($NF,".map")>0 && match($NF,"__j")==0)
    if (!mapped($NF)) { print np": "$NF; np+=1; } }'
