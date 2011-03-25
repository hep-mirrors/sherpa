
SVN_Info.C: SVN_Info.C.in
	@if ! which svn > /dev/null || ! test -d .svn; \
	then cp $(srcdir)/$@ $@.tmp; else \
	cur=$$(echo "/"$(SVNTAG) | sed -e's/[+]/[+]/g'); \
	url=$$(svn info | awk '{ if ($$1=="URL:") { split($$2,a, \
	  "svn/sherpa/"); sub("'$$cur'","",a[2]); print a[2]; } }'); \
	echo -e '#include "ATOOLS/Org/SVN_Info.H"\n' > $@.tmp; \
	echo 'static ATOOLS::SVN_Info initializer' >> $@.tmp; \
	echo '("$(SVNTAG)","'$$url'","'$$(svnversion)'","X");' >> $@.tmp; \
	fi; if test -z $(NOMD5SUM); then sed \
	  -r -e's/".?"\);/"'$$(cat $$(echo *.[CH] | sed 's/'$@'//g') | \
	  $(MD5COMMAND))'");/g' -i $@.tmp; fi; \
	if ! diff $@.tmp $@ > /dev/null 2>&1; \
	  then mv $@.tmp $@; else rm $@.tmp; fi;

.PHONY: SVN_Info.C.in
