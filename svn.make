
SVN_Info.C: SVN_Info.C.in
	@if ! which svn > /dev/null; then exit 0; fi; \
	cur=$$(echo "/"$(SVNTAG) | sed -e's/__/\//g' -e's/[+]/[+]/g'); \
	url=$$(svn info | awk '{ if ($$1=="URL:") { split($$2,a,"svn/sherpa/"); \
	  sub("'$$cur'","",a[2]); print a[2]; } }'); \
	echo -e '#include "ATOOLS/Org/SVN_Info.H"\n' > $(basename $<).tmp; \
	echo 'static ATOOLS::SVN_Info initializer' >> $(basename $<).tmp; \
	echo '("$(SVNTAG)","'$$url'","'$$(svnversion)'");' >> $(basename $<).tmp; \
	if ! diff $(basename $<).tmp $(basename $<) > /dev/null 2>&1; \
	  then mv $(basename $<).tmp $(basename $<); \
	  else rm $(basename $<).tmp; fi;

.PHONY: SVN_Info.C.in
