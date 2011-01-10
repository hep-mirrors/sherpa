
SVN_Info.C: SVN_Info.C.in
	@if ! which svn > /dev/null; then exit 0; fi; \
	cur=$$(echo "/"$(SVNTAG) | sed -e's/_/\//g'); \
	url=$$(svn info | awk '{ if ($$1=="URL:") { \
	  split($$2,a,"svn/sherpa/"); sub("'$$cur'","",a[2]); \
	  gsub("/","\\/",a[2]); print a[2]; } }'); \
	sed -e's/SVNBRANCH/'$$url'/1' -e's/SVNREVISION/'$$(svnversion)'/1' \
	  -e's/SVNTAG/$(SVNTAG)/' < $< > $(basename $<).tmp; \
	if ! diff $(basename $<).tmp $(basename $<) > /dev/null 2>&1; \
	  then mv $(basename $<).tmp $(basename $<); \
	  else rm $(basename $<).tmp; fi;

.PHONY: SVN_Info.C.in
