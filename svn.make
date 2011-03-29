
SVN_Info.C: SVN_Info.C.in
	@if ! which svn > /dev/null || ! test -d $(srcdir)/.svn; then \
		if test -f $(srcdir)/$@; then \
			cp $(srcdir)/$@ $@.tmp; \
		else \
			echo -e '#include "ATOOLS/Org/SVN_Info.H"\n' > $@.tmp; \
			echo 'static ATOOLS::SVN_Info initializer' >> $@.tmp; \
			echo '("$(SVNTAG)","nourlfound","noversionfound","X");' >> $@.tmp; \
		fi; \
	else \
		cur=$$(echo "/"$(SVNTAG) | sed -e's/[+]/[+]/g'); \
		url=$$(svn info $(srcdir) | awk '{ if ($$1=="URL:") { split($$2,a,"svn/sherpa/"); sub("'$$cur'","",a[2]); print a[2]; } }'); \
		echo -e '#include "ATOOLS/Org/SVN_Info.H"\n' > $@.tmp; \
		echo 'static ATOOLS::SVN_Info initializer' >> $@.tmp; \
		echo '("$(SVNTAG)","'$$url'","'$$(svnversion $(srcdir))'","X");' >> $@.tmp; \
	fi

	@if test -z $(NOMD5SUM); then \
		flall="$(addprefix $(srcdir)/,$(filter-out $@,$(SOURCES)) $(HEADERS))"; \
		for file in $$flall; do \
			if test ! -f $$file.in; then \
				fl="$$fl $$file"; \
			fi; \
		done; \
		sed -r -e's/".?"\);/"'$$(cat $$fl | $(MD5COMMAND))'");/g' -i $@.tmp; \
	fi;

	@if ! diff $@.tmp $@ > /dev/null 2>&1; then \
		mv $@.tmp $@; \
	else \
		rm $@.tmp; \
	fi;

.PHONY: SVN_Info.C.in
