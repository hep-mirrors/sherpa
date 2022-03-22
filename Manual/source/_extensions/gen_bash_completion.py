
# -*- coding: utf-8 -*-
"""
A sphinx builder to generate the keyword lists for the sherpa bash compledion.
"""

import codecs
import re
from os import path

from docutils.io import StringOutput
from docutils import nodes
from docutils.core import Publisher

from sphinx.builders import Builder
from sphinx.util.osutil import ensuredir, os_path
from sphinx.environment.adapters.indexentries import IndexEntries

class CompletionBuilder(Builder):
    name = 'completion'
    format = 'text'
    out_suffix = '.index'
    allow_parallel = True
    blacklist = [
        "'PARAMETER: Value'",
        "'Tags: {TAG: Value}'",
        "command line option"
    ]
    brack_re = re.compile('\s*<.*?>\s*')

    def init(self):
        pass

    def get_target_uri(self, docname, typ=None):
        return docname + self.out_suffix

    def get_toctree(self, docname, collapse=True, **kwds):
        if 'includehidden' not in kwds:
            kwds['includehidden'] = False
        toctree = TocTree(self.env).get_toctree_for(docname, self, collapse, **kwds)
        return self.render_partial(toctree)['fragment']

    def get_outdated_docs(self):
        for docname in self.env.found_docs:
            if docname not in self.env.all_docs:
                yield docname
                continue
            targetname = path.join(self.outdir,self.env.doc2path(docname),
                                   self.out_suffix)
            try:
                targetmtime = path.getmtime(targetname)
            except Exception:
                targetmtime = 0
            try:
                srcmtime = path.getmtime(self.env.doc2path(docname))
                if srcmtime > targetmtime:
                    yield docname
            except EnvironmentError:
                # source doesn't exist anymore
                pass

    def prepare_writing(self, docnames):
        pass

    def write_doc(self, docname, doctree):
        pass

    def gen_indices(self):
        genindex = IndexEntries(self.env).create_index(self)
        symbols = []

        for key, entries in genindex:
                    for column in entries:
                        symbols.append(column[0])

        items = [re.sub(self.brack_re, '', sym) \
                 for sym in symbols if sym not in self.blacklist]
        symbols = []
        flags = []

        for sym in items:
            if sym[0] == '-':
                flags.append(sym)
            else:
                symbols.append(sym)

        try:
            completions = path.join(self.outdir, os_path('completion') \
                                + self.out_suffix)
            options = path.join(self.outdir, os_path('options') \
                                + self.out_suffix)

            with codecs.open(completions, 'w', 'utf-8') as f:
                for sym in symbols:
                    f.write('* ' + sym + '\n')

            with codecs.open(options, 'w', 'utf-8') as f:
                f.write(" ".join(flags))

        except (IOError, OSError) as err:
            self.warn("error writing file %s: %s" % (outfilename, err))

    def finish(self):
        self.finish_tasks.add_task(self.gen_indices)


def setup(app):
    app.add_builder(CompletionBuilder)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
