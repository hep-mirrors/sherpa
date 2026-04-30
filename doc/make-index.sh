#!/bin/bash

# this is supposed to run within the CI's "manual-index" job to build the index.html

releases=$(git ls-remote --tags | cut -d '^' -f 1 | cut -d "/" -f 3- | sort -r | uniq); 

cat <<EOF
<!DOCTYPE html>

<html lang="en" data-content_root="./">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" /><meta name="viewport" content="width=device-width, initial-scale=1" />

    <title>Sherpa User Manuals</title>
    <link rel="stylesheet" type="text/css" href="../_static/pygments.css?v=6625fa76" />
    <link rel="stylesheet" type="text/css" href="../_static/alabaster.css?v=7b5ee873" />
    <link rel="stylesheet" type="text/css" href="../_static/css/custom.css?v=ca636962" />
    <script src="../_static/documentation_options.js?v=5929fcd5"></script>
    <script src="../_static/doctools.js?v=9bcbadda"></script>
    <script src="../_static/sphinx_highlight.js?v=dc90522c"></script>
    <link rel="icon" href="../_static/favicon.ico"/>
    <link rel="index" title="Index" href="../genindex.html" />
    <link rel="search" title="Search" href="../search.html" />
    <link rel="next" title="About" href="../monte-carlo.html" />
    <link rel="stylesheet" href="../_static/custom.css" type="text/css" />
    
    
    <meta name="viewport" content="width=device-width, initial-scale=0.9, maximum-scale=0.9" />

  </head>
  
  <body>
    
    <div class="document">
      <div class="documentwrapper">
        <div class="bodywrapper">
          <div class="body" role="main">
            <section id="sherpa-user-manuals">
              <span id="index"></span><h1>Sherpa User Manuals<a class="headerlink" href="#sherpa-user-manuals" title="Link to this heading">¶</a></h1>
              <h2>Released versions</h2>
              <ul class="simple">
EOF

for dir in $releases; do
    echo "                <li><a class='reference external' href='${dir}/index.html'>${dir}</a></li>";
done;

cat <<EOF
              </ul>

              <ul class="simple">
                <li><a class="reference external" href="master/index.html">Latest development branch</a></li>
              </ul>

              <div class="toctree-wrapper compound">
              </div>
            </section>
          </div>
        </div>
      </div>
      
      <div class="sphinxsidebar" role="navigation" aria-label="Main">
        <div class="sphinxsidebarwrapper">
          <p class="logo" style="display: block !important;">
            <a href="../">
              <img class="logo" src="../_static/images/sherpa-logo.png" alt="Logo"/>
            </a>
          </p>
          
          <h3>Navigation</h3>
          <ul>
            <li class="toctree-l1"><a class="reference internal" href="../monte-carlo.html">About</a></li>
            <li class="toctree-l1"><a class="reference internal" href="../changelog.html">Downloads</a></li>
            <li class="toctree-l1"><a class="reference internal" href="../team.html">Sherpa Team</a></li>
            <li class="toctree-l1"><a class="reference internal" href="../publications.html">Publications</a></li>
            <li class="toctree-l1"><a class="reference internal" href="../theses.html">Theses</a></li>
          </ul>
          <hr />
          <ul>
            <li class="toctree-l1"><a href="../sherpa/">User Manual</a></li>
            <li class="toctree-l1"><a href="https://gitlab.com/sherpa-team/sherpa/issues/">Issue Tracker</a></li>
            <li class="toctree-l1"><a href="https://gitlab.com/sherpa-team/sherpa/">Git Repo</a></li>
            <li class="toctree-l1"><a href="https://www.hepforge.org/lists/listinfo/sherpa-announce">Mailing List</a></li>
          </ul>
          
        </div>
      </div>
      <div class="clearer"></div>
    </div>
  </body>
</html>
EOF
