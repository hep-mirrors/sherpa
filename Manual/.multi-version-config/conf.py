import sys
import os

os.chdir("../source")
print(os.path.abspath("."))
sys.path.insert(0, os.path.abspath("."))
from conf import *

# Add multiversion and other config.
extensions += ["sphinx_multiversion"]
templates_path = ["_templates"]
print(extensions)

# for some reason this gets evaluated multiple times, so we have to
# make sure that we do not double apend
extra_sidebars = ["singlemulti.html", "versioning.html"]
for extra in extra_sidebars:
    if extra not in html_sidebars["**"]:
        html_sidebars["**"].append(extra)

# Configure sphinx-multiversion
smv_tag_whitelist = r'^.*$'
smv_branch_whitelist = r"^.*$"
