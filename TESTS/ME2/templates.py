import pkg_resources
from string import Template

mg_template  = Template(pkg_resources.resource_string(__name__, "mg_template.cc"))
sh_template  = Template(pkg_resources.resource_string(__name__, "sh_template.py"))
run_template = Template(pkg_resources.resource_string(__name__, "run_template.dat"))
