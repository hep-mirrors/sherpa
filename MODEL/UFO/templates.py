import pkg_resources
from string import Template

model_template = Template(pkg_resources.resource_string(__name__, "model_template.C"))

lorentz_calc_template = Template(pkg_resources.resource_string(__name__, "lorentz_calc_template.C"))

color_calc_template = Template(pkg_resources.resource_string(__name__, "color_calc_template.C"))

sconstruct_template = Template(pkg_resources.resource_string(__name__, "sconstruct_template"))

run_card_template = Template(pkg_resources.resource_string(__name__, "run_card_template"))
