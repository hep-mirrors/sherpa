import pkg_resources
from string import Template


model_template = Template(pkg_resources.resource_string(__name__, "model_template.C").decode("utf-8"))

lorentz_calc_template = Template(pkg_resources.resource_string(__name__, "lorentz_calc_template.C").decode("utf-8"))

color_calc_template = Template(pkg_resources.resource_string(__name__, "color_calc_template.C").decode("utf-8"))

sconstruct_template = Template(pkg_resources.resource_string(__name__, "sconstruct_template").decode("utf-8"))

run_card_template = Template(pkg_resources.resource_string(__name__, "run_card_template").decode("utf-8"))
