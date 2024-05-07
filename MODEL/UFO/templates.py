import pkgutil
from string import Template

model_template = Template(pkgutil.get_data(__name__, "model_template.C").decode('utf-8'))

lorentz_calc_template = Template(pkgutil.get_data(__name__, "lorentz_calc_template.C").decode('utf-8'))

color_calc_template = Template(pkgutil.get_data(__name__, "color_calc_template.C").decode('utf-8'))

sconstruct_template = Template(pkgutil.get_data(__name__, "sconstruct_template").decode('utf-8'))

run_card_template = Template(pkgutil.get_data(__name__, "run_card_template").decode('utf-8'))
