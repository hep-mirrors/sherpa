from .calculators import calc_parameter


def is_external(parameter):
    return parameter.nature == 'external'


def is_internal(parameter):
    return parameter.nature == 'internal'


def eval_parameter(parameter):
    if is_external(parameter):
        return parameter.name
    return calc_parameter(parameter.value)
