from .calculators import calc_parameter
from .particle import format_particle
from .parameter import is_external, is_internal, eval_parameter
from .vertex import vertex_implementation, chunk_vertices
from string import Template
from pathlib import Path
from .lorentz_writer import _filter_lorentz
import pkgutil
from .message import progress
import time


def order_sort(order):
    if order == 'QCD':
        return (0, order)
    elif order == 'QED':
        return (1, order)
    else:
        return (2, order)


def _physical_particle(particle):
    if particle.GhostNumber != 0:
        return False
    false_values = [0, 0.0, 'False', False]
    for key, value in particle.__dict__.items():
        if 'gold' in key.lower():
            if value in false_values:
                return False
    return True


class ModelWriter:
    def __init__(self, name, ufo, opts):
        self._lorentz_files = []
        self._name = name
        self._ufo = ufo
        self._version = ufo.__version__
        self._opts = opts
        self._params = ''
        self._particles = ''
        self._implementation = ''
        self._calls = ''

    def write(self, path):
        self._path = Path(path) if not isinstance(path, Path) else path
        start = time.time()
        self._write_parameters(self._ufo.all_parameters)
        end = time.time()
        progress(f'Parameter calculation took {end-start} seconds')
        start = time.time()
        self._write_particles(self._ufo.all_particles)
        end = time.time()
        progress(f'Particle calculation took {end-start} seconds')
        start = time.time()
        self._write_couplings(self._ufo.all_couplings)
        end = time.time()
        progress(f'Coulings calculation took {end-start} seconds')
        self._write_vertices(self._ufo.all_vertices, self._ufo.all_orders)
        self._get_lorentz_files(self._ufo.all_lorentz)
        self._write_cmakelists()
        self._write_model()

    def write_run_card(self, filename):
        order_dict = ', '.join([f'{order.name}: Any'
                                for order in self._ufo.all_orders])
        order_statement = f'Order: {{{order_dict}}}'

        all_particles = [p.pdg_code for p in filter(_physical_particle,
                                                    self._ufo.all_particles)]
        all_particles = ', '.join([str(p) for p in all_particles])

        template = pkgutil.get_data(__name__, "Templates/Sherpa.yaml.in")
        template = template.decode('utf-8')
        template = Template(template)
        substitution = {
            'model_name': self._name,
            'order_statement': order_statement,
            'all_particles': all_particles,
            'param_card': f'param_{self._name}.dat',
        }

        with open(filename, 'w') as output:
            output.write(template.safe_substitute(substitution))

    def _write_external_param(self, parameter):
        fmt = 'double {0} = p_dataread->GetEntry<double>("{1}", {2}, 0.0);\n'
        fmt += 'p_constants->insert(std::make_pair(std::string("{0}"), {0}));\n'
        if len(parameter.lhacode) > 2:
            raise ValueError(f'Unknown lhacode format \"{parameter.lhacode}\"')
        lha_indices = ", ".join([str(code) for code in parameter.lhacode])
        return fmt.format(parameter.name, parameter.lhablock, lha_indices)

    def _write_internal_param(self, parameter):
        if parameter.type == 'real':
            fmt = 'double {0} = ToDouble({1});\n'
        else:
            fmt = 'Complex {0} = {1};\n'
        fmt += 'DEBUG_VAR({0});\n'
        value = calc_parameter(parameter.value)

        return fmt.format(parameter.name, value)

    def _write_parameters(self, parameters):
        for parameter in parameters:
            if is_external(parameter):
                self._params += self._write_external_param(parameter)
            elif is_internal(parameter):
                self._params += self._write_internal_param(parameter)
            else:
                raise ValueError(f'Unknown nature {parameter.nature}')

    def _write_particle_info(self, particle):
        return format_particle(particle)

    def _write_particle_mass(self, particle):
        mass = eval_parameter(particle.mass)
        width = eval_parameter(particle.width)

        statement_fmt = 'ATOOLS::Flavour({0}).{1}(ToDouble({2}));\n'
        statement = statement_fmt.format(particle.pdg_code,
                                         'SetWidth', width)
        statement += statement_fmt.format(particle.pdg_code,
                                          'SetMass', mass)
        statement += statement_fmt.format(particle.pdg_code,
                                          'SetHadMass', mass)
        return statement

    def _write_particles(self, particles):
        for particle in particles:
            kfcode = particle.pdg_code
            # Don't explicitly add antiparticles
            if kfcode < 0:
                continue
            self._particles += self._write_particle_info(particle)
            self._params += self._write_particle_mass(particle)

    def _write_couplings(self, couplings):
        fmt = 'p_complexconstants->insert(std::make_pair(std::string("{0}"), {1}));\n'
        fmt += 'DEBUG_VAR((*p_complexconstants)["{0}"]);\n'
        for i, coupling in enumerate(couplings):
            print(f"Writing couplings: {float(i)/len(couplings)*100:0.2f}% finished", end='\r')
            self._params += fmt.format(coupling.name, calc_parameter(coupling.value))

    def _write_vertices(self, vertices, orders):
        hierarchy = sorted([order.name for order in orders], key=order_sort)
        assert (hierarchy[0] == 'QCD' and hierarchy[1] == 'QED')
        func_name = 'void vertices_{}() {{\n'
        for idx, chunk in enumerate(chunk_vertices(vertices, 10)):
            self._calls += f'vertices_{idx}();\n'
            self._implementation += func_name.format(idx)
            self._implementation += vertex_implementation(chunk, hierarchy)

        # Fill the order key getter
        self._index_of_order_key = ""
        for idx, order in enumerate(hierarchy):
            self._index_of_order_key += f'\n      if (key == "{order}") return {idx};'

    def _get_lorentz_files(self, structs):
        structs = filter(lambda struct: _filter_lorentz(struct,
                                                        self._opts["nmax"]),
                         structs)
        for struct in structs:
            self._lorentz_files.append(f'{struct.name}.C')

    def _write_cmakelists(self):
        template = pkgutil.get_data(__name__, "Templates/CMakeLists.txt.in")
        template = template.decode('utf-8')
        template = Template(template)
        substitution = {
            'PROJECT_NAME': self._name,
            'VERSION': self._version,
            'LORENTZ_FILES': '\n  '.join(self._lorentz_files),
            'MODEL_FILES': 'Model.C',
            'SHERPA_ROOT_DIR': self._opts['root_dir'],
            'INSTALL_DIR': self._opts['install_dir'],
            'MODEL_FLAGS': self._opts['model_flags'],
            'LORENTZ_FLAGS': self._opts['lorentz_flags'],
            'COLOR_FILES': '\n  '.join(self._opts['color_files']),

        }
        with open(self._path / 'CMakeLists.txt', 'w') as output:
            output.write(template.safe_substitute(substitution))

    def _write_model(self):
        # TODO: Figure out how to load when used as a package
        template = pkgutil.get_data(__name__, "Templates/model_template.C")
        template = template.decode('utf-8')
        template = Template(template)
        substitution = {
            'model_name': self._name,
            'particle_init': self._particles,
            'param_init': self._params,
            'declarations': self._implementation,
            'calls': self._calls,
            'fill_lorentz_map': '',
            'index_of_order_key': self._index_of_order_key,
        }
        with open(self._path / 'Model.C', 'w') as output:
            output.write(template.safe_substitute(substitution))
