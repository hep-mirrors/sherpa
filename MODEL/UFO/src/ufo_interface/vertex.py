from itertools import groupby
from .calculators import calc_color


def _filter_vertex(vertex):
    for particle in vertex.particles:
        gold_attr = [att for att in particle.__dict__
                     if att.lower().startswith('gold')]
        if len(gold_attr) == 0:
            return False
        falses = [0, 0.0, False]
        if any([particle.__dict__[attr] not in falses for attr in gold_attr]):
            return False
        if hasattr(particle, 'GhostNumber') and particle.GhostNumber != 0:
            return False
        if particle.spin < 0:
            return False
    return True


def chunk_vertices(vertices, n):
    filtered_vertices = list(filter(_filter_vertex, vertices))
    for i in range(0, len(filtered_vertices), n):
        yield filtered_vertices[i:i + n]


def _order_func(value):
    return value[1].order


def _get_orders(coupling, hierachy):
    return [coupling.order.get(name, 0) for name in hierachy]


def vertex_implementation(vertex_list, hierachy):
    indent = '\n  '

    func_impl = ""
    for vert in vertex_list:
        coupling_types = [list(group)
                          for _, group in groupby(vert.couplings.items(),
                                                  key=_order_func)]
        for coupling_type in coupling_types:
            orders = _get_orders(coupling_type[0][1], hierachy)
            func_impl += indent + 'm_v.push_back(Single_Vertex());'
            for part in vert.particles:
                func_impl += (indent + 'm_v.back().AddParticle(ATOOLS::Flavour((kf_code){0}, {1}));'
                              .format(abs(part.pdg_code), 0 if part.pdg_code > 0 else 1))
            for coupling_info in coupling_type:
                color_idx, lorentz_idx = coupling_info[0]
                lorentz = vert.lorentz[lorentz_idx]
                color = calc_color(vert.color[color_idx])
                coupling = coupling_info[1]
                func_impl += (indent +
                              'm_v.back().cpl.push_back(ATOOLS::Kabbala("{0}",'
                              + ' ComplexConstant(std::string("{0}"))));').format(coupling.name)
                func_impl += (indent +
                              'm_v.back().Color.push_back({0});').format(color)
                func_impl += (indent +
                              f'm_v.back().Lorentz.push_back("{lorentz}");')

            func_impl += (indent +
                          f'm_v.back().order.resize({len(hierachy)});')
            for i, order in enumerate(orders):
                func_impl += (indent +
                              f'm_v.back().order[{i}] = {order};')
            func_impl += '\n'

    return func_impl[1:] + '}\n'
