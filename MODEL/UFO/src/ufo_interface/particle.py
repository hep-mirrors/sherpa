PARTICLE_FMT = 'ATOOLS::s_kftable[{0}] = new ATOOLS::Particle_Info('
PARTICLE_FMT += '{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}'


def _is_majorana(particle):
    if particle.spin != 2:
        return False
    return True if particle.name == particle.antiname else False


def _is_self_conjugate(particle):
    if particle.name != particle.antiname:
        return 0
    if _is_majorana(particle):
        return 1
    return -1


def quoted(string):
    return '"' + string + '"'


def tex(string):
    return string.replace(r'\ '.rstrip(), r'\\')


def format_particle(particle):
    if hasattr(particle, 'hadron') and particle.hadron == 1:
        fmt = PARTICLE_FMT + ');\n'
        return fmt.format(particle.pdg_code, 1000.0, 0.0, 0.0,
                          int(3*particle.charge), particle.spin-1,
                          1, 1, quoted(particle.name),
                          quoted(tex(particle.texname)))
    else:
        fmt = PARTICLE_FMT + ', {10}, {11}, {12}, {13}, {14});\n'
        color = 0 if particle.color == 1 else particle.color
        self_conjugate = _is_self_conjugate(particle)
        massive = 0 if particle.mass == 0.0 else 1
        return fmt.format(particle.pdg_code, 1000.0, 0.0, 0.0,
                          int(3*particle.charge), color,
                          particle.spin-1, self_conjugate, 1, 0,
                          massive, quoted(particle.name),
                          quoted(particle.antiname),
                          quoted(tex(particle.texname)),
                          quoted(tex(particle.antitexname)))
