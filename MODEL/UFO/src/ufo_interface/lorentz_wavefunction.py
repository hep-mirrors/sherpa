import sympy
from sympy.tensor.array.expressions import ArraySymbol
from sympy.tensor.array.expressions.array_expressions import ArrayElement
from .lorentz_algebra import METRIC, get_repeated_indices
from .lorentz_algebra import convert_indexed_to_array

WAVEFUNCTIONS = {3: ['j{}', (4,)],
                 -2: ['j{}', (4,)],
                 2: ['j{}', (4,)],
                 -4: ['j{}', (4, 4,)],
                 4: ['j{}', (4, 4,)]}


WAVEFUNCTION_INDEX = {3: ['dummy_L'],
                      -2: ['S'],
                      2: ['S'],
                      -4: ['dummy_L', 'S'],
                      4: ['dummy_L', 'S']}


def contract_wavefunctions(lorentz, spins):
    wavefunctions = []
    is_fermion = True
    for i, spin in enumerate(spins):
        if spin == 1:
            wavefunctions.append(sympy.Symbol(f'j{i}0'))
            continue
        # Assume that second fermion is the bared wavefunction
        sign = 1 if spin % 2 != 0 or is_fermion else -1
        is_fermion = not is_fermion if spin % 2 == 0 else is_fermion
        indices = [sympy.Symbol(f'{label}_{i+1}')
                   for label in WAVEFUNCTION_INDEX[sign*spin]]
        name = WAVEFUNCTIONS[sign*spin][0].format(i)
        shape = WAVEFUNCTIONS[sign*spin][1]
        wavefunction = ArraySymbol(name, shape)
        wavefunctions.append(wavefunction[tuple(indices)])
    rotations = []
    num_ext = len(wavefunctions)
    for i in range(num_ext):
        subwavefunctions = wavefunctions[:i]+wavefunctions[i+1:]
        rotation = {j: (j+(num_ext-1-i)) % num_ext for j in range(num_ext)}
        if hasattr(lorentz, 'args'):
            local = lorentz.copy()
        else:
            local = lorentz
        for subwavefunction in subwavefunctions:
            # Add dummy lorentz index with metric
            if not isinstance(subwavefunction, ArrayElement):
                local *= subwavefunction
                continue
            dummy_indices = [sympy.Symbol(str(index)[6:])
                             for index in subwavefunction.indices
                             if 'L' in str(index)]
            for dummy_index in dummy_indices:
                local *= METRIC[dummy_index,
                                sympy.Symbol(f"dummy_{dummy_index}")]
            local *= subwavefunction
        repeated = get_repeated_indices(local)
        for index in repeated:
            local = sympy.Sum(local, (index, 0, 3))
        local = convert_indexed_to_array(local)

        rotations.append([rotation, local])

    return rotations
