# This file is part of the UFO.
#
# This file contains definitions for functions that
# are extensions of the cmath library, and correspond
# either to functions that are in cmath, but inconvenient
# to access from there (e.g. z.conjugate()),
# or functions that are simply not defined.
#
#
__date__ = "22 July 2010"
__author__ = "claude.duhr@durham.ac.uk"
import cmath
from .object_library import all_functions, Function

#
# shortcuts for functions from cmath
#
complexconjugate = Function(name = 'complexconjugate',
arguments = ('z',),
expression = 'z.conjugate()')

re = Function(name = 're',
arguments = ('z',),
expression = 'z.real')

im = Function(name = 'im',
arguments = ('z',),
expression = 'z.imag')

# New functions (trigonometric)
sec = Function(name = 'sec',
arguments = ('z',),
expression = '1./cmath.cos(z)')

asec = Function(name = 'asec',
arguments = ('z',),
expression = 'cmath.acos(1./z)')

csc = Function(name = 'csc',
arguments = ('z',),
expression = '1./cmath.sin(z)')

acsc = Function(name = 'acsc',
arguments = ('z',),
expression = 'cmath.asin(1./z)')

#==============================================================================
# PION FORM FACTOR FUNCTIONS
#==============================================================================
#
# These functions implement various parametrizations of the pion 
# electromagnetic form factor F_pi(Q^2), which describes the charge 
# distribution in the pion.
#
# Reference: RadioMonteCarLow2 pion-formfactor repository
#           https://gitlab.com/radiomontecarlow2/monte-carlo-results
#

# Constants for pion form factor calculations (in GeV units)
M_RHO = 0.775     # Rho meson mass (primary resonance)
M_RHO_PRIME = 1.465  # Rho' meson mass
M_OMEGA = 0.782   # Omega meson mass
GAMMA_RHO = 0.149 # Rho decay width
GAMMA_RHO_PRIME = 0.400  # Rho' decay width

#
# 1. Vector Dominance Model (VDM) - Single Rho Dominance
#    Simple one-pole parametrization
#
pion_ff_vdm_single = Function(name = 'pion_ff_vdm_single',
arguments = ('Q2',),
expression = '''
M_RHO**2 / (M_RHO**2 - Q2 - 1j*M_RHO*GAMMA_RHO)
''')

#
# 2. Two-Rho Model with Coupling
#    Includes main rho and excited rho' contribution
#
pion_ff_vdm_tworho = Function(name = 'pion_ff_vdm_tworho',
arguments = ('Q2', 'g_rho', 'g_rho_prime'),
expression = '''
(g_rho * M_RHO**2 / (M_RHO**2 - Q2 - 1j*M_RHO*GAMMA_RHO) +
 g_rho_prime * M_RHO_PRIME**2 / (M_RHO_PRIME**2 - Q2 - 1j*M_RHO_PRIME*GAMMA_RHO_PRIME)) / (g_rho + g_rho_prime)
''')

#
# 3. Gounaris-Sakurai Form
#    Energy-dependent width and phase space factors
#    Provides better description near rho pole
#
pion_ff_gounaris_sakurai = Function(name = 'pion_ff_gounaris_sakurai',
arguments = ('Q2'),
expression = '(M_rho**2 / (M_rho**2 - Q2)) * \
(1.0 + (Gamma_rho / (3.0*cmath.pi*(M_rho**2))) * \
 (M_rho**2 * cmath.log(1.0 + 6.0*Q2/(M_rho**2)) - \
  Q2 * (33.0 - 26.0*Q2/M_rho**2) / (1.0 - 0.2*Q2/M_rho**2))) / \
(M_rho**2 - Q2 - 1j*M_rho*Gamma_rho)')

#
# 4. Dispersion Relation Based Form (Lomon parametrization)
#    Connects spacelike and timelike regions via analyticity
#
pion_ff_dispersion = Function(name = 'pion_ff_dispersion',
arguments = ('Q2', 'a0', 'a1', 'a2'),
expression = '''
(a0 + a1*Q2 + a2*Q2**2) / (1.0 + (Q2/M_RHO**2) + (Q2/M_RHO_PRIME**2))
''')

#
# 5. Regge-type asymptotic form
#    Matches perturbative QCD at high Q^2
#
pion_ff_regge = Function(name = 'pion_ff_regge',
arguments = ('Q2', 'Lambda_QCD'),
expression = '''
1.0 / (1.0 + Q2/(4.0*Lambda_QCD**2))**2 * \
(Lambda_QCD**2 / Q2) * cmath.log(1.0 + Q2/(Lambda_QCD**2))**2
''')

#
# 6. Analytic model combining vector dominance and pQCD
#    Smooth interpolation from resonances to asymptotic regime
#
pion_ff_analytic_combined = Function(name = 'pion_ff_analytic_combined',
arguments = ('Q2', 'M_rho', 'Gamma_rho', 'beta', 'Lambda_QCD'),
expression = '''
(M_rho**2 / (M_rho**2 - Q2 - 1j*M_rho*Gamma_rho)) * \
(1.0 / (1.0 + beta*Q2/(4.0*Lambda_QCD**2))**2)
''')

#
# 7. Polynomial form from pion electroproduction fits
#    Empirical parametrization from experimental data
#
pion_ff_poly = Function(name = 'pion_ff_poly',
arguments = ('Q2', 'c0', 'c1', 'c2', 'c3'),
expression = '''
(c0 + c1*Q2 + c2*Q2**2) / (1.0 + c3*Q2**3)
''')

#
# 8. Elastic form factor using rho resonance
#    Standard VMD with energy-dependent width
#
def pion_ff_elastic(Q2, M_rho=M_RHO, Gamma_rho_ref=GAMMA_RHO, Q2_ref=0.1):
    """
    Elastic pion form factor with running width.
    Uses rho resonance with Q2-dependent width.
    
    Args:
        Q2: Four-momentum transfer squared (GeV^2)
        M_rho: Rho meson mass (default: 0.775 GeV)
        Gamma_rho_ref: Reference decay width at Q2_ref
        Q2_ref: Reference momentum for width
    
    Returns:
        Complex form factor value
    """
    # Running width (simplified pipi phase space)
    Gamma_rho = Gamma_rho_ref * (Q2 / Q2_ref) if Q2 > 0 else Gamma_rho_ref
    
    # Breit-Wigner with mass and width
    return M_rho**2 / (M_rho**2 - Q2 - 1j*M_rho*Gamma_rho)

#
# 9. Transition form factor (for gamma* gamma* coupling)
#    Used in hadronic light-by-light scattering
#
pion_ff_transition = Function(name = 'pion_ff_transition',
arguments = ('Q1_sq', 'Q2_sq', 'alpha', 'beta'),
expression = '''
1.0 / ((1.0 + Q1_sq/(0.776**2))**(1.0 + alpha*Q1_sq) * 
       (1.0 + Q2_sq/(0.776**2))**(1.0 + beta*Q2_sq))
''')

#
# 10. Modern high-precision model combining lattice QCD + experiments
#     Includes connected and disconnected contributions
#
pion_ff_lattice_combined = Function(name = 'pion_ff_lattice_combined',
arguments = ('Q2', 'a_connected', 'a_disconnected', 'corr_factor'),
expression = '''
(a_connected + a_disconnected*corr_factor) / \
(1.0 + (Q2/M_RHO**2)*(1.0 + Q2/(2.0*M_RHO_PRIME**2)))
''')

#
# Utility functions for form factor calculations
#

# Calculate Q^2 from four-momentum
Q_squared = Function(name = 'Q_squared',
arguments = ('E_photon', 'theta_angle'),
expression = '''
4.0 * E_photon**2 * cmath.sin(theta_angle/2.0)**2
''')

# Pion radius from form factor (derivative at Q^2=0)
pion_charge_radius = Function(name = 'pion_charge_radius',
arguments = ('dF_dQ2_at_zero',),
expression = '''
cmath.sqrt(-6.0 * dF_dQ2_at_zero)
''')

# Normalization condition: F_pi(0) = 1
pion_ff_normalized = Function(name = 'pion_ff_normalized',
arguments = ('Q2', 'F_unnorm', 'F_at_zero'),
expression = '''
F_unnorm / F_at_zero
''')

# Phase space factor for two-pion production
pion_pipi_phase_space = Function(name = 'pion_pipi_phase_space',
arguments = ('s', 'm_pi'),
expression = '''
cmath.sqrt(1.0 - 4.0*m_pi**2/s) if s > 4.0*m_pi**2 else 0.0+0.0j
''')

# Partial width for rho -> pi pi
rho_pipi_width = Function(name = 'rho_pipi_width',
arguments = ('M_rho', 'm_pi', 'g_rho'),
expression = '''
(g_rho**2 * M_rho) / (24.0*cmath.pi) * \
(cmath.sqrt(1.0 - 4.0*m_pi**2/M_rho**2))**3
''')