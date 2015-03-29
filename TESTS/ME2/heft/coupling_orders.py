# This file was automatically created by FeynRules 2.0.8
# Mathematica version: 8.0 for Linux x86 (64-bit) (February 23, 2011)
# Date: Fri 14 Nov 2014 13:52:48


from object_library import all_orders, CouplingOrder


QCD = CouplingOrder(name = 'QCD',
                    expansion_order = 99,
                    hierarchy = 1)

QED = CouplingOrder(name = 'QED',
                    expansion_order = 99,
                    hierarchy = 2)

HIG = CouplingOrder(name = 'HIG',
                    expansion_order = 1,
                    hierarchy = 3)

HIW = CouplingOrder(name = 'HIW',
                    expansion_order = 1,
                    hierarchy = 4)
