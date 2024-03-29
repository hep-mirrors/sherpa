#!/usr/bin/env python

particle = { 
	'1':'pi0',
        '2':'pi+',
	'3':'eta',
	'4':'K0',
	'5':'K+',
	'6':'eta\'',
	'11':'rho0',
	'12':'rho+',
	'13':'omega',
	'14':'K*0',
	'15':'K*+',
	'16':'phi',
	'21':'D+',
	'22':'D',
	'23':'D_s',
	'24':'eta_c',
	'25':'B+',
	'26':'B',
	'27':'B_s',
	'28':'B_c',
	'29':'eta_b',
	'31':'D*+',
	'32':'D*',
	'33':'D*_s',
	'34':'J/Psi',
	'35':'B*+',
	'36':'B*',
	'37':'B*_s',
	'38':'B*_c',
	'39':'Upsilon1',
	'41':'f_2(1270)',
	'42':'f\'_2(1525)',
	'43':'psi_2S',
	'51':'p',
	'52':'n',
	'53':'Sigma+',
	'54':'Sigma0',
	'55':'Sigma-',
	'56':'Lambda',
	'57':'Xi',
	'58':'Xi-',
	'61':'Delta++',
	'62':'Delta+',
	'63':'Delta0',
	'64':'Delta-',
	'65':'Sigma*+',
	'66':'Sigma*0',
	'67':'Sigma*-',
	'68':'Xi*0',
	'69':'Xi*-',
	'70':'Omega-',
	}

values = {
	'1' :' 9.59    \u00B1 0.33   ',
	'2' :'17.04    \u00B1 0.25   ',
	'3' :' 0.956   \u00B1 0.049  ',
	'4' :' 2.027   \u00B1 0.025  ',
	'5' :' 2.319   \u00B1 0.079  ',
	'6' :' 0.152   \u00B1 0.03   ',
	'11':' 1.295   \u00B1 0.125  ',
	'12':' 2.4     \u00B1 0.43   ',
	'13':' 1.083   \u00B1 0.088  ',
	'14':' 0.761   \u00B1 0.032  ',
	'15':' 0.731   \u00B1 0.058  ',
	'16':' 0.097   \u00B1 0.007  ',
	'21':' 0.184   \u00B1 0.018  ',
	'22':' 0.473   \u00B1 0.026  ',
	'23':' 0.129   \u00B1 0.013  ',
	'24':'                  ',
	'25':'                  ', 
	'26':'                  ',
	'27':'                  ',
	'28':'                  ',
	'29':'                  ',
	'31':' 0.182   \u00B1 0.009  ',
	'32':'                  ',
	'33':' 0.096   \u00B1 0.046  ',
	'34':' 0.00544 \u00B1 0.00029',
	'35':'                  ',
	'36':'                  ',
	'37':'                  ',
	'38':'                  ',
	'39':'                  ',
	'41':' 0.168   \u00B1 0.021  ',
	'42':' 0.02    \u00B1 0.008  ',
	'43':' 0.00229 \u00B1 0.00041',
	'51':' 0.991   \u00B1 0.054  ',
	'52':'                  ',
	'53':' 0.099   \u00B1 0.015  ',
	'54':' 0.074   \u00B1 0.009  ',
	'55':' 0.083   \u00B1 0.011  ',
	'56':' 0.373   \u00B1 0.008  ',
	'57':'                  ',
	'58':' 0.0262  \u00B1 0.001  ',
	'61':' 0.088   \u00B1 0.034  ',
	'62':'                  ',
	'63':'                  ',
	'64':'                  ',
	'65':' 0.0471? \u00B1 0.0046 ',
	'66':'                  ',
	'67':' 0.0471? \u00B1 0.0046 ',
	'68':' 0.0058  \u00B1 0.001  ',
	'69':'                  ',
	'70':' 0.00125 \u00B1 0.00024',
	}

file1 = open("IntermediateHadrons_Multis.dat",'r')
list1 = []
list1 = file1.readlines()
list1.pop(0)

file2 = open("PrimordialHadrons_Multis.dat",'r')
list2 = []
list2 = file2.readlines()
list2.pop(0)

tmpL = []
for itr in particle:
	tmpL.append(int(itr))

tmpL.sort()

for itr in tmpL:
	print(particle[str(itr)].rjust(10), values[str(itr)], list1[int(itr)].split("  ")[1].rjust(10), list2[int(itr)].split("  ")[1].rjust(10))  	
	
	
	
file1.close()
file2.close()
