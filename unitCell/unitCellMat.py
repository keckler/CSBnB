pi = 3.14159265359
AvNum = 6.0221409E23

pinDiameter = 1.2217 #cm
pinPitch = 1.3609 #cm
numRings = 8 #includes central "ring"
cladThickness = 0.6109E-1 #cm
bundleLength = 30. #includes end plates, cm
wirePitch = bundleLength/6 #cm
endPlateThickness = 0.1 #axial thickness of end plate, cm
endPlateRingThickness = 0.5 #thickness of ring connecting pins together, cm
ductThickness = 0.35 #cm
ventLength = 1. #cm
smearDensity = 0.75
axialExpansionLength = 0.04 * bundleLength #cm
enrichment = 0.0072 #from wikipedia
fuelTemp = 600. + 273.15 #K, taken from metal ABR studies
coolTemp = 550. + 273.15 #K, taken from metal ABR studies

#material mass densities
#HT9 taken from 'barrel' function in x_materialcalc.py of ADOPT
#vent taken as SS316 density from 'barre' function in x_materialcalc.py of ADOPT
#fuel taken from section 2.2.1.1 of IAEA's Thermophysical Properties of Materials for Nuclear Engineering 
#Na taken from 'CoolantInletProperties' function in x_coolant.py of ADOPT
mat2dens = {'HT9':7.874 - 3.23E-4*coolTemp,
            'fuel':(19.36E3 - 1.03347*fuelTemp) *1000/100/100/100 * smearDensity,
            'vent':8.804 - 4.209E-4*coolTemp - 3.894E-8*coolTemp**2,
            'Na':(1014. - 0.235*coolTemp)*1000/100/100/100}

#material of each component
component2mat = {'fuel':'fuel',
                 'clad':'HT9',
                 'vent':'vent',
                 'wire':'HT9',
                 'duct':'HT9',
                 'endP':'HT9',
                 'cool':'Na'}

#WEIGHT fractions of each element in core materials
#HT9 taken from ADOPT built-in properties, which are taken from Leibowitz and Blumquist
#vent assumed to be entirely made of iron, as no design info is available
mat2el2frac = {'HT9':{'28':0.005,
                      '24':0.120,
                      '25':0.002,
                      '14':0.0025,
                      '74':0.005,
                      '23':0.005,
                      '6':0.002,
                      '26':0.8485},
               'fuel':{'92':1.0},
               'vent':{'26':1.0},
               'Na':{'11':1.0}}

#isotopic fractions of each element, taken from wikipedia
isoFrac = {'6':{'012':1.0},
           '11':{'023':1.0},
           '14':{'028':0.92223, '029':0.04685, '030':0.03092},
           '23':{'051':0.99750},
           '24':{'050':0.04345, '052':0.83789, '053':0.09501, '054':0.02365},
           '25':{'055':1.0},
           '26':{'054':0.05845, '056':0.91754, '057':0.02119, '058':0.00282},
           '28':{'058':0.680769, '060':0.262231, '061':0.011399, '062':0.036345, '064':0.009256},
           '74':{'182':0.2650, '183':0.1431, '184':0.3064, '186':0.2843},
           '92':{'235':enrichment, '238':1-enrichment}}

#calculate number of pins in assembly based on numRings
numPins = 1
for i in range(2,numRings+1):
    numPins += i*6

#calculate relevant areas
wireDiameter = pinPitch - pinDiameter
wireLength = ((bundleLength - ventLength - 2*endPlateThickness)/wirePitch) * (pinDiameter + wireDiameter) * pi + (bundleLength - ventLength - 2*endPlateThickness)
wireArea = pi * (wireDiameter/2)**2
pinArea = pi * (pinDiameter/2)**2
cladArea = pi * ((pinDiameter/2)**2 - ((pinDiameter/2)-cladThickness)**2)
pinInnerArea = pinArea - cladArea
outerHexFTF = 2.*((numRings-1) * 3.**(1/2.)/2.*(pinDiameter+wireDiameter) + (pinDiameter/2) + wireDiameter + ductThickness)
outerHexArea = 3.**(1/2.) / 2. * outerHexFTF**2.
innerHexFTF = 2.*((numRings-1) * 3.**(1/2.)/2.*(pinDiameter+wireDiameter) + (pinDiameter/2) + wireDiameter)
innerHexArea = 3.**(1/2.) / 2. * innerHexFTF**2.

#calculate volumes of each component based on areas calculated above
outerHexVol = outerHexArea * bundleLength #volume of entire assembly hexagonal region
fuelVol = numPins * pinInnerArea * (bundleLength - ventLength - axialExpansionLength - 2*endPlateThickness - 1*cladThickness)
cladVol = numPins * cladArea * (bundleLength - ventLength - 2*endPlateThickness - 1*cladThickness) + numPins * pi * (pinDiameter/2)**2 * cladThickness
ventVol = pinArea * ventLength
wireVol = wireLength * wireArea
ductVol = (outerHexArea - innerHexArea) * bundleLength
endPlateVol = 2 * (outerHexArea - innerHexArea) * endPlateThickness #JUST AN APPROXIMATION, SHOULD BE CHANGED
coolantVol = outerHexVol - fuelVol - cladVol - ventVol - wireVol - ductVol - endPlateVol

#calculate component volume fractions
componentVolFrac = {'fuel':fuelVol/outerHexVol,
                    'clad':cladVol/outerHexVol,
                    'vent':ventVol/outerHexVol,
                    'wire':wireVol/outerHexVol,
                    'duct':ductVol/outerHexVol,
                    'endP':endPlateVol/outerHexVol,
                    'cool':coolantVol/outerHexVol}

#calculate volume fracs of each unique material
matVolFrac = {}
for component in component2mat.keys():
    mat = component2mat[component]
    if mat in matVolFrac.keys():
        matVolFrac[mat] += componentVolFrac[component]
    else:
        matVolFrac[mat] = 0
        matVolFrac[mat] += componentVolFrac[component]

#calculate number densities of each isotope
isotope2numDens = {}
for mat in matVolFrac.keys():
    for element in mat2el2frac[mat].keys():
        for isotope in isoFrac[element].keys():
            ZAI = element+isotope
            if ZAI in isotope2numDens.keys():
                isotope2numDens[ZAI] += mat2el2frac[mat][element] * isoFrac[element][isotope] * mat2dens[mat] * AvNum * matVolFrac[mat] / float(isotope) / 1E24 #1/barn/cm
            else:
                isotope2numDens[ZAI] = mat2el2frac[mat][element] * isoFrac[element][isotope] * mat2dens[mat] * AvNum * matVolFrac[mat] / float(isotope) / 1E24 #1/barn/cm

fm = open('CSBnB_mburn','w')
fm.write('mat unitCellMat sum burn 1\n')
for isotope in isotope2numDens.keys():
    fm.write('    '+isotope+'.06c '+str(isotope2numDens[isotope])+'\n')