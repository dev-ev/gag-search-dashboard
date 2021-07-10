import pandas as pd
import numpy as np
import csv
import json
import glob, os

MZ_RANGE = (220, 2000)
MAX_CHARGE = 8
SPECIFICITY_MZ_BORDER = 900

DBA_MASS = 129.15175
PROTON_MASS = 1.00728
SO3_MASS = 79.95682
######################### example ##########################
EXAMPLE_LIST = [(1,((458.0610,1),),(),(0,20),'dp2S1 DeltaHexA-GalNAcS'),
            (4, ((268.5053,2),(538.0178,1)), ((458.0610,1),(667.1696,1)),(0,20),'dp2S2 DeltaHexA-GalNAcS2')]
    

def contains_pept(name):
    """Checks if the saccharide name contains the peptide fragment,
        such as EES, GS, SA etc"""
    contains_pept = False
    for pept_stub_name in ('_E', '_G', '_S', '_A'):
        if (pept_stub_name in name) and ('_NRE' not in name):
            contains_pept = True
    return contains_pept

def find_dba_adduct_forms(mass, max_dba_adducts, sulf, max_z):
    all_adduct_forms = []
    if max_dba_adducts > 0:
        dba_adducts = list(range(1, max_dba_adducts + 1))
        dba_adduct_forms = [(((mass + d*DBA_MASS - (sulf - d)*PROTON_MASS)/(sulf - d)),
                             (sulf - d)) for d in dba_adducts]

        dba_adduct_forms_minusOne = []
        if sulf > 2:
            d = max_dba_adducts - 1
            dba_adducts_minusOne = list(range(1, d + 1))
            
            dba_adduct_forms_minusOne = [(((mass + d*DBA_MASS - (sulf - d - 1)*PROTON_MASS)/(sulf - d - 1)),
                                          (sulf - d - 1)) for d in dba_adducts_minusOne]

        all_adduct_forms = dba_adduct_forms + dba_adduct_forms_minusOne
        # Delay the filtering by m/z until later
##        all_adduct_forms = [x for x in all_adduct_forms if ((x[0] >= MZ_RANGE[0]) and (x[0] <= MZ_RANGE[1]))]
        all_adduct_forms = [x for x in all_adduct_forms if (x[1] <= max_z)]
        
    return all_adduct_forms

def find_intact_forms(mass, min_z, max_z):
    charges = list(range(min_z, max_z + 1))
    intact_forms = [((mass - z*PROTON_MASS)/z, z) for z in charges]
    # Delay the filtering by m/z until later
##    intact_forms = [x for x in intact_forms if ((x[0] >= MZ_RANGE[0]) and (x[0] <= MZ_RANGE[1]))]
    return intact_forms

def find_neutral_m_sulfate_loss(mass, sulf, min_z, max_z):
    '''The output contains the (neutral mass with sulfate loss, new min charge, new max charge, new sulfate number)'''
    if (sulf > 1) and (max_z > 1):
        max_sulfate_losses = int(sulf/2)
        masses_sulfate_loss = [((mass-s*SO3_MASS), max(1,int((sulf-s)/2)), min(max_z,(sulf-s)), sulf-s)
                               for s in range(1, max_sulfate_losses + 1)]
        return masses_sulfate_loss
    else:
        return []

def min_max_charges(sulf, hexa, neuac, max_allowed_charge, is_pept=False):

    max_z = min((sulf + hexa + neuac), max_allowed_charge)
    min_z = 1
    
    if is_pept:
        max_z += 1
    
    return (min_z, max_z)
        
if __name__ == "__main__":

    # Use the extension .tab for the input lists to distinguish them from other text files in the folder
    for listname in glob.glob('*.tab'):
        file_suffix = listname[:-4]
        mega_list = []
        with open(f'gag_long_report{file_suffix}.txt', 'w') as ff:
            with open(listname, 'r') as neutral_masses:
                reader = csv.reader(neutral_masses, delimiter='\t')
                for line in reader:
                    name, sulf, hexa, neuac, mass, minrt, maxrt = line
                    is_pept = contains_pept(name)
                    sulf = int(sulf)
                    hexa = int(hexa)
                    neuac = int(neuac)

                    mass = float(mass)
                    minrt = float(minrt)
                    maxrt = float(maxrt)
                    #The maximum charge can be based only on sulfates, or on sulf + hexa + neuac
                    min_z, max_z = min_max_charges(sulf, hexa, neuac, MAX_CHARGE, is_pept)
                    
                    intact_forms = find_intact_forms(mass, min_z, max_z)
                    
                    max_dba_adducts = max(0, (sulf - 1))
                    
                    dba_adduct_forms = find_dba_adduct_forms(mass, max_dba_adducts, sulf, max_z)

                    masses_sulfate_loss = find_neutral_m_sulfate_loss(mass, sulf, min_z, max_z)

                    ff.write(name+' neutral mass '+str(mass)+' with min charge '+str(min_z)+' and max charge '+
                             str(max_z)+' max DBA adducts '+str(max_dba_adducts)+'\n')
                    
                    sulfate_loss_forms = []
                    dba_adduct_and_sulf_loss = []
                    if len(masses_sulfate_loss) > 0:
                        for sulfate_loss_sp in masses_sulfate_loss:
                            mass, min_z, max_z, sulf = sulfate_loss_sp
                            sulfate_loss_forms = sulfate_loss_forms + find_intact_forms(mass, min_z, max_z)
                            max_dba_adducts = max(0, (sulf - 1))
                            dba_adduct_and_sulf_loss = (dba_adduct_and_sulf_loss +
                                                        find_dba_adduct_forms(mass, max_dba_adducts, sulf, max_z))


                    ff.write('Intact forms '+str(intact_forms)+'\nDBA adduct forms '+str(dba_adduct_forms))
                    ff.write('\nSulfate loss neutral masses '+str(masses_sulfate_loss)+'\nSulfate loss forms '+str(sulfate_loss_forms))
                    ff.write('\nForms with sulfate loss and DBA adducts '+str(dba_adduct_and_sulf_loss))
     
                
                    unique_mzs = intact_forms
                    non_unique_mzs = sulfate_loss_forms + dba_adduct_and_sulf_loss
                    for i in dba_adduct_forms:
                        if i[0] >= SPECIFICITY_MZ_BORDER:
                            unique_mzs.append(i)
                        else:
                            non_unique_mzs.append(i)
                    # Calculate the seniority before filtering for m/z to improve consistency
                    # so that the seniority is not changed by the (rather arbitrary) m/z range settings
                    seniority = len(unique_mzs)+len(non_unique_mzs)
                    # Filter the m/z by range
                    unique_mzs = [x for x in unique_mzs if ((x[0] >= MZ_RANGE[0])
                                                            and (x[0] <= MZ_RANGE[1]))]
                    non_unique_mzs = [x for x in non_unique_mzs if ((x[0] >= MZ_RANGE[0])
                                                                    and (x[0] <= MZ_RANGE[1]))]
                    mega_list.append((seniority, unique_mzs, non_unique_mzs, (minrt,maxrt), name))
                    ff.write('\nCombined masses: '+str(len(unique_mzs)+len(non_unique_mzs)))
                    ff.write('\nUnique identifying m/z values are '+str(unique_mzs))
                    ff.write('\nAdditional non-unique m/z values are '+str(non_unique_mzs))
                    ff.write('\n_____________________________________________________________________\n')
        
            

        with open(f'gag_library_{file_suffix}.json','w') as fff:
            json.dump(mega_list,fff)
