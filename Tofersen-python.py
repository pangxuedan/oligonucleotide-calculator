
def create_oligonucleotide_dictionary():
    """
    Create a dictionary with oligonucleotide names and their corresponding molecular weights.
    """
    oligo_dict = {
        # MOE nucleotides
        "MOE G": 419.06646,
        "MOE G P=O": 403.0893,
        "MOE G nucleoside": 341.13353,
        "MOE A": 403.0715,
        "MOE A P=O": 387.09438,
        "MOE A nucleoside": 325.13861,
        "MOE MeU": 394.05997,
        "MOE MeU P=O": 378.08282,
        "MOE MeU nucleoside": 316.12705,
        "MOE MeC": 393.07596,
        "MOE MeC P=O": 377.0988,
        "MOE MeC nucleoside": 315.14304,
        
        # DNA nucleotides
        "dG": 345.02968,
        "dG P=O": 239.05252,
        "dG nucleoside": 267.09675,
        "dA": 329.03476,
        "dA P=O": 313.05761,
        "dA nucleoside": 251.10184,
        "dT": 320.02319,
        "dT P=O": 304.04604,
        "dT nucleoside": 242.09027,
        "dMeC": 319.03918,
        "dMeC P=O": 303.06202,
        "dMeC nucleoside": 241.10626,
        "dC P=O": 289.04637,
        "dC nucleoside": 227.09061,
        
        # NMA nucleotides
        "NMA G": 432.0617,
        "NMA A": 416.06679,
        "NMA G P=O": 416.08455,
        "NMA G nucleoside": 354.12878,
        "NMA MeU": 407.05522,
        "NMA A P=O": 400.08963,
        "NMA A nucleoside": 338.13387,
        "NMA MeC": 406.07121,
        "NMA MeU P=O": 391.07807,
        "NMA MeU nucleoside": 329.1223,
        "NMA MeC P=O": 390.09405,
        "NMA MeC nucleoside": 328.13828,
        
        # m-nucleotides
        "mG": 375.05024,
        "mA": 359.04533,
        "mG P=O": 359.06308,
        "mG nucleoside": 297.10732,
        "mU": 336.01811,
        "mA P=O": 343.06817,
        "mA nucleoside": 281.1124,
        "mC": 335.03409,
        "mU P=O": 320.04095,
        "mU nucleoside": 258.08519,
        "mC P=O": 319.05694,
        "mC nucleoside": 257.10117,
        
        # cET nucleotides
        "cEt G": 387.04024,
        "cEt A": 371.04533,
        "cEt G P=O": 371.06308,
        "cEt G nucleoside": 309.10732,
        "cEt MeU": 362.03376,
        "cEt A P=O": 355.06817,
        "cEt A nucleoside": 293.1124,
        "cEt MeC": 361.04974,
        "cEt MeU P=O": 346.0566,
        "cEt MeU nucleoside": 284.10084,
        "cEt MeC P=O": 345.07259,
        "cEt MeC nucleoside": 283.11682,
        
        # f-nucleotides
        "fG": 363.02025,
        "fA": 347.02534,
        "fG P=O": 347.0431,
        "fG nucleoside": 285.08733,
        "fU": 323.99812,
        "fA P=O": 331.04818,
        "fA nucleoside": 269.09242,
        "fC": 323.01411,
        "fU P=O": 308.02097,
        "fU nucleoside": 246.0652,
        "fC P=O": 307.03695,
        "fC nucleoside": 245.08118,
        
        # Other components
        "AH": 195.04829,
        "AH P=O": 179.07113,
        "GalNAc": 1518.782075
    }
    
    return oligo_dict

def calculate_molecular_weight(sequence):
    """
    Calculate the molecular weight of an oligonucleotide sequence by simple addition.
    
    Parameters:
    sequence (str): A string representing the oligonucleotide sequence
                    Format example: "MOE G-MOE A"
    
    Returns:
    tuple: (total_weight, component_weights, found_all, missing_components)
           total_weight - calculated molecular weight (simple sum of components)
           component_weights - list of (component, weight) tuples
           found_all - boolean indicating if all components were found
           missing_components - list of components not found in the dictionary
    """
    # Get the oligonucleotide dictionary
    oligo_dict = create_oligonucleotide_dictionary()
    
    # Split the sequence by hyphens
    components = [comp.strip() for comp in sequence.split('-')]
    
    # Calculate the total molecular weight by simple addition
    total_weight = 0
    component_weights = []
    found_all = True
    missing_components = []
    
    for component in components:
        if component in oligo_dict:
            weight = oligo_dict[component]
            total_weight += weight
            component_weights.append((component, weight))
        else:
            found_all = False
            missing_components.append(component)
    
    # REMOVED: Water loss correction - now using simple addition only
    # The total_weight is now just the sum of all recognized components
    
    return total_weight, component_weights, found_all, missing_components

def check_po_impurity(components):
    """
    Check if the sequence contains nucleotides that can have P=O impurities.
    
    Parameters:
    components (list): List of component names in the sequence
    
    Returns:
    bool: True if P=O impurity should be calculated, False otherwise
    """
    # Define nucleotides that can have P=O impurities
    po_eligible_nucleotides = [
        "MOE G", "MOE A", "MOE MeU", "MOE MeC",
        "dG", "dA", "dT", "dMeC", 
        "NMA G", "NMA A", "NMA MeU", "NMA MeC",
        "mG", "mA", "mU", "mC",
        "cEt G", "cEt A", "cEt MeU", "cEt MeC",
        "fG", "fA", "fU", "fC"
    ]
    
    # Check if any component in the sequence is eligible for P=O impurity
    for component in components:
        if component in po_eligible_nucleotides:
            return True
    
    return False

def calculate_impurities(sequence, average_mw):
    """
    Calculate the N-1 and N+1 impurities for a given oligonucleotide sequence.
    N+1 impurities: average_mw + component molecular weight
    N-1 impurities: average_mw - component molecular weight
    Duplicated components are counted only once in impurities.
    
    Parameters:
    sequence (str): A string representing the oligonucleotide sequence
                    Format example: "MOE G-MOE A"
    average_mw (float): The average molecular weight from user input
    
    Returns:
    tuple: (n_minus_1, n_plus_1, mobile_phase_impurities)
           n_minus_1 - list of (component, impurity_weight) tuples
           n_plus_1 - list of (component, impurity_weight) tuples
           mobile_phase_impurities - list of mobile phase related impurities
    """
    # Get the oligonucleotide dictionary
    oligo_dict = create_oligonucleotide_dictionary()
    
    # Split the sequence by hyphens
    components = [comp.strip() for comp in sequence.split('-')]
    
    # Check if all components are found
    all_found = all(comp in oligo_dict for comp in components)
    if not all_found:
        return [], [], []
    
    # Get unique components for impurity calculations
    unique_components = []
    for component in components:
        if component not in unique_components and component in oligo_dict:
            unique_components.append(component)
    
    # Calculate N-1 impurities (average_mw - component molecular weight)
    n_minus_1 = []
    for component in unique_components:
        n_minus_1_weight = average_mw - oligo_dict[component]
        n_minus_1.append((component, n_minus_1_weight))
    
    # Calculate N+1 impurities (average_mw + component molecular weight)
    n_plus_1 = []
    for component in unique_components:
        n_plus_1_weight = average_mw + oligo_dict[component]
        n_plus_1.append((component, n_plus_1_weight))
    
    # Check for P=O impurity using the existing method
    if check_po_impurity(components):
        po_impurity_weight = average_mw - 15.977156
        n_plus_1.append(("P=O", po_impurity_weight))
    
    # Check for MOE components that indicate 2-O-CH3 impurity
    moe_components = ["MOE G", "MOE A", "MOE MeU", "MOE MeC", 
                      "MOE G P=O", "MOE A P=O", "MOE MeU P=O", "MOE MeC P=O",
                      "MOE G nucleoside","MOE A nucleoside","MOE MeU nucleoside","MOE MeC nucleoside"]
    
    has_moe_components = any(component in moe_components for component in components)
    
    if has_moe_components:
        # Calculate 2-O-CH3 impurity: x - 44.026215
        moe_impurity_weight = average_mw - 44.026215
        n_plus_1.append(("2-O-CH3", moe_impurity_weight))
    
    # Define lists of nucleotide types to check for debase impurities
    a_nucleotides = ["MOE A", "dA", "NMA A", "mA", "cEt A", "fA", 
                     "MOE A P=O", "dA P=O", "NMA A P=O", "mA P=O", "cEt A P=O", "fA P=O",
                     "MOE A nucleoside", "dA nucleoside", "NMA A nucleoside", "mA nucleoside", 
                     "cEt A nucleoside", "fA nucleoside"]
    
    c_nucleotides = ["MOE MeC", "dMeC", "NMA MeC", "cEt MeC", 
                     "MOE MeC P=O", "dMeC P=O", "NMA MeC P=O", "cEt MeC P=O",
                     "MOE MeC nucleoside", "dMeC nucleoside", "NMA MeC nucleoside", 
                     "cEt MeC nucleoside"]
    
    g_nucleotides = ["MOE G", "dG", "NMA G", "mG", "cEt G", "fG", 
                     "MOE G P=O", "dG P=O", "NMA G P=O", "mG P=O", "cEt G P=O", "fG P=O",
                     "MOE G nucleoside", "dG nucleoside", "NMA G nucleoside", "mG nucleoside", 
                     "cEt G nucleoside", "fG nucleoside"]
    
    t_nucleotides = ["dT", "dT P=O", "dT nucleoside"]
    
    u_nucleotides = ["MOE MeU", "NMA MeU", "cEt MeU", 
                     "MOE MeU P=O", "NMA MeU P=O", "cEt MeU P=O",
                     "MOE MeU nucleoside", "NMA MeU nucleoside", "cEt MeU nucleoside"]
    
    # Check for GalNAc component
    has_galnac = "GalNAc" in components
    
    # Check for debase impurities
    has_a_nucleotides = any(component in a_nucleotides for component in components)
    has_c_nucleotides = any(component in c_nucleotides for component in components)
    has_g_nucleotides = any(component in g_nucleotides for component in components)
    has_t_nucleotides = any(component in t_nucleotides for component in components)
    has_u_nucleotides = any(component in u_nucleotides for component in components)
    
    # Add debase impurities if corresponding nucleotides are found
    if has_a_nucleotides:
        debase_a_weight = average_mw - 117.051755
        n_plus_1.append(("Debase-A", debase_a_weight))
        
        # Add AMPA impurity for A nucleotides
        ampa_weight = average_mw + 98.0368
        n_plus_1.append(("AMPA", ampa_weight))
    
    if has_c_nucleotides:
        debase_c_weight = average_mw - 107.056172
        n_plus_1.append(("Debase-C", debase_c_weight))
        
        # Add OPC impurity for MeC nucleotides
        opc_weight = average_mw + 80.0262
        n_plus_1.append(("OPC", opc_weight))
        
        # Add Deamination impurity for MeC nucleotides
        deamination_weight = average_mw + 0.984
        n_plus_1.append(("Deamination", deamination_weight))
        
        # Add DMT-C-phosphonates impurities for MeC nucleotides
        dmt_thiolation_weight = average_mw + 270.1586
        dmt_oxygen_weight = average_mw + 286.1358
        n_plus_1.append(("DMT-C-phosphonates (Thiolation incomplete)", dmt_thiolation_weight))
        n_plus_1.append(("DMT-C-phosphonates (Oxygen generation incomplete)", dmt_oxygen_weight))
    
    if has_g_nucleotides:
        debase_g_weight = average_mw - 133.04667
        n_plus_1.append(("Debase-G", debase_g_weight))
        
        # Add ADP impurity for G nucleotides
        adp_weight = average_mw + 41.026549
        n_plus_1.append(("ADP", adp_weight))
        
        # Add n+iBu impurity for G nucleotides
        n_ibu_weight = average_mw + 70.0419
        n_plus_1.append(("n+iBu", n_ibu_weight))
        
        # Add IDP impurity for G nucleotides
        idp_weight = average_mw + 69.0579
        n_plus_1.append(("IDP", idp_weight))
    
    if has_t_nucleotides:
        debase_t_weight = average_mw - 108.040188
        n_plus_1.append(("Debase-T", debase_t_weight))
        
        # Add CNET impurity for T nucleotides
        cnet_weight = average_mw + 53.034374
        n_plus_1.append(("CNET", cnet_weight))
    
    if has_u_nucleotides:
        debase_u_weight = average_mw - 108.040188
        n_plus_1.append(("Debase-U", debase_u_weight))
    
    # Add GalNAc-specific impurities if GalNAc is present
    if has_galnac:
        # GalNAc-related impurities
        penta_galnac_weight = average_mw + 923.4951
        n_plus_1.append(("Penta GalNAc", penta_galnac_weight))
        
        tetra_galnac_weight = average_mw + 549.2898
        n_plus_1.append(("Tetra GalNAc", tetra_galnac_weight))
        
        tha_branch_lost_weight = average_mw + 347.2053
        n_plus_1.append(("THA Branch Lost", tha_branch_lost_weight))
        
        hexgal_weight = average_mw - 303.1682
        n_plus_1.append(("HexGal", hexgal_weight))
        
        galnac_h2o_weight = average_mw - 203.0794
        n_plus_1.append(("GalNAc+H2O", galnac_h2o_weight))
        
        n_trisgalnac_weight = average_mw - 1226.6632
        n_plus_1.append(("n-TrisGalNAc", n_trisgalnac_weight))
        
        ahi_ac_weight = average_mw - 1297.7003
        n_plus_1.append(("AHI+Ac", ahi_ac_weight))
        
        n_p_ah_weight = average_mw - 179.0712
        n_plus_1.append(("n+p(AH)", n_p_ah_weight))
        
        n_ade_tea_weight = average_mw - 33.9341
        n_plus_1.append(("n-Ade+TEA", n_ade_tea_weight))
        
        galnac_ac_weight = average_mw - 42.0105
        n_plus_1.append(("GalNAc-Ac", galnac_ac_weight))
        
        aminohexyl_phosphate_weight = average_mw - 1339.7109
        n_plus_1.append(("Aminohexyl Phosphate", aminohexyl_phosphate_weight))
    
    # Add universal impurities (always present)
    # Disulfate/Sulfate impurity for all sequences
    disulfate_weight = average_mw + 15.977156
    n_plus_1.append(("Disulfate/Sulfate", disulfate_weight))
    
    # Uny-TP impurity for all sequences
    uny_tp_weight = average_mw + 276.9811
    n_plus_1.append(("Uny-TP", uny_tp_weight))
    
    # Mobile phase-related impurities (calculated only once, only standard charge state)
    mobile_phase_impurities = []
    
    # Mobile phase pseudo-peaks (always present)
    mobile_phase_1_weight = average_mw + 219.2
    mobile_phase_impurities.append(("Mobile phase false peak 1", mobile_phase_1_weight))
    
    mobile_phase_2_weight = average_mw + 234.8
    mobile_phase_impurities.append(("Mobile phase false peak 2", mobile_phase_2_weight))
    
    mobile_phase_3_weight = average_mw + 255.2
    mobile_phase_impurities.append(("Mobile phase false peak 3", mobile_phase_3_weight))
    
    mobile_phase_4_weight = average_mw + 467.2
    mobile_phase_impurities.append(("Mobile phase false peak 4", mobile_phase_4_weight))
    
    # Tributylamine addition peaks (always present) - MODIFIED: removed "N+" prefix
    n_tbua_weight = average_mw + 185.2143
    mobile_phase_impurities.append(("TBuA", n_tbua_weight))
    
    n_tbua_o_weight = average_mw + 201.2092
    mobile_phase_impurities.append(("TBuA+O", n_tbua_o_weight))
    
    return n_minus_1, n_plus_1, mobile_phase_impurities

def calculate_charge_states(impurities, mobile_phase_impurities, charge, average_mw):
    """
    Calculate molecular weights under different charge states.
    
    Parameters:
    impurities (list): List of (impurity_name, molecular_weight) tuples
    mobile_phase_impurities (list): List of mobile phase related impurities
    charge (int): Charge state
    average_mw (float): Average molecular weight
    
    Returns:
    dict: Dictionary with standard, Na-addition, K-addition, Fe-addition peaks, and mobile phase peaks
    """
    results = {
        'standard': [],
        'na_addition': [],
        'k_addition': [],
        'fe_addition': [],
        'mobile_phase': []
    }
    
    proton_mass = 1.007825
    na_mass = 22.989769
    k_mass = 38.963707
    fe_mass = 55.9349
    
    # For average molecular weight itself
    standard_avg = (average_mw - charge * proton_mass) / charge
    na_avg = (average_mw + na_mass - charge * proton_mass) / charge
    k_avg = (average_mw + k_mass - charge * proton_mass) / charge
    fe_avg = (average_mw + fe_mass - charge * proton_mass) / charge
    
    results['standard'].append(("Average MW", standard_avg))
    results['na_addition'].append(("Average MW", na_avg))
    results['k_addition'].append(("Average MW", k_avg))
    results['fe_addition'].append(("Average MW", fe_avg))
    
    # For each regular impurity (with all charge state variants)
    for impurity_name, mol_weight in impurities:
        # Standard charge state calculation: (x - n*1.007825)/n
        standard_mw = (mol_weight - charge * proton_mass) / charge
        results['standard'].append((impurity_name, standard_mw))
        
        # Na addition peak: (x + 22.989769 - n*1.007825)/n
        na_mw = (mol_weight + na_mass - charge * proton_mass) / charge
        results['na_addition'].append((impurity_name, na_mw))
        
        # K addition peak: (x + 38.963707 - n*1.007825)/n
        k_mw = (mol_weight + k_mass - charge * proton_mass) / charge
        results['k_addition'].append((impurity_name, k_mw))
        
        # Fe addition peak: (x + 55.9349 - n*1.007825)/n
        fe_mw = (mol_weight + fe_mass - charge * proton_mass) / charge
        results['fe_addition'].append((impurity_name, fe_mw))
    
    # For mobile phase impurities (only standard charge state)
    for impurity_name, mol_weight in mobile_phase_impurities:
        standard_mw = (mol_weight - charge * proton_mass) / charge
        results['mobile_phase'].append((impurity_name, standard_mw))
    
    return results

def get_all_results(sequence, average_mw, charge):
    """
    Calculate and organize all results for a given sequence, average MW, and charge.
    
    Parameters:
    sequence (str): Oligonucleotide sequence
    average_mw (float): Average molecular weight
    charge (int): Charge state
    
    Returns:
    dict: Dictionary with all results or None if error occurs
    """
    # Calculate the molecular weight of the sequence (now using simple addition)
    total_weight, component_weights, found_all, missing_components = calculate_molecular_weight(sequence)
    
    if not found_all:
        return {
            'success': False,
            'missing_components': missing_components,
            'oligo_dict': create_oligonucleotide_dictionary()
        }
    
    # Calculate N-1 and N+1 impurities using the modified method
    n_minus_1, n_plus_1, mobile_phase_impurities = calculate_impurities(sequence, average_mw)
    
    # Combine all regular impurities (exclude mobile phase impurities)
    all_impurities = []
    
    # Define the specific order for displaying impurities
    special_impurity_order = [
        "P=O", "2-O-CH3", "CNET", "Debase-A", "Debase-C", "Debase-G",
        "Debase-T", "Debase-U", "Disulfate/Sulfate", "ADP", "AMPA", "OPC", 
        "dGalNAc", "Deamination", "DMT-C-phosphonates (Thiolation incomplete)",
        "DMT-C-phosphonates (Oxygen generation incomplete)", "n+iBu", "IDP", "Uny-TP",
        "Penta GalNAc", "Tetra GalNAc", "THA Branch Lost", "HexGal", "GalNAc+H2O",
        "n-TrisGalNAc", "AHI+Ac", "n+p(AH)", "n-Ade+TEA", "GalNAc-Ac", 
        "Aminohexyl Phosphate"
    ]
    
    # Add special impurities first (in order)
    special_added = set()
    for impurity_type in special_impurity_order:
        for component, weight in n_plus_1:
            if component == impurity_type and component not in special_added:
                all_impurities.append((component, weight))
                special_added.add(component)
    
    # Add N+1 impurities (sequence components that are not special impurities)
    sequence_components_added = set()
    for component, weight in n_plus_1:
        if component not in special_impurity_order and component not in sequence_components_added:
            all_impurities.append((f"N+{component}", weight))
            sequence_components_added.add(component)
    
    # Add N-1 impurities (sequence components)
    for component, weight in n_minus_1:
        if component not in sequence_components_added:
            all_impurities.append((f"N-{component}", weight))
            sequence_components_added.add(component)
        else:
            # Even if the component was already added as N+, we still want to add N-
            all_impurities.append((f"N-{component}", weight))
    
    # Calculate charge state molecular weights
    charge_results = calculate_charge_states(all_impurities, mobile_phase_impurities, charge, average_mw)
    
    return {
        'success': True,
        'sequence': sequence,
        'calculated_mw': total_weight,
        'average_mw': average_mw,
        'charge': charge,
        'component_weights': component_weights,
        'special_impurity_order': special_impurity_order,
        'all_impurities': all_impurities,
        'mobile_phase_impurities': mobile_phase_impurities,
        'n_minus_1': n_minus_1,
        'n_plus_1': n_plus_1,
        'charge_results': charge_results
    }

def format_results(results):
    """
    Format the results for display.
    
    Parameters:
    results (dict): Dictionary with calculation results
    
    Returns:
    str: Formatted results as a string
    """
    if not results['success']:
        output = "\nWarning: The following components were not found in the dictionary:\n"
        for component in results['missing_components']:
            output += f"  - '{component}'\n"
        
        output += "\nPlease check the spelling and try again. Suggested components that might match:\n"
        
        oligo_dict = results['oligo_dict']
        for missing in results['missing_components']:
            found_suggestions = False
            for key in oligo_dict.keys():
                # Look for partial matches
                if missing.lower() in key.lower():
                    if not found_suggestions:
                        output += f"For '{missing}':\n"
                        found_suggestions = True
                    output += f"  - {key}\n"
            if not found_suggestions:
                output += f"No suggestions found for '{missing}'. Type 'list' to see all available terms.\n"
        
        return output
    
    output = f"\nSequence: {results['sequence']}\n"
    output += f"Calculated molecular weight (simple addition): {results['calculated_mw']:.5f}\n"
    output += f"Average molecular weight: {results['average_mw']:.5f}\n"
    output += f"Charge state: {results['charge']}\n"
    
    # Display impurities (molecular weights)
    output += "\nImpurities (molecular weights):\n"
    
    # First display special impurities in their specific order
    for impurity_type in results['special_impurity_order']:
        for name, weight in results['all_impurities']:
            if name == impurity_type:
                output += f"  {name}: {weight:.5f}\n"
    
    # Then display N+/N- impurities (from sequence components)
    for name, weight in results['all_impurities']:
        if name.startswith("N+") or name.startswith("N-"):
            output += f"  {name}: {weight:.5f}\n"
    
    # Display mobile phase impurities separately
    output += "\nMobile phase-related impurities (molecular weights):\n"
    for name, weight in results['mobile_phase_impurities']:
        output += f"  {name}: {weight:.5f}\n"
    
    # Display results for charge states
    output += f"\nMolecular weights under charge state {results['charge']}:\n"
    
    # Standard (deprotonated)
    output += "\nStandard (deprotonated):\n"
    charge = results['charge']
    
    # First display Average MW
    for name, mw in results['charge_results']['standard']:
        if name == "Average MW":
            output += f"  {charge}-{name}: {mw:.5f}\n"
    
    # Then display special impurities in their specific order
    for impurity_type in results['special_impurity_order']:
        for name, mw in results['charge_results']['standard']:
            if name == impurity_type:
                output += f"  {charge}-{name}: {mw:.5f}\n"
    
    # Then display N+/N- impurities
    for name, mw in results['charge_results']['standard']:
        if name.startswith("N+") or name.startswith("N-"):
            output += f"  {charge}-{name}: {mw:.5f}\n"
    
    # Add mobile phase impurities to Standard (deprotonated) section
    output += "\n  Mobile phase-related impurities:\n"
    for name, mw in results['charge_results']['mobile_phase']:
        output += f"    {charge}-{name}: {mw:.5f}\n"
    
    # Na-addition peaks (only Average MW and P=O)
    output += "\nNa-addition peaks:\n"
    
    # Display Average MW
    for name, mw in results['charge_results']['na_addition']:
        if name == "Average MW":
            output += f"  Na-{name}: {mw:.5f}\n"
    
    # Display P=O if present
    for name, mw in results['charge_results']['na_addition']:
        if name == "P=O":
            output += f"  Na-{name}: {mw:.5f}\n"
    
    # K-addition peaks (only Average MW and P=O)
    output += "\nK-addition peaks:\n"
    
    # Display Average MW
    for name, mw in results['charge_results']['k_addition']:
        if name == "Average MW":
            output += f"  K-{name}: {mw:.5f}\n"
    
    # Display P=O if present
    for name, mw in results['charge_results']['k_addition']:
        if name == "P=O":
            output += f"  K-{name}: {mw:.5f}\n"
    
    # Fe-addition peaks (only Average MW and P=O)
    output += "\nFe-addition peaks:\n"
    
    # Display Average MW
    for name, mw in results['charge_results']['fe_addition']:
        if name == "Average MW":
            output += f"  Fe-{name}: {mw:.5f}\n"
    
    # Display P=O if present
    for name, mw in results['charge_results']['fe_addition']:
        if name == "P=O":
            output += f"  Fe-{name}: {mw:.5f}\n"
    
    # Reference section for other Na, K, Fe impurities
    output += "\n" + "="*60 + "\n"
    output += "REFERENCE: Other Na/K/Fe addition impurities\n"
    output += "(For detailed analysis reference only)\n"
    output += "="*60 + "\n"
    
    # Na-addition reference impurities (excluding Average MW and P=O)
    output += "\nOther Na-addition impurities:\n"
    found_na_reference = False
    for name, mw in results['charge_results']['na_addition']:
        if name not in ["Average MW", "P=O"]:
            output += f"  Na-{name}: {mw:.5f}\n"
            found_na_reference = True
    if not found_na_reference:
        output += "  (No additional Na impurities)\n"
    
    # K-addition reference impurities (excluding Average MW and P=O)
    output += "\nOther K-addition impurities:\n"
    found_k_reference = False
    for name, mw in results['charge_results']['k_addition']:
        if name not in ["Average MW", "P=O"]:
            output += f"  K-{name}: {mw:.5f}\n"
            found_k_reference = True
    if not found_k_reference:
        output += "  (No additional K impurities)\n"
    
    # Fe-addition reference impurities (excluding Average MW and P=O)
    output += "\nOther Fe-addition impurities:\n"
    found_fe_reference = False
    for name, mw in results['charge_results']['fe_addition']:
        if name not in ["Average MW", "P=O"]:
            output += f"  Fe-{name}: {mw:.5f}\n"
            found_fe_reference = True
    if not found_fe_reference:
        output += "  (No additional Fe impurities)\n"
    
    return output

def save_results_to_file(filename, results):
    """
    Save results to a file.
    
    Parameters:
    filename (str): Name of the file to save results to
    results (str): Formatted results string
    
    Returns:
    bool: True if successful, False otherwise
    """
    try:
        with open(filename, 'w') as file:
            file.write(results)
        return True
    except Exception as e:
        print(f"Error saving to file: {e}")
        return False

def explain_impurities():
    """Display explanations for different impurity types."""
    print("\nImpurity Type Explanations:")
    print("--------------------------")
    print("N+ Impurities: Represent addition of one component to the sequence")
    print("  Calculated as: Average MW + component molecular weight")
    print("N- Impurities: Represent deletion of one component from the sequence")
    print("  Calculated as: Average MW - component molecular weight")
    print("P=O: Phosphodiester linkage instead of phosphorothioate (calculated as average MW - 15.977156)")
    print("2-O-CH3: Loss of a methoxy group in MOE nucleotides")
    print("Debase-X: Loss of nucleobase X (A, C, G, T, or U)")
    print("AMPA: Addition of aminopropyl group to adenine nucleotides")
    print("OPC: Addition of oxidation product to cytosine nucleotides")
    print("ADP: Addition of deamination product to guanine nucleotides")
    print("CNET: Addition of cyanoethyl group to thymine nucleotides")
    print("Disulfate/Sulfate: Addition of sulfate group")
    print("Deamination: Deamination product for methylcytosine nucleotides")
    print("DMT-C-phosphonates (Thiolation incomplete): Incomplete thiolation during synthesis")
    print("DMT-C-phosphonates (Oxygen generation incomplete): Incomplete oxygen generation during synthesis")
    print("n+iBu: Addition of isobutyl group to guanine nucleotides")
    print("IDP: Isobutyl deamination product for guanine nucleotides")
    print("Uny-TP: Universal impurity present in all sequences")
    print("\nGalNAc-specific impurities (when GalNAc is present):")
    print("Penta GalNAc: Addition of five GalNAc units")
    print("Tetra GalNAc: Addition of four GalNAc units")
    print("THA Branch Lost: Loss of triantennary branching")
    print("HexGal: Loss of hexose galactose")
    print("GalNAc+H2O: GalNAc with water addition")
    print("n-TrisGalNAc: Loss of three GalNAc units")
    print("AHI+Ac: Aminohexyl with acetyl group")
    print("n+p(AH): Addition of phosphorylated aminohexyl")
    print("n-Ade+TEA: Loss of adenine with triethylamine addition")
    print("GalNAc-Ac: GalNAc with loss of acetyl group")
    print("Aminohexyl Phosphate: Aminohexyl phosphate impurity")
    print("\nMobile phase-related impurities (standard charge state only):")
    print("Mobile phase false peaks 1-4: Pseudo-peaks from mobile phase components")
    print("TBuA: Addition of tributylamine (calculated only once)")
    print("TBuA+O: Addition of oxidized tributylamine (calculated only once)")
    

def batch_process(filename):
    """
    Process multiple sequences from a file.
    
    Parameters:
    filename (str): Name of the file containing sequences to process
    """
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        output_filename = filename.split('.')[0] + "_results.txt"
        with open(output_filename, 'w') as outfile:
            for i, line in enumerate(lines):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = [part.strip() for part in line.split(',')]
                if len(parts) != 3:
                    outfile.write(f"Error in line {i+1}: Incorrect format\n")
                    continue
                
                sequence = parts[0]
                try:
                    average_mw = float(parts[1])
                    charge = int(parts[2])
                except ValueError:
                    outfile.write(f"Error in line {i+1}: Invalid numerical values\n")
                    continue
                
                if charge <= 0:
                    outfile.write(f"Error in line {i+1}: Charge must be positive\n")
                    continue
                
                results = get_all_results(sequence, average_mw, charge)
                formatted_results = format_results(results)
                
                outfile.write(f"Results for sequence {i+1}: {sequence}\n")
                outfile.write(formatted_results)
                outfile.write("\n" + "="*80 + "\n\n")
        
        print(f"Batch processing complete. Results saved to {output_filename}")
    except Exception as e:
        print(f"Error during batch processing: {e}")

def display_help():
    """Display help information about the program."""
    print("\nOligonucleotide Molecular Weight Calculator - Help")
    print("================================================")
    print("\nThis program calculates molecular weights and impurities for oligonucleotide sequences.")
    print("\nMain Features:")
    print("1. Calculate molecular weight of a sequence (simple addition of components)")
    print("2. Identify N+1 and N-1 impurities")
    print("3. Calculate molecular weights under different charge states")
    print("4. Mobile phase-related impurities (standard charge state only)")
    print("5. Save results to a file")
    print("6. Batch process multiple sequences")
    print("\nAvailable Commands:")
    print("- 'sequence, average_mw, charge': Calculate results for a sequence")
    print("- 'list': Display all available oligonucleotide components")
    print("- 'save': Save the last calculated results to a file")
    print("- 'help': Display this help information")
    print("- 'explain': Show explanations for different impurity types")
    print("- 'batch filename.txt': Process multiple sequences from a file")
    print("- 'quit': Exit the program")
    
    print("\nInput Format:")
    print("- Sequence must use terms exactly as they appear in the dictionary")
    print("- Components in the sequence must be separated by hyphens")
    print("- All three inputs (sequence, average MW, charge) must be provided, separated by commas")
    
    print("\nExample:")
    print("MOE G-MOE A-dT, 1250.45, 3")
    
    print("\nMODIFIED: Molecular Weight Calculation:")
    print("The calculated molecular weight is now the simple sum of all components")
    print("(No water loss correction is applied)")
    print("For sequence 'MOE G-MOE A', calculated MW = MOE G weight + MOE A weight")
    
    print("\nN+1 and N-1 Impurity Calculations:")
    print("N+1 impurities: Average MW + component molecular weight")
    print("N-1 impurities: Average MW - component molecular weight")
    print("For sequence 'MOE G-MOE A', impurities will be:")
    print("  N+MOE G, N+MOE A, N-MOE G, N-MOE A")
    
    print("\nP=O Calculation:")
    print("P=O impurity is calculated as: Average MW - 15.977156")
    print("P=O impurity is detected when sequence contains any eligible nucleotides")
    
    print("\nNote: Mobile phase-related impurities are calculated only once")
    print("Na/K/Fe addition peaks show only Average MW and P=O in main sections")

def main():
    """
    Main function to interact with the user.
    """
    print("\n=============================================================")
    print("    Oligonucleotide Molecular Weight Calculator v2.8")
    print("          (Modified: Simple Addition Mode)")
    print("=============================================================")
    print("\nIMPORTANT CHANGE: Molecular weight calculation now uses simple addition")
    print("(Water loss correction has been removed)")
    print("\nINPUT INSTRUCTIONS:")
    print("1. Your sequence must use terms EXACTLY as they appear in our dictionary")
    print("   (Type 'list' to see all available oligonucleotides)")
    print("\n2. Create your sequence by combining dictionary terms with hyphens")
    print("   Example: 'MOE G-MOE A-MOE MeC'")
    print("\n3. Please provide ALL THREE of the following inputs, separated by commas:")
    print("   a) Sequence (combination of terms from the dictionary)")
    print("   b) Average molecular weight (a number)")
    print("   c) Charge state (a positive integer)")
    print("\nComplete example: 'MOE G-MOE A-NMA MeC, 1250.45, 3'")
    print("\nAdditional commands:")
    print("- Type 'list' to view the complete dictionary of available terms")
    print("- Type 'save' to save the last results to a file")
    print("- Type 'help' to view help information")
    print("- Type 'explain' to see explanations for impurity types")
    print("- Type 'batch filename.txt' to process multiple sequences from a file")
    print("- Type 'quit' to exit the program")
  
    
    oligo_dict = create_oligonucleotide_dictionary()
    last_results = None
    
    while True:
        print("\n" + "-" * 65)
        user_input = input("Enter: sequence, average MW, charge (or list/save/quit): ")
        
        if user_input.lower() == 'quit':
            print("\nExiting program. Goodbye!")
            break
            
        elif user_input.lower() == 'list':
            print("\nAvailable oligonucleotides in dictionary:")
            
            # Group nucleotides by type for better readability
            types = {
                "MOE": [],
                "DNA (d)": [],
                "NMA": [],
                "m-nucleotides (m)": [],
                "cEt": [],
                "f-nucleotides (f)": [],
                "Other": []
            }
            
            for key, value in oligo_dict.items():
                if key.startswith("MOE"):
                    types["MOE"].append((key, value))
                elif key.startswith("d"):
                    types["DNA (d)"].append((key, value))
                elif key.startswith("NMA"):
                    types["NMA"].append((key, value))
                elif key.startswith("m"):
                    types["m-nucleotides (m)"].append((key, value))
                elif key.startswith("cEt"):
                    types["cEt"].append((key, value))
                elif key.startswith("f"):
                    types["f-nucleotides (f)"].append((key, value))
                else:
                    types["Other"].append((key, value))
            
            # Display grouped nucleotides
            for category, items in types.items():
                if items:
                    print(f"\n{category} nucleotides:")
                    # Sort items within each category
                    for key, value in sorted(items):
                        print(f"  {key}: {value}")
            
            print("\nRemember: Your sequence must be a combination of these exact terms separated by hyphens.")
            
        elif user_input.lower() == 'help':
            display_help()
            
        elif user_input.lower() == 'explain':
            explain_impurities()
            
        elif user_input.lower().startswith('batch '):
            filename = user_input[6:].strip()
            if filename:
                batch_process(filename)
            else:
                print("Error: Please provide a filename for batch processing.")
                print("Example: batch sequences.txt")
                
        elif user_input.lower() == 'save':
            if last_results is None:
                print("No results to save. Please calculate results first.")
            else:
                filename = input("Enter filename to save results (default: oligo_results.txt): ")
                if not filename:
                    filename = "oligo_results.txt"
                
                if save_results_to_file(filename, last_results):
                    print(f"Results successfully saved to {filename}")
                else:
                    print(f"Failed to save results to {filename}")
                    
        else:
            # Parse the input
            parts = [part.strip() for part in user_input.split(',')]
            
            if len(parts) != 3:
                print("\nError: You must provide all three required inputs separated by commas:")
                print("1. Sequence (from dictionary terms)")
                print("2. Average molecular weight")
                print("3. Charge state")
                print("Example: 'MOE G-MOE A, 500, 3'")
                continue
            
            sequence = parts[0]
            
            try:
                average_mw = float(parts[1])
                charge = int(parts[2])
            except ValueError:
                print(f"\nError: Invalid molecular weight or charge. Please enter valid numbers.")
                continue
            
            if charge <= 0:
                print("\nError: Charge must be a positive integer.")
                continue
            
            # Calculate all results
            results = get_all_results(sequence, average_mw, charge)
            
            # Format and display results
            formatted_results = format_results(results)
            print(formatted_results)
            
            # Save the formatted results for potential saving to file later
            last_results = formatted_results

if __name__ == "__main__":
    try:
        # Check for command-line arguments
        import sys
        if len(sys.argv) > 1:
            if sys.argv[1].lower() == 'batch' and len(sys.argv) > 2:
                batch_process(sys.argv[2])
                sys.exit(0)
            elif sys.argv[1].lower() == 'help':
                display_help()
                sys.exit(0)
        
        # No command-line args, run interactive mode
        main()
    except KeyboardInterrupt:
        print("\nProgram interrupted. Exiting.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        print("Please report this error to the developer.")

# Example usage:
if __name__ == "__main__":
    # Test the modified function
    print("Testing modified calculate_molecular_weight function:")
    print("="*50)
    
    # Test case 1
    sequence1 = "MOE G-MOE A-dT"
    result1 = calculate_molecular_weight(sequence1)
    print(f"Sequence: {sequence1}")
    print(f"Total weight (simple addition): {result1[0]:.5f}")
    print(f"Components: {result1[1]}")
    print()
    
    # Test case 2
    sequence2 = "GalNAc-MOE G"
    result2 = calculate_molecular_weight(sequence2)
    print(f"Sequence: {sequence2}")
    print(f"Total weight (simple addition): {result2[0]:.5f}")
    print(f"Components: {result2[1]}")
    print()
    
    print("Note: No water loss correction is applied - weights are simply summed.")
