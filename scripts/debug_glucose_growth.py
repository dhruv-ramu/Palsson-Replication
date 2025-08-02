#!/usr/bin/env python3
"""
Debug Glucose Growth in iJO1366 Model

This script investigates why the iJO1366 model isn't growing on glucose minimal medium.
"""

import cobra

def debug_glucose_growth():
    """Debug glucose growth issues in the iJO1366 model."""
    print("Loading iJO1366 model...")
    model = cobra.io.read_sbml_model("bigg_models/iJO1366.xml")
    
    print(f"\nModel: {model.name}")
    print(f"Reactions: {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    print(f"Genes: {len(model.genes)}")
    
    # Test 1: Default model growth
    print("\n=== Test 1: Default Model Growth ===")
    solution = model.optimize()
    print(f"Default growth rate: {solution.objective_value:.6f} 1/h")
    print(f"Status: {solution.status}")
    
    # Test 2: Check biomass reaction
    print("\n=== Test 2: Biomass Reaction ===")
    biomass_rxn = model.objective
    print(f"Biomass reaction: {biomass_rxn}")
    # Get the actual biomass reaction
    biomass_reactions = [rxn for rxn in model.reactions if 'BIOMASS' in rxn.id]
    if biomass_reactions:
        print(f"Biomass reaction bounds: {biomass_reactions[0].bounds}")
    else:
        print("No biomass reaction found")
    
    # Test 3: Check glucose exchange reaction
    print("\n=== Test 3: Glucose Exchange Reaction ===")
    try:
        glucose_rxn = model.reactions.get_by_id('EX_glc__D_e')
        print(f"Glucose reaction: {glucose_rxn}")
        print(f"Glucose bounds: {glucose_rxn.bounds}")
        print(f"Glucose flux: {solution.fluxes.get('EX_glc__D_e', 'N/A')}")
    except:
        print("Glucose reaction not found")
    
    # Test 4: Check oxygen exchange reaction
    print("\n=== Test 4: Oxygen Exchange Reaction ===")
    try:
        oxygen_rxn = model.reactions.get_by_id('EX_o2_e')
        print(f"Oxygen reaction: {oxygen_rxn}")
        print(f"Oxygen bounds: {oxygen_rxn.bounds}")
        print(f"Oxygen flux: {solution.fluxes.get('EX_o2_e', 'N/A')}")
    except:
        print("Oxygen reaction not found")
    
    # Test 5: Check essential exchanges
    print("\n=== Test 5: Essential Exchange Reactions ===")
    essential_exchanges = ['EX_nh4_e', 'EX_pi_e', 'EX_so4_e', 'EX_k_e', 'EX_na1_e', 
                          'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e', 'EX_fe3_e', 'EX_cl_e']
    
    for ex_id in essential_exchanges:
        try:
            rxn = model.reactions.get_by_id(ex_id)
            flux = solution.fluxes.get(ex_id, 0)
            print(f"{ex_id}: bounds={rxn.bounds}, flux={flux:.6f}")
        except:
            print(f"{ex_id}: Not found")
    
    # Test 6: Try different glucose uptake rates
    print("\n=== Test 6: Different Glucose Uptake Rates ===")
    glucose_rates = [-1, -5, -10, -20]
    
    for rate in glucose_rates:
        model_copy = model.copy()
        try:
            glucose_rxn = model_copy.reactions.get_by_id('EX_glc__D_e')
            glucose_rxn.bounds = (rate, 0)
            solution = model_copy.optimize()
            print(f"Glucose uptake {rate}: growth = {solution.objective_value:.6f} 1/h")
        except:
            print(f"Glucose uptake {rate}: Error")
    
    # Test 7: Check if model can grow with all exchanges open
    print("\n=== Test 7: All Exchanges Open ===")
    model_copy = model.copy()
    for rxn in model_copy.exchanges:
        rxn.bounds = (-1000, 1000)
    
    solution = model_copy.optimize()
    print(f"All exchanges open: growth = {solution.objective_value:.6f} 1/h")
    
    # Test 8: Minimal medium with just glucose and oxygen
    print("\n=== Test 8: Minimal Medium (Glucose + Oxygen) ===")
    model_copy = model.copy()
    
    # Close all exchanges
    for rxn in model_copy.exchanges:
        rxn.bounds = (0, 0)
    
    # Open only glucose and oxygen
    try:
        glucose_rxn = model_copy.reactions.get_by_id('EX_glc__D_e')
        oxygen_rxn = model_copy.reactions.get_by_id('EX_o2_e')
        co2_rxn = model_copy.reactions.get_by_id('EX_co2_e')
        h2o_rxn = model_copy.reactions.get_by_id('EX_h2o_e')
        h_rxn = model_copy.reactions.get_by_id('EX_h_e')
        
        glucose_rxn.bounds = (-10, 0)
        oxygen_rxn.bounds = (-20, 0)
        co2_rxn.bounds = (0, 1000)
        h2o_rxn.bounds = (-1000, 1000)
        h_rxn.bounds = (-1000, 1000)
        
        solution = model_copy.optimize()
        print(f"Minimal medium: growth = {solution.objective_value:.6f} 1/h")
        
    except Exception as e:
        print(f"Minimal medium error: {e}")

if __name__ == "__main__":
    debug_glucose_growth() 