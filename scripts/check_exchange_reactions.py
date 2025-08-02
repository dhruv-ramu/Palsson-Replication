#!/usr/bin/env python3
"""
Check Exchange Reactions in iJO1366 Model

This script lists all available exchange reactions in the iJO1366 model
to help identify the correct reaction IDs for glucose minimal medium setup.
"""

import cobra

def check_exchange_reactions():
    """Check and list exchange reactions in the iJO1366 model."""
    print("Loading iJO1366 model...")
    model = cobra.io.read_sbml_model("bigg_models/iJO1366.xml")
    
    print(f"\nModel: {model.name}")
    print(f"Total reactions: {len(model.reactions)}")
    print(f"Exchange reactions: {len(model.exchanges)}")
    
    print("\nExchange reactions starting with 'EX_':")
    ex_reactions = [rxn for rxn in model.exchanges if rxn.id.startswith('EX_')]
    
    # Group by common patterns
    glucose_patterns = ['glc', 'glucose']
    oxygen_patterns = ['o2', 'oxygen']
    co2_patterns = ['co2', 'carbon']
    water_patterns = ['h2o', 'water']
    proton_patterns = ['h_e', 'proton']
    nitrogen_patterns = ['nh4', 'ammonium', 'nitrogen']
    phosphate_patterns = ['pi', 'phosphate']
    sulfate_patterns = ['so4', 'sulfate']
    potassium_patterns = ['k', 'potassium']
    sodium_patterns = ['na', 'sodium']
    magnesium_patterns = ['mg', 'magnesium']
    calcium_patterns = ['ca', 'calcium']
    iron_patterns = ['fe', 'iron']
    chloride_patterns = ['cl', 'chloride']
    
    patterns = {
        'Glucose': glucose_patterns,
        'Oxygen': oxygen_patterns,
        'CO2': co2_patterns,
        'Water': water_patterns,
        'Proton': proton_patterns,
        'Nitrogen': nitrogen_patterns,
        'Phosphate': phosphate_patterns,
        'Sulfate': sulfate_patterns,
        'Potassium': potassium_patterns,
        'Sodium': sodium_patterns,
        'Magnesium': magnesium_patterns,
        'Calcium': calcium_patterns,
        'Iron': iron_patterns,
        'Chloride': chloride_patterns
    }
    
    for category, pattern_list in patterns.items():
        matches = []
        for rxn in ex_reactions:
            if any(pattern in rxn.id.lower() for pattern in pattern_list):
                matches.append(rxn.id)
        if matches:
            print(f"\n{category} reactions:")
            for match in sorted(matches):
                print(f"  {match}")
    
    print(f"\nAll exchange reactions ({len(ex_reactions)}):")
    for i, rxn in enumerate(sorted(ex_reactions, key=lambda x: x.id)):
        if i < 50:  # Show first 50
            print(f"  {rxn.id}")
        elif i == 50:
            print(f"  ... and {len(ex_reactions) - 50} more")
            break

if __name__ == "__main__":
    check_exchange_reactions() 