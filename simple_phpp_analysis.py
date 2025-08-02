#!/usr/bin/env python3
"""
Simple Phenotype Phase Plane Analysis

A basic implementation to get phenotype phase plane analysis working.
"""

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def test_model_growth():
    """Test if the model can grow with basic constraints."""
    print("Testing model growth capabilities...")
    
    # Load model
    model = cobra.io.read_sbml_model("bigg_models/iJO1366.xml")
    print(f"Model loaded: {model.name}")
    
    # Test basic growth
    solution = model.optimize()
    print(f"Default growth rate: {solution.objective_value:.6f}")
    
    # Test with glucose
    glucose_rxn = model.reactions.get_by_id("EX_glc__D_e")
    glucose_rxn.lower_bound = -10.0
    
    solution = model.optimize()
    print(f"Growth with glucose: {solution.objective_value:.6f}")
    
    return model

def simple_phenotype_phase_plane(model, substrate_rxn, substrate_name):
    """Simple phenotype phase plane analysis."""
    print(f"\nRunning simple phenotype phase plane for {substrate_name}...")
    
    # Define ranges
    substrate_fluxes = np.linspace(0, -10, 11)  # 0 to -10 mmol/gDW/h
    o2_fluxes = np.linspace(0, -15, 16)         # 0 to -15 mmol/gDW/h
    
    results = []
    
    for sub_flux in substrate_fluxes:
        for o2_flux in o2_fluxes:
            # Set bounds
            substrate_rxn.lower_bound = sub_flux
            o2_rxn = model.reactions.get_by_id("EX_o2_e")
            o2_rxn.lower_bound = o2_flux
            
            # Optimize
            solution = model.optimize()
            
            growth_rate = solution.objective_value if solution.status == 'optimal' else 0.0
            
            results.append({
                'substrate_flux': sub_flux,
                'o2_flux': o2_flux,
                'growth_rate': growth_rate
            })
    
    return pd.DataFrame(results)

def plot_results(df, substrate_name):
    """Plot the phenotype phase plane."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Filter for growth
    growth_data = df[df.growth_rate > 1e-6]
    
    if len(growth_data) > 0:
        sc = ax.scatter(growth_data.substrate_flux, growth_data.o2_flux, 
                       c=growth_data.growth_rate, cmap='viridis', s=50)
        plt.colorbar(sc, label='Growth rate (1/h)')
    else:
        ax.scatter(df.substrate_flux, df.o2_flux, c='gray', alpha=0.3, s=10)
    
    ax.set_xlabel(f'{substrate_name.capitalize()} uptake (mmol/gDW/h)')
    ax.set_ylabel('O₂ uptake (mmol/gDW/h)')
    ax.set_title(f'Phenotype Phase Plane: {substrate_name.capitalize()} vs O₂')
    
    plt.tight_layout()
    plt.savefig(f'{substrate_name}_simple_phpp.png', dpi=150, bbox_inches='tight')
    plt.show()

def main():
    """Main function."""
    print("=== Simple Phenotype Phase Plane Analysis ===\n")
    
    # Test model
    model = test_model_growth()
    
    # Test with glucose first (should work)
    print("\n" + "="*50)
    print("GLUCOSE ANALYSIS")
    print("="*50)
    
    glucose_rxn = model.reactions.get_by_id("EX_glc__D_e")
    glucose_results = simple_phenotype_phase_plane(model, glucose_rxn, "glucose")
    
    print(f"Glucose results shape: {glucose_results.shape}")
    print(f"Points with growth: {len(glucose_results[glucose_results.growth_rate > 1e-6])}")
    if len(glucose_results[glucose_results.growth_rate > 1e-6]) > 0:
        print(f"Max growth rate: {glucose_results.growth_rate.max():.6f}")
    
    plot_results(glucose_results, "glucose")
    
    # Test with acetate
    print("\n" + "="*50)
    print("ACETATE ANALYSIS")
    print("="*50)
    
    acetate_rxn = model.reactions.get_by_id("EX_ac_e")
    acetate_results = simple_phenotype_phase_plane(model, acetate_rxn, "acetate")
    
    print(f"Acetate results shape: {acetate_results.shape}")
    print(f"Points with growth: {len(acetate_results[acetate_results.growth_rate > 1e-6])}")
    if len(acetate_results[acetate_results.growth_rate > 1e-6]) > 0:
        print(f"Max growth rate: {acetate_results.growth_rate.max():.6f}")
    
    plot_results(acetate_results, "acetate")
    
    # Save results
    glucose_results.to_csv("glucose_simple_phpp.csv", index=False)
    acetate_results.to_csv("acetate_simple_phpp.csv", index=False)
    
    print("\nAnalysis complete!")
    print("Files saved: glucose_simple_phpp.csv, acetate_simple_phpp.csv")
    print("Plots saved: glucose_simple_phpp.png, acetate_simple_phpp.png")

if __name__ == "__main__":
    main() 