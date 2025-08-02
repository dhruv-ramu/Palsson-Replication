#!/usr/bin/env python3
"""
Reproduction of Edwards et al. 2001 Phenotype Phase Plane Analysis

Based on: Edwards, J. S., Ibarra, R. U., & Palsson, B. O. (2001). 
In silico predictions of Escherichia coli metabolic capabilities are consistent 
with experimental data. Nature Biotechnology, 19(2), 125-130.

This script reproduces the phenotype phase plane analysis for acetate vs oxygen
and succinate vs oxygen uptake as described in the original paper.
"""

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def load_and_setup_model():
    """Load the model and set up for phenotype phase plane analysis."""
    print("Loading E. coli iJO1366 model...")
    model = cobra.io.read_sbml_model("bigg_models/iJO1366.xml")
    print(f"Model loaded: {model.name}")
    print(f"Reactions: {len(model.reactions)}, Metabolites: {len(model.metabolites)}")
    
    # Get key exchange reactions
    acetate_rxn = model.reactions.get_by_id("EX_ac_e")
    succinate_rxn = model.reactions.get_by_id("EX_succ_e")
    oxygen_rxn = model.reactions.get_by_id("EX_o2_e")
    
    print(f"Key reactions: {acetate_rxn.id}, {succinate_rxn.id}, {oxygen_rxn.id}")
    
    return model, acetate_rxn, succinate_rxn, oxygen_rxn

def setup_minimal_medium(model, acetate_rxn, succinate_rxn, oxygen_rxn):
    """Set up minimal medium conditions as described in Edwards et al. 2001."""
    print("\nSetting up minimal medium...")
    
    # Allow essential nutrients for growth
    essential_exchanges = {
        "EX_nh4_e",    # Ammonium
        "EX_pi_e",     # Phosphate  
        "EX_so4_e",    # Sulfate
        "EX_k_e",      # Potassium
        "EX_mg2_e",    # Magnesium
        "EX_fe2_e",    # Iron
        "EX_ca2_e",    # Calcium
        "EX_cl_e",     # Chloride
        "EX_na1_e",    # Sodium
        "EX_h2o_e",    # Water
        "EX_h_e",      # Protons
        "EX_o2_e",     # Oxygen
        "EX_ac_e",     # Acetate
        "EX_succ_e",   # Succinate
    }
    
    # Fix all other carbon sources to zero
    fixed_count = 0
    for rxn in model.exchanges:
        if rxn.id not in essential_exchanges:
            rxn.lower_bound = 0.0
            fixed_count += 1
    
    print(f"Fixed {fixed_count} exchange reactions to zero")
    print(f"Allowed {len(essential_exchanges)} essential exchanges")
    
    return model

def compute_phenotype_phase_plane(model, substrate_rxn, substrate_name, 
                                 substrate_range, o2_range):
    """
    Compute phenotype phase plane for a substrate vs oxygen.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    substrate_rxn : cobra.Reaction
        The substrate exchange reaction
    substrate_name : str
        Name of the substrate
    substrate_range : array
        Range of substrate uptake rates
    o2_range : array
        Range of oxygen uptake rates
    """
    print(f"\nComputing phenotype phase plane for {substrate_name} vs O2...")
    print(f"Substrate range: {substrate_range[0]} to {substrate_range[-1]} mmol/gDW/h")
    print(f"O2 range: {o2_range[0]} to {o2_range[-1]} mmol/gDW/h")
    
    results = []
    total_points = len(substrate_range) * len(o2_range)
    current_point = 0
    
    for sub_flux in substrate_range:
        for o2_flux in o2_range:
            current_point += 1
            if current_point % 50 == 0:
                print(f"Progress: {current_point}/{total_points}")
            
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
    
    df = pd.DataFrame(results)
    
    # Filter for significant growth
    significant_growth = df[df.growth_rate > 1e-6]
    print(f"Total points: {len(df)}")
    print(f"Points with growth: {len(significant_growth)}")
    if len(significant_growth) > 0:
        print(f"Max growth rate: {significant_growth.growth_rate.max():.6f} 1/h")
    
    return df

def find_line_of_optimality(df, threshold=0.99):
    """
    Find the Line of Optimality (LO) for the phenotype phase plane.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Results from phenotype phase plane computation
    threshold : float
        Threshold for considering points as optimal
        
    Returns:
    --------
    tuple : (slope, intercept, lo_points)
    """
    # Filter for significant growth
    growth_data = df[df.growth_rate > 1e-6]
    
    if len(growth_data) == 0:
        print("No significant growth found")
        return 0.0, 0.0, pd.DataFrame()
    
    max_growth = growth_data.growth_rate.max()
    lo_points = growth_data[growth_data.growth_rate >= threshold * max_growth]
    
    print(f"Maximum growth rate: {max_growth:.6f} 1/h")
    print(f"Line of Optimality points: {len(lo_points)}")
    
    if len(lo_points) > 1:
        # Fit linear model: O2 = m * substrate + b
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            lo_points.substrate_flux, lo_points.o2_flux
        )
        
        print(f"LO slope: {slope:.3f} mmol O₂ / mmol substrate")
        print(f"LO intercept: {intercept:.3f} mmol O₂ / gDW/h")
        print(f"R²: {r_value**2:.3f}")
        
        return slope, intercept, lo_points
    else:
        print("Not enough points for Line of Optimality analysis")
        return 0.0, 0.0, lo_points

def create_visualization(df, slope, intercept, lo_points, substrate_name, save_path=None):
    """Create phenotype phase plane visualization."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Filter for growth
    growth_data = df[df.growth_rate > 1e-6]
    
    if len(growth_data) > 0:
        # Scatter plot colored by growth rate
        sc = ax.scatter(growth_data.substrate_flux, growth_data.o2_flux,
                       c=growth_data.growth_rate, cmap='viridis', alpha=0.7, s=30)
        
        # Add Line of Optimality
        if len(lo_points) > 1:
            x_range = np.array([lo_points.substrate_flux.min(), lo_points.substrate_flux.max()])
            y_range = slope * x_range + intercept
            ax.plot(x_range, y_range, 'r--', lw=2, 
                   label=f'LO slope = {slope:.2f}')
            ax.legend()
        
        plt.colorbar(sc, label='Growth rate (1/h)')
    else:
        # Show all points in gray if no growth
        ax.scatter(df.substrate_flux, df.o2_flux, c='gray', alpha=0.3, s=10)
    
    ax.set_xlabel(f'{substrate_name.capitalize()} uptake (mmol/gDW/h)')
    ax.set_ylabel('O₂ uptake (mmol/gDW/h)')
    ax.set_title(f'Phenotype Phase Plane: {substrate_name.capitalize()} vs O₂')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    plt.show()

def compare_with_experimental(slope_sim, intercept_sim, substrate_name):
    """
    Compare simulation results with experimental data from Edwards et al. 2001.
    """
    # Experimental data from Edwards et al. 2001
    experimental_data = {
        'acetate': {
            'slope': 0.67,  # mmol O2 / mmol acetate
            'intercept': -2.0,  # mmol O2 / gDW/h
        },
        'succinate': {
            'slope': 0.50,  # mmol O2 / mmol succinate  
            'intercept': -1.5,  # mmol O2 / gDW/h
        }
    }
    
    exp_data = experimental_data[substrate_name]
    exp_slope = exp_data['slope']
    exp_intercept = exp_data['intercept']
    
    print(f"\n=== COMPARISON WITH EDWARDS ET AL. 2001 ===")
    print(f"Substrate: {substrate_name}")
    print(f"Experimental slope: {exp_slope:.3f} mmol O₂/mmol {substrate_name}")
    print(f"Simulated slope: {slope_sim:.3f} mmol O₂/mmol {substrate_name}")
    print(f"Slope difference: {slope_sim - exp_slope:.3f}")
    print(f"Slope % difference: {((slope_sim - exp_slope) / exp_slope * 100):.1f}%")
    
    print(f"Experimental intercept: {exp_intercept:.3f} mmol O₂/gDW/h")
    print(f"Simulated intercept: {intercept_sim:.3f} mmol O₂/gDW/h")
    print(f"Intercept difference: {intercept_sim - exp_intercept:.3f}")
    print(f"Intercept % difference: {((intercept_sim - exp_intercept) / abs(exp_intercept) * 100):.1f}%")
    
    return {
        'substrate': substrate_name,
        'exp_slope': exp_slope,
        'sim_slope': slope_sim,
        'slope_diff': slope_sim - exp_slope,
        'slope_pct_diff': (slope_sim - exp_slope) / exp_slope * 100,
        'exp_intercept': exp_intercept,
        'sim_intercept': intercept_sim,
        'intercept_diff': intercept_sim - exp_intercept,
        'intercept_pct_diff': (intercept_sim - exp_intercept) / abs(exp_intercept) * 100
    }

def main():
    """Main function to reproduce Edwards et al. 2001 analysis."""
    print("="*60)
    print("REPRODUCTION OF EDWARDS ET AL. 2001")
    print("Phenotype Phase Plane Analysis")
    print("="*60)
    
    # Load and setup model
    model, acetate_rxn, succinate_rxn, oxygen_rxn = load_and_setup_model()
    model = setup_minimal_medium(model, acetate_rxn, succinate_rxn, oxygen_rxn)
    
    # Define ranges based on Edwards et al. 2001
    acetate_range = np.linspace(0, -15, 31)  # 0 to -15 mmol/gDW/h
    succinate_range = np.linspace(0, -15, 31)  # 0 to -15 mmol/gDW/h
    o2_range = np.linspace(0, -20, 41)  # 0 to -20 mmol/gDW/h
    
    results = {}
    
    # Acetate analysis
    print("\n" + "="*50)
    print("ACETATE VS OXYGEN ANALYSIS")
    print("="*50)
    
    acetate_df = compute_phenotype_phase_plane(model, acetate_rxn, "acetate", 
                                              acetate_range, o2_range)
    acetate_slope, acetate_intercept, acetate_lo = find_line_of_optimality(acetate_df)
    
    create_visualization(acetate_df, acetate_slope, acetate_intercept, acetate_lo, 
                        "acetate", "acetate_o2_edwards2001.png")
    
    acetate_comparison = compare_with_experimental(acetate_slope, acetate_intercept, "acetate")
    results['acetate'] = {
        'data': acetate_df,
        'slope': acetate_slope,
        'intercept': acetate_intercept,
        'lo_points': acetate_lo,
        'comparison': acetate_comparison
    }
    
    # Succinate analysis
    print("\n" + "="*50)
    print("SUCCINATE VS OXYGEN ANALYSIS")
    print("="*50)
    
    succinate_df = compute_phenotype_phase_plane(model, succinate_rxn, "succinate", 
                                                succinate_range, o2_range)
    succinate_slope, succinate_intercept, succinate_lo = find_line_of_optimality(succinate_df)
    
    create_visualization(succinate_df, succinate_slope, succinate_intercept, succinate_lo, 
                        "succinate", "succinate_o2_edwards2001.png")
    
    succinate_comparison = compare_with_experimental(succinate_slope, succinate_intercept, "succinate")
    results['succinate'] = {
        'data': succinate_df,
        'slope': succinate_slope,
        'intercept': succinate_intercept,
        'lo_points': succinate_lo,
        'comparison': succinate_comparison
    }
    
    # Save results
    acetate_df.to_csv("acetate_o2_edwards2001.csv", index=False)
    succinate_df.to_csv("succinate_o2_edwards2001.csv", index=False)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    for substrate, result in results.items():
        print(f"\n{substrate.upper()}:")
        print(f"  Simulated slope: {result['slope']:.3f}")
        print(f"  Experimental slope: {result['comparison']['exp_slope']:.3f}")
        print(f"  Slope difference: {result['comparison']['slope_diff']:.3f} ({result['comparison']['slope_pct_diff']:.1f}%)")
        print(f"  Simulated intercept: {result['intercept']:.3f}")
        print(f"  Experimental intercept: {result['comparison']['exp_intercept']:.3f}")
        print(f"  Intercept difference: {result['comparison']['intercept_diff']:.3f} ({result['comparison']['intercept_pct_diff']:.1f}%)")
    
    print(f"\nFiles saved:")
    print(f"  acetate_o2_edwards2001.csv")
    print(f"  succinate_o2_edwards2001.csv")
    print(f"  acetate_o2_edwards2001.png")
    print(f"  succinate_o2_edwards2001.png")
    
    print("\nReproduction complete!")

if __name__ == "__main__":
    main() 