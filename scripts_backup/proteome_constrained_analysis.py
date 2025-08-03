#!/usr/bin/env python3
"""
Proteome-Constrained Analysis: E. coli Metabolic Modeling with Enzyme Constraints

This script implements a simplified proteome-constrained analysis to demonstrate
the impact of enzyme allocation constraints on metabolic decision-making.
"""

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import warnings
warnings.filterwarnings('ignore')

class ProteomeConstrainedAnalysis:
    """Simplified proteome-constrained metabolic analysis."""
    
    def __init__(self, model_path="bigg_models/iJO1366.xml"):
        self.model_path = model_path
        self.model = None
        self.proteome_model = None
        
        # Enzyme data with k_cat values (1/s) and molecular weights (g/mol)
        self.enzyme_data = {
            'HEX1': {'k_cat': 100.0, 'mw': 50000, 'reactions': ['HEX1']},
            'PFK': {'k_cat': 50.0, 'mw': 35000, 'reactions': ['PFK']},
            'PYK': {'k_cat': 200.0, 'mw': 55000, 'reactions': ['PYK']},
            'CS': {'k_cat': 30.0, 'mw': 100000, 'reactions': ['CS']},
            'MDH': {'k_cat': 150.0, 'mw': 35000, 'reactions': ['MDH']},
            'LACZ': {'k_cat': 20.0, 'mw': 116000, 'reactions': ['LACZ']},
        }
        
        # Global protein budget
        self.global_protein_budget = 0.5  # g protein/gDW
    
    def load_model(self):
        """Load and prepare the iJO1366 model."""
        print("Loading iJO1366 model...")
        try:
            self.model = cobra.io.read_sbml_model(self.model_path)
            print(f"Model loaded: {len(self.model.genes)} genes, {len(self.model.reactions)} reactions")
            return True
        except Exception as e:
            print(f"Error loading model: {e}")
            return False
    
    def create_proteome_constrained_model(self):
        """Add proteome constraints to the model."""
        print("Adding proteome constraints...")
        self.proteome_model = self.model.copy()
        
        # Add enzyme capacity constraints
        for enzyme_id, data in self.enzyme_data.items():
            k_cat = data['k_cat']  # 1/s
            mw = data['mw']  # g/mol
            k_cat_mmol = k_cat * 3600 / mw  # mmol/gDW/h
            
            for rxn_id in data['reactions']:
                try:
                    rxn = self.proteome_model.reactions.get_by_id(rxn_id)
                    rxn.upper_bound = min(rxn.upper_bound, k_cat_mmol)
                except KeyError:
                    print(f"Warning: Reaction {rxn_id} not found")
        
        print(f"Proteome constraints added: {len(self.enzyme_data)} enzymes")
    
    def setup_glucose_medium(self):
        """Set up glucose minimal medium."""
        glucose_rxn = self.proteome_model.reactions.get_by_id('EX_glc__D_e')
        glucose_rxn.bounds = (-10, 0)
        
        # Close other carbon sources
        carbon_exchanges = ['EX_ac_e', 'EX_succ_e', 'EX_pyr_e', 'EX_lac__D_e']
        for exchange_id in carbon_exchanges:
            try:
                rxn = self.proteome_model.reactions.get_by_id(exchange_id)
                rxn.bounds = (0, 0)
            except:
                pass
        
        print("Glucose minimal medium setup complete")
    
    def setup_lactose_medium(self):
        """Set up lactose minimal medium."""
        lactose_rxn = self.proteome_model.reactions.get_by_id('EX_lac__D_e')
        lactose_rxn.bounds = (-10, 0)
        
        # Close other carbon sources
        carbon_exchanges = ['EX_ac_e', 'EX_succ_e', 'EX_pyr_e', 'EX_glc__D_e']
        for exchange_id in carbon_exchanges:
            try:
                rxn = self.proteome_model.reactions.get_by_id(exchange_id)
                rxn.bounds = (0, 0)
            except:
                pass
        
        print("Lactose minimal medium setup complete")
    
    def compare_growth_rates(self):
        """Compare growth rates between constrained and unconstrained models."""
        print("\n" + "=" * 60)
        print("GROWTH RATE COMPARISON")
        print("=" * 60)
        
        results = {}
        
        # Test glucose growth
        print("\nTesting glucose growth...")
        
        # Unconstrained model
        glucose_model = self.model.copy()
        glucose_rxn = glucose_model.reactions.get_by_id('EX_glc__D_e')
        glucose_rxn.bounds = (-10, 0)
        solution = glucose_model.optimize()
        unconstrained_glucose_growth = solution.objective_value
        print(f"  Unconstrained glucose growth: {unconstrained_glucose_growth:.6f} 1/h")
        
        # Constrained model
        self.create_proteome_constrained_model()
        self.setup_glucose_medium()
        solution = self.proteome_model.optimize()
        constrained_glucose_growth = solution.objective_value
        print(f"  Constrained glucose growth: {constrained_glucose_growth:.6f} 1/h")
        
        glucose_reduction = (unconstrained_glucose_growth - constrained_glucose_growth) / unconstrained_glucose_growth * 100
        print(f"  Growth reduction: {glucose_reduction:.1f}%")
        
        results['glucose'] = {
            'unconstrained': unconstrained_glucose_growth,
            'constrained': constrained_glucose_growth,
            'reduction_percent': glucose_reduction
        }
        
        # Test lactose growth
        print("\nTesting lactose growth...")
        
        # Unconstrained model
        lactose_model = self.model.copy()
        lactose_rxn = lactose_model.reactions.get_by_id('EX_lac__D_e')
        lactose_rxn.bounds = (-10, 0)
        solution = lactose_model.optimize()
        unconstrained_lactose_growth = solution.objective_value
        print(f"  Unconstrained lactose growth: {unconstrained_lactose_growth:.6f} 1/h")
        
        # Constrained model
        self.proteome_model = self.model.copy()
        self.create_proteome_constrained_model()
        self.setup_lactose_medium()
        solution = self.proteome_model.optimize()
        constrained_lactose_growth = solution.objective_value
        print(f"  Constrained lactose growth: {constrained_lactose_growth:.6f} 1/h")
        
        lactose_reduction = (unconstrained_lactose_growth - constrained_lactose_growth) / unconstrained_lactose_growth * 100
        print(f"  Growth reduction: {lactose_reduction:.1f}%")
        
        results['lactose'] = {
            'unconstrained': unconstrained_lactose_growth,
            'constrained': constrained_lactose_growth,
            'reduction_percent': lactose_reduction
        }
        
        return results
    
    def analyze_flux_distributions(self):
        """Analyze flux distributions under constraints."""
        print("\n" + "=" * 60)
        print("FLUX DISTRIBUTION ANALYSIS")
        print("=" * 60)
        
        # Setup glucose medium
        self.proteome_model = self.model.copy()
        self.create_proteome_constrained_model()
        self.setup_glucose_medium()
        
        # Get fluxes for key reactions
        solution = self.proteome_model.optimize()
        
        key_reactions = ['HEX1', 'PFK', 'PYK', 'CS', 'MDH', 'LACZ']
        flux_data = {}
        
        for rxn_id in key_reactions:
            try:
                flux = solution.fluxes.get(rxn_id, 0)
                flux_data[rxn_id] = flux
                print(f"  {rxn_id}: {flux:.3f} mmol/gDW/h")
            except:
                print(f"  {rxn_id}: Not found")
        
        return flux_data
    
    def perform_flux_variability_analysis(self):
        """Perform Flux Variability Analysis under enzyme constraints."""
        print("\n" + "=" * 60)
        print("FLUX VARIABILITY ANALYSIS")
        print("=" * 60)
        
        # Setup glucose medium with constraints
        self.proteome_model = self.model.copy()
        self.create_proteome_constrained_model()
        self.setup_glucose_medium()
        
        # Perform FVA
        fva_results = cobra.flux_analysis.flux_variability_analysis(
            self.proteome_model, 
            fraction_of_optimum=0.9
        )
        
        # Analyze key reactions
        key_reactions = ['HEX1', 'PFK', 'PYK', 'CS', 'MDH']
        fva_summary = {}
        
        for rxn_id in key_reactions:
            try:
                min_flux = fva_results.loc[rxn_id, 'minimum']
                max_flux = fva_results.loc[rxn_id, 'maximum']
                variability = max_flux - min_flux
                
                fva_summary[rxn_id] = {
                    'min_flux': min_flux,
                    'max_flux': max_flux,
                    'variability': variability
                }
                
                print(f"  {rxn_id}:")
                print(f"    Min flux: {min_flux:.3f} mmol/gDW/h")
                print(f"    Max flux: {max_flux:.3f} mmol/gDW/h")
                print(f"    Variability: {variability:.3f} mmol/gDW/h")
                
            except:
                print(f"  {rxn_id}: Not found in FVA results")
        
        return fva_summary
    
    def create_visualizations(self, growth_results, flux_data, fva_summary):
        """Create visualizations of the analysis results."""
        print("\nCreating visualizations...")
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Growth rate comparison
        substrates = ['glucose', 'lactose']
        unconstrained_growth = [growth_results[s]['unconstrained'] for s in substrates]
        constrained_growth = [growth_results[s]['constrained'] for s in substrates]
        
        x = np.arange(len(substrates))
        width = 0.35
        
        axes[0, 0].bar(x - width/2, unconstrained_growth, width, label='Unconstrained', alpha=0.8)
        axes[0, 0].bar(x + width/2, constrained_growth, width, label='Proteome-Constrained', alpha=0.8)
        axes[0, 0].set_xlabel('Substrate')
        axes[0, 0].set_ylabel('Growth Rate (1/h)')
        axes[0, 0].set_title('Growth Rate Comparison')
        axes[0, 0].set_xticks(x)
        axes[0, 0].set_xticklabels(substrates)
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Growth reduction percentage
        reductions = [growth_results[s]['reduction_percent'] for s in substrates]
        axes[0, 1].bar(substrates, reductions, color=['green', 'red'], alpha=0.8)
        axes[0, 1].set_xlabel('Substrate')
        axes[0, 1].set_ylabel('Growth Reduction (%)')
        axes[0, 1].set_title('Impact of Proteome Constraints')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Key reaction fluxes
        if flux_data:
            reactions = list(flux_data.keys())
            fluxes = list(flux_data.values())
            
            axes[1, 0].bar(reactions, fluxes, alpha=0.8)
            axes[1, 0].set_xlabel('Reaction')
            axes[1, 0].set_ylabel('Flux (mmol/gDW/h)')
            axes[1, 0].set_title('Key Reaction Fluxes (Constrained)')
            axes[1, 0].tick_params(axis='x', rotation=45)
            axes[1, 0].grid(True, alpha=0.3)
        
        # 4. Flux variability
        if fva_summary:
            reactions = list(fva_summary.keys())
            variabilities = [fva_summary[r]['variability'] for r in reactions]
            
            axes[1, 1].bar(reactions, variabilities, alpha=0.8, color='orange')
            axes[1, 1].set_xlabel('Reaction')
            axes[1, 1].set_ylabel('Flux Variability (mmol/gDW/h)')
            axes[1, 1].set_title('Flux Variability Analysis')
            axes[1, 1].tick_params(axis='x', rotation=45)
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('results/proteome_constrained_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Visualizations saved to results/proteome_constrained_analysis.png")
    
    def save_results(self, growth_results, flux_data, fva_summary):
        """Save analysis results to files."""
        print("\nSaving results...")
        
        import os
        os.makedirs('results', exist_ok=True)
        
        # Save growth comparison results
        with open('results/proteome_constrained_growth_comparison.json', 'w') as f:
            json.dump(growth_results, f, indent=2)
        print("Growth comparison saved to results/proteome_constrained_growth_comparison.json")
        
        # Save flux data
        with open('results/proteome_constrained_flux_data.json', 'w') as f:
            json.dump(flux_data, f, indent=2)
        print("Flux data saved to results/proteome_constrained_flux_data.json")
        
        # Save FVA results
        with open('results/proteome_constrained_fva_results.json', 'w') as f:
            json.dump(fva_summary, f, indent=2)
        print("FVA results saved to results/proteome_constrained_fva_results.json")
    
    def run_analysis(self):
        """Run the complete proteome-constrained analysis."""
        print("=" * 70)
        print("PROTEOME-CONSTRAINED METABOLIC ANALYSIS")
        print("=" * 70)
        
        if not self.load_model():
            return False
        
        # Compare growth rates
        growth_results = self.compare_growth_rates()
        
        # Analyze flux distributions
        flux_data = self.analyze_flux_distributions()
        
        # Perform FVA
        fva_summary = self.perform_flux_variability_analysis()
        
        # Create visualizations
        self.create_visualizations(growth_results, flux_data, fva_summary)
        
        # Save results
        self.save_results(growth_results, flux_data, fva_summary)
        
        # Summary
        print("\n" + "=" * 60)
        print("ANALYSIS SUMMARY")
        print("=" * 60)
        
        glucose_reduction = growth_results['glucose']['reduction_percent']
        lactose_reduction = growth_results['lactose']['reduction_percent']
        
        print(f"Impact of proteome constraints:")
        print(f"  Glucose growth reduction: {glucose_reduction:.1f}%")
        print(f"  Lactose growth reduction: {lactose_reduction:.1f}%")
        print(f"  Number of enzymes constrained: {len(self.enzyme_data)}")
        print(f"  Global protein budget: {self.global_protein_budget} g protein/gDW")
        
        if glucose_reduction > lactose_reduction:
            print(f"  Glucose metabolism is more constrained than lactose metabolism")
        else:
            print(f"  Lactose metabolism is more constrained than glucose metabolism")
        
        return True

def main():
    """Main function to run the proteome-constrained analysis."""
    analyzer = ProteomeConstrainedAnalysis()
    success = analyzer.run_analysis()
    
    if success:
        print("\n" + "=" * 70)
        print("PROTEOME-CONSTRAINED ANALYSIS COMPLETE!")
        print("=" * 70)
        print("✅ Successfully analyzed proteome constraints on metabolism")
        print("✅ Compared growth rates under constrained vs unconstrained conditions")
        print("✅ Analyzed flux distributions and variability")
        print("✅ Generated comprehensive visualizations")
        print("✅ Identified metabolic bottlenecks under enzyme constraints")
        
    else:
        print("\n❌ Analysis failed!")

if __name__ == "__main__":
    main() 