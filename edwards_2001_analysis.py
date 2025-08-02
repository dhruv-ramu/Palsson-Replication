#!/usr/bin/env python3
"""
Reproduction of Edwards et al. 2001 Phenotype Phase Plane Analysis

This script reproduces the phenotype phase plane analysis from:
Edwards, J. S., Ibarra, R. U., & Palsson, B. O. (2001). 
In silico predictions of Escherichia coli metabolic capabilities are consistent 
with experimental data. Nature Biotechnology, 19(2), 125-130.

Author: Computational Biology Analysis
Date: 2024
"""

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, List, Optional, Dict
import logging
import warnings
from scipy import stats

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class Edwards2001Reproducer:
    """
    A class to reproduce the Edwards et al. 2001 phenotype phase plane analysis.
    """
    
    def __init__(self, model_path: str):
        """
        Initialize the reproducer with a metabolic model.
        
        Parameters:
        -----------
        model_path : str
            Path to the SBML model file
        """
        self.model_path = model_path
        self.model = None
        self.results = {}
        
        # Experimental data from Edwards et al. 2001
        self.experimental_data = {
            'acetate_o2': {
                'slope': 0.67,  # mmol O2 / mmol acetate
                'intercept': -2.0,  # mmol O2 / gDW/h
                'reference': 'Edwards et al. 2001'
            },
            'succinate_o2': {
                'slope': 0.50,  # mmol O2 / mmol succinate
                'intercept': -1.5,  # mmol O2 / gDW/h
                'reference': 'Edwards et al. 2001'
            }
        }
        
    def load_model(self) -> bool:
        """
        Load the metabolic model from SBML file.
        
        Returns:
        --------
        bool
            True if model loaded successfully, False otherwise
        """
        try:
            logger.info(f"Loading model from {self.model_path}")
            self.model = cobra.io.read_sbml_model(self.model_path)
            logger.info(f"Model loaded successfully: {self.model.name}")
            logger.info(f"Number of reactions: {len(self.model.reactions)}")
            logger.info(f"Number of metabolites: {len(self.model.metabolites)}")
            return True
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            return False
    
    def setup_model_for_analysis(self) -> bool:
        """
        Set up the model for phenotype phase plane analysis.
        
        Returns:
        --------
        bool
            True if setup successful, False otherwise
        """
        try:
            # Get exchange reactions
            self.ac_rxn = self.model.reactions.get_by_id("EX_ac_e")
            self.succ_rxn = self.model.reactions.get_by_id("EX_succ_e")
            self.o2_rxn = self.model.reactions.get_by_id("EX_o2_e")
            
            logger.info(f"Acetate exchange reaction: {self.ac_rxn.id}")
            logger.info(f"Succinate exchange reaction: {self.succ_rxn.id}")
            logger.info(f"Oxygen exchange reaction: {self.o2_rxn.id}")
            
            # Fix all other carbon sources to zero
            fixed_count = 0
            for rxn in self.model.exchanges:
                if rxn.id not in {self.ac_rxn.id, self.succ_rxn.id, self.o2_rxn.id}:
                    rxn.lower_bound = 0.0
                    fixed_count += 1
            
            logger.info(f"Fixed {fixed_count} exchange reactions to zero")
            return True
            
        except Exception as e:
            logger.error(f"Failed to set up model: {e}")
            return False
    
    def compute_phenotype_phase_plane(self, 
                                    substrate_rxn, 
                                    substrate_limits: np.ndarray, 
                                    o2_limits: np.ndarray,
                                    substrate_name: str) -> pd.DataFrame:
        """
        Compute phenotype phase plane for a substrate vs oxygen uptake.
        
        Parameters:
        -----------
        substrate_rxn : cobra.Reaction
            The substrate exchange reaction
        substrate_limits : np.ndarray
            Array of substrate uptake rates to test
        o2_limits : np.ndarray
            Array of oxygen uptake rates to test
        substrate_name : str
            Name of the substrate for logging
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with results of the phenotype phase plane analysis
        """
        logger.info(f"Computing phenotype phase plane for {substrate_name} vs O2...")
        results = []
        total_points = len(substrate_limits) * len(o2_limits)
        current_point = 0
        
        for sub_flux in substrate_limits:
            for o2_flux in o2_limits:
                current_point += 1
                if current_point % 100 == 0:
                    logger.info(f"Progress: {current_point}/{total_points} points processed")
                
                # Set the exchange reaction bounds
                substrate_rxn.lower_bound = sub_flux
                self.o2_rxn.lower_bound = o2_flux
                
                try:
                    # Optimize the model
                    solution = self.model.optimize()
                    
                    if solution.status == 'optimal':
                        # Use a small threshold to distinguish between zero and very small values
                        growth_rate = solution.objective_value if solution.objective_value > 1e-10 else 0.0
                        results.append({
                            'substrate_flux': sub_flux,
                            'o2_flux': o2_flux,
                            'growth_rate': growth_rate,
                            'status': solution.status
                        })
                    else:
                        results.append({
                            'substrate_flux': sub_flux,
                            'o2_flux': o2_flux,
                            'growth_rate': 0.0,
                            'status': solution.status
                        })
                except Exception as e:
                    logger.warning(f"Optimization failed for {substrate_name}={sub_flux}, O2={o2_flux}: {e}")
                    results.append({
                        'substrate_flux': sub_flux,
                        'o2_flux': o2_flux,
                        'growth_rate': 0.0,
                        'status': 'error'
                    })
        
        df = pd.DataFrame(results)
        logger.info(f"Phenotype phase plane computed with shape: {df.shape}")
        
        # Filter out very small growth rates for better analysis
        significant_growth = df[df.growth_rate > 1e-10]
        logger.info(f"Points with significant growth (>1e-10): {len(significant_growth)}")
        
        if len(significant_growth) > 0:
            logger.info(f"Growth rate range: {significant_growth.growth_rate.min():.6f} to {significant_growth.growth_rate.max():.6f} 1/h")
        else:
            logger.warning("No significant growth detected in any condition")
        
        return df
    
    def find_line_of_optimality(self, results_df: pd.DataFrame, threshold: float = 0.99) -> Tuple[float, float, pd.DataFrame]:
        """
        Find the Line of Optimality (LO) for the phenotype phase plane.
        
        Parameters:
        -----------
        results_df : pd.DataFrame
            Results from phenotype phase plane computation
        threshold : float
            Threshold for considering points as optimal (default: 0.99)
            
        Returns:
        --------
        Tuple[float, float, pd.DataFrame]
            Slope, intercept, and points on the line of optimality
        """
        # Filter out zero growth rates
        non_zero_growth = results_df[results_df.growth_rate > 1e-10]
        
        if len(non_zero_growth) == 0:
            logger.warning("No non-zero growth rates found")
            return 0.0, 0.0, pd.DataFrame()
        
        max_growth = non_zero_growth.growth_rate.max()
        lo_points = non_zero_growth[non_zero_growth.growth_rate >= threshold * max_growth]
        
        logger.info(f"Maximum growth rate: {max_growth:.6f} 1/h")
        logger.info(f"Number of LO points: {len(lo_points)}")
        
        if len(lo_points) > 1:
            try:
                # Fit a linear model: O2 = m * substrate + b
                m, b = np.polyfit(lo_points.substrate_flux, lo_points.o2_flux, 1)
                logger.info(f"LO slope: {m:.3f} mmol O₂ / mmol substrate")
                logger.info(f"LO intercept: {b:.3f} mmol O₂ / gDW/h")
                return m, b, lo_points
            except Exception as e:
                logger.error(f"Error fitting line of optimality: {e}")
                return 0.0, 0.0, lo_points
        else:
            logger.warning("Not enough points for Line of Optimality analysis")
            return 0.0, 0.0, lo_points
    
    def create_visualization(self, 
                           results_df: pd.DataFrame,
                           m: float, 
                           b: float, 
                           lo_points: pd.DataFrame,
                           substrate_name: str,
                           save_path: Optional[str] = None) -> None:
        """
        Create and display the phenotype phase plane visualization.
        
        Parameters:
        -----------
        results_df : pd.DataFrame
            Results from phenotype phase plane computation
        m : float
            Slope of the line of optimality
        b : float
            Intercept of the line of optimality
        lo_points : pd.DataFrame
            Points on the line of optimality
        substrate_name : str
            Name of the substrate
        save_path : str, optional
            Path to save the plot (if provided)
        """
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Filter out zero growth rates for better visualization
        plot_data = results_df[results_df.growth_rate > 1e-10]
        
        if len(plot_data) > 0:
            # Scatter plot of all points colored by growth rate
            sc = ax.scatter(
                plot_data.substrate_flux, plot_data.o2_flux,
                c=plot_data.growth_rate, cmap='viridis', alpha=0.7, s=30
            )
            
            # Add the Line of Optimality if we have enough points
            if len(lo_points) > 1:
                ax.plot(
                    lo_points.substrate_flux,
                    m * lo_points.substrate_flux + b,
                    'r--', lw=2, label=f"LO slope = {m:.2f}"
                )
                ax.legend()
        else:
            # If no growth, just show the grid
            ax.scatter(results_df.substrate_flux, results_df.o2_flux, 
                      c='gray', alpha=0.3, s=10, label='No growth')
            ax.legend()
        
        # Plot aesthetics
        ax.set_xlabel(f"{substrate_name.capitalize()} uptake (mmol·gDW⁻¹·h⁻¹)")
        ax.set_ylabel("O₂ uptake (mmol·gDW⁻¹·h⁻¹)")
        ax.set_title(f"Phenotype Phase Plane: {substrate_name.capitalize()} vs. O₂")
        
        if len(plot_data) > 0:
            fig.colorbar(sc, label="Growth rate (1/h)")
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")
        
        plt.show()
    
    def compare_with_experimental(self, 
                                substrate_name: str, 
                                m_sim: float, 
                                b_sim: float) -> Dict:
        """
        Compare simulation results with experimental data.
        
        Parameters:
        -----------
        substrate_name : str
            Name of the substrate ('acetate' or 'succinate')
        m_sim : float
            Simulated slope
        b_sim : float
            Simulated intercept
            
        Returns:
        --------
        Dict
            Comparison results
        """
        key = f"{substrate_name}_o2"
        exp_data = self.experimental_data[key]
        
        m_exp = exp_data['slope']
        b_exp = exp_data['intercept']
        
        # Calculate differences
        slope_diff = m_sim - m_exp
        intercept_diff = b_sim - b_exp
        
        # Calculate percentage differences
        slope_pct_diff = (slope_diff / m_exp) * 100 if m_exp != 0 else float('inf')
        intercept_pct_diff = (intercept_diff / abs(b_exp)) * 100 if b_exp != 0 else float('inf')
        
        comparison = {
            'substrate': substrate_name,
            'experimental_slope': m_exp,
            'simulated_slope': m_sim,
            'slope_difference': slope_diff,
            'slope_percentage_difference': slope_pct_diff,
            'experimental_intercept': b_exp,
            'simulated_intercept': b_sim,
            'intercept_difference': intercept_diff,
            'intercept_percentage_difference': intercept_pct_diff,
            'reference': exp_data['reference']
        }
        
        return comparison
    
    def run_acetate_analysis(self) -> Dict:
        """
        Run phenotype phase plane analysis for acetate vs oxygen.
        
        Returns:
        --------
        Dict
            Analysis results
        """
        logger.info("="*60)
        logger.info("ACETATE VS OXYGEN ANALYSIS")
        logger.info("="*60)
        
        # Define flux limits based on Edwards et al. 2001
        ac_limits = np.linspace(0, -15, 31)   # Acetate uptake: 0 to -15 mmol/gDW/h
        o2_limits = np.linspace(0, -20, 41)   # O2 uptake: 0 to -20 mmol/gDW/h
        
        # Compute phenotype phase plane
        results = self.compute_phenotype_phase_plane(self.ac_rxn, ac_limits, o2_limits, "acetate")
        
        # Find line of optimality
        m, b, lo_points = self.find_line_of_optimality(results)
        
        # Create visualization
        self.create_visualization(results, m, b, lo_points, "acetate", 
                                save_path="acetate_o2_phpp.png")
        
        # Compare with experimental data
        comparison = self.compare_with_experimental("acetate", m, b)
        
        # Store results
        self.results['acetate'] = {
            'data': results,
            'slope': m,
            'intercept': b,
            'lo_points': lo_points,
            'comparison': comparison
        }
        
        return comparison
    
    def run_succinate_analysis(self) -> Dict:
        """
        Run phenotype phase plane analysis for succinate vs oxygen.
        
        Returns:
        --------
        Dict
            Analysis results
        """
        logger.info("="*60)
        logger.info("SUCCINATE VS OXYGEN ANALYSIS")
        logger.info("="*60)
        
        # Define flux limits based on Edwards et al. 2001
        succ_limits = np.linspace(0, -15, 31)   # Succinate uptake: 0 to -15 mmol/gDW/h
        o2_limits = np.linspace(0, -20, 41)     # O2 uptake: 0 to -20 mmol/gDW/h
        
        # Compute phenotype phase plane
        results = self.compute_phenotype_phase_plane(self.succ_rxn, succ_limits, o2_limits, "succinate")
        
        # Find line of optimality
        m, b, lo_points = self.find_line_of_optimality(results)
        
        # Create visualization
        self.create_visualization(results, m, b, lo_points, "succinate", 
                                save_path="succinate_o2_phpp.png")
        
        # Compare with experimental data
        comparison = self.compare_with_experimental("succinate", m, b)
        
        # Store results
        self.results['succinate'] = {
            'data': results,
            'slope': m,
            'intercept': b,
            'lo_points': lo_points,
            'comparison': comparison
        }
        
        return comparison
    
    def save_results(self) -> None:
        """
        Save all results to CSV files.
        """
        for substrate, result in self.results.items():
            filename = f"{substrate}_o2_phpp.csv"
            result['data'].to_csv(filename, index=False)
            logger.info(f"Results saved to {filename}")
    
    def print_summary(self) -> None:
        """
        Print comprehensive summary of all analyses.
        """
        print("\n" + "="*80)
        print("EDWARDS ET AL. 2001 REPRODUCTION SUMMARY")
        print("="*80)
        
        for substrate, result in self.results.items():
            print(f"\n{substrate.upper()} VS OXYGEN ANALYSIS:")
            print("-" * 50)
            
            data = result['data']
            comparison = result['comparison']
            
            print(f"Total points analyzed: {len(data)}")
            significant_growth = data[data.growth_rate > 1e-10]
            print(f"Points with significant growth: {len(significant_growth)}")
            
            if len(significant_growth) > 0:
                print(f"Maximum growth rate: {significant_growth.growth_rate.max():.6f} 1/h")
                print(f"Line of Optimality slope: {result['slope']:.3f} mmol O₂ / mmol {substrate}")
                print(f"Line of Optimality intercept: {result['intercept']:.3f} mmol O₂ / gDW/h")
            
            print(f"\nCOMPARISON WITH EXPERIMENTAL DATA:")
            print(f"Experimental slope: {comparison['experimental_slope']:.3f}")
            print(f"Simulated slope: {comparison['simulated_slope']:.3f}")
            print(f"Slope difference: {comparison['slope_difference']:.3f}")
            print(f"Slope percentage difference: {comparison['slope_percentage_difference']:.1f}%")
            print(f"Experimental intercept: {comparison['experimental_intercept']:.3f}")
            print(f"Simulated intercept: {comparison['simulated_intercept']:.3f}")
            print(f"Intercept difference: {comparison['intercept_difference']:.3f}")
            print(f"Intercept percentage difference: {comparison['intercept_percentage_difference']:.1f}%")


def main():
    """
    Main function to reproduce Edwards et al. 2001 analysis.
    """
    # Configuration
    MODEL_PATH = "bigg_models/iJO1366.xml"
    
    # Initialize reproducer
    reproducer = Edwards2001Reproducer(MODEL_PATH)
    
    # Load model
    if not reproducer.load_model():
        logger.error("Failed to load model. Exiting.")
        return
    
    # Set up model
    if not reproducer.setup_model_for_analysis():
        logger.error("Failed to set up model. Exiting.")
        return
    
    # Run acetate analysis
    acetate_comparison = reproducer.run_acetate_analysis()
    
    # Run succinate analysis
    succinate_comparison = reproducer.run_succinate_analysis()
    
    # Save results
    reproducer.save_results()
    
    # Print summary
    reproducer.print_summary()
    
    logger.info("Edwards et al. 2001 reproduction completed successfully!")


if __name__ == "__main__":
    main() 