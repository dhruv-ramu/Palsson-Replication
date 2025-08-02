#!/usr/bin/env python3
"""
Phenotype Phase Plane Analysis: Acetate vs Oxygen

This script performs phenotype phase plane analysis to examine the relationship 
between acetate uptake and oxygen consumption in E. coli iJO1366.

Author: Computational Biology Analysis
Date: 2024
"""

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, List, Optional
import logging
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PhenotypePhasePlaneAnalyzer:
    """
    A class to perform phenotype phase plane analysis on metabolic models.
    """
    
    def __init__(self, model_path: str):
        """
        Initialize the analyzer with a metabolic model.
        
        Parameters:
        -----------
        model_path : str
            Path to the SBML model file
        """
        self.model_path = model_path
        self.model = None
        self.results = None
        self.ac_rxn = None
        self.o2_rxn = None
        self.succ_rxn = None
        
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
    
    def setup_exchange_reactions(self) -> bool:
        """
        Set up the key exchange reactions for analysis.
        
        Returns:
        --------
        bool
            True if reactions found successfully, False otherwise
        """
        try:
            # Define the key exchange reactions
            self.ac_rxn = self.model.reactions.get_by_id("EX_ac_e")
            self.succ_rxn = self.model.reactions.get_by_id("EX_succ_e")
            self.o2_rxn = self.model.reactions.get_by_id("EX_o2_e")
            
            logger.info(f"Acetate exchange reaction: {self.ac_rxn.id}")
            logger.info(f"Succinate exchange reaction: {self.succ_rxn.id}")
            logger.info(f"Oxygen exchange reaction: {self.o2_rxn.id}")
            return True
        except Exception as e:
            logger.error(f"Failed to set up exchange reactions: {e}")
            return False
    
    def fix_carbon_sources(self) -> int:
        """
        Fix all other carbon sources to zero, keeping only acetate, succinate, and oxygen.
        
        Returns:
        --------
        int
            Number of exchange reactions fixed to zero
        """
        fixed_count = 0
        for rxn in self.model.exchanges:
            if rxn.id not in {self.ac_rxn.id, self.succ_rxn.id, self.o2_rxn.id}:
                rxn.lower_bound = 0.0
                fixed_count += 1
        
        logger.info(f"Fixed {fixed_count} exchange reactions to zero")
        return fixed_count
    
    def test_model_feasibility(self) -> bool:
        """
        Test if the model can grow under the current constraints.
        
        Returns:
        --------
        bool
            True if model can grow, False otherwise
        """
        try:
            # Test with some reasonable bounds
            self.ac_rxn.lower_bound = -10.0
            self.o2_rxn.lower_bound = -10.0
            
            solution = self.model.optimize()
            
            if solution.status == 'optimal' and solution.objective_value > 1e-6:
                logger.info(f"Model can grow with objective value: {solution.objective_value:.6f}")
                return True
            else:
                logger.warning(f"Model cannot grow under current constraints. Status: {solution.status}")
                return False
        except Exception as e:
            logger.error(f"Error testing model feasibility: {e}")
            return False
    
    def compute_phenotype_phase_plane(self, 
                                    ac_limits: np.ndarray, 
                                    o2_limits: np.ndarray) -> pd.DataFrame:
        """
        Compute phenotype phase plane for acetate vs oxygen uptake.
        
        Parameters:
        -----------
        ac_limits : np.ndarray
            Array of acetate uptake rates to test
        o2_limits : np.ndarray
            Array of oxygen uptake rates to test
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with results of the phenotype phase plane analysis
        """
        logger.info("Computing phenotype phase plane...")
        results = []
        total_points = len(ac_limits) * len(o2_limits)
        current_point = 0
        
        for ac_flux in ac_limits:
            for o2_flux in o2_limits:
                current_point += 1
                if current_point % 100 == 0:
                    logger.info(f"Progress: {current_point}/{total_points} points processed")
                
                # Set the exchange reaction bounds
                self.ac_rxn.lower_bound = ac_flux
                self.o2_rxn.lower_bound = o2_flux
                
                try:
                    # Optimize the model
                    solution = self.model.optimize()
                    
                    if solution.status == 'optimal':
                        # Use a small threshold to distinguish between zero and very small values
                        growth_rate = solution.objective_value if solution.objective_value > 1e-10 else 0.0
                        results.append({
                            'EX_ac_e': ac_flux,
                            'EX_o2_e': o2_flux,
                            'growth_rate': growth_rate,
                            'status': solution.status
                        })
                    else:
                        results.append({
                            'EX_ac_e': ac_flux,
                            'EX_o2_e': o2_flux,
                            'growth_rate': 0.0,
                            'status': solution.status
                        })
                except Exception as e:
                    logger.warning(f"Optimization failed for ac={ac_flux}, o2={o2_flux}: {e}")
                    results.append({
                        'EX_ac_e': ac_flux,
                        'EX_o2_e': o2_flux,
                        'growth_rate': 0.0,
                        'status': 'error'
                    })
        
        self.results = pd.DataFrame(results)
        logger.info(f"Phenotype phase plane computed with shape: {self.results.shape}")
        
        # Filter out very small growth rates for better analysis
        significant_growth = self.results[self.results.growth_rate > 1e-10]
        logger.info(f"Points with significant growth (>1e-10): {len(significant_growth)}")
        
        if len(significant_growth) > 0:
            logger.info(f"Growth rate range: {significant_growth.growth_rate.min():.6f} to {significant_growth.growth_rate.max():.6f} 1/h")
        else:
            logger.warning("No significant growth detected in any condition")
        
        return self.results
    
    def find_line_of_optimality(self, threshold: float = 0.99) -> Tuple[float, float, pd.DataFrame]:
        """
        Find the Line of Optimality (LO) for the phenotype phase plane.
        
        Parameters:
        -----------
        threshold : float
            Threshold for considering points as optimal (default: 0.99)
            
        Returns:
        --------
        Tuple[float, float, pd.DataFrame]
            Slope, intercept, and points on the line of optimality
        """
        # Filter out zero growth rates
        non_zero_growth = self.results[self.results.growth_rate > 1e-10]
        
        if len(non_zero_growth) == 0:
            logger.warning("No non-zero growth rates found")
            return 0.0, 0.0, pd.DataFrame()
        
        max_growth = non_zero_growth.growth_rate.max()
        lo_points = non_zero_growth[non_zero_growth.growth_rate >= threshold * max_growth]
        
        logger.info(f"Maximum growth rate: {max_growth:.6f} 1/h")
        logger.info(f"Number of LO points: {len(lo_points)}")
        
        if len(lo_points) > 1:
            try:
                # Fit a linear model: O2 = m * Ac + b
                m_ac_o2, b_ac_o2 = np.polyfit(lo_points.EX_ac_e, lo_points.EX_o2_e, 1)
                logger.info(f"Acetate–O₂ LO slope: {m_ac_o2:.3f} mmol O₂ / mmol acetate")
                logger.info(f"Acetate–O₂ LO intercept: {b_ac_o2:.3f} mmol O₂ / gDW/h")
                return m_ac_o2, b_ac_o2, lo_points
            except Exception as e:
                logger.error(f"Error fitting line of optimality: {e}")
                return 0.0, 0.0, lo_points
        else:
            logger.warning("Not enough points for Line of Optimality analysis")
            return 0.0, 0.0, lo_points
    
    def create_visualization(self, 
                           m_ac_o2: float, 
                           b_ac_o2: float, 
                           lo_points: pd.DataFrame,
                           save_path: Optional[str] = None) -> None:
        """
        Create and display the phenotype phase plane visualization.
        
        Parameters:
        -----------
        m_ac_o2 : float
            Slope of the line of optimality
        b_ac_o2 : float
            Intercept of the line of optimality
        lo_points : pd.DataFrame
            Points on the line of optimality
        save_path : str, optional
            Path to save the plot (if provided)
        """
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Filter out zero growth rates for better visualization
        plot_data = self.results[self.results.growth_rate > 1e-10]
        
        if len(plot_data) > 0:
            # Scatter plot of all points colored by growth rate
            sc = ax.scatter(
                plot_data.EX_ac_e, plot_data.EX_o2_e,
                c=plot_data.growth_rate, cmap='viridis', alpha=0.7, s=30
            )
            
            # Add the Line of Optimality if we have enough points
            if len(lo_points) > 1:
                ax.plot(
                    lo_points.EX_ac_e,
                    m_ac_o2 * lo_points.EX_ac_e + b_ac_o2,
                    'r--', lw=2, label=f"LO slope = {m_ac_o2:.2f}"
                )
                ax.legend()
        else:
            # If no growth, just show the grid
            ax.scatter(self.results.EX_ac_e, self.results.EX_o2_e, 
                      c='gray', alpha=0.3, s=10, label='No growth')
            ax.legend()
        
        # Plot aesthetics
        ax.set_xlabel("Acetate uptake (mmol·gDW⁻¹·h⁻¹)")
        ax.set_ylabel("O₂ uptake (mmol·gDW⁻¹·h⁻¹)")
        ax.set_title("Phenotype Phase Plane: Acetate vs. O₂")
        
        if len(plot_data) > 0:
            fig.colorbar(sc, label="Growth rate (1/h)")
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")
        
        plt.show()
    
    def print_summary_statistics(self, m_ac_o2: float, b_ac_o2: float, lo_points: pd.DataFrame) -> None:
        """
        Print summary statistics of the analysis.
        
        Parameters:
        -----------
        m_ac_o2 : float
            Slope of the line of optimality
        b_ac_o2 : float
            Intercept of the line of optimality
        lo_points : pd.DataFrame
            Points on the line of optimality
        """
        print("\n" + "="*60)
        print("PHENOTYPE PHASE PLANE SUMMARY")
        print("="*60)
        print(f"Total points analyzed: {len(self.results)}")
        
        # Filter for significant growth
        significant_growth = self.results[self.results.growth_rate > 1e-10]
        
        if len(significant_growth) > 0:
            print(f"Points with significant growth (>1e-10): {len(significant_growth)}")
            print(f"Maximum growth rate: {significant_growth.growth_rate.max():.6f} 1/h")
            print(f"Minimum growth rate: {significant_growth.growth_rate.min():.6f} 1/h")
            print(f"Mean growth rate: {significant_growth.growth_rate.mean():.6f} 1/h")
        else:
            print("No significant growth detected in any condition")
        
        if len(lo_points) > 1:
            try:
                r_squared = np.corrcoef(lo_points.EX_ac_e, lo_points.EX_o2_e)[0,1]**2
                print(f"\nLine of Optimality:")
                print(f"  Slope: {m_ac_o2:.3f} mmol O₂ / mmol acetate")
                print(f"  Intercept: {b_ac_o2:.3f} mmol O₂ / gDW/h")
                print(f"  R²: {r_squared:.3f}")
            except:
                print(f"\nLine of Optimality:")
                print(f"  Slope: {m_ac_o2:.3f} mmol O₂ / mmol acetate")
                print(f"  Intercept: {b_ac_o2:.3f} mmol O₂ / gDW/h")
                print(f"  R²: Could not compute")
        else:
            print("\nLine of Optimality: Not enough points for analysis")
    
    def print_best_conditions(self) -> None:
        """
        Print the best growth conditions found.
        """
        # Filter for significant growth
        significant_growth = self.results[self.results.growth_rate > 1e-10]
        
        if len(significant_growth) == 0:
            print("\n" + "="*60)
            print("BEST GROWTH CONDITIONS")
            print("="*60)
            print("No significant growth detected in any condition")
            return
        
        best_conditions = significant_growth.loc[significant_growth.growth_rate.idxmax()]
        
        print("\n" + "="*60)
        print("BEST GROWTH CONDITIONS")
        print("="*60)
        print(f"Acetate uptake: {best_conditions['EX_ac_e']:.2f} mmol/gDW/h")
        print(f"Oxygen uptake: {best_conditions['EX_o2_e']:.2f} mmol/gDW/h")
        print(f"Growth rate: {best_conditions['growth_rate']:.6f} 1/h")
        
        # Show top 5 growth conditions
        print(f"\nTop 5 Growth Conditions:")
        print("-" * 60)
        top_5 = significant_growth.nlargest(5, 'growth_rate')
        print(top_5[['EX_ac_e', 'EX_o2_e', 'growth_rate']].to_string(index=False))
    
    def save_results(self, filepath: str) -> None:
        """
        Save the results to a CSV file.
        
        Parameters:
        -----------
        filepath : str
            Path to save the CSV file
        """
        if self.results is not None:
            self.results.to_csv(filepath, index=False)
            logger.info(f"Results saved to {filepath}")
        else:
            logger.warning("No results to save")


def main():
    """
    Main function to run the phenotype phase plane analysis.
    """
    # Configuration
    MODEL_PATH = "bigg_models/iJO1366.xml"
    AC_LIMITS = np.linspace(0, -20, 21)   # Acetate uptake: 0 to -20 mmol/gDW/h
    O2_LIMITS = np.linspace(0, -30, 31)   # O2 uptake: 0 to -30 mmol/gDW/h
    
    # Initialize analyzer
    analyzer = PhenotypePhasePlaneAnalyzer(MODEL_PATH)
    
    # Load model
    if not analyzer.load_model():
        logger.error("Failed to load model. Exiting.")
        return
    
    # Set up exchange reactions
    if not analyzer.setup_exchange_reactions():
        logger.error("Failed to set up exchange reactions. Exiting.")
        return
    
    # Fix carbon sources
    analyzer.fix_carbon_sources()
    
    # Test model feasibility
    if not analyzer.test_model_feasibility():
        logger.warning("Model may not be able to grow under current constraints")
    
    # Compute phenotype phase plane
    results = analyzer.compute_phenotype_phase_plane(AC_LIMITS, O2_LIMITS)
    
    # Find line of optimality
    m_ac_o2, b_ac_o2, lo_points = analyzer.find_line_of_optimality()
    
    # Create visualization
    analyzer.create_visualization(m_ac_o2, b_ac_o2, lo_points, save_path="phenotype_phase_plane.png")
    
    # Print summary statistics
    analyzer.print_summary_statistics(m_ac_o2, b_ac_o2, lo_points)
    
    # Print best conditions
    analyzer.print_best_conditions()
    
    # Save results
    analyzer.save_results("phpp_ac_o2.csv")
    
    logger.info("Analysis completed successfully!")


if __name__ == "__main__":
    main() 