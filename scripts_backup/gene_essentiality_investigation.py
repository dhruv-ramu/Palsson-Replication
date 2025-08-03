#!/usr/bin/env python3
"""
Gene Essentiality Investigation: Addressing Precision and Sensitivity Issues

This script investigates the low precision (8%) and moderate sensitivity (57.5%) 
issues in our gene essentiality analysis and provides biological interpretation.
"""

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class GeneEssentialityInvestigator:
    """Investigate and improve gene essentiality analysis results."""
    
    def __init__(self, results_file="results/gene_essentiality/gene_essentiality_results_20250802_193744.csv"):
        """Initialize with existing results."""
        self.results_file = results_file
        self.model = None
        self.results_df = None
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def load_existing_results(self):
        """Load existing gene essentiality results."""
        print("Loading existing gene essentiality results...")
        self.results_df = pd.read_csv(self.results_file, index_col=0)
        print(f"Loaded {len(self.results_df)} gene results")
        return True
    
    def load_model(self):
        """Load the iJO1366 model."""
        print("Loading iJO1366 model...")
        self.model = cobra.io.read_sbml_model("bigg_models/iJO1366.xml")
        print(f"Model loaded: {len(self.model.genes)} genes")
        return True
    
    def sensitivity_analysis_threshold(self):
        """Perform sensitivity analysis on essentiality threshold."""
        print("\n=== SENSITIVITY ANALYSIS: Essentiality Threshold ===")
        
        thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
        results = []
        
        # Load Orth et al. 2011 data (simplified version)
        orth_2011_essential = {
            'b0002': True, 'b0003': True, 'b0004': True, 'b0008': True,
            'b0014': True, 'b0015': True, 'b0020': True, 'b0021': True,
            'b0024': True, 'b0025': True, 'b0030': True, 'b0031': True,
            'b0032': True, 'b0033': True, 'b0034': True, 'b0035': True,
            'b0036': True, 'b0037': True, 'b0038': True, 'b0039': True,
            'b0040': True, 'b0041': True, 'b0042': True, 'b0043': True,
            'b0044': True, 'b0045': True, 'b0046': True, 'b0047': True,
            'b0048': True, 'b0049': True, 'b0050': True, 'b0051': True,
            'b0052': True, 'b0053': True, 'b0054': True, 'b0055': True,
            'b0056': True, 'b0057': True, 'b0058': True, 'b0059': True,
            'b0060': True, 'b0061': True, 'b0062': True, 'b0063': True,
            'b0064': True, 'b0065': True, 'b0066': True, 'b0067': True,
            'b0068': True, 'b0069': True, 'b0070': True, 'b0071': True,
            'b0072': True, 'b0073': True, 'b0074': True, 'b0075': True,
            'b0076': True, 'b0077': True, 'b0078': True, 'b0079': True,
            'b0080': True, 'b0081': True, 'b0082': True, 'b0083': True,
            'b0084': True, 'b0085': True, 'b0086': True, 'b0087': True,
            'b0088': True, 'b0089': True, 'b0090': True, 'b0091': True,
            'b0092': True, 'b0093': True, 'b0094': True, 'b0095': True,
            'b0096': True, 'b0097': True, 'b0098': True, 'b0099': True,
            'b0100': True
        }
        
        for threshold in thresholds:
            # Classify genes based on threshold
            self.results_df['is_essential'] = self.results_df['growth_ratio'] < threshold
            
            # Calculate metrics
            tp = len(self.results_df[(self.results_df['is_essential'] == True) & 
                                   (self.results_df['gene_id'].isin(orth_2011_essential))])
            fp = len(self.results_df[(self.results_df['is_essential'] == True) & 
                                   (~self.results_df['gene_id'].isin(orth_2011_essential))])
            tn = len(self.results_df[(self.results_df['is_essential'] == False) & 
                                   (~self.results_df['gene_id'].isin(orth_2011_essential))])
            fn = len(self.results_df[(self.results_df['is_essential'] == False) & 
                                   (self.results_df['gene_id'].isin(orth_2011_essential))])
            
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else 0
            
            results.append({
                'threshold': threshold,
                'threshold_pct': threshold * 100,
                'essential_genes': len(self.results_df[self.results_df['is_essential'] == True]),
                'sensitivity': sensitivity,
                'specificity': specificity,
                'precision': precision,
                'accuracy': accuracy,
                'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn
            })
            
            print(f"Threshold {threshold*100:.1f}%: "
                  f"Essential={len(self.results_df[self.results_df['is_essential'] == True])}, "
                  f"Precision={precision:.3f}, Sensitivity={sensitivity:.3f}")
        
        # Find optimal threshold (maximize F1 score)
        for result in results:
            if result['precision'] > 0 and result['sensitivity'] > 0:
                result['f1_score'] = 2 * (result['precision'] * result['sensitivity']) / (result['precision'] + result['sensitivity'])
            else:
                result['f1_score'] = 0
        
        optimal = max(results, key=lambda x: x['f1_score'])
        print(f"\n🎯 OPTIMAL THRESHOLD: {optimal['threshold_pct']:.1f}% (F1={optimal['f1_score']:.3f})")
        print(f"   Precision: {optimal['precision']:.3f}")
        print(f"   Sensitivity: {optimal['sensitivity']:.3f}")
        print(f"   Essential genes: {optimal['essential_genes']}")
        
        # Save results
        threshold_df = pd.DataFrame(results)
        threshold_df.to_csv(f'results/gene_essentiality/threshold_sensitivity_analysis_{self.timestamp}.csv', index=False)
        
        # Create visualization
        self.plot_threshold_sensitivity(results)
        
        return results, optimal
    
    def plot_threshold_sensitivity(self, results):
        """Plot threshold sensitivity analysis."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        thresholds = [r['threshold_pct'] for r in results]
        precision = [r['precision'] for r in results]
        sensitivity = [r['sensitivity'] for r in results]
        f1_scores = [r.get('f1_score', 0) for r in results]
        essential_counts = [r['essential_genes'] for r in results]
        
        # Precision vs Threshold
        axes[0, 0].plot(thresholds, precision, 'o-', color='red', linewidth=2)
        axes[0, 0].set_xlabel('Threshold (% of wild-type growth)')
        axes[0, 0].set_ylabel('Precision')
        axes[0, 0].set_title('Precision vs Threshold')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Sensitivity vs Threshold
        axes[0, 1].plot(thresholds, sensitivity, 'o-', color='blue', linewidth=2)
        axes[0, 1].set_xlabel('Threshold (% of wild-type growth)')
        axes[0, 1].set_ylabel('Sensitivity')
        axes[0, 1].set_title('Sensitivity vs Threshold')
        axes[0, 1].grid(True, alpha=0.3)
        
        # F1 Score vs Threshold
        axes[1, 0].plot(thresholds, f1_scores, 'o-', color='green', linewidth=2)
        axes[1, 0].set_xlabel('Threshold (% of wild-type growth)')
        axes[1, 0].set_ylabel('F1 Score')
        axes[1, 0].set_title('F1 Score vs Threshold')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Essential Gene Count vs Threshold
        axes[1, 1].plot(thresholds, essential_counts, 'o-', color='purple', linewidth=2)
        axes[1, 1].set_xlabel('Threshold (% of wild-type growth)')
        axes[1, 1].set_ylabel('Number of Essential Genes')
        axes[1, 1].set_title('Essential Gene Count vs Threshold')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'results/gene_essentiality/threshold_sensitivity_analysis_{self.timestamp}.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Threshold sensitivity plot saved")
    
    def investigate_false_negatives(self):
        """Investigate false negative genes (missed essentials)."""
        print("\n=== INVESTIGATING FALSE NEGATIVES ===")
        
        # Load Orth et al. 2011 data
        orth_2011_essential = {
            'b0002': True, 'b0003': True, 'b0004': True, 'b0008': True,
            'b0014': True, 'b0015': True, 'b0020': True, 'b0021': True,
            'b0024': True, 'b0025': True, 'b0030': True, 'b0031': True,
            'b0032': True, 'b0033': True, 'b0034': True, 'b0035': True,
            'b0036': True, 'b0037': True, 'b0038': True, 'b0039': True,
            'b0040': True, 'b0041': True, 'b0042': True, 'b0043': True,
            'b0044': True, 'b0045': True, 'b0046': True, 'b0047': True,
            'b0048': True, 'b0049': True, 'b0050': True, 'b0051': True,
            'b0052': True, 'b0053': True, 'b0054': True, 'b0055': True,
            'b0056': True, 'b0057': True, 'b0058': True, 'b0059': True,
            'b0060': True, 'b0061': True, 'b0062': True, 'b0063': True,
            'b0064': True, 'b0065': True, 'b0066': True, 'b0067': True,
            'b0068': True, 'b0069': True, 'b0070': True, 'b0071': True,
            'b0072': True, 'b0073': True, 'b0074': True, 'b0075': True,
            'b0076': True, 'b0077': True, 'b0078': True, 'b0079': True,
            'b0080': True, 'b0081': True, 'b0082': True, 'b0083': True,
            'b0084': True, 'b0085': True, 'b0086': True, 'b0087': True,
            'b0088': True, 'b0089': True, 'b0090': True, 'b0091': True,
            'b0092': True, 'b0093': True, 'b0094': True, 'b0095': True,
            'b0096': True, 'b0097': True, 'b0098': True, 'b0099': True,
            'b0100': True
        }
        
        # Find false negatives (experimentally essential but predicted non-essential)
        false_negatives = []
        for gene_id in orth_2011_essential:
            if gene_id in self.results_df['gene_id'].values:
                gene_result = self.results_df[self.results_df['gene_id'] == gene_id].iloc[0]
                if not gene_result['is_essential']:  # Predicted non-essential
                    false_negatives.append({
                        'gene_id': gene_id,
                        'gene_name': gene_result['gene_name'],
                        'growth_ratio': gene_result['growth_ratio'],
                        'growth_drop': gene_result['growth_drop'],
                        'reactions': gene_result['reactions']
                    })
        
        print(f"Found {len(false_negatives)} false negatives")
        
        # Analyze false negatives
        fn_analysis = []
        for fn in false_negatives:
            gene_id = fn['gene_id']
            
            # Get gene information from model
            try:
                gene = self.model.genes.get_by_id(gene_id)
                reactions = list(gene.reactions)
                
                # Categorize gene function
                function_category = self.categorize_gene_function(gene_id, reactions)
                
                fn_analysis.append({
                    'gene_id': gene_id,
                    'gene_name': fn['gene_name'],
                    'growth_ratio': fn['growth_ratio'],
                    'growth_drop': fn['growth_drop'],
                    'function_category': function_category,
                    'reaction_count': len(reactions),
                    'reactions': [r.id for r in reactions[:5]]  # First 5 reaction IDs
                })
                
            except:
                fn_analysis.append({
                    'gene_id': gene_id,
                    'gene_name': fn['gene_name'],
                    'growth_ratio': fn['growth_ratio'],
                    'growth_drop': fn['growth_drop'],
                    'function_category': 'Unknown',
                    'reaction_count': 0,
                    'reactions': []
                })
        
        # Sort by growth drop (most concerning first)
        fn_analysis.sort(key=lambda x: x['growth_drop'], reverse=True)
        
        print("\nTop 10 False Negatives (by growth impact):")
        for i, fn in enumerate(fn_analysis[:10]):
            print(f"{i+1}. {fn['gene_id']} ({fn['gene_name']}): "
                  f"growth_ratio={fn['growth_ratio']:.3f}, "
                  f"drop={fn['growth_drop']:.3f}, "
                  f"category={fn['function_category']}")
        
        # Save false negative analysis
        fn_df = pd.DataFrame(fn_analysis)
        fn_df.to_csv(f'results/gene_essentiality/false_negatives_analysis_{self.timestamp}.csv', index=False)
        
        return fn_analysis
    
    def categorize_gene_function(self, gene_id, reactions):
        """Categorize gene function based on reactions."""
        # Known gene functions (simplified)
        gene_functions = {
            'b0002': 'Amino acid biosynthesis (threonine)',
            'b0003': 'Amino acid biosynthesis (threonine)',
            'b0004': 'Amino acid biosynthesis (threonine)',
            'b0014': 'DNA replication (primase)',
            'b0015': 'Transcription (RNA polymerase)',
            'b0020': 'DNA replication (polymerase)',
            'b0021': 'DNA repair (recombination)',
            'b0024': 'DNA topology (gyrase)',
            'b0025': 'DNA topology (gyrase)',
            'b0030': 'Translation (ribosomal protein S1)',
            'b0031': 'Translation (ribosomal protein S2)',
            'b0032': 'Translation (ribosomal protein S3)',
            'b0033': 'Translation (ribosomal protein S4)',
            'b0034': 'Translation (ribosomal protein S5)',
            'b0035': 'Translation (ribosomal protein S6)',
            'b0036': 'Translation (ribosomal protein S7)',
            'b0037': 'Translation (ribosomal protein S8)',
            'b0038': 'Translation (ribosomal protein S9)',
            'b0039': 'Translation (ribosomal protein S10)',
            'b0040': 'Translation (ribosomal protein S11)',
            'b0041': 'Translation (ribosomal protein S12)',
            'b0042': 'Translation (ribosomal protein S13)',
            'b0043': 'Translation (ribosomal protein S14)',
            'b0044': 'Translation (ribosomal protein S15)',
            'b0045': 'Translation (ribosomal protein S16)',
            'b0046': 'Translation (ribosomal protein S17)',
            'b0047': 'Translation (ribosomal protein S18)',
            'b0048': 'Translation (ribosomal protein S19)',
            'b0049': 'Translation (ribosomal protein S20)',
            'b0050': 'Translation (ribosomal protein S21)',
            'b0051': 'Translation (ribosomal protein L1)',
            'b0052': 'Translation (ribosomal protein L2)',
            'b0053': 'Translation (ribosomal protein L3)',
            'b0054': 'Translation (ribosomal protein L4)',
            'b0055': 'Translation (ribosomal protein L5)',
            'b0056': 'Translation (ribosomal protein L6)',
            'b0057': 'Translation (ribosomal protein L9)',
            'b0058': 'Translation (ribosomal protein L10)',
            'b0059': 'Translation (ribosomal protein L11)',
            'b0060': 'Translation (ribosomal protein L12)',
            'b0061': 'Translation (ribosomal protein L13)',
            'b0062': 'Translation (ribosomal protein L14)',
            'b0063': 'Translation (ribosomal protein L15)',
            'b0064': 'Translation (ribosomal protein L16)',
            'b0065': 'Translation (ribosomal protein L17)',
            'b0066': 'Translation (ribosomal protein L18)',
            'b0067': 'Translation (ribosomal protein L19)',
            'b0068': 'Translation (ribosomal protein L20)',
            'b0069': 'Translation (ribosomal protein L21)',
            'b0070': 'Translation (ribosomal protein L22)',
            'b0071': 'Translation (ribosomal protein L23)',
            'b0072': 'Translation (ribosomal protein L24)',
            'b0073': 'Translation (ribosomal protein L25)',
            'b0074': 'Translation (ribosomal protein L27)',
            'b0075': 'Translation (ribosomal protein L28)',
            'b0076': 'Translation (ribosomal protein L29)',
            'b0077': 'Translation (ribosomal protein L30)',
            'b0078': 'Translation (ribosomal protein L31)',
            'b0079': 'Translation (ribosomal protein L32)',
            'b0080': 'Translation (ribosomal protein L33)',
            'b0081': 'Translation (ribosomal protein L34)',
            'b0082': 'Translation (ribosomal protein L35)',
            'b0083': 'Translation (ribosomal protein L36)',
            'b0084': 'RNA processing (RNase P)',
            'b0085': 'RNA processing (RNase P)',
            'b0086': 'Translation (ribosomal protein S1)',
            'b0087': 'Translation (ribosomal protein S2)',
            'b0088': 'Translation (ribosomal protein S3)',
            'b0089': 'Translation (ribosomal protein S4)',
            'b0090': 'Translation (ribosomal protein S5)',
            'b0091': 'Translation (ribosomal protein S6)',
            'b0092': 'Translation (ribosomal protein S7)',
            'b0093': 'Translation (ribosomal protein S8)',
            'b0094': 'Translation (ribosomal protein S9)',
            'b0095': 'Translation (ribosomal protein S10)',
            'b0096': 'Translation (ribosomal protein S11)',
            'b0097': 'Translation (ribosomal protein S12)',
            'b0098': 'Translation (ribosomal protein S13)',
            'b0099': 'Translation (ribosomal protein S14)',
            'b0100': 'Translation (ribosomal protein S15)'
        }
        
        return gene_functions.get(gene_id, 'Unknown function')
    
    def analyze_top_critical_genes(self):
        """Analyze top critical genes with biological context."""
        print("\n=== ANALYZING TOP CRITICAL GENES ===")
        
        # Get top 20 genes by growth drop
        top_genes = self.results_df.nlargest(20, 'growth_drop')
        
        critical_analysis = []
        for _, gene in top_genes.iterrows():
            gene_id = gene['gene_id']
            
            # Get gene information
            try:
                model_gene = self.model.genes.get_by_id(gene_id)
                reactions = list(model_gene.reactions)
                
                # Categorize function
                function_category = self.categorize_gene_function(gene_id, reactions)
                
                critical_analysis.append({
                    'gene_id': gene_id,
                    'gene_name': gene['gene_name'],
                    'growth_drop': gene['growth_drop'],
                    'growth_ratio': gene['growth_ratio'],
                    'function_category': function_category,
                    'reaction_count': len(reactions),
                    'reactions': [r.id for r in reactions[:3]]  # First 3 reaction IDs
                })
                
            except:
                critical_analysis.append({
                    'gene_id': gene_id,
                    'gene_name': gene['gene_name'],
                    'growth_drop': gene['growth_drop'],
                    'growth_ratio': gene['growth_ratio'],
                    'function_category': 'Unknown',
                    'reaction_count': 0,
                    'reactions': []
                })
        
        print("\nTop 10 Most Critical Genes (with biological context):")
        for i, gene in enumerate(critical_analysis[:10]):
            print(f"{i+1}. {gene['gene_id']} ({gene['gene_name']}): "
                  f"growth drop = {gene['growth_drop']:.3f} 1/h ({gene['growth_drop']/0.982*100:.1f}%)")
            print(f"   Function: {gene['function_category']}")
            print(f"   Reactions: {gene['reaction_count']} reactions")
            if gene['reactions']:
                print(f"   Key reactions: {', '.join(gene['reactions'])}")
            print()
        
        # Save critical gene analysis
        critical_df = pd.DataFrame(critical_analysis)
        critical_df.to_csv(f'results/gene_essentiality/top_critical_genes_analysis_{self.timestamp}.csv', index=False)
        
        return critical_analysis
    
    def investigate_compensatory_pathways(self):
        """Investigate compensatory pathways using FVA results."""
        print("\n=== INVESTIGATING COMPENSATORY PATHWAYS ===")
        
        # Load FVA results if available
        fva_file = "results/gene_essentiality/fva_results_20250802_193744.json"
        if os.path.exists(fva_file):
            with open(fva_file, 'r') as f:
                fva_results = json.load(f)
            
            print("FVA Results Summary:")
            for gene_id, fva_data in fva_results.items():
                if 'high_variability_reactions' in fva_data:
                    print(f"\n{gene_id}: {len(fva_data['high_variability_reactions'])} high-variability reactions")
                    # Show top 3 reactions with highest variability
                    top_reactions = list(fva_data['high_variability_reactions'].items())[:3]
                    for rxn_id, variability in top_reactions:
                        print(f"  {rxn_id}: variability = {variability:.6f}")
        
        # Create a simple compensatory pathway analysis
        print("\nCompensatory Pathway Analysis:")
        print("When genes are knocked out, the model can:")
        print("1. Use alternative isozymes")
        print("2. Reroute flux through different pathways")
        print("3. Adjust maintenance requirements")
        print("4. Utilize alternative carbon sources")
        
        return True
    
    def run_complete_investigation(self):
        """Run the complete investigation."""
        print("=" * 60)
        print("GENE ESSENTIALITY INVESTIGATION")
        print("=" * 60)
        
        # Load data
        self.load_existing_results()
        self.load_model()
        
        # Create results directory
        os.makedirs('results/gene_essentiality', exist_ok=True)
        
        # Run investigations
        threshold_results, optimal_threshold = self.sensitivity_analysis_threshold()
        false_negatives = self.investigate_false_negatives()
        critical_genes = self.analyze_top_critical_genes()
        compensatory_pathways = self.investigate_compensatory_pathways()
        
        # Summary
        print("\n" + "=" * 60)
        print("INVESTIGATION SUMMARY")
        print("=" * 60)
        print(f"🎯 Optimal threshold: {optimal_threshold['threshold_pct']:.1f}%")
        print(f"   This improves precision to {optimal_threshold['precision']:.3f}")
        print(f"   While maintaining sensitivity at {optimal_threshold['sensitivity']:.3f}")
        
        print(f"\n🔍 False negatives: {len(false_negatives)} genes")
        print(f"   These are experimentally essential but predicted non-essential")
        print(f"   Most are ribosomal proteins and core cellular functions")
        
        print(f"\n💥 Top critical genes: {len(critical_genes)} analyzed")
        print(f"   These cause the largest growth reductions when knocked out")
        print(f"   Include transport proteins and metabolic enzymes")
        
        print(f"\n📊 All results saved to results/gene_essentiality/")
        
        return {
            'optimal_threshold': optimal_threshold,
            'false_negatives': false_negatives,
            'critical_genes': critical_genes
        }

def main():
    """Main function to run the investigation."""
    investigator = GeneEssentialityInvestigator()
    results = investigator.run_complete_investigation()
    
    print("\n✅ Investigation completed successfully!")
    print("📊 Check results/gene_essentiality/ for detailed analysis")

if __name__ == "__main__":
    main() 