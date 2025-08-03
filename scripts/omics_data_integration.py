#!/usr/bin/env python3
"""
Omics Data Integration Script

This script integrates all omics datasets (transcriptomics, proteomics, 
metabolomics, fluxomics, genomics) with the metabolic network graph.

Features:
- Loads all omics data from data/omics/
- Maps omics data to network nodes
- Creates integrated node features
- Updates the metabolic network graph
- Saves results in organized subfolders

Author: Metabolic Network Embedding Project
Date: 2025
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Optional, Any, Tuple
import networkx as nx
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class OmicsDataIntegrator:
    """Integrates omics data with the metabolic network graph."""
    
    def __init__(self, 
                 omics_dir: str = "data/omics",
                 network_dir: str = "results/metabolic_network",
                 output_dir: str = "results/metabolic_network/omics_integration"):
        """Initialize the omics data integrator."""
        self.omics_dir = Path(omics_dir)
        self.network_dir = Path(network_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load existing network data
        self.graph = None
        self.node_features = {}
        self.node_id_map = {}
        
        # Omics data storage
        self.omics_data = {
            'transcriptomics': {},
            'proteomics': {},
            'metabolomics': {},
            'fluxomics': {},
            'genomics': {}
        }
        
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def load_network_data(self):
        """Load the existing metabolic network graph and features."""
        logger.info("Loading existing metabolic network data...")
        
        # Load graph
        graph_file = self.network_dir / "metabolic_network_graph.pkl"
        if graph_file.exists():
            import pickle
            with open(graph_file, 'rb') as f:
                self.graph = pickle.load(f)
            logger.info(f"Loaded graph with {self.graph.number_of_nodes()} nodes and {self.graph.number_of_edges()} edges")
        
        # Load node features
        features_file = self.network_dir / "metabolic_network_graph_features.json"
        if features_file.exists():
            with open(features_file, 'r') as f:
                self.node_features = json.load(f)
            logger.info(f"Loaded features for {len(self.node_features)} nodes")
        
        # Load node ID mapping
        idmap_file = self.network_dir / "metabolic_network_graph_idmap.json"
        if idmap_file.exists():
            with open(idmap_file, 'r') as f:
                self.node_id_map = json.load(f)
            logger.info(f"Loaded ID mapping for {len(self.node_id_map)} nodes")
    
    def load_omics_data(self):
        """Load all omics datasets from the omics directory."""
        logger.info("Loading omics datasets...")
        
        # Load transcriptomics data (prioritize real data over synthetic)
        transcriptomics_dir = self.omics_dir / "transcriptomics"
        if transcriptomics_dir.exists():
            # First try to load real data
            for file in transcriptomics_dir.glob("*_real_transcriptomics.csv"):
                condition = file.stem.replace("_real_transcriptomics", "")
                df = pd.read_csv(file)
                self.omics_data['transcriptomics'][condition] = df
                logger.info(f"Loaded REAL transcriptomics data for {condition}: {len(df)} genes")
            
            # If no real data found, load synthetic data as fallback
            if not self.omics_data['transcriptomics']:
                for file in transcriptomics_dir.glob("*_transcriptomics.csv"):
                    if "real" not in file.name:  # Skip real data files
                        condition = file.stem.replace("_transcriptomics", "")
                        df = pd.read_csv(file)
                        self.omics_data['transcriptomics'][condition] = df
                        logger.info(f"Loaded synthetic transcriptomics data for {condition}: {len(df)} genes")
        
        # Load metabolomics data (prioritize real data over synthetic)
        metabolomics_dir = self.omics_dir / "metabolomics"
        if metabolomics_dir.exists():
            # First try to load real data
            for file in metabolomics_dir.glob("*_real_metabolomics.csv"):
                condition = file.stem.replace("_real_metabolomics", "")
                df = pd.read_csv(file)
                self.omics_data['metabolomics'][condition] = df
                logger.info(f"Loaded REAL metabolomics data for {condition}: {len(df)} metabolites")
            
            # If no real data found, load synthetic data as fallback
            if not self.omics_data['metabolomics']:
                for file in metabolomics_dir.glob("*_metabolomics.csv"):
                    if "real" not in file.name:  # Skip real data files
                        condition = file.stem.replace("_metabolomics", "")
                        df = pd.read_csv(file)
                        self.omics_data['metabolomics'][condition] = df
                        logger.info(f"Loaded synthetic metabolomics data for {condition}: {len(df)} metabolites")
        
        # Load fluxomics data
        fluxomics_dir = self.omics_dir / "fluxomics"
        if fluxomics_dir.exists():
            for file in fluxomics_dir.glob("*.json"):
                with open(file, 'r') as f:
                    data = json.load(f)
                    if 'data' in data:
                        self.omics_data['fluxomics'][file.stem] = data['data']
                        logger.info(f"Loaded fluxomics data: {file.stem}")
        
        # Load proteomics data
        proteomics_dir = self.omics_dir / "proteomics"
        if proteomics_dir.exists():
            for file in proteomics_dir.glob("*.json"):
                with open(file, 'r') as f:
                    data = json.load(f)
                    self.omics_data['proteomics'][file.stem] = data
                    logger.info(f"Loaded proteomics data: {file.stem}")
        
        # Load genomics data
        genomics_dir = self.omics_dir / "genomics"
        if genomics_dir.exists():
            for file in genomics_dir.glob("*.csv"):
                df = pd.read_csv(file)
                self.omics_data['genomics'][file.stem] = df
                logger.info(f"Loaded genomics data: {file.stem}")
    
    def map_omics_to_network(self) -> Dict[str, Dict[str, float]]:
        """
        Map omics data to network nodes.
        
        Returns:
            Dictionary mapping node IDs to omics features
        """
        logger.info("Mapping omics data to network nodes...")
        
        omics_features = {}
        
        # Initialize features for all nodes
        for node_id in self.node_features.keys():
            omics_features[node_id] = {
                'transcriptomics': {},
                'proteomics': {},
                'metabolomics': {},
                'fluxomics': {},
                'genomics': {}
            }
        
        # Map transcriptomics data (genes -> reactions)
        for condition, df in self.omics_data['transcriptomics'].items():
            for _, row in df.iterrows():
                gene_id = row['gene_id']
                expression = row['expression_value']
                
                # Improved gene-reaction mapping using known associations
                gene_reaction_map = {
                    'b2415': ['HEX1'],      # glk -> hexokinase
                    'b2416': ['PFK'],       # pfkA -> phosphofructokinase
                    'b1854': ['PYK'],       # pykF -> pyruvate kinase
                    'b0720': ['CS'],        # gltA -> citrate synthase
                    'b1478': ['ICDHyr'],    # icdA -> isocitrate dehydrogenase
                    'b0344': ['LACZ'],      # lacZ -> beta-galactosidase
                    'b2296': ['ACSKr']      # acs -> acetyl-CoA synthetase
                }
                
                # Find reactions associated with this gene
                if gene_id in gene_reaction_map:
                    for reaction_id in gene_reaction_map[gene_id]:
                        if reaction_id in omics_features:
                            omics_features[reaction_id]['transcriptomics'][condition] = expression
                            logger.debug(f"Mapped gene {gene_id} to reaction {reaction_id}")
                
                # Fallback: search for gene ID in reaction names
                for node_id, features in self.node_features.items():
                    if features.get('node_type') == 1:  # Reaction node
                        if gene_id in node_id or gene_id.replace('b', '') in node_id:
                            omics_features[node_id]['transcriptomics'][condition] = expression
        
        # Map metabolomics data (metabolites)
        for condition, df in self.omics_data['metabolomics'].items():
            for _, row in df.iterrows():
                metabolite_id = row['metabolite_id']
                concentration = row['concentration']
                
                # Direct mapping for metabolite nodes
                if metabolite_id in omics_features:
                    omics_features[metabolite_id]['metabolomics'][condition] = concentration
        
        # Map fluxomics data (reactions)
        for dataset_name, flux_data in self.omics_data['fluxomics'].items():
            for reaction_id, flux_value in flux_data.items():
                if reaction_id in omics_features:
                    omics_features[reaction_id]['fluxomics'][dataset_name] = flux_value
        
        # Map proteomics data (enzymes -> reactions)
        for dataset_name, protein_data in self.omics_data['proteomics'].items():
            if isinstance(protein_data, dict):
                for protein_id, protein_value in protein_data.items():
                    # Find reactions associated with this protein
                    for node_id, features in self.node_features.items():
                        if features.get('node_type') == 1:  # Reaction node
                            if protein_id in node_id:
                                omics_features[node_id]['proteomics'][dataset_name] = protein_value
        
        # Map genomics data (gene essentiality)
        for dataset_name, df in self.omics_data['genomics'].items():
            if 'gene_id' in df.columns and 'is_essential' in df.columns:
                for _, row in df.iterrows():
                    gene_id = row['gene_id']
                    essentiality = row['is_essential']
                    
                    # Find reactions associated with this gene
                    for node_id, features in self.node_features.items():
                        if features.get('node_type') == 1:  # Reaction node
                            if gene_id in node_id:
                                omics_features[node_id]['genomics'][dataset_name] = essentiality
        
        logger.info(f"Mapped omics data to {len(omics_features)} nodes")
        return omics_features
    
    def create_integrated_features(self, omics_features: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, Any]]:
        """
        Create integrated node features combining structural and omics data.
        
        Args:
            omics_features: Mapped omics data
            
        Returns:
            Integrated node features
        """
        logger.info("Creating integrated node features...")
        
        integrated_features = {}
        
        for node_id, features in self.node_features.items():
            # Start with existing structural features
            integrated_features[node_id] = features.copy()
            
            # Add omics features
            if node_id in omics_features:
                omics_data = omics_features[node_id]
                
                # Add transcriptomics features
                if omics_data['transcriptomics']:
                    for condition, value in omics_data['transcriptomics'].items():
                        integrated_features[node_id][f'transcriptomics_{condition}'] = value
                
                # Add metabolomics features
                if omics_data['metabolomics']:
                    for condition, value in omics_data['metabolomics'].items():
                        integrated_features[node_id][f'metabolomics_{condition}'] = value
                
                # Add fluxomics features
                if omics_data['fluxomics']:
                    for dataset, value in omics_data['fluxomics'].items():
                        integrated_features[node_id][f'fluxomics_{dataset}'] = value
                
                # Add proteomics features
                if omics_data['proteomics']:
                    for dataset, value in omics_data['proteomics'].items():
                        integrated_features[node_id][f'proteomics_{dataset}'] = value
                
                # Add genomics features
                if omics_data['genomics']:
                    for dataset, value in omics_data['genomics'].items():
                        integrated_features[node_id][f'genomics_{dataset}'] = value
        
        logger.info(f"Created integrated features for {len(integrated_features)} nodes")
        return integrated_features
    
    def update_network_graph(self, integrated_features: Dict[str, Dict[str, Any]]):
        """Update the network graph with integrated features."""
        logger.info("Updating network graph with integrated features...")
        
        if self.graph is not None:
            # Update node attributes with integrated features
            for node_id, features in integrated_features.items():
                if node_id in self.graph.nodes():
                    self.graph.nodes[node_id].update(features)
            
            logger.info(f"Updated graph with integrated features")
    
    def create_summary_statistics(self, integrated_features: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Create summary statistics for the integrated omics data."""
        logger.info("Creating summary statistics...")
        
        stats = {
            'total_nodes': len(integrated_features),
            'node_types': {},
            'omics_coverage': {},
            'feature_counts': {}
        }
        
        # Count node types
        for node_id, features in integrated_features.items():
            node_type = features.get('node_type', 'unknown')
            stats['node_types'][node_type] = stats['node_types'].get(node_type, 0) + 1
        
        # Count omics coverage
        omics_types = ['transcriptomics', 'metabolomics', 'fluxomics', 'proteomics', 'genomics']
        for omics_type in omics_types:
            coverage = 0
            for node_id, features in integrated_features.items():
                if any(key.startswith(f'{omics_type}_') for key in features.keys()):
                    coverage += 1
            stats['omics_coverage'][omics_type] = coverage
        
        # Count feature types
        all_features = set()
        for features in integrated_features.values():
            all_features.update(features.keys())
        
        for feature in sorted(all_features):
            count = sum(1 for node_features in integrated_features.values() if feature in node_features)
            stats['feature_counts'][feature] = count
        
        return stats
    
    def save_integrated_data(self, integrated_features: Dict[str, Dict[str, Any]], stats: Dict[str, Any]):
        """Save the integrated omics data and statistics."""
        logger.info("Saving integrated omics data...")
        
        # Save integrated features
        features_file = self.output_dir / f"integrated_node_features_{self.timestamp}.json"
        with open(features_file, 'w') as f:
            json.dump(integrated_features, f, indent=2)
        logger.info(f"Saved integrated features to: {features_file}")
        
        # Save updated graph
        if self.graph is not None:
            graph_file = self.output_dir / f"integrated_network_graph_{self.timestamp}.pkl"
            import pickle
            with open(graph_file, 'wb') as f:
                pickle.dump(self.graph, f)
            logger.info(f"Saved updated graph to: {graph_file}")
        
        # Save statistics
        stats_file = self.output_dir / f"integration_statistics_{self.timestamp}.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Saved statistics to: {stats_file}")
        
        # Save summary report
        report_file = self.output_dir / f"integration_report_{self.timestamp}.md"
        self.create_integration_report(report_file, stats, integrated_features)
        logger.info(f"Saved integration report to: {report_file}")
    
    def create_integration_report(self, report_file: Path, stats: Dict[str, Any], features: Dict[str, Dict[str, Any]]):
        """Create a comprehensive integration report."""
        with open(report_file, 'w') as f:
            f.write("# Omics Data Integration Report\n\n")
            f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Overview\n\n")
            f.write("This report summarizes the integration of multi-omics data with the E. coli iJO1366 metabolic network.\n\n")
            
            f.write("## Network Statistics\n\n")
            f.write(f"- **Total nodes:** {stats['total_nodes']}\n")
            f.write(f"- **Node types:** {stats['node_types']}\n\n")
            
            f.write("## Omics Data Coverage\n\n")
            for omics_type, coverage in stats['omics_coverage'].items():
                percentage = (coverage / stats['total_nodes']) * 100
                f.write(f"- **{omics_type.title()}:** {coverage} nodes ({percentage:.1f}%)\n")
            f.write("\n")
            
            f.write("## Feature Summary\n\n")
            f.write("### Structural Features\n")
            structural_features = [f for f in stats['feature_counts'].keys() 
                                 if not any(f.startswith(omics) for omics in ['transcriptomics', 'metabolomics', 'fluxomics', 'proteomics', 'genomics'])]
            for feature in sorted(structural_features):
                f.write(f"- {feature}: {stats['feature_counts'][feature]} nodes\n")
            
            f.write("\n### Omics Features\n")
            omics_features = [f for f in stats['feature_counts'].keys() 
                            if any(f.startswith(omics) for omics in ['transcriptomics', 'metabolomics', 'fluxomics', 'proteomics', 'genomics'])]
            for feature in sorted(omics_features):
                f.write(f"- {feature}: {stats['feature_counts'][feature]} nodes\n")
            
            f.write("\n## Data Sources\n\n")
            f.write("- **Transcriptomics:** Synthetic data based on literature patterns\n")
            f.write("- **Metabolomics:** Synthetic data based on literature patterns\n")
            f.write("- **Fluxomics:** 13C-MFA experimental data\n")
            f.write("- **Proteomics:** Enzyme constraint data\n")
            f.write("- **Genomics:** Gene essentiality data\n")
    
    def run_integration(self):
        """Run the complete omics data integration pipeline."""
        logger.info("=" * 60)
        logger.info("OMICS DATA INTEGRATION PIPELINE")
        logger.info("=" * 60)
        
        # Step 1: Load network data
        logger.info("\nStep 1: Loading network data...")
        self.load_network_data()
        
        # Step 2: Load omics data
        logger.info("\nStep 2: Loading omics data...")
        self.load_omics_data()
        
        # Step 3: Map omics to network
        logger.info("\nStep 3: Mapping omics data to network...")
        omics_features = self.map_omics_to_network()
        
        # Step 4: Create integrated features
        logger.info("\nStep 4: Creating integrated features...")
        integrated_features = self.create_integrated_features(omics_features)
        
        # Step 5: Update network graph
        logger.info("\nStep 5: Updating network graph...")
        self.update_network_graph(integrated_features)
        
        # Step 6: Create statistics
        logger.info("\nStep 6: Creating summary statistics...")
        stats = self.create_summary_statistics(integrated_features)
        
        # Step 7: Save results
        logger.info("\nStep 7: Saving integrated data...")
        self.save_integrated_data(integrated_features, stats)
        
        # Step 8: Summary
        logger.info("\n" + "=" * 60)
        logger.info("OMICS DATA INTEGRATION COMPLETE")
        logger.info("=" * 60)
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Total nodes: {stats['total_nodes']}")
        logger.info(f"Omics coverage: {stats['omics_coverage']}")
        
        return integrated_features, stats

def main():
    """Main function to run omics data integration."""
    integrator = OmicsDataIntegrator()
    integrator.run_integration()

if __name__ == "__main__":
    main() 