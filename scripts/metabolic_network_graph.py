"""
Metabolic Network Graph Construction Module

This module converts SBML metabolic models to graph structures suitable for
Graph Neural Network processing. It creates nodes for metabolites and reactions,
edges for stoichiometric relationships, and initializes node features.

Author: Metabolic Network Embedding Project
Date: 2025
"""

import json
import networkx as nx
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MetabolicNetworkGraph:
    """
    Converts SBML metabolic models to graph structures for GNN processing.
    
    This class handles:
    - Loading SBML models (JSON format)
    - Creating NetworkX graphs with metabolites and reactions as nodes
    - Adding stoichiometric relationships as edges
    - Initializing node features from experimental data
    - Graph visualization and analysis
    """
    
    def __init__(self, model_path: str):
        """
        Initialize the metabolic network graph constructor.
        
        Args:
            model_path: Path to the SBML model file (JSON format)
        """
        self.model_path = Path(model_path)
        self.model_data = None
        self.graph = None
        self.metabolites = {}
        self.reactions = {}
        self.node_features = {}
        
        # Load the model
        self._load_model()
        
    def _load_model(self):
        """Load the SBML model from JSON file."""
        try:
            with open(self.model_path, 'r') as f:
                self.model_data = json.load(f)
            logger.info(f"Successfully loaded model: {self.model_path}")
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            raise
    
    def build_graph(self) -> nx.DiGraph:
        """
        Build the metabolic network as a directed graph.
        
        Returns:
            NetworkX directed graph representing the metabolic network
        """
        logger.info("Building metabolic network graph...")
        
        # Create directed graph
        self.graph = nx.DiGraph()
        
        # Add metabolite nodes
        self._add_metabolite_nodes()
        
        # Add reaction nodes
        self._add_reaction_nodes()
        
        # Add stoichiometric edges
        self._add_stoichiometric_edges()
        
        # Initialize node features
        self._initialize_node_features()
        
        logger.info(f"Graph built successfully: {self.graph.number_of_nodes()} nodes, "
                   f"{self.graph.number_of_edges()} edges")
        
        return self.graph
    
    def _add_metabolite_nodes(self):
        """Add metabolite nodes to the graph."""
        logger.info("Adding metabolite nodes...")
        
        for metabolite in self.model_data.get('metabolites', []):
            metabolite_id = metabolite['id']
            self.metabolites[metabolite_id] = metabolite
            
            # Add node with metabolite attributes
            self.graph.add_node(
                metabolite_id,
                node_type='metabolite',
                name=metabolite.get('name', ''),
                compartment=metabolite.get('compartment', ''),
                charge=metabolite.get('charge', 0),
                formula=metabolite.get('formula', ''),
                annotations=metabolite.get('annotation', {})
            )
        
        logger.info(f"Added {len(self.metabolites)} metabolite nodes")
    
    def _add_reaction_nodes(self):
        """Add reaction nodes to the graph."""
        logger.info("Adding reaction nodes...")
        
        for reaction in self.model_data.get('reactions', []):
            reaction_id = reaction['id']
            self.reactions[reaction_id] = reaction
            
            # Add node with reaction attributes
            self.graph.add_node(
                reaction_id,
                node_type='reaction',
                name=reaction.get('name', ''),
                subsystem=reaction.get('subsystem', ''),
                lower_bound=reaction.get('lower_bound', -1000),
                upper_bound=reaction.get('upper_bound', 1000),
                gene_reaction_rule=reaction.get('gene_reaction_rule', ''),
                metabolites=reaction.get('metabolites', {})
            )
        
        logger.info(f"Added {len(self.reactions)} reaction nodes")
    
    def _add_stoichiometric_edges(self):
        """Add stoichiometric relationships as edges."""
        logger.info("Adding stoichiometric edges...")
        
        edge_count = 0
        for reaction_id, reaction in self.reactions.items():
            metabolites = reaction.get('metabolites', {})
            
            for metabolite_id, stoichiometry in metabolites.items():
                if metabolite_id in self.metabolites:
                    # Add edge from metabolite to reaction (reactant)
                    if stoichiometry < 0:
                        self.graph.add_edge(
                            metabolite_id, 
                            reaction_id,
                            stoichiometry=abs(stoichiometry),
                            edge_type='reactant'
                        )
                        edge_count += 1
                    
                    # Add edge from reaction to metabolite (product)
                    elif stoichiometry > 0:
                        self.graph.add_edge(
                            reaction_id,
                            metabolite_id,
                            stoichiometry=stoichiometry,
                            edge_type='product'
                        )
                        edge_count += 1
        
        logger.info(f"Added {edge_count} stoichiometric edges")
    
    def _initialize_node_features(self):
        """Initialize node features for GNN processing."""
        logger.info("Initializing node features...")
        
        # Initialize features for all nodes
        for node in self.graph.nodes():
            node_data = self.graph.nodes[node]
            node_type = node_data['node_type']
            
            if node_type == 'metabolite':
                self._initialize_metabolite_features(node, node_data)
            elif node_type == 'reaction':
                self._initialize_reaction_features(node, node_data)
    
    def _initialize_metabolite_features(self, node_id: str, node_data: Dict):
        """Initialize features for metabolite nodes."""
        # Basic metabolite features
        features = {
            'node_type': 0,  # 0 for metabolite, 1 for reaction
            'charge': node_data.get('charge', 0),
            'compartment_encoded': self._encode_compartment(node_data.get('compartment', '')),
            'formula_length': len(node_data.get('formula', '')),
            'has_formula': 1 if node_data.get('formula') else 0,
            'annotation_count': len(node_data.get('annotations', {}))
        }
        
        # Store features
        self.node_features[node_id] = features
    
    def _initialize_reaction_features(self, node_id: str, node_data: Dict):
        """Initialize features for reaction nodes."""
        metabolites = node_data.get('metabolites', {})
        
        # Basic reaction features
        features = {
            'node_type': 1,  # 0 for metabolite, 1 for reaction
            'reactant_count': sum(1 for s in metabolites.values() if s < 0),
            'product_count': sum(1 for s in metabolites.values() if s > 0),
            'total_metabolites': len(metabolites),
            'lower_bound': node_data.get('lower_bound', -1000),
            'upper_bound': node_data.get('upper_bound', 1000),
            'reversible': 1 if node_data.get('lower_bound', 0) < 0 else 0,
            'subsystem_encoded': self._encode_subsystem(node_data.get('subsystem', '')),
            'has_gene_rule': 1 if node_data.get('gene_reaction_rule') else 0
        }
        
        # Store features
        self.node_features[node_id] = features
    
    def _encode_compartment(self, compartment: str) -> int:
        """Encode compartment as integer."""
        compartment_map = {
            'c': 0,  # cytosol
            'e': 1,  # extracellular
            'p': 2,  # periplasm
            'm': 3,  # mitochondria
            'n': 4,  # nucleus
            'g': 5,  # golgi
            'r': 6,  # endoplasmic reticulum
            'x': 7,  # peroxisome
            'v': 8,  # vacuole
            'l': 9   # lysosome
        }
        return compartment_map.get(compartment, 10)  # 10 for unknown
    
    def _encode_subsystem(self, subsystem: str) -> int:
        """Encode subsystem as integer."""
        # Common subsystems in E. coli
        subsystem_map = {
            'Transport, extracellular': 0,
            'Exchange/demand reaction': 1,
            'Glycolysis/gluconeogenesis': 2,
            'Citric acid cycle': 3,
            'Pentose phosphate pathway': 4,
            'Amino acid metabolism': 5,
            'Nucleotide metabolism': 6,
            'Lipid metabolism': 7,
            'Cofactor and vitamin metabolism': 8,
            'Cell envelope biosynthesis': 9,
            'Membrane lipid metabolism': 10,
            'Transport, inner membrane': 11,
            'Transport, outer membrane': 12,
            'Transport, outer membrane porin': 13,
            'Transport, periplasmic': 14,
            'Miscellaneous': 15
        }
        return subsystem_map.get(subsystem, 16)  # 16 for unknown
    
    def get_graph_statistics(self) -> Dict[str, Any]:
        """
        Calculate comprehensive graph statistics.
        
        Returns:
            Dictionary containing graph statistics
        """
        if self.graph is None:
            raise ValueError("Graph not built yet. Call build_graph() first.")
        
        stats = {
            'total_nodes': self.graph.number_of_nodes(),
            'total_edges': self.graph.number_of_edges(),
            'metabolite_nodes': len([n for n, d in self.graph.nodes(data=True) 
                                   if d['node_type'] == 'metabolite']),
            'reaction_nodes': len([n for n, d in self.graph.nodes(data=True) 
                                 if d['node_type'] == 'reaction']),
            'reactant_edges': len([e for e in self.graph.edges(data=True) 
                                 if e[2]['edge_type'] == 'reactant']),
            'product_edges': len([e for e in self.graph.edges(data=True) 
                                if e[2]['edge_type'] == 'product']),
            'density': nx.density(self.graph),
            'is_directed': self.graph.is_directed(),
            'is_connected': nx.is_weakly_connected(self.graph),
            'number_connected_components': nx.number_weakly_connected_components(self.graph),
            'average_clustering': nx.average_clustering(self.graph.to_undirected()),
            'average_shortest_path_length': self._calculate_avg_shortest_path(),
            'degree_distribution': dict(nx.degree(self.graph)),
            'node_features_dim': len(next(iter(self.node_features.values()))) if self.node_features else 0
        }
        
        return stats
    
    def _calculate_avg_shortest_path(self) -> float:
        """Calculate average shortest path length for connected components."""
        try:
            # Get largest connected component
            largest_cc = max(nx.weakly_connected_components(self.graph), key=len)
            subgraph = self.graph.subgraph(largest_cc)
            return nx.average_shortest_path_length(subgraph)
        except:
            return float('inf')
    
    def visualize_network(self, output_path: str = None, max_nodes: int = 100):
        """
        Create a visualization of the metabolic network.
        
        Args:
            output_path: Path to save the visualization
            max_nodes: Maximum number of nodes to display (for performance)
        """
        if self.graph is None:
            raise ValueError("Graph not built yet. Call build_graph() first.")
        
        logger.info("Creating network visualization...")
        
        # Create subgraph with limited nodes for visualization
        if self.graph.number_of_nodes() > max_nodes:
            # Sample nodes from largest connected component
            largest_cc = max(nx.weakly_connected_components(self.graph), key=len)
            sampled_nodes = list(largest_cc)[:max_nodes]
            viz_graph = self.graph.subgraph(sampled_nodes)
        else:
            viz_graph = self.graph
        
        # Create figure
        plt.figure(figsize=(15, 10))
        
        # Layout
        pos = nx.spring_layout(viz_graph, k=1, iterations=50)
        
        # Draw nodes
        metabolite_nodes = [n for n, d in viz_graph.nodes(data=True) 
                          if d['node_type'] == 'metabolite']
        reaction_nodes = [n for n, d in viz_graph.nodes(data=True) 
                         if d['node_type'] == 'reaction']
        
        # Draw metabolites (blue circles)
        nx.draw_networkx_nodes(viz_graph, pos, nodelist=metabolite_nodes,
                             node_color='lightblue', node_size=300, alpha=0.7)
        
        # Draw reactions (red squares)
        nx.draw_networkx_nodes(viz_graph, pos, nodelist=reaction_nodes,
                             node_color='lightcoral', node_size=200, alpha=0.7)
        
        # Draw edges
        nx.draw_networkx_edges(viz_graph, pos, alpha=0.3, arrows=True, arrowsize=10)
        
        # Add labels for a subset of nodes
        labels = {}
        for node in list(viz_graph.nodes())[:20]:  # Show first 20 labels
            node_data = viz_graph.nodes[node]
            if node_data['node_type'] == 'metabolite':
                labels[node] = node_data['name'][:10] if node_data['name'] else node[:10]
            else:
                labels[node] = node_data['name'][:10] if node_data['name'] else node[:10]
        
        nx.draw_networkx_labels(viz_graph, pos, labels, font_size=8)
        
        # Add legend
        plt.scatter([], [], c='lightblue', s=300, alpha=0.7, label='Metabolites')
        plt.scatter([], [], c='lightcoral', s=200, alpha=0.7, label='Reactions')
        plt.legend()
        
        plt.title(f'Metabolic Network Graph\n'
                 f'{viz_graph.number_of_nodes()} nodes, {viz_graph.number_of_edges()} edges')
        plt.axis('off')
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Visualization saved to: {output_path}")
        
        plt.show()
    
    def save_graph(self, output_path: str):
        """
        Save the graph and features to files.
        
        Args:
            output_path: Base path for saving graph files
        """
        if self.graph is None:
            raise ValueError("Graph not built yet. Call build_graph() first.")
        
        output_path = Path(output_path)
        
        # Create a clean graph for export (remove problematic annotations and node names)
        clean_graph = nx.DiGraph()
        node_id_map = {}
        node_counter = 0
        # Add nodes with clean IDs and minimal attributes
        for node in self.graph.nodes():
            node_data = self.graph.nodes[node]
            clean_node_data = {
                'node_type': node_data.get('node_type', 'unknown'),
                'name': node_data.get('name', node),
                'compartment': node_data.get('compartment', 'unknown')
            }
            # Use a sequential node ID for GML compatibility
            clean_node_id = f"n{node_counter}"
            node_id_map[clean_node_id] = node
            clean_graph.add_node(clean_node_id, **clean_node_data)
            node_counter += 1
        # Add edges
        for source, target, edge_data in self.graph.edges(data=True):
            clean_source = [k for k, v in node_id_map.items() if v == source][0]
            clean_target = [k for k, v in node_id_map.items() if v == target][0]
            clean_graph.add_edge(clean_source, clean_target, **edge_data)
        # Save graph as GML in a structured results folder
        output_dir = Path('results/metabolic_network')
        output_dir.mkdir(parents=True, exist_ok=True)
        graph_path = output_dir / (output_path.stem + '.gml')
        nx.write_gml(clean_graph, graph_path)
        logger.info(f"Graph saved to: {graph_path}")
        # Save node ID mapping for interpretability
        idmap_path = output_dir / (output_path.stem + '_idmap.json')
        with open(idmap_path, 'w') as f:
            json.dump(node_id_map, f, indent=2)
        logger.info(f"Node ID map saved to: {idmap_path}")
        # Save original graph as pickle for full data preservation
        import pickle
        pickle_path = output_dir / (output_path.stem + '.pkl')
        with open(pickle_path, 'wb') as f:
            pickle.dump(self.graph, f)
        logger.info(f"Full graph saved to: {pickle_path}")
        # Save node features as JSON
        features_path = output_dir / (output_path.stem + '_features.json')
        with open(features_path, 'w') as f:
            json.dump(self.node_features, f, indent=2)
        logger.info(f"Node features saved to: {features_path}")
        # Save statistics
        stats_path = output_dir / (output_path.stem + '_stats.json')
        with open(stats_path, 'w') as f:
            json.dump(self.get_graph_statistics(), f, indent=2)
        logger.info(f"Graph statistics saved to: {stats_path}")
    
    def get_node_feature_matrix(self) -> Tuple[np.ndarray, List[str]]:
        """
        Get node feature matrix for GNN processing.
        
        Returns:
            Tuple of (feature_matrix, node_ids)
        """
        if not self.node_features:
            raise ValueError("Node features not initialized. Call build_graph() first.")
        
        # Get feature names from first node
        feature_names = list(next(iter(self.node_features.values())).keys())
        
        # Create feature matrix
        node_ids = list(self.node_features.keys())
        feature_matrix = np.zeros((len(node_ids), len(feature_names)))
        
        for i, node_id in enumerate(node_ids):
            features = self.node_features[node_id]
            for j, feature_name in enumerate(feature_names):
                feature_matrix[i, j] = features.get(feature_name, 0)
        
        return feature_matrix, node_ids


def main():
    """Main function for testing the metabolic network graph construction."""
    
    # Path to the iJO1366 model
    model_path = "bigg_models/iJO1366.json"
    
    # Create metabolic network graph
    print("Creating metabolic network graph...")
    metabolic_graph = MetabolicNetworkGraph(model_path)
    
    # Build the graph
    graph = metabolic_graph.build_graph()
    
    # Get statistics
    stats = metabolic_graph.get_graph_statistics()
    print("\nGraph Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    
    # Create visualization
    print("\nCreating visualization...")
    metabolic_graph.visualize_network("results/metabolic_network_visualization.png")
    
    # Save graph
    print("\nSaving graph...")
    metabolic_graph.save_graph("results/metabolic_network_graph")
    
    # Get feature matrix
    feature_matrix, node_ids = metabolic_graph.get_node_feature_matrix()
    print(f"\nFeature matrix shape: {feature_matrix.shape}")
    print(f"Number of nodes: {len(node_ids)}")
    
    print("\nMetabolic network graph construction completed successfully!")


if __name__ == "__main__":
    main() 