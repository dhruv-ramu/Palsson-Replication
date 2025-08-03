#!/usr/bin/env python3
"""
Graph Feature Engineering Module

This module implements comprehensive feature engineering for the metabolic network
as required by Phase 1 Week 3 of the scope.

Features:
- Create node features from multi-omics data
- Engineer edge features for relationships
- Normalize and scale features
- Handle missing data
- Feature importance analysis

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
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.feature_selection import mutual_info_regression, f_regression
import matplotlib.pyplot as plt
import seaborn as sns

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class GraphFeatureEngineer:
    """Comprehensive feature engineering for metabolic network graph."""
    
    def __init__(self, 
                 network_dir: str = "results/metabolic_network",
                 omics_dir: str = "data/omics",
                 output_dir: str = "results/metabolic_network/feature_engineering"):
        """Initialize the feature engineer."""
        self.network_dir = Path(network_dir)
        self.omics_dir = Path(omics_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load existing data
        self.graph = None
        self.node_features = {}
        self.edge_features = {}
        
        # Feature engineering results
        self.processed_node_features = {}
        self.processed_edge_features = {}
        self.feature_importance = {}
        self.scalers = {}
        self.imputers = {}
        
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def load_network_data(self):
        """Load the metabolic network graph and existing features."""
        logger.info("Loading network data for feature engineering...")
        
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
    
    def create_node_features(self) -> Dict[str, np.ndarray]:
        """
        Create comprehensive node features from multi-omics data.
        
        Returns:
            Dictionary of node features
        """
        logger.info("Creating comprehensive node features...")
        
        node_features = {}
        
        for node_id, features in self.node_features.items():
            # Start with structural features
            feature_vector = []
            feature_names = []
            
            # Structural features (already normalized)
            structural_features = [
                'node_type', 'charge', 'compartment_encoded', 'formula_length',
                'has_formula', 'annotation_count', 'has_gene_rule', 'lower_bound',
                'upper_bound', 'reversible', 'subsystem_encoded', 'product_count',
                'reactant_count', 'total_metabolites'
            ]
            
            for feature in structural_features:
                if feature in features:
                    feature_vector.append(float(features[feature]))
                    feature_names.append(feature)
                else:
                    feature_vector.append(0.0)  # Default value
                    feature_names.append(feature)
            
            # Add omics features
            omics_features = [
                'transcriptomics_glucose_minimal', 'transcriptomics_acetate_minimal', 'transcriptomics_lactose_minimal',
                'metabolomics_glucose_minimal', 'metabolomics_acetate_minimal', 'metabolomics_lactose_minimal',
                'fluxomics_glucose_minimal', 'fluxomics_acetate_minimal', 'fluxomics_lactose_minimal',
                'proteomics_glucose_minimal', 'proteomics_acetate_minimal', 'proteomics_lactose_minimal',
                'genomics_essentiality', 'genomics_knockout_effect'
            ]
            
            for feature in omics_features:
                if feature in features:
                    value = features[feature]
                    if isinstance(value, (int, float)) and not np.isnan(value):
                        feature_vector.append(float(value))
                    else:
                        feature_vector.append(0.0)  # Missing value
                    feature_names.append(feature)
                else:
                    feature_vector.append(0.0)  # Missing feature
                    feature_names.append(feature)
            
            # Add derived features
            derived_features = self._create_derived_features(features)
            for feature_name, value in derived_features.items():
                feature_vector.append(value)
                feature_names.append(feature_name)
            
            node_features[node_id] = {
                'features': np.array(feature_vector, dtype=np.float32),
                'feature_names': feature_names
            }
        
        self.processed_node_features = node_features
        logger.info(f"Created features for {len(node_features)} nodes with {len(feature_names)} features each")
        
        return node_features
    
    def _create_derived_features(self, features: Dict[str, Any]) -> Dict[str, float]:
        """Create derived features from existing features."""
        derived = {}
        
        # Metabolic activity indicators
        if 'lower_bound' in features and 'upper_bound' in features:
            lb, ub = features['lower_bound'], features['upper_bound']
            derived['flux_range'] = abs(ub - lb)
            derived['flux_center'] = (lb + ub) / 2
            derived['flux_variability'] = abs(ub - lb) / (abs(ub) + abs(lb) + 1e-8)
        
        # Connectivity features
        if 'product_count' in features and 'reactant_count' in features:
            pc, rc = features['product_count'], features['reactant_count']
            derived['total_connectivity'] = pc + rc
            derived['connectivity_ratio'] = pc / (rc + 1e-8)
            derived['connectivity_balance'] = abs(pc - rc)
        
        # Regulatory complexity
        if 'has_gene_rule' in features and 'annotation_count' in features:
            derived['regulatory_complexity'] = features['has_gene_rule'] * features['annotation_count']
        
        # Metabolic pathway indicators
        if 'subsystem_encoded' in features:
            derived['pathway_specificity'] = features['subsystem_encoded'] / 100.0  # Normalized
        
        return derived
    
    def create_edge_features(self) -> Dict[Tuple[str, str], np.ndarray]:
        """
        Engineer edge features for relationships between nodes.
        
        Returns:
            Dictionary of edge features
        """
        logger.info("Creating edge features...")
        
        edge_features = {}
        
        for edge in self.graph.edges(data=True):
            source, target, edge_data = edge
            
            feature_vector = []
            feature_names = []
            
            # Basic edge features
            feature_vector.extend([
                1.0,  # Edge exists
                edge_data.get('weight', 1.0),  # Edge weight
                edge_data.get('stoichiometry', 1.0),  # Stoichiometric coefficient
            ])
            feature_names.extend(['edge_exists', 'edge_weight', 'stoichiometry'])
            
            # Node type relationship
            source_type = self.node_features.get(source, {}).get('node_type', 0)
            target_type = self.node_features.get(target, {}).get('node_type', 0)
            
            feature_vector.extend([
                source_type,
                target_type,
                1.0 if source_type != target_type else 0.0,  # Cross-type edge
            ])
            feature_names.extend(['source_type', 'target_type', 'cross_type_edge'])
            
            # Metabolic relationship features
            if source_type == 0 and target_type == 1:  # Metabolite -> Reaction
                feature_vector.extend([
                    1.0,  # Metabolite input
                    0.0,  # Metabolite output
                    0.0,  # Reaction input
                    0.0,  # Reaction output
                ])
            elif source_type == 1 and target_type == 0:  # Reaction -> Metabolite
                feature_vector.extend([
                    0.0,  # Metabolite input
                    1.0,  # Metabolite output
                    0.0,  # Reaction input
                    0.0,  # Reaction output
                ])
            else:
                feature_vector.extend([0.0, 0.0, 0.0, 0.0])
            
            feature_names.extend(['metabolite_input', 'metabolite_output', 'reaction_input', 'reaction_output'])
            
            # Pathway relationship
            source_subsystem = self.node_features.get(source, {}).get('subsystem_encoded', 0)
            target_subsystem = self.node_features.get(target, {}).get('subsystem_encoded', 0)
            
            feature_vector.extend([
                source_subsystem,
                target_subsystem,
                1.0 if source_subsystem == target_subsystem else 0.0,  # Same pathway
            ])
            feature_names.extend(['source_subsystem', 'target_subsystem', 'same_pathway'])
            
            # Regulatory relationship
            source_has_gene = self.node_features.get(source, {}).get('has_gene_rule', 0)
            target_has_gene = self.node_features.get(target, {}).get('has_gene_rule', 0)
            
            feature_vector.extend([
                source_has_gene,
                target_has_gene,
                1.0 if source_has_gene and target_has_gene else 0.0,  # Both regulated
            ])
            feature_names.extend(['source_has_gene', 'target_has_gene', 'both_regulated'])
            
            edge_features[(source, target)] = {
                'features': np.array(feature_vector, dtype=np.float32),
                'feature_names': feature_names
            }
        
        self.processed_edge_features = edge_features
        logger.info(f"Created features for {len(edge_features)} edges with {len(feature_names)} features each")
        
        return edge_features
    
    def handle_missing_data(self, features: Dict[str, np.ndarray], method: str = 'knn') -> Dict[str, np.ndarray]:
        """
        Handle missing data in features.
        
        Args:
            features: Dictionary of node/edge features
            method: Imputation method ('mean', 'median', 'knn', 'zero')
            
        Returns:
            Features with missing data imputed
        """
        logger.info(f"Handling missing data using {method} imputation...")
        
        # Convert to matrix
        feature_matrix = []
        node_ids = []
        original_features = {}
        
        for node_id, feature_data in features.items():
            feature_matrix.append(feature_data['features'])
            node_ids.append(node_id)
            original_features[node_id] = feature_data
        
        # Ensure all feature vectors have the same length
        feature_lengths = [len(f) for f in feature_matrix]
        if len(set(feature_lengths)) > 1:
            max_length = max(feature_lengths)
            # Pad shorter vectors with zeros
            for i, features in enumerate(feature_matrix):
                if len(features) < max_length:
                    feature_matrix[i] = np.pad(features, (0, max_length - len(features)), 'constant')
        
        feature_matrix = np.array(feature_matrix)
        
        # Check for missing values
        missing_mask = np.isnan(feature_matrix) | np.isinf(feature_matrix)
        missing_count = np.sum(missing_mask)
        
        if missing_count > 0:
            logger.info(f"Found {missing_count} missing values in feature matrix")
            
            if method == 'mean':
                imputer = SimpleImputer(strategy='mean')
            elif method == 'median':
                imputer = SimpleImputer(strategy='median')
            elif method == 'knn':
                imputer = KNNImputer(n_neighbors=5)
            elif method == 'zero':
                imputer = SimpleImputer(strategy='constant', fill_value=0.0)
            else:
                raise ValueError(f"Unknown imputation method: {method}")
            
            # Impute missing values
            feature_matrix_imputed = imputer.fit_transform(feature_matrix)
            self.imputers[method] = imputer
            
            logger.info(f"Imputed {missing_count} missing values using {method}")
        else:
            feature_matrix_imputed = feature_matrix
            logger.info("No missing values found")
        
        # Convert back to dictionary
        imputed_features = {}
        for i, node_id in enumerate(node_ids):
            imputed_features[node_id] = {
                'features': feature_matrix_imputed[i],
                'feature_names': original_features[node_id]['feature_names']
            }
        
        return imputed_features
    
    def normalize_and_scale_features(self, features: Dict[str, np.ndarray], 
                                   method: str = 'standard') -> Dict[str, np.ndarray]:
        """
        Normalize and scale features.
        
        Args:
            features: Dictionary of features
            method: Scaling method ('standard', 'minmax', 'robust')
            
        Returns:
            Scaled features
        """
        logger.info(f"Normalizing and scaling features using {method} scaling...")
        
        # Convert to matrix
        feature_matrix = []
        node_ids = []
        
        for node_id, feature_data in features.items():
            feature_matrix.append(feature_data['features'])
            node_ids.append(node_id)
        
        # Ensure all feature vectors have the same length
        feature_lengths = [len(f) for f in feature_matrix]
        if len(set(feature_lengths)) > 1:
            max_length = max(feature_lengths)
            # Pad shorter vectors with zeros
            for i, features in enumerate(feature_matrix):
                if len(features) < max_length:
                    feature_matrix[i] = np.pad(features, (0, max_length - len(features)), 'constant')
        
        feature_matrix = np.array(feature_matrix)
        
        # Choose scaler
        if method == 'standard':
            scaler = StandardScaler()
        elif method == 'minmax':
            scaler = MinMaxScaler()
        elif method == 'robust':
            scaler = RobustScaler()
        else:
            raise ValueError(f"Unknown scaling method: {method}")
        
        # Scale features
        feature_matrix_scaled = scaler.fit_transform(feature_matrix)
        self.scalers[method] = scaler
        
        logger.info(f"Scaled features using {method} scaling")
        
        # Convert back to dictionary
        scaled_features = {}
        for i, node_id in enumerate(node_ids):
            scaled_features[node_id] = {
                'features': feature_matrix_scaled[i],
                'feature_names': features[node_id]['feature_names']
            }
        
        return scaled_features
    
    def analyze_feature_importance(self, features: Dict[str, np.ndarray], 
                                 target_feature: str = 'node_type') -> Dict[str, float]:
        """
        Analyze feature importance using mutual information and correlation.
        
        Args:
            features: Dictionary of features
            target_feature: Target feature for importance analysis
            
        Returns:
            Dictionary of feature importance scores
        """
        logger.info(f"Analyzing feature importance for target: {target_feature}")
        
        # Convert to matrix
        feature_matrix = []
        target_values = []
        feature_names = None
        
        for node_id, feature_data in features.items():
            feature_matrix.append(feature_data['features'])
            if feature_names is None:
                feature_names = feature_data['feature_names']
            
            # Get target value
            target_idx = feature_names.index(target_feature) if target_feature in feature_names else 0
            target_values.append(feature_data['features'][target_idx])
        
        feature_matrix = np.array(feature_matrix)
        target_values = np.array(target_values)
        
        # Calculate mutual information
        mi_scores = mutual_info_regression(feature_matrix, target_values, random_state=42)
        
        # Calculate correlation
        correlations = []
        for i in range(feature_matrix.shape[1]):
            corr = np.corrcoef(feature_matrix[:, i], target_values)[0, 1]
            correlations.append(abs(corr) if not np.isnan(corr) else 0.0)
        
        # Combine scores
        importance_scores = {}
        for i, feature_name in enumerate(feature_names):
            importance_scores[feature_name] = {
                'mutual_info': float(mi_scores[i]),
                'correlation': float(correlations[i]),
                'combined_score': float(mi_scores[i] + correlations[i])
            }
        
        # Sort by combined score
        sorted_features = sorted(importance_scores.items(), 
                               key=lambda x: x[1]['combined_score'], reverse=True)
        
        self.feature_importance = dict(sorted_features)
        
        logger.info(f"Analyzed importance for {len(feature_names)} features")
        return self.feature_importance
    
    def create_feature_visualizations(self):
        """Create visualizations for feature analysis."""
        logger.info("Creating feature analysis visualizations...")
        
        # Feature importance plot
        if self.feature_importance:
            top_features = list(self.feature_importance.keys())[:20]
            mi_scores = [self.feature_importance[f]['mutual_info'] for f in top_features]
            corr_scores = [self.feature_importance[f]['correlation'] for f in top_features]
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Mutual information
            ax1.barh(range(len(top_features)), mi_scores)
            ax1.set_yticks(range(len(top_features)))
            ax1.set_yticklabels(top_features)
            ax1.set_xlabel('Mutual Information Score')
            ax1.set_title('Top 20 Features by Mutual Information')
            ax1.invert_yaxis()
            
            # Correlation
            ax2.barh(range(len(top_features)), corr_scores)
            ax2.set_yticks(range(len(top_features)))
            ax2.set_yticklabels(top_features)
            ax2.set_xlabel('Absolute Correlation')
            ax2.set_title('Top 20 Features by Correlation')
            ax2.invert_yaxis()
            
            plt.tight_layout()
            plt.savefig(self.output_dir / f"feature_importance_{self.timestamp}.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        # Feature distribution plots
        if self.processed_node_features:
            feature_matrix = []
            feature_names = None
            
            for node_id, feature_data in self.processed_node_features.items():
                feature_matrix.append(feature_data['features'])
                if feature_names is None:
                    feature_names = feature_data['feature_names']
            
            # Ensure all feature vectors have the same length
            feature_lengths = [len(f) for f in feature_matrix]
            if len(set(feature_lengths)) > 1:
                max_length = max(feature_lengths)
                # Pad shorter vectors with zeros
                for i, features in enumerate(feature_matrix):
                    if len(features) < max_length:
                        feature_matrix[i] = np.pad(features, (0, max_length - len(features)), 'constant')
            
            feature_matrix = np.array(feature_matrix)
            
            # Plot distributions for top features
            top_features = list(self.feature_importance.keys())[:8] if self.feature_importance else feature_names[:8]
            
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            axes = axes.flatten()
            
            for i, feature_name in enumerate(top_features):
                if feature_name in feature_names:
                    feature_idx = feature_names.index(feature_name)
                    axes[i].hist(feature_matrix[:, feature_idx], bins=30, alpha=0.7)
                    axes[i].set_title(feature_name)
                    axes[i].set_xlabel('Value')
                    axes[i].set_ylabel('Frequency')
            
            plt.tight_layout()
            plt.savefig(self.output_dir / f"feature_distributions_{self.timestamp}.png", dpi=300, bbox_inches='tight')
            plt.close()
    
    def save_processed_features(self):
        """Save processed features to files."""
        logger.info("Saving processed features...")
        
        # Save node features
        if self.processed_node_features:
            node_features_data = {}
            for node_id, feature_data in self.processed_node_features.items():
                node_features_data[node_id] = {
                    'features': feature_data['features'].tolist(),
                    'feature_names': feature_data['feature_names']
                }
            
            features_file = self.output_dir / f"processed_node_features_{self.timestamp}.json"
            with open(features_file, 'w') as f:
                json.dump(node_features_data, f, indent=2)
            logger.info(f"Saved node features to: {features_file}")
        
        # Save edge features
        if self.processed_edge_features:
            edge_features_data = {}
            for edge, feature_data in self.processed_edge_features.items():
                edge_key = f"{edge[0]}_{edge[1]}"
                edge_features_data[edge_key] = {
                    'features': feature_data['features'].tolist(),
                    'feature_names': feature_data['feature_names']
                }
            
            edge_file = self.output_dir / f"processed_edge_features_{self.timestamp}.json"
            with open(edge_file, 'w') as f:
                json.dump(edge_features_data, f, indent=2)
            logger.info(f"Saved edge features to: {edge_file}")
        
        # Save feature importance
        if self.feature_importance:
            importance_file = self.output_dir / f"feature_importance_{self.timestamp}.json"
            with open(importance_file, 'w') as f:
                json.dump(self.feature_importance, f, indent=2)
            logger.info(f"Saved feature importance to: {importance_file}")
    
    def run_feature_engineering(self):
        """Run the complete feature engineering pipeline."""
        logger.info("=" * 60)
        logger.info("GRAPH FEATURE ENGINEERING PIPELINE")
        logger.info("=" * 60)
        
        # Step 1: Load network data
        logger.info("\nStep 1: Loading network data...")
        self.load_network_data()
        
        # Step 2: Create node features
        logger.info("\nStep 2: Creating node features...")
        node_features = self.create_node_features()
        
        # Step 3: Create edge features
        logger.info("\nStep 3: Creating edge features...")
        edge_features = self.create_edge_features()
        
        # Step 4: Handle missing data
        logger.info("\nStep 4: Handling missing data...")
        node_features_imputed = self.handle_missing_data(node_features, method='knn')
        edge_features_imputed = self.handle_missing_data(edge_features, method='mean')
        
        # Step 5: Normalize and scale features
        logger.info("\nStep 5: Normalizing and scaling features...")
        node_features_scaled = self.normalize_and_scale_features(node_features_imputed, method='standard')
        edge_features_scaled = self.normalize_and_scale_features(edge_features_imputed, method='minmax')
        
        # Step 6: Analyze feature importance
        logger.info("\nStep 6: Analyzing feature importance...")
        feature_importance = self.analyze_feature_importance(node_features_scaled)
        
        # Step 7: Create visualizations
        logger.info("\nStep 7: Creating visualizations...")
        self.create_feature_visualizations()
        
        # Step 8: Save results
        logger.info("\nStep 8: Saving results...")
        self.save_processed_features()
        
        # Step 9: Summary
        logger.info("\n" + "=" * 60)
        logger.info("FEATURE ENGINEERING COMPLETE")
        logger.info("=" * 60)
        logger.info(f"Node features: {len(node_features_scaled)} nodes")
        logger.info(f"Edge features: {len(edge_features_scaled)} edges")
        logger.info(f"Feature importance analyzed for {len(feature_importance)} features")
        logger.info(f"Output directory: {self.output_dir}")
        
        return node_features_scaled, edge_features_scaled, feature_importance

def main():
    """Main function to run feature engineering."""
    engineer = GraphFeatureEngineer()
    engineer.run_feature_engineering()

if __name__ == "__main__":
    main() 