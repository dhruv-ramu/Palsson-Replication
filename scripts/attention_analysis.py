#!/usr/bin/env python3
"""
Attention Pattern Analysis for Metabolic Network Embedding

This module implements attention pattern analysis for Phase 2 Week 2,
including comprehensive analysis capabilities as specified in scope.md.

Features:
- Attention pattern analysis
- Node relationship analysis
- Condition-specific patterns
- Metabolic pathway attention analysis
- Statistical attention analysis

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
import json
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import logging
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA, NMF
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from scipy.stats import pearsonr, spearmanr, entropy
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AttentionAnalyzer:
    """
    Comprehensive attention pattern analyzer for metabolic networks.
    
    Features:
    - Attention pattern clustering
    - Node relationship analysis
    - Condition-specific pattern analysis
    - Metabolic pathway attention analysis
    - Statistical attention analysis
    """
    
    def __init__(self, output_dir: str = "results/metabolic_network/attention_analysis"):
        """Initialize the attention analyzer."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        logger.info(f"Initialized AttentionAnalyzer with output directory: {self.output_dir}")
    
    def analyze_attention_patterns(self, 
                                 attention_weights: torch.Tensor,
                                 node_features: torch.Tensor,
                                 node_types: Optional[torch.Tensor] = None,
                                 node_labels: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Analyze attention patterns comprehensively.
        
        Args:
            attention_weights: Attention weight tensor
            node_features: Node feature tensor
            node_types: Node type indicators (optional)
            node_labels: Node labels (optional)
            
        Returns:
            Dictionary containing analysis results
        """
        logger.info("Analyzing attention patterns...")
        
        # Process attention weights
        if attention_weights.dim() == 4:
            attention_weights = attention_weights.squeeze(0)
        if attention_weights.dim() == 3:
            attention_weights = attention_weights.mean(dim=0)
        
        attention_matrix = attention_weights.detach().cpu().numpy()
        features = node_features.detach().cpu().numpy()
        
        # Initialize results dictionary
        results = {
            'attention_matrix_shape': attention_matrix.shape,
            'node_features_shape': features.shape,
            'analysis_timestamp': self.timestamp
        }
        
        # 1. Basic attention statistics
        results['basic_stats'] = self._compute_basic_attention_stats(attention_matrix)
        
        # 2. Attention clustering analysis
        results['clustering'] = self._analyze_attention_clustering(attention_matrix)
        
        # 3. Node relationship analysis
        results['node_relationships'] = self._analyze_node_relationships(attention_matrix, features)
        
        # 4. Node type analysis (if available)
        if node_types is not None:
            results['node_type_analysis'] = self._analyze_node_type_patterns(
                attention_matrix, node_types.detach().cpu().numpy()
            )
        
        # 5. Feature-attention correlation analysis
        results['feature_correlation'] = self._analyze_feature_attention_correlation(
            attention_matrix, features
        )
        
        # 6. Attention centrality analysis
        results['centrality'] = self._analyze_attention_centrality(attention_matrix)
        
        # 7. Attention sparsity analysis
        results['sparsity'] = self._analyze_attention_sparsity(attention_matrix)
        
        # 8. Attention entropy analysis
        results['entropy'] = self._analyze_attention_entropy(attention_matrix)
        
        return results
    
    def _compute_basic_attention_stats(self, attention_matrix: np.ndarray) -> Dict[str, float]:
        """Compute basic attention statistics."""
        return {
            'total_attention': float(attention_matrix.sum()),
            'mean_attention': float(attention_matrix.mean()),
            'std_attention': float(attention_matrix.std()),
            'max_attention': float(attention_matrix.max()),
            'min_attention': float(attention_matrix.min()),
            'median_attention': float(np.median(attention_matrix)),
            'attention_range': float(attention_matrix.max() - attention_matrix.min()),
            'attention_variance': float(attention_matrix.var())
        }
    
    def _analyze_attention_clustering(self, attention_matrix: np.ndarray) -> Dict[str, Any]:
        """Analyze attention patterns using clustering."""
        logger.info("Performing attention clustering analysis...")
        
        # Reshape attention matrix for clustering
        num_nodes = attention_matrix.shape[0]
        attention_vectors = attention_matrix.reshape(num_nodes, -1)
        
        # PCA for dimensionality reduction
        pca = PCA(n_components=min(10, attention_vectors.shape[1]))
        attention_pca = pca.fit_transform(attention_vectors)
        
        # K-means clustering
        kmeans = KMeans(n_clusters=min(5, num_nodes // 10), random_state=42)
        kmeans_labels = kmeans.fit_predict(attention_pca)
        
        # DBSCAN clustering
        dbscan = DBSCAN(eps=0.5, min_samples=5)
        dbscan_labels = dbscan.fit_predict(attention_pca)
        
        # Compute clustering metrics
        kmeans_silhouette = silhouette_score(attention_pca, kmeans_labels) if len(set(kmeans_labels)) > 1 else 0
        kmeans_calinski = calinski_harabasz_score(attention_pca, kmeans_labels) if len(set(kmeans_labels)) > 1 else 0
        
        return {
            'pca_explained_variance_ratio': pca.explained_variance_ratio_.tolist(),
            'kmeans_labels': kmeans_labels.tolist(),
            'dbscan_labels': dbscan_labels.tolist(),
            'kmeans_silhouette_score': float(kmeans_silhouette),
            'kmeans_calinski_harabasz_score': float(kmeans_calinski),
            'num_kmeans_clusters': len(set(kmeans_labels)),
            'num_dbscan_clusters': len(set(dbscan_labels)) - (1 if -1 in dbscan_labels else 0)
        }
    
    def _analyze_node_relationships(self, attention_matrix: np.ndarray, features: np.ndarray) -> Dict[str, Any]:
        """Analyze node relationships based on attention patterns."""
        logger.info("Analyzing node relationships...")
        
        # Compute attention-based similarity matrix
        attention_similarity = np.corrcoef(attention_matrix)
        
        # Find highly connected node pairs
        threshold = np.percentile(attention_matrix, 95)
        high_attention_pairs = np.where(attention_matrix > threshold)
        
        # Compute feature similarity
        feature_similarity = np.corrcoef(features)
        
        # Correlation between attention and feature similarity
        attention_feature_corr = np.corrcoef(attention_matrix.flatten(), feature_similarity.flatten())[0, 1]
        
        return {
            'attention_similarity_mean': float(attention_similarity.mean()),
            'attention_similarity_std': float(attention_similarity.std()),
            'feature_similarity_mean': float(feature_similarity.mean()),
            'feature_similarity_std': float(feature_similarity.std()),
            'attention_feature_correlation': float(attention_feature_corr),
            'high_attention_threshold': float(threshold),
            'num_high_attention_pairs': len(high_attention_pairs[0]),
            'high_attention_pairs': list(zip(high_attention_pairs[0].tolist(), high_attention_pairs[1].tolist()))
        }
    
    def _analyze_node_type_patterns(self, attention_matrix: np.ndarray, node_types: np.ndarray) -> Dict[str, Any]:
        """Analyze attention patterns by node type."""
        logger.info("Analyzing node type patterns...")
        
        # Separate metabolites and reactions
        metabolite_mask = node_types[:, 0] > 0.5
        reaction_mask = ~metabolite_mask
        
        # Compute attention statistics by node type
        metabolite_attention = attention_matrix[metabolite_mask]
        reaction_attention = attention_matrix[reaction_mask]
        
        # Cross-type attention
        metabolite_to_reaction = attention_matrix[metabolite_mask][:, reaction_mask]
        reaction_to_metabolite = attention_matrix[reaction_mask][:, metabolite_mask]
        
        return {
            'num_metabolites': int(metabolite_mask.sum()),
            'num_reactions': int(reaction_mask.sum()),
            'metabolite_attention_mean': float(metabolite_attention.mean()),
            'reaction_attention_mean': float(reaction_attention.mean()),
            'metabolite_attention_std': float(metabolite_attention.std()),
            'reaction_attention_std': float(reaction_attention.std()),
            'metabolite_to_reaction_mean': float(metabolite_to_reaction.mean()),
            'reaction_to_metabolite_mean': float(reaction_to_metabolite.mean()),
            'metabolite_to_reaction_std': float(metabolite_to_reaction.std()),
            'reaction_to_metabolite_std': float(reaction_to_metabolite.std())
        }
    
    def _analyze_feature_attention_correlation(self, attention_matrix: np.ndarray, features: np.ndarray) -> Dict[str, Any]:
        """Analyze correlation between node features and attention patterns."""
        logger.info("Analyzing feature-attention correlations...")
        
        # Compute node-level attention statistics
        node_attention_in = attention_matrix.sum(axis=0)
        node_attention_out = attention_matrix.sum(axis=1)
        
        # Compute correlations with each feature
        feature_correlations = {}
        for i in range(features.shape[1]):
            feature = features[:, i]
            
            # Pearson correlation
            pearson_in, p_in = pearsonr(feature, node_attention_in)
            pearson_out, p_out = pearsonr(feature, node_attention_out)
            
            # Spearman correlation
            spearman_in, sp_in = spearmanr(feature, node_attention_in)
            spearman_out, sp_out = spearmanr(feature, node_attention_out)
            
            feature_correlations[f'feature_{i}'] = {
                'pearson_incoming': float(pearson_in),
                'pearson_incoming_pvalue': float(p_in),
                'pearson_outgoing': float(pearson_out),
                'pearson_outgoing_pvalue': float(p_out),
                'spearman_incoming': float(spearman_in),
                'spearman_incoming_pvalue': float(sp_in),
                'spearman_outgoing': float(spearman_out),
                'spearman_outgoing_pvalue': float(sp_out)
            }
        
        return {
            'feature_correlations': feature_correlations,
            'num_features': features.shape[1]
        }
    
    def _analyze_attention_centrality(self, attention_matrix: np.ndarray) -> Dict[str, Any]:
        """Analyze attention centrality patterns."""
        logger.info("Analyzing attention centrality...")
        
        # Compute centrality measures
        node_attention_in = attention_matrix.sum(axis=0)
        node_attention_out = attention_matrix.sum(axis=1)
        
        # Find most central nodes
        top_incoming = np.argsort(node_attention_in)[-10:]
        top_outgoing = np.argsort(node_attention_out)[-10:]
        
        # Compute centrality statistics
        return {
            'incoming_attention_mean': float(node_attention_in.mean()),
            'incoming_attention_std': float(node_attention_in.std()),
            'outgoing_attention_mean': float(node_attention_out.mean()),
            'outgoing_attention_std': float(node_attention_out.std()),
            'top_incoming_nodes': top_incoming.tolist(),
            'top_outgoing_nodes': top_outgoing.tolist(),
            'top_incoming_scores': node_attention_in[top_incoming].tolist(),
            'top_outgoing_scores': node_attention_out[top_outgoing].tolist(),
            'centrality_correlation': float(np.corrcoef(node_attention_in, node_attention_out)[0, 1])
        }
    
    def _analyze_attention_sparsity(self, attention_matrix: np.ndarray) -> Dict[str, Any]:
        """Analyze attention sparsity patterns."""
        logger.info("Analyzing attention sparsity...")
        
        # Compute sparsity at different thresholds
        thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
        sparsity_at_thresholds = {}
        
        for threshold in thresholds:
            sparsity = (attention_matrix < threshold).sum() / attention_matrix.size
            sparsity_at_thresholds[f'threshold_{threshold}'] = float(sparsity)
        
        # Compute effective sparsity (using entropy)
        attention_entropy = entropy(attention_matrix.flatten() + 1e-10)
        max_entropy = entropy(np.ones_like(attention_matrix.flatten()) / attention_matrix.size)
        effective_sparsity = 1 - (attention_entropy / max_entropy)
        
        return {
            'sparsity_at_thresholds': sparsity_at_thresholds,
            'effective_sparsity': float(effective_sparsity),
            'attention_entropy': float(attention_entropy),
            'max_entropy': float(max_entropy)
        }
    
    def _analyze_attention_entropy(self, attention_matrix: np.ndarray) -> Dict[str, Any]:
        """Analyze attention entropy patterns."""
        logger.info("Analyzing attention entropy...")
        
        # Compute entropy for each node's attention distribution
        node_entropies = []
        for i in range(attention_matrix.shape[0]):
            node_dist = attention_matrix[i, :]
            node_entropy = entropy(node_dist + 1e-10)
            node_entropies.append(node_entropy)
        
        node_entropies = np.array(node_entropies)
        
        return {
            'mean_node_entropy': float(node_entropies.mean()),
            'std_node_entropy': float(node_entropies.std()),
            'min_node_entropy': float(node_entropies.min()),
            'max_node_entropy': float(node_entropies.max()),
            'node_entropies': node_entropies.tolist(),
            'high_entropy_nodes': np.argsort(node_entropies)[-10:].tolist(),
            'low_entropy_nodes': np.argsort(node_entropies)[:10].tolist()
        }
    
    def analyze_condition_comparison(self, 
                                   glucose_attention: torch.Tensor,
                                   acetate_attention: torch.Tensor,
                                   lactose_attention: torch.Tensor) -> Dict[str, Any]:
        """
        Analyze attention patterns across different growth conditions.
        
        Args:
            glucose_attention: Attention weights for glucose condition
            acetate_attention: Attention weights for acetate condition
            lactose_attention: Attention weights for lactose condition
            
        Returns:
            Dictionary containing condition comparison analysis
        """
        logger.info("Analyzing attention patterns across conditions...")
        
        # Process attention weights
        def process_attention(attn):
            if attn.dim() == 4:
                attn = attn.squeeze(0)
            if attn.dim() == 3:
                attn = attn.mean(dim=0)
            return attn.detach().cpu().numpy()
        
        glucose_matrix = process_attention(glucose_attention)
        acetate_matrix = process_attention(acetate_attention)
        lactose_matrix = process_attention(lactose_attention)
        
        # Initialize results
        results = {
            'condition_comparison': {},
            'condition_similarities': {},
            'condition_differences': {}
        }
        
        # Analyze each condition
        conditions = {
            'glucose': glucose_matrix,
            'acetate': acetate_matrix,
            'lactose': lactose_matrix
        }
        
        for condition_name, attention_matrix in conditions.items():
            results['condition_comparison'][condition_name] = self._compute_basic_attention_stats(attention_matrix)
        
        # Compute condition similarities
        condition_names = list(conditions.keys())
        for i, name1 in enumerate(condition_names):
            for j, name2 in enumerate(condition_names[i+1:], i+1):
                matrix1 = conditions[name1]
                matrix2 = conditions[name2]
                
                # Correlation similarity
                correlation = np.corrcoef(matrix1.flatten(), matrix2.flatten())[0, 1]
                
                # Structural similarity (using cosine similarity)
                cosine_sim = np.dot(matrix1.flatten(), matrix2.flatten()) / (
                    np.linalg.norm(matrix1.flatten()) * np.linalg.norm(matrix2.flatten())
                )
                
                # Mean absolute difference
                mean_diff = np.mean(np.abs(matrix1 - matrix2))
                
                pair_name = f"{name1}_vs_{name2}"
                results['condition_similarities'][pair_name] = {
                    'correlation': float(correlation),
                    'cosine_similarity': float(cosine_sim),
                    'mean_absolute_difference': float(mean_diff)
                }
        
        # Find condition-specific patterns
        for condition_name, attention_matrix in conditions.items():
            # Find unique patterns in this condition
            other_matrices = [m for name, m in conditions.items() if name != condition_name]
            other_mean = np.mean(other_matrices, axis=0)
            
            # Compute condition-specific attention
            condition_specific = attention_matrix - other_mean
            condition_specific_positive = np.where(condition_specific > 0, condition_specific, 0)
            
            results['condition_differences'][condition_name] = {
                'mean_condition_specific_attention': float(condition_specific.mean()),
                'max_condition_specific_attention': float(condition_specific.max()),
                'num_condition_specific_connections': int((condition_specific > 0.01).sum()),
                'condition_specific_attention_sum': float(condition_specific_positive.sum())
            }
        
        return results
    
    def analyze_attention_evolution(self, 
                                  attention_weights_list: List[torch.Tensor],
                                  layer_names: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Analyze attention evolution across different layers.
        
        Args:
            attention_weights_list: List of attention weights from different layers
            layer_names: Names of the layers
            
        Returns:
            Dictionary containing evolution analysis
        """
        logger.info("Analyzing attention evolution across layers...")
        
        if layer_names is None:
            layer_names = [f'Layer_{i}' for i in range(len(attention_weights_list))]
        
        # Process attention weights
        processed_matrices = []
        for attn in attention_weights_list:
            if attn.dim() == 4:
                attn = attn.squeeze(0)
            if attn.dim() == 3:
                attn = attn.mean(dim=0)
            processed_matrices.append(attn.detach().cpu().numpy())
        
        # Initialize results
        results = {
            'layer_analysis': {},
            'evolution_patterns': {},
            'layer_similarities': {}
        }
        
        # Analyze each layer
        for i, (matrix, layer_name) in enumerate(zip(processed_matrices, layer_names)):
            results['layer_analysis'][layer_name] = self._compute_basic_attention_stats(matrix)
        
        # Compute layer similarities
        for i in range(len(processed_matrices) - 1):
            matrix1 = processed_matrices[i]
            matrix2 = processed_matrices[i + 1]
            
            # Correlation between consecutive layers
            correlation = np.corrcoef(matrix1.flatten(), matrix2.flatten())[0, 1]
            
            # Attention change magnitude
            change_magnitude = np.mean(np.abs(matrix2 - matrix1))
            
            pair_name = f"{layer_names[i]}_to_{layer_names[i+1]}"
            results['layer_similarities'][pair_name] = {
                'correlation': float(correlation),
                'change_magnitude': float(change_magnitude),
                'relative_change': float(change_magnitude / (matrix1.mean() + 1e-8))
            }
        
        # Analyze evolution patterns
        # 1. Attention concentration evolution
        attention_concentration = []
        for matrix in processed_matrices:
            # Compute Gini coefficient-like measure
            sorted_attention = np.sort(matrix.flatten())
            n = len(sorted_attention)
            concentration = 2 * np.sum((np.arange(1, n + 1) * sorted_attention)) / (n * np.sum(sorted_attention)) - (n + 1) / n
            attention_concentration.append(concentration)
        
        results['evolution_patterns']['attention_concentration'] = attention_concentration
        
        # 2. Attention sparsity evolution
        attention_sparsity = []
        for matrix in processed_matrices:
            sparsity = (matrix < 0.01).sum() / matrix.size
            attention_sparsity.append(sparsity)
        
        results['evolution_patterns']['attention_sparsity'] = attention_sparsity
        
        # 3. Attention entropy evolution
        attention_entropy = []
        for matrix in processed_matrices:
            entropy_val = entropy(matrix.flatten() + 1e-10)
            attention_entropy.append(entropy_val)
        
        results['evolution_patterns']['attention_entropy'] = attention_entropy
        
        return results
    
    def save_analysis_results(self, results: Dict[str, Any], filename: Optional[str] = None) -> str:
        """
        Save analysis results to JSON file.
        
        Args:
            results: Analysis results dictionary
            filename: Output filename (optional)
            
        Returns:
            Path to saved file
        """
        if filename is None:
            filename = f"attention_analysis_{self.timestamp}.json"
        
        filepath = self.output_dir / filename
        
        with open(filepath, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Analysis results saved to: {filepath}")
        return str(filepath)
    
    def generate_analysis_summary(self, results: Dict[str, Any]) -> str:
        """
        Generate a human-readable summary of analysis results.
        
        Args:
            results: Analysis results dictionary
            
        Returns:
            Summary string
        """
        summary = []
        summary.append("=" * 60)
        summary.append("ATTENTION PATTERN ANALYSIS SUMMARY")
        summary.append("=" * 60)
        
        # Basic statistics
        if 'basic_stats' in results:
            stats = results['basic_stats']
            summary.append(f"\nBasic Attention Statistics:")
            summary.append(f"  - Mean attention: {stats['mean_attention']:.4f}")
            summary.append(f"  - Std attention: {stats['std_attention']:.4f}")
            summary.append(f"  - Attention range: {stats['attention_range']:.4f}")
        
        # Clustering results
        if 'clustering' in results:
            clustering = results['clustering']
            summary.append(f"\nClustering Analysis:")
            summary.append(f"  - K-means clusters: {clustering['num_kmeans_clusters']}")
            summary.append(f"  - K-means silhouette score: {clustering['kmeans_silhouette_score']:.4f}")
            summary.append(f"  - DBSCAN clusters: {clustering['num_dbscan_clusters']}")
        
        # Centrality analysis
        if 'centrality' in results:
            centrality = results['centrality']
            summary.append(f"\nCentrality Analysis:")
            summary.append(f"  - Mean incoming attention: {centrality['incoming_attention_mean']:.4f}")
            summary.append(f"  - Mean outgoing attention: {centrality['outgoing_attention_mean']:.4f}")
            summary.append(f"  - Centrality correlation: {centrality['centrality_correlation']:.4f}")
        
        # Sparsity analysis
        if 'sparsity' in results:
            sparsity = results['sparsity']
            summary.append(f"\nSparsity Analysis:")
            summary.append(f"  - Effective sparsity: {sparsity['effective_sparsity']:.4f}")
            summary.append(f"  - Sparsity at 0.01 threshold: {sparsity['sparsity_at_thresholds']['threshold_0.01']:.4f}")
        
        # Node type analysis
        if 'node_type_analysis' in results:
            node_type = results['node_type_analysis']
            summary.append(f"\nNode Type Analysis:")
            summary.append(f"  - Metabolites: {node_type['num_metabolites']}")
            summary.append(f"  - Reactions: {node_type['num_reactions']}")
            summary.append(f"  - Metabolite attention mean: {node_type['metabolite_attention_mean']:.4f}")
            summary.append(f"  - Reaction attention mean: {node_type['reaction_attention_mean']:.4f}")
        
        summary.append("\n" + "=" * 60)
        
        return "\n".join(summary)

def test_attention_analysis():
    """Test the attention analysis tools."""
    logger.info("Testing attention analysis tools...")
    
    # Create test data
    num_nodes = 100
    num_heads = 8
    
    # Test attention weights
    attention_weights = torch.randn(1, num_heads, num_nodes, num_nodes)
    attention_weights = F.softmax(attention_weights, dim=-1)
    
    # Test node features
    node_features = torch.randn(num_nodes, 64)
    
    # Test node types
    node_types = torch.randint(0, 2, (num_nodes, 2)).float()
    
    # Initialize analyzer
    analyzer = AttentionAnalyzer()
    
    # Test attention pattern analysis
    logger.info("Testing attention pattern analysis...")
    results = analyzer.analyze_attention_patterns(attention_weights, node_features, node_types)
    
    # Test condition comparison
    logger.info("Testing condition comparison analysis...")
    glucose_attn = torch.randn(1, num_heads, num_nodes, num_nodes)
    acetate_attn = torch.randn(1, num_heads, num_nodes, num_nodes)
    lactose_attn = torch.randn(1, num_heads, num_nodes, num_nodes)
    
    for attn in [glucose_attn, acetate_attn, lactose_attn]:
        attn = F.softmax(attn, dim=-1)
    
    condition_results = analyzer.analyze_condition_comparison(glucose_attn, acetate_attn, lactose_attn)
    
    # Test attention evolution
    logger.info("Testing attention evolution analysis...")
    attention_weights_list = [attention_weights, attention_weights * 1.5, attention_weights * 0.5]
    evolution_results = analyzer.analyze_attention_evolution(attention_weights_list)
    
    # Save results
    logger.info("Saving analysis results...")
    analyzer.save_analysis_results(results, "test_attention_analysis.json")
    
    # Generate summary
    logger.info("Generating analysis summary...")
    summary = analyzer.generate_analysis_summary(results)
    logger.info(f"\n{summary}")
    
    logger.info("All attention analysis tests completed successfully!")

if __name__ == "__main__":
    test_attention_analysis() 