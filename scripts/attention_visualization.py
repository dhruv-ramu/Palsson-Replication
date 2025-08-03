#!/usr/bin/env python3
"""
Attention Visualization Tools for Metabolic Network Embedding

This module implements attention visualization tools for Phase 2 Week 2,
including attention weight plotting and analysis capabilities as specified in scope.md.

Features:
- Attention weight heatmaps
- Node relationship visualization
- Condition comparison plots
- Attention pattern analysis
- Interactive attention exploration

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn.functional as F
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
import json
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import logging
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AttentionVisualizer:
    """
    Comprehensive attention visualization toolkit for metabolic networks.
    
    Features:
    - Attention weight heatmaps
    - Node relationship graphs
    - Condition comparison visualizations
    - Attention pattern analysis
    - Interactive plots
    """
    
    def __init__(self, output_dir: str = "results/metabolic_network/attention_visualization"):
        """Initialize the attention visualizer."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Set up plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
        logger.info(f"Initialized AttentionVisualizer with output directory: {self.output_dir}")
    
    def plot_attention_heatmap(self, 
                              attention_weights: torch.Tensor,
                              node_labels: Optional[List[str]] = None,
                              title: str = "Attention Weights Heatmap",
                              save_path: Optional[str] = None) -> None:
        """
        Plot attention weights as a heatmap.
        
        Args:
            attention_weights: Attention weight tensor [num_heads, num_nodes, num_nodes]
            node_labels: Node labels for x and y axes
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating attention heatmap: {title}")
        
        # Convert to numpy and average across heads if multiple heads
        if attention_weights.dim() == 4:  # [batch, heads, nodes, nodes]
            attention_weights = attention_weights.squeeze(0)  # Remove batch dimension
        
        if attention_weights.dim() == 3:  # [heads, nodes, nodes]
            attention_weights = attention_weights.mean(dim=0)  # Average across heads
        
        attention_matrix = attention_weights.detach().cpu().numpy()
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Create heatmap
        im = ax.imshow(attention_matrix, cmap='viridis', aspect='auto')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Attention Weight', rotation=270, labelpad=20)
        
        # Set labels
        if node_labels is not None:
            # Sample labels if too many nodes
            if len(node_labels) > 50:
                step = len(node_labels) // 50
                sampled_labels = node_labels[::step]
                sampled_indices = np.arange(0, len(node_labels), step)
                ax.set_xticks(sampled_indices)
                ax.set_yticks(sampled_indices)
                ax.set_xticklabels(sampled_labels, rotation=45, ha='right')
                ax.set_yticklabels(sampled_labels)
            else:
                ax.set_xticks(range(len(node_labels)))
                ax.set_yticks(range(len(node_labels)))
                ax.set_xticklabels(node_labels, rotation=45, ha='right')
                ax.set_yticklabels(node_labels)
        else:
            ax.set_xlabel('Target Nodes')
            ax.set_ylabel('Source Nodes')
        
        # Set title
        ax.set_title(title, fontsize=16, fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"attention_heatmap_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Attention heatmap saved to: {save_path}")
    
    def plot_node_attention_network(self, 
                                   attention_weights: torch.Tensor,
                                   node_features: torch.Tensor,
                                   node_types: Optional[torch.Tensor] = None,
                                   top_k: int = 50,
                                   title: str = "Node Attention Network",
                                   save_path: Optional[str] = None) -> None:
        """
        Plot attention relationships as a network graph.
        
        Args:
            attention_weights: Attention weight tensor
            node_features: Node feature tensor
            node_types: Node type indicators (optional)
            top_k: Number of top attention connections to show
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating node attention network: {title}")
        
        # Process attention weights
        if attention_weights.dim() == 4:
            attention_weights = attention_weights.squeeze(0)
        if attention_weights.dim() == 3:
            attention_weights = attention_weights.mean(dim=0)
        
        attention_matrix = attention_weights.detach().cpu().numpy()
        
        # Create network graph
        G = nx.DiGraph()
        
        # Add nodes
        num_nodes = attention_matrix.shape[0]
        for i in range(num_nodes):
            node_type = "metabolite" if node_types is None or node_types[i, 0] > 0.5 else "reaction"
            G.add_node(i, type=node_type)
        
        # Add edges (top-k attention connections)
        flat_attention = attention_matrix.flatten()
        top_indices = np.argsort(flat_attention)[-top_k:]
        
        for idx in top_indices:
            source = idx // num_nodes
            target = idx % num_nodes
            weight = attention_matrix[source, target]
            if weight > 0.01:  # Only add significant connections
                G.add_edge(source, target, weight=weight)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(15, 12))
        
        # Position nodes using spring layout
        pos = nx.spring_layout(G, k=1, iterations=50)
        
        # Draw nodes
        node_colors = ['red' if G.nodes[node]['type'] == 'metabolite' else 'blue' for node in G.nodes()]
        node_sizes = [1000 if G.nodes[node]['type'] == 'metabolite' else 800 for node in G.nodes()]
        
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.7, ax=ax)
        
        # Draw edges
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        edge_alphas = [min(1.0, w * 10) for w in edge_weights]  # Scale for visibility
        
        nx.draw_networkx_edges(G, pos, alpha=edge_alphas, edge_color='gray', 
                              width=[w * 5 for w in edge_weights], ax=ax)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.7, label='Metabolites'),
            Patch(facecolor='blue', alpha=0.7, label='Reactions')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Set title
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.axis('off')
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"attention_network_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Attention network saved to: {save_path}")
    
    def plot_condition_comparison(self, 
                                 glucose_attention: torch.Tensor,
                                 acetate_attention: torch.Tensor,
                                 lactose_attention: torch.Tensor,
                                 node_labels: Optional[List[str]] = None,
                                 title: str = "Condition Comparison",
                                 save_path: Optional[str] = None) -> None:
        """
        Plot attention comparison across different growth conditions.
        
        Args:
            glucose_attention: Attention weights for glucose condition
            acetate_attention: Attention weights for acetate condition
            lactose_attention: Attention weights for lactose condition
            node_labels: Node labels
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating condition comparison plot: {title}")
        
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
        
        # Create subplots
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        
        conditions = [
            (glucose_matrix, 'Glucose', 'viridis'),
            (acetate_matrix, 'Acetate', 'plasma'),
            (lactose_matrix, 'Lactose', 'inferno')
        ]
        
        for i, (matrix, condition_name, cmap) in enumerate(conditions):
            im = axes[i].imshow(matrix, cmap=cmap, aspect='auto')
            axes[i].set_title(f'{condition_name} Condition', fontsize=14, fontweight='bold')
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=axes[i])
            cbar.set_label('Attention Weight', rotation=270, labelpad=20)
            
            # Set labels
            if node_labels is not None and len(node_labels) <= 30:
                axes[i].set_xticks(range(len(node_labels)))
                axes[i].set_yticks(range(len(node_labels)))
                axes[i].set_xticklabels(node_labels, rotation=45, ha='right', fontsize=8)
                axes[i].set_yticklabels(node_labels, fontsize=8)
            else:
                axes[i].set_xlabel('Target Nodes')
                axes[i].set_ylabel('Source Nodes')
        
        # Set main title
        fig.suptitle(title, fontsize=16, fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"condition_comparison_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Condition comparison saved to: {save_path}")
    
    def plot_attention_patterns(self, 
                               attention_weights: torch.Tensor,
                               node_features: torch.Tensor,
                               node_types: Optional[torch.Tensor] = None,
                               title: str = "Attention Patterns Analysis",
                               save_path: Optional[str] = None) -> None:
        """
        Analyze and visualize attention patterns.
        
        Args:
            attention_weights: Attention weight tensor
            node_features: Node feature tensor
            node_types: Node type indicators (optional)
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating attention patterns analysis: {title}")
        
        # Process attention weights
        if attention_weights.dim() == 4:
            attention_weights = attention_weights.squeeze(0)
        if attention_weights.dim() == 3:
            attention_weights = attention_weights.mean(dim=0)
        
        attention_matrix = attention_weights.detach().cpu().numpy()
        features = node_features.detach().cpu().numpy()
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Attention distribution
        axes[0, 0].hist(attention_matrix.flatten(), bins=50, alpha=0.7, color='skyblue')
        axes[0, 0].set_title('Attention Weight Distribution', fontweight='bold')
        axes[0, 0].set_xlabel('Attention Weight')
        axes[0, 0].set_ylabel('Frequency')
        
        # 2. Node attention centrality
        node_attention_in = attention_matrix.sum(axis=0)
        node_attention_out = attention_matrix.sum(axis=1)
        
        axes[0, 1].scatter(node_attention_in, node_attention_out, alpha=0.6, color='orange')
        axes[0, 1].set_title('Node Attention Centrality', fontweight='bold')
        axes[0, 1].set_xlabel('Incoming Attention')
        axes[0, 1].set_ylabel('Outgoing Attention')
        
        # 3. Feature-attention correlation
        if features.shape[1] > 1:
            # Use PCA for high-dimensional features
            pca = PCA(n_components=1)
            features_pca = pca.fit_transform(features).flatten()
        else:
            features_pca = features.flatten()
        
        axes[1, 0].scatter(features_pca, node_attention_in, alpha=0.6, color='green')
        axes[1, 0].set_title('Feature-Attention Correlation', fontweight='bold')
        axes[1, 0].set_xlabel('Principal Component 1')
        axes[1, 0].set_ylabel('Incoming Attention')
        
        # 4. Node type attention analysis
        if node_types is not None:
            node_types_np = node_types.detach().cpu().numpy()
            metabolite_mask = node_types_np[:, 0] > 0.5
            reaction_mask = ~metabolite_mask
            
            metabolite_attention = node_attention_in[metabolite_mask]
            reaction_attention = node_attention_in[reaction_mask]
            
            axes[1, 1].boxplot([metabolite_attention, reaction_attention], 
                              labels=['Metabolites', 'Reactions'])
            axes[1, 1].set_title('Attention by Node Type', fontweight='bold')
            axes[1, 1].set_ylabel('Incoming Attention')
        else:
            # Fallback: attention sparsity
            attention_sparsity = (attention_matrix < 0.01).sum() / attention_matrix.size
            axes[1, 1].text(0.5, 0.5, f'Attention Sparsity: {attention_sparsity:.3f}', 
                           ha='center', va='center', transform=axes[1, 1].transAxes,
                           fontsize=14, fontweight='bold')
            axes[1, 1].set_title('Attention Sparsity', fontweight='bold')
            axes[1, 1].axis('off')
        
        # Set main title
        fig.suptitle(title, fontsize=16, fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"attention_patterns_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Attention patterns analysis saved to: {save_path}")
    
    def create_interactive_attention_plot(self, 
                                        attention_weights: torch.Tensor,
                                        node_labels: Optional[List[str]] = None,
                                        title: str = "Interactive Attention Visualization") -> None:
        """
        Create an interactive attention visualization using Plotly.
        
        Args:
            attention_weights: Attention weight tensor
            node_labels: Node labels
            title: Plot title
        """
        logger.info(f"Creating interactive attention plot: {title}")
        
        # Process attention weights
        if attention_weights.dim() == 4:
            attention_weights = attention_weights.squeeze(0)
        if attention_weights.dim() == 3:
            attention_weights = attention_weights.mean(dim=0)
        
        attention_matrix = attention_weights.detach().cpu().numpy()
        
        # Create interactive heatmap
        if node_labels is None:
            node_labels = [f'Node_{i}' for i in range(attention_matrix.shape[0])]
        
        # Sample labels if too many nodes
        if len(node_labels) > 100:
            step = len(node_labels) // 100
            sampled_labels = node_labels[::step]
            sampled_matrix = attention_matrix[::step, ::step]
        else:
            sampled_labels = node_labels
            sampled_matrix = attention_matrix
        
        # Create Plotly figure
        fig = go.Figure(data=go.Heatmap(
            z=sampled_matrix,
            x=sampled_labels,
            y=sampled_labels,
            colorscale='Viridis',
            colorbar=dict(title='Attention Weight')
        ))
        
        fig.update_layout(
            title=title,
            xaxis_title='Target Nodes',
            yaxis_title='Source Nodes',
            width=800,
            height=600
        )
        
        # Save interactive plot
        save_path = self.output_dir / f"interactive_attention_{self.timestamp}.html"
        fig.write_html(str(save_path))
        
        logger.info(f"Interactive attention plot saved to: {save_path}")
    
    def plot_attention_evolution(self, 
                                attention_weights_list: List[torch.Tensor],
                                layer_names: Optional[List[str]] = None,
                                title: str = "Attention Evolution Across Layers",
                                save_path: Optional[str] = None) -> None:
        """
        Plot attention evolution across different layers.
        
        Args:
            attention_weights_list: List of attention weights from different layers
            layer_names: Names of the layers
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating attention evolution plot: {title}")
        
        if layer_names is None:
            layer_names = [f'Layer_{i}' for i in range(len(attention_weights_list))]
        
        # Process attention weights
        processed_weights = []
        for attn in attention_weights_list:
            if attn.dim() == 4:
                attn = attn.squeeze(0)
            if attn.dim() == 3:
                attn = attn.mean(dim=0)
            processed_weights.append(attn.detach().cpu().numpy())
        
        # Create subplots
        num_layers = len(processed_weights)
        fig, axes = plt.subplots(1, num_layers, figsize=(5*num_layers, 5))
        
        if num_layers == 1:
            axes = [axes]
        
        for i, (weights, layer_name) in enumerate(zip(processed_weights, layer_names)):
            im = axes[i].imshow(weights, cmap='viridis', aspect='auto')
            axes[i].set_title(f'{layer_name}', fontweight='bold')
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=axes[i])
            cbar.set_label('Attention Weight', rotation=270, labelpad=20)
            
            if i == 0:
                axes[i].set_ylabel('Source Nodes')
            axes[i].set_xlabel('Target Nodes')
        
        # Set main title
        fig.suptitle(title, fontsize=16, fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"attention_evolution_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Attention evolution plot saved to: {save_path}")
    
    def generate_attention_report(self, 
                                attention_weights: torch.Tensor,
                                node_features: torch.Tensor,
                                node_types: Optional[torch.Tensor] = None,
                                save_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Generate a comprehensive attention analysis report.
        
        Args:
            attention_weights: Attention weight tensor
            node_features: Node feature tensor
            node_types: Node type indicators (optional)
            save_path: Path to save the report
            
        Returns:
            Dictionary containing attention statistics
        """
        logger.info("Generating attention analysis report...")
        
        # Process attention weights
        if attention_weights.dim() == 4:
            attention_weights = attention_weights.squeeze(0)
        if attention_weights.dim() == 3:
            attention_weights = attention_weights.mean(dim=0)
        
        attention_matrix = attention_weights.detach().cpu().numpy()
        
        # Calculate attention statistics
        stats = {
            'total_attention': attention_matrix.sum(),
            'mean_attention': attention_matrix.mean(),
            'std_attention': attention_matrix.std(),
            'max_attention': attention_matrix.max(),
            'min_attention': attention_matrix.min(),
            'attention_sparsity': (attention_matrix < 0.01).sum() / attention_matrix.size,
            'attention_entropy': -np.sum(attention_matrix * np.log(attention_matrix + 1e-10)),
            'num_nodes': attention_matrix.shape[0],
            'num_connections': (attention_matrix > 0.01).sum()
        }
        
        # Node-level statistics
        node_attention_in = attention_matrix.sum(axis=0)
        node_attention_out = attention_matrix.sum(axis=1)
        
        stats.update({
            'mean_incoming_attention': node_attention_in.mean(),
            'std_incoming_attention': node_attention_in.std(),
            'mean_outgoing_attention': node_attention_out.mean(),
            'std_outgoing_attention': node_attention_out.std(),
            'top_attended_nodes': np.argsort(node_attention_in)[-10:].tolist(),
            'top_attending_nodes': np.argsort(node_attention_out)[-10:].tolist()
        })
        
        # Node type analysis
        if node_types is not None:
            node_types_np = node_types.detach().cpu().numpy()
            metabolite_mask = node_types_np[:, 0] > 0.5
            reaction_mask = ~metabolite_mask
            
            metabolite_attention = node_attention_in[metabolite_mask]
            reaction_attention = node_attention_in[reaction_mask]
            
            stats.update({
                'metabolite_mean_attention': metabolite_attention.mean(),
                'reaction_mean_attention': reaction_attention.mean(),
                'metabolite_std_attention': metabolite_attention.std(),
                'reaction_std_attention': reaction_attention.std(),
                'num_metabolites': metabolite_mask.sum(),
                'num_reactions': reaction_mask.sum()
            })
        
        # Save report
        if save_path is None:
            save_path = self.output_dir / f"attention_report_{self.timestamp}.json"
        else:
            save_path = self.output_dir / save_path
        
        with open(save_path, 'w') as f:
            json.dump(stats, f, indent=2)
        
        logger.info(f"Attention report saved to: {save_path}")
        logger.info(f"Attention statistics: {stats}")
        
        return stats

def test_attention_visualization():
    """Test the attention visualization tools."""
    logger.info("Testing attention visualization tools...")
    
    # Create test data
    num_nodes = 50
    num_heads = 8
    
    # Test attention weights
    attention_weights = torch.randn(1, num_heads, num_nodes, num_nodes)
    attention_weights = F.softmax(attention_weights, dim=-1)
    
    # Test node features
    node_features = torch.randn(num_nodes, 64)
    
    # Test node types
    node_types = torch.randint(0, 2, (num_nodes, 2)).float()
    
    # Initialize visualizer
    visualizer = AttentionVisualizer()
    
    # Test different visualization functions
    logger.info("Testing attention heatmap...")
    visualizer.plot_attention_heatmap(attention_weights, title="Test Attention Heatmap")
    
    logger.info("Testing node attention network...")
    visualizer.plot_node_attention_network(attention_weights, node_features, node_types, 
                                         title="Test Attention Network")
    
    logger.info("Testing condition comparison...")
    glucose_attn = torch.randn(1, num_heads, num_nodes, num_nodes)
    acetate_attn = torch.randn(1, num_heads, num_nodes, num_nodes)
    lactose_attn = torch.randn(1, num_heads, num_nodes, num_nodes)
    
    for attn in [glucose_attn, acetate_attn, lactose_attn]:
        attn = F.softmax(attn, dim=-1)
    
    visualizer.plot_condition_comparison(glucose_attn, acetate_attn, lactose_attn,
                                       title="Test Condition Comparison")
    
    logger.info("Testing attention patterns...")
    visualizer.plot_attention_patterns(attention_weights, node_features, node_types,
                                     title="Test Attention Patterns")
    
    logger.info("Testing interactive plot...")
    visualizer.create_interactive_attention_plot(attention_weights, title="Test Interactive Plot")
    
    logger.info("Testing attention evolution...")
    attention_weights_list = [attention_weights, attention_weights * 1.5, attention_weights * 0.5]
    layer_names = ['Layer 1', 'Layer 2', 'Layer 3']
    visualizer.plot_attention_evolution(attention_weights_list, layer_names,
                                      title="Test Attention Evolution")
    
    logger.info("Testing attention report...")
    stats = visualizer.generate_attention_report(attention_weights, node_features, node_types)
    
    logger.info("All attention visualization tests completed successfully!")

if __name__ == "__main__":
    test_attention_visualization() 