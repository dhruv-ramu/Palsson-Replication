"""
Network Visualization Module

Advanced visualization tools for metabolic networks, including interactive plots,
attention weight visualization, and cross-condition comparison plots.

Author: Metabolic Network Embedding Project
Date: 2025
"""

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Dict, List, Tuple, Optional, Any
import json
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class NetworkVisualizer:
    """
    Advanced visualization tools for metabolic networks.
    
    Features:
    - Interactive network plots with Plotly
    - Attention weight visualization
    - Cross-condition comparison plots
    - 3D network representations
    - Statistical analysis plots
    """
    
    def __init__(self, graph: nx.DiGraph, node_features: Dict = None):
        """
        Initialize the network visualizer.
        
        Args:
            graph: NetworkX graph representing the metabolic network
            node_features: Dictionary of node features
        """
        self.graph = graph
        self.node_features = node_features or {}
        self.node_ids = list(self.graph.nodes())
        
        # Set up color schemes
        self.setup_color_schemes()
    
    def setup_color_schemes(self):
        """Set up color schemes for different node types and attributes."""
        self.colors = {
            'metabolite': '#3498db',  # Blue
            'reaction': '#e74c3c',    # Red
            'reactant_edge': '#2ecc71',  # Green
            'product_edge': '#f39c12',   # Orange
            'attention_high': '#e74c3c',  # Red for high attention
            'attention_medium': '#f39c12', # Orange for medium attention
            'attention_low': '#3498db'     # Blue for low attention
        }
        
        # Compartment colors
        self.compartment_colors = {
            'c': '#3498db',  # cytosol - blue
            'e': '#e74c3c',  # extracellular - red
            'p': '#2ecc71',  # periplasm - green
            'm': '#f39c12',  # mitochondria - orange
            'n': '#9b59b6',  # nucleus - purple
            'g': '#1abc9c',  # golgi - teal
            'r': '#34495e',  # endoplasmic reticulum - dark blue
            'x': '#e67e22',  # peroxisome - dark orange
            'v': '#95a5a6',  # vacuole - gray
            'l': '#16a085'   # lysosome - dark teal
        }
    
    def create_interactive_network(self, output_path: str = None, 
                                 max_nodes: int = 200,
                                 attention_weights: Dict = None):
        """
        Create an interactive network visualization using Plotly.
        
        Args:
            output_path: Path to save the HTML file
            max_nodes: Maximum number of nodes to display
            attention_weights: Dictionary of attention weights for nodes
        """
        logger.info("Creating interactive network visualization...")
        
        # Create subgraph for visualization
        if self.graph.number_of_nodes() > max_nodes:
            # Sample from largest connected component
            largest_cc = max(nx.weakly_connected_components(self.graph), key=len)
            sampled_nodes = list(largest_cc)[:max_nodes]
            viz_graph = self.graph.subgraph(sampled_nodes)
        else:
            viz_graph = self.graph
        
        # Get layout
        pos = nx.spring_layout(viz_graph, k=2, iterations=50)
        
        # Prepare node data
        node_x = []
        node_y = []
        node_text = []
        node_color = []
        node_size = []
        node_type = []
        
        for node in viz_graph.nodes():
            node_data = viz_graph.nodes[node]
            x, y = pos[node]
            
            node_x.append(x)
            node_y.append(y)
            
            # Node label
            name = node_data.get('name', node)
            node_type_label = node_data.get('node_type', 'unknown')
            node_text.append(f"{name}<br>Type: {node_type_label}<br>ID: {node}")
            
            # Node color based on type
            if node_data.get('node_type') == 'metabolite':
                node_color.append(self.colors['metabolite'])
                node_size.append(15)
                node_type.append('Metabolite')
            else:
                node_color.append(self.colors['reaction'])
                node_size.append(10)
                node_type.append('Reaction')
            
            # Apply attention weights if provided
            if attention_weights and node in attention_weights:
                weight = attention_weights[node]
                if weight > 0.7:
                    node_color[-1] = self.colors['attention_high']
                elif weight > 0.3:
                    node_color[-1] = self.colors['attention_medium']
                else:
                    node_color[-1] = self.colors['attention_low']
                node_size[-1] = max(5, node_size[-1] * weight * 2)
        
        # Prepare edge data
        edge_x = []
        edge_y = []
        edge_color = []
        
        for edge in viz_graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            
            edge_data = viz_graph.edges[edge]
            if edge_data.get('edge_type') == 'reactant':
                edge_color.extend([self.colors['reactant_edge']] * 3)
            else:
                edge_color.extend([self.colors['product_edge']] * 3)
        
        # Create the plot
        fig = go.Figure()
        
        # Add edges
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='gray'),
            hoverinfo='none',
            mode='lines',
            showlegend=False
        ))
        
        # Add nodes
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            text=node_text,
            textposition="top center",
            marker=dict(
                size=node_size,
                color=node_color,
                line=dict(width=2, color='white'),
                opacity=0.8
            ),
            showlegend=False
        ))
        
        # Update layout
        fig.update_layout(
            title=f'Interactive Metabolic Network<br>'
                  f'<sub>{viz_graph.number_of_nodes()} nodes, {viz_graph.number_of_edges()} edges</sub>',
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white'
        )
        
        # Add legend
        legend_traces = []
        for node_type, color in [('Metabolites', self.colors['metabolite']), 
                                ('Reactions', self.colors['reaction'])]:
            legend_traces.append(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=10, color=color),
                name=node_type,
                showlegend=True
            ))
        
        fig.add_traces(legend_traces)
        
        if output_path:
            fig.write_html(output_path)
            logger.info(f"Interactive network saved to: {output_path}")
        
        return fig
    
    def plot_attention_weights(self, attention_weights: Dict[str, float], 
                             output_path: str = None,
                             top_k: int = 50):
        """
        Create a visualization of attention weights.
        
        Args:
            attention_weights: Dictionary mapping node IDs to attention weights
            output_path: Path to save the plot
            top_k: Number of top attention nodes to show
        """
        logger.info("Creating attention weights visualization...")
        
        # Sort nodes by attention weight
        sorted_nodes = sorted(attention_weights.items(), 
                            key=lambda x: x[1], reverse=True)[:top_k]
        
        node_ids, weights = zip(*sorted_nodes)
        
        # Get node names
        node_names = []
        node_types = []
        for node_id in node_ids:
            if node_id in self.graph.nodes():
                node_data = self.graph.nodes[node_id]
                name = node_data.get('name', node_id)
                node_names.append(name[:20])  # Truncate long names
                node_types.append(node_data.get('node_type', 'unknown'))
            else:
                node_names.append(node_id[:20])
                node_types.append('unknown')
        
        # Create DataFrame
        df = pd.DataFrame({
            'Node': node_names,
            'Attention_Weight': weights,
            'Node_Type': node_types,
            'Node_ID': node_ids
        })
        
        # Create the plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12))
        
        # Bar plot of top attention weights
        colors = [self.colors['attention_high'] if w > 0.7 else 
                 self.colors['attention_medium'] if w > 0.3 else 
                 self.colors['attention_low'] for w in weights]
        
        bars = ax1.bar(range(len(weights)), weights, color=colors, alpha=0.7)
        ax1.set_xlabel('Node Rank')
        ax1.set_ylabel('Attention Weight')
        ax1.set_title(f'Top {top_k} Nodes by Attention Weight')
        ax1.set_xticks(range(len(weights)))
        ax1.set_xticklabels([f"{i+1}" for i in range(len(weights))], rotation=45)
        
        # Add value labels on bars
        for i, (bar, weight) in enumerate(zip(bars, weights)):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{weight:.3f}', ha='center', va='bottom', fontsize=8)
        
        # Distribution plot
        ax2.hist(weights, bins=20, alpha=0.7, color=self.colors['attention_medium'])
        ax2.set_xlabel('Attention Weight')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Distribution of Attention Weights')
        ax2.axvline(np.mean(weights), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(weights):.3f}')
        ax2.axvline(np.median(weights), color='green', linestyle='--', 
                   label=f'Median: {np.median(weights):.3f}')
        ax2.legend()
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Attention weights plot saved to: {output_path}")
        
        plt.show()
        
        return fig, df
    
    def plot_cross_condition_comparison(self, condition_data: Dict[str, Dict], 
                                      output_path: str = None):
        """
        Create comparison plots across different growth conditions.
        
        Args:
            condition_data: Dictionary mapping condition names to data
            output_path: Path to save the plot
        """
        logger.info("Creating cross-condition comparison plot...")
        
        conditions = list(condition_data.keys())
        n_conditions = len(conditions)
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()
        
        # Plot 1: Growth rates comparison
        growth_rates = [condition_data[c].get('growth_rate', 0) for c in conditions]
        bars1 = axes[0].bar(conditions, growth_rates, color=self.colors['metabolite'], alpha=0.7)
        axes[0].set_title('Growth Rates Across Conditions')
        axes[0].set_ylabel('Growth Rate (1/h)')
        axes[0].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, rate in zip(bars1, growth_rates):
            axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{rate:.2f}', ha='center', va='bottom')
        
        # Plot 2: Network density comparison
        densities = [condition_data[c].get('network_density', 0) for c in conditions]
        bars2 = axes[1].bar(conditions, densities, color=self.colors['reaction'], alpha=0.7)
        axes[1].set_title('Network Density Across Conditions')
        axes[1].set_ylabel('Network Density')
        axes[1].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, density in zip(bars2, densities):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.0001,
                        f'{density:.4f}', ha='center', va='bottom')
        
        # Plot 3: Node count comparison
        node_counts = [condition_data[c].get('total_nodes', 0) for c in conditions]
        bars3 = axes[2].bar(conditions, node_counts, color=self.colors['attention_medium'], alpha=0.7)
        axes[2].set_title('Total Nodes Across Conditions')
        axes[2].set_ylabel('Number of Nodes')
        axes[2].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, count in zip(bars3, node_counts):
            axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                        f'{count}', ha='center', va='bottom')
        
        # Plot 4: Edge count comparison
        edge_counts = [condition_data[c].get('total_edges', 0) for c in conditions]
        bars4 = axes[3].bar(conditions, edge_counts, color=self.colors['attention_low'], alpha=0.7)
        axes[3].set_title('Total Edges Across Conditions')
        axes[3].set_ylabel('Number of Edges')
        axes[3].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, count in zip(bars4, edge_counts):
            axes[3].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
                        f'{count}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Cross-condition comparison saved to: {output_path}")
        
        plt.show()
        
        return fig
    
    def create_3d_network(self, output_path: str = None, max_nodes: int = 100):
        """
        Create a 3D network visualization.
        
        Args:
            output_path: Path to save the HTML file
            max_nodes: Maximum number of nodes to display
        """
        logger.info("Creating 3D network visualization...")
        
        # Create subgraph for visualization
        if self.graph.number_of_nodes() > max_nodes:
            largest_cc = max(nx.weakly_connected_components(self.graph), key=len)
            sampled_nodes = list(largest_cc)[:max_nodes]
            viz_graph = self.graph.subgraph(sampled_nodes)
        else:
            viz_graph = self.graph
        
        # Get 3D layout
        pos = nx.spring_layout(viz_graph, k=2, iterations=50, dim=3)
        
        # Prepare node data
        node_x = []
        node_y = []
        node_z = []
        node_color = []
        node_size = []
        node_text = []
        
        for node in viz_graph.nodes():
            node_data = viz_graph.nodes[node]
            x, y, z = pos[node]
            
            node_x.append(x)
            node_y.append(y)
            node_z.append(z)
            
            # Node label
            name = node_data.get('name', node)
            node_type = node_data.get('node_type', 'unknown')
            node_text.append(f"{name}<br>Type: {node_type}")
            
            # Node color and size
            if node_data.get('node_type') == 'metabolite':
                node_color.append(self.colors['metabolite'])
                node_size.append(8)
            else:
                node_color.append(self.colors['reaction'])
                node_size.append(5)
        
        # Prepare edge data
        edge_x = []
        edge_y = []
        edge_z = []
        
        for edge in viz_graph.edges():
            x0, y0, z0 = pos[edge[0]]
            x1, y1, z1 = pos[edge[1]]
            
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_z.extend([z0, z1, None])
        
        # Create 3D plot
        fig = go.Figure()
        
        # Add edges
        fig.add_trace(go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='gray', width=1),
            hoverinfo='none',
            showlegend=False
        ))
        
        # Add nodes
        fig.add_trace(go.Scatter3d(
            x=node_x, y=node_y, z=node_z,
            mode='markers',
            marker=dict(
                size=node_size,
                color=node_color,
                opacity=0.8,
                line=dict(color='white', width=1)
            ),
            text=node_text,
            hoverinfo='text',
            showlegend=False
        ))
        
        # Update layout
        fig.update_layout(
            title=f'3D Metabolic Network<br>'
                  f'<sub>{viz_graph.number_of_nodes()} nodes, {viz_graph.number_of_edges()} edges</sub>',
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Z',
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.5)
                )
            ),
            margin=dict(l=0, r=0, b=0, t=30),
            showlegend=False
        )
        
        if output_path:
            fig.write_html(output_path)
            logger.info(f"3D network saved to: {output_path}")
        
        return fig
    
    def plot_network_statistics(self, output_path: str = None):
        """
        Create comprehensive network statistics plots.
        
        Args:
            output_path: Path to save the plot
        """
        logger.info("Creating network statistics plots...")
        
        # Calculate statistics
        stats = {
            'total_nodes': self.graph.number_of_nodes(),
            'total_edges': self.graph.number_of_edges(),
            'metabolite_nodes': len([n for n, d in self.graph.nodes(data=True) 
                                   if d['node_type'] == 'metabolite']),
            'reaction_nodes': len([n for n, d in self.graph.nodes(data=True) 
                                 if d['node_type'] == 'reaction']),
            'density': nx.density(self.graph),
            'connected_components': nx.number_weakly_connected_components(self.graph),
            'average_clustering': nx.average_clustering(self.graph.to_undirected())
        }
        
        # Degree distribution
        degrees = dict(nx.degree(self.graph))
        degree_values = list(degrees.values())
        
        # Create subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # Plot 1: Node type distribution
        node_types = [stats['metabolite_nodes'], stats['reaction_nodes']]
        labels = ['Metabolites', 'Reactions']
        colors = [self.colors['metabolite'], self.colors['reaction']]
        axes[0].pie(node_types, labels=labels, colors=colors, autopct='%1.1f%%')
        axes[0].set_title('Node Type Distribution')
        
        # Plot 2: Degree distribution histogram
        axes[1].hist(degree_values, bins=30, alpha=0.7, color=self.colors['attention_medium'])
        axes[1].set_xlabel('Degree')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('Degree Distribution')
        axes[1].axvline(np.mean(degree_values), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(degree_values):.1f}')
        axes[1].legend()
        
        # Plot 3: Degree distribution log-log
        degree_counts = pd.Series(degree_values).value_counts().sort_index()
        axes[2].loglog(degree_counts.index, degree_counts.values, 'o-', 
                      color=self.colors['attention_high'])
        axes[2].set_xlabel('Degree (log)')
        axes[2].set_ylabel('Count (log)')
        axes[2].set_title('Degree Distribution (Log-Log)')
        
        # Plot 4: Network statistics bar plot
        stat_names = ['Total Nodes', 'Total Edges', 'Density', 'Components', 'Clustering']
        stat_values = [stats['total_nodes'], stats['total_edges'], 
                      stats['density'], stats['connected_components'], 
                      stats['average_clustering']]
        bars = axes[3].bar(stat_names, stat_values, color=self.colors['attention_low'], alpha=0.7)
        axes[3].set_title('Network Statistics')
        axes[3].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, value in zip(bars, stat_values):
            axes[3].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(stat_values)*0.01,
                        f'{value:.3f}', ha='center', va='bottom')
        
        # Plot 5: Compartment distribution (if available)
        compartments = {}
        for node, data in self.graph.nodes(data=True):
            if data.get('node_type') == 'metabolite':
                comp = data.get('compartment', 'unknown')
                compartments[comp] = compartments.get(comp, 0) + 1
        
        if compartments:
            comp_names = list(compartments.keys())
            comp_values = list(compartments.values())
            comp_colors = [self.compartment_colors.get(c, '#95a5a6') for c in comp_names]
            axes[4].bar(comp_names, comp_values, color=comp_colors, alpha=0.7)
            axes[4].set_title('Metabolite Distribution by Compartment')
            axes[4].tick_params(axis='x', rotation=45)
        else:
            axes[4].text(0.5, 0.5, 'No compartment data available', 
                        ha='center', va='center', transform=axes[4].transAxes)
            axes[4].set_title('Compartment Distribution')
        
        # Plot 6: Subsystem distribution (if available)
        subsystems = {}
        for node, data in self.graph.nodes(data=True):
            if data.get('node_type') == 'reaction':
                subsystem = data.get('subsystem', 'unknown')
                subsystems[subsystem] = subsystems.get(subsystem, 0) + 1
        
        if subsystems:
            # Show top 10 subsystems
            top_subsystems = sorted(subsystems.items(), key=lambda x: x[1], reverse=True)[:10]
            sub_names, sub_values = zip(*top_subsystems)
            axes[5].barh(range(len(sub_names)), sub_values, color=self.colors['reaction'], alpha=0.7)
            axes[5].set_yticks(range(len(sub_names)))
            axes[5].set_yticklabels([name[:20] for name in sub_names])
            axes[5].set_title('Top 10 Reaction Subsystems')
        else:
            axes[5].text(0.5, 0.5, 'No subsystem data available', 
                        ha='center', va='center', transform=axes[5].transAxes)
            axes[5].set_title('Subsystem Distribution')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Network statistics saved to: {output_path}")
        
        plt.show()
        
        return fig, stats


def main():
    """Main function for testing the network visualization."""
    
    # This would typically be called after creating a metabolic network graph
    print("Network visualization module created successfully!")
    print("Use this module with a MetabolicNetworkGraph instance to create visualizations.")


if __name__ == "__main__":
    main() 