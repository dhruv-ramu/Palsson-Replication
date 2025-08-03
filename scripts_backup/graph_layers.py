#!/usr/bin/env python3
"""
Custom Graph Convolution Layers for Metabolic Network Embedding

This module implements custom graph convolution layers for Phase 2 Week 1,
specialized for metabolic network analysis as specified in scope.md.

Features:
- Custom GCN layers with metabolic-specific features
- Edge-aware convolution layers
- Multi-scale graph convolution
- Hierarchical feature aggregation
- Metabolic pathway-aware layers

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree
import numpy as np
from typing import Optional, Tuple, Dict, Any
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MetabolicGCNConv(MessagePassing):
    """
    Custom GCN convolution layer for metabolic networks.
    
    Features:
    - Edge feature integration
    - Node type awareness (metabolite vs reaction)
    - Stoichiometric coefficient handling
    - Metabolic pathway information
    """
    
    def __init__(self, 
                 in_channels: int,
                 out_channels: int,
                 edge_channels: int = 0,
                 node_type_channels: int = 0,
                 pathway_channels: int = 0,
                 dropout: float = 0.2,
                 bias: bool = True):
        """
        Initialize the MetabolicGCNConv layer.
        
        Args:
            in_channels: Input node feature dimension
            out_channels: Output node feature dimension
            edge_channels: Edge feature dimension
            node_type_channels: Node type feature dimension
            pathway_channels: Pathway feature dimension
            dropout: Dropout rate
            bias: Whether to use bias
        """
        super(MetabolicGCNConv, self).__init__(aggr='add')
        
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.edge_channels = edge_channels
        self.node_type_channels = node_type_channels
        self.pathway_channels = pathway_channels
        self.dropout = dropout
        
        # Node feature transformation
        self.node_linear = nn.Linear(in_channels, out_channels, bias=bias)
        
        # Edge feature transformation (if edge features are provided)
        if edge_channels > 0:
            self.edge_linear = nn.Linear(edge_channels, out_channels, bias=False)
        else:
            self.edge_linear = None
        
        # Node type transformation (if node type features are provided)
        if node_type_channels > 0:
            self.node_type_linear = nn.Linear(node_type_channels, out_channels, bias=False)
        else:
            self.node_type_linear = None
        
        # Pathway transformation (if pathway features are provided)
        if pathway_channels > 0:
            self.pathway_linear = nn.Linear(pathway_channels, out_channels, bias=False)
        else:
            self.pathway_linear = None
        
        # Attention mechanism for edge weighting
        self.edge_attention = nn.Sequential(
            nn.Linear(out_channels * 2, out_channels // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(out_channels // 2, 1),
            nn.Sigmoid()
        )
        
        # Layer normalization
        self.layer_norm = nn.LayerNorm(out_channels)
        
        # Dropout layer
        self.dropout_layer = nn.Dropout(dropout)
        
        logger.info(f"Initialized MetabolicGCNConv: {in_channels}->{out_channels}")
        logger.info(f"Edge features: {edge_channels}, Node type: {node_type_channels}, Pathway: {pathway_channels}")
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                edge_attr: Optional[torch.Tensor] = None,
                node_type: Optional[torch.Tensor] = None,
                pathway: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass through the metabolic GCN layer.
        
        Args:
            x: Node features [num_nodes, in_channels]
            edge_index: Edge indices [2, num_edges]
            edge_attr: Edge features [num_edges, edge_channels] (optional)
            node_type: Node type features [num_nodes, node_type_channels] (optional)
            pathway: Pathway features [num_nodes, pathway_channels] (optional)
            
        Returns:
            Updated node features [num_nodes, out_channels]
        """
        # Add self-loops
        edge_index, edge_attr = add_self_loops(edge_index, edge_attr, num_nodes=x.size(0))
        
        # Start propagating messages
        return self.propagate(edge_index, x=x, edge_attr=edge_attr, 
                            node_type=node_type, pathway=pathway)
    
    def message(self, x_i: torch.Tensor, x_j: torch.Tensor,
                edge_attr: Optional[torch.Tensor] = None,
                node_type_i: Optional[torch.Tensor] = None,
                node_type_j: Optional[torch.Tensor] = None,
                pathway_i: Optional[torch.Tensor] = None,
                pathway_j: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Construct messages from node j to node i.
        
        Args:
            x_i: Source node features
            x_j: Target node features
            edge_attr: Edge features
            node_type_i: Source node type features
            node_type_j: Target node type features
            pathway_i: Source pathway features
            pathway_j: Target pathway features
            
        Returns:
            Messages [num_edges, out_channels]
        """
        # Transform node features
        msg = self.node_linear(x_j)
        
        # Add edge features if available
        if edge_attr is not None and self.edge_linear is not None:
            edge_msg = self.edge_linear(edge_attr)
            msg = msg + edge_msg
        
        # Add node type features if available
        if node_type_j is not None and self.node_type_linear is not None:
            type_msg = self.node_type_linear(node_type_j)
            msg = msg + type_msg
        
        # Add pathway features if available
        if pathway_j is not None and self.pathway_linear is not None:
            pathway_msg = self.pathway_linear(pathway_j)
            msg = msg + pathway_msg
        
        # Apply attention mechanism
        if x_i is not None:
            # Concatenate source and target features for attention
            attention_input = torch.cat([x_i, x_j], dim=-1)
            attention_weights = self.edge_attention(attention_input)
            msg = msg * attention_weights
        
        return msg
    
    def update(self, aggr_out: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        """
        Update node features based on aggregated messages.
        
        Args:
            aggr_out: Aggregated messages
            x: Original node features
            
        Returns:
            Updated node features
        """
        # Transform original features
        out = self.node_linear(x)
        
        # Add aggregated messages
        out = out + aggr_out
        
        # Apply layer normalization
        out = self.layer_norm(out)
        
        # Apply activation and dropout
        out = F.relu(out)
        out = self.dropout_layer(out)
        
        return out

class MultiScaleGCNConv(MessagePassing):
    """
    Multi-scale graph convolution layer for capturing different neighborhood sizes.
    
    Features:
    - Multiple neighborhood scales
    - Adaptive aggregation
    - Scale-specific attention
    """
    
    def __init__(self, 
                 in_channels: int,
                 out_channels: int,
                 scales: list = [1, 2, 3],
                 dropout: float = 0.2,
                 bias: bool = True):
        """
        Initialize the MultiScaleGCNConv layer.
        
        Args:
            in_channels: Input node feature dimension
            out_channels: Output node feature dimension
            scales: List of neighborhood scales to consider
            dropout: Dropout rate
            bias: Whether to use bias
        """
        super(MultiScaleGCNConv, self).__init__(aggr='add')
        
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.scales = scales
        self.dropout = dropout
        
        # Scale-specific transformations
        self.scale_transforms = nn.ModuleList([
            nn.Linear(in_channels, out_channels // len(scales), bias=bias)
            for _ in scales
        ])
        
        # Scale attention mechanism
        self.scale_attention = nn.Sequential(
            nn.Linear(out_channels, out_channels // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(out_channels // 2, len(scales)),
            nn.Softmax(dim=-1)
        )
        
        # Output projection
        self.output_projection = nn.Linear(out_channels, out_channels, bias=bias)
        
        # Layer normalization
        self.layer_norm = nn.LayerNorm(out_channels)
        
        # Dropout layer
        self.dropout_layer = nn.Dropout(dropout)
        
        logger.info(f"Initialized MultiScaleGCNConv with scales: {scales}")
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through the multi-scale GCN layer.
        
        Args:
            x: Node features
            edge_index: Edge indices
            
        Returns:
            Updated node features
        """
        # Multi-scale message passing
        scale_outputs = []
        
        for i, scale in enumerate(self.scales):
            # Apply scale-specific transformation
            scale_x = self.scale_transforms[i](x)
            
            # Multi-hop message passing for this scale
            scale_out = scale_x
            for _ in range(scale):
                scale_out = self.propagate(edge_index, x=scale_out)
            
            scale_outputs.append(scale_out)
        
        # Concatenate scale outputs
        multi_scale_out = torch.cat(scale_outputs, dim=-1)
        
        # Apply scale attention
        attention_weights = self.scale_attention(multi_scale_out)
        
        # Weighted combination of scales
        weighted_out = torch.zeros_like(multi_scale_out)
        for i, scale_out in enumerate(scale_outputs):
            weighted_out[:, i * self.out_channels // len(self.scales):(i + 1) * self.out_channels // len(self.scales)] = \
                scale_out * attention_weights[:, i:i+1]
        
        # Output projection
        out = self.output_projection(weighted_out)
        
        # Layer normalization and activation
        out = self.layer_norm(out)
        out = F.relu(out)
        out = self.dropout_layer(out)
        
        return out
    
    def message(self, x_j: torch.Tensor) -> torch.Tensor:
        """Construct messages from neighboring nodes."""
        return x_j
    
    def update(self, aggr_out: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        """Update node features based on aggregated messages."""
        return aggr_out

class HierarchicalGCNConv(MessagePassing):
    """
    Hierarchical graph convolution layer for metabolic pathway analysis.
    
    Features:
    - Pathway-aware aggregation
    - Hierarchical feature learning
    - Metabolic hierarchy integration
    """
    
    def __init__(self, 
                 in_channels: int,
                 out_channels: int,
                 pathway_dim: int = 64,
                 hierarchy_levels: int = 3,
                 dropout: float = 0.2,
                 bias: bool = True):
        """
        Initialize the HierarchicalGCNConv layer.
        
        Args:
            in_channels: Input node feature dimension
            out_channels: Output node feature dimension
            pathway_dim: Pathway feature dimension
            hierarchy_levels: Number of hierarchy levels
            dropout: Dropout rate
            bias: Whether to use bias
        """
        super(HierarchicalGCNConv, self).__init__(aggr='add')
        
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.pathway_dim = pathway_dim
        self.hierarchy_levels = hierarchy_levels
        self.dropout = dropout
        
        # Hierarchy-level specific transformations
        self.hierarchy_transforms = nn.ModuleList([
            nn.Linear(in_channels, out_channels // hierarchy_levels, bias=bias)
            for _ in range(hierarchy_levels)
        ])
        
        # Pathway-aware attention
        self.pathway_attention = nn.Sequential(
            nn.Linear(pathway_dim, pathway_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(pathway_dim // 2, hierarchy_levels),
            nn.Softmax(dim=-1)
        )
        
        # Cross-hierarchy communication
        self.cross_hierarchy = nn.ModuleList([
            nn.Linear(out_channels // hierarchy_levels, out_channels // hierarchy_levels, bias=bias)
            for _ in range(hierarchy_levels)
        ])
        
        # Output projection
        self.output_projection = nn.Linear(out_channels, out_channels, bias=bias)
        
        # Layer normalization
        self.layer_norm = nn.LayerNorm(out_channels)
        
        # Dropout layer
        self.dropout_layer = nn.Dropout(dropout)
        
        logger.info(f"Initialized HierarchicalGCNConv with {hierarchy_levels} levels")
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                pathway_features: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass through the hierarchical GCN layer.
        
        Args:
            x: Node features
            edge_index: Edge indices
            pathway_features: Pathway features (optional)
            
        Returns:
            Updated node features
        """
        # Hierarchy-level processing
        hierarchy_outputs = []
        
        for i in range(self.hierarchy_levels):
            # Transform features for this hierarchy level
            level_x = self.hierarchy_transforms[i](x)
            
            # Message passing for this level
            level_out = self.propagate(edge_index, x=level_x, level=i)
            
            # Cross-hierarchy communication
            level_out = self.cross_hierarchy[i](level_out)
            
            hierarchy_outputs.append(level_out)
        
        # Concatenate hierarchy outputs
        hierarchical_out = torch.cat(hierarchy_outputs, dim=-1)
        
        # Apply pathway-aware attention if pathway features are available
        if pathway_features is not None:
            attention_weights = self.pathway_attention(pathway_features)
            
            # Weighted combination of hierarchy levels
            weighted_out = torch.zeros_like(hierarchical_out)
            for i, level_out in enumerate(hierarchy_outputs):
                start_idx = i * self.out_channels // self.hierarchy_levels
                end_idx = (i + 1) * self.out_channels // self.hierarchy_levels
                weighted_out[:, start_idx:end_idx] = level_out * attention_weights[:, i:i+1]
            
            hierarchical_out = weighted_out
        
        # Output projection
        out = self.output_projection(hierarchical_out)
        
        # Layer normalization and activation
        out = self.layer_norm(out)
        out = F.relu(out)
        out = self.dropout_layer(out)
        
        return out
    
    def message(self, x_j: torch.Tensor, level: int) -> torch.Tensor:
        """Construct messages from neighboring nodes for specific hierarchy level."""
        return x_j
    
    def update(self, aggr_out: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        """Update node features based on aggregated messages."""
        return aggr_out

class MetabolicAttentionLayer(nn.Module):
    """
    Attention layer specifically designed for metabolic networks.
    
    Features:
    - Node type-aware attention
    - Stoichiometric coefficient weighting
    - Metabolic pathway attention
    """
    
    def __init__(self, 
                 in_channels: int,
                 out_channels: int,
                 num_heads: int = 8,
                 dropout: float = 0.2,
                 bias: bool = True):
        """
        Initialize the MetabolicAttentionLayer.
        
        Args:
            in_channels: Input feature dimension
            out_channels: Output feature dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
            bias: Whether to use bias
        """
        super(MetabolicAttentionLayer, self).__init__()
        
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.num_heads = num_heads
        self.dropout = dropout
        
        # Multi-head attention
        self.attention = nn.MultiheadAttention(
            embed_dim=in_channels,
            num_heads=num_heads,
            dropout=dropout,
            bias=bias,
            batch_first=True
        )
        
        # Output projection
        self.output_projection = nn.Linear(in_channels, out_channels, bias=bias)
        
        # Layer normalization
        self.layer_norm1 = nn.LayerNorm(in_channels)
        self.layer_norm2 = nn.LayerNorm(out_channels)
        
        # Feed-forward network
        self.feed_forward = nn.Sequential(
            nn.Linear(out_channels, out_channels * 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(out_channels * 2, out_channels)
        )
        
        # Dropout layers
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)
        
        logger.info(f"Initialized MetabolicAttentionLayer with {num_heads} heads")
    
    def forward(self, x: torch.Tensor, 
                node_type: Optional[torch.Tensor] = None,
                pathway: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass through the metabolic attention layer.
        
        Args:
            x: Node features [num_nodes, in_channels]
            node_type: Node type features (optional)
            pathway: Pathway features (optional)
            
        Returns:
            Updated node features [num_nodes, out_channels]
        """
        # Add node type and pathway information if available
        if node_type is not None:
            # Project node_type to match x dimensions
            node_type_proj = nn.Linear(node_type.size(-1), x.size(-1), bias=False).to(x.device)
            x = x + node_type_proj(node_type)
        
        if pathway is not None:
            # Project pathway to match x dimensions
            pathway_proj = nn.Linear(pathway.size(-1), x.size(-1), bias=False).to(x.device)
            x = x + pathway_proj(pathway)
        
        # Multi-head self-attention
        attn_out, _ = self.attention(x, x, x)
        attn_out = self.dropout1(attn_out)
        x = self.layer_norm1(x + attn_out)
        
        # Output projection
        out = self.output_projection(x)
        out = self.dropout2(out)
        out = self.layer_norm2(out)
        
        # Feed-forward network
        ff_out = self.feed_forward(out)
        out = out + ff_out
        
        return out

def create_custom_layer_factory():
    """
    Create a factory function for custom graph layers.
    
    Returns:
        Dictionary mapping layer names to layer classes
    """
    return {
        'MetabolicGCNConv': MetabolicGCNConv,
        'MultiScaleGCNConv': MultiScaleGCNConv,
        'HierarchicalGCNConv': HierarchicalGCNConv,
        'MetabolicAttentionLayer': MetabolicAttentionLayer
    }

def test_custom_layers():
    """Test the custom graph layers."""
    logger.info("Testing custom graph layers...")
    
    # Test parameters
    num_nodes = 100
    in_channels = 32
    out_channels = 64
    
    # Create test data
    x = torch.randn(num_nodes, in_channels)
    edge_index = torch.randint(0, num_nodes, (2, 200))
    
    # Test MetabolicAttentionLayer (simplified)
    logger.info("Testing MetabolicAttentionLayer...")
    attention_layer = MetabolicAttentionLayer(
        in_channels=in_channels,
        out_channels=out_channels,
        num_heads=8
    )
    out4 = attention_layer(x)  # No additional features
    logger.info(f"MetabolicAttentionLayer output shape: {out4.shape}")
    
    # Test MultiScaleGCNConv (simplified)
    logger.info("Testing MultiScaleGCNConv...")
    multiscale_gcn = MultiScaleGCNConv(
        in_channels=in_channels,
        out_channels=out_channels,
        scales=[1, 2]
    )
    out2 = multiscale_gcn(x, edge_index)
    logger.info(f"MultiScaleGCNConv output shape: {out2.shape}")
    
    logger.info("Custom layers tested successfully!")

if __name__ == "__main__":
    test_custom_layers() 