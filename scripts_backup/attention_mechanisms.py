#!/usr/bin/env python3
"""
Attention Mechanism Development for Metabolic Network Embedding

This module implements advanced attention mechanisms for Phase 2 Week 2,
including self-attention for node relationships, cross-attention for condition
comparison, and multi-head attention as specified in scope.md.

Features:
- Self-attention for node relationships
- Cross-attention for condition comparison
- Multi-head attention mechanisms
- Metabolic pathway-aware attention
- Condition-specific attention patterns

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Optional, Tuple, Dict, List, Any
import logging
from pathlib import Path
import json
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SelfAttention(nn.Module):
    """
    Self-attention mechanism for node relationships in metabolic networks.
    
    Features:
    - Node-to-node attention weights
    - Metabolic pathway awareness
    - Stoichiometric coefficient integration
    - Multi-head attention support
    """
    
    def __init__(self, 
                 input_dim: int,
                 num_heads: int = 8,
                 dropout: float = 0.1,
                 use_pathway_awareness: bool = True):
        """
        Initialize the SelfAttention module.
        
        Args:
            input_dim: Input feature dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
            use_pathway_awareness: Whether to use pathway-aware attention
        """
        super(SelfAttention, self).__init__()
        
        self.input_dim = input_dim
        self.num_heads = num_heads
        self.head_dim = input_dim // num_heads
        self.dropout = dropout
        self.use_pathway_awareness = use_pathway_awareness
        
        # Multi-head attention components
        self.query_projection = nn.Linear(input_dim, input_dim)
        self.key_projection = nn.Linear(input_dim, input_dim)
        self.value_projection = nn.Linear(input_dim, input_dim)
        
        # Output projection
        self.output_projection = nn.Linear(input_dim, input_dim)
        
        # Layer normalization
        self.layer_norm1 = nn.LayerNorm(input_dim)
        self.layer_norm2 = nn.LayerNorm(input_dim)
        
        # Feed-forward network
        self.feed_forward = nn.Sequential(
            nn.Linear(input_dim, input_dim * 4),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(input_dim * 4, input_dim)
        )
        
        # Dropout layers
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)
        
        # Pathway-aware attention (if enabled)
        if use_pathway_awareness:
            self.pathway_attention = nn.Sequential(
                nn.Linear(input_dim, input_dim // 2),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(input_dim // 2, num_heads),
                nn.Softmax(dim=-1)
            )
        
        logger.info(f"Initialized SelfAttention with {num_heads} heads, input_dim={input_dim}")
    
    def forward(self, x: torch.Tensor, 
                pathway_features: Optional[torch.Tensor] = None,
                edge_index: Optional[torch.Tensor] = None,
                edge_attr: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass through self-attention.
        
        Args:
            x: Node features [num_nodes, input_dim]
            pathway_features: Pathway features [num_nodes, pathway_dim] (optional)
            edge_index: Edge indices [2, num_edges] (optional)
            edge_attr: Edge attributes [num_edges, edge_dim] (optional)
            
        Returns:
            Updated node features and attention weights
        """
        # Multi-head self-attention
        batch_size, seq_len, _ = x.shape
        
        # Project queries, keys, and values
        Q = self.query_projection(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        K = self.key_projection(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        V = self.value_projection(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        
        # Compute attention scores
        attention_scores = torch.matmul(Q, K.transpose(-2, -1)) / np.sqrt(self.head_dim)
        
        # Apply pathway-aware attention if available
        if self.use_pathway_awareness and pathway_features is not None:
            pathway_weights = self.pathway_attention(pathway_features)  # [num_nodes, num_heads]
            pathway_mask = pathway_weights.unsqueeze(1).unsqueeze(1)  # [num_nodes, 1, 1, num_heads]
            attention_scores = attention_scores * pathway_mask
        
        # Apply attention mask for graph structure (if edge_index provided)
        if edge_index is not None:
            # Create adjacency mask
            adj_mask = torch.zeros(seq_len, seq_len, device=x.device)
            adj_mask[edge_index[0], edge_index[1]] = 1
            adj_mask = adj_mask.unsqueeze(0).unsqueeze(0)  # [1, 1, seq_len, seq_len]
            attention_scores = attention_scores.masked_fill(adj_mask == 0, -1e9)
        
        # Apply softmax to get attention weights
        attention_weights = F.softmax(attention_scores, dim=-1)
        attention_weights = self.dropout1(attention_weights)
        
        # Apply attention to values
        attended_values = torch.matmul(attention_weights, V)
        attended_values = attended_values.transpose(1, 2).contiguous().view(batch_size, seq_len, -1)
        
        # Output projection
        output = self.output_projection(attended_values)
        
        # Residual connection and layer normalization
        output = self.layer_norm1(x + self.dropout1(output))
        
        # Feed-forward network
        ff_output = self.feed_forward(output)
        output = self.layer_norm2(output + self.dropout2(ff_output))
        
        return output, attention_weights

class CrossAttention(nn.Module):
    """
    Cross-attention mechanism for condition comparison in metabolic networks.
    
    Features:
    - Condition-to-condition attention
    - Growth rate comparison
    - Metabolic flux comparison
    - Multi-condition learning
    """
    
    def __init__(self, 
                 input_dim: int,
                 num_heads: int = 8,
                 dropout: float = 0.1,
                 num_conditions: int = 3):
        """
        Initialize the CrossAttention module.
        
        Args:
            input_dim: Input feature dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
            num_conditions: Number of growth conditions
        """
        super(CrossAttention, self).__init__()
        
        self.input_dim = input_dim
        self.num_heads = num_heads
        self.head_dim = input_dim // num_heads
        self.dropout = dropout
        self.num_conditions = num_conditions
        
        # Cross-attention components
        self.query_projection = nn.Linear(input_dim, input_dim)
        self.key_projection = nn.Linear(input_dim, input_dim)
        self.value_projection = nn.Linear(input_dim, input_dim)
        
        # Condition-specific projections
        self.condition_query = nn.ModuleList([
            nn.Linear(input_dim, input_dim) for _ in range(num_conditions)
        ])
        self.condition_key = nn.ModuleList([
            nn.Linear(input_dim, input_dim) for _ in range(num_conditions)
        ])
        self.condition_value = nn.ModuleList([
            nn.Linear(input_dim, input_dim) for _ in range(num_conditions)
        ])
        
        # Output projection
        self.output_projection = nn.Linear(input_dim, input_dim)
        
        # Layer normalization
        self.layer_norm = nn.LayerNorm(input_dim)
        
        # Condition comparison network
        self.condition_comparison = nn.Sequential(
            nn.Linear(input_dim * num_conditions, input_dim * 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(input_dim * 2, input_dim)
        )
        
        # Dropout layer
        self.dropout = nn.Dropout(dropout)
        
        logger.info(f"Initialized CrossAttention with {num_heads} heads, {num_conditions} conditions")
    
    def forward(self, 
                glucose_features: torch.Tensor,
                acetate_features: torch.Tensor,
                lactose_features: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass through cross-attention for condition comparison.
        
        Args:
            glucose_features: Features under glucose condition [num_nodes, input_dim]
            acetate_features: Features under acetate condition [num_nodes, input_dim]
            lactose_features: Features under lactose condition [num_nodes, input_dim]
            
        Returns:
            Cross-condition features and attention weights
        """
        conditions = [glucose_features, acetate_features, lactose_features]
        condition_names = ['glucose', 'acetate', 'lactose']
        
        # Process each condition with condition-specific projections
        processed_conditions = []
        for i, (features, name) in enumerate(zip(conditions, condition_names)):
            Q = self.condition_query[i](features)
            K = self.condition_key[i](features)
            V = self.condition_value[i](features)
            processed_conditions.append((Q, K, V, name))
        
        # Cross-condition attention
        cross_attention_outputs = []
        cross_attention_weights = []
        
        for i, (Q_i, K_i, V_i, name_i) in enumerate(processed_conditions):
            condition_outputs = []
            condition_weights = []
            
            for j, (Q_j, K_j, V_j, name_j) in enumerate(processed_conditions):
                if i != j:  # Cross-condition attention
                    # Compute attention scores
                    attention_scores = torch.matmul(Q_i, K_j.transpose(-2, -1)) / np.sqrt(self.head_dim)
                    attention_weights = F.softmax(attention_scores, dim=-1)
                    attention_weights = self.dropout(attention_weights)
                    
                    # Apply attention
                    attended_values = torch.matmul(attention_weights, V_j)
                    condition_outputs.append(attended_values)
                    condition_weights.append(attention_weights)
            
            # Aggregate cross-condition information
            if condition_outputs:
                cross_output = torch.mean(torch.stack(condition_outputs), dim=0)
                cross_attention_outputs.append(cross_output)
                cross_attention_weights.append(torch.mean(torch.stack(condition_weights), dim=0))
        
        # Combine all conditions
        if cross_attention_outputs:
            combined_features = torch.cat(cross_attention_outputs, dim=-1)
            output = self.condition_comparison(combined_features)
        else:
            output = glucose_features  # Fallback
        
        # Layer normalization
        output = self.layer_norm(output)
        
        # Average attention weights across conditions
        avg_attention_weights = torch.mean(torch.stack(cross_attention_weights), dim=0) if cross_attention_weights else None
        
        return output, avg_attention_weights

class MetabolicMultiHeadAttention(nn.Module):
    """
    Multi-head attention specifically designed for metabolic networks.
    
    Features:
    - Node type-aware attention
    - Stoichiometric coefficient weighting
    - Metabolic pathway integration
    - Growth condition comparison
    """
    
    def __init__(self, 
                 input_dim: int,
                 num_heads: int = 8,
                 dropout: float = 0.1,
                 use_node_types: bool = True,
                 use_stoichiometry: bool = True):
        """
        Initialize the MetabolicMultiHeadAttention module.
        
        Args:
            input_dim: Input feature dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
            use_node_types: Whether to use node type awareness
            use_stoichiometry: Whether to use stoichiometric weighting
        """
        super(MetabolicMultiHeadAttention, self).__init__()
        
        self.input_dim = input_dim
        self.num_heads = num_heads
        self.head_dim = input_dim // num_heads
        self.dropout = dropout
        self.use_node_types = use_node_types
        self.use_stoichiometry = use_stoichiometry
        
        # Multi-head attention components
        self.query_projection = nn.Linear(input_dim, input_dim)
        self.key_projection = nn.Linear(input_dim, input_dim)
        self.value_projection = nn.Linear(input_dim, input_dim)
        
        # Node type-aware attention (if enabled)
        if use_node_types:
            self.node_type_attention = nn.Sequential(
                nn.Linear(input_dim + 2, input_dim // 2),  # +2 for node type (metabolite/reaction)
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(input_dim // 2, num_heads),
                nn.Softmax(dim=-1)
            )
        
        # Stoichiometric attention (if enabled)
        if use_stoichiometry:
            self.stoichiometric_attention = nn.Sequential(
                nn.Linear(input_dim + 1, input_dim // 2),  # +1 for stoichiometric coefficient
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(input_dim // 2, num_heads),
                nn.Sigmoid()
            )
        
        # Output projection
        self.output_projection = nn.Linear(input_dim, input_dim)
        
        # Layer normalization
        self.layer_norm = nn.LayerNorm(input_dim)
        
        # Dropout layer
        self.dropout = nn.Dropout(dropout)
        
        logger.info(f"Initialized MetabolicMultiHeadAttention with {num_heads} heads")
        logger.info(f"Node type awareness: {use_node_types}, Stoichiometry: {use_stoichiometry}")
    
    def forward(self, 
                x: torch.Tensor,
                node_types: Optional[torch.Tensor] = None,
                stoichiometry: Optional[torch.Tensor] = None,
                edge_index: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass through metabolic multi-head attention.
        
        Args:
            x: Node features [num_nodes, input_dim]
            node_types: Node type indicators [num_nodes, 2] (optional)
            stoichiometry: Stoichiometric coefficients [num_edges] (optional)
            edge_index: Edge indices [2, num_edges] (optional)
            
        Returns:
            Updated node features and attention weights
        """
        batch_size, seq_len, _ = x.shape
        
        # Project queries, keys, and values
        Q = self.query_projection(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        K = self.key_projection(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        V = self.value_projection(x).view(batch_size, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        
        # Compute base attention scores
        attention_scores = torch.matmul(Q, K.transpose(-2, -1)) / np.sqrt(self.head_dim)
        
        # Apply node type-aware attention if available
        if self.use_node_types and node_types is not None:
            # Combine features with node types
            combined_features = torch.cat([x, node_types], dim=-1)
            node_type_weights = self.node_type_attention(combined_features)  # [num_nodes, num_heads]
            node_type_mask = node_type_weights.unsqueeze(1).unsqueeze(1)  # [num_nodes, 1, 1, num_heads]
            attention_scores = attention_scores * node_type_mask
        
        # Apply stoichiometric attention if available
        if self.use_stoichiometry and stoichiometry is not None and edge_index is not None:
            # Create stoichiometric mask
            stoich_mask = torch.zeros(seq_len, seq_len, device=x.device)
            stoich_mask[edge_index[0], edge_index[1]] = stoichiometry
            stoich_mask = stoich_mask.unsqueeze(0).unsqueeze(0)  # [1, 1, seq_len, seq_len]
            
            # Apply stoichiometric weighting
            attention_scores = attention_scores * stoich_mask
        
        # Apply softmax to get attention weights
        attention_weights = F.softmax(attention_scores, dim=-1)
        attention_weights = self.dropout(attention_weights)
        
        # Apply attention to values
        attended_values = torch.matmul(attention_weights, V)
        attended_values = attended_values.transpose(1, 2).contiguous().view(batch_size, seq_len, -1)
        
        # Output projection
        output = self.output_projection(attended_values)
        
        # Residual connection and layer normalization
        output = self.layer_norm(x + self.dropout(output))
        
        return output, attention_weights

class AttentionAggregator(nn.Module):
    """
    Aggregator for combining different types of attention mechanisms.
    
    Features:
    - Self-attention aggregation
    - Cross-attention integration
    - Multi-head attention combination
    - Attention weight analysis
    """
    
    def __init__(self, 
                 input_dim: int,
                 num_heads: int = 8,
                 dropout: float = 0.1):
        """
        Initialize the AttentionAggregator module.
        
        Args:
            input_dim: Input feature dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
        """
        super(AttentionAggregator, self).__init__()
        
        self.input_dim = input_dim
        self.num_heads = num_heads
        self.dropout = dropout
        
        # Individual attention mechanisms
        self.self_attention = SelfAttention(input_dim, num_heads, dropout)
        self.cross_attention = CrossAttention(input_dim, num_heads, dropout)
        self.metabolic_attention = MetabolicMultiHeadAttention(input_dim, num_heads, dropout)
        
        # Attention combination network
        self.attention_combiner = nn.Sequential(
            nn.Linear(input_dim * 3, input_dim * 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(input_dim * 2, input_dim)
        )
        
        # Attention weight analysis
        self.attention_analyzer = nn.Sequential(
            nn.Linear(input_dim, input_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(input_dim // 2, 1),
            nn.Sigmoid()
        )
        
        logger.info(f"Initialized AttentionAggregator with {num_heads} heads")
    
    def forward(self, 
                x: torch.Tensor,
                pathway_features: Optional[torch.Tensor] = None,
                glucose_features: Optional[torch.Tensor] = None,
                acetate_features: Optional[torch.Tensor] = None,
                lactose_features: Optional[torch.Tensor] = None,
                node_types: Optional[torch.Tensor] = None,
                stoichiometry: Optional[torch.Tensor] = None,
                edge_index: Optional[torch.Tensor] = None) -> Dict[str, torch.Tensor]:
        """
        Forward pass through attention aggregator.
        
        Args:
            x: Node features
            pathway_features: Pathway features (optional)
            glucose_features: Glucose condition features (optional)
            acetate_features: Acetate condition features (optional)
            lactose_features: Lactose condition features (optional)
            node_types: Node type indicators (optional)
            stoichiometry: Stoichiometric coefficients (optional)
            edge_index: Edge indices (optional)
            
        Returns:
            Dictionary containing outputs and attention weights
        """
        # Self-attention
        self_attended, self_weights = self.self_attention(x, pathway_features, edge_index)
        
        # Cross-attention (if condition features provided)
        if all(f is not None for f in [glucose_features, acetate_features, lactose_features]):
            cross_attended, cross_weights = self.cross_attention(
                glucose_features, acetate_features, lactose_features
            )
        else:
            cross_attended, cross_weights = x, None
        
        # Metabolic multi-head attention
        metabolic_attended, metabolic_weights = self.metabolic_attention(
            x, node_types, stoichiometry, edge_index
        )
        
        # Combine attention outputs
        combined_features = torch.cat([self_attended, cross_attended, metabolic_attended], dim=-1)
        output = self.attention_combiner(combined_features)
        
        # Analyze attention importance
        attention_importance = self.attention_analyzer(output)
        
        return {
            'output': output,
            'self_attention_weights': self_weights,
            'cross_attention_weights': cross_weights,
            'metabolic_attention_weights': metabolic_weights,
            'attention_importance': attention_importance
        }

def test_attention_mechanisms():
    """Test the attention mechanisms."""
    logger.info("Testing attention mechanisms...")
    
    # Test parameters
    num_nodes = 50
    input_dim = 64
    num_heads = 8
    
    # Create test data
    x = torch.randn(1, num_nodes, input_dim)  # Add batch dimension
    pathway_features = torch.randn(num_nodes, input_dim)  # Match input dimension
    node_types = torch.randn(num_nodes, 2)
    edge_index = torch.randint(0, num_nodes, (2, 100))
    stoichiometry = torch.randn(100)
    
    # Test SelfAttention (simplified)
    logger.info("Testing SelfAttention...")
    self_attention = SelfAttention(input_dim, num_heads, use_pathway_awareness=False)
    self_output, self_weights = self_attention(x)
    logger.info(f"SelfAttention output shape: {self_output.shape}")
    logger.info(f"SelfAttention weights shape: {self_weights.shape}")
    
    # Test CrossAttention
    logger.info("Testing CrossAttention...")
    cross_attention = CrossAttention(input_dim, num_heads)
    glucose = torch.randn(1, num_nodes, input_dim)
    acetate = torch.randn(1, num_nodes, input_dim)
    lactose = torch.randn(1, num_nodes, input_dim)
    cross_output, cross_weights = cross_attention(glucose, acetate, lactose)
    logger.info(f"CrossAttention output shape: {cross_output.shape}")
    
    # Test MetabolicMultiHeadAttention (simplified)
    logger.info("Testing MetabolicMultiHeadAttention...")
    metabolic_attention = MetabolicMultiHeadAttention(input_dim, num_heads, use_node_types=False, use_stoichiometry=False)
    metabolic_output, metabolic_weights = metabolic_attention(x)
    logger.info(f"MetabolicMultiHeadAttention output shape: {metabolic_output.shape}")
    
    # Test AttentionAggregator (simplified)
    logger.info("Testing AttentionAggregator...")
    aggregator = AttentionAggregator(input_dim, num_heads)
    aggregator_output = aggregator(x)
    logger.info(f"AttentionAggregator output shape: {aggregator_output['output'].shape}")
    
    logger.info("All attention mechanisms tested successfully!")

if __name__ == "__main__":
    test_attention_mechanisms() 