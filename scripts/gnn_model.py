#!/usr/bin/env python3
"""
Graph Neural Network Implementation for Metabolic Network Embedding

This module implements the core GNN architecture for Phase 2 Week 1,
including GCN layers, message passing framework, hierarchical feature learning,
and multi-layer architecture as specified in scope.md.

Features:
- GCN layers using PyTorch Geometric
- Message passing framework
- Hierarchical feature learning
- Multi-layer architecture
- Attention mechanisms
- Multi-task learning support

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv, GraphConv, global_mean_pool, global_max_pool
from torch_geometric.data import Data, DataLoader
import numpy as np
import json
from pathlib import Path
import logging
from typing import Dict, List, Optional, Tuple, Any
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MetabolicGNN(nn.Module):
    """
    Core GNN implementation for metabolic network embedding.
    
    Architecture:
    - Input layer: Node feature processing
    - GCN layers: Graph convolution with residual connections
    - Attention layer: Self-attention for node relationships
    - Output layer: Multi-task prediction
    """
    
    def __init__(self, 
                 input_dim: int,
                 hidden_dim: int = 128,
                 output_dim: int = 64,
                 num_layers: int = 3,
                 dropout: float = 0.2,
                 use_attention: bool = True,
                 use_residual: bool = True):
        """
        Initialize the MetabolicGNN model.
        
        Args:
            input_dim: Dimension of input node features
            hidden_dim: Dimension of hidden layers
            output_dim: Dimension of output embeddings
            num_layers: Number of GCN layers
            dropout: Dropout rate
            use_attention: Whether to use attention mechanism
            use_residual: Whether to use residual connections
        """
        super(MetabolicGNN, self).__init__()
        
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.num_layers = num_layers
        self.dropout = dropout
        self.use_attention = use_attention
        self.use_residual = use_residual
        
        # Input projection layer
        self.input_projection = nn.Linear(input_dim, hidden_dim)
        
        # GCN layers
        self.gcn_layers = nn.ModuleList()
        for i in range(num_layers):
            if i == 0:
                self.gcn_layers.append(GCNConv(hidden_dim, hidden_dim))
            else:
                self.gcn_layers.append(GCNConv(hidden_dim, hidden_dim))
        
        # Attention layer (if enabled)
        if use_attention:
            self.attention = GATConv(hidden_dim, hidden_dim, heads=4, dropout=dropout)
            self.attention_out = nn.Linear(hidden_dim * 4, hidden_dim)
        
        # Output projection
        self.output_projection = nn.Linear(hidden_dim, output_dim)
        
        # Layer normalization
        self.layer_norms = nn.ModuleList([
            nn.LayerNorm(hidden_dim) for _ in range(num_layers)
        ])
        
        # Dropout layers
        self.dropout_layers = nn.ModuleList([
            nn.Dropout(dropout) for _ in range(num_layers)
        ])
        
        logger.info(f"Initialized MetabolicGNN with {input_dim}->{hidden_dim}->{output_dim}")
        logger.info(f"Layers: {num_layers}, Attention: {use_attention}, Residual: {use_residual}")
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, 
                batch: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass through the GNN.
        
        Args:
            x: Node features [num_nodes, input_dim]
            edge_index: Edge indices [2, num_edges]
            batch: Batch assignment [num_nodes] (optional)
            
        Returns:
            Node embeddings [num_nodes, output_dim]
        """
        # Input projection
        x = self.input_projection(x)
        x = F.relu(x)
        x = self.dropout_layers[0](x)
        
        # GCN layers with residual connections
        for i, (gcn_layer, layer_norm, dropout_layer) in enumerate(
            zip(self.gcn_layers, self.layer_norms, self.dropout_layers[1:])
        ):
            # GCN convolution
            h = gcn_layer(x, edge_index)
            h = layer_norm(h)
            h = F.relu(h)
            h = dropout_layer(h)
            
            # Residual connection
            if self.use_residual and i > 0:
                x = x + h
            else:
                x = h
        
        # Attention mechanism (if enabled)
        if self.use_attention:
            x = self.attention(x, edge_index)
            x = self.attention_out(x)
            x = F.relu(x)
            x = self.dropout_layers[-1](x)
        
        # Output projection
        x = self.output_projection(x)
        
        return x
    
    def get_embeddings(self, x: torch.Tensor, edge_index: torch.Tensor,
                      batch: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Get node embeddings without final projection.
        
        Args:
            x: Node features
            edge_index: Edge indices
            batch: Batch assignment (optional)
            
        Returns:
            Node embeddings before final projection
        """
        # Input projection
        x = self.input_projection(x)
        x = F.relu(x)
        x = self.dropout_layers[0](x)
        
        # GCN layers
        for i, (gcn_layer, layer_norm, dropout_layer) in enumerate(
            zip(self.gcn_layers, self.layer_norms, self.dropout_layers[1:])
        ):
            h = gcn_layer(x, edge_index)
            h = layer_norm(h)
            h = F.relu(h)
            h = dropout_layer(h)
            
            if self.use_residual and i > 0:
                x = x + h
            else:
                x = h
        
        # Attention mechanism
        if self.use_attention:
            x = self.attention(x, edge_index)
            x = self.attention_out(x)
            x = F.relu(x)
            x = self.dropout_layers[-1](x)
        
        return x

class MultiTaskMetabolicGNN(nn.Module):
    """
    Multi-task GNN for metabolic network analysis.
    
    Supports multiple prediction tasks:
    - Node classification (metabolite vs reaction)
    - Node regression (flux prediction)
    - Graph-level prediction (growth rate)
    """
    
    def __init__(self, 
                 input_dim: int,
                 hidden_dim: int = 128,
                 embedding_dim: int = 64,
                 num_layers: int = 3,
                 dropout: float = 0.2,
                 num_node_classes: int = 2,
                 num_regression_targets: int = 1):
        """
        Initialize the MultiTaskMetabolicGNN model.
        
        Args:
            input_dim: Dimension of input node features
            hidden_dim: Dimension of hidden layers
            embedding_dim: Dimension of node embeddings
            num_layers: Number of GCN layers
            dropout: Dropout rate
            num_node_classes: Number of node classification classes
            num_regression_targets: Number of regression targets
        """
        super(MultiTaskMetabolicGNN, self).__init__()
        
        # Core GNN encoder
        self.gnn_encoder = MetabolicGNN(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            output_dim=embedding_dim,
            num_layers=num_layers,
            dropout=dropout,
            use_attention=True,
            use_residual=True
        )
        
        # Task-specific heads
        self.node_classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, num_node_classes)
        )
        
        self.node_regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, num_regression_targets)
        )
        
        # Graph-level prediction (global pooling)
        self.graph_predictor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
        
        logger.info(f"Initialized MultiTaskMetabolicGNN with {num_node_classes} classes and {num_regression_targets} regression targets")
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                batch: Optional[torch.Tensor] = None) -> Dict[str, torch.Tensor]:
        """
        Forward pass with multi-task predictions.
        
        Args:
            x: Node features
            edge_index: Edge indices
            batch: Batch assignment (optional)
            
        Returns:
            Dictionary of predictions for different tasks
        """
        # Get node embeddings
        node_embeddings = self.gnn_encoder.get_embeddings(x, edge_index, batch)
        
        # Node-level predictions
        node_class_logits = self.node_classifier(node_embeddings)
        node_regression = self.node_regressor(node_embeddings)
        
        # Graph-level predictions (if batch is provided)
        graph_predictions = None
        if batch is not None:
            # Global pooling
            graph_embeddings = global_mean_pool(node_embeddings, batch)
            graph_predictions = self.graph_predictor(graph_embeddings)
        
        return {
            'node_embeddings': node_embeddings,
            'node_class_logits': node_class_logits,
            'node_regression': node_regression,
            'graph_predictions': graph_predictions
        }

class GraphDataProcessor:
    """Process graph data for GNN training."""
    
    def __init__(self, feature_dir: str = "results/metabolic_network/feature_engineering"):
        """Initialize the data processor."""
        self.feature_dir = Path(feature_dir)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def load_processed_features(self) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """Load processed node and edge features."""
        logger.info("Loading processed features...")
        
        # Find latest feature files
        node_files = list(self.feature_dir.glob("processed_node_features_*.json"))
        edge_files = list(self.feature_dir.glob("processed_edge_features_*.json"))
        
        if not node_files or not edge_files:
            raise FileNotFoundError("Processed feature files not found")
        
        # Load latest files
        latest_node_file = max(node_files, key=lambda x: x.stat().st_mtime)
        latest_edge_file = max(edge_files, key=lambda x: x.stat().st_mtime)
        
        logger.info(f"Loading node features from: {latest_node_file}")
        logger.info(f"Loading edge features from: {latest_edge_file}")
        
        with open(latest_node_file, 'r') as f:
            node_features = json.load(f)
        
        with open(latest_edge_file, 'r') as f:
            edge_features = json.load(f)
        
        # Convert to numpy arrays
        node_features_processed = {}
        for node_id, data in node_features.items():
            node_features_processed[node_id] = {
                'features': np.array(data['features'], dtype=np.float32),
                'feature_names': data['feature_names']
            }
        
        edge_features_processed = {}
        for edge_id, data in edge_features.items():
            edge_features_processed[edge_id] = {
                'features': np.array(data['features'], dtype=np.float32),
                'feature_names': data['feature_names']
            }
        
        logger.info(f"Loaded {len(node_features_processed)} node features")
        logger.info(f"Loaded {len(edge_features_processed)} edge features")
        
        return node_features_processed, edge_features_processed
    
    def create_pytorch_geometric_data(self, 
                                    node_features: Dict[str, np.ndarray],
                                    edge_features: Dict[str, np.ndarray],
                                    graph_file: str = "results/metabolic_network/metabolic_network_graph.pkl") -> Data:
        """
        Create PyTorch Geometric Data object from features.
        
        Args:
            node_features: Processed node features
            edge_features: Processed edge features
            graph_file: Path to network graph file
            
        Returns:
            PyTorch Geometric Data object
        """
        logger.info("Creating PyTorch Geometric Data object...")
        
        # Load graph
        import pickle
        with open(graph_file, 'rb') as f:
            graph = pickle.load(f)
        
        # Create node feature matrix
        node_ids = list(node_features.keys())
        node_id_to_idx = {node_id: idx for idx, node_id in enumerate(node_ids)}
        
        # Get feature dimensions and find maximum
        feature_dims = [len(node_features[node_id]['features']) for node_id in node_ids]
        max_feature_dim = max(feature_dims)
        min_feature_dim = min(feature_dims)
        
        if max_feature_dim != min_feature_dim:
            logger.warning(f"Feature dimensions vary: min={min_feature_dim}, max={max_feature_dim}")
            logger.warning("Padding shorter feature vectors with zeros")
        
        node_feature_matrix = np.zeros((len(node_ids), max_feature_dim), dtype=np.float32)
        
        for idx, node_id in enumerate(node_ids):
            features = node_features[node_id]['features']
            if len(features) < max_feature_dim:
                # Pad with zeros
                padded_features = np.pad(features, (0, max_feature_dim - len(features)), 'constant')
                node_feature_matrix[idx] = padded_features
            else:
                node_feature_matrix[idx] = features
        
        # Create edge index matrix
        edge_list = []
        edge_attr_list = []
        
        for edge in graph.edges():
            source, target = edge
            if source in node_id_to_idx and target in node_id_to_idx:
                source_idx = node_id_to_idx[source]
                target_idx = node_id_to_idx[target]
                
                edge_list.append([source_idx, target_idx])
                
                # Get edge features if available
                edge_key = f"{source}_{target}"
                if edge_key in edge_features:
                    edge_attr_list.append(edge_features[edge_key]['features'])
                else:
                    # Default edge features
                    edge_attr_list.append(np.zeros(16, dtype=np.float32))
        
        edge_index = np.array(edge_list, dtype=np.int64).T
        edge_attr = np.array(edge_attr_list, dtype=np.float32)
        
        # Create PyTorch Geometric Data object
        data = Data(
            x=torch.tensor(node_feature_matrix, dtype=torch.float32),
            edge_index=torch.tensor(edge_index, dtype=torch.long),
            edge_attr=torch.tensor(edge_attr, dtype=torch.float32)
        )
        
        logger.info(f"Created Data object with {data.x.shape[0]} nodes and {data.edge_index.shape[1]} edges")
        logger.info(f"Node features: {data.x.shape[1]}, Edge features: {data.edge_attr.shape[1]}")
        
        return data

class GNNTrainer:
    """Trainer for GNN models."""
    
    def __init__(self, model: nn.Module, device: str = 'cpu'):
        """Initialize the trainer."""
        self.model = model.to(device)
        self.device = device
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def train_step(self, data: Data, optimizer: torch.optim.Optimizer) -> Dict[str, float]:
        """Single training step."""
        self.model.train()
        optimizer.zero_grad()
        
        # Move data to device
        data = data.to(self.device)
        
        # Forward pass
        outputs = self.model(data.x, data.edge_index, data.batch)
        
        # Calculate losses (placeholder - customize based on tasks)
        loss = torch.tensor(0.0, device=self.device)
        
        # Backward pass
        loss.backward()
        optimizer.step()
        
        return {'loss': loss.item()}
    
    def validate_step(self, data: Data) -> Dict[str, float]:
        """Single validation step."""
        self.model.eval()
        
        with torch.no_grad():
            data = data.to(self.device)
            outputs = self.model(data.x, data.edge_index, data.batch)
            
            # Calculate metrics (placeholder)
            loss = torch.tensor(0.0, device=self.device)
            
        return {'val_loss': loss.item()}

def create_model_architecture_visualization(model: nn.Module, 
                                          output_path: str = "results/metabolic_network/model_architecture.png"):
    """Create visualization of the model architecture."""
    logger.info("Creating model architecture visualization...")
    
    # Create a simple architecture diagram
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define layer positions
    layers = [
        ('Input Features', 0),
        ('Input Projection', 1),
        ('GCN Layer 1', 2),
        ('GCN Layer 2', 3),
        ('GCN Layer 3', 4),
        ('Attention Layer', 5),
        ('Output Projection', 6),
        ('Node Embeddings', 7)
    ]
    
    # Draw layers
    for layer_name, y_pos in layers:
        ax.add_patch(plt.Rectangle((0.1, y_pos), 0.8, 0.6, 
                                 facecolor='lightblue', edgecolor='black'))
        ax.text(0.5, y_pos + 0.3, layer_name, ha='center', va='center', fontsize=10)
    
    # Draw connections
    for i in range(len(layers) - 1):
        ax.arrow(0.5, layers[i][1] + 0.6, 0, 0.4, head_width=0.02, head_length=0.1, fc='black')
    
    # Add annotations
    ax.text(0.5, -0.5, 'Metabolic GNN Architecture', ha='center', va='center', fontsize=14, fontweight='bold')
    ax.text(0.5, -0.8, f'Input: {model.input_dim} → Hidden: {model.hidden_dim} → Output: {model.output_dim}', 
            ha='center', va='center', fontsize=10)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Architecture visualization saved to: {output_path}")

def main():
    """Main function to demonstrate GNN implementation."""
    logger.info("=" * 60)
    logger.info("PHASE 2 WEEK 1: GRAPH NEURAL NETWORK IMPLEMENTATION")
    logger.info("=" * 60)
    
    # Step 1: Load processed features
    logger.info("\nStep 1: Loading processed features...")
    processor = GraphDataProcessor()
    node_features, edge_features = processor.load_processed_features()
    
    # Step 2: Create PyTorch Geometric data
    logger.info("\nStep 2: Creating PyTorch Geometric data...")
    data = processor.create_pytorch_geometric_data(node_features, edge_features)
    
    # Step 3: Initialize GNN model
    logger.info("\nStep 3: Initializing GNN model...")
    input_dim = data.x.shape[1]
    model = MetabolicGNN(
        input_dim=input_dim,
        hidden_dim=128,
        output_dim=64,
        num_layers=3,
        dropout=0.2,
        use_attention=True,
        use_residual=True
    )
    
    # Step 4: Create multi-task model
    logger.info("\nStep 4: Creating multi-task model...")
    multi_task_model = MultiTaskMetabolicGNN(
        input_dim=input_dim,
        hidden_dim=128,
        embedding_dim=64,
        num_layers=3,
        dropout=0.2,
        num_node_classes=2,
        num_regression_targets=1
    )
    
    # Step 5: Test forward pass
    logger.info("\nStep 5: Testing forward pass...")
    with torch.no_grad():
        # Test basic GNN
        embeddings = model(data.x, data.edge_index)
        logger.info(f"Basic GNN output shape: {embeddings.shape}")
        
        # Test multi-task GNN
        outputs = multi_task_model(data.x, data.edge_index)
        logger.info(f"Multi-task outputs:")
        for key, value in outputs.items():
            if value is not None:
                logger.info(f"  {key}: {value.shape}")
    
    # Step 6: Create architecture visualization
    logger.info("\nStep 6: Creating architecture visualization...")
    create_model_architecture_visualization(model)
    
    # Step 7: Save model configuration
    logger.info("\nStep 7: Saving model configuration...")
    config = {
        'model_type': 'MetabolicGNN',
        'input_dim': input_dim,
        'hidden_dim': 128,
        'output_dim': 64,
        'num_layers': 3,
        'dropout': 0.2,
        'use_attention': True,
        'use_residual': True,
        'num_nodes': data.x.shape[0],
        'num_edges': data.edge_index.shape[1],
        'node_features': data.x.shape[1],
        'edge_features': data.edge_attr.shape[1] if hasattr(data, 'edge_attr') else 0
    }
    
    config_path = Path("results/metabolic_network/gnn_config.json")
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    logger.info(f"Model configuration saved to: {config_path}")
    
    # Step 8: Summary
    logger.info("\n" + "=" * 60)
    logger.info("GNN IMPLEMENTATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Input dimension: {input_dim}")
    logger.info(f"Hidden dimension: {128}")
    logger.info(f"Output dimension: {64}")
    logger.info(f"Number of layers: {3}")
    logger.info(f"Attention mechanism: Enabled")
    logger.info(f"Residual connections: Enabled")
    logger.info(f"Multi-task support: Implemented")
    logger.info(f"Data compatibility: PyTorch Geometric")
    
    return model, multi_task_model, data

if __name__ == "__main__":
    main() 