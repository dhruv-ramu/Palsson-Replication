#!/usr/bin/env python3
"""
Robust Data Handling for Metabolic Network Embedding

This module provides a complete rewrite of data handling with:
- Consistent data structures
- Proper validation
- Robust batching
- Type safety
- Error handling

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import numpy as np
from typing import Dict, List, Optional, Tuple, Any, Union
import logging
from pathlib import Path
import json
from datetime import datetime
from collections import defaultdict
import warnings
from dataclasses import dataclass
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import networkx as nx

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class MetabolicGraph:
    """
    Consistent data structure for metabolic network graphs.
    
    This ensures all components work with the same data format.
    """
    # Node features: [num_nodes, feature_dim]
    node_features: torch.Tensor
    
    # Edge indices: [2, num_edges]
    edge_indices: torch.Tensor
    
    # Edge features: [num_edges, edge_feature_dim] (optional)
    edge_features: Optional[torch.Tensor] = None
    
    # Node types: [num_nodes] (optional)
    node_types: Optional[torch.Tensor] = None
    
    # Graph-level features: [graph_feature_dim] (optional)
    graph_features: Optional[torch.Tensor] = None
    
    # Metadata
    num_nodes: int = 0
    num_edges: int = 0
    feature_dim: int = 0
    edge_feature_dim: int = 0
    
    def __post_init__(self):
        """Validate and set metadata after initialization."""
        if self.node_features.dim() != 2:
            raise ValueError(f"node_features must be 2D, got {self.node_features.dim()}D")
        
        if self.edge_indices.dim() != 2 or self.edge_indices.size(0) != 2:
            raise ValueError(f"edge_indices must be [2, num_edges], got {self.edge_indices.shape}")
        
        # Set metadata
        self.num_nodes = self.node_features.size(0)
        self.num_edges = self.edge_indices.size(1)
        self.feature_dim = self.node_features.size(1)
        
        # Validate edge indices
        if self.edge_indices.max() >= self.num_nodes:
            raise ValueError(f"Edge indices exceed node count: max={self.edge_indices.max()}, nodes={self.num_nodes}")
        
        # Validate edge features if provided
        if self.edge_features is not None:
            if self.edge_features.size(0) != self.num_edges:
                raise ValueError(f"Edge features count mismatch: {self.edge_features.size(0)} vs {self.num_edges}")
            self.edge_feature_dim = self.edge_features.size(1)
        
        # Validate node types if provided
        if self.node_types is not None:
            if self.node_types.size(0) != self.num_nodes:
                raise ValueError(f"Node types count mismatch: {self.node_types.size(0)} vs {self.num_nodes}")
    
    def to_device(self, device: torch.device) -> 'MetabolicGraph':
        """Move all tensors to specified device."""
        return MetabolicGraph(
            node_features=self.node_features.to(device),
            edge_indices=self.edge_indices.to(device),
            edge_features=self.edge_features.to(device) if self.edge_features is not None else None,
            node_types=self.node_types.to(device) if self.node_types is not None else None,
            graph_features=self.graph_features.to(device) if self.graph_features is not None else None,
            num_nodes=self.num_nodes,
            num_edges=self.num_edges,
            feature_dim=self.feature_dim,
            edge_feature_dim=self.edge_feature_dim
        )
    
    def clone(self) -> 'MetabolicGraph':
        """Create a deep copy of the graph."""
        return MetabolicGraph(
            node_features=self.node_features.clone(),
            edge_indices=self.edge_indices.clone(),
            edge_features=self.edge_features.clone() if self.edge_features is not None else None,
            node_types=self.node_types.clone() if self.node_types is not None else None,
            graph_features=self.graph_features.clone() if self.graph_features is not None else None,
            num_nodes=self.num_nodes,
            num_edges=self.num_edges,
            feature_dim=self.feature_dim,
            edge_feature_dim=self.edge_feature_dim
        )

@dataclass
class MetabolicTargets:
    """
    Consistent structure for target values.
    """
    # Node-level targets
    node_classification: Optional[torch.Tensor] = None  # [num_nodes, num_classes]
    node_regression: Optional[torch.Tensor] = None      # [num_nodes]
    
    # Graph-level targets
    graph_classification: Optional[torch.Tensor] = None  # [num_classes]
    graph_regression: Optional[torch.Tensor] = None      # [1]
    
    # Condition-specific targets
    condition_prediction: Optional[torch.Tensor] = None  # [num_conditions]
    pathway_analysis: Optional[torch.Tensor] = None      # [num_pathways]
    growth_rate: Optional[torch.Tensor] = None           # [1]
    
    def __post_init__(self):
        """Validate target structure."""
        targets = []
        if self.node_classification is not None:
            targets.append(('node_classification', self.node_classification))
        if self.node_regression is not None:
            targets.append(('node_regression', self.node_regression))
        if self.graph_classification is not None:
            targets.append(('graph_classification', self.graph_classification))
        if self.graph_regression is not None:
            targets.append(('graph_regression', self.graph_regression))
        if self.condition_prediction is not None:
            targets.append(('condition_prediction', self.condition_prediction))
        if self.pathway_analysis is not None:
            targets.append(('pathway_analysis', self.pathway_analysis))
        if self.growth_rate is not None:
            targets.append(('growth_rate', self.growth_rate))
        
        if not targets:
            raise ValueError("At least one target must be provided")
    
    def to_device(self, device: torch.device) -> 'MetabolicTargets':
        """Move all tensors to specified device."""
        return MetabolicTargets(
            node_classification=self.node_classification.to(device) if self.node_classification is not None else None,
            node_regression=self.node_regression.to(device) if self.node_regression is not None else None,
            graph_classification=self.graph_classification.to(device) if self.graph_classification is not None else None,
            graph_regression=self.graph_regression.to(device) if self.graph_regression is not None else None,
            condition_prediction=self.condition_prediction.to(device) if self.condition_prediction is not None else None,
            pathway_analysis=self.pathway_analysis.to(device) if self.pathway_analysis is not None else None,
            growth_rate=self.growth_rate.to(device) if self.growth_rate is not None else None
        )

class MetabolicDataset(Dataset):
    """
    Robust dataset for metabolic network data.
    
    Features:
    - Consistent data structures
    - Proper validation
    - Type safety
    - Error handling
    """
    
    def __init__(self, 
                 graphs: List[MetabolicGraph],
                 targets: List[MetabolicTargets],
                 condition_labels: Optional[List[int]] = None,
                 transform: Optional[callable] = None):
        """
        Initialize the dataset.
        
        Args:
            graphs: List of MetabolicGraph objects
            targets: List of MetabolicTargets objects
            condition_labels: List of condition labels
            transform: Optional data transformation
        """
        if len(graphs) != len(targets):
            raise ValueError(f"Number of graphs ({len(graphs)}) must match number of targets ({len(targets)})")
        
        self.graphs = graphs
        self.targets = targets
        self.condition_labels = condition_labels
        self.transform = transform
        
        # Validate condition labels
        if self.condition_labels is not None and len(self.condition_labels) != len(graphs):
            raise ValueError(f"Number of condition labels ({len(self.condition_labels)}) must match number of graphs ({len(graphs)})")
        
        logger.info(f"Initialized MetabolicDataset with {len(graphs)} samples")
        
        # Log dataset statistics
        self._log_statistics()
    
    def _log_statistics(self):
        """Log dataset statistics for debugging."""
        num_nodes = [g.num_nodes for g in self.graphs]
        num_edges = [g.num_edges for g in self.graphs]
        feature_dims = [g.feature_dim for g in self.graphs]
        
        logger.info(f"Dataset statistics:")
        logger.info(f"  Total samples: {len(self.graphs)}")
        logger.info(f"  Nodes: min={min(num_nodes)}, max={max(num_nodes)}, mean={np.mean(num_nodes):.1f}")
        logger.info(f"  Edges: min={min(num_edges)}, max={max(num_edges)}, mean={np.mean(num_edges):.1f}")
        logger.info(f"  Feature dims: {set(feature_dims)}")
        
        if self.condition_labels is not None:
            unique_conditions = set(self.condition_labels)
            logger.info(f"  Conditions: {unique_conditions}")
    
    def __len__(self):
        return len(self.graphs)
    
    def __getitem__(self, idx):
        """Get a single sample."""
        graph = self.graphs[idx]
        targets = self.targets[idx]
        
        sample = {
            'graph': graph,
            'targets': targets
        }
        
        if self.condition_labels is not None:
            sample['condition_label'] = self.condition_labels[idx]
        
        if self.transform:
            sample = self.transform(sample)
        
        return sample

class DataPreprocessor:
    """
    Robust data preprocessor with proper validation and error handling.
    """
    
    def __init__(self, 
                 feature_scaler: str = 'standard',
                 handle_missing: str = 'mean',
                 normalize_graphs: bool = True):
        """
        Initialize the preprocessor.
        
        Args:
            feature_scaler: Type of feature scaling ('standard', 'minmax', 'none')
            handle_missing: Strategy for missing data ('mean', 'median', 'drop')
            normalize_graphs: Whether to normalize graph structures
        """
        self.feature_scaler = feature_scaler
        self.handle_missing = handle_missing
        self.normalize_graphs = normalize_graphs
        
        # Initialize scalers
        if feature_scaler == 'standard':
            self.scaler = StandardScaler()
        elif feature_scaler == 'minmax':
            self.scaler = MinMaxScaler()
        else:
            self.scaler = None
        
        self.is_fitted = False
        
        logger.info(f"Initialized DataPreprocessor with {feature_scaler} scaling and {handle_missing} missing data handling")
    
    def fit_transform(self, 
                     graphs: List[MetabolicGraph],
                     targets: List[MetabolicTargets]) -> Tuple[List[MetabolicGraph], List[MetabolicTargets]]:
        """
        Fit the preprocessor and transform the data.
        
        Args:
            graphs: List of input graphs
            targets: List of input targets
            
        Returns:
            Preprocessed graphs and targets
        """
        logger.info("Starting data preprocessing...")
        
        # Validate input
        if len(graphs) != len(targets):
            raise ValueError(f"Number of graphs ({len(graphs)}) must match number of targets ({len(targets)})")
        
        # Handle missing data
        processed_graphs = []
        for i, graph in enumerate(graphs):
            try:
                processed_graph = self._process_graph(graph)
                processed_graphs.append(processed_graph)
            except Exception as e:
                logger.error(f"Error processing graph {i}: {str(e)}")
                raise
        
        # Scale features if needed
        if self.scaler is not None:
            processed_graphs = self._scale_features(processed_graphs)
            self.is_fitted = True
        
        # Normalize graph structures if needed
        if self.normalize_graphs:
            processed_graphs = self._normalize_graphs(processed_graphs)
        
        # Validate targets
        processed_targets = []
        for i, target in enumerate(targets):
            try:
                processed_target = self._process_targets(target)
                processed_targets.append(processed_target)
            except Exception as e:
                logger.error(f"Error processing targets {i}: {str(e)}")
                raise
        
        logger.info(f"Data preprocessing completed. Processed {len(processed_graphs)} graphs")
        
        return processed_graphs, processed_targets
    
    def _process_graph(self, graph: MetabolicGraph) -> MetabolicGraph:
        """Process a single graph."""
        # Handle missing data in node features
        node_features = self._handle_missing_data(graph.node_features)
        
        # Create new graph with processed features
        return MetabolicGraph(
            node_features=node_features,
            edge_indices=graph.edge_indices,
            edge_features=graph.edge_features,
            node_types=graph.node_types,
            graph_features=graph.graph_features
        )
    
    def _handle_missing_data(self, data: torch.Tensor) -> torch.Tensor:
        """Handle missing data in tensor."""
        if torch.isnan(data).any():
            logger.warning(f"Found {torch.isnan(data).sum().item()} missing values in data")
            
            if self.handle_missing == 'mean':
                data_clean = torch.nan_to_num(data, nan=torch.nanmean(data))
            elif self.handle_missing == 'median':
                data_clean = torch.nan_to_num(data, nan=torch.nanmedian(data))
            elif self.handle_missing == 'drop':
                # Drop rows with any missing values
                mask = ~torch.isnan(data).any(dim=1)
                data_clean = data[mask]
                logger.warning(f"Dropped {len(data) - len(data_clean)} samples with missing values")
            else:
                raise ValueError(f"Unknown missing data strategy: {self.handle_missing}")
        else:
            data_clean = data
        
        return data_clean
    
    def _scale_features(self, graphs: List[MetabolicGraph]) -> List[MetabolicGraph]:
        """Scale node features across all graphs."""
        # Collect all features
        all_features = torch.cat([g.node_features for g in graphs], dim=0)
        
        # Fit and transform
        scaled_features = torch.tensor(self.scaler.fit_transform(all_features), dtype=torch.float32)
        
        # Split back to individual graphs
        processed_graphs = []
        start_idx = 0
        for graph in graphs:
            end_idx = start_idx + graph.num_nodes
            scaled_graph_features = scaled_features[start_idx:end_idx]
            
            processed_graph = MetabolicGraph(
                node_features=scaled_graph_features,
                edge_indices=graph.edge_indices,
                edge_features=graph.edge_features,
                node_types=graph.node_types,
                graph_features=graph.graph_features
            )
            processed_graphs.append(processed_graph)
            start_idx = end_idx
        
        return processed_graphs
    
    def _normalize_graphs(self, graphs: List[MetabolicGraph]) -> List[MetabolicGraph]:
        """Normalize graph structures (optional)."""
        # For now, just return the graphs as-is
        # This could include edge normalization, node degree normalization, etc.
        return graphs
    
    def _process_targets(self, targets: MetabolicTargets) -> MetabolicTargets:
        """Process target values."""
        # For now, just return targets as-is
        # This could include target scaling, encoding, etc.
        return targets

class GraphBatch:
    """
    Proper batching for graph data.
    
    Handles variable-sized graphs correctly.
    """
    
    def __init__(self, graphs: List[MetabolicGraph]):
        """Initialize batch from list of graphs."""
        self.graphs = graphs
        self.batch_size = len(graphs)
        
        # Create batch indices
        self.batch_indices = self._create_batch_indices()
        
        # Combine node features
        self.node_features = torch.cat([g.node_features for g in graphs], dim=0)
        
        # Combine edge indices with offset
        self.edge_indices = self._combine_edge_indices()
        
        # Combine edge features if available
        self.edge_features = self._combine_edge_features()
        
        # Combine node types if available
        self.node_types = self._combine_node_types()
        
        # Combine graph features if available
        self.graph_features = self._combine_graph_features()
        
        # Set metadata
        self.num_nodes = self.node_features.size(0)
        self.num_edges = self.edge_indices.size(1)
        self.feature_dim = self.node_features.size(1)
    
    def _create_batch_indices(self) -> torch.Tensor:
        """Create batch indices for each node."""
        batch_indices = []
        for i, graph in enumerate(self.graphs):
            batch_indices.extend([i] * graph.num_nodes)
        return torch.tensor(batch_indices, dtype=torch.long)
    
    def _combine_edge_indices(self) -> torch.Tensor:
        """Combine edge indices with proper offset."""
        edge_indices = []
        node_offset = 0
        
        for graph in self.graphs:
            # Offset edge indices
            offset_edges = graph.edge_indices + node_offset
            edge_indices.append(offset_edges)
            node_offset += graph.num_nodes
        
        return torch.cat(edge_indices, dim=1)
    
    def _combine_edge_features(self) -> Optional[torch.Tensor]:
        """Combine edge features if available."""
        edge_features = []
        for graph in self.graphs:
            if graph.edge_features is not None:
                edge_features.append(graph.edge_features)
            else:
                # Create default edge features if missing
                default_features = torch.zeros(graph.num_edges, 1, dtype=torch.float32)
                edge_features.append(default_features)
        
        if edge_features:
            return torch.cat(edge_features, dim=0)
        return None
    
    def _combine_node_types(self) -> Optional[torch.Tensor]:
        """Combine node types if available."""
        node_types = []
        for graph in self.graphs:
            if graph.node_types is not None:
                node_types.append(graph.node_types)
            else:
                # Create default node types if missing
                default_types = torch.zeros(graph.num_nodes, dtype=torch.long)
                node_types.append(default_types)
        
        if node_types:
            return torch.cat(node_types, dim=0)
        return None
    
    def _combine_graph_features(self) -> Optional[torch.Tensor]:
        """Combine graph features if available."""
        graph_features = []
        for graph in self.graphs:
            if graph.graph_features is not None:
                graph_features.append(graph.graph_features)
            else:
                # Create default graph features if missing
                default_features = torch.zeros(1, dtype=torch.float32)
                graph_features.append(default_features)
        
        if graph_features:
            return torch.stack(graph_features, dim=0)
        return None
    
    def to_device(self, device: torch.device) -> 'GraphBatch':
        """Move batch to specified device."""
        batch = GraphBatch.__new__(GraphBatch)
        batch.graphs = self.graphs
        batch.batch_size = self.batch_size
        batch.batch_indices = self.batch_indices.to(device)
        batch.node_features = self.node_features.to(device)
        batch.edge_indices = self.edge_indices.to(device)
        batch.edge_features = self.edge_features.to(device) if self.edge_features is not None else None
        batch.node_types = self.node_types.to(device) if self.node_types is not None else None
        batch.graph_features = self.graph_features.to(device) if self.graph_features is not None else None
        batch.num_nodes = self.num_nodes
        batch.num_edges = self.num_edges
        batch.feature_dim = self.feature_dim
        return batch

def collate_metabolic_data(batch: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Custom collate function for metabolic data.
    
    Handles variable-sized graphs and multiple target types.
    """
    graphs = [item['graph'] for item in batch]
    targets = [item['targets'] for item in batch]
    
    # Create graph batch
    graph_batch = GraphBatch(graphs)
    
    # Combine targets
    combined_targets = _combine_targets(targets)
    
    # Combine condition labels if available
    condition_labels = None
    if 'condition_label' in batch[0]:
        condition_labels = torch.tensor([item['condition_label'] for item in batch], dtype=torch.long)
    
    return {
        'graph_batch': graph_batch,
        'targets': combined_targets,
        'condition_labels': condition_labels
    }

def _combine_targets(targets: List[MetabolicTargets]) -> MetabolicTargets:
    """Combine multiple target sets into a single target set."""
    # Initialize combined targets
    combined = {}
    
    # Node-level targets
    if targets[0].node_classification is not None:
        node_class_tensors = [t.node_classification for t in targets if t.node_classification is not None]
        if node_class_tensors:
            combined['node_classification'] = torch.cat(node_class_tensors, dim=0)
    
    if targets[0].node_regression is not None:
        node_reg_tensors = [t.node_regression for t in targets if t.node_regression is not None]
        if node_reg_tensors:
            combined['node_regression'] = torch.cat(node_reg_tensors, dim=0)
    
    # Graph-level targets
    if targets[0].graph_classification is not None:
        graph_class_tensors = [t.graph_classification for t in targets if t.graph_classification is not None]
        if graph_class_tensors:
            combined['graph_classification'] = torch.stack(graph_class_tensors, dim=0)
    
    if targets[0].graph_regression is not None:
        graph_reg_tensors = [t.graph_regression for t in targets if t.graph_regression is not None]
        if graph_reg_tensors:
            combined['graph_regression'] = torch.cat(graph_reg_tensors, dim=0)
    
    # Other targets
    if targets[0].condition_prediction is not None:
        cond_tensors = [t.condition_prediction for t in targets if t.condition_prediction is not None]
        if cond_tensors:
            combined['condition_prediction'] = torch.stack(cond_tensors, dim=0)
    
    if targets[0].pathway_analysis is not None:
        pathway_tensors = [t.pathway_analysis for t in targets if t.pathway_analysis is not None]
        if pathway_tensors:
            combined['pathway_analysis'] = torch.stack(pathway_tensors, dim=0)
    
    if targets[0].growth_rate is not None:
        growth_tensors = [t.growth_rate for t in targets if t.growth_rate is not None]
        if growth_tensors:
            combined['growth_rate'] = torch.cat(growth_tensors, dim=0)
    
    return MetabolicTargets(**combined)

def create_synthetic_data(num_samples: int = 100,
                         num_nodes_range: Tuple[int, int] = (20, 50),
                         feature_dim: int = 35,
                         edge_probability: float = 0.1) -> Tuple[List[MetabolicGraph], List[MetabolicTargets]]:
    """
    Create synthetic metabolic network data for testing.
    
    Args:
        num_samples: Number of graphs to create
        num_nodes_range: Range for number of nodes per graph
        feature_dim: Dimension of node features
        edge_probability: Probability of edge between any two nodes
        
    Returns:
        List of graphs and targets
    """
    graphs = []
    targets = []
    
    for i in range(num_samples):
        # Random number of nodes
        num_nodes = np.random.randint(num_nodes_range[0], num_nodes_range[1] + 1)
        
        # Create random node features
        node_features = torch.randn(num_nodes, feature_dim, dtype=torch.float32)
        
        # Create random edges
        edge_list = []
        for u in range(num_nodes):
            for v in range(u + 1, num_nodes):
                if np.random.random() < edge_probability:
                    edge_list.append([u, v])
                    edge_list.append([v, u])  # Undirected graph
        
        if not edge_list:
            # Ensure at least some edges
            edge_list = [[0, 1], [1, 0]]
        
        edge_indices = torch.tensor(edge_list, dtype=torch.long).t()
        
        # Create graph
        graph = MetabolicGraph(
            node_features=node_features,
            edge_indices=edge_indices
        )
        graphs.append(graph)
        
        # Create targets
        target = MetabolicTargets(
            node_classification=torch.randint(0, 2, (num_nodes,), dtype=torch.long),
            node_regression=torch.randn(num_nodes, dtype=torch.float32),
            graph_classification=torch.randint(0, 3, (1,), dtype=torch.long),
            graph_regression=torch.randn(1, dtype=torch.float32),
            condition_prediction=torch.randn(3, dtype=torch.float32),
            pathway_analysis=torch.randn(10, dtype=torch.float32),
            growth_rate=torch.randn(1, dtype=torch.float32)
        )
        targets.append(target)
    
    logger.info(f"Created synthetic data: {num_samples} graphs")
    return graphs, targets

def test_data_handling():
    """Test the data handling system thoroughly."""
    logger.info("Testing data handling system...")
    
    # Create synthetic data
    graphs, targets = create_synthetic_data(num_samples=50)
    
    # Test preprocessing
    logger.info("Testing data preprocessing...")
    try:
        preprocessor = DataPreprocessor()
        processed_graphs, processed_targets = preprocessor.fit_transform(graphs, targets)
        logger.info("✅ Data preprocessing successful")
    except Exception as e:
        logger.error(f"❌ Data preprocessing failed: {str(e)}")
        return
    
    # Test dataset creation
    logger.info("Testing dataset creation...")
    try:
        dataset = MetabolicDataset(processed_graphs, processed_targets)
        logger.info("✅ Dataset creation successful")
    except Exception as e:
        logger.error(f"❌ Dataset creation failed: {str(e)}")
        return
    
    # Test data loading
    logger.info("Testing data loading...")
    try:
        dataloader = DataLoader(dataset, batch_size=8, shuffle=True, collate_fn=collate_metabolic_data)
        
        for batch_idx, batch in enumerate(dataloader):
            graph_batch = batch['graph_batch']
            targets = batch['targets']
            
            logger.info(f"Batch {batch_idx}: {graph_batch.batch_size} graphs, "
                       f"{graph_batch.num_nodes} total nodes, {graph_batch.num_edges} total edges")
            
            # Test device movement
            device = torch.device('cpu')
            graph_batch_device = graph_batch.to_device(device)
            targets_device = targets.to_device(device)
            
            if batch_idx >= 2:  # Test first few batches
                break
        
        logger.info("✅ Data loading successful")
    except Exception as e:
        logger.error(f"❌ Data loading failed: {str(e)}")
        return
    
    # Test edge case: empty dataset
    logger.info("Testing edge cases...")
    try:
        empty_dataset = MetabolicDataset([], [])
        logger.info("✅ Empty dataset handling successful")
    except ValueError as e:
        logger.info(f"✅ Empty dataset properly rejected: {str(e)}")
    
    logger.info("🎉 Data handling test completed successfully!")

if __name__ == "__main__":
    test_data_handling() 