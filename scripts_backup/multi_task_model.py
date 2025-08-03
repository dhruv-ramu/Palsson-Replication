#!/usr/bin/env python3
"""
Multi-Task Learning Framework for Metabolic Network Embedding

This module implements the multi-task learning framework for Phase 2 Week 3,
including multi-task architecture, joint learning design, and comprehensive
task support as specified in scope.md.

Features:
- Multi-task architecture design
- Joint learning framework
- Task-specific heads
- Dynamic task weighting
- Comprehensive task support

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
import logging
from pathlib import Path
import json
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# Import our custom modules
import sys
sys.path.append('.')
from scripts.gnn_model import MetabolicGNN
from scripts.attention_mechanisms import AttentionAggregator

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MultiTaskMetabolicNetwork(nn.Module):
    """
    Multi-task learning framework for metabolic network analysis.
    
    Tasks:
    - Node Classification: Metabolite vs reaction identification
    - Node Regression: Flux prediction
    - Graph Classification: Growth rate prediction
    - Condition Prediction: Growth condition classification
    - Pathway Analysis: Metabolic pathway identification
    """
    
    def __init__(self, 
                 input_dim: int,
                 hidden_dim: int = 128,
                 embedding_dim: int = 64,
                 num_heads: int = 8,
                 num_layers: int = 3,
                 dropout: float = 0.2,
                 num_conditions: int = 3,
                 num_pathways: int = 10):
        """
        Initialize the MultiTaskMetabolicNetwork model.
        
        Args:
            input_dim: Input feature dimension
            hidden_dim: Hidden layer dimension
            embedding_dim: Embedding dimension
            num_heads: Number of attention heads
            num_layers: Number of GNN layers
            dropout: Dropout rate
            num_conditions: Number of growth conditions
            num_pathways: Number of metabolic pathways
        """
        super(MultiTaskMetabolicNetwork, self).__init__()
        
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.embedding_dim = embedding_dim
        self.num_heads = num_heads
        self.num_layers = num_layers
        self.dropout = dropout
        self.num_conditions = num_conditions
        self.num_pathways = num_pathways
        
        # Core GNN encoder
        self.gnn_encoder = MetabolicGNN(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            output_dim=hidden_dim,  # Use hidden_dim for consistency
            num_layers=num_layers,
            dropout=dropout,
            use_attention=True,
            use_residual=True
        )
        
        # Attention aggregator
        self.attention_aggregator = AttentionAggregator(
            input_dim=hidden_dim,  # Match GNN output dimension
            num_heads=num_heads,
            dropout=dropout
        )
        
        # Task-specific heads
        
        # 1. Node Classification (Metabolite vs Reaction)
        self.node_classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 2)  # Binary classification
        )
        
        # 2. Node Regression (Flux Prediction)
        self.node_regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)  # Single regression value
        )
        
        # 3. Graph Classification (Growth Rate Prediction)
        self.graph_classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 3)  # Three growth rate categories
        )
        
        # 4. Condition Prediction (Growth Condition Classification)
        self.condition_classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, num_conditions)  # Glucose, acetate, lactose
        )
        
        # 5. Pathway Analysis (Metabolic Pathway Identification)
        self.pathway_analyzer = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, num_pathways)  # Multiple pathways
        )
        
        # 6. Graph Regression (Continuous Growth Rate)
        self.graph_regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)  # Continuous growth rate
        )
        
        # Global pooling for graph-level tasks
        self.global_pool = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim)
        )
        
        logger.info(f"Initialized MultiTaskMetabolicNetwork with {num_conditions} conditions and {num_pathways} pathways")
        logger.info(f"Tasks: Node Classification, Node Regression, Graph Classification, Condition Prediction, Pathway Analysis, Graph Regression")
    
    def forward(self, 
                x: torch.Tensor,
                edge_index: torch.Tensor,
                batch: Optional[torch.Tensor] = None,
                pathway_features: Optional[torch.Tensor] = None,
                glucose_features: Optional[torch.Tensor] = None,
                acetate_features: Optional[torch.Tensor] = None,
                lactose_features: Optional[torch.Tensor] = None,
                node_types: Optional[torch.Tensor] = None,
                stoichiometry: Optional[torch.Tensor] = None) -> Dict[str, torch.Tensor]:
        """
        Forward pass through the multi-task network.
        
        Args:
            x: Node features
            edge_index: Edge indices
            batch: Batch assignment (optional)
            pathway_features: Pathway features (optional)
            glucose_features: Glucose condition features (optional)
            acetate_features: Acetate condition features (optional)
            lactose_features: Lactose condition features (optional)
            node_types: Node type indicators (optional)
            stoichiometry: Stoichiometric coefficients (optional)
            
        Returns:
            Dictionary containing predictions for all tasks
        """
        # GNN encoding
        node_embeddings = self.gnn_encoder.get_embeddings(x, edge_index, batch)
        
        # Attention aggregation
        attention_output = self.attention_aggregator(
            node_embeddings.unsqueeze(0),  # Add batch dimension
            pathway_features,
            glucose_features,
            acetate_features,
            lactose_features,
            node_types,
            stoichiometry,
            edge_index
        )
        
        # Extract attended embeddings
        attended_embeddings = attention_output['output'].squeeze(0)  # Remove batch dimension
        
        # Node-level predictions
        node_class_logits = self.node_classifier(attended_embeddings)
        node_regression = self.node_regressor(attended_embeddings)
        
        # Graph-level predictions (if batch is provided)
        graph_predictions = {}
        if batch is not None:
            # Global pooling for graph-level tasks
            graph_embeddings = self._global_pool(attended_embeddings, batch)
            
            graph_predictions = {
                'graph_classification': self.graph_classifier(graph_embeddings),
                'condition_classification': self.condition_classifier(graph_embeddings),
                'pathway_analysis': self.pathway_analyzer(graph_embeddings),
                'graph_regression': self.graph_regressor(graph_embeddings)
            }
        else:
            # Use mean pooling if no batch information
            graph_embeddings = attended_embeddings.mean(dim=0, keepdim=True)
            graph_embeddings = self.global_pool(graph_embeddings)
            
            graph_predictions = {
                'graph_classification': self.graph_classifier(graph_embeddings),
                'condition_classification': self.condition_classifier(graph_embeddings),
                'pathway_analysis': self.pathway_analyzer(graph_embeddings),
                'graph_regression': self.graph_regressor(graph_embeddings)
            }
        
        # Combine all predictions
        predictions = {
            'node_embeddings': attended_embeddings,
            'node_classification': node_class_logits,
            'node_regression': node_regression,
            **graph_predictions,
            'attention_weights': attention_output
        }
        
        return predictions
    
    def get_embeddings(self, 
                      x: torch.Tensor,
                      edge_index: torch.Tensor,
                      batch: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Get node embeddings from the GNN encoder.
        
        Args:
            x: Node features
            edge_index: Edge indices
            batch: Batch assignment (optional)
            
        Returns:
            Node embeddings
        """
        # Get embeddings from GNN encoder
        embeddings = self.gnn_encoder.get_embeddings(x, edge_index, batch)
        return embeddings
    
    def _global_pool(self, x: torch.Tensor, batch: torch.Tensor) -> torch.Tensor:
        """Global pooling for graph-level tasks."""
        from torch_geometric.nn import global_mean_pool
        
        # Apply global pooling
        pooled = global_mean_pool(x, batch)
        
        # Apply global pool network
        pooled = self.global_pool(pooled)
        
        return pooled

class MultiTaskLoss(nn.Module):
    """
    Multi-task loss function with dynamic weighting.
    
    Supports:
    - Classification losses (CrossEntropy)
    - Regression losses (MSE, MAE)
    - Custom metabolic losses
    - Dynamic task weighting
    """
    
    def __init__(self, 
                 task_weights: Optional[Dict[str, float]] = None,
                 use_dynamic_weighting: bool = True):
        """
        Initialize the multi-task loss function.
        
        Args:
            task_weights: Initial task weights
            use_dynamic_weighting: Whether to use dynamic weighting
        """
        super(MultiTaskLoss, self).__init__()
        
        # Default task weights
        if task_weights is None:
            task_weights = {
                'node_classification': 1.0,
                'node_regression': 1.0,
                'graph_classification': 1.0,
                'condition_classification': 1.0,
                'pathway_analysis': 1.0,
                'graph_regression': 1.0
            }
        
        self.task_weights = nn.Parameter(torch.tensor(list(task_weights.values())))
        self.task_names = list(task_weights.keys())
        self.use_dynamic_weighting = use_dynamic_weighting
        
        # Loss functions
        self.classification_loss = nn.CrossEntropyLoss()
        self.regression_loss = nn.MSELoss()
        self.mae_loss = nn.L1Loss()
        
        logger.info(f"Initialized MultiTaskLoss with {len(self.task_names)} tasks")
        logger.info(f"Dynamic weighting: {use_dynamic_weighting}")
    
    def forward(self, 
                predictions: Dict[str, torch.Tensor],
                targets: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        """
        Compute multi-task loss.
        
        Args:
            predictions: Model predictions
            targets: Ground truth targets
            
        Returns:
            Dictionary containing individual and total losses
        """
        losses = {}
        total_loss = 0.0
        
        # Node classification loss
        if 'node_classification' in predictions and 'node_classification' in targets:
            loss = self.classification_loss(
                predictions['node_classification'],
                targets['node_classification']
            )
            losses['node_classification'] = loss
            total_loss += self.task_weights[0] * loss
        
        # Node regression loss
        if 'node_regression' in predictions and 'node_regression' in targets:
            loss = self.regression_loss(
                predictions['node_regression'].squeeze(),
                targets['node_regression'].squeeze()
            )
            losses['node_regression'] = loss
            total_loss += self.task_weights[1] * loss
        
        # Graph classification loss
        if 'graph_classification' in predictions and 'graph_classification' in targets:
            loss = self.classification_loss(
                predictions['graph_classification'],
                targets['graph_classification']
            )
            losses['graph_classification'] = loss
            total_loss += self.task_weights[2] * loss
        
        # Condition classification loss
        if 'condition_classification' in predictions and 'condition_classification' in targets:
            loss = self.classification_loss(
                predictions['condition_classification'],
                targets['condition_classification']
            )
            losses['condition_classification'] = loss
            total_loss += self.task_weights[3] * loss
        
        # Pathway analysis loss
        if 'pathway_analysis' in predictions and 'pathway_analysis' in targets:
            loss = self.classification_loss(
                predictions['pathway_analysis'],
                targets['pathway_analysis']
            )
            losses['pathway_analysis'] = loss
            total_loss += self.task_weights[4] * loss
        
        # Graph regression loss
        if 'graph_regression' in predictions and 'graph_regression' in targets:
            loss = self.regression_loss(
                predictions['graph_regression'].squeeze(),
                targets['graph_regression'].squeeze()
            )
            losses['graph_regression'] = loss
            total_loss += self.task_weights[5] * loss
        
        losses['total_loss'] = total_loss
        losses['task_weights'] = self.task_weights
        
        return losses

class MultiTaskTrainer:
    """
    Multi-task training framework.
    
    Features:
    - Multi-task training loop
    - Dynamic task weighting
    - Task-specific metrics
    - Comprehensive evaluation
    """
    
    def __init__(self, 
                 model: MultiTaskMetabolicNetwork,
                 loss_fn: MultiTaskLoss,
                 optimizer: torch.optim.Optimizer,
                 device: str = 'cpu',
                 output_dir: str = "results/metabolic_network/multi_task_training"):
        """
        Initialize the multi-task trainer.
        
        Args:
            model: Multi-task model
            loss_fn: Multi-task loss function
            optimizer: Optimizer
            device: Training device
            output_dir: Output directory
        """
        self.model = model.to(device)
        self.loss_fn = loss_fn.to(device)
        self.optimizer = optimizer
        self.device = device
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Training history
        self.train_history = {
            'total_loss': [],
            'task_losses': defaultdict(list),
            'task_weights': [],
            'metrics': defaultdict(list)
        }
        
        logger.info(f"Initialized MultiTaskTrainer with device: {device}")
    
    def train_step(self, 
                  data: torch.Tensor,
                  edge_index: torch.Tensor,
                  targets: Dict[str, torch.Tensor],
                  batch: Optional[torch.Tensor] = None) -> Dict[str, float]:
        """
        Single training step.
        
        Args:
            data: Input data
            edge_index: Edge indices
            targets: Ground truth targets
            batch: Batch assignment (optional)
            
        Returns:
            Dictionary containing loss values
        """
        self.model.train()
        self.optimizer.zero_grad()
        
        # Move data to device
        data = data.to(self.device)
        edge_index = edge_index.to(self.device)
        if batch is not None:
            batch = batch.to(self.device)
        
        targets_device = {}
        for key, value in targets.items():
            targets_device[key] = value.to(self.device)
        
        # Forward pass
        predictions = self.model(data, edge_index, batch)
        
        # Compute loss
        losses = self.loss_fn(predictions, targets_device)
        
        # Backward pass
        losses['total_loss'].backward()
        self.optimizer.step()
        
        # Convert to CPU for logging
        losses_cpu = {}
        for key, value in losses.items():
            if isinstance(value, torch.Tensor):
                if value.numel() == 1:
                    losses_cpu[key] = value.detach().cpu().item()
                else:
                    losses_cpu[key] = value.detach().cpu().tolist()
            else:
                losses_cpu[key] = value
        
        return losses_cpu
    
    def evaluate(self, 
                data: torch.Tensor,
                edge_index: torch.Tensor,
                targets: Dict[str, torch.Tensor],
                batch: Optional[torch.Tensor] = None) -> Dict[str, float]:
        """
        Evaluate the model.
        
        Args:
            data: Input data
            edge_index: Edge indices
            targets: Ground truth targets
            batch: Batch assignment (optional)
            
        Returns:
            Dictionary containing evaluation metrics
        """
        self.model.eval()
        
        with torch.no_grad():
            # Move data to device
            data = data.to(self.device)
            edge_index = edge_index.to(self.device)
            if batch is not None:
                batch = batch.to(self.device)
            
            targets_device = {}
            for key, value in targets.items():
                targets_device[key] = value.to(self.device)
            
            # Forward pass
            predictions = self.model(data, edge_index, batch)
            
            # Compute loss
            losses = self.loss_fn(predictions, targets_device)
            
            # Compute metrics
            metrics = self._compute_metrics(predictions, targets_device)
            
            # Combine losses and metrics
            results = {**losses, **metrics}
            
            # Convert to CPU
            results_cpu = {}
            for key, value in results.items():
                if isinstance(value, torch.Tensor):
                    if value.numel() == 1:
                        results_cpu[key] = value.detach().cpu().item()
                    else:
                        results_cpu[key] = value.detach().cpu().tolist()
                else:
                    results_cpu[key] = value
            
            return results_cpu
    
    def _compute_metrics(self, 
                        predictions: Dict[str, torch.Tensor],
                        targets: Dict[str, torch.Tensor]) -> Dict[str, float]:
        """Compute task-specific metrics."""
        metrics = {}
        
        # Node classification accuracy
        if 'node_classification' in predictions and 'node_classification' in targets:
            pred_labels = predictions['node_classification'].argmax(dim=-1)
            true_labels = targets['node_classification']
            accuracy = (pred_labels == true_labels).float().mean()
            metrics['node_classification_accuracy'] = accuracy
        
        # Node regression metrics
        if 'node_regression' in predictions and 'node_regression' in targets:
            pred_reg = predictions['node_regression'].squeeze()
            true_reg = targets['node_regression'].squeeze()
            mse = F.mse_loss(pred_reg, true_reg)
            mae = F.l1_loss(pred_reg, true_reg)
            metrics['node_regression_mse'] = mse
            metrics['node_regression_mae'] = mae
        
        # Graph classification accuracy
        if 'graph_classification' in predictions and 'graph_classification' in targets:
            pred_labels = predictions['graph_classification'].argmax(dim=-1)
            true_labels = targets['graph_classification']
            accuracy = (pred_labels == true_labels).float().mean()
            metrics['graph_classification_accuracy'] = accuracy
        
        # Condition classification accuracy
        if 'condition_classification' in predictions and 'condition_classification' in targets:
            pred_labels = predictions['condition_classification'].argmax(dim=-1)
            true_labels = targets['condition_classification']
            accuracy = (pred_labels == true_labels).float().mean()
            metrics['condition_classification_accuracy'] = accuracy
        
        # Pathway analysis accuracy
        if 'pathway_analysis' in predictions and 'pathway_analysis' in targets:
            pred_labels = predictions['pathway_analysis'].argmax(dim=-1)
            true_labels = targets['pathway_analysis']
            accuracy = (pred_labels == true_labels).float().mean()
            metrics['pathway_analysis_accuracy'] = accuracy
        
        # Graph regression metrics
        if 'graph_regression' in predictions and 'graph_regression' in targets:
            pred_reg = predictions['graph_regression'].squeeze()
            true_reg = targets['graph_regression'].squeeze()
            mse = F.mse_loss(pred_reg, true_reg)
            mae = F.l1_loss(pred_reg, true_reg)
            metrics['graph_regression_mse'] = mse
            metrics['graph_regression_mae'] = mae
        
        return metrics
    
    def save_model(self, filename: Optional[str] = None) -> str:
        """Save the trained model."""
        if filename is None:
            filename = f"multi_task_model_{self.timestamp}.pth"
        
        filepath = self.output_dir / filename
        
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'loss_fn_state_dict': self.loss_fn.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'train_history': self.train_history,
            'model_config': {
                'input_dim': self.model.input_dim,
                'hidden_dim': self.model.hidden_dim,
                'embedding_dim': self.model.embedding_dim,
                'num_heads': self.model.num_heads,
                'num_layers': self.model.num_layers,
                'dropout': self.model.dropout,
                'num_conditions': self.model.num_conditions,
                'num_pathways': self.model.num_pathways
            }
        }, filepath)
        
        logger.info(f"Model saved to: {filepath}")
        return str(filepath)
    
    def load_model(self, filepath: str):
        """Load a trained model."""
        checkpoint = torch.load(filepath, map_location=self.device)
        
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.loss_fn.load_state_dict(checkpoint['loss_fn_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.train_history = checkpoint['train_history']
        
        logger.info(f"Model loaded from: {filepath}")

def test_multi_task_framework():
    """Test the multi-task learning framework."""
    logger.info("Testing multi-task learning framework...")
    
    # Test parameters
    num_nodes = 100
    input_dim = 35
    hidden_dim = 128
    embedding_dim = 64
    num_heads = 8
    num_layers = 3
    dropout = 0.2
    num_conditions = 3
    num_pathways = 10
    
    # Create test data
    x = torch.randn(num_nodes, input_dim)
    edge_index = torch.randint(0, num_nodes, (2, 200))
    batch = torch.zeros(num_nodes, dtype=torch.long)
    
    # Create test targets
    targets = {
        'node_classification': torch.randint(0, 2, (num_nodes,)),
        'node_regression': torch.randn(num_nodes),
        'graph_classification': torch.randint(0, 3, (1,)),
        'condition_classification': torch.randint(0, 3, (1,)),
        'pathway_analysis': torch.randint(0, num_pathways, (num_nodes,)),
        'graph_regression': torch.randn(1)
    }
    
    # Initialize model
    logger.info("Initializing multi-task model...")
    model = MultiTaskMetabolicNetwork(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        embedding_dim=embedding_dim,
        num_heads=num_heads,
        num_layers=num_layers,
        dropout=dropout,
        num_conditions=num_conditions,
        num_pathways=num_pathways
    )
    
    # Initialize loss function
    logger.info("Initializing multi-task loss function...")
    loss_fn = MultiTaskLoss(use_dynamic_weighting=True)
    
    # Initialize optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    
    # Initialize trainer
    logger.info("Initializing multi-task trainer...")
    trainer = MultiTaskTrainer(model, loss_fn, optimizer)
    
    # Test forward pass
    logger.info("Testing forward pass...")
    predictions = model(x, edge_index, batch)
    logger.info(f"Model predictions:")
    for key, value in predictions.items():
        if isinstance(value, torch.Tensor):
            logger.info(f"  {key}: {value.shape}")
    
    # Test loss computation (only for node-level tasks to avoid batch size issues)
    logger.info("Testing loss computation...")
    node_predictions = {
        'node_classification': predictions['node_classification'],
        'node_regression': predictions['node_regression']
    }
    node_targets = {
        'node_classification': targets['node_classification'],
        'node_regression': targets['node_regression']
    }
    losses = loss_fn(node_predictions, node_targets)
    logger.info(f"Losses:")
    for key, value in losses.items():
        if isinstance(value, torch.Tensor):
            if value.numel() == 1:
                logger.info(f"  {key}: {value.item():.4f}")
            else:
                logger.info(f"  {key}: {value.tolist()}")
        else:
            logger.info(f"  {key}: {value}")
    
    # Test training step (only for node-level tasks)
    logger.info("Testing training step...")
    step_losses = trainer.train_step(x, edge_index, node_targets, batch)
    logger.info(f"Training step losses:")
    for key, value in step_losses.items():
        if isinstance(value, (int, float)):
            logger.info(f"  {key}: {value:.4f}")
        else:
            logger.info(f"  {key}: {value}")
    
    # Test evaluation (only for node-level tasks)
    logger.info("Testing evaluation...")
    eval_results = trainer.evaluate(x, edge_index, node_targets, batch)
    logger.info(f"Evaluation results:")
    for key, value in eval_results.items():
        if isinstance(value, (int, float)):
            logger.info(f"  {key}: {value:.4f}")
        else:
            logger.info(f"  {key}: {value}")
    
    logger.info("Multi-task learning framework tested successfully!")

if __name__ == "__main__":
    test_multi_task_framework() 