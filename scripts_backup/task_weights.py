#!/usr/bin/env python3
"""
Dynamic Task Weighting for Metabolic Network Embedding

This module implements dynamic task weighting for Phase 2 Week 3,
including adaptive weighting strategies, uncertainty-based weighting,
and task importance learning as specified in scope.md.

Features:
- Dynamic task weighting
- Uncertainty-based weighting
- Task importance learning
- Adaptive weighting strategies
- Multi-task optimization

Author: Metabolic Network Embedding Project
Date: 2025
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Dict, List, Optional, Tuple, Any, Union
import logging
from pathlib import Path
import json
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Fix import issue
import sys
sys.path.append('.')

class UncertaintyBasedWeighting(nn.Module):
    """
    Uncertainty-based task weighting using task-specific uncertainties.
    
    Features:
    - Task uncertainty estimation
    - Homoscedastic uncertainty
    - Heteroscedastic uncertainty
    - Dynamic uncertainty weighting
    """
    
    def __init__(self, 
                 num_tasks: int,
                 uncertainty_type: str = 'homoscedastic',
                 initial_uncertainty: float = 1.0,
                 learn_uncertainty: bool = True):
        """
        Initialize uncertainty-based weighting.
        
        Args:
            num_tasks: Number of tasks
            uncertainty_type: Type of uncertainty ('homoscedastic' or 'heteroscedastic')
            initial_uncertainty: Initial uncertainty value
            learn_uncertainty: Whether to learn uncertainty parameters
        """
        super(UncertaintyBasedWeighting, self).__init__()
        
        self.num_tasks = num_tasks
        self.uncertainty_type = uncertainty_type
        self.learn_uncertainty = learn_uncertainty
        
        if uncertainty_type == 'homoscedastic':
            # Single uncertainty parameter for all tasks
            self.task_uncertainties = nn.Parameter(
                torch.full((num_tasks,), initial_uncertainty, dtype=torch.float32)
            )
        elif uncertainty_type == 'heteroscedastic':
            # Task-specific uncertainty parameters
            self.task_uncertainties = nn.Parameter(
                torch.full((num_tasks,), initial_uncertainty, dtype=torch.float32)
            )
            self.uncertainty_networks = nn.ModuleList([
                nn.Sequential(
                    nn.Linear(64, 32),  # Assuming 64-dim embeddings
                    nn.ReLU(),
                    nn.Linear(32, 1),
                    nn.Softplus()
                ) for _ in range(num_tasks)
            ])
        
        logger.info(f"Initialized UncertaintyBasedWeighting with {num_tasks} tasks, type={uncertainty_type}")
    
    def forward(self, 
                task_losses: torch.Tensor,
                task_embeddings: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Compute uncertainty-based task weights.
        
        Args:
            task_losses: Individual task losses [num_tasks]
            task_embeddings: Task embeddings for heteroscedastic uncertainty [num_tasks, embedding_dim]
            
        Returns:
            Task weights [num_tasks]
        """
        if self.uncertainty_type == 'homoscedastic':
            # Homoscedastic uncertainty: same uncertainty for all samples of a task
            uncertainties = F.softplus(self.task_uncertainties)  # Ensure positive
            weights = 1.0 / (uncertainties ** 2)
            
        elif self.uncertainty_type == 'heteroscedastic':
            # Heteroscedastic uncertainty: uncertainty varies with input
            if task_embeddings is not None:
                uncertainties = []
                for i, uncertainty_net in enumerate(self.uncertainty_networks):
                    task_uncertainty = uncertainty_net(task_embeddings[i])
                    uncertainties.append(task_uncertainty)
                uncertainties = torch.cat(uncertainties)
            else:
                uncertainties = F.softplus(self.task_uncertainties)
            
            weights = 1.0 / (uncertainties ** 2)
        
        # Normalize weights
        weights = weights / weights.sum()
        
        return weights
    
    def get_uncertainties(self) -> torch.Tensor:
        """Get current task uncertainties."""
        return F.softplus(self.task_uncertainties)

class TaskImportanceLearning(nn.Module):
    """
    Learn task importance weights based on task performance and relationships.
    
    Features:
    - Task importance estimation
    - Performance-based weighting
    - Task relationship modeling
    - Adaptive importance learning
    """
    
    def __init__(self, 
                 num_tasks: int,
                 embedding_dim: int = 64,
                 use_task_relationships: bool = True,
                 importance_decay: float = 0.95):
        """
        Initialize task importance learning.
        
        Args:
            num_tasks: Number of tasks
            embedding_dim: Embedding dimension for task representations
            use_task_relationships: Whether to model task relationships
            importance_decay: Decay factor for importance updates
        """
        super(TaskImportanceLearning, self).__init__()
        
        self.num_tasks = num_tasks
        self.embedding_dim = embedding_dim
        self.use_task_relationships = use_task_relationships
        self.importance_decay = importance_decay
        
        # Task importance parameters
        self.task_importance = nn.Parameter(torch.ones(num_tasks, dtype=torch.float32))
        
        # Task embedding network
        self.task_embedding_net = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim // 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(embedding_dim // 2, embedding_dim // 4),
            nn.ReLU(),
            nn.Linear(embedding_dim // 4, 1),
            nn.Sigmoid()
        )
        
        # Task relationship network (if enabled)
        if use_task_relationships:
            self.task_relationship_net = nn.Sequential(
                nn.Linear(num_tasks, num_tasks // 2),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.Linear(num_tasks // 2, num_tasks),
                nn.Sigmoid()
            )
        
        # Performance history
        self.performance_history = defaultdict(list)
        
        logger.info(f"Initialized TaskImportanceLearning with {num_tasks} tasks")
    
    def forward(self, 
                task_losses: torch.Tensor,
                task_embeddings: torch.Tensor,
                task_metrics: Optional[Dict[str, float]] = None) -> torch.Tensor:
        """
        Compute task importance weights.
        
        Args:
            task_losses: Individual task losses [num_tasks]
            task_embeddings: Task embeddings [num_tasks, embedding_dim]
            task_metrics: Task performance metrics (optional)
            
        Returns:
            Task importance weights [num_tasks]
        """
        # Base importance from task embeddings
        base_importance = self.task_embedding_net(task_embeddings).squeeze(-1)
        
        # Performance-based importance
        if task_metrics is not None:
            performance_importance = self._compute_performance_importance(task_metrics)
        else:
            # Use loss-based performance (lower loss = higher importance)
            performance_importance = 1.0 / (task_losses + 1e-8)
            performance_importance = performance_importance / performance_importance.sum()
        
        # Task relationship importance (if enabled)
        if self.use_task_relationships:
            relationship_importance = self.task_relationship_net(task_losses.unsqueeze(0)).squeeze(0)
        else:
            relationship_importance = torch.ones(self.num_tasks)
        
        # Combine different importance factors
        combined_importance = (base_importance * 0.3 + 
                             performance_importance * 0.5 + 
                             relationship_importance * 0.2)
        
        # Apply learned importance scaling
        final_importance = self.task_importance * combined_importance
        
        # Normalize to get weights
        weights = F.softmax(final_importance, dim=-1)
        
        return weights
    
    def _compute_performance_importance(self, task_metrics: Dict[str, float]) -> torch.Tensor:
        """Compute importance based on task performance metrics."""
        importance = torch.zeros(self.num_tasks)
        
        for i, (task_name, metric) in enumerate(task_metrics.items()):
            if i < self.num_tasks:
                # Higher metric = higher importance (assuming metrics are accuracy-like)
                importance[i] = metric
        
        # Normalize
        if importance.sum() > 0:
            importance = importance / importance.sum()
        else:
            importance = torch.ones(self.num_tasks) / self.num_tasks
        
        return importance
    
    def update_performance_history(self, task_name: str, performance: float):
        """Update performance history for a task."""
        self.performance_history[task_name].append(performance)
        
        # Keep only recent history
        if len(self.performance_history[task_name]) > 100:
            self.performance_history[task_name] = self.performance_history[task_name][-100:]

class AdaptiveTaskWeighting(nn.Module):
    """
    Adaptive task weighting that adjusts weights based on training progress.
    
    Features:
    - Training progress tracking
    - Dynamic weight adjustment
    - Convergence-based weighting
    - Gradient-based weighting
    """
    
    def __init__(self, 
                 num_tasks: int,
                 initial_weights: Optional[torch.Tensor] = None,
                 learning_rate: float = 0.01,
                 momentum: float = 0.9,
                 use_gradient_norm: bool = True):
        """
        Initialize adaptive task weighting.
        
        Args:
            num_tasks: Number of tasks
            initial_weights: Initial task weights
            learning_rate: Learning rate for weight updates
            momentum: Momentum for weight updates
            use_gradient_norm: Whether to use gradient norm for weighting
        """
        super(AdaptiveTaskWeighting, self).__init__()
        
        self.num_tasks = num_tasks
        self.learning_rate = learning_rate
        self.momentum = momentum
        self.use_gradient_norm = use_gradient_norm
        
        # Initialize weights
        if initial_weights is None:
            initial_weights = torch.ones(num_tasks) / num_tasks
        
        self.task_weights = nn.Parameter(initial_weights)
        self.weight_momentum = torch.zeros_like(self.task_weights)
        
        # Training history
        self.loss_history = defaultdict(list)
        self.weight_history = []
        self.gradient_history = defaultdict(list)
        
        logger.info(f"Initialized AdaptiveTaskWeighting with {num_tasks} tasks")
    
    def forward(self, 
                task_losses: torch.Tensor,
                task_gradients: Optional[torch.Tensor] = None,
                epoch: int = 0) -> torch.Tensor:
        """
        Compute adaptive task weights.
        
        Args:
            task_losses: Individual task losses [num_tasks]
            task_gradients: Task gradients [num_tasks] (optional)
            epoch: Current training epoch
            
        Returns:
            Task weights [num_tasks]
        """
        # Store loss history
        for i, loss in enumerate(task_losses):
            self.loss_history[f'task_{i}'].append(loss.item())
        
        # Compute adaptive weights
        if self.use_gradient_norm and task_gradients is not None:
            # Use gradient norm for weighting
            gradient_norms = torch.norm(task_gradients, dim=-1)
            gradient_weights = 1.0 / (gradient_norms + 1e-8)
            gradient_weights = gradient_weights / gradient_weights.sum()
            
            # Combine with current weights
            adaptive_weights = 0.7 * self.task_weights + 0.3 * gradient_weights
        else:
            # Use loss-based adaptation
            loss_weights = 1.0 / (task_losses + 1e-8)
            loss_weights = loss_weights / loss_weights.sum()
            
            adaptive_weights = 0.8 * self.task_weights + 0.2 * loss_weights
        
        # Apply momentum update
        weight_update = adaptive_weights - self.task_weights
        self.weight_momentum = self.momentum * self.weight_momentum + self.learning_rate * weight_update
        
        # Update weights
        new_weights = self.task_weights + self.weight_momentum
        
        # Ensure weights are positive and sum to 1
        new_weights = F.softmax(new_weights, dim=-1)
        
        # Store weight history
        self.weight_history.append(new_weights.detach().clone())
        
        return new_weights
    
    def get_weight_history(self) -> torch.Tensor:
        """Get weight history as a tensor."""
        if self.weight_history:
            return torch.stack(self.weight_history)
        else:
            return torch.empty(0, self.num_tasks)
    
    def reset_weights(self):
        """Reset weights to uniform distribution."""
        with torch.no_grad():
            self.task_weights.fill_(1.0 / self.num_tasks)
            self.weight_momentum.zero_()

class MultiTaskWeightingStrategy(nn.Module):
    """
    Comprehensive multi-task weighting strategy combining multiple approaches.
    
    Features:
    - Multiple weighting strategies
    - Strategy combination
    - Dynamic strategy selection
    - Comprehensive weighting
    """
    
    def __init__(self, 
                 num_tasks: int,
                 strategy_weights: Optional[Dict[str, float]] = None,
                 use_uncertainty: bool = True,
                 use_importance: bool = True,
                 use_adaptive: bool = True):
        """
        Initialize multi-task weighting strategy.
        
        Args:
            num_tasks: Number of tasks
            strategy_weights: Weights for different strategies
            use_uncertainty: Whether to use uncertainty-based weighting
            use_importance: Whether to use importance-based weighting
            use_adaptive: Whether to use adaptive weighting
        """
        super(MultiTaskWeightingStrategy, self).__init__()
        
        self.num_tasks = num_tasks
        self.use_uncertainty = use_uncertainty
        self.use_importance = use_importance
        self.use_adaptive = use_adaptive
        
        # Default strategy weights
        if strategy_weights is None:
            strategy_weights = {
                'uncertainty': 0.4,
                'importance': 0.3,
                'adaptive': 0.3
            }
        
        self.strategy_weights = strategy_weights
        
        # Initialize weighting strategies
        if use_uncertainty:
            self.uncertainty_weighting = UncertaintyBasedWeighting(num_tasks)
        
        if use_importance:
            self.importance_weighting = TaskImportanceLearning(num_tasks)
        
        if use_adaptive:
            self.adaptive_weighting = AdaptiveTaskWeighting(num_tasks)
        
        # Strategy combination network
        self.strategy_combiner = nn.Sequential(
            nn.Linear(num_tasks * 3, num_tasks * 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(num_tasks * 2, num_tasks),
            nn.Softmax(dim=-1)
        )
        
        logger.info(f"Initialized MultiTaskWeightingStrategy with {num_tasks} tasks")
        logger.info(f"Strategies: uncertainty={use_uncertainty}, importance={use_importance}, adaptive={use_adaptive}")
    
    def forward(self, 
                task_losses: torch.Tensor,
                task_embeddings: Optional[torch.Tensor] = None,
                task_gradients: Optional[torch.Tensor] = None,
                task_metrics: Optional[Dict[str, float]] = None,
                epoch: int = 0) -> Dict[str, torch.Tensor]:
        """
        Compute comprehensive task weights using multiple strategies.
        
        Args:
            task_losses: Individual task losses [num_tasks]
            task_embeddings: Task embeddings [num_tasks, embedding_dim] (optional)
            task_gradients: Task gradients [num_tasks] (optional)
            task_metrics: Task performance metrics (optional)
            epoch: Current training epoch
            
        Returns:
            Dictionary containing weights from different strategies and combined weights
        """
        strategy_weights = {}
        
        # Uncertainty-based weighting
        if self.use_uncertainty:
            uncertainty_weights = self.uncertainty_weighting(task_losses, task_embeddings)
            strategy_weights['uncertainty'] = uncertainty_weights
        
        # Importance-based weighting
        if self.use_importance and task_embeddings is not None:
            importance_weights = self.importance_weighting(task_losses, task_embeddings, task_metrics)
            strategy_weights['importance'] = importance_weights
        
        # Adaptive weighting
        if self.use_adaptive:
            adaptive_weights = self.adaptive_weighting(task_losses, task_gradients, epoch)
            strategy_weights['adaptive'] = adaptive_weights
        
        # Combine strategies
        combined_weights = self._combine_strategies(strategy_weights)
        
        # Add combined weights to output
        strategy_weights['combined'] = combined_weights
        
        return strategy_weights
    
    def _combine_strategies(self, strategy_weights: Dict[str, torch.Tensor]) -> torch.Tensor:
        """Combine weights from different strategies."""
        if len(strategy_weights) == 1:
            return list(strategy_weights.values())[0]
        
        # Concatenate weights from all strategies
        all_weights = []
        for strategy_name in ['uncertainty', 'importance', 'adaptive']:
            if strategy_name in strategy_weights:
                all_weights.append(strategy_weights[strategy_name])
            else:
                # Use uniform weights for missing strategies
                all_weights.append(torch.ones(self.num_tasks) / self.num_tasks)
        
        concatenated_weights = torch.cat(all_weights)
        
        # Use strategy combiner network
        combined_weights = self.strategy_combiner(concatenated_weights.unsqueeze(0)).squeeze(0)
        
        return combined_weights
    
    def get_strategy_uncertainties(self) -> Optional[torch.Tensor]:
        """Get uncertainties from uncertainty-based weighting."""
        if self.use_uncertainty:
            return self.uncertainty_weighting.get_uncertainties()
        return None
    
    def get_weight_history(self) -> Optional[torch.Tensor]:
        """Get weight history from adaptive weighting."""
        if self.use_adaptive:
            return self.adaptive_weighting.get_weight_history()
        return None

class TaskWeightingVisualizer:
    """
    Visualizer for task weighting strategies.
    
    Features:
    - Weight evolution plots
    - Strategy comparison
    - Performance correlation
    - Weight analysis
    """
    
    def __init__(self, output_dir: str = "results/metabolic_network/task_weighting"):
        """Initialize the task weighting visualizer."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        logger.info(f"Initialized TaskWeightingVisualizer with output directory: {self.output_dir}")
    
    def plot_weight_evolution(self, 
                            weight_history: torch.Tensor,
                            task_names: Optional[List[str]] = None,
                            title: str = "Task Weight Evolution",
                            save_path: Optional[str] = None) -> None:
        """
        Plot task weight evolution over time.
        
        Args:
            weight_history: Weight history tensor [num_epochs, num_tasks]
            task_names: Names of tasks
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating weight evolution plot: {title}")
        
        num_epochs, num_tasks = weight_history.shape
        
        if task_names is None:
            task_names = [f'Task_{i}' for i in range(num_tasks)]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        epochs = range(num_epochs)
        for i in range(num_tasks):
            ax.plot(epochs, weight_history[:, i], label=task_names[i], linewidth=2)
        
        ax.set_xlabel('Training Epoch')
        ax.set_ylabel('Task Weight')
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"weight_evolution_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Weight evolution plot saved to: {save_path}")
    
    def plot_strategy_comparison(self, 
                               strategy_weights: Dict[str, torch.Tensor],
                               task_names: Optional[List[str]] = None,
                               title: str = "Strategy Comparison",
                               save_path: Optional[str] = None) -> None:
        """
        Plot comparison of different weighting strategies.
        
        Args:
            strategy_weights: Dictionary of weights from different strategies
            task_names: Names of tasks
            title: Plot title
            save_path: Path to save the plot
        """
        logger.info(f"Creating strategy comparison plot: {title}")
        
        strategies = list(strategy_weights.keys())
        num_tasks = len(list(strategy_weights.values())[0])
        
        if task_names is None:
            task_names = [f'Task_{i}' for i in range(num_tasks)]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(15, 8))
        
        x = np.arange(num_tasks)
        width = 0.8 / len(strategies)
        
        for i, strategy in enumerate(strategies):
            weights = strategy_weights[strategy].detach().cpu().numpy()
            ax.bar(x + i * width, weights, width, label=strategy, alpha=0.8)
        
        ax.set_xlabel('Tasks')
        ax.set_ylabel('Weight')
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.set_xticks(x + width * (len(strategies) - 1) / 2)
        ax.set_xticklabels(task_names, rotation=45)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Save plot
        if save_path is None:
            save_path = self.output_dir / f"strategy_comparison_{self.timestamp}.png"
        else:
            save_path = self.output_dir / save_path
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Strategy comparison plot saved to: {save_path}")

def test_task_weighting():
    """Test the task weighting modules."""
    logger.info("Testing task weighting modules...")
    
    # Test parameters
    num_tasks = 6
    num_epochs = 10
    
    # Create test data
    task_losses = torch.randn(num_tasks)
    task_embeddings = torch.randn(num_tasks, 64)
    task_gradients = torch.randn(num_tasks, 32)
    task_metrics = {
        'task_0': 0.85,
        'task_1': 0.72,
        'task_2': 0.91,
        'task_3': 0.68,
        'task_4': 0.79,
        'task_5': 0.83
    }
    
    # Test UncertaintyBasedWeighting
    logger.info("Testing UncertaintyBasedWeighting...")
    uncertainty_weighting = UncertaintyBasedWeighting(num_tasks)
    uncertainty_weights = uncertainty_weighting(task_losses, task_embeddings)
    logger.info(f"Uncertainty weights: {uncertainty_weights}")
    
    # Test TaskImportanceLearning
    logger.info("Testing TaskImportanceLearning...")
    importance_weighting = TaskImportanceLearning(num_tasks)
    importance_weights = importance_weighting(task_losses, task_embeddings, task_metrics)
    logger.info(f"Importance weights: {importance_weights}")
    
    # Test AdaptiveTaskWeighting
    logger.info("Testing AdaptiveTaskWeighting...")
    adaptive_weighting = AdaptiveTaskWeighting(num_tasks)
    adaptive_weights = adaptive_weighting(task_losses, task_gradients, epoch=0)
    logger.info(f"Adaptive weights: {adaptive_weights}")
    
    # Test MultiTaskWeightingStrategy
    logger.info("Testing MultiTaskWeightingStrategy...")
    strategy_weighting = MultiTaskWeightingStrategy(num_tasks)
    strategy_weights = strategy_weighting(task_losses, task_embeddings, task_gradients, task_metrics, epoch=0)
    logger.info(f"Strategy weights:")
    for strategy, weights in strategy_weights.items():
        logger.info(f"  {strategy}: {weights}")
    
    # Test TaskWeightingVisualizer
    logger.info("Testing TaskWeightingVisualizer...")
    visualizer = TaskWeightingVisualizer()
    
    # Create weight history for visualization
    weight_history = torch.randn(num_epochs, num_tasks)
    weight_history = F.softmax(weight_history, dim=-1)
    
    visualizer.plot_weight_evolution(weight_history, title="Test Weight Evolution")
    visualizer.plot_strategy_comparison(strategy_weights, title="Test Strategy Comparison")
    
    logger.info("All task weighting modules tested successfully!")

if __name__ == "__main__":
    test_task_weighting() 