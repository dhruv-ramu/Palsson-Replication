#!/usr/bin/env python3
"""
Custom Loss Functions for Metabolic Network Embedding

This module implements custom loss functions for Phase 2 Week 3,
including task-specific objectives, metabolic-specific losses, and
advanced loss functions as specified in scope.md.

Features:
- Task-specific loss functions
- Metabolic-specific losses
- Advanced loss functions
- Loss function combinations
- Custom loss metrics

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

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Fix import issue
import sys
sys.path.append('.')

class MetabolicClassificationLoss(nn.Module):
    """
    Custom classification loss for metabolic network tasks.
    
    Features:
    - Node type classification (metabolite vs reaction)
    - Growth condition classification
    - Pathway classification
    - Focal loss for imbalanced classes
    - Label smoothing
    """
    
    def __init__(self, 
                 alpha: float = 1.0,
                 gamma: float = 2.0,
                 label_smoothing: float = 0.1,
                 class_weights: Optional[torch.Tensor] = None):
        """
        Initialize the metabolic classification loss.
        
        Args:
            alpha: Focal loss alpha parameter
            gamma: Focal loss gamma parameter
            label_smoothing: Label smoothing factor
            class_weights: Class weights for imbalanced data
        """
        super(MetabolicClassificationLoss, self).__init__()
        
        self.alpha = alpha
        self.gamma = gamma
        self.label_smoothing = label_smoothing
        self.class_weights = class_weights
        
        # Standard cross-entropy loss
        self.ce_loss = nn.CrossEntropyLoss(
            weight=class_weights,
            label_smoothing=label_smoothing
        )
        
        logger.info(f"Initialized MetabolicClassificationLoss with alpha={alpha}, gamma={gamma}")
    
    def forward(self, 
                predictions: torch.Tensor,
                targets: torch.Tensor,
                use_focal: bool = True) -> torch.Tensor:
        """
        Compute metabolic classification loss.
        
        Args:
            predictions: Model predictions [batch_size, num_classes]
            targets: Ground truth labels [batch_size]
            use_focal: Whether to use focal loss
            
        Returns:
            Loss value
        """
        if use_focal:
            return self._focal_loss(predictions, targets)
        else:
            return self.ce_loss(predictions, targets)
    
    def _focal_loss(self, predictions: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        """Compute focal loss for imbalanced classification."""
        # Apply softmax to get probabilities
        probs = F.softmax(predictions, dim=-1)
        
        # Get probability of correct class
        batch_size = predictions.size(0)
        probs_correct = probs[torch.arange(batch_size), targets]
        
        # Compute focal loss
        focal_weight = (1 - probs_correct) ** self.gamma
        
        # Apply alpha weighting if class weights are provided
        if self.class_weights is not None:
            alpha_weight = self.class_weights[targets]
            focal_weight = alpha_weight * focal_weight
        
        # Compute cross-entropy
        ce_loss = F.cross_entropy(predictions, targets, reduction='none')
        
        # Apply focal weighting
        focal_loss = focal_weight * ce_loss
        
        return focal_loss.mean()

class MetabolicRegressionLoss(nn.Module):
    """
    Custom regression loss for metabolic network tasks.
    
    Features:
    - Flux prediction regression
    - Growth rate regression
    - Huber loss for robustness
    - Relative error support
    """
    
    def __init__(self, 
                 delta: float = 1.0,
                 use_huber: bool = True):
        """
        Initialize the metabolic regression loss.
        
        Args:
            delta: Huber loss delta parameter
            use_huber: Whether to use Huber loss
        """
        super(MetabolicRegressionLoss, self).__init__()
        
        self.delta = delta
        self.use_huber = use_huber
        
        # Standard MSE loss
        self.mse_loss = nn.MSELoss()
        
        logger.info(f"Initialized MetabolicRegressionLoss with delta={delta}, use_huber={use_huber}")
    
    def forward(self, 
                predictions: torch.Tensor,
                targets: torch.Tensor,
                use_huber: bool = True,
                use_relative: bool = False) -> torch.Tensor:
        """
        Compute metabolic regression loss.
        
        Args:
            predictions: Model predictions
            targets: Ground truth values
            use_huber: Whether to use Huber loss
            use_relative: Whether to use relative error
            
        Returns:
            Loss value
        """
        if use_huber:
            return self._huber_loss(predictions, targets, use_relative)
        else:
            return self.mse_loss(predictions, targets)
    
    def _huber_loss(self, 
                   predictions: torch.Tensor,
                   targets: torch.Tensor,
                   use_relative: bool = False) -> torch.Tensor:
        """Compute Huber loss with optional relative error."""
        if use_relative:
            # Relative error: |pred - target| / (|target| + epsilon)
            epsilon = 1e-8
            relative_error = torch.abs(predictions - targets) / (torch.abs(targets) + epsilon)
            return F.huber_loss(relative_error, torch.zeros_like(relative_error), delta=self.delta)
        else:
            return F.huber_loss(predictions, targets, delta=self.delta)

class StoichiometricConsistencyLoss(nn.Module):
    """
    Loss function to enforce stoichiometric consistency in metabolic networks.
    
    Features:
    - Mass balance constraints
    - Stoichiometric coefficient consistency
    - Reaction direction constraints
    """
    
    def __init__(self, 
                 stoichiometric_matrix: torch.Tensor,
                 reaction_directions: Optional[torch.Tensor] = None,
                 weight: float = 1.0):
        """
        Initialize the stoichiometric consistency loss.
        
        Args:
            stoichiometric_matrix: Stoichiometric coefficient matrix
            reaction_directions: Reaction direction indicators (1 for forward, -1 for reverse)
            weight: Loss weight
        """
        super(StoichiometricConsistencyLoss, self).__init__()
        
        self.stoichiometric_matrix = stoichiometric_matrix
        self.reaction_directions = reaction_directions
        self.weight = weight
        
        logger.info(f"Initialized StoichiometricConsistencyLoss with weight={weight}")
    
    def forward(self, 
                flux_predictions: torch.Tensor,
                metabolite_predictions: torch.Tensor) -> torch.Tensor:
        """
        Compute stoichiometric consistency loss.
        
        Args:
            flux_predictions: Predicted reaction fluxes
            metabolite_predictions: Predicted metabolite concentrations
            
        Returns:
            Loss value
        """
        # Mass balance constraint: S * v = 0
        mass_balance = torch.matmul(self.stoichiometric_matrix, flux_predictions)
        mass_balance_loss = F.mse_loss(mass_balance, torch.zeros_like(mass_balance))
        
        # Reaction direction constraints
        direction_loss = 0.0
        if self.reaction_directions is not None:
            # Enforce reaction directions
            direction_violations = torch.relu(-flux_predictions * self.reaction_directions)
            direction_loss = direction_violations.mean()
        
        # Metabolite concentration constraints (non-negative)
        concentration_loss = torch.relu(-metabolite_predictions).mean()
        
        total_loss = mass_balance_loss + direction_loss + concentration_loss
        
        return self.weight * total_loss

class MetabolicPathwayLoss(nn.Module):
    """
    Loss function for metabolic pathway analysis.
    
    Features:
    - Pathway coherence
    - Metabolic flux distribution
    - Pathway-specific constraints
    """
    
    def __init__(self, 
                 pathway_adjacency: torch.Tensor,
                 pathway_weights: Optional[torch.Tensor] = None,
                 coherence_weight: float = 1.0):
        """
        Initialize the metabolic pathway loss.
        
        Args:
            pathway_adjacency: Pathway adjacency matrix
            pathway_weights: Pathway importance weights
            coherence_weight: Pathway coherence weight
        """
        super(MetabolicPathwayLoss, self).__init__()
        
        self.pathway_adjacency = pathway_adjacency
        self.pathway_weights = pathway_weights
        self.coherence_weight = coherence_weight
        
        logger.info(f"Initialized MetabolicPathwayLoss with coherence_weight={coherence_weight}")
    
    def forward(self, 
                pathway_predictions: torch.Tensor,
                node_embeddings: torch.Tensor) -> torch.Tensor:
        """
        Compute metabolic pathway loss.
        
        Args:
            pathway_predictions: Predicted pathway assignments
            node_embeddings: Node embeddings
            
        Returns:
            Loss value
        """
        # Pathway coherence: similar nodes should have similar pathway assignments
        pathway_similarity = torch.matmul(pathway_predictions, pathway_predictions.t())
        embedding_similarity = F.cosine_similarity(node_embeddings.unsqueeze(1), 
                                                 node_embeddings.unsqueeze(0), dim=2)
        
        coherence_loss = F.mse_loss(pathway_similarity, embedding_similarity)
        
        # Pathway structure consistency
        structure_loss = 0.0
        if self.pathway_adjacency is not None:
            # Nodes connected in pathway should have similar pathway assignments
            structure_similarity = torch.matmul(self.pathway_adjacency, pathway_predictions)
            structure_loss = F.mse_loss(pathway_predictions, structure_similarity)
        
        # Pathway weight regularization
        weight_loss = 0.0
        if self.pathway_weights is not None:
            weight_loss = F.mse_loss(pathway_predictions, self.pathway_weights)
        
        total_loss = (self.coherence_weight * coherence_loss + 
                     structure_loss + weight_loss)
        
        return total_loss

class GrowthRateLoss(nn.Module):
    """
    Loss function for growth rate prediction.
    
    Features:
    - Growth rate constraints
    - Condition-specific growth rates
    - Biological plausibility
    """
    
    def __init__(self, 
                 min_growth_rate: float = 0.0,
                 max_growth_rate: float = 2.0,
                 condition_growth_rates: Optional[Dict[str, float]] = None):
        """
        Initialize the growth rate loss.
        
        Args:
            min_growth_rate: Minimum biologically plausible growth rate
            max_growth_rate: Maximum biologically plausible growth rate
            condition_growth_rates: Expected growth rates for different conditions
        """
        super(GrowthRateLoss, self).__init__()
        
        self.min_growth_rate = min_growth_rate
        self.max_growth_rate = max_growth_rate
        self.condition_growth_rates = condition_growth_rates or {
            'glucose': 0.85,
            'acetate': 0.42,
            'lactose': 0.35
        }
        
        logger.info(f"Initialized GrowthRateLoss with range [{min_growth_rate}, {max_growth_rate}]")
    
    def forward(self, 
                growth_predictions: torch.Tensor,
                condition_predictions: torch.Tensor,
                targets: torch.Tensor) -> torch.Tensor:
        """
        Compute growth rate loss.
        
        Args:
            growth_predictions: Predicted growth rates
            condition_predictions: Predicted growth conditions
            targets: Target growth rates
            
        Returns:
            Loss value
        """
        # Basic MSE loss
        mse_loss = F.mse_loss(growth_predictions, targets)
        
        # Growth rate bounds
        lower_bound_loss = torch.relu(self.min_growth_rate - growth_predictions).mean()
        upper_bound_loss = torch.relu(growth_predictions - self.max_growth_rate).mean()
        
        # Condition-specific growth rate constraints
        condition_loss = 0.0
        for i, condition in enumerate(['glucose', 'acetate', 'lactose']):
            condition_mask = (condition_predictions.argmax(dim=-1) == i)
            if condition_mask.sum() > 0:
                expected_rate = self.condition_growth_rates[condition]
                condition_predictions_masked = growth_predictions[condition_mask]
                condition_loss += F.mse_loss(condition_predictions_masked, 
                                           torch.full_like(condition_predictions_masked, expected_rate))
        
        total_loss = mse_loss + lower_bound_loss + upper_bound_loss + condition_loss
        
        return total_loss

class AttentionConsistencyLoss(nn.Module):
    """
    Loss function to enforce attention consistency.
    
    Features:
    - Attention weight regularization
    - Attention sparsity
    - Attention consistency across conditions
    """
    
    def __init__(self, 
                 sparsity_weight: float = 0.1,
                 consistency_weight: float = 0.5,
                 entropy_weight: float = 0.1):
        """
        Initialize the attention consistency loss.
        
        Args:
            sparsity_weight: Weight for attention sparsity
            consistency_weight: Weight for attention consistency
            entropy_weight: Weight for attention entropy
        """
        super(AttentionConsistencyLoss, self).__init__()
        
        self.sparsity_weight = sparsity_weight
        self.consistency_weight = consistency_weight
        self.entropy_weight = entropy_weight
        
        logger.info(f"Initialized AttentionConsistencyLoss with sparsity={sparsity_weight}, consistency={consistency_weight}")
    
    def forward(self, 
                attention_weights: torch.Tensor,
                condition_attention_weights: Optional[Dict[str, torch.Tensor]] = None) -> torch.Tensor:
        """
        Compute attention consistency loss.
        
        Args:
            attention_weights: Main attention weights
            condition_attention_weights: Attention weights for different conditions
            
        Returns:
            Loss value
        """
        # Attention sparsity (L1 regularization)
        sparsity_loss = torch.norm(attention_weights, p=1)
        
        # Attention entropy (encourage diverse attention)
        attention_probs = F.softmax(attention_weights, dim=-1)
        entropy = -torch.sum(attention_probs * torch.log(attention_probs + 1e-8), dim=-1)
        entropy_loss = -entropy.mean()  # Negative because we want to maximize entropy
        
        # Attention consistency across conditions
        consistency_loss = 0.0
        if condition_attention_weights is not None:
            condition_weights = list(condition_attention_weights.values())
            for i in range(len(condition_weights)):
                for j in range(i + 1, len(condition_weights)):
                    consistency_loss += F.mse_loss(condition_weights[i], condition_weights[j])
        
        total_loss = (self.sparsity_weight * sparsity_loss + 
                     self.entropy_weight * entropy_loss + 
                     self.consistency_weight * consistency_loss)
        
        return total_loss

class CombinedMetabolicLoss(nn.Module):
    """
    Combined loss function for all metabolic tasks.
    
    Features:
    - Multi-task loss combination
    - Dynamic loss weighting
    - Task-specific loss functions
    """
    
    def __init__(self, 
                 task_weights: Optional[Dict[str, float]] = None,
                 use_dynamic_weighting: bool = True):
        """
        Initialize the combined metabolic loss.
        
        Args:
            task_weights: Weights for different tasks
            use_dynamic_weighting: Whether to use dynamic weighting
        """
        super(CombinedMetabolicLoss, self).__init__()
        
        # Default task weights
        if task_weights is None:
            task_weights = {
                'node_classification': 1.0,
                'node_regression': 1.0,
                'graph_classification': 1.0,
                'condition_classification': 1.0,
                'pathway_analysis': 1.0,
                'graph_regression': 1.0,
                'stoichiometric_consistency': 0.5,
                'metabolic_pathway': 0.3,
                'growth_rate': 1.0,
                'attention_consistency': 0.1
            }
        
        self.task_weights = nn.Parameter(torch.tensor(list(task_weights.values())))
        self.task_names = list(task_weights.keys())
        self.use_dynamic_weighting = use_dynamic_weighting
        
        # Initialize individual loss functions
        self.node_classification_loss = MetabolicClassificationLoss()
        self.node_regression_loss = MetabolicRegressionLoss()
        self.graph_classification_loss = MetabolicClassificationLoss()
        self.condition_classification_loss = MetabolicClassificationLoss()
        self.pathway_analysis_loss = MetabolicClassificationLoss()
        self.graph_regression_loss = MetabolicRegressionLoss()
        self.stoichiometric_loss = StoichiometricConsistencyLoss(torch.randn(10, 10))
        self.metabolic_pathway_loss = MetabolicPathwayLoss(torch.randn(10, 10))
        self.growth_rate_loss = GrowthRateLoss()
        self.attention_consistency_loss = AttentionConsistencyLoss()
        
        logger.info(f"Initialized CombinedMetabolicLoss with {len(self.task_names)} tasks")
    
    def forward(self, 
                predictions: Dict[str, torch.Tensor],
                targets: Union[Dict[str, torch.Tensor], Any],  # Can be dict or MetabolicTargets
                additional_data: Optional[Dict[str, torch.Tensor]] = None) -> Dict[str, torch.Tensor]:
        """
        Compute combined metabolic loss.
        
        Args:
            predictions: Model predictions
            targets: Ground truth targets (dict or MetabolicTargets)
            additional_data: Additional data for specialized losses
            
        Returns:
            Dictionary containing individual and total losses
        """
        # Convert MetabolicTargets to dictionary if needed
        if hasattr(targets, 'node_classification'):  # It's a MetabolicTargets dataclass
            targets_dict = {}
            if targets.node_classification is not None:
                targets_dict['node_classification'] = targets.node_classification
            if targets.node_regression is not None:
                targets_dict['node_regression'] = targets.node_regression
            if targets.graph_classification is not None:
                targets_dict['graph_classification'] = targets.graph_classification
            if targets.graph_regression is not None:
                targets_dict['graph_regression'] = targets.graph_regression
            if targets.condition_prediction is not None:
                targets_dict['condition_prediction'] = targets.condition_prediction
            if targets.pathway_analysis is not None:
                targets_dict['pathway_analysis'] = targets.pathway_analysis
            if targets.growth_rate is not None:
                targets_dict['growth_rate'] = targets.growth_rate
        else:
            targets_dict = targets  # It's already a dictionary
        
        losses = {}
        total_loss = 0.0
        
        # Node classification loss
        if 'node_classification' in predictions and 'node_classification' in targets_dict:
            loss = self.node_classification_loss(
                predictions['node_classification'],
                targets['node_classification']
            )
            losses['node_classification'] = loss
            total_loss += self.task_weights[0] * loss
        
        # Node regression loss
        if 'node_regression' in predictions and 'node_regression' in targets_dict:
            loss = self.node_regression_loss(
                predictions['node_regression'],
                targets_dict['node_regression']
            )
            losses['node_regression'] = loss
            total_loss += self.task_weights[1] * loss
        
        # Graph classification loss
        if 'graph_classification' in predictions and 'graph_classification' in targets_dict:
            loss = self.graph_classification_loss(
                predictions['graph_classification'],
                targets_dict['graph_classification']
            )
            losses['graph_classification'] = loss
            total_loss += self.task_weights[2] * loss
        
        # Condition classification loss
        if 'condition_classification' in predictions and 'condition_classification' in targets_dict:
            loss = self.condition_classification_loss(
                predictions['condition_classification'],
                targets_dict['condition_classification']
            )
            losses['condition_classification'] = loss
            total_loss += self.task_weights[3] * loss
        
        # Pathway analysis loss
        if 'pathway_analysis' in predictions and 'pathway_analysis' in targets_dict:
            loss = self.pathway_analysis_loss(
                predictions['pathway_analysis'],
                targets_dict['pathway_analysis']
            )
            losses['pathway_analysis'] = loss
            total_loss += self.task_weights[4] * loss
        
        # Graph regression loss
        if 'graph_regression' in predictions and 'graph_regression' in targets_dict:
            loss = self.graph_regression_loss(
                predictions['graph_regression'],
                targets_dict['graph_regression']
            )
            losses['graph_regression'] = loss
            total_loss += self.task_weights[5] * loss
        
        # Specialized metabolic losses (if additional data provided)
        if additional_data is not None:
            # Stoichiometric consistency loss
            if 'flux_predictions' in predictions and 'metabolite_predictions' in predictions:
                loss = self.stoichiometric_loss(
                    predictions['flux_predictions'],
                    predictions['metabolite_predictions']
                )
                losses['stoichiometric_consistency'] = loss
                total_loss += self.task_weights[6] * loss
            
            # Metabolic pathway loss
            if 'pathway_predictions' in predictions and 'node_embeddings' in predictions:
                loss = self.metabolic_pathway_loss(
                    predictions['pathway_predictions'],
                    predictions['node_embeddings']
                )
                losses['metabolic_pathway'] = loss
                total_loss += self.task_weights[7] * loss
            
            # Growth rate loss
            if ('growth_predictions' in predictions and 
                'condition_predictions' in predictions and 
                'growth_rate' in targets_dict):
                loss = self.growth_rate_loss(
                    predictions['growth_predictions'],
                    predictions['condition_predictions'],
                    targets_dict['growth_rate']
                )
                losses['growth_rate'] = loss
                total_loss += self.task_weights[8] * loss
            
            # Attention consistency loss
            if 'attention_weights' in predictions:
                condition_attention = additional_data.get('condition_attention_weights')
                loss = self.attention_consistency_loss(
                    predictions['attention_weights'],
                    condition_attention
                )
                losses['attention_consistency'] = loss
                total_loss += self.task_weights[9] * loss
        
        losses['total_loss'] = total_loss
        losses['task_weights'] = self.task_weights
        
        return losses

def test_custom_loss_functions():
    """Test the custom loss functions."""
    logger.info("Testing custom loss functions...")
    
    # Test parameters
    batch_size = 32
    num_classes = 3
    num_nodes = 100
    
    # Create test data
    predictions = torch.randn(batch_size, num_classes)
    targets = torch.randint(0, num_classes, (batch_size,))
    regression_predictions = torch.randn(batch_size)
    regression_targets = torch.randn(batch_size)
    
    # Test MetabolicClassificationLoss
    logger.info("Testing MetabolicClassificationLoss...")
    classification_loss = MetabolicClassificationLoss()
    loss = classification_loss(predictions, targets)
    logger.info(f"Classification loss: {loss.item():.4f}")
    
    # Test MetabolicRegressionLoss
    logger.info("Testing MetabolicRegressionLoss...")
    regression_loss = MetabolicRegressionLoss()
    loss = regression_loss(regression_predictions, regression_targets)
    logger.info(f"Regression loss: {loss.item():.4f}")
    
    # Test StoichiometricConsistencyLoss
    logger.info("Testing StoichiometricConsistencyLoss...")
    stoichiometric_matrix = torch.randn(50, 100)
    stoichiometric_loss = StoichiometricConsistencyLoss(stoichiometric_matrix)
    flux_predictions = torch.randn(100)
    metabolite_predictions = torch.randn(50)
    loss = stoichiometric_loss(flux_predictions, metabolite_predictions)
    logger.info(f"Stoichiometric loss: {loss.item():.4f}")
    
    # Test MetabolicPathwayLoss
    logger.info("Testing MetabolicPathwayLoss...")
    pathway_adjacency = torch.randn(10, 10)
    pathway_loss = MetabolicPathwayLoss(pathway_adjacency)
    pathway_predictions = torch.randn(10, 5)
    node_embeddings = torch.randn(10, 64)
    loss = pathway_loss(pathway_predictions, node_embeddings)
    logger.info(f"Pathway loss: {loss.item():.4f}")
    
    # Test GrowthRateLoss
    logger.info("Testing GrowthRateLoss...")
    growth_loss = GrowthRateLoss()
    growth_predictions = torch.randn(batch_size)
    condition_predictions = torch.randn(batch_size, 3)
    growth_targets = torch.randn(batch_size)
    loss = growth_loss(growth_predictions, condition_predictions, growth_targets)
    logger.info(f"Growth rate loss: {loss.item():.4f}")
    
    # Test AttentionConsistencyLoss
    logger.info("Testing AttentionConsistencyLoss...")
    attention_loss = AttentionConsistencyLoss()
    attention_weights = torch.randn(num_nodes, num_nodes)
    loss = attention_loss(attention_weights)
    logger.info(f"Attention consistency loss: {loss.item():.4f}")
    
    # Test CombinedMetabolicLoss
    logger.info("Testing CombinedMetabolicLoss...")
    combined_loss = CombinedMetabolicLoss()
    
    # Create test predictions and targets
    test_predictions = {
        'node_classification': torch.randn(num_nodes, 2),
        'node_regression': torch.randn(num_nodes),
        'graph_classification': torch.randn(1, 3),
        'condition_classification': torch.randn(1, 3),
        'pathway_analysis': torch.randn(num_nodes, 5),
        'graph_regression': torch.randn(1)
    }
    
    test_targets = {
        'node_classification': torch.randint(0, 2, (num_nodes,)),
        'node_regression': torch.randn(num_nodes),
        'graph_classification': torch.randint(0, 3, (1,)),
        'condition_classification': torch.randint(0, 3, (1,)),
        'pathway_analysis': torch.randint(0, 5, (num_nodes,)),
        'graph_regression': torch.randn(1)
    }
    
    losses = combined_loss(test_predictions, test_targets)
    logger.info(f"Combined loss: {losses['total_loss'].item():.4f}")
    
    logger.info("All custom loss functions tested successfully!")

if __name__ == "__main__":
    test_custom_loss_functions() 