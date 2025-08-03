#!/usr/bin/env python3
"""
Domain Adaptation for Metabolic Network Embedding

This module implements domain adaptation techniques for Phase 3 Week 1,
including domain adaptation methods, adversarial training, and
domain-invariant feature learning as specified in scope.md.

Features:
- Domain adaptation techniques
- Adversarial training
- Domain-invariant feature learning
- Maximum Mean Discrepancy (MMD)
- Domain confusion loss

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

# Import our custom modules
import sys
sys.path.append('.')
from scripts.multi_task_model import MultiTaskMetabolicNetwork, MultiTaskTrainer
from scripts.loss_functions import CombinedMetabolicLoss

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DomainDiscriminator(nn.Module):
    """
    Domain discriminator for adversarial domain adaptation.
    
    Features:
    - Domain classification
    - Gradient reversal
    - Adversarial training
    """
    
    def __init__(self, 
                 input_dim: int,
                 hidden_dim: int = 128,
                 num_domains: int = 2,
                 dropout: float = 0.2):
        """
        Initialize the domain discriminator.
        
        Args:
            input_dim: Input feature dimension
            hidden_dim: Hidden layer dimension
            num_domains: Number of domains to discriminate
            dropout: Dropout rate
        """
        super(DomainDiscriminator, self).__init__()
        
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.num_domains = num_domains
        self.dropout = dropout
        
        # Domain discriminator network
        self.domain_classifier = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, num_domains)
        )
        
        logger.info(f"Initialized DomainDiscriminator with {num_domains} domains")
    
    def forward(self, features: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through domain discriminator.
        
        Args:
            features: Input features [batch_size, input_dim]
            
        Returns:
            Domain predictions [batch_size, num_domains]
        """
        return self.domain_classifier(features)

class GradientReversalLayer(nn.Module):
    """
    Gradient reversal layer for adversarial training.
    
    This layer reverses the gradient during backpropagation,
    allowing adversarial training of the feature extractor.
    """
    
    def __init__(self, alpha: float = 1.0):
        """
        Initialize the gradient reversal layer.
        
        Args:
            alpha: Gradient reversal coefficient
        """
        super(GradientReversalLayer, self).__init__()
        self.alpha = alpha
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass (identity function).
        
        Args:
            x: Input tensor
            
        Returns:
            Same tensor (identity function)
        """
        return x
    
    def backward(self, grad_output: torch.Tensor) -> torch.Tensor:
        """
        Backward pass with gradient reversal.
        
        Args:
            grad_output: Gradient from next layer
            
        Returns:
            Reversed gradient
        """
        return -self.alpha * grad_output

class AdversarialDomainAdaptation(nn.Module):
    """
    Adversarial domain adaptation with gradient reversal.
    
    Features:
    - Adversarial training
    - Gradient reversal
    - Domain-invariant features
    """
    
    def __init__(self, 
                 feature_extractor: nn.Module,
                 task_classifier: nn.Module,
                 domain_discriminator: DomainDiscriminator,
                 alpha: float = 1.0):
        """
        Initialize adversarial domain adaptation.
        
        Args:
            feature_extractor: Feature extraction network
            task_classifier: Task-specific classifier
            domain_discriminator: Domain discriminator
            alpha: Gradient reversal coefficient
        """
        super(AdversarialDomainAdaptation, self).__init__()
        
        self.feature_extractor = feature_extractor
        self.task_classifier = task_classifier
        self.domain_discriminator = domain_discriminator
        self.gradient_reversal = GradientReversalLayer(alpha)
        
        logger.info("Initialized AdversarialDomainAdaptation")
    
    def forward(self, 
                x: torch.Tensor,
                edge_index: torch.Tensor,
                batch: Optional[torch.Tensor] = None) -> Dict[str, torch.Tensor]:
        """
        Forward pass through adversarial domain adaptation.
        
        Args:
            x: Input features
            edge_index: Edge indices
            batch: Batch assignment (optional)
            
        Returns:
            Dictionary containing task predictions and domain predictions
        """
        # Extract features
        features = self.feature_extractor(x, edge_index, batch)
        
        # Task predictions
        task_predictions = self.task_classifier(features)
        
        # Domain predictions (with gradient reversal)
        reversed_features = self.gradient_reversal(features)
        domain_predictions = self.domain_discriminator(reversed_features)
        
        return {
            'task_predictions': task_predictions,
            'domain_predictions': domain_predictions,
            'features': features
        }

class MaximumMeanDiscrepancy(nn.Module):
    """
    Maximum Mean Discrepancy (MMD) for domain adaptation.
    
    Features:
    - MMD computation
    - Kernel-based distance
    - Domain-invariant learning
    """
    
    def __init__(self, 
                 kernel_mul: float = 2.0,
                 kernel_num: int = 5,
                 fix_sigma: Optional[float] = None):
        """
        Initialize MMD.
        
        Args:
            kernel_mul: Kernel multiplier
            kernel_num: Number of kernels
            fix_sigma: Fixed sigma value (optional)
        """
        super(MaximumMeanDiscrepancy, self).__init__()
        
        self.kernel_mul = kernel_mul
        self.kernel_num = kernel_num
        self.fix_sigma = fix_sigma
        
        logger.info(f"Initialized MMD with {kernel_num} kernels")
    
    def forward(self, 
                source_features: torch.Tensor,
                target_features: torch.Tensor) -> torch.Tensor:
        """
        Compute MMD between source and target features.
        
        Args:
            source_features: Source domain features
            target_features: Target domain features
            
        Returns:
            MMD loss
        """
        batch_size = source_features.size(0)
        
        # Compute kernel matrices
        kernels = self._guassian_kernel(source_features, target_features)
        
        # Compute MMD
        loss = 0
        for kernel in kernels:
            s1, s2 = kernel[:batch_size, :batch_size], kernel[batch_size:, batch_size:]
            s12 = kernel[:batch_size, batch_size:]
            loss += torch.mean(s1) + torch.mean(s2) - 2 * torch.mean(s12)
        
        return loss
    
    def _guassian_kernel(self, 
                        source_features: torch.Tensor,
                        target_features: torch.Tensor) -> List[torch.Tensor]:
        """
        Compute Gaussian kernel matrices.
        
        Args:
            source_features: Source domain features
            target_features: Target domain features
            
        Returns:
            List of kernel matrices
        """
        n_samples = source_features.size(0) + target_features.size(0)
        total_features = torch.cat([source_features, target_features], dim=0)
        
        total_features0 = total_features.unsqueeze(0).expand(total_features.size(0), total_features.size(0), total_features.size(1))
        total_features1 = total_features.unsqueeze(1).expand(total_features.size(0), total_features.size(0), total_features.size(1))
        
        L2_distance = ((total_features0 - total_features1) ** 2).sum(2)
        
        if self.fix_sigma:
            bandwidth = self.fix_sigma
        else:
            bandwidth = torch.sum(L2_distance.data) / (n_samples ** 2 - n_samples)
        
        bandwidth /= self.kernel_mul ** (self.kernel_num // 2)
        bandwidth_list = [bandwidth * (self.kernel_mul ** i) for i in range(self.kernel_num)]
        
        kernel_val = [torch.exp(-L2_distance / bandwidth_temp) for bandwidth_temp in bandwidth_list]
        
        return kernel_val

class DomainConfusionLoss(nn.Module):
    """
    Domain confusion loss for domain-invariant feature learning.
    
    Features:
    - Domain confusion
    - Entropy maximization
    - Domain-invariant features
    """
    
    def __init__(self, 
                 num_domains: int = 2,
                 temperature: float = 1.0):
        """
        Initialize domain confusion loss.
        
        Args:
            num_domains: Number of domains
            temperature: Temperature for softmax
        """
        super(DomainConfusionLoss, self).__init__()
        
        self.num_domains = num_domains
        self.temperature = temperature
        
        logger.info(f"Initialized DomainConfusionLoss with {num_domains} domains")
    
    def forward(self, domain_predictions: torch.Tensor) -> torch.Tensor:
        """
        Compute domain confusion loss.
        
        Args:
            domain_predictions: Domain predictions [batch_size, num_domains]
            
        Returns:
            Domain confusion loss
        """
        # Apply temperature scaling
        scaled_predictions = domain_predictions / self.temperature
        
        # Compute softmax probabilities
        domain_probs = F.softmax(scaled_predictions, dim=1)
        
        # Compute entropy (higher entropy = more confusion = better)
        entropy = -torch.sum(domain_probs * torch.log(domain_probs + 1e-8), dim=1)
        
        # Return negative entropy (we want to maximize entropy)
        return -torch.mean(entropy)

class DomainAdaptationTrainer:
    """
    Trainer for domain adaptation.
    
    Features:
    - Adversarial training
    - MMD-based training
    - Domain confusion training
    - Multi-objective optimization
    """
    
    def __init__(self, 
                 model: nn.Module,
                 domain_discriminator: DomainDiscriminator,
                 output_dir: str = "results/metabolic_network/domain_adaptation"):
        """
        Initialize domain adaptation trainer.
        
        Args:
            model: Main model
            domain_discriminator: Domain discriminator
            output_dir: Output directory
        """
        self.model = model
        self.domain_discriminator = domain_discriminator
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Initialize loss functions
        self.task_loss = CombinedMetabolicLoss()
        self.domain_loss = nn.CrossEntropyLoss()
        self.mmd_loss = MaximumMeanDiscrepancy()
        self.confusion_loss = DomainConfusionLoss()
        
        # Training history
        self.training_history = {
            'task_losses': [],
            'domain_losses': [],
            'mmd_losses': [],
            'confusion_losses': [],
            'total_losses': []
        }
        
        logger.info(f"Initialized DomainAdaptationTrainer with output directory: {self.output_dir}")
    
    def train_adversarial(self, 
                         source_data: Dict[str, torch.Tensor],
                         target_data: Dict[str, torch.Tensor],
                         num_epochs: int = 100,
                         learning_rate: float = 0.001,
                         lambda_adversarial: float = 0.1) -> Dict[str, Any]:
        """
        Train with adversarial domain adaptation.
        
        Args:
            source_data: Source domain data
            target_data: Target domain data
            num_epochs: Number of training epochs
            learning_rate: Learning rate
            lambda_adversarial: Weight for adversarial loss
            
        Returns:
            Training history
        """
        logger.info("Starting adversarial domain adaptation training...")
        
        # Initialize optimizers
        model_optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        discriminator_optimizer = torch.optim.Adam(self.domain_discriminator.parameters(), lr=learning_rate)
        
        # Training loop
        for epoch in range(num_epochs):
            # Train domain discriminator
            discriminator_loss = self._train_discriminator(source_data, target_data, discriminator_optimizer)
            
            # Train main model (with adversarial loss)
            model_losses = self._train_model_adversarial(source_data, target_data, model_optimizer, lambda_adversarial)
            
            # Log progress
            if epoch % 10 == 0:
                logger.info(f"Epoch {epoch}/{num_epochs}: "
                          f"Task Loss: {model_losses['task_loss']:.4f}, "
                          f"Domain Loss: {discriminator_loss:.4f}, "
                          f"Adversarial Loss: {model_losses['adversarial_loss']:.4f}")
            
            # Store history
            self.training_history['task_losses'].append(model_losses['task_loss'])
            self.training_history['domain_losses'].append(discriminator_loss)
            self.training_history['total_losses'].append(model_losses['total_loss'])
        
        return self.training_history
    
    def train_mmd(self, 
                  source_data: Dict[str, torch.Tensor],
                  target_data: Dict[str, torch.Tensor],
                  num_epochs: int = 100,
                  learning_rate: float = 0.001,
                  lambda_mmd: float = 0.1) -> Dict[str, Any]:
        """
        Train with MMD-based domain adaptation.
        
        Args:
            source_data: Source domain data
            target_data: Target domain data
            num_epochs: Number of training epochs
            learning_rate: Learning rate
            lambda_mmd: Weight for MMD loss
            
        Returns:
            Training history
        """
        logger.info("Starting MMD-based domain adaptation training...")
        
        # Initialize optimizer
        optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        
        # Training loop
        for epoch in range(num_epochs):
            # Forward pass
            source_features = self.model.get_embeddings(source_data['node_features'], 
                                                       source_data['edge_index'])
            target_features = self.model.get_embeddings(target_data['node_features'], 
                                                       target_data['edge_index'])
            
            # Compute losses
            task_loss = self.task_loss(self.model(source_data['node_features'], 
                                                 source_data['edge_index']), 
                                     source_data['targets'])
            mmd_loss = self.mmd_loss(source_features, target_features)
            
            total_loss = task_loss['total_loss'] + lambda_mmd * mmd_loss
            
            # Backward pass
            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()
            
            # Log progress
            if epoch % 10 == 0:
                logger.info(f"Epoch {epoch}/{num_epochs}: "
                          f"Task Loss: {task_loss['total_loss']:.4f}, "
                          f"MMD Loss: {mmd_loss:.4f}")
            
            # Store history
            self.training_history['task_losses'].append(task_loss['total_loss'].item())
            self.training_history['mmd_losses'].append(mmd_loss.item())
            self.training_history['total_losses'].append(total_loss.item())
        
        return self.training_history
    
    def train_confusion(self, 
                       source_data: Dict[str, torch.Tensor],
                       target_data: Dict[str, torch.Tensor],
                       num_epochs: int = 100,
                       learning_rate: float = 0.001,
                       lambda_confusion: float = 0.1) -> Dict[str, Any]:
        """
        Train with domain confusion loss.
        
        Args:
            source_data: Source domain data
            target_data: Target domain data
            num_epochs: Number of training epochs
            learning_rate: Learning rate
            lambda_confusion: Weight for confusion loss
            
        Returns:
            Training history
        """
        logger.info("Starting domain confusion training...")
        
        # Initialize optimizer
        optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        
        # Training loop
        for epoch in range(num_epochs):
            # Forward pass
            source_output = self.model(source_data['node_features'], 
                                     source_data['edge_index'])
            target_output = self.model(target_data['node_features'], 
                                     target_data['edge_index'])
            
            # Compute losses
            task_loss = self.task_loss(source_output, source_data['targets'])
            
            # Get domain predictions from domain discriminator
            source_features = self.model.get_embeddings(source_data['node_features'], 
                                                       source_data['edge_index'])
            target_features = self.model.get_embeddings(target_data['node_features'], 
                                                       target_data['edge_index'])
            
            source_domain_pred = self.domain_discriminator(source_features)
            target_domain_pred = self.domain_discriminator(target_features)
            
            # Combine domain predictions for confusion loss
            combined_domain_pred = torch.cat([source_domain_pred, target_domain_pred], dim=0)
            confusion_loss = self.confusion_loss(combined_domain_pred)
            
            total_loss = task_loss['total_loss'] + lambda_confusion * confusion_loss
            
            # Backward pass
            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()
            
            # Log progress
            if epoch % 10 == 0:
                logger.info(f"Epoch {epoch}/{num_epochs}: "
                          f"Task Loss: {task_loss['total_loss']:.4f}, "
                          f"Confusion Loss: {confusion_loss:.4f}")
            
            # Store history
            self.training_history['task_losses'].append(task_loss['total_loss'].item())
            self.training_history['confusion_losses'].append(confusion_loss.item())
            self.training_history['total_losses'].append(total_loss.item())
        
        return self.training_history
    
    def _train_discriminator(self, 
                           source_data: Dict[str, torch.Tensor],
                           target_data: Dict[str, torch.Tensor],
                           optimizer: torch.optim.Optimizer) -> float:
        """Train domain discriminator."""
        self.domain_discriminator.train()
        
        # Source domain (label 0)
        source_features = self.model.get_embeddings(source_data['node_features'], 
                                                   source_data['edge_index'])
        source_domain_pred = self.domain_discriminator(source_features)
        source_domain_loss = self.domain_loss(source_domain_pred, 
                                            torch.zeros(source_features.size(0), dtype=torch.long))
        
        # Target domain (label 1)
        target_features = self.model.get_embeddings(target_data['node_features'], 
                                                   target_data['edge_index'])
        target_domain_pred = self.domain_discriminator(target_features)
        target_domain_loss = self.domain_loss(target_domain_pred, 
                                            torch.ones(target_features.size(0), dtype=torch.long))
        
        # Total discriminator loss
        discriminator_loss = source_domain_loss + target_domain_loss
        
        # Backward pass
        optimizer.zero_grad()
        discriminator_loss.backward()
        optimizer.step()
        
        return discriminator_loss.item()
    
    def _train_model_adversarial(self, 
                               source_data: Dict[str, torch.Tensor],
                               target_data: Dict[str, torch.Tensor],
                               optimizer: torch.optim.Optimizer,
                               lambda_adversarial: float) -> Dict[str, float]:
        """Train main model with adversarial loss."""
        self.model.train()
        
        # Task loss on source domain
        source_output = self.model(source_data['node_features'], 
                                 source_data['edge_index'])
        task_loss = self.task_loss(source_output, source_data['targets'])
        
        # Adversarial loss (domain discriminator should be confused)
        source_features = self.model.get_embeddings(source_data['node_features'], 
                                                   source_data['edge_index'])
        target_features = self.model.get_embeddings(target_data['node_features'], 
                                                   target_data['edge_index'])
        
        source_domain_pred = self.domain_discriminator(source_features)
        target_domain_pred = self.domain_discriminator(target_features)
        
        # We want the discriminator to be confused (uniform predictions)
        uniform_target = torch.ones_like(source_domain_pred) / source_domain_pred.size(1)
        adversarial_loss = (F.mse_loss(F.softmax(source_domain_pred, dim=1), uniform_target) +
                           F.mse_loss(F.softmax(target_domain_pred, dim=1), uniform_target))
        
        total_loss = task_loss['total_loss'] + lambda_adversarial * adversarial_loss
        
        # Backward pass
        optimizer.zero_grad()
        total_loss.backward()
        optimizer.step()
        
        return {
            'task_loss': task_loss['total_loss'].item(),
            'adversarial_loss': adversarial_loss.item(),
            'total_loss': total_loss.item()
        }
    
    def visualize_domain_adaptation(self):
        """Create visualizations for domain adaptation results."""
        logger.info("Creating domain adaptation visualizations...")
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Task loss over time
        if self.training_history['task_losses']:
            axes[0, 0].plot(self.training_history['task_losses'])
            axes[0, 0].set_title('Task Loss Over Time')
            axes[0, 0].set_xlabel('Epoch')
            axes[0, 0].set_ylabel('Task Loss')
            axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Domain loss over time
        if self.training_history['domain_losses']:
            axes[0, 1].plot(self.training_history['domain_losses'])
            axes[0, 1].set_title('Domain Loss Over Time')
            axes[0, 1].set_xlabel('Epoch')
            axes[0, 1].set_ylabel('Domain Loss')
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. MMD loss over time
        if self.training_history['mmd_losses']:
            axes[1, 0].plot(self.training_history['mmd_losses'])
            axes[1, 0].set_title('MMD Loss Over Time')
            axes[1, 0].set_xlabel('Epoch')
            axes[1, 0].set_ylabel('MMD Loss')
            axes[1, 0].grid(True, alpha=0.3)
        
        # 4. Total loss over time
        if self.training_history['total_losses']:
            axes[1, 1].plot(self.training_history['total_losses'])
            axes[1, 1].set_title('Total Loss Over Time')
            axes[1, 1].set_xlabel('Epoch')
            axes[1, 1].set_ylabel('Total Loss')
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save visualization
        viz_path = self.output_dir / f"domain_adaptation_visualization_{self.timestamp}.png"
        plt.savefig(viz_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Domain adaptation visualization saved to: {viz_path}")

def test_domain_adaptation():
    """Test the domain adaptation modules."""
    logger.info("Testing domain adaptation modules...")
    
    # Model configuration
    model_config = {
        'input_dim': 35,
        'hidden_dim': 128,
        'embedding_dim': 64,
        'num_heads': 8,
        'num_layers': 3,
        'dropout': 0.2,
        'num_conditions': 3,
        'num_pathways': 10
    }
    
    # Initialize model and domain discriminator
    model = MultiTaskMetabolicNetwork(**model_config)
    domain_discriminator = DomainDiscriminator(input_dim=128, num_domains=2)
    
    # Test MMD
    logger.info("Testing MMD...")
    mmd = MaximumMeanDiscrepancy()
    source_features = torch.randn(50, 128)
    target_features = torch.randn(50, 128)
    mmd_loss = mmd(source_features, target_features)
    logger.info(f"MMD loss: {mmd_loss.item():.4f}")
    
    # Test domain confusion loss
    logger.info("Testing domain confusion loss...")
    confusion_loss = DomainConfusionLoss(num_domains=2)
    domain_predictions = torch.randn(100, 2)
    confusion_loss_value = confusion_loss(domain_predictions)
    logger.info(f"Domain confusion loss: {confusion_loss_value.item():.4f}")
    
    # Test domain adaptation trainer
    logger.info("Testing domain adaptation trainer...")
    trainer = DomainAdaptationTrainer(model, domain_discriminator)
    
    # Create test data
    source_data = {
        'node_features': torch.randn(100, 35),
        'edge_index': torch.randint(0, 100, (2, 200)),
        'targets': {
            'node_classification': torch.randint(0, 2, (100,)),
            'node_regression': torch.randn(100)
        }
    }
    
    target_data = {
        'node_features': torch.randn(80, 35),
        'edge_index': torch.randint(0, 80, (2, 150)),
        'targets': {
            'node_classification': torch.randint(0, 2, (80,)),
            'node_regression': torch.randn(80)
        }
    }
    
    # Test MMD training
    logger.info("Testing MMD training...")
    mmd_history = trainer.train_mmd(source_data, target_data, num_epochs=5)
    
    # Test confusion training
    logger.info("Testing confusion training...")
    confusion_history = trainer.train_confusion(source_data, target_data, num_epochs=5)
    
    # Create visualizations
    logger.info("Creating visualizations...")
    trainer.visualize_domain_adaptation()
    
    logger.info("Domain adaptation modules tested successfully!")

if __name__ == "__main__":
    test_domain_adaptation() 