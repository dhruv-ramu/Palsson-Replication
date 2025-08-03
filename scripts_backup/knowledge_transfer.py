#!/usr/bin/env python3
"""
Knowledge Transfer for Metabolic Network Embedding

This module implements knowledge transfer techniques for Phase 3 Week 1,
including cross-condition knowledge transfer, knowledge distillation,
and transfer learning strategies as specified in scope.md.

Features:
- Cross-condition knowledge transfer
- Knowledge distillation
- Teacher-student learning
- Feature alignment
- Progressive knowledge transfer

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
import copy

# Import our custom modules
import sys
sys.path.append('.')
from scripts.multi_task_model import MultiTaskMetabolicNetwork, MultiTaskTrainer
from scripts.loss_functions import CombinedMetabolicLoss

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class KnowledgeDistillationLoss(nn.Module):
    """
    Knowledge distillation loss for teacher-student learning.
    
    Features:
    - Soft target distillation
    - Hard target learning
    - Temperature scaling
    - Feature distillation
    """
    
    def __init__(self, 
                 temperature: float = 4.0,
                 alpha: float = 0.7,
                 beta: float = 0.3):
        """
        Initialize knowledge distillation loss.
        
        Args:
            temperature: Temperature for softmax scaling
            alpha: Weight for soft target loss
            beta: Weight for hard target loss
        """
        super(KnowledgeDistillationLoss, self).__init__()
        
        self.temperature = temperature
        self.alpha = alpha
        self.beta = beta
        
        # Standard cross-entropy loss for hard targets
        self.hard_loss = nn.CrossEntropyLoss()
        
        # KL divergence loss for soft targets
        self.soft_loss = nn.KLDivLoss(reduction='batchmean')
        
        logger.info(f"Initialized KnowledgeDistillationLoss with temperature={temperature}, alpha={alpha}, beta={beta}")
    
    def forward(self, 
                student_logits: torch.Tensor,
                teacher_logits: torch.Tensor,
                hard_targets: torch.Tensor) -> torch.Tensor:
        """
        Compute knowledge distillation loss.
        
        Args:
            student_logits: Student model logits
            teacher_logits: Teacher model logits
            hard_targets: Hard target labels
            
        Returns:
            Combined distillation loss
        """
        # Soft target loss (KL divergence)
        student_probs = F.log_softmax(student_logits / self.temperature, dim=1)
        teacher_probs = F.softmax(teacher_logits / self.temperature, dim=1)
        soft_loss = self.soft_loss(student_probs, teacher_probs) * (self.temperature ** 2)
        
        # Hard target loss (cross-entropy)
        hard_loss = self.hard_loss(student_logits, hard_targets)
        
        # Combined loss
        total_loss = self.alpha * soft_loss + self.beta * hard_loss
        
        return total_loss

class FeatureAlignmentLoss(nn.Module):
    """
    Feature alignment loss for knowledge transfer.
    
    Features:
    - Feature space alignment
    - L2 distance minimization
    - Cosine similarity maximization
    - Multi-scale feature alignment
    """
    
    def __init__(self, 
                 alignment_type: str = 'l2',
                 normalize_features: bool = True):
        """
        Initialize feature alignment loss.
        
        Args:
            alignment_type: Type of alignment ('l2', 'cosine', 'both')
            normalize_features: Whether to normalize features
        """
        super(FeatureAlignmentLoss, self).__init__()
        
        self.alignment_type = alignment_type
        self.normalize_features = normalize_features
        
        logger.info(f"Initialized FeatureAlignmentLoss with type={alignment_type}")
    
    def forward(self, 
                student_features: torch.Tensor,
                teacher_features: torch.Tensor) -> torch.Tensor:
        """
        Compute feature alignment loss.
        
        Args:
            student_features: Student model features
            teacher_features: Teacher model features
            
        Returns:
            Feature alignment loss
        """
        if self.normalize_features:
            student_features = F.normalize(student_features, p=2, dim=1)
            teacher_features = F.normalize(teacher_features, p=2, dim=1)
        
        if self.alignment_type == 'l2':
            # L2 distance minimization
            loss = F.mse_loss(student_features, teacher_features)
        
        elif self.alignment_type == 'cosine':
            # Cosine similarity maximization
            cosine_sim = F.cosine_similarity(student_features, teacher_features, dim=1)
            loss = 1 - torch.mean(cosine_sim)
        
        elif self.alignment_type == 'both':
            # Both L2 and cosine
            l2_loss = F.mse_loss(student_features, teacher_features)
            cosine_sim = F.cosine_similarity(student_features, teacher_features, dim=1)
            cosine_loss = 1 - torch.mean(cosine_sim)
            loss = l2_loss + cosine_loss
        
        else:
            raise ValueError(f"Unknown alignment type: {self.alignment_type}")
        
        return loss

class CrossConditionKnowledgeTransfer:
    """
    Cross-condition knowledge transfer framework.
    
    Features:
    - Cross-condition transfer
    - Progressive knowledge transfer
    - Multi-task knowledge transfer
    - Condition-specific adaptation
    """
    
    def __init__(self, 
                 model_config: Dict[str, Any],
                 output_dir: str = "results/metabolic_network/knowledge_transfer"):
        """
        Initialize cross-condition knowledge transfer.
        
        Args:
            model_config: Model configuration
            output_dir: Output directory
        """
        self.model_config = model_config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Initialize teacher and student models
        self.teacher_model = MultiTaskMetabolicNetwork(**model_config)
        self.student_model = MultiTaskMetabolicNetwork(**model_config)
        
        # Knowledge transfer history
        self.transfer_history = {
            'distillation': {},
            'feature_alignment': {},
            'progressive_transfer': {}
        }
        
        logger.info(f"Initialized CrossConditionKnowledgeTransfer with output directory: {self.output_dir}")
    
    def setup_teacher(self, teacher_model_path: str):
        """
        Set up teacher model from pre-trained weights.
        
        Args:
            teacher_model_path: Path to teacher model checkpoint
        """
        checkpoint = torch.load(teacher_model_path)
        self.teacher_model.load_state_dict(checkpoint['model_state_dict'])
        self.teacher_model.eval()  # Set to evaluation mode
        
        logger.info(f"Teacher model loaded from: {teacher_model_path}")
    
    def knowledge_distillation(self, 
                             target_data: Dict[str, torch.Tensor],
                             num_epochs: int = 100,
                             learning_rate: float = 0.001,
                             temperature: float = 4.0,
                             alpha: float = 0.7,
                             beta: float = 0.3) -> Dict[str, Any]:
        """
        Perform knowledge distillation from teacher to student.
        
        Args:
            target_data: Target condition data
            num_epochs: Number of training epochs
            learning_rate: Learning rate
            temperature: Distillation temperature
            alpha: Weight for soft target loss
            beta: Weight for hard target loss
            
        Returns:
            Distillation history
        """
        logger.info("Starting knowledge distillation...")
        
        # Initialize loss functions and optimizer
        distillation_loss = KnowledgeDistillationLoss(temperature, alpha, beta)
        optimizer = torch.optim.Adam(self.student_model.parameters(), lr=learning_rate)
        
        # Training history
        distillation_history = {
            'epoch_losses': [],
            'soft_losses': [],
            'hard_losses': [],
            'best_loss': float('inf'),
            'best_model_state': None
        }
        
        for epoch in range(num_epochs):
            # Forward pass
            student_output = self.student_model(target_data['node_features'], 
                                              target_data['edge_index'])
            
            with torch.no_grad():
                teacher_output = self.teacher_model(target_data['node_features'], 
                                                  target_data['edge_index'])
            
            # Compute distillation loss
            total_loss = distillation_loss(student_output['node_classification'],
                                         teacher_output['node_classification'],
                                         target_data['targets']['node_classification'])
            
            # Backward pass
            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()
            
            # Log progress
            if epoch % 10 == 0:
                logger.info(f"Epoch {epoch}/{num_epochs}: Distillation Loss: {total_loss.item():.4f}")
            
            # Store history
            distillation_history['epoch_losses'].append(total_loss.item())
            
            # Save best model
            if total_loss.item() < distillation_history['best_loss']:
                distillation_history['best_loss'] = total_loss.item()
                distillation_history['best_model_state'] = copy.deepcopy(self.student_model.state_dict())
        
        # Save distilled model
        checkpoint_path = self.output_dir / f"distilled_model_{self.timestamp}.pth"
        torch.save({
            'model_state_dict': distillation_history['best_model_state'],
            'model_config': self.model_config,
            'distillation_history': distillation_history,
            'temperature': temperature,
            'alpha': alpha,
            'beta': beta
        }, checkpoint_path)
        logger.info(f"Distilled model saved to: {checkpoint_path}")
        
        # Store in transfer history
        self.transfer_history['distillation'] = distillation_history
        
        return distillation_history
    
    def feature_alignment_transfer(self, 
                                 target_data: Dict[str, torch.Tensor],
                                 num_epochs: int = 100,
                                 learning_rate: float = 0.001,
                                 alignment_weight: float = 0.5,
                                 alignment_type: str = 'both') -> Dict[str, Any]:
        """
        Perform feature alignment-based knowledge transfer.
        
        Args:
            target_data: Target condition data
            num_epochs: Number of training epochs
            learning_rate: Learning rate
            alignment_weight: Weight for feature alignment loss
            alignment_type: Type of feature alignment
            
        Returns:
            Feature alignment history
        """
        logger.info("Starting feature alignment transfer...")
        
        # Initialize loss functions and optimizer
        task_loss = CombinedMetabolicLoss()
        alignment_loss = FeatureAlignmentLoss(alignment_type=alignment_type)
        optimizer = torch.optim.Adam(self.student_model.parameters(), lr=learning_rate)
        
        # Training history
        alignment_history = {
            'epoch_losses': [],
            'task_losses': [],
            'alignment_losses': [],
            'best_loss': float('inf'),
            'best_model_state': None
        }
        
        for epoch in range(num_epochs):
            # Forward pass
            student_output = self.student_model(target_data['node_features'], 
                                              target_data['edge_index'])
            
            with torch.no_grad():
                teacher_output = self.teacher_model(target_data['node_features'], 
                                                  target_data['edge_index'])
            
            # Compute losses
            task_loss_value = task_loss(student_output, target_data['targets'])
            alignment_loss_value = alignment_loss(student_output['node_embeddings'],
                                                teacher_output['node_embeddings'])
            
            total_loss = task_loss_value['total_loss'] + alignment_weight * alignment_loss_value
            
            # Backward pass
            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()
            
            # Log progress
            if epoch % 10 == 0:
                logger.info(f"Epoch {epoch}/{num_epochs}: "
                          f"Task Loss: {task_loss_value['total_loss']:.4f}, "
                          f"Alignment Loss: {alignment_loss_value:.4f}")
            
            # Store history
            alignment_history['epoch_losses'].append(total_loss.item())
            alignment_history['task_losses'].append(task_loss_value['total_loss'].item())
            alignment_history['alignment_losses'].append(alignment_loss_value.item())
            
            # Save best model
            if total_loss.item() < alignment_history['best_loss']:
                alignment_history['best_loss'] = total_loss.item()
                alignment_history['best_model_state'] = copy.deepcopy(self.student_model.state_dict())
        
        # Save aligned model
        checkpoint_path = self.output_dir / f"aligned_model_{self.timestamp}.pth"
        torch.save({
            'model_state_dict': alignment_history['best_model_state'],
            'model_config': self.model_config,
            'alignment_history': alignment_history,
            'alignment_weight': alignment_weight,
            'alignment_type': alignment_type
        }, checkpoint_path)
        logger.info(f"Aligned model saved to: {checkpoint_path}")
        
        # Store in transfer history
        self.transfer_history['feature_alignment'] = alignment_history
        
        return alignment_history
    
    def progressive_knowledge_transfer(self, 
                                    target_data: Dict[str, torch.Tensor],
                                    num_epochs_per_stage: int = 50,
                                    learning_rates: List[float] = [0.001, 0.0001, 0.00001],
                                    transfer_strategies: List[str] = ['distillation', 'alignment', 'combined']) -> Dict[str, Any]:
        """
        Perform progressive knowledge transfer through multiple stages.
        
        Args:
            target_data: Target condition data
            num_epochs_per_stage: Number of epochs per stage
            learning_rates: Learning rates for each stage
            transfer_strategies: Transfer strategies for each stage
            
        Returns:
            Progressive transfer history
        """
        logger.info("Starting progressive knowledge transfer...")
        
        progressive_history = {
            'stage_results': [],
            'overall_losses': [],
            'best_loss': float('inf'),
            'best_model_state': None
        }
        
        current_model = copy.deepcopy(self.student_model)
        
        for stage, (lr, strategy) in enumerate(zip(learning_rates, transfer_strategies)):
            logger.info(f"Progressive transfer stage {stage + 1}/{len(transfer_strategies)}: {strategy} with lr={lr}")
            
            # Initialize optimizer
            optimizer = torch.optim.Adam(current_model.parameters(), lr=lr)
            
            # Initialize loss functions based on strategy
            if strategy == 'distillation':
                distillation_loss = KnowledgeDistillationLoss()
                task_loss = CombinedMetabolicLoss()
            elif strategy == 'alignment':
                alignment_loss = FeatureAlignmentLoss()
                task_loss = CombinedMetabolicLoss()
            elif strategy == 'combined':
                distillation_loss = KnowledgeDistillationLoss()
                alignment_loss = FeatureAlignmentLoss()
                task_loss = CombinedMetabolicLoss()
            
            stage_losses = []
            
            for epoch in range(num_epochs_per_stage):
                # Forward pass
                student_output = current_model(target_data['node_features'], 
                                            target_data['edge_index'])
                
                with torch.no_grad():
                    teacher_output = self.teacher_model(target_data['node_features'], 
                                                      target_data['edge_index'])
                
                # Compute loss based on strategy
                if strategy == 'distillation':
                    total_loss = distillation_loss(student_output['node_classification'],
                                                 teacher_output['node_classification'],
                                                 target_data['targets']['node_classification'])
                
                elif strategy == 'alignment':
                    task_loss_value = task_loss(student_output, target_data['targets'])
                    alignment_loss_value = alignment_loss(student_output['node_embeddings'],
                                                        teacher_output['node_embeddings'])
                    total_loss = task_loss_value['total_loss'] + 0.5 * alignment_loss_value
                
                elif strategy == 'combined':
                    distillation_loss_value = distillation_loss(student_output['node_classification'],
                                                             teacher_output['node_classification'],
                                                             target_data['targets']['node_classification'])
                    task_loss_value = task_loss(student_output, target_data['targets'])
                    alignment_loss_value = alignment_loss(student_output['node_embeddings'],
                                                        teacher_output['node_embeddings'])
                    total_loss = (0.4 * distillation_loss_value + 
                                0.4 * task_loss_value['total_loss'] + 
                                0.2 * alignment_loss_value)
                
                # Backward pass
                optimizer.zero_grad()
                total_loss.backward()
                optimizer.step()
                
                stage_losses.append(total_loss.item())
                
                if epoch % 10 == 0:
                    logger.info(f"Stage {stage + 1}, Epoch {epoch}/{num_epochs_per_stage}: Loss: {total_loss.item():.4f}")
            
            # Store stage results
            stage_result = {
                'strategy': strategy,
                'learning_rate': lr,
                'final_loss': stage_losses[-1],
                'losses': stage_losses
            }
            progressive_history['stage_results'].append(stage_result)
            progressive_history['overall_losses'].extend(stage_losses)
            
            # Update best model
            if stage_losses[-1] < progressive_history['best_loss']:
                progressive_history['best_loss'] = stage_losses[-1]
                progressive_history['best_model_state'] = copy.deepcopy(current_model.state_dict())
        
        # Save progressive model
        checkpoint_path = self.output_dir / f"progressive_model_{self.timestamp}.pth"
        torch.save({
            'model_state_dict': progressive_history['best_model_state'],
            'model_config': self.model_config,
            'progressive_history': progressive_history,
            'transfer_strategies': transfer_strategies,
            'learning_rates': learning_rates
        }, checkpoint_path)
        logger.info(f"Progressive model saved to: {checkpoint_path}")
        
        # Store in transfer history
        self.transfer_history['progressive_transfer'] = progressive_history
        
        return progressive_history
    
    def compare_transfer_methods(self) -> Dict[str, Any]:
        """Compare different knowledge transfer methods."""
        logger.info("Comparing knowledge transfer methods...")
        
        comparison = {
            'distillation': self.transfer_history.get('distillation', {}),
            'feature_alignment': self.transfer_history.get('feature_alignment', {}),
            'progressive_transfer': self.transfer_history.get('progressive_transfer', {})
        }
        
        # Calculate performance metrics
        performance_summary = {}
        
        if comparison['distillation']:
            performance_summary['distillation'] = {
                'best_loss': comparison['distillation'].get('best_loss', float('inf')),
                'final_loss': comparison['distillation'].get('epoch_losses', [float('inf')])[-1] if comparison['distillation'].get('epoch_losses') else float('inf')
            }
        
        if comparison['feature_alignment']:
            performance_summary['feature_alignment'] = {
                'best_loss': comparison['feature_alignment'].get('best_loss', float('inf')),
                'final_loss': comparison['feature_alignment'].get('epoch_losses', [float('inf')])[-1] if comparison['feature_alignment'].get('epoch_losses') else float('inf')
            }
        
        if comparison['progressive_transfer']:
            performance_summary['progressive_transfer'] = {
                'best_loss': comparison['progressive_transfer'].get('best_loss', float('inf')),
                'final_loss': comparison['progressive_transfer'].get('overall_losses', [float('inf')])[-1] if comparison['progressive_transfer'].get('overall_losses') else float('inf')
            }
        
        # Save comparison results
        comparison_path = self.output_dir / f"transfer_methods_comparison_{self.timestamp}.json"
        with open(comparison_path, 'w') as f:
            json.dump(performance_summary, f, indent=2)
        
        logger.info(f"Transfer methods comparison saved to: {comparison_path}")
        logger.info(f"Performance summary: {performance_summary}")
        
        return performance_summary
    
    def visualize_knowledge_transfer(self):
        """Create visualizations for knowledge transfer results."""
        logger.info("Creating knowledge transfer visualizations...")
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Distillation loss curve
        if self.transfer_history.get('distillation', {}).get('epoch_losses'):
            axes[0, 0].plot(self.transfer_history['distillation']['epoch_losses'])
            axes[0, 0].set_title('Knowledge Distillation Loss')
            axes[0, 0].set_xlabel('Epoch')
            axes[0, 0].set_ylabel('Loss')
            axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Feature alignment loss curves
        if self.transfer_history.get('feature_alignment'):
            alignment_history = self.transfer_history['feature_alignment']
            if alignment_history.get('task_losses'):
                axes[0, 1].plot(alignment_history['task_losses'], label='Task Loss')
            if alignment_history.get('alignment_losses'):
                axes[0, 1].plot(alignment_history['alignment_losses'], label='Alignment Loss')
            axes[0, 1].set_title('Feature Alignment Losses')
            axes[0, 1].set_xlabel('Epoch')
            axes[0, 1].set_ylabel('Loss')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Progressive transfer stages
        if self.transfer_history.get('progressive_transfer', {}).get('stage_results'):
            stage_results = self.transfer_history['progressive_transfer']['stage_results']
            strategies = [result['strategy'] for result in stage_results]
            final_losses = [result['final_loss'] for result in stage_results]
            
            axes[1, 0].bar(strategies, final_losses)
            axes[1, 0].set_title('Progressive Transfer Stage Results')
            axes[1, 0].set_ylabel('Final Loss')
            axes[1, 0].tick_params(axis='x', rotation=45)
        
        # 4. Method comparison
        performance_data = self.compare_transfer_methods()
        if performance_data:
            methods = list(performance_data.keys())
            best_losses = [performance_data[method]['best_loss'] for method in methods]
            
            axes[1, 1].bar(methods, best_losses)
            axes[1, 1].set_title('Transfer Methods Comparison')
            axes[1, 1].set_ylabel('Best Loss')
            axes[1, 1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # Save visualization
        viz_path = self.output_dir / f"knowledge_transfer_visualization_{self.timestamp}.png"
        plt.savefig(viz_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Knowledge transfer visualization saved to: {viz_path}")

def test_knowledge_transfer():
    """Test the knowledge transfer modules."""
    logger.info("Testing knowledge transfer modules...")
    
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
    
    # Initialize knowledge transfer framework
    knowledge_transfer = CrossConditionKnowledgeTransfer(model_config)
    
    # Create test data
    target_data = {
        'node_features': torch.randn(100, 35),
        'edge_index': torch.randint(0, 100, (2, 200)),
        'targets': {
            'node_classification': torch.randint(0, 2, (100,)),
            'node_regression': torch.randn(100)
        }
    }
    
    # Test knowledge distillation
    logger.info("Testing knowledge distillation...")
    distillation_history = knowledge_transfer.knowledge_distillation(
        target_data, num_epochs=5, learning_rate=0.001
    )
    
    # Test feature alignment transfer
    logger.info("Testing feature alignment transfer...")
    alignment_history = knowledge_transfer.feature_alignment_transfer(
        target_data, num_epochs=5, learning_rate=0.001
    )
    
    # Test progressive knowledge transfer
    logger.info("Testing progressive knowledge transfer...")
    progressive_history = knowledge_transfer.progressive_knowledge_transfer(
        target_data, 
        num_epochs_per_stage=3,
        learning_rates=[0.001, 0.0001],
        transfer_strategies=['distillation', 'alignment']
    )
    
    # Compare transfer methods
    logger.info("Comparing transfer methods...")
    performance_comparison = knowledge_transfer.compare_transfer_methods()
    
    # Create visualizations
    logger.info("Creating visualizations...")
    knowledge_transfer.visualize_knowledge_transfer()
    
    logger.info("Knowledge transfer modules tested successfully!")

if __name__ == "__main__":
    test_knowledge_transfer() 