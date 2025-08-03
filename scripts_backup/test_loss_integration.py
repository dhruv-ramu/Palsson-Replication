#!/usr/bin/env python3
"""
Test Loss Integration with MetabolicTargets

Simple test to verify that the loss function works correctly with the new data structures.
"""

import torch
import sys
sys.path.append('.')

from scripts.data_handling import MetabolicTargets, create_synthetic_data
from scripts.loss_functions import CombinedMetabolicLoss
from scripts.multi_task_model import MultiTaskMetabolicNetwork

def test_loss_integration():
    """Test that loss function works with MetabolicTargets."""
    print("Testing loss function integration...")
    
    # Create synthetic data
    graphs, targets = create_synthetic_data(num_samples=5)
    
    # Create model
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
    model = MultiTaskMetabolicNetwork(**model_config)
    
    # Create loss function
    loss_fn = CombinedMetabolicLoss()
    
    # Test with a single sample
    graph = graphs[0]
    target = targets[0]
    
    print(f"Graph: {graph.num_nodes} nodes, {graph.num_edges} edges")
    print(f"Target has: node_classification={target.node_classification is not None}, "
          f"node_regression={target.node_regression is not None}")
    
    # Forward pass
    predictions = model(graph.node_features, graph.edge_indices)
    print(f"Predictions keys: {list(predictions.keys())}")
    
    # Test loss computation
    try:
        loss = loss_fn(predictions, target)
        print(f"✅ Loss computation successful!")
        print(f"Loss keys: {list(loss.keys())}")
        print(f"Total loss: {loss['total_loss'].item():.4f}")
    except Exception as e:
        print(f"❌ Loss computation failed: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_loss_integration() 