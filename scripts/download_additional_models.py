#!/usr/bin/env python3
"""
Download Additional E. coli Models for Comparison
"""

import requests
import gzip
import os

def download_model(model_name):
    """Download a model from BiGG database."""
    print(f"Downloading {model_name}...")
    
    # Create bigg_models directory if it doesn't exist
    os.makedirs('bigg_models', exist_ok=True)
    
    # Download JSON format
    json_url = f"http://bigg.ucsd.edu/static/models/{model_name}.json"
    json_path = f"bigg_models/{model_name}.json"
    
    try:
        response = requests.get(json_url)
        if response.status_code == 200:
            with open(json_path, 'wb') as f:
                f.write(response.content)
            print(f"  ✅ {model_name}.json downloaded")
            return True
        else:
            print(f"  ❌ {model_name}.json not found (status: {response.status_code})")
            return False
    except Exception as e:
        print(f"  ❌ Error downloading {model_name}.json: {e}")
        return False

def main():
    """Download additional E. coli models."""
    print("=" * 60)
    print("DOWNLOADING ADDITIONAL E. COLI MODELS")
    print("=" * 60)
    
    # List of models to download
    models = [
        'iAF1260',  # Alternative model
        'iJR904',   # Another alternative
        'iML1515',  # Latest model
        'iEC1372'   # Another alternative
    ]
    
    success_count = 0
    
    for model in models:
        if download_model(model):
            success_count += 1
    
    print(f"\n✅ Downloaded {success_count}/{len(models)} additional models")
    print("Models are ready for comprehensive comparison!")

if __name__ == "__main__":
    main() 