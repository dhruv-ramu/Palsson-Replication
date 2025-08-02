import os
import gzip
import shutil
import urllib.request

# Create output directory
os.makedirs("bigg_models", exist_ok=True)

# Define URLs and output paths
downloads = {
    "iJO1366.xml.gz": "http://bigg.ucsd.edu/static/models/iJO1366.xml.gz",
    "iJO1366.json.gz": "http://bigg.ucsd.edu/static/models/iJO1366.json.gz"
}

for filename_gz, url in downloads.items():
    output_gz = os.path.join("bigg_models", filename_gz)
    output_file = output_gz.replace(".gz", "")

    # Download
    print(f"Downloading {filename_gz}...")
    urllib.request.urlretrieve(url, output_gz)

    # Decompress
    print(f"Decompressing {filename_gz}...")
    with gzip.open(output_gz, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    print(f"Saved to: {output_file}")

print("\n✅ All files downloaded and extracted successfully.")
