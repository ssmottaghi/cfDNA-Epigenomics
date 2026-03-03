import pandas as pd
import os

input_folder = "Data/raw_bed"
output_folder = "Data/processed_bed"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def data_preprocess(data):
    """
    Filters for high-score peaks and removes duplicates based on genomic start coordinates.
    """
    # Select specific columns of interest
    preprocessed_data = data[["chrom", "chromStart", "chromEnd", "score", "signalValue", "pValue"]]
    
    # Quality filter: Keep only peaks with a score greater than 100
    preprocessed_data = preprocessed_data[preprocessed_data["score"] > 100]
    
    # Sorting to ensure the highest pValue/score is kept during deduplication
    preprocessed_data = preprocessed_data.sort_values(by=['chrom', 'chromStart', "score", "pValue"])
    
    # Remove duplicates based on chromosome and start position
    preprocessed_data = preprocessed_data.drop_duplicates(subset=["chrom", "chromStart"], keep='last', ignore_index=True)
    
    return preprocessed_data

# Iterate through files and process
files = [f for f in os.listdir(input_folder) if f.endswith(".bed")]

for file in files:
    file_path = os.path.join(input_folder, file)
    
    # Read the narrowPeak/BED format (typically 10 columns)
    data = pd.read_table(file_path, header=0, names=[
        "chrom", "chromStart", "chromEnd", "name", "score",
        "strand", "signalValue", "pValue", "qValue", "peak"
    ])
    
    # Apply preprocessing
    cleaned_data = data_preprocess(data)
    
    # Save output with 'n.bed' suffix in the output folder
    output_filename = file[:-4] + "n.bed"
    cleaned_data.to_csv(os.path.join(output_folder, output_filename), 
                        header=False, index=False, sep='\t')

print(f"Successfully preprocessed {len(files)} files.")
