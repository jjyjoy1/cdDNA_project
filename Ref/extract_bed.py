import pandas as pd

# Load the CSV file
csv_file_path = "/home/jjy/Documents/GitHubRepos/cfDNA_cfRNA_project/config/cfDNA/cfDNA_amplicons.csv"
df = pd.read_csv(csv_file_path, sep=',')

# Display the first few rows to check the structure
df.head()

# Extract chromosome, start, and end positions from the "COORDINATES" column
bed_data = []

for coord in df["COORDINATES"].dropna():
    parts = coord.split(":")
    genome_version = parts[0]  # Should be "GRCH37"
    chrom = parts[1].split('chr')[-1]  # Chromosome number
    start, end = parts[2].split("-")  # Start and end positions
    strand = parts[3]
    strand_symbol = '-' if strand == '-1' else "+"
    # Convert to BED format (0-based start position)
    bed_data.append([f"chr{chrom}", int(start) - 1, int(end), strand_symbol])

# Create a DataFrame for BED format
bed_df = pd.DataFrame(bed_data, columns=["Chromosome", "Start", "End", "Strand"])
bed_df['Chromosome'] = bed_df['Chromosome'].str.replace('chr', '')
# Save the BED file
bed_file_path = "/home/jjy/Documents/GitHubRepos/cfDNA_cfRNA_project/Ref/cfDNA_amplicons.bed"
bed_df.to_csv(bed_file_path, sep="\t", index=False, header=False)




