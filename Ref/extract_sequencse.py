import sys
from Bio import SeqIO
from Bio.Seq import Seq

def parse_coordinates(coordinate_str):
    # Parse a line like "GRCH37:14:105246510-105246605:-1"
    parts = coordinate_str.strip().split(':')
    chrom = parts[1]
    start_end = parts[2].split('-')
    start = int(start_end[0]) - 1  # BED is 0-based
    end = int(start_end[1])
    strand = '-' if parts[3] == '-1' else '+'
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'name': coordinate_str.strip()
    }

def extract_sequences(coordinates_file, reference_fasta, output_fasta):
    # Load the reference genome
    genome = SeqIO.index(reference_fasta, 'fasta')
    
    # Parse coordinates
    with open(coordinates_file, 'r') as f:
        lines = f.readlines()[1:]  # Skip header line
    
    # Process each coordinate
    with open(output_fasta, 'w') as out:
        for line in lines:
            if not line.strip():
                continue
            region = parse_coordinates(line)
            chrom_seq = genome[region['chrom']].seq
            
            # Extract sequence
            seq = chrom_seq[region['start']:region['end']]
            
            # Reverse complement if strand is '-'
            if region['strand'] == '-':
                seq = seq.reverse_complement()
            
            # Write to FASTA
            out.write(f">{region['name']}\n{seq}\n")



if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python extract_sequences.py <coordinates.txt> <reference.fasta> <output.fasta>")
        sys.exit(1)
    
    coordinates_file = sys.argv[1]
    reference_fasta = sys.argv[2]
    output_fasta = sys.argv[3]
    
    extract_sequences(coordinates_file, reference_fasta, output_fasta)
