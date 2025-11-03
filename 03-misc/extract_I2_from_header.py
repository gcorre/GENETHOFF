#!/usr/bin/env python3
"""
Extract BC and UMI from reads header if it exists
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
import sys


def extract_index_from_header(header):
    """
    Extract index after '+' in header in R2 FASTQ file.
    
    Args:
        header: Description part of the FASTQ header
    
    Returns:
        Index if found
    """
    idx=header.find('+')
    if idx != -1:
        index = header[idx + 1:]
        return index
    else:
        print("No index found. Cannot continue processing this sample")
        sys.exit(1)


def process_fastq_gz(input_file, output_file):
    """
    Take R2 fastq file and extract BC2 and UMI sequences.
    
    Args:
        input_file: path to input FASTQ.gz
        output_file: path to output FASTQ.gz
    
    Returns:
        Number of treated reads
    """
    nb_with_index = 0
    nb_without_index = 0
    total_reads = 0
    
    # Compte d'abord le nombre total de reads
    print(f"Counting reads in {input_file}...")
    with gzip.open(input_file, 'rt') as handle:
        total_reads = sum(1 for _ in handle) // 4
    
    print(f"Total: {total_reads} reads to process\n")
    print(f"Writting output file...")
    
    # Traite et écrit les reads au fur et à mesure
    with gzip.open(input_file, 'rt') as input_handle, \
         gzip.open(output_file, 'wt') as output_handle:
        
        for i, record in enumerate(SeqIO.parse(input_handle, "fastq"), 1):
            # Affiche la progression tous les 1000 reads
            if i % 1000 == 0:
                percentage = (i / total_reads) * 100
                print(f"  Treated: {i}/{total_reads} reads ({percentage:.1f}%)...", end='\r')
            
            # Extrait l'index du header
            index_seq = extract_index_from_header(record.description)


            if index_seq:
                # Crée un nouveau record avec l'index comme séquence
                new_record = SeqRecord(
                    Seq(index_seq),
                    id=record.id,
                    name=record.name,
                    description=record.description,
                    letter_annotations={'phred_quality': [40] * len(index_seq)}
                )
                SeqIO.write(new_record, output_handle, "fastq")
                nb_with_index += 1
            else:
                # Écrit le record original si pas d'index trouvé
                SeqIO.write(record, output_handle, "fastq")
                nb_without_index += 1
    
    print(f"  Treated: {total_reads} reads (100.0%) - Completed !  ")
    
    print(f"\n Extraction completed:")
    print(f"  - {nb_with_index} reads with index")
    print(f"  - {nb_without_index} reads without index")
    print(f"  - Total: {total_reads} reads")
    
    return total_reads


if __name__ == "__main__":
    # Vérification des arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.fastq.gz> <output.fastq.gz>")
        print("\nExemple:")
        print("  python script.py input.fastq.gz output.fastq.gz")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Traitement du fichier
    try:
        process_fastq_gz(input_file, output_file)
    except FileNotFoundError:
        print(f"Erreur: File {input_file} doesn't exists.")
        sys.exit(1)
    except Exception as e:
        print(f"Error while processing: {e}")
        sys.exit(1)