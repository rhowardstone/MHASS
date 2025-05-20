#!/usr/bin/env python3

import os
from pathlib import Path
from tqdm import tqdm

def combine_fastqs(input_dir, output_fastq):
    """Combine FASTQ files and add template prefix to headers."""
    input_dir = Path(input_dir)
    total_records = 0
    total_files = 0
    
    # First, count total files to process for progress display
    subdirs = [d for d in input_dir.iterdir() if d.is_dir()]
    total_subdirs = len(subdirs)
    
    print(f"Found {total_subdirs} template directories to process")
    
    with open(output_fastq, 'w') as outfile:
        for subdir in tqdm(sorted(subdirs), desc="Combining FASTQ files"):
            template_name = subdir.name
            fastq_file = subdir / "ccs.fastq"
            
            if not fastq_file.exists():
                print(f"[Warning] Missing FASTQ: {fastq_file}")
                continue
            
            total_files += 1
            records_in_file = 0
            
            with open(fastq_file, 'r') as f:
                while True:
                    header = f.readline().strip()
                    if not header:
                        break  # end of file
                    
                    seq = f.readline().strip()
                    plus = f.readline().strip()
                    qual = f.readline().strip()
                    
                    # Rewrite header to prefix with template name
                    if header.startswith('@'):
                        new_header = f"@{template_name}/{header[1:]}"
                    else:
                        new_header = f"{template_name}/{header}"
                    
                    # Write to output with proper FASTQ format
                    outfile.write(f"{new_header}\n{seq}\n{plus}\n{qual}\n")
                    
                    records_in_file += 1
                    total_records += 1
    
    print(f"[Done] Combined {total_records} records from {total_files} FASTQ files into: {output_fastq}")
    return output_fastq
