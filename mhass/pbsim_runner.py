#!/usr/bin/env python3

import subprocess
import re
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def extract_passes(template_file):
    """Extract the number of passes from the template filename."""
    match = re.search(r"_np(\d+)", template_file.name)
    return int(match.group(1)) if match else None

def run_chunk(chunk_id, fasta_chunk, output_base, pbsim_path, qshmm_path):
    """Run PBSIM on a chunk of template files."""
    results = []
    for template_file in fasta_chunk:
        np_value = extract_passes(template_file)
        if np_value is None:
            print(f"Warning: Could not extract np value from {template_file.name}")
            continue
            
        template = template_file.stem
        out_dir = output_base / template
        out_dir.mkdir(parents=True, exist_ok=True)

        # Run PBSIM
        pbsim_cmd = [
            str(pbsim_path),
            "--strategy", "templ",
            "--template", str(template_file),
            "--method", "qshmm",
            "--qshmm", str(qshmm_path),
            "--accuracy-mean", "0.85",
            "--pass-num", str(np_value),
            "--prefix", str(out_dir / "sim")
        ]
        
        try:
            # Capture the output instead of discarding it
            result = subprocess.run(pbsim_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, text=True)
            if result.returncode != 0:
                print(f"Error running PBSIM on {template_file}")
                print(f"Command: {' '.join(pbsim_cmd)}")
                print(f"Return code: {result.returncode}")
                print(f"Error output: {result.stderr}")
                continue
        except Exception as e:
            print(f"Exception when running PBSIM on {template_file}: {str(e)}")
            continue
            
        # Run CCS
        ccs_cmd = ['ccs', '-j', "1", str(out_dir / "sim.bam"), str(out_dir / "ccs.fastq")]
        try:
            subprocess.run(ccs_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError:
            print(f"Error running CCS on {out_dir / 'sim.bam'}")
            continue
            
        results.append(template_file.name)
    
    return results

def chunked(lst, n_chunks):
    """Split a list into approximately equal chunks."""
    if n_chunks <= 0:
        return [lst]
    
    avg = max(1, len(lst) // n_chunks)
    chunks = [lst[i * avg: (i + 1) * avg] for i in range(n_chunks)]
    
    # Handle any remaining items
    remainder = lst[n_chunks * avg:]
    for i, extra in enumerate(remainder):
        chunks[i % len(chunks)].append(extra)
    
    # Remove empty chunks
    return [chunk for chunk in chunks if chunk]

def run_pbsim(template_dir, output_base, pbsim_path, qshmm_path, threads):
    """Run PBSIM on template files in parallel."""
    template_dir = Path(template_dir)
    output_base = Path(output_base)
    output_base.mkdir(parents=True, exist_ok=True)
    
    # Ensure PBSIM and QSHMM paths are valid
    pbsim_path = Path(pbsim_path)
    qshmm_path = Path(qshmm_path)
    
    # If pbsim is not found at the specified path, try to find it in the PATH
    if not pbsim_path.exists():
        pbsim_in_path = shutil.which("pbsim")
        if pbsim_in_path:
            pbsim_path = Path(pbsim_in_path)
            print(f"Using pbsim from PATH: {pbsim_path}")
        else:
            print(f"Warning: PBSIM not found at {pbsim_path}. Simulation may fail.")
    
    # If QSHMM is not found, look for it in standard locations
    if not qshmm_path.exists():
        # Try common locations
        common_locations = [
            Path("/usr/local/share/pbsim3/data/QSHMM-RSII.model"),
            Path("/usr/share/pbsim3/data/QSHMM-RSII.model"),
            Path("./pbsim3/data/QSHMM-RSII.model")
        ]
        
        for loc in common_locations:
            if loc.exists():
                qshmm_path = loc
                print(f"Using QSHMM model from: {qshmm_path}")
                break
        else:
            print(f"Warning: QSHMM model not found at {qshmm_path}. Simulation may fail.")
    
    # Find all template FASTA files
    fastas = sorted(template_dir.glob("*.fasta"))
    if not fastas:
        print(f"Error: No template FASTA files found in {template_dir}")
        return None
    
    # Create chunks for parallel processing
    num_chunks = min(len(fastas), threads * 5)  # Create more chunks than threads for better load balancing
    chunks = chunked(fastas, num_chunks)
    print(f"Total template files: {len(fastas)}, split into {len(chunks)} chunks across {threads} threads")
    
    # Process chunks in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(run_chunk, i, chunk, output_base, pbsim_path, qshmm_path): i
            for i, chunk in enumerate(chunks)
        }
        
        total_processed = 0
        with tqdm(total=len(fastas), desc="Simulating reads") as pbar:
            for future in as_completed(futures):
                result = future.result()
                total_processed += len(result)
                pbar.update(len(result))
    
    print(f"Successfully processed {total_processed}/{len(fastas)} template files")
    return output_base
