#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This programme is designed to obtain consensus sequences based on reference fasta and vcf files with a specified frequency threshold

The program implementation is based on three main steps:
    1) Obtaining the desired N specified fragments, L specified lengths
    2) Filtering the VCF file by allele frequency
    3) Obtaining consensus for the sample with applied mutations
"""

import argparse
import sys
import os
import pickle
import re
import random
import tempfile
from typing import List, Tuple, Dict, Union



def create_simple_index(fasta_path: str) -> Dict[str, Tuple[int, int]]:
    """
    Creates index: {chrom: (start_pos, length)}
    """
    index = {}
    with open(fasta_path, 'rb') as f:
        chrom = None
        start_pos = 0
        length = 0
        pos = 0
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                if chrom:
                    index[chrom] = (start_pos, length)
                break
            if line.startswith(b'>'):
                if chrom:
                    index[chrom] = (start_pos, length)
                chrom = line[1:].strip().decode().split()[0]
                start_pos = pos + len(line)
                length = 0
            else:
                length += len(line.strip())
    index_file = fasta_path + '.idx'
    with open(index_file, 'wb') as f:
        pickle.dump(index, f)
    return index


def load_or_create_index(fasta_path: str) -> Dict[str, Tuple[int, int]]:
    """
    Checks whether there is an index for fasta. If not, creates one
    """
    index_file = fasta_path + '.idx'
    if os.path.exists(index_file):
        fasta_mtime = os.path.getmtime(fasta_path)
        index_mtime = os.path.getmtime(index_file)
        if index_mtime >= fasta_mtime:
            with open(index_file, 'rb') as f:
                return pickle.load(f)
    return create_simple_index(fasta_path)


def get_sequence(fasta_path: str, index: Dict[str, Tuple[int, int]], chrom: str, start: int, end: int) -> Union[str, None]:
    """
    The function extracts a sequence from FASTA
    """
    if chrom not in index:
        return None
    start_pos, chrom_length = index[chrom]
    start = max(1, start)
    end = min(chrom_length, end)
    if start > end:
        return ""
    with open(fasta_path, 'rb') as f:
        f.seek(start_pos + start - 1)
        raw_bytes_to_read = (end - start + 1) * 2 
        seq_bytes = f.read(raw_bytes_to_read)
    cleaned_seq = seq_bytes.decode().replace('\n', '')
    result_seq = cleaned_seq[:end - start + 1]

    return result_seq.upper()


def get_random_fragments(fasta_path: str, fragment_length: int, num_fragments: int) -> Dict[str, str]:
    """
    Generates random fragments of a specified length from a FASTA file
    """
    index = load_or_create_index(fasta_path)
    fragments = {}
    chroms = [chrom for chrom, (_, length) in index.items() if length >= fragment_length]
    if not chroms:
        raise ValueError("There are no chroms long enough to extract fragments")
    chrom_weights = [index[chrom][1] for chrom in chroms]
    
    attempts = 0
    max_attempts = num_fragments * 42 
    
    while len(fragments) < num_fragments and attempts < max_attempts:
        attempts += 1
        chrom = random.choices(chroms, weights=chrom_weights, k=1)[0]
        _, chrom_length = index[chrom]
        if chrom_length >= fragment_length:
            start = random.randint(1, chrom_length - fragment_length + 1)
            end = start + fragment_length - 1
            seq = get_sequence(fasta_path, index, chrom, start, end)
            if seq and len(seq) == fragment_length:
                fragment_name = f"{chrom}:{start}-{end}"
                if fragment_name not in fragments:
                    fragments[fragment_name] = seq
    if len(fragments) < num_fragments:
         print(f"Warning: {num_fragments} fragments requested, {len(fragments)} generated. The FASTA file may be too small or contain many repetitive regions.", file=sys.stderr)
    return fragments


def calculate_internal_af(fields: List[str]) -> float:
    """
    The function for calculating internal frequency in vcf. The AD or GT values are used for the calculation
    """
    format_field = fields[8]
    samples = fields[9:]
    if not samples or not format_field:
        return 0.0
    format_fields = format_field.split(':')
    try:
        gt_index = format_fields.index('GT')
    except ValueError:
        return 0.0
    alt_alleles = fields[4].split(',') if fields[4] != '.' else []
    num_alt_alleles = len(alt_alleles)
    if num_alt_alleles == 0:
        return 0.0
    alt_allele_count = 0
    total_alleles = 0
    for sample_data in samples:
        if sample_data == '.' or sample_data == './.' or sample_data == '.|.':
            continue
        sample_fields = sample_data.split(':')
        if gt_index >= len(sample_fields):
            continue
        gt = sample_fields[gt_index]
        if gt == '.':
            continue
        alleles = []
        i = 0
        while i < len(gt):
            if gt[i].isdigit():
                num_start = i
                while i < len(gt) and gt[i].isdigit():
                    i += 1
                alleles.append(gt[num_start:i])
            else:
                i += 1
        for allele in alleles:
            if allele.isdigit():
                allele_num = int(allele)
                total_alleles += 1
                if allele_num > 0:
                    alt_allele_count += 1
    if total_alleles == 0:
        return 0.0
    return float(alt_allele_count) / float(total_alleles)


def filter_vcf_file(vcf_path: str, fragments: Dict[str, str], min_af: float, filtered_vcf_path: str) -> None:
    """
    Filters VCF files based on fragments and AF
    """
    valid_chromosomes = set(fragment_key.split(':')[0] for fragment_key in fragments.keys())
    input_file = open(vcf_path, 'r')
    output_file = open(filtered_vcf_path, 'w')
    vcf_chromosomes = set()
    try:
        for line in input_file:
            if line.startswith('#'):
                output_file.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom = fields[0]
            vcf_chromosomes.add(chrom)
            pos = int(fields[1])
            info_field = fields[7]
            if chrom not in valid_chromosomes:
                continue
            af_match = re.search(r'AF=([^;]+)', info_field)
            if af_match:
                af_value = float(af_match.group(1))
            else:
                af_value = calculate_internal_af(fields)
                if info_field == '.':
                    info_field = 'AF=' + str(af_value)
                else:
                    info_field = info_field + ';AF=' + str(af_value)
                fields[7] = info_field
            if af_value >= min_af:
                output_file.write('\t'.join(fields) + '\n')
    finally:
        input_file.close()
        output_file.close()
        common_chroms = valid_chromosomes.intersection(vcf_chromosomes)
        if vcf_chromosomes and valid_chromosomes and not common_chroms:
            error_message = "Ops! The VCF file does not match the FASTA file. Please check the files"
            print(f"Error: {error_message}")
            print(f"Chromosomes in FASTA: {sorted(valid_chromosomes)}")
            print(f"Chromosomes in VCF: {sorted(vcf_chromosomes)}")
            raise ValueError(error_message)


def get_sample_names_from_vcf(vcf_path: str) -> List[str]:
    """
    The function extracts sample names from VCF file
    """
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                if len(parts) > 9:
                    return parts[9:]
                else:
                    return []
    return []


def extract_allele_frequencies(info_field: str, num_alt_alleles: int) -> List[float]:
    """
    Extracts allele frequencies from the INFO field
    """
    frequencies = []
    info_parts = info_field.split(';')
    af_found = False
    for part in info_parts:
        if part.startswith('AF='):
            af_found = True
            af_value = part[3:]
            if af_value:
                freq_strings = af_value.split(',')
                try:
                    frequencies = [float(f) for f in freq_strings if f]
                except ValueError:
                    pass
            break
    if not af_found:            #Note: with the current logic of the code, this block will never be used because this function uses an already filtered VCF, which will always have AF due to the use of def calculate_internal_af. However, I thought it advisable to leave this block, realising that it will never be executed in the current logic of the programme. 
        ac_match = None
        an_value = None
        for part in info_parts:
            if part.startswith('AC='):
                ac_match = part[3:]
            elif part.startswith('AN='):
                try:
                    an_value = int(part[3:])
                except ValueError:
                    pass
        if ac_match and an_value:
            ac_values = ac_match.split(',')
            try:
                counts = [int(c) for c in ac_values if c]
                frequencies = [count / an_value for count in counts]
            except ValueError:
                pass
    if not frequencies:
        frequencies = [1.0 / (num_alt_alleles + 1)] * num_alt_alleles if num_alt_alleles > 0 else []
    while len(frequencies) < num_alt_alleles:
        frequencies.append(0.0)
    return frequencies


def parse_genotype_with_priority(sample_data: str, format_fields: list, ref: str,
                               alt_alleles: list, allele_frequencies: list) -> str:
    """
    Parse genotype of sample with priority selection of alleles
    """
    if sample_data == '.' or sample_data == './.' or sample_data == '.|.':
        return ref
    sample_fields = sample_data.split(':')
    gt_index = -1
    for i, field in enumerate(format_fields):
        if field == 'GT':
            gt_index = i
            break
    if gt_index == -1 or gt_index >= len(sample_fields):
        return ref
    gt = sample_fields[gt_index]
    if gt == '.':
        return ref
    alleles = []
    i = 0
    while i < len(gt):
        if gt[i].isdigit():
            num_start = i
            while i < len(gt) and gt[i].isdigit():
                i += 1
            alleles.append(gt[num_start:i])
        else:
            i += 1
    if not alleles:
        return ref
    if len(alleles) == 2 and alleles[0] == alleles[1]:
        allele_num = alleles[0]
        if allele_num != '0' and allele_num.isdigit() and int(allele_num) > 0:
            allele_index = int(allele_num) - 1
            if allele_index < len(alt_alleles):
                return alt_alleles[allele_index]
        elif allele_num == '0':
            return ref
    else:
        valid_alt_alleles = []
        valid_alt_frequencies = []
        for allele_num in alleles:
            if allele_num != '0' and allele_num.isdigit() and int(allele_num) > 0:
                allele_index = int(allele_num) - 1
                if allele_index < len(alt_alleles):
                    valid_alt_alleles.append(alt_alleles[allele_index])
                    if allele_index < len(allele_frequencies):
                        valid_alt_frequencies.append(allele_frequencies[allele_index])
                    else:
                        valid_alt_frequencies.append(0.0)
        if valid_alt_alleles and valid_alt_frequencies:
            max_freq_index = 0
            for i in range(1, len(valid_alt_frequencies)):
                if valid_alt_frequencies[i] > valid_alt_frequencies[max_freq_index]:
                    max_freq_index = i
            return valid_alt_alleles[max_freq_index]
        elif valid_alt_alleles:
            return valid_alt_alleles[0]
        else:
            return ref
    return ref


def parse_vcf_mutations_with_allele_freq(vcf_path: str) -> Dict[str, List[Tuple[int, str, List[str], List[str], List[float]]]]:
    """
    Parses VCF file and returns a dictionary of mutations
    mutations[chrom] = [(pos, ref, alt_alleles, sample_genotypes, frequencies), ...]
    """
    mutations = {}
    sample_names = get_sample_names_from_vcf(vcf_path)
    if not sample_names:
        sample_names = ["sample"]
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 + len(sample_names):
                continue
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_alleles = fields[4].split(',') if fields[4] != '.' else []
            info = fields[7]
            format_field = fields[8]
            samples = fields[9:]
            format_fields = format_field.split(':')
            frequencies = extract_allele_frequencies(info, len(alt_alleles))
            sample_genotypes = []     
            for sample_data in samples:
                sample_fields = sample_data.split(':')
                try:
                    gt_index = format_fields.index('GT')
                    if gt_index < len(sample_fields):
                        gt_str = sample_fields[gt_index]
                    else:
                        gt_str = '.'
                except ValueError:
                    gt_str = '.'
                sample_genotypes.append(gt_str)
            
            if chrom not in mutations:
                mutations[chrom] = []
            mutations[chrom].append((pos, ref, alt_alleles, sample_genotypes, frequencies))
    return mutations


def find_overlapping_mutations(mutations, chrom, start, end):
    """
    Finds mutations that overlap with target fragment
    """
    overlapping = []
    if chrom not in mutations:
        return overlapping
    for mutation in mutations[chrom]:
        pos, ref, _, _, _ = mutation
        mut_start = pos
        mut_end = pos + len(ref) - 1
        if start <= mut_end and end >= mut_start:
            overlapping.append(mutation)
    return overlapping



def apply_all_mutations_from_single_sample_indel(
    fragment_seq: str,
    fragment_start: int,
    fragment_end: int,
    mutations: List[Tuple[int, str, List[str], List[str], List[float]]],
    chosen_sample_idx: int 
) -> str: 
    """
    Applies all crossing mutations from a single sample with phase preservation
    """
 
    seq_list = list(fragment_seq.upper())
    sorted_mutations = sorted(mutations, key=lambda x: x[0])
    haplotype_choice_position = None 

    offset = 0

    for pos, ref_allele, alt_alleles, sample_genotypes, allele_frequencies in sorted_mutations:
        if chosen_sample_idx >= len(sample_genotypes):
            continue

        chosen_allele_seq = None

        gt_str = sample_genotypes[chosen_sample_idx]
        
        if haplotype_choice_position is None:
            if gt_str == '.' or ('/' not in gt_str and '|' not in gt_str):
                gt_alleles = ['.', '.']
            else:
                separator = '|' if '|' in gt_str else '/'
                gt_alleles = gt_str.split(separator)

            if len(gt_alleles) >= 2 and len(set(gt_alleles)) > 1 and '.' not in gt_alleles:
                if gt_alleles[1] != '0':
                    haplotype_choice_position = 1
                elif gt_alleles[0] != '0':
                    haplotype_choice_position = 0
                else:
                    haplotype_choice_position = 0
        
        if gt_str == '.' or ('/' not in gt_str and '|' not in gt_str):
            gt_alleles_current = ['.', '.']
        else:
            separator = '|' if '|' in gt_str else '/'
            gt_alleles_current = gt_str.split(separator)
            
        if not gt_alleles_current or len(gt_alleles_current) < 2:
             continue

        allele_num_to_use = None 


        if haplotype_choice_position is not None and haplotype_choice_position < len(gt_alleles_current):
            allele_num_to_use = gt_alleles_current[haplotype_choice_position]
        else:
            chosen_allele_seq = parse_genotype_with_priority(gt_str, ['GT'], ref_allele, alt_alleles, allele_frequencies)
            if not chosen_allele_seq or chosen_allele_seq == '.':
                continue
            chosen_allele_seq = chosen_allele_seq.upper()
            ref_allele = ref_allele.upper()
            if ref_allele == chosen_allele_seq:
                continue

        if chosen_allele_seq is None:
            if allele_num_to_use is None:
                 allele_num_to_use = gt_alleles_current[0] if gt_alleles_current else '0'
            
            if allele_num_to_use == '.':
                chosen_allele_seq = '.'
            elif not allele_num_to_use.isdigit():
                chosen_allele_seq = ref_allele
            else:
                idx = int(allele_num_to_use)
                if idx == 0:
                    chosen_allele_seq = ref_allele
                elif 1 <= idx <= len(alt_alleles):
                    chosen_allele_seq = alt_alleles[idx - 1]
                else:
                    chosen_allele_seq = ref_allele

        if not chosen_allele_seq or chosen_allele_seq == '.':
            continue
        chosen_allele_seq = chosen_allele_seq.upper()
        ref_allele = ref_allele.upper()
        if ref_allele == chosen_allele_seq:
            continue

        relative_pos = pos - fragment_start + offset

        ref_len = len(ref_allele)
        alt_len = len(chosen_allele_seq)

        if 0 <= relative_pos < len(seq_list) and ref_len > 0:
            seq_list[relative_pos:relative_pos + ref_len] = list(chosen_allele_seq)
            offset += alt_len - ref_len

    return ''.join(seq_list)


def create_mutated_fasta(fasta_path: str,
                         vcf_path: str,
                         output_fasta_path: str,
                         fragment_length: int,
                         num_fragments: int,
                         min_af: float):
    """
    Creates a FASTA file where each fragment contains mutations from the VCF file.
    This function connects the entire pipeline
    """

    print(f'[1/3] Generating {num_fragments} random fragments with a length of {fragment_length}...')
    fragments = get_random_fragments(fasta_path, fragment_length, num_fragments)
    print(f'{len(fragments)} fragments generated.')
    print(f'[2/3] Filtering VCF file {vcf_path} with threshold AF>={min_af}...')
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as tmp_vcf:
        output_path = tmp_vcf.name
    try:
        filter_vcf_file(vcf_path, fragments, min_af, output_path)
        print(f' VCF file filtered and saved to a temporary file.')
        print(f'[3/3] Creating a FASTA file with mutations {output_fasta_path}...')
        mutations = parse_vcf_mutations_with_allele_freq(output_path)
        sample_names = get_sample_names_from_vcf(output_path)
        
        if not sample_names:
            sample_names = ["sample"]

        with open(output_fasta_path, 'w') as output_file:
            for fragment_name, fragment_seq in fragments.items():
                chrom, start_end = fragment_name.split(':', 1)
                start, end = map(int, start_end.split('-'))
                overlapping_mutations = find_overlapping_mutations(mutations, chrom, start, end)
                chosen_sample_idx = random.randint(0, len(sample_names) - 1)
                chosen_sample_name = sample_names[chosen_sample_idx] if chosen_sample_idx < len(sample_names) else f"sample_{chosen_sample_idx}"
                
                if overlapping_mutations:
                    try:
                        mutated_seq = apply_all_mutations_from_single_sample_indel(
                            fragment_seq, start, end, overlapping_mutations, chosen_sample_idx
                        )
                    except Exception as e:
                        print(f"Error applying mutations to fragment {fragment_name}: {e}")
                        mutated_seq = fragment_seq
                else:
                    mutated_seq = fragment_seq

                header = f">{chrom}:{start}-{end}_{chosen_sample_name}"
                output_file.write(header + '\n')
                for i in range(0, len(mutated_seq), 42):
                    output_file.write(mutated_seq[i:i+42] + '\n')
        print(f" Done! FASTA file with mutations saved: {output_fasta_path}")

    finally:
        if os.path.exists(output_path):
            os.remove(output_path)

def main():
    parser = argparse.ArgumentParser(
        description="A program for obtaining consensus sequences for samples from a given VCF file based on a reference FASTA and a specified frequency threshold.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--fasta', '-f', required=True, help='The path to the original FASTA file.')
    parser.add_argument('--vcf', '-v', required=True, help='The path to the original VCF file.')
    parser.add_argument('--output', '-o', default='./consensus.fa' , help='The path to the output FASTA file with mutations.')
    parser.add_argument('--fragment_length', '-l', type=int, default=100,
                        help='Desired length of each fragment (default: 100)')
    parser.add_argument('--num_fragments', '-n', type=int, default=100,
                        help='Desired number of fragments to generate (default: 100)')
    parser.add_argument('--min_af', '-a', type=float, default=0.01,
                        help='Desired allele frequency threshold for VCF filtering (default: 0.01)')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Ops! FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.vcf):
        print(f"Ops! VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    if args.fragment_length <= 0:
        print(f"Ops! The fragment length must be a positive number. Received:{args.fragment_length}", file=sys.stderr)
        sys.exit(1)
    if args.num_fragments <= 0:
        print(f"Ops! The number of fragments must be a positive number. Received: {args.num_fragments}", file=sys.stderr)
        sys.exit(1)
    if not (0.0 <= args.min_af <= 1.0):
        print(f"Ops! The minimum allele frequency must be between 0.0 and 1.0. Received: {args.min_af}", file=sys.stderr)
        sys.exit(1)

    try:
        create_mutated_fasta(
            fasta_path=args.fasta,
            vcf_path=args.vcf,
            output_fasta_path=args.output,
            fragment_length=args.fragment_length,
            num_fragments=args.num_fragments,
            min_af=args.min_af
        )
        print("Done it!")
        print(f"Result saved in: {args.output}")
    except Exception as e:
        print(f"Runtime error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
