# GCF program for your multi-vcf file.

CLI program for generating FASTA files with mutations taken from a VCF file.

This script combines three main steps:

1.  **Random Fragment Generation:** Extracts random fragments of a specified length from the reference FASTA file.
2.  **VCF File Filtering:** Selects variants from the VCF file that overlap with the generated fragments and have an allele frequency (AF) above a given threshold.
3.  **Mutated FASTA Creation:** For each generated fragment, a single sample is randomly chosen from the filtered VCF file, and that sample's mutations are applied to the fragment's sequence. The result is written to a new FASTA file.

## Features

*   **Phase Preservation (Haplotype):** When applying mutations for the selected sample, the script attempts to preserve phase. It finds the first suitable heterozygous mutation and selects one "side" of the genotype (e.g., for `0/1`, either `0` or `1` is chosen). The same "side" is used for all subsequent mutations in that fragment from the same sample.
*   **Fallback Allele Selection Mechanism:** If the initial mutations are homozygous or have missing data (`.`), a fallback method (`parse_genotype_with_priority`) is used, which selects an allele based on frequency or priority.
*   **FASTA Indexing:** An index file (`<fasta_file>.idx`) is created and used for faster sequence access.
*   **Allele Frequency Handling:** If the `AF` field is missing in the VCF, the frequency is calculated based on sample genotypes (`GT`).

## Dependencies

*   **Python 3.6 or higher**
*  **No third-party libraries are required.** Only standard modules are used: `argparse`, `sys`, `os`, `pickle`, `re`, `random`, `tempfile`, `typing`

## Installation

1.  Ensure Python 3.6 or higher is installed.
2.  Copy the `get_consensus.py` file to the desired directory.

## Usage

```bash
python get_consensus.py --fasta FASTA --vcf VCF [--output OUTPUT] [--fragment_length FRAGMENT_LENGTH] [--num_fragments NUM_FRAGMENTS] [--min_af MIN_AF]
```
**Notes:** The program uses default parameters for --fragment_length (-l = 100), --num_fragments (-n = 100), --min_af (-a = 0.01).
## Graphical visualization of the program's algorithm
<pre>
+--------------------------------------------------+
|                 Input Data                       |
|             FASTA file (reference)               |
|             VCF file (variants)                  |
|             Command-line Parameters              |
|  (fragment length, number of fragments, min AF)  |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
|                   main()                         |
|       Parses command-line arguments              |
|       Validates input data and parameters        |
|       Calls create_mutated_fasta()               |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
|         create_mutated_fasta()                   |
|  Main function orchestrating the pipeline:       |
|        1. Generate Fragments                     |
|        2. Filter VCF                             |
|        3. Create Mutated FASTA                   |
+--------------------------------------------------+
         |           |                    |
         |           |                    +--------------------------------------+
         |           |                                                           |
         v           v                                                           v
+-----------------------------+    +---------------------------------+     +-----------------------------+
| 1. get_random_fragments()   |    | 2. filter_vcf_file()            |     | 3. parse_vcf_mutations_...()|
|                             |    |                                 |     |                             |
|    load_or_create_index()   |--> |   calculate_internal_af()       |     | get_sample_names_from_vcf() |
|    create_simple_index()    |    |   extract_allele_frequencies()  |     | extract_allele_frequencies()|
|    get_sequence()           |    |                                 |     |                             |
+-----------------------------+    +---------------------------------+     +-----------------------------+
         |                                   |                                       |
         |                                   |                                       |
         v                                   |                                       v
+-----------------------------+              |                     +---------------------------------+
|   load_or_create_index()    |              |                     |   find_overlapping_mutations()  |
|   create_simple_index()     |              |                     |                                 |
|   get_sequence()            |              |                     |   (in loop for each fragment)   |
+-----------------------------+              |                     +---------------------------------+
                                             |                                       |
                                             |                                       v
                                             |              +---------------------------------------+
                                             |              | (in loop for each fragment)           |
                                             |              | random.randint() (sample selection)   |
                                             |              | apply_all_mutations_from_...()        |
                                             |              | parse_genotype_with_priority()        |
                                             |              |  (fallback allele determination       |
                                             |              |      if phasing logic fails)          |
                                             |              +---------------------------------------+
                                             |                                       |
                                             |                                       |
                                             +---------------------------------------+
                                             |
                                             v
                              +---------------------------------+
                              | Write header and modified       |
                              | sequence to output FASTA file   |
                              +---------------------------------+
                                             |
                                             v
                              +---------------------------------+
                              |    Output Mutated FASTA File    |
                              +---------------------------------+
  </pre>
