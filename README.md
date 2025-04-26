# hmm_bigwigs

A Python tool to perform Hidden Markov Model (HMM) analysis on genomic data stored in bigWig files. This tool uses the [pomegranate](https://github.com/jmschrei/pomegranate) library to identify genomic states based on signal intensity patterns.

## Installation

Install using pip:

```bash
pip install hmm_bigwigs
```

## Requirements

- numpy
- pandas
- bioframe
- pomegranate
- pybbi

## Usage

```bash
bigwig_hmm.py -i <input.bw> -g <genome> -n <num_states> -o <output_prefix>
```

### Arguments

- `-i, --inputfile` (required): Path to input bigWig file
- `-g, --genome`: Reference genome (e.g., hg38, mm10)
- `-n, --num_states`: Number of HMM states to identify (default: 2)
- `-o, --outputfile`: Output file prefix (will generate files ending in "_state_HMM.bed")
- `--cmap`: Matplotlib colormap name for state visualization (default: 'coolwarm')

### Output Files

The tool generates several BED format files:
- `{prefix}_n_state_HMM_colored.bed`: Complete segmentation with all states
- `{prefix}_0_state.bed`: Regions assigned to state 0
- `{prefix}_1_state.bed`: Regions assigned to state 1
- `{prefix}_2_state.bed`: Regions assigned to state 2 (if n_states > 2)

Each BED file contains the following columns:
- chromosome
- start position
- end position
- state assignment
- score (0)
- strand (.)
- RGB color values

## Example

```bash
bigwig_hmm.py -i H3K27me3.bw -g hg38 -n 3 -o H3K27me3_domains
```

This command will:
1. Read the H3K27me3.bw bigWig file
2. Segment the genome into 3 states
3. Generate color-coded BED files for visualization

## Features

- Automatic handling of chromosome selection (excludes chrM, chrX, chrY by default)
- Intelligent state assignment based on signal intensity
- Color-coded output for easy visualization in genome browsers
- Efficient merging of adjacent regions with identical states

## License

[Add your license information here]

## Citation

If you use this tool in your research, please cite:
[Add citation information here]

## Contact

- Author: George Spracklin
- Issues: https://github.com/gspracklin/hmm_bigwigs/issues
