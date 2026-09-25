# genetic-variant-interpreter
A Python tool for detecting SNPs and analyzing their effects on RNA and protein sequences.

## Features

- Accepts coding DNA sequences through manual input or FASTA files
- Detects single-nucleotide polymorphisms (SNPs)
- Converts coding DNA sequences into RNA
- Translates RNA sequences into amino acid sequences
- Supports reading frames 1, 2, and 3
- Classifies mutations as synonymous, missense, nonsense, or stop-loss
- Compares reference and experimental sequences
- Includes automated tests using pytest

## Requirements

- Python 3
- tabulate
- pytest

## Installation

Install the required packages:

pip install -r requirements.txt

## Usage

Run the program with:

python genetic_variant_interpreter.py

The program will ask whether you want to enter sequences manually or provide FASTA files.

## Testing

Run the automated tests with:

pytest
