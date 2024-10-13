# NucleoTools

NucleoTools is a toolkit for filtering sequenced DNA reads. It also provides functionality for converting DNA and RNA sequences. Additionally, this script contains two functions for processing biological data: converting multiline sequences to single-line format in FASTA and extracting protein names that best match the database from BLAST results.

## Contents
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)

## Description

### 1. Filtering
The `filter_fastq` function allows you to filter DNA sequences based on several parameters:
- Average read quality (`quality_threshold`)
- GC content of the read (`gc_bounds`)
- Length of the read (`length_bounds`)

### 2. DNA and RNA Conversion
The `run_dna_rna_tools` function allows you to perform the following actions on DNA and RNA sequences:
- Transcribe a coding DNA strand to mRNA (`transcribe`)
- Write the sequence in reverse order (`reverse`)
- Create a complementary sequence (`complement`)
- Write the complementary sequence in reverse order (`reverse_complement`)

### 3. Converting Multiline Sequences to Single Line
The `convert_multiline_fasta_to_oneline` function reads the specified FASTA file, where sequences (DNA/RNA/proteins, etc.) may be split across multiple lines, and saves the data to a new FASTA file where each sequence is contained in a single line.

### 4. Extracting Protein Names for Best Matches
The `parse_blast_output` function extracts protein names from a text file with BLAST results that correspond to the best match and saves the names in alphabetical order to a new file.

## Installation

1. Clone the repository:
    ```bash
    git clone git@github.com:Polina1010123/dna-rna-tools-and-fastq-filtrator.git
    ```

2. Navigate to the project directory:
    ```bash
    cd path-to-project-directory
    ```

3. Set the correct path to the `modules` directory in `fastq_filtrator.py` for successful imports.

4. Set the correct paths for reading and writing your data files in `fastq_filtrator.py`.

## Usage

1. **Filtering**
   
   In `fastq_filtrator.py`, specify the path to the FASTQ file you want to filter in the `input_fastq` argument of the `filter_fastq` function, and the path to the file where you want to save the filtered sequences in the `output_fastq` argument.
    
    ```python
    def filter_fastq(
        input_fastq: str = "c:/Users/kozlo/Downloads/example_fastq.fastq",
        quality_threshold: float = 20,
        gc_bounds: float | tuple[float, float] = (20, 80),
        length_bounds: int | tuple[float, float] = (80),
        output_fastq: str = "c:/Users/kozlo/Downloads/filtered/filtered_sequences.fastq",
        ) -> dict[str, tuple[str, str]]:
    ```
   
   In the same file, set the parameters for filtering: `quality_threshold`, `gc_bounds`, and `length_bounds`. 
   Note that in the `quality_threshold` parameter you can specify only a single number, while in `gc_bounds` and `length_bounds` you can provide either a range of values or a single number. In the latter case, reads with a GC content or length greater than the specified number will be removed. 
   In the example above, reads with an average quality score above 20, GC content from 20 to 80%, and length less than 80 bp will be saved.

2. **DNA and RNA Conversion**

   In `fastq_filtrator.py`, specify the sequence for conversion and the type of conversion in the `run_dna_rna_tools` function.
   
   Examples:
    ```python
    # Transcribe a coding DNA strand to mRNA
    result = run_dna_rna_tools("ATG", "transcribe")
    print(result)
    ```
    ```python
    # Write the sequence in reverse order
    result = run_dna_rna_tools("ATG", "reverse")
    print(result)
    ```

3. **Converting Multiline Sequences to Single Line**

    In the `convert_multiline_fasta_to_oneline` function, specify the path to the input FASTA file (argument `input_fasta`) containing multiline sequences, as well as the path to the output FASTA file (argument `output_fasta`) where single-line sequences will be written.

4. **Extracting Protein Names for Best Matches**

    In the `parse_blast_output` function, specify the path to the input text file with BLAST results (argument `input_file`), and the path to the output text file where the best matches will be written (argument `output_file`).
