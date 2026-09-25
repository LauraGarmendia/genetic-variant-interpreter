import sys
from tabulate import tabulate



sequences = {}

def main():
    print("")
    print("---- GENETIC VARIANT INTERPRETER ----")
    print("")
    print("")
    print("How would you like to provide the sequences?")
    print("")
    print("1. Enter sequences manually")
    print("2. Upload two FASTA files")
    print("3. Upload one FASTA file containing two sequences")
    print("")

    try:
        choice = int(input("Please enter the number corresponding to your choice: ").rstrip())
        print("")

        if choice == 1:
            sequences, reading_frame = option_1()
        elif choice == 2:
            sequences, reading_frame = option_2()
        elif choice == 3:
            sequences, reading_frame = option_3()
        else:
            print("Invalid input. Please enter a valid number.")
            print("")
            return

    except ValueError:
        print("Invalid input. Please enter a valid number.")
        print("")
        return

    SNPs = analyze_SNP(sequences)
    print("")
    analyze_RNA(sequences)
    print("")
    analyze_protein(sequences, reading_frame)
    print("")
    classify_mutations(sequences, reading_frame, SNPs)





def option_1():

    experimental = input("Please enter the experimental coding DNA sequence: ").rstrip().upper()
    print("")

    if experimental == "" or not set(experimental.upper()).issubset({"A", "C", "T", "G"}):
        print("Experimental sequence is not valid.")
        print("")
        sys.exit(1)

    reference = input("Please enter the reference coding DNA sequence: ").rstrip().upper()
    print("")

    if reference == "" or not set(reference.upper()).issubset({"A", "C", "T", "G"}):
        print("Reference sequence is not valid.")
        print("")
        sys.exit(1)

    sequences["experimental"] = experimental
    sequences["reference"] = reference

    reading_frame = input("Enter the desired reading frame (1,2,3): ").rstrip()

    if reading_frame not in ["1", "2", "3"]:
        print("Invalid reading frame. Reading frame must be 1, 2, or 3.")
        print("")
        sys.exit(1)

    return sequences, reading_frame






def option_2():

    experimental_file = input("Enter the experimental FASTA file name: ").rstrip()

    try:
        with open(experimental_file, "r") as f:
            experimental = f.read()
    except:
        print("")
        print("Invalid experimental file name.")
        print("")
        sys.exit(1)

    experimental_lines = experimental.splitlines()
    sequences["experimental"] = "".join(experimental_lines[1:]).upper()

    reference_file = input("Enter the reference FASTA file name: ").rstrip()

    try:
        with open(reference_file, "r") as f:
            reference = f.read()
    except:
        print("")
        print("Invalid reference file name.")
        print("")
        sys.exit(1)

    reference_lines = reference.splitlines()
    sequences["reference"] = "".join(reference_lines[1:]).upper()

    if sequences["experimental"] == "" or not set(sequences["experimental"]).issubset({"A", "C", "T", "G"}):
        print("Experimental sequence is not valid.")
        print("")
        sys.exit(1)

    if sequences["reference"] == "" or not set(sequences["reference"]).issubset({"A", "C", "T", "G"}):
        print("Reference sequence is not valid.")
        print("")
        sys.exit(1)

    reading_frame = input("Enter the desired reading frame (1,2,3): ")
    print("")

    if reading_frame not in ["1","2","3"]:
        print("Invalid reading frame. Reading frame must be 1, 2, or 3.")
        print("")
        sys.exit(1)

    return sequences, reading_frame






def option_3():
    fasta_file = input("Enter the FASTA file name: ").rstrip()

    try:
        with open(fasta_file, "r") as f:
            file = f.read()
    except:
        print("")
        print("Invalid file name.")
        print("")
        sys.exit(1)

    reading_frame = input("Enter the desired reading frame (1,2,3): ")
    print("")

    if reading_frame not in ["1","2","3"]:
        print("Invalid reading frame. Reading frame must be 1, 2, or 3.")
        print("")
        sys.exit(1)

    seq_list = file.split(">")
    seq_list.remove('')

    if len(seq_list) != 2:
        print("The FASTA file must contain exactly two sequences.")
        print("")
        sys.exit(1)

    experimental = seq_list[0].splitlines()
    reference = seq_list[1].splitlines()

    sequences["experimental"] = "".join(experimental[1:]).upper()
    sequences["reference"] = "".join(reference[1:]).upper()

    if sequences["experimental"] == "" or not set(sequences["experimental"]).issubset({"A", "C", "T", "G"}):
        print("Experimental sequence is not valid.")
        print("")
        sys.exit(1)

    if sequences["reference"] == "" or not set(sequences["reference"]).issubset({"A", "C", "T", "G"}):
        print("Reference sequence is not valid.")
        print("")
        sys.exit(1)

    return sequences, reading_frame






def analyze_SNP(sequences):

    experimental = sequences["experimental"]
    reference = sequences["reference"]

    if len(experimental) != len(reference):
        print("Insertions and deletions are not currently supported.")
        print("")
        sys.exit(1)

    SNPs = []

    for position, (ref, exp) in enumerate(zip(reference, experimental), start = 1):
        if exp != ref:
            SNPs.append([position, ref, exp])
    print("")
    print("--- Summary of found SNPs ---")
    print("")
    print(tabulate(SNPs, headers=["Position", "Reference", "Experimental"], tablefmt = "github"))
    print("")

    return SNPs




def analyze_RNA(sequences):

    experimental = sequences["experimental"]
    reference = sequences["reference"]
    ref_rna = reference.replace("T", "U")
    exp_rna = experimental.replace("T", "U")

    print("--- RNA Sequence ---")
    print("")
    print(f"Reference RNA sequence: {ref_rna}")
    print(f"Experimental RNA sequence: {exp_rna}")
    print("")

    return ref_rna, exp_rna


def analyze_protein(sequences, reading_frame):

    experimental = sequences["experimental"]
    reference = sequences["reference"]
    exp_rna = experimental.replace("T", "U")
    ref_rna = reference.replace("T", "U")

    stop_codon = ["UAA", "UAG", "UGA"]

    codon_table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UAU": "Y", "UAC": "Y",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "UGU": "C", "UGC": "C", "UGG": "W",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }


    ref_protein = []
    exp_protein = []


    for i in range(int(reading_frame)-1, len(ref_rna), 3):

        codon = ref_rna[i:i+3]

        if len(codon) < 3:
            break
        elif codon in stop_codon:
            break
        elif codon in codon_table:
            ref_protein.append(codon_table[codon])

    ref_protein_sequence = "".join(ref_protein)


    for i in range(int(reading_frame)-1, len(exp_rna), 3):

        codon = exp_rna[i:i+3]
        if len(codon) < 3:
            break
        elif codon in stop_codon:
            break
        elif codon in codon_table:
            exp_protein.append(codon_table[codon])

    exp_protein_sequence = "".join(exp_protein)

    print("--- Protein Sequence ---")
    print("")
    print(f"Reference protein sequence: {ref_protein_sequence}")
    print(f"Experimental protein sequence: {exp_protein_sequence}")
    print("")

    return ref_protein_sequence, exp_protein_sequence



def classify_mutations(sequences, reading_frame, SNPs):


    stop_codon = ["UAA", "UAG", "UGA"]

    codon_table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UAU": "Y", "UAC": "Y",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "UGU": "C", "UGC": "C", "UGG": "W",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", "UAA": "STOP",
"UAG": "STOP",
"UGA": "STOP"
    }


    experimental = sequences["experimental"]
    reference = sequences["reference"]
    exp_rna = experimental.replace("T", "U")
    ref_rna = reference.replace("T", "U")
    reading_frame = int(reading_frame)

    mutations = []

    for SNP in SNPs:
        position = SNP[0]
        ref_nucleotide = SNP[1]
        exp_nucleotide = SNP[2]

        index_position = position - 1

        codon_start = index_position - ((index_position - (reading_frame - 1)) % 3)

        ref_codon = reference[codon_start : codon_start + 3]
        exp_codon = experimental[codon_start : codon_start + 3]

        ref_rna_codon = ref_codon.replace("T", "U")
        exp_rna_codon = exp_codon.replace("T", "U")

        if len(ref_rna_codon) < 3 or len(exp_rna_codon) < 3:
            continue

        ref_aminoacid = codon_table[ref_rna_codon]
        exp_aminoacid = codon_table[exp_rna_codon]

        if ref_aminoacid == exp_aminoacid:
            mutation_type = "Synonymous"
        elif ref_aminoacid != exp_aminoacid:
            if exp_aminoacid == "STOP":
                mutation_type = "Nonsense"
            elif ref_aminoacid == "STOP":
                mutation_type = "Stop-loss"
            else:
                mutation_type = "Missense"


        mutations.append((position, ref_nucleotide + " --> " + exp_nucleotide, ref_rna_codon + " --> " + exp_rna_codon, ref_aminoacid + " --> " + exp_aminoacid, mutation_type))

    print("")
    print("--- Summary of found mutations ---")
    print("")
    print(tabulate(mutations, headers=["Position", "DNA change", "RNA change", "Aminoacid change", "Mutation type"], tablefmt = "github"))
    print("")

    return mutations




if __name__ == "__main__":
    main()
