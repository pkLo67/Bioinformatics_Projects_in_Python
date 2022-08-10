# The OOP-based program for open reading frame analysis
# Use E_coli_genomic_DNA_seqs.fasta as an input for program testing
# ORF_analysis_result.fasta is the output from the program

class ORF_analysis_program:
    
    # Initiate the class and define the class attributes for fasta, orf_size_threshold,
    # and output filename
    def __init__(self):
        
        import os
        filename = input("Enter the fasta filename here:").strip()
        
        while filename not in os.listdir():
            print("The entered file name does not exist. Please try again.")
            filename = input("Re-enter the fasta filename here:").strip()
        
        with open(filename) as file:
            fasta = file.readlines()
        
        error = True
        
        while error:
            try:
                orf_size_threshold = int(input("Set the minimal size of ORF (default: 50 bp): "))
                if orf_size_threshold < 9:
                    print("Your entered ORF size is less than 9 bp. Please try again.")
                    error = True
                else:
                    error = False
            except:
                print("No entry detected. Please re-enter the minimal size of ORF.")
        
        output_filename = input("Enter the output filename for ORF analysis: ").strip()
        
        self.fasta = fasta
        self.orf_size_threshold = orf_size_threshold
        self.output_filename = output_filename

    # Define a class method to generate a dictionary to store ID and sequence 
    # extracted from the fasta
    def sequence_dictionary(self):
        seq_dict = {}
        for line in self.fasta:
            if ">" in line:
                seq = ""
                items = line.split()
                seq_ID = items[0] + "|" + " ".join([items[1], items[2], items[3]]) + "|"
            elif line == "\n":
                seq_dict[seq_ID] = seq
            else:
                seq += line.strip().upper().replace(" ", "")

        self.seq_dict = seq_dict
        return self.seq_dict
    
    # Define a class method to convert the DNA sequence into its reverse complement sequence
    def reverse_complement(self, sequence):
        nt_complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        sequence_c = ""

        for nt in sequence:
            sequence_c += nt_complement_dict.get(nt)

        sequence_rc = sequence_c[::-1]
        return sequence_rc
    
    # Define a class method to split the DNA sequence into 3-nt codons from three different frames
    def three_frames_codon_arrays(self, sequence):
        frame_1 = []
        frame_2 = []
        frame_3 = []
        three_frames = [frame_1, frame_2, frame_3]

        for i in range(3):
            for n in range(i, len(sequence)-2, 3):
                codon = sequence[n:n+3]
                three_frames[i].append(codon)
            
        return three_frames  
    
    
    # Define a class method to generate codon arrays from 6 different reading frames 
    # (3 frames from each forward strand and reverse-complement strand sequences)
    def six_frames_codon_arrays(self, sequence):
        sequence_rc = self.reverse_complement(sequence)
        frames_1 = self.three_frames_codon_arrays(sequence)
        frames_2 = self.three_frames_codon_arrays(sequence_rc)
        six_frames = frames_1 + frames_2
        six_frames_dict = {}
        for i in range(6):
            six_frames_dict[i+1] = six_frames[i]
        return six_frames_dict
    
    
    # Define a class method to identify the open reading frames from a DNA sequence
    def ORF_finder(self, frame_array, orf_size_threshold):
        import re

        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        orf = ""
        orf_ls = []

        for codon in frame_array:
            if codon == start_codon:
                if re.search(r"^ATG", orf):
                    orf += codon + " "
                else:
                    orf = ""
                    orf += codon + " "
            elif codon not in stop_codons:
                orf += codon + " "
            else:
                orf += codon
                orf_ls.append(orf)
                orf = ""

        orf_list = []
        for orf in orf_ls:
            if re.search(r"^ATG", orf):
                orf_list.append(orf)

        orf_list_select = []
        for orf in orf_list:
            if len(orf.replace(" ", "")) >= orf_size_threshold:
                orf_list_select.append(orf)

        return orf_list_select   
    
    
    # Define a class method to generate a dictionary to store all ORFs found from a sequence
    def ORF_dictionary(self, sequence, orf_size_threshold):
        orf_dict = {}
        six_frames_dict = self.six_frames_codon_arrays(sequence)
        for frame_num, frame_array in six_frames_dict.items():
            orf_list_select = self.ORF_finder(frame_array, orf_size_threshold)
            orf_dict[frame_num] = orf_list_select
        return orf_dict
    
    
    # Define a class method to generate a dictionary to store all ORFs found from multiple sequences
    def ORF_dictionary_multi_seqs(self, orf_size_threshold=50):
        orf_dict_multi_seq = {}
        for ID, sequence in self.seq_dict.items():
            orf_dict = self.ORF_dictionary(sequence, orf_size_threshold)
            orf_dict_multi_seq[ID] = [sequence, orf_dict]
        return orf_dict_multi_seq
    
    
    # Define a function to return gencode dictionary
    def gencode(self):
        codons_dict = {
                       'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                       'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                       'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                       'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                       'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                       'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                       'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                       'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                       'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                       'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                       'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                       'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                       'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                       'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                       'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                       'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        
        self.codons_dict = codons_dict        
        return self.codons_dict
    
    
    # Define a function to translate ORF into a protein sequence
    def protein_translation(self, orf):
        codons_dict = self.gencode()
        codons = orf.split()
        protein = ""
        for codon in codons:
            if codons_dict.get(codon) != "_":
                protein += codons_dict.get(codon)
            else:
                break
        return protein
    
    
    # Define a class method to write identified ORFs from all six reading frames into a file
    def ORF_analysis_output(self, orf_dictionary):
        import re
        output = open(self.output_filename, "w")
        for ID, seq_orf_dict in orf_dictionary.items():
            sequence = seq_orf_dict[0]
            orf_dict_ = seq_orf_dict[1]
            for num, orfs in orf_dict_.items():
                if num < 4:
                    for orf in orfs:
                        orf_ns = orf.replace(" ", "")
                        length = len(orf_ns)
                        start = re.search(orf_ns, sequence).start() + 1
                        end = re.search(orf_ns, sequence).end()
                        title = f"{ID}FRAME={num} POS={start} LEN={length}"
                        output.write(title + "\n")
                        for n in range(0, len(orf), 60):
                            output.write(orf[n:n+60] + "\n")
                        output.write("\n")
                        
                        protein = self.protein_translation(orf)
                        protein_id = f">The translated protein sequence, {len(protein)} amino acids."
                        output.write(f"{protein_id}\n")
                        for a in range(0, len(protein), 60):
                            output.write(protein[a:a+60] + "\n")
                        output.write("\n")
                        
                else:
                    for orf in orfs:
                        orf_ns = orf.replace(" ", "")
                        length = len(orf_ns)
                        sequence_rc = self.reverse_complement(sequence)
                        start = re.search(orf_ns, sequence_rc).start() + 1
                        end = re.search(orf_ns, sequence_rc).end()
                        title = f"{ID}FRAME={num} POS={-(start)} LEN={length}"
                        output.write(title + "\n")
                        for n in range(0, len(orf), 60):
                            output.write(orf[n:n+60] + "\n")
                        output.write("\n")
                        
                        protein = self.protein_translation(orf)
                        protein_id = f">The translated protein sequence, {len(protein)} amino acids."
                        output.write(f"{protein_id}\n")
                        for a in range(0, len(protein), 60):
                            output.write(protein[a:a+60] + "\n")
                        output.write("\n")
                        
        output.close()
        
    
    # Define a main() method to execute the entire ORF analysis program
    def ORF_analysis_main(self):
        
        print("The program can identify open reading frames from multiple DNA sequences")
        print("-" * 80)
        
        self.sequence_dictionary()
        orf_dict_multi_seq = self.ORF_dictionary_multi_seqs(self.orf_size_threshold) 
        self.ORF_analysis_output(orf_dict_multi_seq)
        
        print(f"The ORF analysis is done and the output result has been written into {self.output_filename}.")
        

# Assign the program class and call the main method for ORF 
ORF_analysis = ORF_analysis_program()
ORF_analysis.ORF_analysis_main()

