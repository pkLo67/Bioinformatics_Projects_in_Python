## The Program for Analysis of Restriction Enzyme Cut Sites in a DNA Sequence (OOP Code)

class RE_cut_analysis_program:
    
    def __init__(self):
        
        import os
        # Use the while loop to validate whether the entered file name exists in the working directory
        seq_file_name = input("Enter the file name of the DNA sequence file: ").strip()
        while seq_file_name not in os.listdir():
            print("The entered file name is not found. Please try again.")
            seq_file_name = input("Re-enter the sequence file name: ").strip()
        
        RE_file_name = input("Enter the file name with the restriction enzyme list: ")
        while RE_file_name not in os.listdir():
            print("The entered file name is not found. Please try again.")
            RE_file_name = input("Re-enter the RE list file name: ").strip()

        RE_name = input("Enter a restriction enzyme name for analysis: ").strip()
                
        self.seq_file_name = seq_file_name
        self.RE_file_name = RE_file_name
        self.RE_name = RE_name
    
    # Define a function to extract the DNA sequence from the fasta file for the downstream analysis
    def sequence_assembler(self):
        file = open(self.seq_file_name)

        for line in file:
            if line.startswith(">"):
                seq = ""
            else:
                seq += line.strip().upper()

        file.close()
        
        self.sequence = seq

        return self.sequence
    
    # Define a function to create a dictionary to store RE names and their cut sites
    def RE_dictionary(self):

        import re
        RE_dict = {}

        # Create a dictionary for nucleotide symbols that are used in RE recognition sites
        nt_codes = {"A":"A", "C":"C", "G":"G", "T":"T",
                    "R":"[AG]", "Y":"[CT]", "M":"[AC]",
                    "K":"[GT]", "S":"[CG]", "W":"[AT]",
                    "B":"[CGT]", "D":"[AGT]", "H":"[ACT]",
                    "V":"[ACG]", "N":"[ACGT]"}

        # Use a for loop to store RE names as keys and their recognition sites as values
        # in a dictionary. The cut positions and lengths of RE recognition sites are also
        # stored as parts of values.

        f = open(self.RE_file_name)
        for line in f:
            items = line.strip().split()
            RE_name = items[1].upper()
            RE_site = items[0]
            RE_dict[RE_name] = [RE_site.replace("/", ""), RE_site.index("/"), len(RE_site)-1]

        f.close()    

        # Use a for loop to convert RE recognition sites into regular expression patterns
        for RE_key, RE_site in RE_dict.items():
            for nt_key, nt_value in nt_codes.items():
                RE_site[0] = RE_site[0].replace(nt_key, nt_codes[nt_key])
                RE_dict[RE_key] = RE_site
        
        self.RE_dict = RE_dict
        return self.RE_dict
    
    
    # Define a function to identify RE cut sites in a sequence
    def RE_cut_analysis(self):
        import re

        RE_dict = self.RE_dictionary()
        dna_seq = self.sequence_assembler()

        while self.RE_name.upper() not in RE_dict.keys():
            print("Your entered RE name is invalid and not found.")
            self.RE_name = input("Please re-enter the RE name: ").strip()
        RE_name = self.RE_name
        # Use the re.finditer function to identify RE sites in the DNA sequence
        RE_sites = []
        all_starts = []
        all_ends = []
        all_cuts = [0]
        for match in re.finditer(RE_dict[RE_name.upper()][0], dna_seq):
            RE_sites.append(match.group(0))
            all_starts.append(match.start()+1)
            all_ends.append(match.end())
            all_cuts.append(match.start() + RE_dict[RE_name.upper()][1])

        all_cuts.append(len(dna_seq))

        print(f"      The Identified {RE_name} Cut Sites in the Analyzed Sequence:")
        print("-" * 70)
        print(f"       {RE_name} cut site              Start position       End position")
        print("-----------------------------     --------------       -------------")
        for RE_site, start, end in zip(RE_sites, all_starts, all_ends):
            len_1 = round((30-len(RE_site))/2)
            len_2 = 37-len(RE_site)-len_1
            print(" " * len_1 +  RE_site + " " * len_2 + str(start) +\
                  " " * (60-len(RE_site)-len(str(start))-len_1-len_2) + str(end))
        print("\nNote: The program displays the start and end nucleotide positions of RE cut sites.")

        # Use a for loop to calculate DNA fragment sizes after RE digestion
        print("\nThe predicted DNA fragment sizes after RE digestion:")
        for i in range(1,len(all_cuts)):
            current_cut_position = all_cuts[i]
            prior_cut_position = all_cuts[i-1]
            fragment_size = current_cut_position - prior_cut_position
            print(f"The size of fragment #{i} is {fragment_size}")


# Assign the program class to a class object and initiate the program
RE_program = RE_cut_analysis_program()
print(RE_program.seq_file_name)
print(RE_program.RE_file_name)
print(RE_program.RE_name)


# Call the sequence_assembler() method and print the generated sequence
RE_program.sequence_assembler()
print(RE_program.sequence)


# Call the RE_dictionary() method and print items in the RE dictiionary
RE_program.RE_dictionary()
for item in RE_program.RE_dict.items():
    print(item)


# Call the RE_cut_analysis() method to execute the RE cut analysis program
RE_program.RE_cut_analysis()


# Check whether the RE_cut_analysis() method can validate the entered RE name
RE_program.RE_name = "ScaXXX"
RE_program.RE_cut_analysis()

