## Blast Alignment Analysis Using Biopython

"""The program for Blast analysis using Biopython is based on Biopython Tutorial and Cookbook (http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)"""

## The IL-6 cDNA sequence (IL-6_cds_seq.fasta) is used for the demonstration of this program.


class Blast_alignment_analysis_program:
    
    def __init__(self):
        
        import os
        
        # Use a while loop to validate the entered file name
        filename = input("Enter the fasta file name for Blast analysis: ").strip()
        while filename not in os.listdir():
            print("The file name does not exist. Please try again.")
            filename = input("Re-enter the fasta file name for Blast analysis: ").strip()
        
        # Use a while loop to validate the entered blast program
        blast_program = input("Enter the Blast program (blastn, blastp, blastx, tblastn, tblastx): ").strip()
        while blast_program not in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']:
            print("The provided Blast program is invalid. Please try again.")
            blast_program = input("Re-enter the Blast program (blastn, blastp, blastx, tblastn, tblastx): ").strip()
        
        # Use a while loop to validate the entered database
        link = "https://ncbi.github.io/blast-cloud/blastdb/available-blastdbs.html"
        print(f"The information for choosing a database can be found in the following link:\n{link}")
        db_type = input("Enter a database for Blast analysis: ").strip()
        db_list = ['nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', 'pdbaa', 'pdbnt']
        
        while db_type not in db_list:
            print("The provided database is invalid. Please try again.")
            db_type = input("Re-enter a database for Blast analysis: ").strip()
        
        output_filename = input("Enter the output file name (filename.xml) for saving Blast analysis result: ").strip()
        
        self.filename = filename
        self.blast_program = blast_program
        self.db_type = db_type
        self.output_filename = output_filename
    
    # Define a class method for Blast alignment analysis
    def Blast_alignment_analysis(self):
        
        # This function performs Blast alignment analysis for a query sequence
        # The Blast alignment result is exported as a XML file as the output
        from Bio.Blast import NCBIWWW
        from Bio import SeqIO
        
        record = SeqIO.read(self.filename, format="fasta")
        
        result_handle = NCBIWWW.qblast(self.blast_program, self.db_type, record.format("fasta"))
        print("Blast alignment analysis is complete.")
        
        with open(self.output_filename, "w") as out_handle:
            out_handle.write(result_handle.read())
        
        result_handle.close()
    
    
    # Define a class method to display the Blast alignment analysis results
    def Blast_result_displayer(self):
        
        import os
        from Bio.Blast import NCBIXML
        
        option = input(f"Is {self.output_filename} opened for displaying Blast results (yes/no)? ").strip().upper()
        
        if option == "YES":
            xml_file = self.output_filename
        
        else:
            xml_file = input("Enter the Blast result file (XML) for display: ").strip()
            while xml_file not in listdir():
                print("The entered file name does not exist. Please try again.")
                xml_file = input("Re-enter the Blast result file (XML) for display: ").strip()
        
        display_no = int(input("Enter the number of top alignments from Blast analysis: "))
        E_Value_Threshold = float(input("Enter the E-value threshold: "))
        
        with open(xml_file) as result_handle:
            blast_record = NCBIXML.read(result_handle)
        
        count = 0
        for alignment in list(blast_record.alignments)[:display_no]:
            count += 1
            for hsp in alignment.hsps:
                if hsp.expect < E_Value_Threshold:
                    print(f"\n****Alignment #{count}****")
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e-value:", hsp.expect)
                    for i in range(0, len(hsp.query), 70):
                        print(hsp.query[i:i+70])
                        print(hsp.match[i:i+70])
                        print(hsp.sbjct[i:i+70])


# Assign the called program class to an object and check created object attributes
blast_analysis = Blast_alignment_analysis_program()
print(blast_analysis.filename)
print(blast_analysis.blast_program)
print(blast_analysis.db_type)
print(blast_analysis.output_filename)


# Call the class method Blast_alignment_analysis() to perform Blast alignment analysis
blast_analysis.Blast_alignment_analysis()


# Call the class method Blast_result_displayer() to show the Blast analysis results
blast_analysis.Blast_result_displayer()

