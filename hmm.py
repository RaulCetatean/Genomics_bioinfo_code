#hmm

fasta = open("chr22.fa", "r")


def get_genome(fasta):
        genome = ""
        for line in fasta:
                if line[0] != ">" or "N":
                        genome += line[:].upper()

        return genome 
print(get_genome)
