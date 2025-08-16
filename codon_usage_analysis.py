from Bio import SeqIO

# opening files here
uniprot_blast = open("SK926\SK926.blast_swissprot", "r")
fasta = open("SK926\SK926.fna", "r")
dat_file = open("SK926\swissprot_db.dat", "r")

# saving the phyla of each domain here
bacteria_phyla = []
archaea_phyla = []
eukaryota_phyla = []
viruses_phyla = []

# list of dictionaries that connect phylum, organism ID, and contig ID together
list_of_bacteria_phyla_dict = []
list_of_eukaryota_phyla_dict = []
list_of_archaea_phyla_dict = []
list_of_viruses_phyla_dict = []

# the dictionaries that contain phylum, organism ID, and contig ID
bacteria_phyla_dict = {}
eukaryota_phyla_dict = {}
archaea_phyla_dict = {}
viruses_phyla_dict = {}

# dictionary where key = organism_id and value = contig_id
dict_contig = {}
for line in uniprot_blast.readlines()[:len(uniprot_blast.readlines()) - 2]: # the -2 because there are 2 empty lines at the end
    line = line.strip().split("\t")
    organism_id = line[1]
    contig_id = line[0]
    dict_contig[organism_id] = contig_id

# parsing the phyla of each domain
for line in dat_file.readlines()[:len(dat_file.readlines()) - 3]: # the -3 because there are 3 empty lines at the end
    line = line.strip()
    info = line.split("\t") # splitting over the tab to create 3 elements, the third of which contains taxonomy information and first contains organism ID
    Organism_ID = info[0]
    taxa = info[2] # retrieving taxonomy information
    taxa_detailed = taxa.split(";") # this is to separate the domain, class, phylum, etc.
    if (len(taxa_detailed) >= 3): # to make sure the taxonomy information contains all information
        phylum = taxa_detailed[2].strip()

        if taxa_detailed[0].strip() == "Bacteria": # if domain is bacteria, save the phylum (if it was not saved before) to its respective list (bacteria_phyla_

            if phylum not in bacteria_phyla:
                bacteria_phyla.append(phylum)

            # Connecting phylum, Organism ID, and Contig ID
            bacteria_phyla_dict["Phylum"] = phylum
            bacteria_phyla_dict["Organism_ID"] = Organism_ID
            bacteria_phyla_dict["Contig_ID"] = dict_contig[Organism_ID]
            list_of_bacteria_phyla_dict.append(bacteria_phyla_dict.copy())

        elif taxa_detailed[0].strip() == "Eukaryota": # so on for the rest of domains

            if phylum not in eukaryota_phyla:
                eukaryota_phyla.append(phylum)

            eukaryota_phyla_dict["Phylum"] = phylum
            eukaryota_phyla_dict["Organism_ID"] = Organism_ID
            eukaryota_phyla_dict["Contig_ID"] = dict_contig[Organism_ID]
            list_of_eukaryota_phyla_dict.append(eukaryota_phyla_dict.copy())

        elif taxa_detailed[0].strip() == "Archaea":

            if phylum not in archaea_phyla:
                archaea_phyla.append(phylum)

            archaea_phyla_dict["Phylum"] = phylum
            archaea_phyla_dict["Organism_ID"] = Organism_ID
            archaea_phyla_dict["Contig_ID"] = dict_contig[Organism_ID]
            list_of_archaea_phyla_dict.append(archaea_phyla_dict.copy())

        elif taxa_detailed[0].strip() == "Viruses":

            if phylum not in viruses_phyla:
                viruses_phyla.append(phylum)

            viruses_phyla_dict["Phylum"] = phylum
            viruses_phyla_dict["Organism_ID"] = Organism_ID
            viruses_phyla_dict["Contig_ID"] = dict_contig[Organism_ID]
            list_of_viruses_phyla_dict.append(viruses_phyla_dict.copy())

    else: # if phylum name is not found, ignore this organism and move on
        continue

def connect_contig_with_seq(): # function that makes a dictionary where the key is the contig id and the value is the sequence
    dict = {}
    # I want to parse only the first N lines, so I create a new file that has N lines from the original fasta file. N in this example is 2000
    new_file = open("new.fasta", "a")
    original_file = open("SK926\SK926.fna", "r").readlines()
    for line in original_file[:2000]:
        new_file.write(line)
    for record in SeqIO.parse("new.fasta", "fasta"): # I parse with Biopython package
        dict[record.id] = str(record.seq)
    return dict

def get_contig_id_from_phylum(phylum): # function that takes a phylum and gives me all the contig IDs of that phylum
    phylum_contig_IDs = []
    if phylum in bacteria_phyla:
        for dict in list_of_bacteria_phyla_dict:
            if dict["Phylum"] == phylum:
                phylum_contig_IDs.append(dict["Contig_ID"])
    elif phylum in eukaryota_phyla:
        for dict in list_of_eukaryota_phyla_dict:
            if dict["Phylum"] == phylum:
                phylum_contig_IDs.append(dict["Contig_ID"])
    elif phylum in archaea_phyla:
        for dict in list_of_archaea_phyla_dict:
            if dict["Phylum"] == phylum:
                phylum_contig_IDs.append(dict["Contig_ID"])
    elif phylum in viruses_phyla:
        for dict in list_of_viruses_phyla_dict:
            if dict["Phylum"] == phylum:
                phylum_contig_IDs.append(dict["Contig_ID"])
    else:
        pass

    return phylum_contig_IDs

def get_sequence_from_phylum(phylum): # function that takes a phylum and gives me all the sequences in that phylum --> by using the above two functions
    ids = get_contig_id_from_phylum(phylum)
    dict_seqs = connect_contig_with_seq()
    list_of_seq = []
    for id in ids:
        list_of_seq.append(dict_seqs.get(id, "N/A"))
    list_of_seq = [i for i in list_of_seq if i != "N/A"]

    return list_of_seq

# creating files to save the sequences of each phylum
for phylum in bacteria_phyla:
    f_name = "Sequences\Bacteria/" + phylum + ".txt"
    f = open(f_name, "a")
    f.writelines(seq + '\n' for seq in get_sequence_from_phylum(phylum)) # writing the sequences of bacterial phyla, each in its own line
for phylum in eukaryota_phyla:
    f_name = "Sequences\Eukaryota/" + phylum + ".txt"
    f = open(f_name, "a")
    f.writelines(seq + '\n' for seq in get_sequence_from_phylum(phylum)) # writing the sequences of eykaryotic phyla, each in its own line
for phylum in archaea_phyla:
    f_name = "Sequences\Archaea/" + phylum + ".txt"
    f = open(f_name, "a")
    f.writelines(seq + '\n' for seq in get_sequence_from_phylum(phylum)) # writing the sequences of archaeal, each in its own line
for phylum in viruses_phyla:
    f_name = "Sequences\Viruses/" + phylum + ".txt"
    f = open(f_name, "a")
    f.writelines(seq + '\n' for seq in get_sequence_from_phylum(phylum)) # writing the sequences of viral phyla, each in its own line

# dictionary where the key is the amino acid and the value is a list of all the corresponding codons
aa_dict = {
"A": ["GCG", "GCA", "GCT", "GCC"],
"C": ["TGT", "TGC"], "D": ["GAT", "GAC"], "E": ["GAG", "GAA"], "F":  ["TTT", "TTC"],
"G": ["GGG", "GGA", "GGT", "GGC"], "H": ["CAT", "CAC"], "I": ["ATA", "ATT", "ATC"],
"K": ["AAG", "AAA"], "L": ["TTG", "TTA", "CTG", "CTA", "CTT", "CTC"], "M": ["ATG"],
"N": ["AAT", "AAC"], "P": ["CCG", "CCA", "CCT", "CCC"], "Q": ["CAG", "CAA"],
"R": ["AGG", "AGA", "CGG", "CGA", "CGT", "CGC"], "S": ["AGT", "AGC", "TCG", "TCA", "TCT", "TCC"],
"T": ["ACG", "ACA", "ACT", "ACC"], "V": ["GTG", "GTA", "GTT", "GTC"], "W": ["TGG"],
"Y": ["TAT", "TAC"], "STOP": ["TGA", "TAG", "TAA"]}

# dictionary where the key is codon and the value is the corresponding amino acid
table = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
             "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
             "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP",
             "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
             "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
             "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
             "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
             "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
             "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
             "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
             "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
             "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
             "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
             "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
             "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
             "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}


# I will make only one comparison between bacilli and gammaproteobacteria as an example, but the logic shall be the same for all the other pairs
# I did that because there will be NC2 files for each domain, where N is the number of phyla. This is extremely a lot.
f1 = open("Sequences\Bacteria/""Bacilli.txt", "r")
f2 = open("Sequences\Bacteria/""Gammaproteobacteria.txt", "r")

dect_1 = {} # saving the codon count of phylum number one here
dect_2 = {} # saving the codon count of phylum number two here

def get_aa_from_codon(codon): # takes a codon and returns the corresponding amino acid
    return table[codon]
def get_codons_from_aa(aa): # takes an amino acid and returns the corresponding codon(s)
    return aa_dict[aa]

for line in f1.readlines(): # getting the codon counts of phylum number 1
    line = line.rstrip()
    for codon in table:
        if codon in dect_1:
            dect_1[codon] += line.count(codon)
        else:
            dect_1[codon] = line.count(codon)

for line in f2.readlines(): # getting the codon counts of phylum number 2
    line = line.rstrip()
    for codon in table:
        if codon in dect_2:
            dect_2[codon] += line.count(codon)
        else:
            dect_2[codon] = line.count(codon)

output_name = "Bacilli-Gammaproteobacteria.txt" # name of output file
output_f = open(output_name, "w")
output_f.write("Bacilli\tGammaproteobacteria\n") # headers to indicate the position of each phylum in the file, so that we know which codon usage is for which

for codon in table:
    # the next few lines are to know the division factor for each codon.
    # The general formula for codon usage is: codon_count / sum_of_all_codon_counts_that_correspond_to_the_same_amino_acid
    # division factor is the denominator in the formula above --> sum_of_all... etc.
    aa = get_aa_from_codon(codon)
    codons = get_codons_from_aa(aa)
    division_factor_1 = 0
    division_factor_2 = 0
    for c in codons:
        division_factor_1 += dect_1[c]
        division_factor_2 += dect_2[c]
        # writing the codon usage in the file, where the codon usage is codon count / division factor
        # if we want percentage we can simply multiply by 100
    output_f.write(codon+"\t"+str(dect_1[codon] / division_factor_1)+"\t"+str(dect_2[codon] / division_factor_2)+"\n")

# AND WE'RE DONE!!! TOOK ME THE WHOLE DAY BUT IT'S WORTH IT.