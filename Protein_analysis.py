#!/usr/bin/python3
import os,shutil,subprocess,sys,re,glob,fileinput
#Importing all the libraries I use
try:
	import ete3
except ModuleNotFoundError:
	print('ete3 not installed, installing it now...')
	os.system('pip install ete3')
	import ete3
try:
	import pandas as pd
except ModuleNotFoundError:
	print('pandas not installed, installing it now...')
	os.system('pip install pandas')
	import pandas as pd
try:
	import numpy as np
except ModuleNotFoundError:
	print('numpy not installed, installing it now...')
	os.system('pip install numpy')
	import numpy as np
try:
	import argparse
except ModuleNotFoundError:
	print('argparse not installed, installing it now...')
	os.system('pip install argparse')
	import argparse
from ete3 import Tree
from collections import Counter

def get_sequences():
	#use clustalo to align, and downloading clustalo's infered tree, in this case it is a neighbor-joining tree based on the alignmentc
	print('Aligning sequences + building phylogenetic tree...')
	os.system('clustalo --threads=10 -i protein_sequences/input.fa -o protein_sequences/aligned_sequences.fa --guidetree-out=clustalo_tree') # I don't think the --threads is helping much

#Adding arguments and a small description that is displayed if the user types -h
def new_arg(name,abreviation,description):
	parser.add_argument('--'+name,'-'+abreviation, type=str, help=description)

def make_acc_dict():
	acc_list=[]
	species=[]
	acc_spp = {}
	for line in open("acc_spp.txt", "r"):#Tab sepparated between acc and spp
		acc_list.append(line.split('\t')[0])
		species.append(line.split('\t')[1].strip('\n'))#Getting rid of the newline character
		acc_spp[line.split('\t')[0]] = line.split('\t')[1].strip('\n')
	non_redundant_species = set(species)
	return non_redundant_species,acc_spp,species,acc_list

parser=argparse.ArgumentParser(
    description='''proteins.py performs an analysis based on user input.\n\
The user can provide a Taxonomic Group and a Protein Family of their choice.\n''',
    epilog="""Please message <my.email> to report any bugs.""")
new_arg('taxonomic_group','t','Taxonomic group of your choice')
new_arg('protein_family','p','Protein family of your choice.')
new_arg('dir','d','Directory where program will be executed')
new_arg('winsize','w','Window size of conservation plot and hydropathy plot')
new_arg('partial','partial','y if you want to include partial sequences, anything else will exclude them')


args=parser.parse_args()#Parsing arguments and assigning them to variables
prot_fam = args.protein_family 
tax_grp = args.taxonomic_group
wd = args.dir
winsize = args.winsize
partial = args.partial


#If no arguments are provided, the user is asked to enter information that is then assigned to the variables
if wd is None:
	wd = input('Your current working directory is {}, please type the path where you want this analysis made:\n\
(the default is the current working directory, you can just press enter to stay in the current directory):\
	\n\t'.format(os.getcwd()))
if tax_grp is None:
	tax_grp = input("No taxonomic group was selected. Please enter the taxonomic group you want to investigate:\n\t")
if prot_fam is None:
	prot_fam = input("No protein family was selected. Please enter the desired protein family:\n\t")
if partial is None:
	partial = input('Do you want to include partial sequences? Press y if yes, any other key will discard partial sequences:\n\t')
if partial == 'y':
	partial = ''
else:
	partial = ' NOT PARTIAL'
#Changing directory to the one that was specified, with error traps that default to the current dir
try:
	os.chdir(wd)
except FileNotFoundError:
	print("Directory: {} does not exist, defaulting to current directory".format(wd))
except NotADirectoryError:
	print("{} is not a directory, defaulting to current directory".format(wd))
os.system("mkdir protein_sequences && mkdir protein_sequences/individual_sequences && mkdir prosite_motifs reports")

print('Fetching Accession values for {} proteins in {}...'.format(prot_fam,tax_grp))
print('Please quit now if this is not the desired query')

#Downloading species name and acc values, making a dictionary using the accession value as the key
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" \
| efetch -format gbc | xtract -pattern INSDSeq -element INSDSeq_accession-version -element INSDSeq_organism \
| cut -f1,2 >  acc_spp.txt')
non_redundant_species,acc_spp,species,acc_list = make_acc_dict()
print('There are {} protein sequences, and {} species that match this query.'.format(len(acc_list),len(non_redundant_species)))

#read the accession numbers to get the sequences, and print a warning if longer than 1000, or if there are 0 sequences, or 1 sequence

if len(acc_list) == 0:
	print('There are no protein sequences that match your query. Please restart with different input arguments.')
	print('Exiting analysis')
	sys.exit()
elif len(acc_list) == 1:
	print('There is only one protein that matches your query, this program can\'t use a single sequence as input. \
Please restart with different input arguments.')
	print('Exiting analysis')
	sys.exit()
if len(acc_list) > 1000:
	while True:
		limit = input('There are more than 1000 sequences in this analysis, do you want to set 1000 as the limit?\n \
Press y if you want the limit, n if you want to download all the sequences\n')
		if limit == 'y':
			max_seq = 1000
			break
		if limit == 'n':
			max_seq = None
			break
		else:
			choice = input('That was not an option...\nPlease select y or n:\t')
if len(acc_list) <= 1000:
	max_seq = None
try:
	if str(max_seq).isnumeric():
		print('Fetching 1000 Protein Sequences...')
		os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" \
		| efetch -format gbc -stop {}| xtract -pattern INSDSeq -element INSDSeq_accession-version -element INSDSeq_organism \
		| cut -f1,2 >  acc_spp.txt'.format(max_seq))
		non_redundant_species,acc_spp,species,acc_list = make_acc_dict()#Need to amke it again, or values will be wrong
		os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" \
		| efetch -format fasta -stop {} > protein_sequences/all_sequences.fa'.format(max_seq))
	else:
		print('Fetching all Protein Sequences...')
		os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" \
		| efetch -format fasta > protein_sequences/all_sequences.fa')

except AttributeError: #In casd there is an error, this will not happen
	print('Fetching all Protein Sequences...')
	os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" \
	| efetch -format fasta > protein_sequences/all_sequences.fa')


#Check the formats of the fasta file, to make sure that they will be analysed by pullseq. At this point, any sequence that doesn't start with the accession value will not be recognised by pullseq

for line in open("protein_sequences/all_sequences.fa", "r"):
	if line.startswith('>'):
		if any(acc in line for acc in acc_list) != True: #If the accession value is not in the header of the sequence
			print('WARNING! Unrecognised accession value\nThe following sequence:\n{}\n\
				is either not in fasta format, or is in a format that is not supported in this program\n'.format(line))

#Asking the user if they want the most/least common species, all of them, or a pick and choose approach
#I then select the accession values of those sequences and add them to a pullseq loop
#By doing the pullseq -g in a loop, I can also use the (annoying) sequences that start with the db name like sp|
#all sequences
while True:
	choice = input('\nDo you want to filter your sequences? Please select one of the following options, or press any other key to exit:\n\n\t\
-Use all sequences for the analysis (Will take a long time for large datasets) (Press 1)\n\n\t\
-Use sequences from the most common species (Press 2)\n\n\t\
-Use sequences from the least common species (Press 3)\n\n\t\
-Cycle one by one and chose which species you want (Not recommended for large datasets) (Press 4)\n\n\t')
	if choice == '1' :
		print('All sequences will be used for the analysis. \nThe alignment may take a while, depending on the amount of sequences')
		os.system('cp protein_sequences/all_sequences.fa protein_sequences/input.fa ')
		get_sequences()
		filtered_dict = acc_spp
		break
	#print the most common species by using the argument Counter
	elif choice == '2' :
		while True:
			spp_no = input('These are the 10 most common species, please select how many species you want to use:\n\t{}\n'\
					.format(Counter(species).most_common(10)))
			if str(spp_no).isnumeric():
				spp_no=int(spp_no)
				break
			else:
				print('please enter a number...')
		if int(spp_no) > 10:
			print('The number you selected is greater than the displayed number, Chosing 10...')
			spp_no = 10
		
		selected_species = []
		numbers = np.linspace(1,spp_no,spp_no)
		for x in numbers:
			selected_species.append(Counter(species).most_common(10)[int(x)-1][0]) #x-1 is the position of the [species, freq] in the list, and [0] is for the spp name
		
		filtered_acc = [] #Going back to the acc_spp dicctionary to extract the accession values of the chosen species
		for species in selected_species:
			keys = [key for key,value in acc_spp.items() if value == species]
			for key in keys:
				filtered_acc.append(key)
		filtered_dict = {}
		for acc in filtered_acc: #Make a new dict for the pd df later as the final report
			filtered_dict[acc] = acc_spp[acc]
		#pullseq loop to extract the selected species
		for acc in filtered_acc:
			os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -g {} \
				>> protein_sequences/input.fa'.format(acc))
		get_sequences()
		break
	#Same approach as above, only with the least common species
	elif choice == '3':
		while True:
			spp_no = input('These are the 10 least common species, please select how many species you want to use:\n\t{}\n'\
				.format(Counter(species).most_common()[-11:-1]))
			if str(spp_no).isnumeric():
				spp_no=int(spp_no)
				break
			else:
				print('please enter a number...')
		if int(spp_no) > 10:
			print('The number you selected is greater than the displayed number, Chosing 10...')
			spp_no = 10
		selected_species = []
		numbers = np.linspace(1,spp_no,spp_no)
		for x in numbers:
			selected_species.append(Counter(species).most_common()[-int(x)][0])
		filtered_acc = []
		for species in selected_species:
			keys = [key for key,value in acc_spp.items() if value == species]
			for key in keys:
				filtered_acc.append(key)
		filtered_dict = {}
		for acc in filtered_acc: #Make a new dict for the pd df later as the final report
			filtered_dict[acc] = acc_spp[acc]
		for acc in filtered_acc:
			os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -g {} \
			>> protein_sequences/input.fa'.format(acc))
		get_sequences()
		break
	#Pick and choose approach
	elif choice == '4':
		filtered_dict = {}
		for current_species in non_redundant_species:
			if input('There are {} occurrances of {} in this dataset. Do you want to include the sequences in the analysis?\n \
Press y if yes, any other key will discard the sequence(s)\n\t'\
			.format(species.count(current_species),current_species)) =='y':
				
				keys = [key for key,value in acc_spp.items() if value == current_species]
				for accession in keys:
					filtered_dict[accession]=acc_spp[accession]
					os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -g {} \
					>> protein_sequences/input.fa'.format(accession))
		get_sequences()
		break

#If user selects another option...
	else:
		print('That wasn\'t one of the options...\nPlease choose one of the options')
	


#Making tree file, everyone likes a tree
#Now I need to change all acc values in the tree file for the format acc_spp, as a tree that only has the accession value is not very informative
try:
	newtree = open("newick_tree.txt", "a")
	for line in open("clustalo_tree", "r"):
		if line.startswith('(') or line.startswith(')') or line.startswith(',') or line.startswith(';'):
			newtree.write(line)
		elif '|' in line[2:4]:#Those that have NCBI identifiers in the headers
			acc = line.split('|')[1] #many sequences are found in this format, if a sequence is downloaded with another format, such as tr, it will work, but the whole header will be there, instead of the acc and spp
			spp = acc_spp.get(acc)
			spp = spp.replace(' ','_')
			newtree.write(acc + '_' + spp +':' + line.split(':')[1])
		else:
			acc = line.split(':')[0]
			spp = acc_spp.get(acc)
			spp = spp.replace(' ','_')
			newtree.write(acc + '_' + spp +':' + line.split(':')[1]) 
	newtree.close()
	sys.stdout = open("tree.txt", "w") #tree is printed to stdout, so I send it to a file
	t = Tree("newick_tree.txt", format=1) #same version as the clustalo one, but I have formatted the seq labels to acc_spp
	print(t)
	sys.stdout.close()
	#reopening stdout, I need it to print messages
	sys.stdout = open("/dev/stdout", "w")
	#If winsize was not specified as argument, and adding warnings for large and small winsizes
	print('Done. File newick_tree.txt contains Guide tree written in Newick format.\n\
	This tree was generated from the alignment, using a neighbor-joining approach.\n \
	tree.txt contains a visual representation of the infered tree.') 
except:#If an error occurs when building the tree, the analysis can still proceed.
	sys.stdout = open("/dev/stdout", "w")
	print('Sequences aligned correctly, however an error was encountered when creating the newick tree.\n\
The tree can only be found in newick format.')
if winsize is None:
		winsize = input('Please specify a window size for the conservation plot (4 is recommended). \n\
The same window size will be used if you decide to perform more analysis:\n\t')
while True:
	if winsize.isnumeric() == True:
		break
	else:
		winsize = input('That was not a number...\nPlease choose a window size:\t')

if int(winsize) > 30:
	print('WARNING:\nThis window size is large, conservation plot might not be accurate if sequences are short')
if int(winsize) < 4:
	print('WARNING:\nThis window size is low, conservation plot might be noisy')

print('Producing conservation plot. Graph will appear on screen, but will also be saved as plotcon.svg')
#plotcon for svg and for screen output
os.system('plotcon -sequences protein_sequences/aligned_sequences.fa -winsize {} -graph x11 -graph svg \
	&& plotcon -sequences protein_sequences/aligned_sequences.fa -winsize {} -graph x11'.format(winsize,winsize))
print('Done')
print('Searching for known protein motifs in the query sequences...')

#to scan for motifs, file needs to be split into individual sequences
inputseq = open("protein_sequences/input.fa", "r")
content = inputseq.read()
allsequences = content.split('>')
inputseq.close()
for sequence in allsequences[1:]:
	acc = sequence.split()[0] #retrieve the accession value
	if '|' in sequence[2:4]:#in those sequences that contain NCBI identifiers
		acc = acc.split('|')[1]
	current_file = open('protein_sequences/individual_sequences/' + acc + '.fa','w')
	current_file.write('>' + sequence)
	current_file.close()
#Scan prosite
#loop through files
directory = 'protein_sequences/individual_sequences/'
for filename in os.listdir(directory):
	file = os.path.join(directory, filename)
	os.system('patmatmotifs -sequence {} -sid1 {} -outfile prosite_motifs/{}.patmatmotifs -auto'\
	.format(file,file[39:],file[39:])) #file[39:] is the accession number

#Extract motifs, and write them down to a file
files = glob.glob('prosite_motifs/*')
for file in files:
	motifs = open("reports/motifs.txt", "a")
	elements = []
	with open(file,'r') as f:
		for line in f:
			pattern1 = 'Sequence'
			if re.search( pattern1, line):
				elements.append(file.split('/')[1][:-16]) #the accession value, will be wrong if starts with sp| etc. Patmatmotifs can't deal with fasta headers in this format, acc name will be wrong
			pattern2 = 'Motif'
			if re.search(pattern2,line):
				elements.append(line.split()[2]) # the motif
	for element in elements:
		motifs.write(element + '\t')
	motifs.write('\n')
print('Done. The motifs for each sequence can be found in reports/motifs.txt')

#Option to perform tmap analysis
if input("Do you want to predict and plot possible transmembrane segments in the sequences (tmap, small datasets might result in no hits)? If yes, press y:\n\t") =='y':
	os.system('tmap -sequences protein_sequences/aligned_sequences.fa -graph svg -outfile reports/report.tmap -auto && tmap -sequences protein_sequences/aligned_sequences.fa -graph x11 -outfile reports/report.tmap -auto')
	print('Done, report can be found in reports/report.tmap')
#iep analysis
if input('Would you like to calculate the isoelectric point of the chosen proteins? If yes, press y:\n\t') == 'y':
	os.system('iep -sequence protein_sequences/aligned_sequences.fa -outfile reports/report.iep')
#hmoment analysis
if input('Would you like to calculate the hydrophobic moment for the chosen proteins? If yes, press y\n\t') == 'y':
	os.system('hmoment -seqall protein_sequences/aligned_sequences.fa  -outfile reports/report.hmoment ')
#pepwindowall analysis
if input('Would you like to plot a superimposed hydropathy plot for the chosen sequences (not recommended for large datasets)? If yes, press y\n\t') == 'y':
	os.system('pepwindowall -sequences protein_sequences/aligned_sequences.fa -window {} -gxtitle="Residue Position" -gytitle="Hydropathy" -graph svg && pepwindowall -sequences protein_sequences/aligned_sequences.fa -window {} -gxtitle="Residue Position" -gytitle="Hydropathy" -graph x11 '.format(winsize,winsize))

#Summarising all the results, with dictionaries
acc_motif = {}
for line in open("reports/motifs.txt", "r"):
	line = line.strip() # remove newline
	motifs = ','.join(line.split('\t')[1:]) #Joining all the motifd that were found into comma-separated string
	acc_motif[line.split('\t')[0]] = motifs

#And making a dataframe with the results, columns will be empty if tests were not performed

df = pd.DataFrame(list(filtered_dict.items()))
df.columns = ['Accession_number','Species']
df["Motifs"] = ""
df["Isoelectric_point"] = ""
df["Max_Hmoment"] = ""
for index, row in df.iterrows():
	df.loc[index,'Motifs'] = acc_motif[df.loc[index,'Accession_number']]
	if df.loc[index,'Motifs'] == '':
		df.loc[index,'Motifs'] = 'NaN'

#Making a dict with accession, iep so that I can add it to the df at the correct places
try:
	acc_iep = {}
	for line in open("reports/report.iep","r"):
		if line.startswith('IEP'):
			acc = line.split(' ')[2]
		elif line.startswith('Isoelectric'):#The way the emboss reports are formated
			iep = line.strip().split(' ')[3]
			acc_iep[acc] = iep
	for index, row in df.iterrows():
		try:
			df.loc[index,'Isoelectric_point'] = acc_iep[df.loc[index,'Accession_number']] # for every row, I find the iep in the dict by using tht same rows acc as key
		except KeyError:
			df.loc[index,'Isoelectric_point'] = 'NaN' #For some reason emboss programs mess up the sequence id when the fasta format contains NCBI identifiers, i.e, 'sp|'
except FileNotFoundError:
		print('Isoelectric point was not calculated, skipping from final report')#No report means that theuser didn't want to use this program
		df['Isoelectric_point'] = 'NaN'
#same concept as before, this time with the hmoment instead
try:
	acc_hmoment = {}
	for line in open('reports/report.hmoment','r'):
		if line.startswith('HMOMENT'):
			acc = line.split(' ')[2]
		elif line.startswith('Window'):
			hmoment = line.strip().split(' ')[-1]
			acc_hmoment[acc] = hmoment
	for index, row in df.iterrows():
		try:
			df.loc[index,'Max_Hmoment'] = acc_hmoment[df.loc[index,'Accession_number']]
		except KeyError:
			df.loc[index,'Max_Hmoment'] = 'NaN'

except FileNotFoundError:
		print('Hydrophobic moment was not calculated, skipping from final report')
		df['Max_Hmoment'] = 'NaN'
#Save the df
df.to_csv('output.tsv',sep='\t')
os.system('rm -fr acc_spp.txt && rm -fr clustalo_tree') # Unwanted files
