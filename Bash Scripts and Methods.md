## Membrane protein identification & collection of their protein sequences
#### To eliminate TM=0 & SP=Y & sort by no. of TM helices:
```bash
awk '$2!==0 && $3!=="Y"' *phobiuspredict.txt > *editedphobius.txt
sort -gk 2n -o *editedphobius.txt *editedphobius.txt
```
#### To get list of accession numbers for membrane proteins in text file:
```bash
awk -F"\ " '{print $1}' *editedphobius.txt > *TMproteins.txt
```
#### Remove lcl prefix & SEQUENCE_ID title:
```bash
sed -i 's/lcl|//;s/SEQUENCE_ID//g' *TMproteins.txt; sed 1d *TMproteins.txt
```
#### Make & run grep command to write list of sequence headers to parse in Python script to pair headers and sequences:
```bash
awk -v d="|" '{s=(NR==1?s:s d)$0}END{print s}' *TMproteins.txt > command.sh
sed -i 's/|/\\|/g' command.sh; sed -i "s/.*/grep '&' *.fa/" command.sh
chmod +x command.sh
./command.sh > TM_AN_proteins.txt
```
Python script `protein-sequences.ipynb` used to make a dictionary of the protein sequences, and create a text file header-sequence-pairs.txt using the list of sequence headers containing only the sequences of membrane proteins.
*This was initially opened with Jupyter notebooks, input filenames edited and run in the browser before using Atom text editor to edit input filenames and the command runipy was used (with Miniconda installed in Linux Terminal):*
```bash
runipy protein-sequences.ipynb
mv header-sequence-pairs.txt Bacteria_Name_translated_cds_mps.fa
```
#### Check number of lines equal in both files (in case grep command selected multiple sequence headers for one sequence ID):
```bash
cat TM_AN_proteins.txt | wc -l; cat *TMproteins.txt | wc -l
```

# Transmembrane protein information (number of helices, protein names)
#### To get number of proteins with each number of predicted TM helices:
```bash
awk '$2==1' *editedphobius.txt | wc -l
awk '$2==2' *editedphobius.txt | wc -l
# Repeated for all values of column 2 (up to max number of predicted TM helices)
```
#### To get a list of the protein sequence IDs and their protein names (according to genome annotation): 
```bash
sed 's/ \[/;[/g' TM_AN_proteins.txt > list_of_proteins.txt
names="$(awk -F";" '{print $1,$4}' list_of_proteins.txt | grep "protein=")"
names2="$(awk -F";" '{print $1,$5}' list_of_proteins.txt | grep "protein=")"
echo "$(echo -e "$names\n$names2" | sort -gk 1n)" > list_of_proteins.txt
```

# Identification of human non-homologous membrane proteins 
#### To make human membrane proteome local BLAST database from membrane protein sequences fasta file:
```bash
makeblastdb -in GRCh38.p13_translated_cds_mps.fa -out blastdb -parse_seqids -dbtype prot
```
#### To run local BLASTp search against human membrane proteome:
```bash
blastp -db blastdb -query Bacteria_Name_translated_cds_mps.fa -outfmt 0 -out Bacteria_name_results.txt -num_threads 4
blastp -db blastdb -query Bacteria_Name_translated_cds_mps.fa -outfmt 6 -out Bacteria_name_results_parse.txt -num_threads 4
```
#### To filter BLASTp search results by E Value (>0.001) & Percent Identity (<35%):
```bash
# For number of NO protein hits:
grep -c "No hits" *_results.txt

# For number of total protein hits:
awk -F"\t" '{print $1}' *_results_parse.txt | sort | uniq | wc -l

# To report how many homologous proteins and list homologous hits in text file:
awk -F"\t" '$3>35 && $11<1e-3' *_results_parse.txt > homologous_hits.txt
grep -f homologous_hits.txt *_results_parse.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l

# To collect list of all the non-homologous proteins (incl. no hits):
grep -v "$(awk -F"\t" '{print $1}' homologous_hits.txt | sort | uniq)" *TMproteins.txt | sort > non_hom_proteins.txt

# For total number of non-homologous proteins (incl. no hits):
echo "$(grep -v "$(awk -F"\t" '{print $1}' homologous_hits.txt)" *_results_parse.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l) + $(grep -c "No hit" *_results.txt)" | bc
```

# Identification of essential human non-homologous membrane proteins
#### To make essential gene local BLAST database from DEG amino acid sequences file:
```bash
makeblastdb -in DEG.aa -out DEG10 -parse_seqids -dbtype prot
```
#### To run local BLASTp search against Database of Essential Genes (DEG):
```bash
blastp -db DEG10 -query non_hom_proteins.fa -outfmt 6 -out non_hom_deg_results.txt -num_threads 4
```
#### To check DEG Blastp results:
•	awk -F"\t" '$3>35 && $11<1e-3' non_hom_deg_results* (for 9 matched species)
•	awk -F "\t" '$3>35 && $11<1e-4 && $12>100' non_hom_deg_results.txt | awk -F "\t" '{print $1}' | sort | uniq | wc -l
•	awk -F "\t" '$3>35 && $11<1e-4 && $12>100' non_hom_deg_results.txt | awk -F "\t" '{print $1}' | sort | uniq > unique_essential_proteins.txt
#### To copy python script and GCF fa to next folder and edit python script to include name of proteins (e.g. 6-span-TMproteins_uniq.txt or i-i-6-span...etc):
•	cp protein-sequences.ipynb GCF*_translated_cds.fa ../6-span-proteins/; cd ../6-span-proteins/; ls -l; atom protein-sequences.ipynb
#### To get fasta sequences for sequence alignment:
•	runipy protein-sequences.ipynb; mv header-sequence-pairs.txt ./6-span-proteins_uniq.fa

Protein Clusters
•	cd-hit -i combined-unique_essential_proteins.fa -o combined-unique_essential_proteins_cdhit_0.6 -c 0.6 -n 4 -d 0 -g 1
•	cd-hit -i combined-unique_essential_proteins.fa -o combined-unique_essential_proteins_cdhit_0.8 -c 0.8 -n 5 -d 0 -g 1
•	cd-hit -i combined-unique_essential_proteins.fa -o combined-unique_essential_proteins_cdhit_0.9 -c 0.9 -n 5 -d 0 -g 1
•	make_multi_seq combined-unique_essential_proteins.fa combined-unique_essential_proteins_cdhit_0.9.clstr multi-seq_0.9 5
•	sed 's/^$/-|/g' clustrs_AN.txt
•	csplit --digits=2 --quiet --prefix=outfile infile "/-|/+1" "{*}"

