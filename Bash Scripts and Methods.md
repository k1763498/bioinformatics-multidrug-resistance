## Membrane protein identification & collection of their protein sequences
To eliminate TM=0 & SP=Y & sort by no. of TM helices:
```bash
awk '$2!==0 && $3!=="Y"' *phobiuspredict.txt > *editedphobius.txt
sort -gk 2n -o *editedphobius.txt *editedphobius.txt
```
To get list of accession numbers for membrane proteins in text file:
```bash
awk -F"\ " '{print $1}' *editedphobius.txt > *TMproteins.txt
```
Remove lcl prefix & SEQUENCE_ID title:
```bash
sed -i 's/lcl|//;s/SEQUENCE_ID//g' *TMproteins.txt; sed 1d *TMproteins.txt
```
Make & run grep command to write list of sequence headers to parse in Python script to pair headers and sequences:
```bash
awk -v d="|" '{s=(NR==1?s:s d)$0}END{print s}' *TMproteins.txt > command.sh
sed -i 's/|/\\|/g' command.sh; sed -i "s/.*/grep '&' *.fa/" command.sh
chmod +x command.sh
./command.sh > TM_AN_proteins.txt
```
Python script used to make a dictionary of the protein sequences, and create a text file header-sequence-pairs.txt using the list of sequence headers containing only the sequences of membrane proteins.
This was initially opened with Jupyter notebooks, input filenames edited and run in the browser before using Atom text editor to edit input filenames and the command runipy was used (with Miniconda installed in Linux Terminal):
•	runipy protein-sequences.ipynb
•	mv header-sequence-pairs.txt Bacteria_Name_translated_cds_mps.fa
Check number of lines equal in both files (in case grep command selected multiple sequence headers for one sequence ID):
•	cat TM_AN_proteins.txt | wc -l; cat *TMproteins.txt | wc -l
