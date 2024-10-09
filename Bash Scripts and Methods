Membrane protein identification & collection of their protein sequences
To eliminate TM=0 & SP=Y & sort by no. of TM helices:
•	awk '$2!==0 && $3!=="Y"' *phobiuspredict.txt > *editedphobius.txt
•	sort -gk 2n -o *editedphobius.txt *editedphobius.txt
To get list of accession numbers for membrane proteins in text file:
•	awk -F"\ " '{print $1}' *editedphobius.txt > *TMproteins.txt
Remove lcl prefix & SEQUENCE_ID title:
•	sed -i 's/lcl|//;s/SEQUENCE_ID//g' *TMproteins.txt; sed 1d *TMproteins.txt
Make & run grep command to write list of sequence headers to parse in Python script to pair headers and sequences:
•	awk -v d="|" '{s=(NR==1?s:s d)$0}END{print s}' *TMproteins.txt > command.sh
•	sed -i 's/|/\\|/g' command.sh; sed -i "s/.*/grep '&' *.fa/" command.sh
•	chmod +x command.sh
•	./command.sh > TM_AN_proteins.txt
