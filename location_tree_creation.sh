#! /bin/bash

##############
#### Help ####
##############
Help()
{
	# Display Help
	echo "This script create a phylogenetic tree by taking as input a location name, the name of the used amplicon, a taxonomy and a sequence table."
	echo
	echo "Syntax: location_tree_creation.sh [-t|s|a|u|f|r|g|l|o|h]"
	echo
	echo "Options:"
	echo "-n	name of location you want to study (mandatory)"
	echo "-t	taxonomy table (mandatory)"
	echo "-s	occurence table (mandatory)"
	echo "-u	unwanted taxa list"
	echo "-a	amplicon name (mandatory)"
	echo "-f	forward primer sequence"
	echo "-r	reverse primer sequence"
	echo "-g	gene name (if different than amplicon name)"
	echo "-l	maximum length of the amplicon (advice: length of the amplicon for outgroup + 100bp)"
	echo "-o	outgroup fasta file already formatted"
	echo "-h	display this help message."
}

######################
#### Main program ####
######################

###############################
#### Process input options ####
###############################

while getopts hn:t:s:a:f:r:g:l:o: option
do
	case $option in
		(h) # display Help
			Help
			exit;;
		(n) # enter location name
			location=$OPTARG;;
		(t) # enter taxonomy table name
			tax_table=$OPTARG;;
		(s) # enter sequence table name
			occ_table=$OPTARG;;
		(u) # enter unwanted taxa list
			unwanted_taxa=$$OPTARG;;
		(a) # enter amplicon name
			amplicon=$OPTARG;;
		(f) # forward primer sequence
			primer_forward=$OPTARG;;
		(r) # reverse primer sequence
			primer_reverse=$OPTARG;;
		(g) # gene name
			gene_name=$OPTARG;;
		(l) # maximum length of the amplicon
			length_amp=$OPTARG;;
		(o) # outgroup
			outgroup=$OPTARG;;
		(\?) # Invalid option
			echo -e "\nError: Invalid option\n\n"
			Help
			exit;;
	esac
done

# Load conda environment
source activate tree-creation

echo -e "Checking the parameters given by the user to be sure that all needed informations are provided.\n"

# Checking if family name available

if [ "$location" == "" ]
then
        echo -e "Location name not detected. Please provide a location name to start the analysis."
        exit 1
else
        echo -e "Location name detected: $location"
fi


# Checking if taxonomy table available

if [ "$tax_table" == "" ]
then
	echo -e "Taxonomy table not detected. Please provide a file containing taxonomy information of the ASVs."
	exit 1
else
	echo -e "Taxonomy table detected: $tax_table"
fi

# Checking if sequence table available

if [ "$occ_table" == "" ]
then
        echo -e "Occurence table not detected. Please provide a file containing the occurence of each ASV."
	exit 1
else
        echo -e "Sequence table detected: $occ_table"
fi

# Checking if unwanted taxa list available

if [ "$unwanted_taxa" == "" ]
then
        echo -e "You didn't provide the Unwanted taxa list, the default list will be used to exclude some taxa from the phylogenetic tree. If you want to use your own list please provide a file name with the -u parameter."
		unwanted_taxa=removed_organisms.txt
else
        echo -e "Sequence table detected: $unwanted_taxa"
fi

# Checking if amplicon name provided and if primer sequences available if primer sequence not in list

if [ "$amplicon" == "" ]
then
	echo -e "No amplicon name detected. Please provide the name of the amplicon you used to start the analysis."
	exit 1
else
	if [ "$primer_forward" == "" ] && [ "$primer_reverse" == "" ]
	then
		grep_amp=`grep -o -P "^$amplicon\t" list_primers.tsv`
		if [ -z $grep_amp ]
		then
			echo -e "You didn't provide the forward and reverse primers of your amplicon and the amplicon name you provided is not in the default list. Please provide the primer sequences or use an amplicon name from the list."
			exit 1
		else
			echo -e "You didn't provide the forward and reverse primers of your amplicon but the amplicon name you provided is in the default list. Default primers from the list will be used for the trimming of the reference sequences.\n"
			primer_forward=`grep "$amplicon" list_primers.tsv | cut -f2`
			primer_reverse=`grep "$amplicon" list_primers.tsv | cut -f3`
			if [ "$gene_name" == "" ]
			then
				echo -e "You didn't provide the real gene name but the amplicon name you provided is in the default list. Default gene name from the list will be used to retrieve the reference sequences.\n"
				gene_name=`grep "$amplicon" list_primers.tsv | cut -f4`
			else
				echo -e "Gene name detected. The following gene name will be used to retrieve the reference sequences: $gene_name"
			fi
			if [ "$length_amp" == "" ]
                        then
                                echo -e "You didn't provide the maximum length of the amplicon but the amplicon name you provided is in the default list. Default maximum length from the list will be used to remove too long sequences before alignment.\n"
                                length_amp=`grep "$amplicon" list_primers.tsv | cut -f5 | tr -d '\r'`
                        else
                                echo -e "Maximum length of the amplicon detected. The following length will be used to remove too long sequences before the alignment: $length_amp"
                        fi
			if [ "$outgroup" == "" ]
                        then
                                echo -e "You didn't provide the fasta file of the outgroup but the amplicon name you provided is in the default list. Default outgroup from the list will be used.\n"
				outgroup=./outgroup_sequences/petromyzon_marinus_${amplicon}.fa
                        else
                                echo -e "Fasta file of the outgroup detected. The following outgroup will be used: $outgroup"
                        fi
		fi
	elif ["$primer_forward" == "" ] || [ "$primer_reverse" == "" ]
	then
		grep_amp=`grep -o -P "^$amplicon\t" list_primers.tsv`
                if [ -z $grep_amp ]
                then
                        echo -e "You only provided one primer for your amplicon and the amplicon name you provided is not in the default list. Please provide both forward and reverse primer sequences or use an amplicon name from the list."
                        exit 1
                else
                        echo -e "You only provided one primer for your amplicon but the amplicon name you provided is in the default list. Default primers from the list will be used for the trimming of the reference sequences.\n"
                        primer_forward=`grep "$amplicon" list_primers.tsv | cut -f2`
                        primer_reverse=`grep "$amplicon" list_primers.tsv | cut -f3`
			if [ "$gene_name" == "" ]
                        then
                                echo -e "You didn't provide the real gene name but the amplicon name you provided is in the default list. Default gene name from the list will be used to retrieve the reference sequences.\n"
                                gene_name=`grep "$amplicon" list_primers.tsv | cut -f4`
                        else
                                echo -e "Gene name detected. The following gene name will be used to retrieve the reference sequences: $gene_name"
			fi
			if [ "$length_amp" == "" ]
			then
				echo -e "You didn't provide the maximum length of the amplicon but the amplicon name you provided is in the default list. Default maximum length from the list will be used to remove too long sequences before alignment.\n"
				length_amp=`grep "$amplicon" list_primers.tsv | cut -f5 | tr -d '\r'`
			else
				echo -e "Maximum length of the amplicon detected. The following length will be used to remove too long sequences before the alignment: $length_amp"
			fi
			if [ "$outgroup" == "" ]
                        then
                                echo -e "You didn't provide the fasta file of the outgroup but the amplicon name you provided is in the default list. Default outgroup from the list will be used.\n"
                                outgroup=./outgroup_sequences/petromyzon_marinus_${amplicon}.fa
                        else
                                echo -e "Fasta file of the outgroup detected. The following outgroup will be used: $outgroup"
                        fi
                fi
	else
		echo -e "The following forward ${primer_forward} and reverse ${primer_reverse} sequences will be used to trim the NCBI reference sequences.\n"
		if [ "$gene_name" == "" ]
		then
			if [ -z $grep_amp ]
			then
				echo -e "You didn't provide the real gene name and the amplicon name you provided is not in the default list. Please provide a gene name to be able to retrieve the reference sequences.\n"
				exit 1
			else
                                echo -e "You didn't provide the real gene name but the amplicon name you provided is in the default list. Default gene name from the list will be used to retrieve the reference sequences.\n"
                                gene_name=`grep "$amplicon" list_primers.tsv | cut -f4`
			fi
		else
			echo -e "Gene name detected. The following gene name will be used to retrieve the reference sequences: $gene_name\n"
		fi
		if [ "$length_amp" == "" ]
		then
			if [ -z $grep_amp ]
			then
				echo -e "You didn't provide the maximum length of the amplicon and the amplicon name you provided is not in the default list. Please provide a maximum length of the amplicon to be able to remove too long sequences before alignment.\n"
				exit 1
			else
				echo -e "You didn't provide the maximum length of the amplicon but the amplicon name you provided is in the default list. Default maximum length from the list will be used to remove too long sequences before alignment.\n"
				length_amp=`grep "$amplicon" list_primers.tsv | cut -f5 | tr -d '\r'`
			fi
		else
			echo -e "Maximum length of the amplicon detected. The following length will be used to remove too long sequences before the alignment: $length_amp\n"
		fi
		if [ "$outgroup" == "" ]
                then
                        if [ -z $grep_amp ]
                        then
                                echo -e "You didn't provide the fasta file of the outgroup and the amplicon name you provided is not in the default list. Please provide the outgroup file to be able to plot your tree.\n"
                                exit 1
                        else
                                echo -e "You didn't provide the fasta file of the outgroup but the amplicon name you provided is in the default list. Default outgroup from the list will be used.\n"
                                outgroup=./outgroup_sequences/petromyzon_marinus_${amplicon}.fa
                        fi
                else
                        echo -e "Fasta file of the outgroup detected. The following outgroup will be used: $outgroup"
                fi
	fi
fi

echo -e "Starting tree creation for the $location location\n"

mkdir $location
cd $location

# Create file containing all ASV sequences of this sample
echo -e "Starting to retrieve all the ASV sequences for this location...\n----------\n"

	## Get column numbers of interesting fields
	DNA_col=`awk -v RS='\t' '/DNA_sequence/{print NR; exit}' ../${tax_table}`
	spe_col=`awk -v RS='\t' '/scientificName/{print NR; exit}' ../${tax_table}`
	rank_col=`awk -v RS='\t' '/taxonRank/{print NR; exit}' ../${tax_table}`

	## Retrieve ASV sequences
	grep "^asv" ../${tax_table} | sed 's/\t/|/g' | while read line; do asv=`echo $line | cut -d'|' -f1`; species=`echo $line | cut -d'|' -f${spe_col} | sed 's/ /_/g'`; echo -e ">${asv}_${species}" >> ${location}_ASV.fa; seq=`echo $line | cut -d'|' -f${DNA_col}`; echo $seq >> ${location}_ASV.fa; done

	## Remove duplicated sequences
	seqkit rmdup -s < ${location}_ASV.fa > ${location}_ASV_nodup.fa

	## Transform sequences to get DNA on one line for each ASV
	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_ASV_nodup.fa > ${location}_ASV_nodup_correct.fa

# Remove the ASVs identified as human
echo -e "Remove ASVs identified as human...\n----------\n"

	while read line
	do
		grep_res=`echo $line | grep -v "asv." | head -n1 | cut -d' ' -f1`
		if [ -z $grep_res ] # if lines contains asv. do the following
		then
			species=`echo $line | cut -d'_' -f2-` # get species name
			name=`echo $line | cut -d'>' -f2` # get the complete sequence name
			if [[ $species == "Homo_sapiens" || $species == "Incertae_sedis" ]] # check if species is human or assigned incertain
			then
				num=`grep -n "$name" ${location}_ASV_nodup_correct.fa | cut -d":" -f1` # if it is retrieve line number of the asv sequence
				other=$(( $num + 1 ))
				sed -i ''$num','$other'd' ${location}_ASV_nodup_correct.fa # and remove these lines
			fi
		fi
	done < ${location}_ASV_nodup_correct.fa

	## Remove temporary files
	rm ${location}_ASV.fa ${location}_ASV_nodup.fa
	mv ${location}_ASV_nodup_correct.fa ${location}_ASV.fa

# Retrieve NCBI sequences to add known species to the tree
echo -e "Starting to retrieve the NCBI sequences corresponding to the genus or species identified with blast and vsearch...\n----------\n"

	## Retrieve the accession numbers
	touch temp_spe_file.txt
	while read line
	do
		grep_res=`echo $line | grep -v "asv."`
		if [ -z $grep_res ]
		then
			asv=`echo $line | cut -d'_' -f1 | cut -d'>' -f2`
			rank=`grep -P "^$asv\t" ../${tax_table} | cut -f${rank_col}`
			genus=`echo $line | cut -d'_' -f2`
			grep_gen=`grep "$genus" temp_spe_file.txt | head -n1 | cut -d' ' -f1`
			if [ -z $grep_gen ]
			then
				if [ "$rank" == "genus" ]
				then
					grep "$genus " ../all_nt_db_acc.txt | grep -E 'genome|${gene_name}' | grep "mitochon" | sort -u -t' ' -k3 | cut -d' ' -f1 | cut -d'>' -f2 >> ${location}_accessions.txt
					echo $genus >> temp_spe_file.txt
				elif [ "$rank" == "species" ]
				then
					# spe=`echo $line | cut -d'_' -f2- | sed 's/_/ /g'`
					grep "$genus " ../all_nt_db_acc.txt | grep -E 'genome|${gene_name}' | grep "mitochon" | sort -u -t' ' -k3 | cut -d' ' -f1 | cut -d'>' -f2 >> ${location}_accessions.txt
					echo $genus >> temp_spe_file.txt
				fi
			fi
		fi
	done < ${location}_ASV.fa

	## Remove temporary files
	rm temp_spe_file.txt

	## Download sequences of these accession numbers
	echo -e "Downloading the sequences from NCBI...\n----------\n"
	sort ${location}_accessions.txt | uniq > ${location}_acc_nodup.txt
	mkdir ncbi_acc_files
	while read line
	do
		acc=`echo $line`
		ncbi-acc-download --format fasta --out ncbi_acc_files/${acc}.fasta $acc
	done < ${location}_acc_nodup.txt
	cat ncbi_acc_files/*.fasta > ${location}_organisms.fasta

	## Cleaning temporary files
	rm -r ncbi_acc_files/
	rm ${location}_accessions.txt

	## Put all fasta sequences on one line
	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_organisms.fasta > ${location}_organisms_1l.fa

	## Clean
	rm ${location}_organisms.fasta
	mv ${location}_organisms_1l.fa ${location}_organisms.fa

# Remove unwanted taxa and put all sequences (organisms and asv) in one file
echo -e "Removing unwanted taxa like Hominidae and Putting all the sequences in one file...\n----------\n"

	## Remove all unwanted sequences
	grep -f ../$unwanted_taxa ${location}_organisms.fa | sed 's/^>//g' > idlist.txt
	seqkit grep -vi --pattern-file idlist.txt -n ${location}_organisms.fa > ${location}_organisms_clean.fa
	sed -i 's/PREDICTED: /PREDICTED:/g' ${location}_organisms_clean.fa
	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_organisms_clean.fa > ${location}_organisms_clean_1l.fa

	## Remove temporary files
	rm idlist.txt ${location}_organisms_clean.fa
	mv ${location}_organisms_clean_1l.fa ${location}_organisms_clean.fa

	## Rename the sequences to keep only the organism genus and species
	grep "^>" ${location}_organisms_clean.fa > names_organisms.txt
	while read line
	do
		name=`echo $line | cut -d' ' -f2-3`
		sed -i "s|$line|>${name}|" ${location}_organisms_clean.fa > sed_output.txt 2>&1
	done < names_organisms.txt

	## Add asv sequences to the same file and remove temporary files
	cat ${location}_ASV.fa >> ${location}_organisms_final.fa
	cat ${location}_organisms_clean.fa >> ${location}_organisms_final.fa
	rm names_organisms.txt

# Apply cutadapt on all the sequences
echo -e "Applying cutadapt to the sequences in order to keep only the fragment of interest...\n----------\n"

	## Detect forward primer
	cutadapt -g $primer_forward -o ${location}_final-short.fa -e 4 -j 1 ${location}_organisms_final.fa > ${location}_cutadapt1.log 2>&1

    ## Detect reverse primer
    cutadapt -a $primer_reverse -o ${location}_final-short2.fa -e 4 -j 1 ${location}_final-short.fa > ${location}_cutadapt2.log 2>&1

	## Clean
    rm ${location}_final-short.fa

# Refine sequence file to make it clean
echo -e "Refining the sequences and changing sequence's names...\n----------\n"

	# remove too long sequences which were not cutted with cutadapt (no amplicon in seq?)
	seqkit seq -M $length_amp ${location}_final-short2.fa | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > ${location}_final-short2_correct.fa

	# remove duplicated organisms
	echo -e "Removing duplicate sequences and duplicated names...\n----------\n"
	seqkit rmdup -n < ${location}_final-short2_correct.fa > ${location}_final-short_nodup.fa

	## Put the sequences on one line each
	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_final-short_nodup.fa > ${location}_final-short2_nodup_correct.fa

	## Remove temporary files
	rm ${location}_final-short2.fa ${location}_final-short_nodup.fa
	mv ${location}_final-short2_nodup_correct.fa ${location}_final-short.fa

# Add the outgroup to the file
	echo -e "Adding the outgroup...\n----------\n"
	cat ../$outgroup >> ${location}_final-short.fa

# Infer the phylogenetic tree
	## Align the sequences with Mafft
	echo -e "Aligning the sequences with Mafft...\n----------\n"
	mafft --thread -1 --preservecase --adjustdirection ${location}_final-short.fa > ${location}.aln.fa 2> ${location}_mafft.log

	## Replace spaces by underscores in the organisms names
	sed -i 's/ /_/g' ${location}.aln.fa

	## Removed special characters in the names
	sed -i 's/://g' ${location}.aln.fa
	sed -i 's/_R_asv./asv./g' ${location}.aln.fa
	sed -i 's/_R_//g' ${location}.aln.fa

	## Curate the alignments with Gblocks
	echo -e "Refining the alignment with Gblocks...\n----------\n"
	Gblocks ${location}.aln.fa -t=d -b4=5 -b5=h -e=.gb > ${location}_gblock.log 2>&1

	## Create the tree with Fasttree
	echo -e "Creating the phylogenetic tree with Fasttree...\n----------\n"
	FastTreeMP -gtr -nt ${location}.aln.fa.gb > ${location}.nwk 2>${location}_FastTree.log

# Create the file to retrieve the number of sampling site for each ASV
	echo -e "Creating the table file with quantitative information on the ASV (number of site detecting the ASV)...\n----------\n"

	## Retrieve all the sequence names
	grep "^>" ${location}.aln.fa.gb | sed 's/>//g' > ${location}_quant_names.tsv

	## For each ASV retrieve the relative abundance
	### occurence table file: ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Occurrence_table.tsv
	while read line
	do
		asv=`echo $line | cut -d'_' -f1`
		sites=`grep "${asv}_EE" ../${occ_table} | wc -l`
		echo "$sites" >> ${location}_quant_num.tsv
	done < ${location}_quant_names.tsv

	# Assemble names and percentage in the same file and give headers to the table
	paste ${location}_quant_names.tsv ${location}_quant_num.tsv > ${location}_quant.tsv
	sed -i '1itaxa\tquantSites' ${location}_quant.tsv

	# Remove temporary files
	rm ${location}_quant_names.tsv ${location}_quant_num.tsv

# Produce a file with the tree visualization
	echo -e "Creating the tree image...\n----------\n"

	## Change to conda environment for tree construction
	conda activate R-tree-creation

	## call the R script
	path=`pwd`
	Rscript ../plot_tree_location.r ${path} ${location}.nwk ${location}_quant.tsv ${location} ${location}_root.pdf

	## Remove temporary files
	rm Rplots.pdf

echo -e "Analysis finished! You can find the tree image in the folder ${location}/ under the name ${location}_root.png"
