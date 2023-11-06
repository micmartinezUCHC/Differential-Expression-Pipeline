#!/usr/bin/bash

#--------Differential Expression Pipeline (DEPl) version 3.0
#-----Mike Martinez
#---Rosenberg Lab UCHC

#-----Print pipeline banner
figlet -f doom 'DEPl v3.0'
echo "Mike Martinez 2023"
echo "Rosenberg Lab UCHC"

#-----Prompt the user to select a runMode
echo "Select Run Mode:"
echo "1: Single Factor Design DESeq2"
echo "2: Single Factor Design DESeq2 + ClusterProfiler GSEA"
read mode
clear

#-----Ensure run mode is of either 1 or 3
if [ "$mode" != "1" ] && [ "$mode" != "2" ]; then
	echo "Error: Invalid run mode selection. Must be of either 1 or 2"
	exit 1
fi

#-----Set base output directory, date stamped
resultDir="/Users/mikemartinez/Desktop/PREVENT/Pipeline_Results/$(date +'%B_%d_%Y_%H:%M')"

#-----Check if the baseoutput directory exists
if [ -d "$resultDir" ]; then
	echo "$resultDir already exists"
else
	mkdir -p "$resultDir"
	echo "Directory $resultDir successfully created"
fi

#-----Prompt the user to define the parameters for the run
echo "Enter the absolute path to the raw counts file: "
read counts
clear

#-----Check that the file path specified exists and is not empty
if [ -e "$counts" ]; then
	echo ""
else
	echo "Error Counts file does not exist. Please check the path."
	exit 2
fi

echo "Enter a name for the comparison: "
read name
clear

echo "Enter reference level for comparison: "
read ref
clear

echo "Enter number of observations of reference level: "
read nRef
clear

echo "Enter numerator level for comparison: "
read treatment
clear

echo "Enter number of observations of numerator level: "
read nTreatment
clear

echo "rno or mmu (Do not input response with any capital letters)"
read organism
clear

#-----Check that a valid organism code was specified. 
if [ "$organism" != "rno" ] && [ "$organism" != "mmu" ]; then
	echo "Error: Invalid organism code. DEG pipeline v3.0 is compatible with either mmu or rno"
	exit 3
fi

#-----Det DESeq2 results output folder within the resultDir
DESeq2OP="$resultDir/$name"


echo -e "Running DEPl in Mode $mode \n \n"
echo -e "Parameters:\n Run Name: $name \n Denominator (reference) level: $ref \n # Observations denominator: $nRef \n Numerator level: $treatment \n # Observations numerator: $nTreatment \n Organism: $organism"
echo "#------------------------------#"

#-----Run RScripts dependent on mode selection
Rscript automated_DESeq2_v3.R --counts "$counts" --name "$name" --ref "$ref" --n "$nRef" --treatment "$treatment" --m "$nTreatment" --out "$DESeq2OP" --org "$organism" 


if [ "$mode" = "2" ]; then

	#-----Set the variables for the GSEA analysis
	DEG_file="${DESeq2OP}/DESeq2_Results/${name}_All_DEGs.csv"
	echo -e "\n"
	echo -e "#----- Beginning GSEA -----#\n"

	#-----Check that the DEG file exists and is not empty
	if [ -e "$DEG_file" ]; then
		if [ -s "$DEG_file" ]; then
			echo "DEG file exists and is not empty!"
		else
			echo "DEG file exists, but is empty."
		fi
	else
		echo "DEG does not exist"
		exit 4
	fi

	#-----GSEA output folder
	gseaOP="$DESeq2OP/GSEA_Results" 
	echo "$gseaOP"

	#-----Run the GSEA Rscript
	Rscript automated_GSEA_v3.R --DEGs "$DEG_file" --name "$name" --ref "$ref" --treatment "$treatment" --out "$gseaOP" --org "$organism" 
fi

echo -e "\n \n \n DEPl completed!"










 
