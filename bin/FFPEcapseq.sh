#!/bin/bash
#SBATCH -t 48:00:00 -N 1 -n 2 -J FFPECAPSEQ --account=gertz-kp --partition=gertz-kp

source "$HOME"/FFPEcapseq/bin/arg_parser.sh

bash "$HOME"/FFPEcapseq/bin/requirment_check.sh


[[ $? -eq 1 ]] && exit

usage () {
cat << help_message
USAGE:  ${0##*/} -c common_name [OPTIONAL ARGUMENTS]

Argument options:

Mandatory:
-c	Common name that all experiments share
	or a space separated list of fastqs belonging
	to the same experiment

-g 	eneter a genome name.
	options hg19, hg38, mm10, Scerv

Optional:
-u	UMIs present, Default:True, 
	set to False if no UMIs

-p	# of cores to use for parallel processes
	type max for all physical cores, else enter 
	number of cores to use.  Default:1

-v	Set to True to see what the program is doing in 
	more detail.  Default:""

-get	y to get Genomes from FFPEcapseq repository 

help_message
}

get_flags verbose~-v umi~-u:True genome=-g common=-c get_genomes~-get cores~-p:1

if [[ $verbose ]]
	then 
	echo "Setting verbose mode"
	verbose=True
	set -v 
fi

if [[ ! $SLURM_CPUS_ON_NODE ]]
    then
    echo submitting sbatch $0 -c $common -u $umi -p $cores -g $genome 
	sbatch $0 -c $common -u $umi -p $cores -g $genome
    exit
fi


if [[ $(uname -a) =~ Darwin ]]
	then 
	available_cores=$(system_profiler SPHardwareDataType | awk '$0~"Total Number of Cores:"{print $NF}')
	else
	available_cores=$(lscpu | awk '$0~"Core\\(s\\) per socket:"{C=$NF}$0~"Socket\\(s\\):"{S=$NF}END{print S*C}')
fi

if [[ $(echo $cores | awk '{print tolower($1)}') == "max" || $cores -gt $available_cores ]]
	then 
	cores=$available_cores
elif [[ $cores -gt 1 && $cores -lt $available_cores ]]
	then cores=$cores
fi

cat << settings

Running with the following settings:

Fastq Files: 
$(ls $common*)

Genome: $genome

Cores: $cores

UMI: $umi

Verbose: $verbose

settings
[[ $verbose ]] && verbose="-v"
echo common="$common" > .settings
echo cores="$cores" >> .settings
echo genome="$genome" >> .settings
echo umi="$umi" >> .settings
echo verbose="$verbose" >> .settings


startdir=$(pwd)
while [[ $get_genomes == "Y" || $get_genomes == "y" ]]
        do PS3='Select Genome: '
        cd ~/FFPEcapseq || (echo 'No ~/FFPEcaseq folder!'; exit 1)
        list=($(curl -s http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/Genomes))
	echo -e "\nGenome Options:"
	select x in "${list[@]}"
		do if [[ "$x" == "0" ]]
			then break
		elif [[ -z $x ]]
			then echo "Invalid Option"
		else 
			if [[ ! -d $x ]]
                                then curl -O http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/$x\.tar.gz
                                tar -xzf $x.tar.gz
                                rm $x.tar.gz
				break
                        else
                                echo -e "\n$x is already present.\nIf you want to replace, please delete ~/FFPEcapseq/$x before trying to install.\n"
				break
                        fi
		fi
	done
	read -p "Do you want to get more bowtie indexes? (y/n):" get_genomes
        if [[ $get_genomes != "y" && $get_genomes != "Y" ]]
                then cd $startdir || (echo "No folder $startdir !"; exit 1)
        fi
done

if [[ ! -z $get_genomes ]]
	then echo 
	read -p "Continue with FFPEcapseq pipeline?(y/n):" response
	if [[ $repsone != "Y" || $response != "y" ]]
		then exit
	fi
fi

if [[ ! -e ~/FFPEcapseq/$genome ]]
	then indexes=$(curl -s http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/)
	if [[ $indexes =~ $genome\.tar.gz ]]
		then read -p "You don't have that genome in $HOME/FFPEcapseq.  Would you like to download it? (Y/N):" response
		if [[ $response == "Y" || $response == "y" ]]
			then startdir=$(pwd)
			cd ~/FFPEcapseq || (echo 'No ~/FFPEcaseq folder!'; exit 1)
			curl -O http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/$genome\.tar.gz
			tar -xzf $genome\.tar.gz
			rm $genome\.tar.gz
			cd $startdir || (echo "No folder $startdir !"; exit 1)
		else
			echo "You will need to build genome and mRNA indexes for your genome and put/link them to your $HOME/FFPEcapseq folder"
			echo "See README.txt for more information"
			exit 1
		fi
	else
		echo -e "I have not built the $genome indexes.  You will need to build genome and mRNA indexes for your genome and put/link them to your $HOME/FFPEcapseq folder.\n\nSee e README.txt for more information on indexes\n\nAvailable Genomes are:"
		curl -s http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/Genomes
		exit 1
	fi
fi

###############################
###  pretreatment of Fastqs ###
###############################

echo unzipping, renaming, extracting UMIs etc...
namechange=$(bash $verbose ~/FFPEcapseq/bin/pretreat.sh $common $umi $cores)
if [[ $? -ne 0 ]]
	then 
	echo pretreat failed
	exit 1
fi

echo pretreat finished


[[ -z $namechange ]] || common=$namechange\_$common
echo New common name=$common

###############################
###   Alignment of fastqs   ###
###############################
echo aligning > /dev/stderr
for x in $common*.fastq*
	do 
	name=$(echo ${x##*/}| sed 's/\..*//g')
	if [[ ! -e $name/$name\_sorted.sam || ! -e $name/$name\_sorted_refseq.sam ]]
		then echo missing some alignements for $name, running bowtie
		bash $verbose ~/FFPEcapseq/bin/bowtie.sh $x $genome $cores $umi
		if [[ $? -ne 0 ]]
			then
			echo pretreat failed
			exit 1
		fi
	fi
done

################################
###   awk analysis of sams   ###
################################
[[ -e commands.txt ]] && rm commands.txt
for x in $common*
	do if [[ -d $x ]]
		then echo "bash $verbose ~/FFPEcapseq/bin/awkanalysis.sh $x $genome $umi" >> commands.txt 
		if [[ $? -ne 0 ]]
			then
			echo pretreat failed > /dev/stderr
			exit 1
		fi
	fi
done

~/FFPEcapseq/bin/queue commands.txt $cores


rm commands.txt
echo -e '\a'
echo "making report"

#########################
###  Generate Report  ###
#########################

n=0

for x in $common*
do
        if [[ -d $x ]]
        then
                if [[ $n == 0 ]]
                        then
                        awk -v samp=$x -v ORS="\t" -v FS="=" 'BEGIN{print "SAMPLE"}$0~"="{print $1}' ./$x/$x\_summary.txt > $common\_final_summary.tab
                        echo "" >> $common\_final_summary.tab
                fi
        awk -v samp=$x -v ORS="\t" -v FS="=" 'BEGIN{print samp}$0~"="{print $2}' ./$x/$x\_summary.txt >> $common\_final_summary.tab
        echo "" >> $common\_final_summary.tab
        n=$n+1
        fi
done

if [[ $(uname) == "Linux" ]]
        then ~/FFPEcapseq/bin/sort -Vk 1,1 $common\_final_summary.tab > temp; mv temp $common\_final_summary.tab
elif [[ $(uname) == "Darwin" ]]
        then ~/FFPEcapseq/bin/gsort -Vk 1,1 $common\_final_summary.tab > temp; mv temp $common\_final_summary.tab
else
        sort -Vk 1,1 $common\_final_summary.tab > temp; mv temp $common\_final_summary.tab
fi



[[ -z $namechange ]] && nc="" || nc=$namechange\_
awk '
FNR==NR{
        a["'$nc'"$1]=$2"\t"$3"\t"$4
        next
        }
FNR!=NR{
        if (FNR==1)
                print $0"Concatemer %\tRev Primer %\tBoth %"
        else if ($1 in a)
                print $a[$1]
        else
                print $0
}
' Concatamer_summary.tab $common\_final_summary.tab > temp

mv temp $common\_final_summary.tab

for x in $common*
        do
        if [[ -d $x ]]
                then echo $x/$x\_geneexpression.tab >> geneexpressionlist.txt
        fi
done

echo FFPEcapseq pipeline Finished

echo Final report:
echo "$(pwd)"/"$common"\_final_summary.tab


