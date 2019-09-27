#!/bin/bash

#if [[ $(uname) == Darwin ]]
#then
#	split="$HOME/FFPEcapseq/bin/gsplit"
#else 
#	split=$(command -v split)
#fi

if [[ -z "$1" ]]
then 
	echo "Requires an argument" > /dev/stderr
	exit 1
fi

common="$1"
umi="$2"
cores="$3"

for x in "$common"*.*
	do if [[ "$x" =~ txt ]]
		then 
		echo mv "$x" "${x/.txt/.fastq}" > /dev/stderr
		mv "$x" "${x/.txt/.fastq}"

	fi
done

for x in "$common"*.*
	do if [[ "$x" =~ .gz$ ]] 
			then 
			if [[ ! -e Concatamer_summary.tab ]] || [[ $umi == True && ! -e UMI_"$x" && ! -e UMI_${x/.gz/} ]] || [[ ! $umi == True && ! -e Trimmed_"$x" && ! -e Trimmed_${x/.gz/} ]]
			then
			echo "Unzipping $x" > /dev/stderr
			if [[ -e ${x/.gz/} ]]
				then continue
			fi
			[[ $(type pigz 2>/dev/null) ]] && pigz -d -p "$cores" -c "$x" > "${x/.gz/}" || gunzip -c "$x" > "${x/.gz/}"
		else
			echo "Fastq already modified for $x" > /dev/stderr
		fi
	fi
done

if [[ ! -e Concatamer_summary.tab ]]
	then
	echo "Counting concatamers and rev primer" > /dev/stderr
	python ~/FFPEcapseq/bin/MP_concatamer.py "$cores" "$common"*.fastq 
	sort -Vk 1,1 Concatamer_summary.tab  > tmp; mv tmp Concatamer_summary.tab

	awk '
	BEGIN{
		PROCINFO["sorted_in"]="@val_num_desc"
	}
	{
		a[$1"_"$2"_"$3]+=$4
	}
	END{
		n=0
		for (x in a) {
			print x,a[x]
			n++
			if (n>1000) 
				break
		}
	}' Contamination_sequences.txt > tmp; mv tmp Contamination_sequences.txt

	echo "Contamination Done" > /dev/stderr
else 
	echo "Contamination calculations already done" > /dev/stderr
fi




retvalue=""

if [[ $umi == False ]]
	then
	umi_trim () {
					awk '
						NR%4==1{
							name=$1
						}
						NR%4==2{
							seq=substr($1,10)
							sub("^G{0,2}","",seq)
						}
						NR%4==3{
							p=$1
						}
						NR%4==0{
							seq_q=substr($1,(1+length($1)-length(seq)))
							printf "%s\n%s\n%s\n%s\n",name,seq,p,seq_q > "Trimmed_"FILENAME}END{print FILENME" done"> "/dev/stderr"
						}' "$1"
				}
	echo "UMI == False, Trimming AGGG" > /dev/stderr
	for x in "$common"*.fastq
		do 
		if [[ "$x" == ^Trimmed ]]
			then continue	
		else 
			echo trimming AGGG > /dev/stderr
			#pieces=$(echo "$cores" | awk '{if ($1<4) print 1; else print $1 - $1%4}')
			#lines=$(wc -l $x | awk -v p="$pieces" '{print $1/p}')
			#"$split" -l "$lines" "$x" split_
			awk '
                    BEGIN{
                        n=0
                        core='$cores'
                    }
                    ENDFILE{ if (FNR==NR)
                        cores=core-1
                        f=sprintf("%d",(FNR/4)/cores)
                    }
                    FNR==NR{next}
                    {
                        print $x > "split_"n
                    }
                    FNR%(f*4)==0{
                        n++
                    }
            ' $x $x
			for chunk in split_*
				do
					echo umi_trim "$chunk" >> commands.txt
			done
			export -f umi_trim
			~/FFPEcapseq/bin/queue commands.txt "$cores"
			cat Trimmed_split_* >> Trimmed_"$x"
			rm Trimmed_split*
			rm commands.txt
			rm split*
			retvalue="Trimmed"
		fi
	done
	for x in "$common"*.fastq
		do \rm -r "$x"
	done

elif [[ -z $umi || $umi == True ]]
	then 
	echo "UMI == True, trimming and compressing" > /dev/stderr
	for x in *"$common"*.fastq
		do
		if  [[ "$x" =~ ^UMI || -e "UMI_$x" ]]
			then 
			if [[ ! $common =~ ^UMI ]]
			then 
				retvalue="UMI"
			fi
			echo "Already Treated $x" > /dev/stderr
			continue 
		#This step assumes a 9bp UMI and an AGGG on 5' end
		else 
			echo "UMI demultiplexing..." > /dev/stderr
			[[ -e commands.txt ]] && rm commands.txt
			awk '
					BEGIN{
						n=0
						core='$cores'
					}
					ENDFILE{ if (FNR==NR)
						cores=core-1
						f=sprintf("%d",(FNR/4)/cores)
					}
					FNR==NR{next}
					{
						print $x > "split_"n
					}
					FNR%(f*4)==0{
						n++
					}
			' $x $x
			#pieces=$(echo "$cores" | awk '{if ($1<4) print 1; else print $1 - $1%4}')
			#lines=$(wc -l $x | awk -v p="$pieces" '{printf "%i",$1/p}')
			#"$split" -l "$lines" "$x" split_
			for chunk in split_*
			do
				umi_trim () {
					awk '
						NR%4==1{
							name=$1
							}
						NR%4==2{
							seq=substr($1,14) #Assumes a 9bp UMI and AGGG on 5 prime end
							sub("^G{0,2}","",seq) #strips out extra Gs from Template switching
							UMI=substr($1,0,9) #Assumes a 9bp UMI
							}
						NR%4==3{
							p=$1
							}
						NR%4==0{
							seq_q=substr($1,(1+length($1)-length(seq)))
							umi_q=substr($1,0,9) #Assumes a 9bp UMI
							printf "%s_%s_%s\n%s\n%s\n%s\n",name,UMI,umi_q,seq,p,seq_q > "UMI_"FILENAME}END{print FILENAME" done" > "/dev/stderr"
						}' "$1"
				}
				echo umi_trim "$chunk" >> commands.txt
			done
			export -f umi_trim
			~/FFPEcapseq/bin/queue commands.txt "$cores"
			cat UMI_split_* >> UMI_"$x" 
			rm UMI_split_*
			rm commands.txt
			rm split_*
			retvalue="UMI"
		fi
	done
	for x in "$common"*.fastq;
		do [[ ! $common =~ ^UMI ]] && rm -r "$x"
	done
fi

echo "Name change=$retvalue" > /dev/stderr

echo $retvalue

