#!/bin/bash
source ~/FFPEcapseq/bin/arg_parser.sh

usage (){
		echo -e "\nUsage:"
		echo `basename $0`" -f file.sam/bam [  -t output-type: (SAM or BAM) -o outfilename ] " 
}

get_flags infile=-f outfile~-o outtype~-t 


tempfile=`mktemp`
outtype=${outtype:="SAM"}
outfile=${outfile:="/dev/stdout"}
ftype=`file $infile`
if [[ $ftype =~ "/dev/stdin" ]]
then 
	echo -e "\nERROR!:\nThere are problems using this program with a pipe.\nPlease feed a file in through the -f flag"  >&2
	exit 2
fi

if [[ $ftype =~ "data" ]]
then 
	samtools view -h $infile  > $tempfile
else 
	cp $infile $tempfile
fi


awk '
function compare(str1,str2,max){
	split(str1,s1,"")
	split(str2,s2,"")
	if (length(s1)!=length(s2))
		return 0
	n=0
	for (x=1;x<=length(s1);x++){
		if (s1[x] != s2[x])
			n++
		if (n > max ) 
			return 0
		}
	return 1
}
BEGIN{
	FS="\t"
	OFS="\t"
	PROCINFO["sorted_in"]="@ind_str_asc"
}
NR==1{print}
$1~"^@"{
	print
	next
}

{	
	split($1,id,"_")
	data[$3,$4][id[2]]=$0
}
END{
	for (loci in data){
		UMI=""
		num=1
		for (umi in data[loci]){
			if (!(compare(umi,UMI,2))){
				sub("\t","_"num"\t",data[loci][umi])
				print data[loci][umi]
				UMI=umi
			}
			else num++
		}
	}
}
' $tempfile | samtools sort -O $outtype > $outfile 

rm $tempfile

