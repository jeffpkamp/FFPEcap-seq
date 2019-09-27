#!/bin/bash

ccheck=""
name="$1"
genome="$2"
UMI="$3"

#gets complexity (positions / million), unique molecules (UMI + position), and total aligned reads.

eval $(awk -v name="$name" -v OFMT=%0.1f -v CONVFMT=%0.1f -F "\t" '
	BEGIN{comp=0}
	$1!~"@"{
		if($2!=4){
			split($1,b,"_")
			RNAs[b[2]" "$3" "$4]++
			n++
			if (k<1000001){
				a[$3,$4]
				k++
			}
		}
		t++
	}
	END{
		for (x in RNAs)
			print x,RNAs[x] > name"/"name"_SeqperMolecule.tab"
		print "Unique_RNAs="length(RNAs)
		print "complexity="length(a)
		print "compscore="length(a)*100/k
		print "nreads="n
	}' ./"$name"/"$name"\_sorted.sam)




#get alignment percentage from bowtie output

eval $(awk '
	$0~"# reads that failed to align"{
		split($0,a,"(")
		sub("%)","",a[2])
		print "gperaligned="100-a[2]
	}' ./"$name"/alignout )


#reads refseq alignment to get the number of rRNA reads, total aligned refseq reads, Unique Molecules, alignment precentage

eval $(awk -v OFMT=%0.1f -v name="$name" -v CONVFMT=%0.1f -v FS="\t" '
	FNR==NR{rRNA_labels[$0];next}
	$1!~"@"{
	if ($3 != "*" && $3 != "NR_003286" && $3 != "NR_003287"){
		split($1,a,"_")
		b[a[2],$3,$4]
		gene[$3][a[2]]
		total++
		}
	if ($3 in rRNA_labels) 
		r++ 
	if($2!=4) 
		n++
	t++
	}
	END{
		print "rRNA="r 
		print "rnreads="n
		print "refperaligned="n/t*100
		print "Unique_mRNAs="length(b)
		for (g in gene){
			for (umi in gene[g])
				umi_total++  
			printf "%s\t%0.2f\n",g,umi_total > name"/"name"_geneexpression.tab"
			umi_total=0
			} 
	}' ~/FFPEcapseq/"$genome"/"$genome"\_rRNA_genes.tab ./"$name"/"$name"\_sorted_refseq.sam)


refperaligned=$(awk '
	$0~"reads with at least one reported alignment:"{
		sub(".*\\(","",$0)
		sub("%\\)","",$0)
		print $1
	}
' "$name"/refout)


rRNAali=$(echo  | awk -v r="$rRNA" -v n="$nreads" 'END{print r/n*100}')
rRNAref=$(echo  | awk -v r="$rRNA" -v n="$rnreads" 'END{print r/n*100}')




if [[ $UMI == "True" ]]
then 


#gets complexity/million for collapsed reads
	eval $(awk -v OFTMT=%0.1f -v CONVFMT=%0.1f -v FS="\t" '
		BEGIN{comp=0}
		$1!~"@"{
			if($2!=4){
				split($1,b,"_")
				RNAs[b[2]" "$3" "$4]++
				n++
				if (k<1000001){
					a[$3,$4]
					k++
				}
			}
			t++
		}
		END{
			print "col_complexity="length(a)
			print "col_compscore="length(a)*100/k
		}' ./"$name"/"$name"\_sorted_collapsed.sam)

#Calculate Gene Coverage, ie where the transcript aligns on the gene body

awk  '
	BEGIN{
		OFS="\t"
		PROCINFO["sorted_in"]="@val_num_asc"
		print "Gene","Len","pos","% pos 5 end"
	}
	$1~"@SQ"{
		sub("^SN:","",$2)
		sub("^LN:","",$3)
		a[$2]=$3
	}
	$1!~"^@"{
		if ($3 in a) 
			b[$3][$1]=$4
		}
	END{
		for (gene in b) 
			for (read in b[gene]) 
				printf "%s\t%d\t%d\t%d\n",gene,a[gene],b[gene][read],b[gene][read]/a[gene]*100
	}' ./"$name"/"$name"\_sorted_collapsed_refseq.sam > ./"$name"/"$name"\_collapsed_coverage.tab	


UMI_per_million=$(awk -F"\t" '$1!~"^@" && $2 != 4 {split($1,a,"_");print a[2]"_"$3"_"$4}' "$name"/"$name"\_sorted.sam | python ~/FFPEcapseq/bin/UMI_capture_per_million.py)

awk  '
	BEGIN{
		OFS="\t"
	}
	FNR==NR{
		rRNA_label[$0]
		next
	}
	NR>1 && (!($1 in rRNA_label))  {
		a[$4]++
		c++
		if ($4<11) 
			b++
	}
	END{
		for (i in a) 
			print i,a[i]/c*100
	}' ~/FFPEcapseq/"$genome"/"$genome"\_rRNA_genes.tab ./"$name"/"$name"\_collapsed_coverage.tab | sort -nk 1 > ./"$name"/"$name"\_collapsed_coveragesummary.tab

under10=$(awk '$1<10{b+=$2}END{print b}' ./"$name"/"$name"\_collapsed_coveragesummary.tab)
under100=$(awk '$1<100{b+=$2}END{print b}' ./"$name"/"$name"\_collapsed_coveragesummary.tab)


else

awk  '
	BEGIN{
		OFS="\t"
		PROCINFO["sorted_in"]="@val_num_asc"
		print "Gene","Len","pos","% pos 5 end"
	}
	$1~"@SQ"{
		sub("^SN:","",$2)
		sub("^LN:","",$3)
		a[$2]=$3
	}
	$1!~"^@"{
		if ($3 in a) 
			b[$3][$1]=$4
		}
	END{
		for (gene in b) 
			for (read in b[gene]) 
				printf "%s\t%d\t%d\t%d\n",gene,a[gene],b[gene][read],b[gene][read]/a[gene]*100
	}' ./"$name"/"$name"\_sorted_refseq.sam > ./"$name"/"$name"\_coverage.tab	



awk  '
	BEGIN{
		OFS="\t"
	}
	FNR==NR{
		rRNA_label[$0]
		next
	}
	NR>1 && (!($1 in rRNA_label))  {
		a[$4]++
		c++
		if ($4<11) 
			b++
	}
	END{
		for (i in a) 
			print i,a[i]/c*100
	}' ~/FFPEcapseq/"$genome"/"$genome"\_rRNA_genes.tab ./"$name"/"$name"\_coverage.tab | sort -nk 1 > ./"$name"/"$name"\_coveragesummary.tab

under10=$(awk '$1<10{b+=$2}END{print b}' ./"$name"/"$name"\_coveragesummary.tab)
under100=$(awk '$1<100{b+=$2}END{print b}' ./"$name"/"$name"\_coveragesummary.tab)

	UMI_per_million="NA"
	col_compscore="NA"
	col_complexity="NA"
fi



cat "$name"/alignout > "$name"/"$name"\_summary.txt
cat "$name"/refout >> "$name"/"$name"\_summary.txt
cat << output >> $name/$name\_summary.txt
Aligned Genomic Reads = $nreads
Unique Reads per M$ccheck = $complexity
Aligned Ref Reads = $rnreads
rRNA Reads = $rRNA
Genomic aligned % = $gperaligned
Ref Aligne % = $refperaligned
rRNA % Ref Reads  = $rRNAref
Unique gRNA Molecules = $Unique_RNAs
Unique mRNA Molecules = $Unique_mRNAs
Complexity Score (%) = $compscore
UMI Location Score(%) = $col_compscore
%10 Coverage Score$startposcheck = $under10
UMI per million = $UMI_per_million
output

#[[ ! -d ./"$name"/fastqcreport ]]  && mkdir ./"$name"/fastqcreport
#~/bin/FastQC/fastqc "$name".fastq 
#mv "$name"*.html "$name"
#rm *.zip
