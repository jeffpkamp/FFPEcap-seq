curl -O  http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/FFPEcapseq.tar
tar -xf FFPEcapseq.tar
[[ "$(pwd)" != "$HOME" ]] && mv FFPEcapseq "$HOME"/FFPEcapseq
rm FFPEcapseq.tar

echo "Set up FFPEcapseq in your home directory."

sleep 1

echo "Checking for required files and programs"

bash "$HOME"/FFPEcapseq/bin/requirment_check.sh

get_genomes="Y"

while [[ $get_genomes == "Y" || $get_genomes == "y" ]]
        do PS3='Select Genome: '
        cd ~/FFPEcapseq || (echo 'No ~/FFPEcaseq folder!'; exit 1)
        list=($(curl -s http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/Genomes))
        list=(quit ${list[@]})
	echo -e "\nDownload Indexes to run FFPEcapseq\nGenome Options:"
        select x in "${list[@]}"
                do 
		if [[ "$x" == "quit" ]]
                	then break
                elif [[ -z $x ]]
                        then echo "Invalid Option"
                else
                        if [[ ! -d $x ]]
                                then
				curl -O http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes/"$x".tar.gz
                                tar -xzf "$x".tar.gz
                                rm "$x".tar.gz
                                break
                        else
                                echo -e "\n$x is already present.\nIf you want to replace, please delete ~/FFPEcapseq/$x before trying to install.\n"
                                break
                        fi
                fi
        done
        read -p "Do you want to get any other bowtie indexes? (y/n):" get_genomes
        if [[ $get_genomes != "y" && $get_genomes != "Y" ]]
                then cd "$startdir"
        fi
done

echo

echo "To recheck requirements:"
echo "~/FFPEcapseq/bin/requirment_check.sh again to make sure all required programs are installed"

echo

echo "To Run Program:" 
echo "~/FFPEcapseq/bin/FFPEcapseq.sh to start program"

echo

echo "For more information read ~/FFPEcapseq/README.txt"
