#!/bin/bash
n=0

### check python
if [[ $(command -v python) == "" ]]
        then n=$((n + 1))
	echo -e "#################\n###  WARNING  ###\n#################\nPython2.7 required for program, please install and/or add to path"
        if [[ $(uname) == Darwin ]]
		then cat << 'help'

You can install this with homebrew
get homebrew with:
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
then run 
brew install python2.7

help
	elif [[ $(uname) == Linux ]]
		then cat << 'help'

If you have apt, you can install this with:
sudo apt-get install python2.7

help
	fi
fi
	
python -c 'import psutil' 2>/dev/null  || {
	echo -e "Python Module psutil required\nInstall psutil with:\necho \"pip install psutil\""
	n=$((n+1))
}
	


### check awk
if [[ $(command -v awk) == "" ]]
        then n=$((n+1))
	echo -e "#################\n###  WARNING  ###\n#################\nawk is required for program, please install and/or add to path"
        if [[ $(uname) == Darwin ]]
                then cat << 'help'

You can install this with homebrew
get homebrew with
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
then run
brew install gawk

help
        elif [[ $(uname) == Linux ]]
                then cat << help

If you have apt, you can install this with:
sudo apt-get install gawk

help
        fi
fi


### check samtools
if [[ $(command -v samtools) == "" ]]
        then 
	n=$((n+1))
	echo -e "#################\n###  WARNING  ###\n#################\nsamtools required for program, please install and/or add to path!!!"

cat << 'help'
######################
####    Warning   ####
######################

samtools version 1.5 or greater is required for program.  Please install and/or add to path

Follow these instructions to get the latest samtools:

cd ~
curl -O https://astuteinternet.dl.sourceforge.net/project/samtools/samtools/1.9/samtools-1.9.tar.bz2
tar -xjf samtools-1.9.tar.bz2
cd samtools-1.9/
./configure
make
sudo make install 
cd ~
rm samtools-1.9

now restart terminal

help


elif [[ $(samtools 2>&1| awk '$1~"Version:"{if ($2 !~ "1.[5-9]" ) {print 1}}') ]]
	then 
	n=$((n+1))
	samtools 2>&1 | awk '$1~"Version:"{if ($2 !~ "1.[5-9]" ) {print "\n\n#################\n###  WARNING  ###\n#################\n\nCurrent Samtools Version is "$2; print "samtools version needs to be updated to at least 1.9\n" }}'
	cat << 'help'
Follow these instructions to get the latest samtools:

cd ~
curl -O https://astuteinternet.dl.sourceforge.net/project/samtools/samtools/1.9/samtools-1.9.tar.bz2
tar -xjf samtools-1.9.tar.bz2
cd samtools-1.9/
./configure
make
sudo make install
cd ~
rm samtools-1.9

now restart terminal

help
fi



### check bowtie
if [[ $(command -v bowtie) == "" ]]
	then n=$((n+1))
	if [[ $(uname) == Darwin ]]
	then
cat << 'help'
######################
####    Warning   ####
######################

bowtie is required for program.  Please install and/or add to path.

get homebrew with:

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"


install bowtie with:

brew tap brewsci/bio 
brew install bowtie

homebrew (macOS)



or install from the source code here: 
https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2

help

	elif [[ $(uname) == Linux ]]
	then
cat << 'help'
######################
####    Warning   ####
######################

bowtie is required for program.  Please install and/or add to path.

Bowtie can be installed using apt-get (debian) or installed from the source code here: 
https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2

help
	fi

fi

if [[ $n -eq 0 ]]
	then echo 
	echo "All required programs are present!"  
	echo
	exit 0
else 
	exit 1
fi

