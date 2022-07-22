# napSeeker

napSeeker: A computational software for identifying non-capped RNA from NAP-seq data.

Overview:
---------

napSeeker is a software to identify non-capped RNA by judging the position of special adapters and calculating the coverage of the candidate non-capped RNA from NAP-seq data. 

Usage:
---------

Usage:  ncapSeeker [options] --fa <genome seq> --fai <fai file> --input <BAM alignments><BR>
ncapSeeker: for discovering novel napRNAs from NAP-seq data<BR>
[options]<BR>
-v/--verbose                   : verbose information<BR>
-V/--version                   : napSeeker version<BR>
-h/--help                      : help informations<BR>
-n/--norm                      : normalized reads to the locus number<BR>
-N/--nano                      : nanopore sequencing<BR>
-f/--fold                      : rnafold the sequence with length<500nt<BR>
-k/--keep-dup                  : keep duplication, deault is false<BR>
-P/--pair                      : input is paired-end format<BR>
-s/--strand                    : strand-specific<BR>
-i/--input <string>            : input file<BAM format><BR>
-o/--outfile <string>          : output file<BR>
-t/--min-tag <double>          : minimum tag number for each start or end pos, default>=5.0 reads<BR>
-c/--cov <double>              : minimum coverage in each contig, default>=2 reads<BR>
-r/--rpm <double>              : minimum rpm value for each psi, default>=0.05<BR>
-l/--min-len <int>             : minimum length of reads, default=15<BR>
-L/--max-len <int>             : maximum length of reads, default=1000<BR>
-m <int>                       : minimum cluster length default=20<BR>
-M <int>                       : maximum cluster length, default=10000<BR>
-F <double>                    : minimum fold-change between start/end and mean reads[default=2]<BR>
-e <double>                    : maximum fold-change between start and end reads[default=5]<BR>


Installation:<BR>
---------

Typical install time: within 5 min.
Download napSeeker-1.0.tar.gz from https://github.com/junhong-huang/napSeeker/releases ; unpack it, and make:<BR>
tar -xzvf napSeeker-1.0.tar.gz<BR>
cd napSeeker-1.0<BR>
make<BR>

System requirements:
---------

Operating system: napSeeker is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.<BR>
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers.<BR>
Libraries and other installation requirements: napSeeker includes other software libraries: BamTools, cdflib, alglibsrc, RNAfoldLib, and RNAshapesLib library package. All will automatically compile during napSeeker installation process.<BR>
By default, napSeeker does not require any additional libraries to be installed by you.<BR>

Prerequisites:<BR>
---------

Dependencies: The input of napSeeker is BAM file. So you need the read mapper STAR or other mappers<BR>
You can get the most fresh versions:<BR>
(1)  STAR: https://github.com/alexdobin/STAR<BR>
(2)  Samtools: http://www.htslib.org/<BR><BR>You need to have the reference genome, fai file, and  STAR indexes.<BR>You can constructed these datasets by yourself using following steps:<BR>
As an example, let's assume you use human genome (version hg38).<BR>
(1)  Genome:<BR>
mkdir genome<BR>
wget -c 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'<BR>
gzip -d hg38.fa.gz<BR>
cd ..<BR>
(2) Build the genome index and align reads to genome:<BR>
STAR --runMode genomeGenerate --genomeDir ./starIndex --genomeFastaFiles ./genome/hg38.fa --sjdbGTFfile gencode.v30.annotation.gtf --sjdbOverhang 100<BR>
(3)build the fai index:<BR>
samtools faidx hg38.fa<BR>
(4)  Align reads to genome using STAR:<BR>
--genomeLoad NoSharedMemory --limitBAMsortRAM 60000000000 --alignEndsType EndToEnd --outFilterType BySJout --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMatchNmin 20 --outFilterMatchNminOverLread 0.8 --seedSearchStartLmax 15 --seedSearchStartLmaxOverLread 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 20 --alignSJDBoverhangMin 10 --outSAMtype BAM Unsorted --outSAMmode Full --outSAMattributes All --outSAMunmapped None --outSAMorder Paired --outSAMprimaryFlag AllBestScore --outSAMreadID Standard --outReadsUnmapped Fastx --limitOutSJcollapsed 5000000 --alignEndsProtrude 150 ConcordantPair --readFilesCommand zcat<BR>

run napSeeker:
---------

napSeeker --fa hg38.fa --fai hg38.fa.fai --input NAP-seq.bam \> napSeeker_candidate_napRNAs.txt<BR>

Output:
---------

#chrom	start	end	name	score	strand	contigLen	clusterLen	readNum	startReadNum	endReadNum	startPvalue	endPvalue	startFold	endFold	upFold	downFold	up20ntFold	down20ntFold	clusterSeq	cdScore	structure	shape	mfe<BR>
chr8	123434450	123434610	ncapSeeker_chr8_123434450_123434610	48763	+	160	4070	48762.50758	22038.50758	26724	-160427	-198624	117.5569	162.15689	272.25286	125.3105	37.45919	16.31894	ACCCCATCTCTACTAAAAATACAAAAATTAGCTGGGTGTGGTGGCGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGTGGGAAAATCACTCGAACCCCGGAGGTGGAGGTTGCAGTGAGCTGGGATTGTGCCACTGTGCTCCAGCCTGGGTGACAGAG	0	..............................((((((((((((((((((......((((((((..(((((....))))).........(((((..(((((((....))).))))..)))))))))))))))))))))))))).)))))(((.....)))..	[[][]][]	-60.4<BR>

Note: # is comment line<BR>

Acknowledgements:
---------

Thanks a lot to everyone who contributed to the public code (e.g. BamTools, Samtools) used by napSeeker.<BR>

Contact :
---------

*****************************************************************************************<BR>
 \*	napSeeker - A computational software for identifying non-capped RNA from NAP-seq data.<BR>
 \*<BR>
 \*	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn><BR>
 \* <BR>
 \*	RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
 \*	<BR>
 \*  Create date: 02/07/2022<BR>
 \*  <BR>
 \*  last modified time: 02/07/2022<BR>
 ****************************************************************************************<BR>