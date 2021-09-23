# fdfBLAST_version2021


================

The current working version fdfBLAST_version2021 is an update of fdfBLAST_1.20120704.pl. The new version is executed by a shell script, instead of the original interface input.

If you have any questions, please contact yuann@big.ac.cn.

================





## Install:


PERL Modules Math::BigFloat BioPerl GD

HMMER3 and PFAM-A

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz

## Running fdfBLAST 

### 1. Run BLAST+

makeblastdb for each genome and blastp alignment for pairwise genome.

main parameters:

num_threads(running cpu), evalue(set according to species,ex:1e-70), num_alignments(ex:5)


for i in `ls *.fas`;do 

	./ncbi-blast-2.12.0+/bin/makeblastdb  -in $i -dbtype prot  -parse_seqids  -out $i

done


for i in `cat genomes_lst`;do

	for j in `cat genomes_lst`;do 
	
      	./ncbi-blast-2.12.0+/bin/blastp -db $i.fas -num_threads 10  -outfmt 0 -query $j.fas -out  g2gc/"$j"_"$i".bpo
      
	done
	
done


### 2. Generate reciprocal hits

perl step3.pl $genomo_dir $run_dir $pValue

#ex：

#perl step3_hits.pl /xtdisk/liufan_group/yuanna/fdf/fdfBLAST-master/4genomes ./ 1e-10

### 3. Select best hits

Run this step for every pairwise genome.

#ex:

#perl split_3part_by_2pair_alignBest.pl $gene_hit_dir  GenomeA genomeB $output_dir

### 4. PFAM domains and draw domains

cd runDir

#### 1) First run step5.pl and  generate file: composite_list.csv and split_list.csv.

touch output.parse

perl step5.pl output.parse  ./results  ./
<br>
#### 2) Sort and uniq file: composite_list.csv and split_list.csv

cat composite_list_sorted.csv split_list_sorted_uniq.csv |sort|uniq  > all_domains.csv
<br>
#### 3) Extract seq from genome file

perl  -pe 's[,][]g' all_domains.csv > ex_lst

##### #first makeblastdb

../../ncbi-blast-2.12.0+/bin/makeblastdb -in GenomeA.fas -dbtype prot  -parse_seqids  -out GenomeA

../../ncbi-blast-2.12.0+/bin/makeblastdb -in GenomeB.fas -dbtype prot  -parse_seqids  -out GenomeB

../../ncbi-blast-2.12.0+/bin/makeblastdb -in GenomeC.fas -dbtype prot  -parse_seqids  -out GenomeC

##### #extract sequences from genome using ex_lst

../../../ncbi-blast-2.12.0+/bin/blastdbcmd -db  ../../3genomes/GenomeA  -entry_batch ex_lst   -out ex_GenomeA.fa 

../../../ncbi-blast-2.12.0+/bin/blastdbcmd -db  ../../3genomes/GenomeB  -entry_batch ex_lst   -out ex_GenomeB.fa 

../../../ncbi-blast-2.12.0+/bin/blastdbcmd -db  ../../3genomes/GenomeC  -entry_batch ex_lst   -out ex_GenomeC.fa
<br>
####  4) Run Pfam and parse

cat ex_*.fa >all_domain_seqs.fas

../../../pfam/hmmer3-3/bin/hmmscan  --cpu 6 --noali  ../../../pfam/Pfam-A.hmm   all_domain_seqs.fas >output.out 

perl hmm3scan_parser.pl -i output.out -o output.parse
<br>
#### 5) Run step5.pl again using domains file（output.parse） 

perl step5.pl output.parse  ./results ./


