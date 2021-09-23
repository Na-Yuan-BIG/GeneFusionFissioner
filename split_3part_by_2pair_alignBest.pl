#!/usr/bin/perl -w
#Filename:
#Author:	Na Yuan
#EMail:		yuann@big.ac.cn
#Modified:	
#Description: 
my $version=1.00;
#use strict;
my $HIT_LIMIT=250;
my $GENE_HITS_DIR=$ARGV[0];
my $OUT_DIR=$ARGV[3];
#example run:
#for GenomeC 
#perl split_3part_by_2pair_align.pl gene_hit_dir  GenomeC GenomeA output_dir
#perl split_3part_by_2pair_align.pl gene_hit_dir  GenomeC GenomeB output_dir
#perl split_3part_by_2pair_align.pl gene_hit_dir  GenomeC genomeD output_dir

#perl split_3part_by_2pair_align.pl gene_hit_dir  GenomeC genomeD output_dir
#opening initial/genomeD_GenomeC.csv 

my $tfile=$GENE_HITS_DIR."/initial/".$ARGV[2]."_".$ARGV[1].".csv";
open IN2, "$tfile";
my %rev;
while(<IN2>){
	next if(/^0,/);
	my @r=split(/,/,$_);
	$rev{$r[1]}=$_;
}
close IN2;


#opening initial/GenomeA_genomeD.csv
my $qfile=$GENE_HITS_DIR."/initial/".$ARGV[1]."_".$ARGV[2].".csv";
my $non_file=$OUT_DIR."/no_hits_".$ARGV[2].".csv";
my $dup_file=$OUT_DIR."/duplications_".$ARGV[1]."_fixed.csv";
my $dif_file=$OUT_DIR."/differentials_".$ARGV[1]."_fixed.csv";
open IN1, "$qfile";
open NON,">>$non_file";
open DUP,">>$dup_file";
open DIF,">>$dif_file";
my %dupid;
while(<IN1>){
	my @t=split(/,/,$_);
	next if($t[0]>$HIT_LIMIT);
	if($t[0] ==0 ){
		print NON $t[1],",\n";
	}elsif($t[0] ==1){	
		if(defined $rev{$t[3]}){
			my @s=split(/,/,$rev{$t[3]});
			if($t[1] eq $s[3]){
				print DUP "0;",$t[1],";",$t[2],";",$t[3],";",$t[4],";",$t[5],";",$t[6],";\n";
			}
		}
	}else{
		my $num=0;
		for ( my $k=1;$k<=$t[0];$k++){
			my $lx=3+($k-1)*6;
			if(defined $rev{$t[$lx]}){
				my @s2=split(/,/,$rev{$t[$lx]});
				if($t[1] eq $s2[3]){
					print DIF $num,";",$t[1],";",$t[2],";",$t[$lx],";",$t[$lx+1],";",$t[$lx+2],";",$t[$lx+3],";",$t[$lx+4],";",$t[$lx+5],";\n";
					$dupid{$t[1]}=1;
					if(defined $dupid{$t[1]}){
						$num++;
					}
				}
			}	
		}
	}
}
close IN1;
close NON;
close DUP;
close DIF;
