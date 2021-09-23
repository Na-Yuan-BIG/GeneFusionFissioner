#!/usr/bin/perl
#Filename:
#Author:	Na Yuan
#EMail:		yuann@big.ac.cn
#Modified:	
#Description: 

#使用说明：
#cd /xtdisk/liufan_group/yuanna/fdf/fdfBLAST-master/run/1626318028
#perl step3.pl /xtdisk/liufan_group/yuanna/fdf/fdfBLAST-master/4genomes ./ 1e-10

use Bio::SearchIO;     # Bioperl for input/output/persing of BLAST etc
use Cwd;               # Gets pathname of current working directory
use File::Basename;    # Remove path information and extract 8.3 filenames
use GD;                # Creates PNG images

# use GD::SVG;               # Creates SVG images - uncomment if you want SVG output
use Math::BigFloat;        # Arbitrary size floating point math package (handles e-values)
use feature qw(switch);    # This is the given replacement of deprecated switch statement
use Time::Local;           # For time elapsed when running different stages


use Getopt::Long;


my $version=1.00;

my $GENOME_DIR=$ARGV[0];
my $G2GC_DIR="$GENOME_DIR/g2gc";
my $RUN_DIR=$ARGV[1];
my $lower_limit = "0";
$upper_limit = $ARGV[2];

my $GENE_HITS_DIR     = "$RUN_DIR/gene_hits";
my $GENE_HITS_INITIAL = "$GENE_HITS_DIR/initial";
my $GENE_HITS_LOOKUP  = "$GENE_HITS_DIR/lookup";



@GENOME_LIST = get_genomes();
my @g2gc_file_names = glob "$G2GC_DIR/*.bpo";
run_gene_hits_helper();
run_gene_hits( \@g2gc_file_names, $upper_limit, $lower_limit );
generate_lookup_tables(@GENOME_LIST);



sub run_gene_hits {

    my @g2gc_files  = @{ $_[0] };
    my $upper_limit = $_[1];
    my $lower_limit = $_[2];

    #print "XXX\n@g2gc_files\n$#g2gc_files\n$upper_limit\n$lower_limit\nXXX\n";

    for ( my $i = 0 ; $i <= $#g2gc_files ; $i++ ) {
        my $current = $i + 1;
        print "\n$current of " . ( $#g2gc_files + 1 );

        # Assign file name in iterative loop
        my $in_file = "$g2gc_files[$i]";

        my ( $in_filename, $dir, $ext ) = fileparse( $in_file, '\..*' );

        print "    $in_filename";

        # X = Blast sequence number
        my $x_pos = 0;

        undef(@comparison_length);
        #my @comparison_length = $EMPTY;

        my $search = new Bio::SearchIO(
            '-format' => 'blast',
            '-file'   => $in_file
        );
        while ( my $result = $search->next_result ) {
            $y_pos                   = 0;
            $query_accession         = $result->query_accession;
            $query_length            = $result->query_length;
            $query_genome[$x_pos][0] = "$query_accession";
            $query_genome[$x_pos][1] = "$query_length";
            while ( my $hit = $result->next_hit ) {
                $hit_accession    = $hit->accession;
                $hit_length       = $hit->length;
                $hit_significance = $hit->significance;

                while ( my $hsp = $hit->next_hsp ) {
                    $percent_identity = $hsp->percent_identity;

                    # Range of query that hits the subject gene
                    @query_range = $hsp->range('query');
                    #####@hit_range   = $hsp->range('hit');

                    # We don't really need 13 decimal places, round to none..
                    $percent_identity = sprintf( "%.0f", $percent_identity );
                }

                $comparison_genome[$x_pos][$y_pos] =
                  "$hit_accession,$hit_length,$hit_significance,$percent_identity,$query_range[0],$query_range[1],";
                $y_pos++;
            }

            push( @comparison_length, $y_pos );
            $x_pos++;
        }

        open my $out_fh, '>', "$GENE_HITS_INITIAL/$in_filename.csv";
        ##### Surely we don't need this twice? 
        # I don't know Guy, what were you trying to do?
        $query_genome_length = $x_pos;

        #Append hit numbers to query_genome array
        for ( my $x = 0 ; $x < $query_genome_length ; $x++ ) {


            my $comparison_genome_length = $comparison_length[$x];
            my $count                    = 0;

            for ( my $d = 0 ; $d < $comparison_genome_length ; $d++ ) {
                $evalue = "$comparison_genome[$x][$d]";
                $evalue =~ m/(.*?)\,(\d*)\,(.*?)\,(.*?)\,/;
                $evalue = $3;
                my $value = evaluate( $evalue );
                $evalue = $value;
                if (   ( $evalue <= $upper_limit )
                    && ( $evalue >= $lower_limit ) ) {
                    $count++;
                }
            }
            $query_genome[$x][2] = $count;
        }

        for ( my $x = 0 ; $x < $query_genome_length ; $x++ ) {

            print $out_fh "$query_genome[$x][2],$query_genome[$x][0],$query_genome[$x][1],";

            my $comparison_genome_length2 = $query_genome[$x][2];

            for ( $y = 0 ; $y < $comparison_genome_length2 ; $y++ ) {

                $evalue = "$comparison_genome[$x][$y]";
                $evalue =~ m/(.*?)\,(\d*)\,(.*?)\,(.*?)\,/;
                $evalue = $3;
                my $value = evaluate( $evalue );
                $evalue = $value;

                ## 0 is the lower limit as the lower the E-value,
                ## or the closer it is to zero, the more "significant" the match is
                ## therefore
                if (   ( $evalue <= $upper_limit )
                    && ( $evalue >= $lower_limit ) )
                {
                    print $out_fh "$comparison_genome[$x][$y]";
                }
            }
            print $out_fh "\n";
        }
        close($out_fh);

        undef(@comparison_genome);
        undef(@query_genome);
    }
    print "\nFinished\n";

}

# A method to create tables of the recorded hits for each genome group
sub generate_lookup_tables {

    my @genomes = @_;

    print "Generating Lookup Tables";

    # Check to see if the initial hit lists exist
    if ( -e $GENE_HITS_INITIAL && -d $GENE_HITS_INITIAL ) {

        # Edited gene hit lists get their own directory
        mkdir( $GENE_HITS_LOOKUP, 0755 );

        # Internal counter
        $run = 0;

        # Two for loops to iterate through each gene hits initial file
        # and add the values to a 2d array for output to file
        for ( $i = 0 ; $i <= $#genomes ; $i++ ) {
            for ( $j = 0 ; $j <= $#genomes ; $j++ ) {

                ( $file_i, $dir, $ext ) = fileparse( $genomes[$i], '\..*' );
                ( $file_j, $dir, $ext ) = fileparse( $genomes[$j], '\..*' );
                open my $in_fh, '<', "$GENE_HITS_INITIAL/$file_i\_$file_j\.csv";

                # We're using while with an iterator instead of foreach to
                # read in the file line by line adding each column to the 2d array
                $x = 0;
                while (<$in_fh>) {
                    my ($line) = $_;
                    chomp($line);
                    @temp = split( /,/, $line );
                    $AoA[$x][$j] = $temp[0];
                    $x++;
                }

                # Then we add up the values in each column for each file
                # using a for loop
                $total_in_array_column = 0;
                for my $i ( 0 .. $#AoA ) {
                    $total_in_array_column += $AoA[$i][$j];
                }

                # Then those values are added plainly to an array,
                # the sum of all values is calculated and if they are
                # equal to 0 then there are no hits, else we continue...
                push( @check_for_zero, $total_in_array_column );
                $total_of_array_columns += $_ for @check_for_zero;
                if ( $total_of_array_columns == 0 ) {
                    print "\nThere are no hits at all, please try some different e-values...\n";
                    last;
                }
                ###
            }    # End second (internal) for loop
            ##close($in_fh);
            open my $out_fh, '>', "$GENE_HITS_LOOKUP/$file_i\.csv";

            # Here we use two for loops to iterate through the 2D array
            # for $i from 0 to length of @AoA
            for my $i ( 0 .. $#AoA ) {

                # A reference to the row at position $i
                $aref = $AoA[$i];

                # Length of row
                $n          = @$aref - 1;
                $line_total = 0;

                # Second for, $j from 0 to length of row
                for my $j ( 0 .. $n ) {

                    # Sum the individual numbers of each row
                    $line_total = $line_total + $AoA[$i][$j];

                    # Print out each element of the row
                    print $out_fh "$AoA[$i][$j],";
                }
                print $out_fh "\n";
            }

            # Increment run
            $run = $run + 1;

            # Reset array, so length is not kept at largest
            # @AoA = ();
            undef(@AoA);
            print ".";
        }
        ##close($out_fh);
    }
    print "\n";
}
sub evaluate {
    my $value = shift;

    # Perl method bstr will only change '1e' not 'e' to decimal, therefore prefix value with '1'
    if ( $value =~ m/^e/ ) {
        $value = "1" . $value;
        $value = Math::BigFloat->bstr($value);
    }
    if ( $value =~ m/e/ ) {
        $value = Math::BigFloat->bstr($value);
    }
    return $value;
}

sub run_gene_hits_helper {

    # Create the gene_hits output folder, if unable then on error, quit.
    # This could probably be handled better - overwrite? sub directory? etc
    mkdir( $GENE_HITS_DIR, 0755 )
      || die "gene_hits has already been make, please delete gene_hits dir and then run.\n";

    mkdir( $GENE_HITS_INITIAL, 0755 );
}

sub get_genomes {
    my @file_names = glob "$GENOME_DIR/*.fas";
    foreach my $file (@file_names) {
        $file = fileparse($file);
    }
    return @file_names;
}