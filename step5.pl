#!/usr/bin/perl
#Filename:
#Author:	Na Yuan
#EMail:		yuann@big.ac.cn
#Modified:	
#Description: 


# Import Modules
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
my $lower_ratio = "0.1";
my $higher_ratio = "1.0";

my  $domain_in = $ARGV[0] ;
my  $dir = $ARGV[1] ;
my $GENE_HITS_DIFFERENTIALS="$ARGV[1]/differentials";
my $RUN_DIR=$ARGV[2] ;

&menu_five;

sub menu_five {
    my $menu_choice = "Identify Fusions";
    read_domain_file();
    fusion_scan();
}

sub fusion_scan {

    print "\nPerforming Fusion Scans\n";

    @fusion_filenames = glob "$GENE_HITS_DIFFERENTIALS/*fixed.csv";
    $fusion_number    = @fusion_filenames;

    print "Number of files = $fusion_number\n";

    for ( my $i = 0 ; $i < $fusion_number ; $i++ ) {

        undef(@fusion);
        open my $fusion_fh, '<', "$fusion_filenames[$i]";
        $input_file = fileparse( $fusion_filenames[$i] );
        print "Opening $fusion_filenames[$i]\n";

        while (<$fusion_fh>) {
            my ($line) = $_;
            chomp($line);
            push( @fusion, "$line" );
        }

        $fusion_array_length = @fusion;

        for ( my $j = 0 ; $j < $fusion_array_length ; $j++ ) {

            @array_line = split( /;/, $fusion[$j] );
            $line_number = $array_line[0];

            if ( $line_number == 0 ) {

                $query_accession = $array_line[1];
                $query_length    = $array_line[2];

                $subject_accession       = $array_line[3];
                $subject_length          = $array_line[4];
                $subject_evalue          = $array_line[5];
                $subject_hit_range_start = $array_line[7];
                $subject_hit_range_end   = $array_line[8];
                $match_length            = $subject_hit_range_end - $subject_hit_range_start;
                push( @subject_hits,
                    "$subject_accession,$subject_hit_range_start,$subject_hit_range_end,$subject_length" )
                  if $match_length > 50 && $subject_length > 50;
                #####
                print "\nQ $query_accession -> $subject_accession\n"
                  if $DEBUG >= 1;
                #####

                for ( my $k = $j + 1 ; $k < $fusion_array_length ; $k++ ) {

                    @next_array_line = split( /;/, $fusion[$k] );
                    $next_line_number = $next_array_line[0];

                    if ( $next_line_number != 0 ) {

                        $next_query_accession = $next_array_line[1];

                        $next_subject_accession       = $next_array_line[3];
                        $next_subject_length          = $next_array_line[4];
                        $next_subject_evalue          = $next_array_line[5];
                        $next_subject_hit_range_start = $next_array_line[7];
                        $next_subject_hit_range_end   = $next_array_line[8];

                        # > 50 should be user selectable limit for now it is hard-coded
                        # to exclude all sequences with <=50 bases.
                        # Another limit is Matched Length (end - start) this should also be
                        # more than 50.
                        $next_match_length = $next_subject_hit_range_end - $next_subject_hit_range_start;

                        if (   $next_line_number != 0
                            && $query_accession eq $next_query_accession
                            && $next_subject_length > 50
                            && $next_match_length > 50 )
                        {
                            print
"\n$next_line_number ne 0 && $query_accession eq $next_query_accession && $next_subject_length > 50 && $next_match_length > 50\n"
                              if $DEBUG >= 1;
                            push( @subject_hits,
"$next_subject_accession,$next_subject_hit_range_start,$next_subject_hit_range_end,$next_subject_length,$next_subject_evalue"
                            );
                        }
                    }
                    else {
                        last;
                    }
                }
                &rank_sort;
                undef(@subject_hits);
            }
        }
        close($fusion_fh);
    }
}

sub rank_sort {
    $subject_hit_length = @subject_hits;

    # Get fusion ORF length and half
    $fusionORF = $query_length;
    $half_fusion_length = sprintf( "%.0f", $query_length / 2 );
    print "Fusion Length - $fusionORF / 2 = $half_fusion_length\n"
      if $DEBUG >= 1;

    for ( $q = 0 ; $q < $subject_hit_length ; $q++ ) {
        @unfused = split( /,/, $subject_hits[$q] );

        #print "Q = $q\n";
        if ( $unfused[1] <= $half_fusion_length ) {
            $pos1 = "L";
        }
        elsif ( $unfused[1] >= $half_fusion_length ) {
            $pos1 = "R";
        }
        if ( $unfused[2] <= $half_fusion_length ) {
            $pos2 = "L";
        }
        elsif ( $unfused[2] >= $half_fusion_length ) {
            $pos2 = "R";
        }

        #print "$pos1 - $unfused[1] , $unfused[2] - $pos2\n";
        if ( "$pos1$pos2" eq "LL" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is Leftmost\n";
            push( @leftmost, "$subject_hits[$q]" );
        }
        elsif ( "$pos1$pos2" eq "RR" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is Rightmost\n";
            push( @rightmost, "$subject_hits[$q]" );
        }
        elsif ( "$pos1$pos2" eq "LR" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is more left\n";
            push( @middles, "$subject_hits[$q]" );
        }
        elsif ( "$pos1$pos2" eq "RL" ) {

            #print "$pos1$pos2\t$unfused[0] ORF is more right\n";
            push( @middles, "$subject_hits[$q]" );
        }
    }

    $left_num   = @leftmost;
    $right_num  = @rightmost;
    $middle_num = @middles;
    print "L:$left_num\tR:$right_num\tM:$middle_num\n" if $DEBUG >= 1;

    if ( $left_num == 0 && $right_num == 0 && $middle_num == 0 ) {

        # Discard
    }
    elsif ($left_num == 0 && $right_num == 0
        || $right_num == 0 && $middle_num == 0
        || $left_num == 0  && $middle_num == 0 )
    {

        # Discard
    }
    elsif ( $left_num == 0 && $right_num >= 1 && $middle_num >= 1 ) {

        &ignore_orthologues;
        print "\t\tXX $continue XX\n" if $DEBUG >= 1;
        if ( $continue eq "yes" ) {
            &middle_and_right;
        }
    }
    elsif ( $right_num == 0 && $left_num >= 1 && $middle_num >= 1 ) {

        &ignore_orthologues;
        print "\t\tXX $continue XX\n" if $DEBUG >= 1;
        if ( $continue eq "yes" ) {
            &middle_and_left;
        }

    }
    elsif ( $middle_num == 0 && $left_num >= 1 && $right_num >= 1 ) {

        print "011\n" if $DEBUG >= 1;
        &left_and_right;

    }
    elsif ( $left_num >= 1 && $right_num >= 1 && $middle_num >= 1 ) {

        #print "\nXX - All - XX\n";
        &ignore_orthologues;
        print "\t\tXX $continue XX\n" if $DEBUG >= 1;
        if ( $continue eq "yes" ) {
            print "111\n" if $DEBUG >= 1;
            &left_and_right;
            &middle_and_left;
            &middle_and_right;
        }
    }
    undef(@middles);
    undef(@leftmost);
    undef(@rightmost);
}
sub left_and_right {
    print "\nLR\n" if $DEBUG >= 1;
    for ( my $a = 0 ; $a < $left_num ; $a++ ) {
        print "A:$a\n" if $DEBUG >= 1;
        @unfused_one = split( /,/, $leftmost[$a] );
        my $one_end = $unfused_one[2];

        for ( my $b = 0 ; $b < $right_num ; $b++ ) {
            print "\tB:$b\n" if $DEBUG >= 1;
            @unfused_two = split( /,/, $rightmost[$b] );
            my $two_start = $unfused_two[1];

            print "\t$query_accession -> $unfused_one[0] + $unfused_two[0]\n"
              if $DEBUG >= 1;

            $ratio = $one_end / $two_start;
            $ratio = sprintf( "%.2f", $ratio );
            print "\tRlr: $one_end / $two_start = $ratio\n" if $DEBUG >= 1;

            if ( $ratio >= $lower_ratio && $ratio <= $higher_ratio ) {
                undef(@subject_hits);
                push( @subject_hits, "$unfused_one[0],$unfused_one[1],$unfused_one[2],$unfused_one[3]" );
                push( @subject_hits, "$unfused_two[0],$unfused_two[1],$unfused_two[2],$unfused_two[3]" );
                $subject_hit_length = @subject_hits;
                &generate_image;
            }
        }
    }
}

sub middle_and_left {
    print "\nML\n" if $DEBUG >= 1;
    for ( my $a = 0 ; $a < $left_num ; $a++ ) {
        print "A:$a\n" if $DEBUG >= 1;
        @unfused_one = split( /,/, $leftmost[$a] );

        my $one_end = $unfused_one[2];

        for ( my $b = 0 ; $b < $middle_num ; $b++ ) {
            print "\tB:$b\n" if $DEBUG >= 1;
            @unfused_two = split( /,/, $middles[$b] );

            my $two_start = $unfused_two[1];
            $two_match_length = $unfused_two[2] - $unfused_two[1];
            print "$query_accession -> $unfused_one[0] + $unfused_two[0]\n"
              if $DEBUG >= 1;

            $ratio = $one_end / $two_start;
            $ratio = sprintf( "%.2f", $ratio );

            if ( $two_start <= $one_end ) {
                print "\t$two_start <= $one_end Overlap!\n" if $DEBUG >= 1;
            }
            else {
                print "\tRml: $one_end / $two_start = $ratio\n"
                  if $DEBUG >= 1;
                if ( $ratio >= $lower_ratio && $ratio <= $higher_ratio ) {
                    undef(@subject_hits);
                    push( @subject_hits, "$unfused_one[0],$unfused_one[1],$unfused_one[2],$unfused_one[3]" );
                    push( @subject_hits, "$unfused_two[0],$unfused_two[1],$unfused_two[2],$unfused_two[3]" );
                    $subject_hit_length = @subject_hits;
                    &generate_image;
                }
            }
        }
    }
}

sub middle_and_right {
    print "\nMR\n" if $DEBUG >= 1;
    for ( my $a = 0 ; $a < $middle_num ; $a++ ) {
        print "A:$a\n" if $DEBUG >= 1;
        @unfused_one = split( /,/, $middles[$a] );

        my $one_end = $unfused_one[2];

        for ( my $b = 0 ; $b < $right_num ; $b++ ) {
            print "\tB:$a\n" if $DEBUG >= 1;
            @unfused_two = split( /,/, $rightmost[$b] );

            my $two_start = $unfused_two[1];
            $two_match_length = $unfused_two[2] - $unfused_two[1];
            print "$query_accession -> $unfused_one[0] + $unfused_two[0]\n"
              if $DEBUG >= 1;

            $ratio = $one_end / $two_start;
            $ratio = sprintf( "%.2f", $ratio );

            if ( $one_end >= $two_start ) {
                print "\t$one_end >= $two_start Overlap!\n" if $DEBUG >= 1;
            }
            else {
                print "\tRmr: $one_end / $two_start = $ratio\n"
                  if $DEBUG >= 1;
                if ( $ratio >= $lower_ratio && $ratio <= $higher_ratio ) {
                    undef(@subject_hits);
                    push( @subject_hits, "$unfused_one[0],$unfused_one[1],$unfused_one[2],$unfused_one[3]" );
                    push( @subject_hits, "$unfused_two[0],$unfused_two[1],$unfused_two[2],$unfused_two[3]" );
                    $subject_hit_length = @subject_hits;
                    &generate_image;
                }
            }
        }
    }
}


sub generate_image {

    $query_acc_size = ( length $query_accession );
    push( @names_length, $query_acc_size );

    $number_of_hits = $subject_hit_length;

    $unfused_one_length = ( length $unfused_one[0] );
    push( @names_length, $unfused_one_length );
    $unfused_two_length = ( length $unfused_two[0] );
    push( @names_length, $unfused_two_length );

    # Get longest gene accession and * 6 for gene accession length in 'pixels'
    my $highest = 0;
    for (@names_length) {
        $highest = $_ if $_ > $highest;
    }
    $name_length = $highest * 6;

    # Image specific
    $padding_left  = 100;
    $padding_right = 100;
    $padding       = $padding_left + $padding_right;

    # place to start drawing sequences from
    $left_pos = $padding_left + $name_length;

    # Image Width(length?)
    $image_length = $query_length + $padding + $name_length;
    if ( $image_length < 500 ) {
        $image_length += 250;
    }

    # number of hits * 50 + 50 for query + 20 for padding
    $image_height = ( $number_of_hits * 50 ) + 50 + 20;

    # set positions
    $font_vpos      = 2;
    $bar_vpos       = 20;
    $accession_vpos = 15;
    $start          = 0;
    $end            = 0;

    # create a new image
    $im = new GD::Image( $image_length, $image_height, 1 );    # Not paletted - allows for alpha transparency of domains

    # allocate some specific colors
    $white = $im->colorAllocate( 255, 255, 255 );
    $black = $im->colorAllocate( 0,   0,   0 );
    $red   = $im->colorAllocate( 255, 0,   0 );
    $blue  = $im->colorAllocate( 0,   21,  181 );

    # make the background non-transparent and interlaced
    $im->transparent(-1);                                      # no transparency
    $im->interlaced('true');
    $im->alphaBlending(1);

    # White background
    $im->filledRectangle( 0, 0, $image_length, $image_height, $white );

    # Put a black frame around the picture
    $im->rectangle( 0, 0, $image_length - 1, $image_height - 1, $black );

    $actual_length = "";
    &draw_query;
    &draw_subjects;
    &watermark;

    @unfused_one = split( /,/, $subject_hits[0] );
    @unfused_two = split( /,/, $subject_hits[1] );

    # We only really need 1dp for the folders, 2dp ratio is printed in image...
    $ratio = sprintf( "%.1f", $ratio );

    $query_accession_dir = "$GENE_HITS_DIFFERENTIALS/$query_accession";

    &print_comp_list("$query_accession,\n");
    &print_split_list("$unfused_one[0],\n$unfused_two[0],\n");

    if ( -e $query_accession_dir && -d $query_accession_dir ) {

        $ratio_dir = "$query_accession_dir/$ratio";
        if ( -e $ratio_dir && -d $ratio_dir ) {

            open my $out_fh, '>', "$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
        else {
            mkdir( $ratio_dir, 0755 );
            open $out_fh, '>', "$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
    }
    else {

        mkdir( $query_accession_dir, 0755 );
        $ratio_dir = "$query_accession_dir/$ratio";
        if ( -e $ratio_dir && -d $ratio_dir ) {

            open $out_fh, '>', "$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
        else {
            mkdir( $ratio_dir, 0755 );
            open $out_fh, '>', "$ratio_dir/$query_accession\_\_$unfused_one[0]\_\_$unfused_two[0].png";
            binmode $out_fh;
            print $out_fh $im->png(9);
            close($out_fh);
        }
    }
}

sub draw_query {
    $bar_length = $query_length;
    $accession  = $query_accession;
    &scale_bars;

    #
    $end       = $query_length;
    $end_scale = $query_length;
    &length_bars;

    #
    &bar($blue);
    &accession_name($blue);
    &draw_domain;
}

sub draw_subjects {

    for ( $c = 0 ; $c < $number_of_hits ; $c++ ) {

        @temp          = split( /,/, $subject_hits[$c] );
        $left_pos      = $left_pos + $temp[1];
        $bar_vpos      = $bar_vpos + 50;
        $bar_length    = ( $temp[2] - $temp[1] );
        $actual_length = $temp[3];
        $font_vpos     = $font_vpos + 50;

        #
        $colour_ratio = $bar_length / $actual_length;
        $colour_ratio = sprintf( "%.2f", $colour_ratio );

        &get_colour;

        #
        &bar($ratio_colour);

        #
        $start     = $temp[1];
        $end       = $bar_length;
        $end_scale = $end + $start;
        &scale_bars;

        #
        &length_bars;

        #
        $accession = $temp[0];
        $left_pos  = $left_pos - $temp[1];
        ##$left_pos = $padding_left + $name_length;
        $accession_vpos = $accession_vpos + 50;
        &accession_name($ratio_colour);
        ###############################
        # I think this section fixes the second domain out of sync bug!
        if ( $c == 0 ) {

            #print "***C*** = $c\n";
            &draw_domain;
        }
        elsif ( $c == 1 ) {

            #$left_pos = $padding_left + $name_length + 400;
            $left_pos = $padding_left + $temp[1] + $name_length;

            #print "---C--- = $c\n";
            &draw_domain;
        }
        ###############################
    }
}

sub draw_domain {

    for ( my $i = 0 ; $i <= $#domains ; $i++ ) {

        my @temp = split( /\t/, $domains[$i] );

        if ( $accession eq $temp[0] ) {
            print "$accession = $temp[0]\n";

            my @temp2 = split( /\~/, $temp[2] );
            for ( my $j = 0 ; $j <= $#temp2 ; $j++ ) {

                $var = $j * 3 + 3;

                $domain       = $temp2[$j];
                $domain_start = $temp[$var];
                $domain_end   = $temp[ $var + 1 ];
                $domain_type  = $temp[ $var + 2 ];
                chomp($domain_type);    # Remove end of line, causing last domain to be ignored.

                print "\t$#temp2 - $j - $domain\t$domain_start\t$domain_end\t$domain_type\n";

                &domain_colour($domain);
                if ( $domain_type eq ".." ) {
                    &domain_start_open;
                    &domain_middle;
                    &domain_end_open;
                }
                elsif ( $domain_type eq "[." ) {
                    &domain_start_closed;
                    &domain_middle;
                    &domain_end_open;
                }
                elsif ( $domain_type eq ".]" ) {
                    &domain_start_open;
                    &domain_middle;
                    &domain_end_closed;
                }
                elsif ( $domain_type eq "[]" ) {
                    &domain_start_closed;
                    &domain_middle;
                    &domain_end_closed;
                }
            }
        }
    }
}

sub domain_end_open {

    # make a polygon
    $hpos = $left_pos + $domain_end - 12;

    $poly = new GD::Polygon;
    $poly->addPt( $hpos,     $bar_vpos );
    $poly->addPt( $hpos + 4, $bar_vpos - 4 );
    $poly->addPt( $hpos,     $bar_vpos - 4 );
    $poly->addPt( $hpos,     $bar_vpos );
    $poly->addPt( $hpos,     $bar_vpos + 8 );
    $poly->addPt( $hpos + 4, $bar_vpos + 4 );
    $poly->addPt( $hpos,     $bar_vpos );

    # draw the polygon, filling it with a color
    $im->filledPolygon( $poly, $dom_colour );
}

sub domain_start_open {

    # make a polygon
    $hpos = $left_pos + $domain_start + 9;

    $poly = new GD::Polygon;
    $poly->addPt( $hpos + 4, $bar_vpos - 4 );    #40
    $poly->addPt( $hpos,     $bar_vpos - 4 );    #00
    $poly->addPt( $hpos + 4, $bar_vpos - 4 );    #44
    $poly->addPt( $hpos,     $bar_vpos );        #08
    $poly->addPt( $hpos + 4, $bar_vpos + 4 );    #412
    $poly->addPt( $hpos,     $bar_vpos + 8 );    #016
    $poly->addPt( $hpos + 4, $bar_vpos + 8 );    #416

    # draw the polygon, filling it with a color
    $im->filledPolygon( $poly, $dom_colour );

}

sub domain_end_closed {

    $im->setAntiAliased($dom_colour);
    $im->filledArc( $left_pos + $domain_end - 13, $bar_vpos + 2, 13, 13, 270, 90, $dom_colour, gdArc );
}

sub domain_start_closed {

    $im->filledArc( $left_pos + $domain_start + 14, $bar_vpos + 2, 13, 13, 90, 270, $dom_colour, gdArc );
    $im->setAntiAliased($dom_colour);
}

sub domain_middle {
    $im->filledRectangle(
        $left_pos + $domain_start + 14,
        $bar_vpos - 4,
        $left_pos + $domain_end - 13,
        $bar_vpos + 8, $dom_colour
    );    # + 14/-13 for caps
    $domain_name_length = length($domain) * 3;
    $dark_grey = $im->colorAllocate( 83, 83, 83 );
    if ( $domain_name_length < ( $domain_end - $domain_start ) ) {
        $im->string( gdSmallFont,
            $left_pos + $domain_start + ( ( ( $domain_end - $domain_start ) / 2 ) - $domain_name_length ) - 1,
            $bar_vpos - 4 - 1,
            "$domain", $dark_grey
        );    # shadow
        $im->string( gdSmallFont,
            $left_pos + $domain_start + ( ( ( $domain_end - $domain_start ) / 2 ) - $domain_name_length ),
            $bar_vpos - 4,
            "$domain", $white
        );
    }
}

sub domain_colour {
    my $temp = shift;

    # Add the ascii values of the domain text to generate
    # colours specific to each name...
    my @domain = unpack( "C*", "$temp" );
    $sum = eval { join '+', @domain };
    $sum = $sum / 2;

    my $r = int( ( ( $sum * 10 ) - 40 ) - 240 );
    my $g = int( ( ( $sum * 10 ) - 40 ) - 120 );
    my $b = int( ( ( $sum * 10 ) - 40 ) - 80 );
    my $a = "45";

    #$dom_colour = $im->colorAllocate( $r, $g, $b );
    $dom_colour = $im->colorAllocateAlpha( $r, $g, $b, $a );

    return $dom_colour;
}

sub scale_bars {
    $scale_quotient = int( $bar_length / 100 );

    #Large Scale Bars
    for ( my $a = 0 ; $a <= $scale_quotient ; $a++ ) {
        $scale = $a * 100;
        $im->string( gdSmallFont, $scale + $left_pos, $font_vpos, "$scale", $black );
        $im->rectangle( $scale + $left_pos, $bar_vpos - 5, $scale + $left_pos + 1, $bar_vpos - 1, $black );
    }
}

sub length_bars {

    $end_padded = $end + $left_pos;

    # Start
    $im->string( gdSmallFont, $left_pos, $font_vpos + 30, "$start", $red );
    $im->rectangle( $left_pos, $bar_vpos + 6, $left_pos + 1, $bar_vpos + 10, $red );

    #End
    $im->string( gdSmallFont, $end_padded, $font_vpos + 30, "$end_scale", $red );
    $im->rectangle( $end_padded - 1, $bar_vpos + 6, $end_padded, $bar_vpos + 10, $red );

    #Length End
    #$im->string( gdSmallFont, $end_padded, $font_vpos, "$bar_length", $black );
    #$im->rectangle( $end_padded - 1, $bar_vpos - 5, $end_padded, $bar_vpos - 1, $black );

}

sub bar {
    $colour = shift;
    $im->filledRectangle( $left_pos, $bar_vpos, $bar_length + $left_pos, $bar_vpos + 5, $colour );
}

sub accession_name {
    my $colour = shift;
    if ( $actual_length eq "" ) {
        $im->string( gdSmallFont, $padding_left - 50, $accession_vpos,      "$accession",    $colour );
        $im->string( gdTinyFont,  $padding_left - 50, $accession_vpos + 12, "L=$bar_length", $colour );
    }
    else {
        $im->string( gdSmallFont, $padding_left - 50, $accession_vpos,      "$accession",        $colour );
        $im->string( gdTinyFont,  $padding_left - 50, $accession_vpos + 20, "ML=$bar_length",    $colour );
        $im->string( gdTinyFont,  $padding_left - 50, $accession_vpos + 12, "L =$actual_length", $colour );
    }
}

sub watermark {

    $dark_grey  = $im->colorAllocate( 83,  83,  83 );
    $light_grey = $im->colorAllocate( 192, 192, 192 );

    #$dark_red   = $im->colorAllocate( 146, 84,  83 );
    #$dark_green = $im->colorAllocate( 129, 158, 107 );

    $darker_grey  = $im->colorAllocate( 52,  50,  51 );
    $darker_red   = $im->colorAllocate( 123, 59,  59 );
    $purple       = $im->colorAllocate( 96,  84,  112 );
    $dark_blue    = $im->colorAllocate( 33,  135, 204 );
    $darker_green = $im->colorAllocate( 43,  166, 22 );
    $pink         = $im->colorAllocate( 255, 0,   255 );

    #$logo_left_x = $image_length - $padding_right;
    $logo_left_y = $image_height - 20;

    $fdfBLAST = "fdfBLAST v$VERSION";

    $im->string( gdSmallFont, $image_length - 100, $logo_left_y + 7, "$fdfBLAST", $light_grey );

    $im->string( gdSmallFont, $image_length - 180, $logo_left_y + 7, "Ratio = $ratio", $light_grey );

    $im->string( gdTinyFont, $image_length - 425, $logo_left_y + 10, "Length Match %", $light_grey );

    $im->filledRectangle( $image_length - 350, $logo_left_y + 10, $image_length - 320, $logo_left_y + 17,
        $darker_grey );
    $im->filledRectangle( $image_length - 320, $logo_left_y + 10, $image_length - 290, $logo_left_y + 17, $darker_red );
    $im->filledRectangle( $image_length - 290, $logo_left_y + 10, $image_length - 260, $logo_left_y + 17, $purple );
    $im->filledRectangle( $image_length - 260, $logo_left_y + 10, $image_length - 230, $logo_left_y + 17, $dark_blue );
    $im->filledRectangle(
        $image_length - 230,
        $logo_left_y + 10,
        $image_length - 200,
        $logo_left_y + 17,
        $darker_green
    );
    $im->filledRectangle( $image_length - 200, $logo_left_y + 10, $image_length - 190, $logo_left_y + 17, $pink );

    $im->string( gdTinyFont, $image_length - 343, $logo_left_y + 10, "<40",    $light_grey );
    $im->string( gdTinyFont, $image_length - 317, $logo_left_y + 10, "40-60",  $light_grey );
    $im->string( gdTinyFont, $image_length - 287, $logo_left_y + 10, "60-70",  $light_grey );
    $im->string( gdTinyFont, $image_length - 257, $logo_left_y + 10, "70-80",  $light_grey );
    $im->string( gdTinyFont, $image_length - 229, $logo_left_y + 10, "80-100", $light_grey );
    $im->string( gdTinyFont, $image_length - 197, $logo_left_y + 10, "!",      $light_grey );
}

sub get_colour {

    if ( $colour_ratio le "0.4" ) {

        # Black
        my $r = int( 42 - ( 40 - ( 10 * ( 10 * 1 ) ) ) );
        my $g = int( 40 - ( 40 - ( 10 * ( 10 * 1 ) ) ) );
        my $b = int( 41 - ( 40 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Black $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.4" && $colour_ratio le "0.6" ) {

        # Red
        my $r = int( 123 - ( 60 - ( 10 * ( 10 * 1 ) ) ) );
        my $g = int( 59 -  ( 60 - ( 10 * ( 10 * 1 ) ) ) );
        my $b = int( 59 -  ( 60 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Red $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.6" && $colour_ratio le "0.7" ) {

        # Purple
        my $r = int( 96 -  ( 70 - ( 10 * ( 10 * 1 ) ) ) );
        my $g = int( 84 -  ( 70 - ( 10 * ( 10 * 1 ) ) ) );
        my $b = int( 112 - ( 70 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Purple $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.7" && $colour_ratio le "0.8" ) {

        # Blue
        my $r = int( 33 -  ( 80 - ( 10 * ( 10 * 1 ) ) ) );
        my $g = int( 135 - ( 80 - ( 10 * ( 10 * 1 ) ) ) );
        my $b = int( 204 - ( 80 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Blue $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    elsif ( $colour_ratio gt "0.8" && $colour_ratio le "1.0" ) {

        # Green
        my $r = int( 43 -  ( 100 - ( 10 * ( 10 * 1 ) ) ) );
        my $g = int( 166 - ( 100 - ( 10 * ( 10 * 1 ) ) ) );
        my $b = int( 22 -  ( 100 - ( 10 * ( 10 * 1 ) ) ) );

        #print "Green $r, $g, $b\n";
        $ratio_colour = $im->colorAllocate( $r, $g, $b );
    }
    else {

        # Bright Pink
        $ratio_colour = $im->colorAllocate( 255, 0, 255 );
    }
    return $ratio_colour;
}

sub ignore_orthologues {

    print "Ignoring Potential Orthologues\n" if $DEBUG >= 1;
    $continue = "yes";

    for ( my $a = 0 ; $a < $middle_num ; $a++ ) {
        @unfused_middles = split( /,/, $middles[$a] );
        my $middle_length = $unfused_middles[3];
        print "\t$query_accession - $middle_length\n" if $DEBUG >= 1;
        $length_ratio = $middle_length / $query_length;

        if ( $length_ratio le "0.95" ) {

            #$continue = "yes";
            push( @cont, "yes" );
            print "\t\t$unfused_middles[0] - $middle_length / $query_length = $length_ratio - yes\n"
              if $DEBUG >= 1;
        }
        elsif ( $length_ratio gt "0.95" ) {

            push( @cont, "no" );
            print "\t\t$unfused_middles[0] - $middle_length / $query_length = $length_ratio - no\n"
              if $DEBUG >= 1;
        }

    }
    foreach my $item (@cont) {
        print $item . "," if $DEBUG >= 1;
    }
    print "\n" if $DEBUG >= 1;

    $continue = "no" if ( grep { /^no$/ } @cont );

    undef(@cont);
    return $continue;
}

sub print_comp_list {
    my $message = shift;
    open my $log_fh, '>>', "$RUN_DIR/composite_list.csv";
    print $log_fh "$message";
    close($log_fh);
}

sub print_split_list {
    my $message = shift;
    open my $log_fh, '>>', "$RUN_DIR/split_list.csv";
    print $log_fh "$message";
    close($log_fh);
}


sub read_domain_file {
    open my $domains_fh, '<', $domain_in;
    while (<$domains_fh>) {
        $line = $_;
        push( @domains, $line );
    }
    close($domains_fh);
}

sub dhms {

    my %months = (
        Jan => 1,
        Feb => 2,
        Mar => 3,
        Apr => 4,
        May => 5,
        Jun => 6,
        Jul => 7,
        Aug => 8,
        Sep => 9,
        Oct => 10,
        Nov => 11,
        Dec => 12
    );

    my $beginning = $start_time;
    my $end       = $end_time;

    my @b = split( /\s+|:/sm, $beginning );
    my @e = split( /\s+|:/sm, $end );

    my $b = timelocal( $b[5], $b[4], $b[3], $b[2], $months{ $b[1] } - 1, $b[-1] );
    my $e = timelocal( $e[5], $e[4], $e[3], $e[2], $months{ $e[1] } - 1, $e[-1] );

    my $elapsed = $e - $b;

    @parts = gmtime($elapsed);

    return @parts;

}
