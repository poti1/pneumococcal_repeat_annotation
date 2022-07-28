#!/bin/env perl

use v5.26;
use strict;
use warnings FATAL => "all";
use Bio::SeqIO;
use Bio::SeqI;


=head1 NAME

    pneumococcal_repeat_annotation

=head1 SYNOPSIS

    perl pneumococcal_repeat_annotation.pl my.fasta

    pneumococcal_repeat_annotation.pl Fasta_Files/ERR2090225.fasta

=head1 DESCRIPTION

Script to process repeat sequences to realign overlapping elements
and identify disrupted repeat elements.

=cut


##############################################################################
#                                Main
##############################################################################

my $DEBUG = 0;
my $n     = 0;
my %start;
my %end;
my %type;
my %score;
my %strand;
my $r = 0;
my %BOX;
my $boxstart;
my $boxend;
my $boxnum = 0;
my $boxstrand;
my %boxstart;
my %boxstrand;
my %boxend;
my $A_start;
my $A_end;
my $A_score = 0;
my $A_hitstrand;
my $A_hitstart;
my $A_hitend;
my $B_start;
my $B_end;
my $B_score = 0;
my $B_hitstrand;
my $B_hitstart;
my $B_hitend;
my $in_box = 0;
my $start;
my $end;
my $hitscore  = 0;
my $hitstart  = 0;
my $hitend    = 0;
my $hitstrand = 0;
my %brokencoords;
my %brokenscore;
my %brokendisrupt;
my $hmmsearch = "/home/sreeram/hmmer-1.8.4/hmmls -c -t";
$hmmsearch = "hmmls -c -t";

my $repeat_hmm_dir = "Repeat_HMMs";

my $repeats = {
    boxA => {
        cutoff     => 30,
        hmm        => "$repeat_hmm_dir/boxA.hmm",    # paths of HMM for analysis
        input_file => "boxA.rep",    # output files from HMM analysis
    },
    boxB => {
        cutoff     => 14,
        hmm        => "$repeat_hmm_dir/boxB.hmm",
        input_file => "boxB.rep",
    },
    boxC => {
        cutoff     => 28,
        hmm        => "$repeat_hmm_dir/boxC.hmm",
        input_file => "boxC.rep",
    },
    RUP => {
        cutoff     => 47,
        hmm        => "$repeat_hmm_dir/RUP.hmm",
        input_file => "RUP.rep",
    },
    SPRITE => {
        cutoff     => 66,
        hmm        => "$repeat_hmm_dir/SPRITE.hmm",
        input_file => "SPRITE.rep",
    },
};

run_sequences( $repeats ) if not $DEBUG;

foreach my $genome ( @ARGV ) {

    build_hmm_structure( $repeats, $genome );

    print STDERR "Identified repeat units in sequence $genome\n";
    identify_composite_box_elements();

    find_box_elements_with_no_starting_boxA_module();

    print STDERR "Identified composite BOX elements\n";
    delete( $BOX{$boxnum} );

    my %overlaps;
    identify_overlapping_repeat_elements( \%overlaps ) if defined $start;
    print STDERR "Realigning overlapping elements\n";

    my $fasta    = Bio::SeqIO->new( -file => $genome, -format => "fasta" );
    my $sequence = $fasta->next_seq;

    say "About to process overlaps";
    rescan_overlapping_repeats( \%overlaps, $sequence );

    print STDERR "Disrupted repeat elements identified\n";
    write_table( $genome );

    undef( %start );
    undef( %end );
    undef( %type );
    undef( %score );
    undef( %strand );
    undef( %BOX );
    undef( %boxstart );
    undef( %boxstrand );
    undef( %boxend );
    undef( %brokencoords );
    undef( %brokenscore );
    undef( %brokendisrupt );

}

close OUT;

# cleanup();  # TODO: put cleanup back.
print STDERR "All Done\n";


##############################################################################
#                          Functions
##############################################################################

sub run_sequences {
    my ( $repeats ) = @_;

    foreach my $seq ( @ARGV ) {
        foreach my $repeat ( sort keys %$repeats ) {
            my $node = $repeats->{$repeat};

            print STDERR "Searching sequence $seq for repeat $repeat...";
            system(
"$hmmsearch $node->{cutoff} $node->{hmm} $seq > $seq.$node->{input_file}"
            );
            print STDERR "done\n";
        }
    }
}

# parse information from HMM output
sub build_hmm_structure {
    my ( $repeats, $genome ) = @_;

    my $is_score = qr{ ^ \s* ( [\d.]+ ) }x;
    my $is_start = qr{ \b f \b \s* : \s* ( \d+ ) }x;
    my $is_end   = qr{ \b t \b \s* : \s* ( \d+ ) }x;

    foreach my $repeat ( sort keys %$repeats ) {
        my $node = $repeats->{$repeat};
        my $file = "$genome.$node->{input_file}";
        open my $fh, "<", $file or die "Cannot open input file: $file\n";

        foreach ( <$fh> ) {
            chomp;
            say "\nline: [$_]" if $DEBUG >= 2;

            my ( $score ) = /$is_score/;
            next if not $score;

            my ( $start ) = /$is_start/;
            if ( !$start ) {
                print STDERR "Missing start value for: $genome line $.\n";
                next;
            }

            my ( $end ) = /$is_end/;
            if ( !$end ) {
                print STDERR "Missing end value for: $genome line $.\n";
                next;
            }

            if ( $score > $node->{cutoff} ) {
                $start{$r}  = $start;
                $end{$r}    = $end;
                $type{$r}   = $repeat;
                $score{$r}  = $score;
                $strand{$r} = ( $end > $start ) ? 1 : -1;

                # say dumper {
                #     _00_r      => $r,
                #     _21_start  => $start{$r},
                #   # _22_startH => \%start,
                #     _24_end    => $end{$r},
                #   # _25_endH   => \%end,
                #     _31_type   => $type{$r},
                #     _32_score  => $score{$r},
                #     _33_strand => $strand{$r},
                # } if $DEBUG;
            }
            $r++;
        }

        close $fh;
    }
}

# identify composite BOX elements from box modules
sub identify_composite_box_elements {
    foreach my $repeata ( sort keys %start ) {
        if ( $type{$repeata} =~ /boxA/ ) {
            $boxstrand = $strand{$repeata};
            $boxstart  = $start{$repeata};
            @{ $BOX{$boxnum} } = "$repeata";
            if ( $boxstrand == 1 ) {
                $boxend = $end{$repeata} +
                  5;    # five base leeway for finding the adjacent boxB element
                foreach
                  my $repeatb ( sort { $start{$a} <=> $start{$b} } keys %start )
                {
                    if (
                        $type{$repeatb} =~ /boxB/
                        && (   $start{$repeatb} <= $boxend
                            && $start{$repeatb} >= $boxstart )
                        && $strand{$repeatb} == $boxstrand
                      )
                    {
                        $boxend = $end{$repeatb} + 1;
                        push( @{ $BOX{$boxnum} }, $repeatb );
                    }
                }
                foreach
                  my $repeatc ( sort { $start{$a} <=> $start{$b} } keys %start )
                {
                    if (
                        $type{$repeatc} =~ /boxC/
                        && (   $start{$repeatc} <= $boxend
                            && $start{$repeatc} >= $boxstart )
                        && $strand{$repeatc} == $boxstrand
                      )
                    {
                        $boxend = $end{$repeatc} + 1;
                        push( @{ $BOX{$boxnum} }, $repeatc );
                    }
                }
                if ( $#{ $BOX{$boxnum} } > 0 ) {
                    $boxstart{$boxnum}  = $boxstart;
                    $boxend{$boxnum}    = $boxend - 1;
                    $boxstrand{$boxnum} = $boxstrand;
                    $boxnum++;
                }
            }
            else {
                $boxend = $end{$repeata} - 5;
                foreach
                  my $repeatb ( sort { $start{$b} <=> $start{$a} } keys %start )
                {
                    if (
                        $type{$repeatb} =~ /boxB/
                        && (   $start{$repeatb} >= $boxend
                            && $start{$repeatb} <= $boxstart )
                        && $strand{$repeatb} == $boxstrand
                      )
                    {
                        $boxend = $end{$repeatb} - 1;
                        push( @{ $BOX{$boxnum} }, $repeatb );
                    }
                }
                foreach
                  my $repeatc ( sort { $start{$b} <=> $start{$a} } keys %start )
                {
                    if (
                        $type{$repeatc} =~ /boxC/
                        && (   $start{$repeatc} >= $boxend
                            && $start{$repeatc} <= $boxstart )
                        && $strand{$repeatc} == $boxstrand
                      )
                    {
                        $boxend = $end{$repeatc} - 1;
                        push( @{ $BOX{$boxnum} }, $repeatc );
                    }
                }
                if ( $#{ $BOX{$boxnum} } > 0 ) {
                    $boxstart{$boxnum}  = $boxstart;
                    $boxend{$boxnum}    = $boxend + 1;
                    $boxstrand{$boxnum} = $boxstrand;
                    $boxnum++;
                }
            }
        }
    }
}

# find BC BOX elements with no starting boxA module
sub find_box_elements_with_no_starting_boxA_module {
    $in_box = 0;
    foreach my $repeatc ( sort keys %start ) {
        if ( $type{$repeatc} =~ /boxC/ ) {
            foreach $boxnum ( sort keys %BOX ) {
                if ( grep( /$repeatc/, @{ $BOX{$boxnum} } ) ) {
                    $in_box = 1;
                }
            }
            if ( $in_box == 0 ) {
                $boxend    = $end{$repeatc};
                $boxstrand = $strand{$repeatc};
                @{ $BOX{$boxnum} } = "$repeatc";
                if ( $boxstrand == -1 ) {
                    $boxstart = $start{$repeatc} + 1;
                    foreach my $repeatb (
                        sort { $start{$a} <=> $start{$b} }
                        keys %start
                      )
                    {
                        if (
                            $type{$repeatb} =~ /boxB/
                            && (   $end{$repeatb} <= $boxstart
                                && $start{$repeatb} >= $boxstart )
                            && $strand{$repeatb} == $boxstrand
                          )
                        {
                            $boxstart = $start{$repeatb} + 1;
                            push( @{ $BOX{$boxnum} }, $repeatb );
                        }
                    }
                    if ( $#{ $BOX{$boxnum} } > 0 ) {
                        $boxstart{$boxnum}  = $boxstart - 1;
                        $boxend{$boxnum}    = $boxend;
                        $boxstrand{$boxnum} = $boxstrand;
                        $boxnum++;
                    }
                }
                else {
                    $boxstart = $start{$repeatc} - 1;
                    foreach my $repeatb (
                        sort { $start{$b} <=> $start{$a} }
                        keys %start
                      )
                    {
                        if (
                            $type{$repeatb} =~ /boxB/
                            && (   $end{$repeatb} >= $boxstart
                                && $start{$repeatb} <= $boxstart )
                            && $strand{$repeatb} == $boxstrand
                          )
                        {
                            $boxstart = $start{$repeatb} - 1;
                            push( @{ $BOX{$boxnum} }, $repeatb );
                        }
                    }
                    if ( $#{ $BOX{$boxnum} } > 0 ) {
                        $boxstart{$boxnum}  = $boxstart + 1;
                        $boxend{$boxnum}    = $boxend;
                        $boxstrand{$boxnum} = $boxstrand;
                        $boxnum++;
                    }
                }
            }
            $in_box = 0;
        }
    }
}

# identify overlapping repeat elements, except BOX modules
sub identify_overlapping_repeat_elements {
    my ( $overlaps ) = @_;

    foreach my $repeata ( sort keys %start ) {
        foreach my $repeatb ( sort keys %start ) {
            if (
                (
                       $start{$repeatb} <= $start{$repeata}
                    && $end{$repeatb} >= $start{$repeata}
                )
                || (   $start{$repeatb} <= $end{$repeata}
                    && $end{$repeatb} >= $end{$repeata} )
              )
            {
                unless (
                       ( $type{$repeata} =~ /box/ && $type{$repeatb} =~ /box/ )
                    || ( $repeata == $repeatb ) )
                {
                    my @unsorted = ( $repeata, $repeatb );
                    my @sorted   = sort( @unsorted );
                    $overlaps->{ $sorted[0] } = $sorted[1];
                }
            }
        }
    }
}

# rescan regions around overlapping repeats to look for disruptions
sub rescan_overlapping_repeats {
    my ( $overlaps, $sequence ) = @_;

    foreach my $repeata ( sort keys %$overlaps ) {
        boundaries( $repeata );
        $A_start = $start;
        $A_end   = $end;
        my $upper_bound = $end + 251;
        my $lower_bound = $start - 251;
        my $lower_seq   = $sequence->subseq( $lower_bound, $start - 1 );
        my $upper_seq   = $sequence->subseq( $end + 1,     $upper_bound );
        my $total_seq   = "$lower_seq" . "$upper_seq";
        print_fasta( "first.seq", $total_seq );
        system(
"$hmmsearch $repeats->{$type{$overlaps->{$repeata}}}->{hmmh} first.seq > first.out"
        );
        hmm_results( "first.out" );
        $A_score     = $hitscore;
        $A_hitstrand = $hitstrand;
        $A_hitstart  = $hitstart;
        $A_hitend    = $hitend;
        boundaries( $overlaps->{$repeata} );
        $B_start     = $start;
        $B_end       = $end;
        $upper_bound = $end + 251;
        $lower_bound = $start - 251;
        $lower_seq   = $sequence->subseq( $lower_bound, $start - 1 );
        $upper_seq   = $sequence->subseq( $end + 1,     $upper_bound );
        $total_seq   = "$lower_seq" . "$upper_seq";
        print_fasta( "second.seq", $total_seq );
        system(
"$hmmsearch $repeats->{$type{$repeata}}->{hmm} second.seq > second.out"
        );
        hmm_results( "second.out" );
        $B_score     = $hitscore;
        $B_hitstrand = $hitstrand;
        $B_hitstart  = $hitstart;
        $B_hitend    = $hitend;
        my $A_inc = ( $A_score - $score{ $overlaps->{$repeata} } ) /
          $score{ $overlaps->{$repeata} };
        my $B_inc = ( $B_score - $score{$repeata} ) / $score{$repeata};

        if ( $A_inc >= 0 || $B_inc >= 0 ) {    # identify disruption events
            my $firstpoint;
            my $secondpoint;
            my $thirdpoint;
            my $fourthpoint;
            if ( $A_inc > $B_inc ) {
                $start{ $overlaps->{$repeata} } = "BROKEN";
                if ( $A_hitstrand == 1 ) {
                    $secondpoint = $A_start - 1;
                    $thirdpoint  = $A_end + 1;
                    $firstpoint  = $A_start - 252 + $A_hitstart;
                    $fourthpoint = $A_end - 251 + $A_hitend;
                    @{ $brokencoords{ $overlaps->{$repeata} } } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{$repeata} =~ /box/ ) {
                        $brokendisrupt{ $overlaps->{$repeata} } = "BOX";
                    }
                    else {
                        $brokendisrupt{ $overlaps->{$repeata} } =
                          $type{$repeata};
                    }
                    $brokenscore{ $overlaps->{$repeata} } = $A_score;
                }
                else {
                    $secondpoint = $A_start - 1;
                    $thirdpoint  = $A_end + 1;
                    $firstpoint  = $A_start - 252 + $A_hitend;
                    $fourthpoint = $A_end - 251 + $A_hitstart;
                    @{ $brokencoords{ $overlaps->{$repeata} } } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{$repeata} =~ /box/ ) {
                        $brokendisrupt{ $overlaps->{$repeata} } = "BOX";
                    }
                    else {
                        $brokendisrupt{ $overlaps->{$repeata} } =
                          $type{$repeata};
                    }
                    $brokenscore{ $overlaps->{$repeata} } = $A_score;
                }
            }
            elsif ( $B_inc > $A_inc ) {
                $start{$repeata} = "BROKEN";
                if ( $B_hitstrand == 1 ) {
                    $secondpoint = $B_start - 1;
                    $thirdpoint  = $B_end + 1;
                    $firstpoint  = $B_start - 252 + $B_hitstart;
                    $fourthpoint = $B_end - 251 + $B_hitend;
                    @{ $brokencoords{$repeata} } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{ $overlaps->{$repeata} } =~ /box/ ) {
                        $brokendisrupt{$repeata} = "BOX";
                    }
                    else {
                        $brokendisrupt{$repeata} =
                          $type{ $overlaps->{$repeata} };
                    }
                    $brokenscore{$repeata} = $B_score;
                }
                else {
                    $secondpoint = $B_start - 1;
                    $thirdpoint  = $B_end + 1;
                    $firstpoint  = $B_start - 252 + $B_hitend;
                    $fourthpoint = $B_end - 251 + $B_hitstart;
                    @{ $brokencoords{$repeata} } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{ $overlaps->{$repeata} } =~ /box/ ) {
                        $brokendisrupt{$repeata} = "BOX";
                    }
                    else {
                        $brokendisrupt{$repeata} =
                          $type{ $overlaps->{$repeata} };
                    }
                    $brokenscore{$repeata} = $B_score;
                }
            }
        }

    }
}

sub write_table {
    my ( $genome ) = @_;

    open my $fh, ">", "$genome.repeats.tab" or die "Cannot open output file\n";
    print STDERR "Printing output files\n";

    foreach my $repeat ( sort keys %start ) {
        my $coords = $brokencoords{$repeat};

        if ( $start{$repeat} eq "BROKEN" ) {
            my $type = ( $strand{$repeat} == 1 ) ? "order" : "complement";
            print $fh <<~"OUTPUT";
            FT   repeat_unit     $type($coords->[0]..$coords->[1],$coords->[2]..$coords->[3])
            FT                   /colour=2
            FT                   /label=$type{$repeat}
            FT                   /note=Detected using HMMER $hmmsearch; appears to have been disrupted through $brokendisrupt{$repeat} insertion
            FT                   /note=Initial match of score $score{$repeat} to model $type{$repeat}; realignment score of $brokenscore{$repeat}
            OUTPUT
        }
        else {
            my $range = "$start{$repeat}..$end{$repeat}";
            my $type_and_range =
              ( $strand{$repeat} == 1 ) ? $range : "complement($range)";
            print $fh <<~"OUTPUT";
            FT   repeat_unit     $type_and_range
            FT                   /colour=2
            FT                   /label=$type{$repeat}
            FT                   /note=Detected using HMMER $hmmsearch; match of score $score{$repeat} to model $type{$repeat}
            OUTPUT
        }
    }

    foreach my $boxnum ( sort keys %BOX ) {
        my $range = "$boxstart{$boxnum}..$boxend{$boxnum}";
        my $type_and_range =
          ( $boxstrand{$boxnum} == 1 ) ? $range : "complement($range)";
        print $fh <<~"OUTPUT";
        FT   repeat_unit     $type_and_range
        FT                   /colour=4
        FT                   /note=Composite BOX element
        FT                   /label=BOX
        OUTPUT
    }

    close $fh;
}

# subroutine for identifying repeat boundaries, esp BOX elements
sub boundaries {
    my ( $rep ) = @_;

    if ( $type{$rep} =~ /box/ ) {
        foreach $boxnum ( sort keys %BOX ) {
            foreach my $module ( @{ $BOX{$boxnum} } ) {
                if ( $module == $rep ) {
                    my @unsorted = ( "$boxstart{$boxnum}", "$boxend{$boxnum}" );
                    my @sorted   = sort( @unsorted );
                    $start  = $sorted[0];
                    $end    = $sorted[1];
                    $in_box = 1;
                }
            }
        }
    }
    if ( $in_box == 0 ) {
        my @unsorted = ( "$start{$rep}", "$end{$rep}" );
        my @sorted   = sort( @unsorted );
        $start = $sorted[0];
        $end   = $sorted[1];
    }
    $in_box = 0;
    return ( $start, $end );
}

# subroutine for printing sequence in a suitable format for /home/sreeram/hmmer-1.8.4/hmmls
sub print_fasta {
    my $filename   = shift;
    my $dna_string = shift;
    open OUT, "> $filename";
    print OUT ">$filename\n";
    my $offset = 0;
    while ( $offset < 440 ) {
        my $line = substr( $dna_string, $offset, 60 );
        print OUT "$line\n";
        $offset += 60;
    }
    my $line = substr( $dna_string, $offset, ( 500 - $offset ) );
    print OUT "$line\n";
    close OUT;
}

# subroutine for picking the top hit from the /home/sreeram/hmmer-1.8.4/hmmls results
sub hmm_results {
    my $filename = shift;
    open HMM, $filename, or die "Cannot open HMM file\n";
    $hitscore = 0;
    foreach ( <HMM> ) {
        my @data = split( /\s+/, $_ );
        if ( substr( $_, 0, 2 ) =~ /\d\d/ ) {
            if ( $data[0] > $hitscore ) {
                if ( $data[3] < $data[5] && $data[3] <= 250 && $data[5] >= 250 )
                {
                    $hitstrand = 1;
                    $hitscore  = $data[0];
                    $hitstart  = $data[3];
                    $hitend    = $data[5];
                }
                elsif ($data[5] < $data[3]
                    && $data[3] >= 250
                    && $data[5] <= 250 )
                {
                    $hitstrand = -1;
                    $hitscore  = $data[0];
                    $hitstart  = $data[3];
                    $hitend    = $data[5];
                }
            }
        }
    }
    close HMM;
    return ( $hitstart, $hitend, $hitscore, $hitstrand );
}

sub cleanup {
    my @files_to_remove =
      glob "first.out first.seq second.out second.seq *.rep";

    for ( @files_to_remove ) {
        if ( -e ) {
            unlink or warn "Cannot remove '$_': $!";
        }
    }
}

