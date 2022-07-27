#! /software/bin/perl -w

=head1 NAME

    pneumococcal_repeat_annotation

=head1 SYNOPSIS

    perl pneumococcal_repeat_annotation.pl my.fasta

=head1 DESCRIPTION

Script to process repeat sequences to realign overlapping elements
and identify disrupted repeat elements.

=cut

use v5.26;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqI;

#TODO: Remove this debug code !!!
use feature    qw(say);
use Mojo::Util qw(dumper);

my $r = 0;
my %BOX;
my $boxstartScalar;
my $boxEndScalar;
my $boxnum = 0;
my $boxstrand;
my %boxStartHash;
my %boxStrandHash;
my %boxEndHash;
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
my $startScalar;
my $endScalar;
my $hitscore  = 0;
my $hitstart  = 0;
my $hitend    = 0;
my $hitstrand = 0;
my %brokencoords;
my %brokenscore;
my %brokendisrupt;
my $hmmls_script = "/home/sreeram/hmmer-1.8.4/hmmls";
$hmmls_script = "hmmsearch";


########################################################################
#                             Main
########################################################################

my $options = {
    boxA => {
        cutoff     => 30,
        hmm        => "./Repeat_HMMs/boxA.hmm", # paths of HMM for analysis
        input_file => "boxA.rep",               # output files from HMM analysis
    },
    boxB => {
        cutoff     => 14,
        hmm        => "./Repeat_HMMs/boxB.hmm",
        input_file => "boxB.rep",
    },
    boxC => {
        cutoff     => 28,
        hmm        => "./Repeat_HMMs/boxC.hmm",
        input_file => "boxC.rep",
    },
    RUP => {
        cutoff     => 47,
        hmm        => "./Repeat_HMMs/RUP.hmm",
        input_file => "RUP.rep",
    },
    SPRITE => {
        cutoff     => 66,
        hmm        => "./Repeat_HMMs/SPRITE.hmm",
        input_file => "SPRITE.rep",
    },

};

sequence_search( $options );

for my $genome ( @ARGV ) {

    my $hmm_spec = parse_hmms( $genome, $options );

    print STDERR "Identified repeat units in sequence $genome\n";

    identify_composite_box_elements( $hmm_spec );

    find_box_elemets_with_no_starting_boxA_module( $hmm_spec );

    print STDERR "Identified composite BOX elements\n";
    delete( $BOX{$boxnum} );

    my %overlaps;

    if ( defined $startScalar ) {
        identify_overlapping_repeat_elements_except_box_modules( $hmm_spec,
            \%overlaps );
    }

    print STDERR "Realigning overlapping elements\n";

    my $fasta    = Bio::SeqIO->new( -file => $genome, -format => "fasta" );
    my $sequence = $fasta->next_seq;

    rescan_regions_around_overlapping_repeats_to_find_disruptions( $hmm_spec,
        \%overlaps, $sequence );

    print STDERR "Disrupted repeat elements identified\n";
    print STDERR "Printing output files\n";

    save_results( $hmm_spec, \%brokencoords, $genome, \%BOX );

    undef( %BOX );
    undef( %boxStartHash );
    undef( %boxStrandHash );
    undef( %boxEndHash );
    undef( %brokencoords );
    undef( %brokenscore );
    undef( %brokendisrupt );

}

cleanup();
print STDERR "Done\n";


########################################################################
#                           Functions
########################################################################

sub sequence_search {
    my ( $options ) = @_;

    for my $seq ( @ARGV ) {
        for my $repeat ( sort keys %$options ) {
            my $node = $options->{$repeat};

            print STDERR "Searching sequence $seq for repeat $repeat...";
            system(
"$hmmls_script -t $node->{cutoff} $node->{hmm} $seq > $seq.$node->{input_file}"
            );
            print STDERR "done\n";
        }
    }
}

sub parse_hmms {
    my ( $genome, $options ) = @_;
    my $hmm_spec = {};

    for my $repeat ( sort keys %$options ) { # parse information from HMM output
        open my $fh, "$genome.$options->{$repeat}{input_file}"
          or die "Cannot open input file\n";
        for ( <$fh> ) {
            chomp;
            if ( substr( $_, 0, 2 ) =~ /\d\d/ ) {
                my @array = split( /\s+/, $_ );
                my @start = split( /:/,   $array[2] );
                my @end   = split( /:/,   $array[3] );
                if ( $array[0] > $options->{$repeat}{cutoff} ) {
                    $hmm_spec->{$r} = {
                        type   => $repeat,
                        score  => $array[0],
                        start  => $start[1],
                        end    => $end[1],
                        strand => $end[1] > $start[1] ? 1 : -1,
                    };
                }
                $r++;
            }
        }
    }

    $hmm_spec;
}

# identify composite BOX elements from box modules
sub identify_composite_box_elements {
    my ( $hmm_spec ) = @_;

    for my $repeata ( sort keys %$hmm_spec ) {
        my $nodeA = $hmm_spec->{$repeata};
        if ( $nodeA->{type} =~ /boxA/ ) {
            $boxstrand      = $nodeA->{strand};
            $boxstartScalar = $nodeA->{start};
            @{ $BOX{$boxnum} } = $repeata;

            if ( $boxstrand == 1 ) {

                # five base leeway for finding the adjacent boxB element
                $boxEndScalar = $nodeA->{end} + 5;
                for my $repeatb (
                    sort { $hmm_spec->{$a}{start} <=> $hmm_spec->{$b}{start} }
                    keys %$hmm_spec
                  )
                {
                    if (
                        $hmm_spec->{$repeatb}->{type} =~ /boxB/
                        && (   $hmm_spec->{$repeatb}{start} <= $boxEndScalar
                            && $hmm_spec->{$repeatb}{start} >= $boxstartScalar )
                        && $hmm_spec->{$repeatb}->{strand} == $boxstrand
                      )
                    {
                        $boxEndScalar = $hmm_spec->{$repeatb}->{end} + 1;
                        push( @{ $BOX{$boxnum} }, $repeatb );
                    }
                }
                for my $repeatc (
                    sort { $hmm_spec->{$a}{start} <=> $hmm_spec->{$b}{start} }
                    keys %$hmm_spec
                  )
                {
                    if (
                        $hmm_spec->{$repeatc}->{type} =~ /boxC/
                        && (   $hmm_spec->{$repeatc}{start} <= $boxEndScalar
                            && $hmm_spec->{$repeatc}{start} >= $boxstartScalar )
                        && $hmm_spec->{$repeatc}->{strand} == $boxstrand
                      )
                    {
                        $boxEndScalar = $hmm_spec->{$repeatc}->{end} + 1;
                        push( @{ $BOX{$boxnum} }, $repeatc );
                    }
                }
                if ( $#{ $BOX{$boxnum} } > 0 ) {
                    $boxStartHash{$boxnum}  = $boxstartScalar;
                    $boxEndHash{$boxnum}    = $boxEndScalar - 1;
                    $boxStrandHash{$boxnum} = $boxstrand;
                    $boxnum++;
                }
            }
            else {
                $boxEndScalar = $nodeA->{end} - 5;
                for my $repeatb (
                    sort { $hmm_spec->{$b}{start} <=> $hmm_spec->{$a}{start} }
                    keys %$hmm_spec
                  )
                {
                    if (
                        $hmm_spec->{$repeatb}->{type} =~ /boxB/
                        && (   $hmm_spec->{$repeatb}{start} >= $boxEndScalar
                            && $hmm_spec->{$repeatb}{start} <= $boxstartScalar )
                        && $hmm_spec->{$repeatb}->{strand} == $boxstrand
                      )
                    {
                        $boxEndScalar = $hmm_spec->{$repeatb}->{end} - 1;
                        push( @{ $BOX{$boxnum} }, $repeatb );
                    }
                }
                for my $repeatc (
                    sort { $hmm_spec->{$b}{start} <=> $hmm_spec->{$a}{start} }
                    keys %$hmm_spec
                  )
                {
                    if (
                        $hmm_spec->{$repeatc}->{type} =~ /boxC/
                        && (   $hmm_spec->{$repeatc}{start} >= $boxEndScalar
                            && $hmm_spec->{$repeatc}{start} <= $boxstartScalar )
                        && $hmm_spec->{$repeatc}->{strand} == $boxstrand
                      )
                    {
                        $boxEndScalar = $hmm_spec->{$repeatc}->{end} - 1;
                        push( @{ $BOX{$boxnum} }, $repeatc );
                    }
                }
                if ( $#{ $BOX{$boxnum} } > 0 ) {
                    $boxStartHash{$boxnum}  = $boxstartScalar;
                    $boxEndHash{$boxnum}    = $boxEndScalar + 1;
                    $boxStrandHash{$boxnum} = $boxstrand;
                    $boxnum++;
                }
            }
        }
    }
}

sub find_box_elemets_with_no_starting_boxA_module {
    my ( $hmm_spec ) = @_;

    $in_box = 0;

    # find BC BOX elements with no starting boxA module
    for my $repeatc ( sort keys %$hmm_spec ) {
        if ( $hmm_spec->{$repeatc}->{type} =~ /boxC/ ) {
            for $boxnum ( sort keys %BOX ) {
                if ( grep( /$repeatc/, @{ $BOX{$boxnum} } ) ) {
                    $in_box = 1;
                }
            }
            if ( $in_box == 0 ) {
                $boxEndScalar = $hmm_spec->{$repeatc}->{end};
                $boxstrand    = $hmm_spec->{$repeatc}->{strand};
                @{ $BOX{$boxnum} } = "$repeatc";
                if ( $boxstrand == -1 ) {
                    $boxstartScalar = $hmm_spec->{$repeatc}{start} + 1;
                    for my $repeatb (
                        sort {
                            $hmm_spec->{$a}{start} <=> $hmm_spec->{$b}{start}
                        }
                        keys %$hmm_spec
                      )
                    {
                        if (
                            $hmm_spec->{$repeatb}->{type} =~ /boxB/
                            && ( $hmm_spec->{$repeatb}->{end} <= $boxstartScalar
                                && $hmm_spec->{$repeatb}{start} >=
                                $boxstartScalar )
                            && $hmm_spec->{$repeatb}->{strand} == $boxstrand
                          )
                        {
                            $boxstartScalar = $hmm_spec->{$repeatb}{start} + 1;
                            push( @{ $BOX{$boxnum} }, $repeatb );
                        }
                    }
                    if ( $#{ $BOX{$boxnum} } > 0 ) {
                        $boxStartHash{$boxnum}  = $boxstartScalar - 1;
                        $boxEndHash{$boxnum}    = $boxEndScalar;
                        $boxStrandHash{$boxnum} = $boxstrand;
                        $boxnum++;
                    }
                }
                else {
                    $boxstartScalar = $hmm_spec->{$repeatc}{start} - 1;
                    for my $repeatb (
                        sort {
                            $hmm_spec->{$b}{start} <=> $hmm_spec->{$a}{start}
                        }
                        keys %$hmm_spec
                      )
                    {
                        if (
                            $hmm_spec->{$repeatb}->{type} =~ /boxB/
                            && ( $hmm_spec->{$repeatb}->{end} >= $boxstartScalar
                                && $hmm_spec->{$repeatb}{start} <=
                                $boxstartScalar )
                            && $hmm_spec->{$repeatb}->{strand} == $boxstrand
                          )
                        {
                            $boxstartScalar = $hmm_spec->{$repeatb}{start} - 1;
                            push( @{ $BOX{$boxnum} }, $repeatb );
                        }
                    }
                    if ( $#{ $BOX{$boxnum} } > 0 ) {
                        $boxStartHash{$boxnum}  = $boxstartScalar + 1;
                        $boxEndHash{$boxnum}    = $boxEndScalar;
                        $boxStrandHash{$boxnum} = $boxstrand;
                        $boxnum++;
                    }
                }
            }
            $in_box = 0;
        }
    }
}

# identify overlapping repeat elements, except BOX modules
sub identify_overlapping_repeat_elements_except_box_modules {
    my ( $hmm_spec, $overlaps ) = @_;

    for my $repeata ( sort keys %$hmm_spec ) {
        my $nodeA = $hmm_spec->{$repeata};
        for my $repeatb ( sort keys %$hmm_spec ) {
            if (
                (
                       $hmm_spec->{$repeatb}{start} <= $nodeA->{start}
                    && $hmm_spec->{$repeatb}->{end} >= $nodeA->{start}
                )
                || (   $hmm_spec->{$repeatb}{start} <= $nodeA->{end}
                    && $hmm_spec->{$repeatb}->{end} >= $nodeA->{end} )
              )
            {
                unless (
                    (
                           $nodeA->{type} =~ /box/
                        && $hmm_spec->{$repeatb}->{type} =~ /box/
                    )
                    || ( $repeata == $repeatb )
                  )
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
sub rescan_regions_around_overlapping_repeats_to_find_disruptions {
    my ( $hmm_spec, $overlaps, $sequence ) = @_;

    for my $repeata ( sort keys %$overlaps ) {

        boundaries( $repeata, $hmm_spec );
        $A_start = $startScalar;
        $A_end   = $endScalar;
        my $upper_bound = $endScalar + 251;
        my $lower_bound = $startScalar - 251;
        my $lower_seq   = $sequence->subseq( $lower_bound,   $startScalar - 1 );
        my $upper_seq   = $sequence->subseq( $endScalar + 1, $upper_bound );
        my $total_seq   = "$lower_seq" . "$upper_seq";
        print_fasta( "first.seq", $total_seq );
        system(
"$hmmls_script $options->{$hmm_spec->{$overlaps->{$repeata}}->{type}}->{hmm} first.seq > first.out"
        );
        hmm_results( "first.out" );
        $A_score     = $hitscore;
        $A_hitstrand = $hitstrand;
        $A_hitstart  = $hitstart;
        $A_hitend    = $hitend;
        boundaries( $overlaps->{$repeata}, $hmm_spec );
        $B_start     = $startScalar;
        $B_end       = $endScalar;
        $upper_bound = $endScalar + 251;
        $lower_bound = $startScalar - 251;
        $lower_seq   = $sequence->subseq( $lower_bound,   $startScalar - 1 );
        $upper_seq   = $sequence->subseq( $endScalar + 1, $upper_bound );
        $total_seq   = "$lower_seq" . "$upper_seq";
        print_fasta( "second.seq", $total_seq );
        system(
"$hmmls_script $options->{$hmm_spec->{$repeata}->{type}}->{hmm} second.seq > second.out"
        );
        hmm_results( "second.out" );
        $B_score     = $hitscore;
        $B_hitstrand = $hitstrand;
        $B_hitstart  = $hitstart;
        $B_hitend    = $hitend;
        my $A_inc =
          ( $A_score - $hmm_spec->{ $overlaps->{$repeata} }->{score} ) /
          $hmm_spec->{ $overlaps->{$repeata} }->{score};
        my $B_inc = ( $B_score - $hmm_spec->{$repeata}->{score} ) /
          $hmm_spec->{$repeata}->{score};

        if ( $A_inc >= 0 || $B_inc >= 0 ) {    # identify disruption events
            my $firstpoint;
            my $secondpoint;
            my $thirdpoint;
            my $fourthpoint;
            if ( $A_inc > $B_inc ) {
                $hmm_spec->{ $overlaps->{$repeata} }->{start} = "BROKEN";
                if ( $A_hitstrand == 1 ) {
                    $secondpoint = $A_start - 1;
                    $thirdpoint  = $A_end + 1;
                    $firstpoint  = $A_start - 252 + $A_hitstart;
                    $fourthpoint = $A_end - 251 + $A_hitend;
                    @{ $brokencoords{ $overlaps->{$repeata} } } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $hmm_spec->{$repeata}->{type} =~ /box/ ) {
                        $brokendisrupt{ $overlaps->{$repeata} } = "BOX";
                    }
                    else {
                        $brokendisrupt{ $overlaps->{$repeata} } =
                          $hmm_spec->{$repeata}->{type};
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
                    if ( $hmm_spec->{$repeata}->{type} =~ /box/ ) {
                        $brokendisrupt{ $overlaps->{$repeata} } = "BOX";
                    }
                    else {
                        $brokendisrupt{ $overlaps->{$repeata} } =
                          $hmm_spec->{$repeata}->{type};
                    }
                    $brokenscore{ $overlaps->{$repeata} } = $A_score;
                }
            }
            elsif ( $B_inc > $A_inc ) {
                $hmm_spec->{$repeata}{start} = "BROKEN";
                if ( $B_hitstrand == 1 ) {
                    $secondpoint = $B_start - 1;
                    $thirdpoint  = $B_end + 1;
                    $firstpoint  = $B_start - 252 + $B_hitstart;
                    $fourthpoint = $B_end - 251 + $B_hitend;
                    @{ $brokencoords{$repeata} } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $hmm_spec->{ $overlaps->{$repeata} }->{type} =~ /box/ )
                    {
                        $brokendisrupt{$repeata} = "BOX";
                    }
                    else {
                        $brokendisrupt{$repeata} =
                          $hmm_spec->{ $overlaps->{$repeata} }->{type};
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
                    if ( $hmm_spec->{ $overlaps->{$repeata} }->{type} =~ /box/ )
                    {
                        $brokendisrupt{$repeata} = "BOX";
                    }
                    else {
                        $brokendisrupt{$repeata} =
                          $hmm_spec->{ $overlaps->{$repeata} }->{type};
                    }
                    $brokenscore{$repeata} = $B_score;
                }
            }
        }

    }
}

sub save_results {
    my ( $hmm_spec, $brokencoords, $genome, $BOX ) = @_;

    my $file = "$genome.repeats.tab";
    open my $fh, ">", $file or die "Cannot open output file: $file\n";

    # Broken regions.
    for my $repeat ( sort keys %$hmm_spec ) {
        my $node = $hmm_spec->{$repeat};

        if ( $node->{start} eq "BROKEN" ) {
            my $coords = $brokencoords->{$repeat};
            my $type   = $node->{strand} == 1 ? "order" : "complement";

            print $fh <<~OUTPUT;
            FT   repeat_unit     $type($coords->[0]..$coords->[1],$coords->[2]..$coords->[3])
            FT                   /colour=2
            FT                   /label=$node->{type}
            FT                   /note=Detected using HMMER $hmmls_script; appears to have been disrupted through $brokendisrupt{$repeat} insertion
            FT                   /note=Initial match of score $node->{score} to model $node->{type}; realignment score of $brokenscore{$repeat}
            OUTPUT
        }
        else {
            my $type = $node->{strand} == 1 ? "" : "complement";

            print $fh <<~OUTPUT;
            FT   repeat_unit     $node->{start}..$node->{end}
            FT                   /colour=2
            FT                   /label=$node->{type}
            FT                   /note=Detected using HMMER $hmmls_script; match of score $node->{score} to model $node->{type}
            OUTPUT
        }
    }

    for my $boxnum ( sort keys %$BOX ) {
        my $type = $boxStrandHash{$boxnum} == 1 ? "" : "complement";
        print $fh <<~OUTPUT;
        FT   repeat_unit     $type($boxStartHash{$boxnum}..$boxEndHash{$boxnum})
        FT                   /colour=4
        FT                   /note=Composite BOX element
        FT                   /label=BOX
        OUTPUT
    }

    close $fh;
}

# subroutine for identifying repeat boundaries, esp BOX elements
sub boundaries {
    my ( $rep, $hmm_spec ) = @_;
    my $node = $hmm_spec->{$rep};

    if ( $node->{type} =~ /box/ ) {
        for $boxnum ( sort keys %BOX ) {
            for my $module ( @{ $BOX{$boxnum} } ) {
                if ( $module eq $rep ) {
                    my @unsorted =
                      ( $boxStartHash{$boxnum}, $boxEndHash{$boxnum} );
                    my @sorted = sort( @unsorted );
                    $startScalar = $sorted[0];
                    $endScalar   = $sorted[1];
                    $in_box      = 1;
                }
            }
        }
    }

    if ( not $in_box ) {
        my @unsorted = ( $node->{start}, $node->{end} );
        my @sorted   = sort( @unsorted );
        $startScalar = $sorted[0];
        $endScalar   = $sorted[1];
    }

    $in_box = 0;

    return ( $startScalar, $endScalar );
}

# subroutine for printing sequence in a suitable format for /home/sreeram/hmmer-1.8.4/hmmls
sub print_fasta {
    my $filename   = shift;
    my $dna_string = shift;

    open my $fh, ">", $filename;
    print $fh ">$filename\n";

    my $offset = 0;
    while ( $offset < 440 ) {
        my $line = substr( $dna_string, $offset, 60 );
        print $fh "$line\n";
        $offset += 60;
    }
    my $line = substr( $dna_string, $offset, ( 500 - $offset ) );
    print $fh "$line\n";

    close $fh;
}

# subroutine for picking the top hit from the /home/sreeram/hmmer-1.8.4/hmmls results
sub hmm_results {
    my $filename = shift;
    open HMM, $filename, or die "Cannot open HMM file\n";
    $hitscore = 0;
    for ( <HMM> ) {
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


