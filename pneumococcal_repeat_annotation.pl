#! /software/bin/perl -w

use v5.26;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqI;

my %startHash;
my %endHash;
my %type;
my %score;
my %strand;
my $r = 0;
my %BOX;
my $boxstartScalar;
my $boxEndScalar;
my $boxnum = 0;
my $boxstrand;
my %boxStartHash;
my %boxstrand;
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

my @repeats =
  ( "boxA", "boxB", "boxC", "RUP", "SPRITE" );    # repeats being analysed
my %cutoffs = (    # cutoff scores for HMM analysis
    "boxA"   => 30,
    "boxB"   => 14,
    "boxC"   => 28,
    "RUP"    => 47,
    "SPRITE" => 66,
);
my %hmms = (
    "boxA"   => "./Repeat_HMMs/boxA.hmm",     # paths of HMM for analysis
    "boxB"   => "./Repeat_HMMs/boxB.hmm",
    "boxC"   => "./Repeat_HMMs/boxC.hmm",
    "RUP"    => "./Repeat_HMMs/RUP.hmm",
    "SPRITE" => "./Repeat_HMMs/SPRITE.hmm",
);
my %input_files = (
    "boxA"   => "boxA.rep",                   # output files from HMM analysis
    "boxB"   => "boxB.rep",
    "boxC"   => "boxC.rep",
    "RUP"    => "RUP.rep",
    "SPRITE" => "SPRITE.rep",
);

foreach my $seq ( @ARGV ) {
    foreach my $repeat ( @repeats ) {
        print STDERR "Searching sequence $seq for repeat $repeat...";
        system(
"/home/sreeram/hmmer-1.8.4/hmmls -c -t $cutoffs{$repeat} $hmms{$repeat} $seq > $seq.$input_files{$repeat}"
        );
        print STDERR "done\n";
    }
}

foreach my $genome ( @ARGV ) {
    foreach my $repeat ( sort keys %input_files )
    {    # parse information from HMM output
        open IN, "$genome.$input_files{$repeat}"
          or die print "Cannot open input file\n";
        foreach ( <IN> ) {
            chomp;
            if ( substr( $_, 0, 2 ) =~ /\d\d/ ) {
                my @array      = split( /\s+/, $_ );
                my @startArray = split( /:/,   $array[2] );
                my @finish     = split( /:/,   $array[3] );
                if ( $array[0] > $cutoffs{$repeat} ) {
                    $startHash{$r} = $startArray[1];
                    $endHash{$r}   = $finish[1];
                    $type{$r}      = $repeat;
                    $score{$r}     = $array[0];
                    if ( $endHash{$r} > $startHash{$r} ) {
                        $strand{$r} = 1;
                    }
                    else {
                        $strand{$r} = -1;
                    }
                }
                $r++;
            }
        }
    }
    print STDERR "Identified repeat units in sequence $genome\n";
    foreach my $repeata ( sort keys %startHash )
    {    # identify composite BOX elements from box modules
        if ( $type{$repeata} =~ /boxA/ ) {
            $boxstrand      = $strand{$repeata};
            $boxstartScalar = $startHash{$repeata};
            @{ $BOX{$boxnum} } = "$repeata";
            if ( $boxstrand == 1 ) {
                $boxEndScalar = $endHash{$repeata} +
                  5;    # five base leeway for finding the adjacent boxB element
                foreach my $repeatb (
                    sort { $startHash{$a} <=> $startHash{$b} }
                    keys %startHash
                  )
                {
                    if (
                        $type{$repeatb} =~ /boxB/
                        && (   $startHash{$repeatb} <= $boxEndScalar
                            && $startHash{$repeatb} >= $boxstartScalar )
                        && $strand{$repeatb} == $boxstrand
                      )
                    {
                        $boxEndScalar = $endHash{$repeatb} + 1;
                        push( @{ $BOX{$boxnum} }, $repeatb );
                    }
                }
                foreach my $repeatc (
                    sort { $startHash{$a} <=> $startHash{$b} }
                    keys %startHash
                  )
                {
                    if (
                        $type{$repeatc} =~ /boxC/
                        && (   $startHash{$repeatc} <= $boxEndScalar
                            && $startHash{$repeatc} >= $boxstartScalar )
                        && $strand{$repeatc} == $boxstrand
                      )
                    {
                        $boxEndScalar = $endHash{$repeatc} + 1;
                        push( @{ $BOX{$boxnum} }, $repeatc );
                    }
                }
                if ( $#{ $BOX{$boxnum} } > 0 ) {
                    $boxStartHash{$boxnum} = $boxstartScalar;
                    $boxEndHash{$boxnum}   = $boxEndScalar - 1;
                    $boxstrand{$boxnum}    = $boxstrand;
                    $boxnum++;
                }
            }
            else {
                $boxEndScalar = $endHash{$repeata} - 5;
                foreach my $repeatb (
                    sort { $startHash{$b} <=> $startHash{$a} }
                    keys %startHash
                  )
                {
                    if (
                        $type{$repeatb} =~ /boxB/
                        && (   $startHash{$repeatb} >= $boxEndScalar
                            && $startHash{$repeatb} <= $boxstartScalar )
                        && $strand{$repeatb} == $boxstrand
                      )
                    {
                        $boxEndScalar = $endHash{$repeatb} - 1;
                        push( @{ $BOX{$boxnum} }, $repeatb );
                    }
                }
                foreach my $repeatc (
                    sort { $startHash{$b} <=> $startHash{$a} }
                    keys %startHash
                  )
                {
                    if (
                        $type{$repeatc} =~ /boxC/
                        && (   $startHash{$repeatc} >= $boxEndScalar
                            && $startHash{$repeatc} <= $boxstartScalar )
                        && $strand{$repeatc} == $boxstrand
                      )
                    {
                        $boxEndScalar = $endHash{$repeatc} - 1;
                        push( @{ $BOX{$boxnum} }, $repeatc );
                    }
                }
                if ( $#{ $BOX{$boxnum} } > 0 ) {
                    $boxStartHash{$boxnum} = $boxstartScalar;
                    $boxEndHash{$boxnum}   = $boxEndScalar + 1;
                    $boxstrand{$boxnum}    = $boxstrand;
                    $boxnum++;
                }
            }
        }
    }
    $in_box = 0;
    foreach my $repeatc ( sort keys %startHash )
    {    # find BC BOX elements with no starting boxA module
        if ( $type{$repeatc} =~ /boxC/ ) {
            foreach $boxnum ( sort keys %BOX ) {
                if ( grep( /$repeatc/, @{ $BOX{$boxnum} } ) ) {
                    $in_box = 1;
                }
            }
            if ( $in_box == 0 ) {
                $boxEndScalar = $endHash{$repeatc};
                $boxstrand    = $strand{$repeatc};
                @{ $BOX{$boxnum} } = "$repeatc";
                if ( $boxstrand == -1 ) {
                    $boxstartScalar = $startHash{$repeatc} + 1;
                    foreach my $repeatb (
                        sort { $startHash{$a} <=> $startHash{$b} }
                        keys %startHash
                      )
                    {
                        if (
                            $type{$repeatb} =~ /boxB/
                            && (   $endHash{$repeatb} <= $boxstartScalar
                                && $startHash{$repeatb} >= $boxstartScalar )
                            && $strand{$repeatb} == $boxstrand
                          )
                        {
                            $boxstartScalar = $startHash{$repeatb} + 1;
                            push( @{ $BOX{$boxnum} }, $repeatb );
                        }
                    }
                    if ( $#{ $BOX{$boxnum} } > 0 ) {
                        $boxStartHash{$boxnum} = $boxstartScalar - 1;
                        $boxEndHash{$boxnum}   = $boxEndScalar;
                        $boxstrand{$boxnum}    = $boxstrand;
                        $boxnum++;
                    }
                }
                else {
                    $boxstartScalar = $startHash{$repeatc} - 1;
                    foreach my $repeatb (
                        sort { $startHash{$b} <=> $startHash{$a} }
                        keys %startHash
                      )
                    {
                        if (
                            $type{$repeatb} =~ /boxB/
                            && (   $endHash{$repeatb} >= $boxstartScalar
                                && $startHash{$repeatb} <= $boxstartScalar )
                            && $strand{$repeatb} == $boxstrand
                          )
                        {
                            $boxstartScalar = $startHash{$repeatb} - 1;
                            push( @{ $BOX{$boxnum} }, $repeatb );
                        }
                    }
                    if ( $#{ $BOX{$boxnum} } > 0 ) {
                        $boxStartHash{$boxnum} = $boxstartScalar + 1;
                        $boxEndHash{$boxnum}   = $boxEndScalar;
                        $boxstrand{$boxnum}    = $boxstrand;
                        $boxnum++;
                    }
                }
            }
            $in_box = 0;
        }
    }
    print STDERR "Identified composite BOX elements\n";
    delete( $BOX{$boxnum} );

    my %overlaps;

    if ( defined $startScalar ) {
        foreach my $repeata ( sort keys %startHash )
        {    # identify overlapping repeat elements, except BOX modules
            foreach my $repeatb ( sort keys %startHash ) {
                if (
                    (
                           $startHash{$repeatb} <= $startHash{$repeata}
                        && $endHash{$repeatb} >= $startHash{$repeata}
                    )
                    || (   $startHash{$repeatb} <= $endHash{$repeata}
                        && $endHash{$repeatb} >= $endHash{$repeata} )
                  )
                {
                    unless (
                        (
                               $type{$repeata} =~ /box/
                            && $type{$repeatb} =~ /box/
                        )
                        || ( $repeata == $repeatb )
                      )
                    {
                        my @unsorted = ( $repeata, $repeatb );
                        my @sorted   = sort( @unsorted );
                        $overlaps{ $sorted[0] } = $sorted[1];
                    }
                }
            }
        }
    }

    print STDERR "Realigning overlapping elements\n";

    my $fasta    = Bio::SeqIO->new( -file => $genome, -format => "fasta" );
    my $sequence = $fasta->next_seq;

    foreach my $repeata ( sort keys %overlaps )
    {    # rescan regions around overlapping repeats to look for disruptions
        boundaries( $repeata );
        $A_start = $startScalar;
        $A_end   = $endScalar;
        my $upper_bound = $endScalar + 251;
        my $lower_bound = $startScalar - 251;
        my $lower_seq   = $sequence->subseq( $lower_bound,   $startScalar - 1 );
        my $upper_seq   = $sequence->subseq( $endScalar + 1, $upper_bound );
        my $total_seq   = "$lower_seq" . "$upper_seq";
        print_fasta( "first.seq", $total_seq );
        system(
"/home/sreeram/hmmer-1.8.4/hmmls -c $hmms{$type{$overlaps{$repeata}}} first.seq > first.out"
        );
        hmm_results( "first.out" );
        $A_score     = $hitscore;
        $A_hitstrand = $hitstrand;
        $A_hitstart  = $hitstart;
        $A_hitend    = $hitend;
        boundaries( $overlaps{$repeata} );
        $B_start     = $startScalar;
        $B_end       = $endScalar;
        $upper_bound = $endScalar + 251;
        $lower_bound = $startScalar - 251;
        $lower_seq   = $sequence->subseq( $lower_bound,   $startScalar - 1 );
        $upper_seq   = $sequence->subseq( $endScalar + 1, $upper_bound );
        $total_seq   = "$lower_seq" . "$upper_seq";
        print_fasta( "second.seq", $total_seq );
        system(
"/home/sreeram/hmmer-1.8.4/hmmls -c $hmms{$type{$repeata}} second.seq > second.out"
        );
        hmm_results( "second.out" );
        $B_score     = $hitscore;
        $B_hitstrand = $hitstrand;
        $B_hitstart  = $hitstart;
        $B_hitend    = $hitend;
        my $A_inc = ( $A_score - $score{ $overlaps{$repeata} } ) /
          $score{ $overlaps{$repeata} };
        my $B_inc = ( $B_score - $score{$repeata} ) / $score{$repeata};

        if ( $A_inc >= 0 || $B_inc >= 0 ) {    # identify disruption events
            my $firstpoint;
            my $secondpoint;
            my $thirdpoint;
            my $fourthpoint;
            if ( $A_inc > $B_inc ) {
                $startHash{ $overlaps{$repeata} } = "BROKEN";
                if ( $A_hitstrand == 1 ) {
                    $secondpoint = $A_start - 1;
                    $thirdpoint  = $A_end + 1;
                    $firstpoint  = $A_start - 252 + $A_hitstart;
                    $fourthpoint = $A_end - 251 + $A_hitend;
                    @{ $brokencoords{ $overlaps{$repeata} } } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{$repeata} =~ /box/ ) {
                        $brokendisrupt{ $overlaps{$repeata} } = "BOX";
                    }
                    else {
                        $brokendisrupt{ $overlaps{$repeata} } = $type{$repeata};
                    }
                    $brokenscore{ $overlaps{$repeata} } = $A_score;
                }
                else {
                    $secondpoint = $A_start - 1;
                    $thirdpoint  = $A_end + 1;
                    $firstpoint  = $A_start - 252 + $A_hitend;
                    $fourthpoint = $A_end - 251 + $A_hitstart;
                    @{ $brokencoords{ $overlaps{$repeata} } } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{$repeata} =~ /box/ ) {
                        $brokendisrupt{ $overlaps{$repeata} } = "BOX";
                    }
                    else {
                        $brokendisrupt{ $overlaps{$repeata} } = $type{$repeata};
                    }
                    $brokenscore{ $overlaps{$repeata} } = $A_score;
                }
            }
            elsif ( $B_inc > $A_inc ) {
                $startHash{$repeata} = "BROKEN";
                if ( $B_hitstrand == 1 ) {
                    $secondpoint = $B_start - 1;
                    $thirdpoint  = $B_end + 1;
                    $firstpoint  = $B_start - 252 + $B_hitstart;
                    $fourthpoint = $B_end - 251 + $B_hitend;
                    @{ $brokencoords{$repeata} } =
                      ( $firstpoint, $secondpoint, $thirdpoint, $fourthpoint );
                    if ( $type{ $overlaps{$repeata} } =~ /box/ ) {
                        $brokendisrupt{$repeata} = "BOX";
                    }
                    else {
                        $brokendisrupt{$repeata} = $type{ $overlaps{$repeata} };
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
                    if ( $type{ $overlaps{$repeata} } =~ /box/ ) {
                        $brokendisrupt{$repeata} = "BOX";
                    }
                    else {
                        $brokendisrupt{$repeata} = $type{ $overlaps{$repeata} };
                    }
                    $brokenscore{$repeata} = $B_score;
                }
            }
        }

    }

    print STDERR "Disrupted repeat elements identified\n";

    open OUT, "> $genome.repeats.tab"
      or die print STDERR "Cannot open output file\n";
    print STDERR "Printing output files\n";

    foreach my $repeat ( sort keys %startHash ) {
        if ( $startHash{$repeat} eq "BROKEN" ) {
            if ( $strand{$repeat} == 1 ) {
                print OUT
"FT   repeat_unit     order(${$brokencoords{$repeat}}[0]..${$brokencoords{$repeat}}[1],${$brokencoords{$repeat}}[2]..${$brokencoords{$repeat}}[3])\n";
            }
            elsif ( $strand{$repeat} == -1 ) {
                print OUT
"FT   repeat_unit     complement(${$brokencoords{$repeat}}[0]..${$brokencoords{$repeat}}[1],${$brokencoords{$repeat}}[2]..${$brokencoords{$repeat}}[3])\n";
            }
            print OUT "FT                   /colour=2\n";
            print OUT "FT                   /label=$type{$repeat}\n";
            print OUT
"FT                   /note=Detected using HMMER /home/sreeram/hmmer-1.8.4/hmmls; appears to have been disrupted through $brokendisrupt{$repeat} insertion\n";
            print OUT
"FT                   /note=Initial match of score $score{$repeat} to model $type{$repeat}; realignment score of $brokenscore{$repeat}\n";
        }
        else {
            if ( $strand{$repeat} == 1 ) {
                print OUT
"FT   repeat_unit     $startHash{$repeat}..$endHash{$repeat}\n";
            }
            elsif ( $strand{$repeat} == -1 ) {
                print OUT
"FT   repeat_unit     complement($startHash{$repeat}..$endHash{$repeat})\n";
            }
            print OUT "FT                   /colour=2\n";
            print OUT "FT                   /label=$type{$repeat}\n";
            print OUT
"FT                   /note=Detected using HMMER /home/sreeram/hmmer-1.8.4/hmmls; match of score $score{$repeat} to model $type{$repeat}\n";
        }
    }
    foreach my $boxnum ( sort keys %BOX ) {
        if ( $boxstrand{$boxnum} == 1 ) {
            print OUT
"FT   repeat_unit     $boxStartHash{$boxnum}..$boxEndHash{$boxnum}\n";
        }
        else {
            print OUT
"FT   repeat_unit     complement($boxStartHash{$boxnum}..$boxEndHash{$boxnum})\n";
        }
        print OUT "FT                   /colour=4\n";
        print OUT "FT                   /note=Composite BOX element\n";
        print OUT "FT                   /label=BOX\n";
    }

    undef( %startHash );
    undef( %endHash );
    undef( %type );
    undef( %score );
    undef( %strand );
    undef( %BOX );
    undef( %boxStartHash );
    undef( %boxstrand );
    undef( %boxEndHash );
    undef( %brokencoords );
    undef( %brokenscore );
    undef( %brokendisrupt );

}

close OUT;
system( "rm first.out first.seq second.out second.seq *.rep;" );
print STDERR "Done\n";

sub boundaries
{    # subroutine for identifying repeat boundaries, esp BOX elements
    my $rep = shift;
    if ( $type{$rep} =~ /box/ ) {
        foreach $boxnum ( sort keys %BOX ) {
            foreach my $module ( @{ $BOX{$boxnum} } ) {
                if ( $module == $rep ) {
                    my @unsorted =
                      ( "$boxStartHash{$boxnum}", "$boxEndHash{$boxnum}" );
                    my @sorted = sort( @unsorted );
                    $startScalar = $sorted[0];
                    $endScalar   = $sorted[1];
                    $in_box      = 1;
                }
            }
        }
    }
    if ( $in_box == 0 ) {
        my @unsorted = ( "$startHash{$rep}", "$endHash{$rep}" );
        my @sorted   = sort( @unsorted );
        $startScalar = $sorted[0];
        $endScalar   = $sorted[1];
    }
    $in_box = 0;
    return ( $startScalar, $endScalar );
}

sub print_fasta
{ # subroutine for printing sequence in a suitable format for /home/sreeram/hmmer-1.8.4/hmmls
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

sub hmm_results
{ # subroutine for picking the top hit from the /home/sreeram/hmmer-1.8.4/hmmls results
    my $filename = shift;
    open HMM, $filename, or die print STDERR "Cannot open HMM file\n";
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
