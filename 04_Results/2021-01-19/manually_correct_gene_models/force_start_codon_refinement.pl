#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Fasta_reader;
use Nuc_translator;
use Carp;
use Data::Dumper;
use List::Util qw (max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use PWM;
use File::Basename;


my $atg_pwm_pos = 20;
my $adj_dist = 30;
my $adj_pct = 15;

my $usage = <<__EOUSAGE__;

#####################################################################################
#
#  --transcripts <string>   transcripts.fasta file targeted by transdecoder
#
#  --gff3_file <string>     target gff3 file
#
#  optional:
#
#  --adj_dist <int>         distance allowed for start adjustment (default: $adj_dist)
# 
#  --adj_pct <int>          pecentage of orf length for examining start adjustment (default: $adj_pct)
#
#  --debug                  verbose
#
#######################################################################################



__EOUSAGE__

    ;


my $transcripts_file;
my $gff3_file;
my $DEBUG = 0;
my $workdir;

&GetOptions('transcripts=s' => \$transcripts_file,
            'gff3_file=s' => \$gff3_file,
            'adj_dist=i' => \$adj_dist,
            'atg_pct=i' => \$adj_pct,
            'debug' => \$DEBUG,
    );


unless ($transcripts_file && $gff3_file) {
    die $usage;
}

if ($adj_pct > 30 || $adj_pct < 0) {
    die "Error, --adj_pct is out of range...  must be between 0 and 30 ";
}

main: {

    print STDERR "-reading transcripts: $transcripts_file\n" if $DEBUG;
    my $fasta_reader = new Fasta_reader($transcripts_file);

    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    print STDERR "-parsing orf annotations: $gff3_file\n" if $DEBUG;
    my $gene_obj_indexer_href = {};

    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

    my $num_starts_revised = 0;
    
    foreach my $transcript_acc (sort keys %$asmbl_id_to_gene_list_href) {

        print STDERR "-processing: $transcript_acc\n" if $DEBUG;
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$transcript_acc}};
        

        my $transcript_seq = uc $seqs{$transcript_acc};

        foreach my $gene_id (@gene_ids) {

            my $gene_obj = $gene_obj_indexer_href->{$gene_id};
            
            my $revised_start_flag = &refine_start_codon_position($transcript_acc, $gene_id,
                                                                  $gene_obj, $transcript_seq);

            if ($revised_start_flag) {
                $num_starts_revised++;

                ## update naming convention based on now having an updated start codon.
                if ($gene_obj->{com_name} =~ /internal/) {
                    $gene_obj->{com_name} =~ s/internal/3prime_partial/;
                }
                elsif ($gene_obj->{com_name} =~ /5prime_partial/) {
                    $gene_obj->{com_name} =~ s/5prime_partial/complete/;
                }
            } 
            print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
            

        }
    }
    
    print STDERR "-number of revised start positions: $num_starts_revised\n";

    exit(0);
    
}


####
sub refine_start_codon_position {
    my ($transcript_acc, $gene_id,
        $gene_obj, $transcript_seq) = @_;
    
    my $revised_start_flag = 0;
    
    my $orient = $gene_obj->get_orientation();
    
    my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_model_span();
    my $orf_len = $rend - $lend + 1;
    if ($orf_len % 3 != 0) {
        die "Error, $orf_len is not mod 3 " . $gene_obj->toString();
    }
    
    my $orig_start_pos = $lend;
    
    my $start_pos = $lend;
    if ($orient eq '-') {
        $transcript_seq = &reverse_complement($transcript_seq);
        $start_pos = length($transcript_seq) - $rend + 1;
    }


    # only work on 5' partials
    my $start_index = $start_pos - 1; # zero based

    if (substr($transcript_seq, $start_index, 3) eq "ATG") {
        return(0);
    }
    
    my @alt_starts;


    
    my $max_search_pos = max($start_index + $adj_dist, $start_index + int($adj_pct * $orf_len / 100));
    
    my $best_alt_start = undef;

    while ($transcript_seq =~ /(ATG)/g) {
        my $pos = $-[0];
        if ($pos > $max_search_pos) { last; } # too far
        if ($pos > $start_index 
            && 
            ($pos - $start_index) % 3  == 0) { # in frame start

            $best_alt_start = $pos; 
	    last; # Found first alt-start, end loop
        }
    }

    unless ($best_alt_start) {
        return($revised_start_flag);
    }
    
    
    if ($best_alt_start) {
        $best_alt_start++; # make 1-based coord
        my $new_start = $best_alt_start;
        if ($orient eq '-') {
            $new_start = length($transcript_seq) - $best_alt_start + 1;
        }
        
        my ($exon_obj) = $gene_obj->get_exons();
        my $cds_obj = $exon_obj->get_CDS_obj();
        $cds_obj->{end5} = $new_start;
        
        $gene_obj->refine_gene_object();

        $revised_start_flag = 1;
        
        print STDERR "# refined start codon: $orig_start_pos -> $new_start\n" if $DEBUG;
        print "# refined start codon: $orig_start_pos -> $new_start\n";
    }

    return($revised_start_flag);
    
}

