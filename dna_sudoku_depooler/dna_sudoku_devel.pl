#!/usr/bin/perl
use 5.10.0;
use strict;
use warnings;

use File::Basename;
use File::Spec;

#---OPTIONS---
our $path_sep = File::Spec->catfile('', '');

our $help;
our $threads_number = 2;
our $config_file; # = "test.txt";                                               #debugging
our $gatk_location = catfile(catdir(dirname(__FILE__), "opt"), "GenomeAnalysisTK.jar"); 
our $bt2script_location = catfile(catdir(dirname(__FILE__), "opt"), "bt2alignment.sh");
our $depooler_location = catfile(catdir(dirname(__FILE__), "opt"), "SudokuDePooler-0.1.1-SNAPSHOT-jar-with-dependencies.jar"); 
our ($output_folder, $prefix, $opt_file) = ('run', 'run');
our ($scheme, $fq_file, $dmp_folder, $algn_folder, $bam_file, $snp_folder, $vcf_file, $depool_folder);
our ($demult_seq, @trim5_seq, @trim3_seq);
our ($bt2index, $reference, $regions);
our ($trimert, $trimqual, $trimlen) = (0.2, 25, 40);
our ($organism_ploidy, $poolsize);
our $monitor;

our %OPT = (
    HELP1           => 'help', 
    HELP2           => 'h',
    THREADS         => 'n',
    CFG_FILE        => 'cfg',
    GATK            => 'gatk',
    SCHEME          => 'scheme',
    DEMULT_FOLDER   => 'dmpf',
    FQ_FILE         => 'fq',
    BAM_FILE        => 'bam',
    VCF_FILE        => 'vcf',
    PREFIX          => 'prefix',
    OUTPUT          => 'output',
    DEMULT          => 'demult',
    TR_SEQ5         => 'trseq5',
    TR_SEQ3         => 'trseq3',
    TR_ERROR        => 'trer',
    TR_QUAL         => 'trqual',
    READS_LEN       => 'readlen',
    BT2_INDEX       => 'bt2index',
    REFERENCE       => 'reference',
    REGIONS         => 'regions',
    ORGANISM_PLOIDY => 'org_ploidy',
    POOL_SIZE       => 'pool_size',
    GRAFIC_MONITOR  => 'gm',
    );

#---LAUNCH---
&parse_options(@ARGV);
&print_help()   if ($help);
&parse_options(&read_config_file()) and &parse_options(@ARGV)   if ($config_file);
&create_output_folder();
$opt_file = catfile($output_folder, ($prefix ? "$prefix\_" : '') . "opt.txt");

if      ($vcf_file)     {&do_depooling();}
elsif   ($bam_file)     {&do_snpcalling();}
elsif   ($dmp_folder)   {&do_mapping();}
else                    {&do_demultiplexing();}

#---PROCEDURES---
sub do_demultiplexing {
    die "No fasta source file"              if not defined($fq_file) or not -e $fq_file;
    die "No demultiplexing sequences file"  if not defined($demult_seq) or not -e $demult_seq;
    $dmp_folder = &create_folder("dmp-reads");
    
    say "***Demultiplexing***";
    
    my $exit_code = system "cutadapt -g file:$demult_seq -o ${dmp_folder}${path_sep}\{name}.fq $fq_file \\
    >> $dmp_folder$path_sep\\dmp-report.txt";
    
    &do_reads_trimming() if not $exit_code;
}

sub do_reads_trimming {
    my @pools_list = &get_pools_list();

    say "***Reads Trimming***";

    foreach (@pools_list) {
        my $poolId = $_;
        my $raw_fq = catfile($dmp_folder, $_ . '.fq');
        die "No file with reads $raw_fq"    if not -e $raw_fq;
        my $tr_fq = catfile($dmp_folder, $_ . '_tr.fq');
        
        my $command;
        foreach (@trim5_seq) {
            my $source = $command ? '-' : $raw_fq;
            my $part = fileparse($_, qr/\.[^.]*/);
            my $report_file = catfile($dmp_folder, $poolId . '_' . $part .'_tr-rep.txt');
            $command .= "cutadapt -g file:$_ -e $trimert --discard-untrimmed $source 2> $report_file | ";
        }
        foreach (@trim3_seq) {
            my $source = $command ? '-' : $raw_fq;
            my $part = fileparse($_, qr/\.[^.]*/);
            my $report_file = catfile($dmp_folder, $poolId . '_' . $part .'_tr-rep.txt');
            $command .= "cutadapt -a file:$_ -e $trimert --discard-untrimmed $source 2> $report_file | ";
        }
        {
            my $source = $command ? '-' : $raw_fq;
            my $report_file = catfile($dmp_folder, $poolId . '_filt-rep.txt');
            $command .= "cutadapt -q $trimqual,$trimqual -m $trimlen -o $tr_fq $source > $report_file";
        }
        
        say $poolId;
        my $exit_code = system $command;
        die if $exit_code;
    }

    &write_demultiplexing_options() and &do_mapping();
}

sub do_mapping {
    die "No demultiplexing sequences file"  if not defined($demult_seq) or not -e $demult_seq;
    die "No bowtie2 index"                  if not defined($bt2index) or not -e "$bt2index.1.bt2";
    my @pools_list = &get_pools_list();
    $algn_folder = &create_folder("alignments");
    $bam_file = catfile($algn_folder, ($prefix ? "${prefix}_" : '') . 'merged.bam');
    my $merged_command = "samtools merge -f $bam_file ";

    say "***Reads Mapping***";

    foreach (@pools_list) {
        my $tr_fq = catfile($dmp_folder, $_ . '_tr.fq');
        my $poolbam = catfile($algn_folder, ($prefix ? "${prefix}_" : '') . $_);
        say $_;
        my $exit_code = system "$bt2script_location $poolbam $bt2index $tr_fq $_ $threads_number";
        die if $exit_code;
        $merged_command .= "${poolbam}.bam ";
    }

    say "***Alignments Merging***";

    my $exit_code = system $merged_command;
    die if $exit_code;
    $exit_code = system "samtools index $bam_file";

    &write_mapping_options() and &do_snpcalling() if not $exit_code;
}

sub do_snpcalling {
    die "No reference"          if not defined($reference) or not -e "$reference";
    die "No regions"            if not defined($regions) or not -e "$regions";
    die "No organism ploidy"    if not defined($organism_ploidy);
    die "No pool size"          if not defined($poolsize);
    $snp_folder = &create_folder("snp");
    $vcf_file = catfile($snp_folder, ($prefix ? "${prefix}_" : '') . 'snp.vcf');
    my $snp_calling_report = catfile($snp_folder, ($prefix ? "${prefix}_" : '') . 'report.txt');

    say "***SNP-calling***";
    my $poolploidy = $organism_ploidy * $poolsize;
    my $command = "java -jar $gatk_location -T HaplotypeCaller -gt_mode DISCOVERY "
    . "-R $reference --activeRegionIn $regions -I $bam_file --sample_ploidy $poolploidy "
    . "-o $vcf_file -nct $threads_number 2> $snp_calling_report";
    my $exit_code = system $command;

    &write_snpcalling_options() and &do_depooling() if not $exit_code;
}

sub do_depooling {
    die "No pooling scheme"     if not defined($scheme) or not -e $scheme;
    die "No vcf file"           if not defined($vcf_file) or not -e $vcf_file;
    $depool_folder = &create_folder("depooling");
    
    say "***De-pooling***";
    my $command = "java -jar $depooler_location -sch $scheme -vcf $vcf_file ";
    $command .= $monitor
        ? "-m" 
        : "> " . catfile($depool_folder, ($prefix ? "${prefix}_" : '') . 'depooling.txt'); 
    my $exit_code = system $command;

    &write_depooling_options() and say "Finished successfully" if not $exit_code;
}

#---OPTIONS FUNCTIONS---
sub parse_options {
    use Getopt::Long qw(GetOptionsFromArray);
    
    my $opt_parsing = GetOptionsFromArray(\@_,
        "$OPT{HELP1}|$OPT{HELP2}"   => \$help,
        "$OPT{THREADS}=i"           => \$threads_number,
        "$OPT{CFG_FILE}=s"          => \$config_file,
        "$OPT{GATK}=s"              => \$gatk_location,
        "$OPT{SCHEME}=s"            => \$scheme,
        "$OPT{FQ_FILE}=s"           => \$fq_file,
        "$OPT{DEMULT_FOLDER}=s"     => \$dmp_folder,
        "$OPT{BAM_FILE}=s"          => \$bam_file,
        "$OPT{VCF_FILE}=s"          => \$vcf_file,
        "$OPT{PREFIX}=s"            => \$prefix,
        "$OPT{OUTPUT}=s"            => \$output_folder,
        "$OPT{DEMULT}=s"            => \$demult_seq,
        "$OPT{TR_SEQ5}=s"           => \@trim5_seq,
        "$OPT{TR_SEQ3}=s"           => \@trim3_seq,
        "$OPT{TR_ERROR}=f"          => \$trimert,
        "$OPT{TR_QUAL}=i"           => \$trimqual, 
        "$OPT{READS_LEN}=i"         => \$trimlen,
        "$OPT{BT2_INDEX}=s"         => \$bt2index,
        "$OPT{REFERENCE}=s"         => \$reference, 
        "$OPT{REGIONS}=s"           => \$regions,
        "$OPT{ORGANISM_PLOIDY}=i"   => \$organism_ploidy,
        "$OPT{POOL_SIZE}=i"         => \$poolsize,
        "$OPT{GRAFIC_MONITOR}"      => \$monitor,
    );
    die "Unknown option"  if (not $opt_parsing);
    @trim5_seq = split(/,/,join(',',@trim5_seq));
    @trim3_seq = split(/,/,join(',',@trim3_seq));
}

sub write_demultiplexing_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{FQ_FILE},          $fq_file                if $fq_file;
    printf OPTF "-%s\t%s\n",    $OPT{TR_SEQ5},          join(',',@trim5_seq)    if @trim5_seq;
    printf OPTF "-%s\t%s\n",    $OPT{TR_SEQ3},          join(',',@trim3_seq)    if @trim3_seq;
    printf OPTF "-%s\t%s\n",    $OPT{TR_ERROR},         $trimert;
    printf OPTF "-%s\t%s\n",    $OPT{TR_QUAL},          $trimqual;
    printf OPTF "-%s\t%s\n",    $OPT{READS_LEN},        $trimlen;
    close OPTF;
}

sub write_mapping_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{DEMULT},           $demult_seq;
    printf OPTF "-%s\t%s\n",    $OPT{DEMULT_FOLDER},    $dmp_folder;
    printf OPTF "-%s\t%s\n",    $OPT{BT2_INDEX},        $bt2index;
    close OPTF;
}

sub write_snpcalling_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{BAM_FILE},         $bam_file;
    printf OPTF "-%s\t%s\n",    $OPT{REFERENCE},        $reference;
    printf OPTF "-%s\t%s\n",    $OPT{REGIONS},          $regions;
    printf OPTF "-%s\t%s\n",    $OPT{ORGANISM_PLOIDY},  $organism_ploidy;
    printf OPTF "-%s\t%s\n",    $OPT{POOL_SIZE},        $poolsize;
    printf OPTF "-%s\t%s\n",    $OPT{GATK},             $gatk_location;
    close OPTF;
}

sub write_depooling_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{VCF_FILE},         $vcf_file;
    printf OPTF "-%s\t%s\n",    $OPT{SCHEME},           $scheme;
    printf OPTF "-%s\n",        $OPT{GRAFIC_MONITOR},   $monitor    if $monitor;
    printf OPTF "-%s\t%s\n",    $OPT{PREFIX},           $prefix     if $prefix;
    close OPTF;
}

#---FUNCTIONS---
sub print_help {
    say "i'll help you";
    die;
}

sub read_config_file {
    say "Configuration file is '$config_file'";
    my @array;
    if (open CONFIG, "<", $config_file) {
        while (<CONFIG>) {
            my @tokens = split( /\s+/, $_ );
            if (@tokens) {
                $tokens[0] =~ s/(^\w)/-$1/;
#                say $tokens[0];
                push (@array, @tokens);
            }
        }
        close CONFIG;
        @array;
    } else {
        die "Cannot read contig file .";
    };
}

sub get_pools_list {
    my @pools_list;
    if (open FILE, "<", $demult_seq) {
        while (<FILE>) {
            if (s/^>//) {
                push @pools_list, ((split(/\s+/, $_))[0]);
            }
        }
        close FILE;
        @pools_list;
    } else {
        die "Cannot read pool identificator file .";
    }
}

sub create_output_folder {
    my $folder = canonpath($output_folder);
    my $index = 0;
    while(-d $folder) {
        $folder = canonpath($output_folder) . '_' . ++$index;
    }
    die "Cannot create folder $folder"  if not make_path ($folder);
    $output_folder = $folder;
}
    
sub create_folder {
    use File::Path qw(make_path);
    use File::Spec::Functions;
    
    my $folder = catdir(canonpath($output_folder), $prefix ? "$prefix\_$_[0]" : $_[0]);
    if (not -d $folder) {
        my $s = make_path ($folder);
        die "Cannot create folder $folder"  if not $s;
    } else {
#        my $st = (stat($folder))[2] & 060;                                      #debug
#        say $st;                                                                #debug
        die "Cannot write to $dmp_folder"   if ((stat($folder))[2] & 060) != 48; 
    }
    $folder;
}
