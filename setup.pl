#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use Sys::Hostname;
use Data::Dumper;

my $host = hostname;
my $base_path = $Bin;

my $tools_path = "$base_path/src/main/resources/tools/";
system("mkdir -p $tools_path");
system("rm -rf $tools_path/*");

# required perl modules
my @perl_modules=('CPAN::DistnameInfo','Try::Tiny','Config::Simple', 'Parallel::ForkManager', 'Log::Log4perl', 'Getopt::Long', 'Text::CSV', 'Text::Unidecode');

# required binaries
my %tools = ('esearch', '-', 'efetch', '-', 'xtract' , '-', 'fastx_collapser', '-', 'samtools', '-','kraken2', '-');

my $out_conf= "$base_path/conf.txt";

my $trimmomaticURL = "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip";
my $blastnURL = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz";
my $sratoolkitURL = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.5/sratoolkit.2.10.5-ubuntu64.tar.gz";
my $bowtie2URL = "https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip";
my $fastqcURL = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip";
my $kraken_viral_db = "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220908.tar.gz";

# edirect utils: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/install_ubuntu.txt
# samtools: http://www.htslib.org/download/


sub check_exists_command { 
    my $check = `sh -c 'command -v $_[0]'`; 
    return $check;
}

sub install_perl_module {
    my ($module) = (@_);
    if ($module eq '') {
        print "No Module name provided. Installer will exit.\n";
        exit();
    }
    eval {
        ( my $toload = "$module.pm") =~ s/::/\//g;
        require $toload;
        $module->import();
    };
    if($@) {
        print "Installing $module from CPAN: \n";
        system("yes | perl -MCPAN -e 'install $module'");
    }
    else{
        print "\n\nModule $module exists. Skipping installation\n";
    }
}

foreach my $mod(@perl_modules){
    install_perl_module($mod);
}

my $wget_path = `which wget`;
chomp $wget_path;
if ($wget_path eq ''){
    print "ERROR : Unable to find wget in the path. Please install wget utility in your system before starting the analysis\n";
}

# # setup trimmomatic
my $trimmomatic_path = `which trimmomatic`;
chomp $trimmomatic_path;
my $trimmomatic_adapter_PE = '<Please provide a valid path for the paired end adapter fasta file>';
my $trimmomatic_adapter_SE = '<Please provide a valid path for the single end adapter fasta file>';

if ($trimmomatic_path !~ /trimmomatic$/){
    system ("wget $trimmomaticURL -P $tools_path/");
    system ("unzip -o $tools_path/Trimmomatic-0.39.zip -d $tools_path/");
    system ("rm $tools_path/Trimmomatic-0.39.zip");
    if (-s "$tools_path/Trimmomatic-0.39/trimmomatic-0.39.jar"){
        $trimmomatic_path = "$tools_path/Trimmomatic-0.39/trimmomatic-0.39.jar";
        $trimmomatic_adapter_PE = "$tools_path/Trimmomatic-0.39/adapters/TruSeq2-PE.fa";
        $trimmomatic_adapter_SE = "$tools_path/Trimmomatic-0.39/adapters/TruSeq2-SE.fa";
    }
    else {
        $trimmomatic_path = '<Please provide a valid path for the executable jar file>';
        print "ERROR : Unable to setup Trimmomatic. Please provide a valid path for trimmomatic and the adapters in the config file: $out_conf\n";
    }
}
else {
    # $trimmomatic_path = "trimmomatic";
    my @tmp_path= split(/\//,$trimmomatic_path);
    pop @tmp_path; pop @tmp_path;
    $trimmomatic_path = join("/",@tmp_path).'/share/trimmomatic/trimmomatic.jar';
    if (! -s "$trimmomatic_path"){
        $trimmomatic_path = '<Please provide a valid path for the executable jar file>';
        print "ERROR : Unable to setup Trimmomatic. Please provide a valid path for trimmomatic and the adapters in the config file: $out_conf\n";
    }
    if (-s "$base_path/src/main/resources/adapters/TruSeq2-PE.fa"){
        $trimmomatic_adapter_PE = "$base_path/src/main/resources/adapters/TruSeq2-PE.fa";
    }
    if (-s "$base_path/src/main/resources/adapters/TruSeq2-SE.fa"){
        $trimmomatic_adapter_SE = "$base_path/src/main/resources/adapters/TruSeq2-SE.fa";
    }
}

# # setup blastn
my $blastn_path = `which blastn`;
my $makeblastdb_path = `which makeblastdb`;
chomp $blastn_path;
chomp $makeblastdb_path;
if ($blastn_path eq ''){
    system ("wget $blastnURL -P $tools_path/");
    system ("tar -xvf $tools_path/ncbi-blast-2.10.0+-x64-linux.tar.gz -C $tools_path/");
    system ("rm $tools_path/ncbi-blast-2.10.0+-x64-linux.tar.gz");
    if (! check_exists_command "$tools_path/ncbi-blast-2.10.0+/bin/blastn"){
        print "ERROR : Unable to setup ncbi-blast-2.10.0+-x64-linux. Please provide a valid path for NCBI blastn in the config file: $out_conf\n";
        $blastn_path = '<Please provide a valid path for the executable binary>';
    }
    else {
        $blastn_path = "$tools_path/ncbi-blast-2.10.0+/bin/blastn";
        $makeblastdb_path = "$tools_path/ncbi-blast-2.10.0+/bin/makeblastdb";
    }
}


# # setup sratoolkit
my $fastqdump_path = `which fastq-dump`;
chomp $fastqdump_path;
if ($fastqdump_path eq ""){
    system ("wget $sratoolkitURL -P $tools_path/");
    system ("tar -xvf $tools_path/sratoolkit.2.10.5-ubuntu64.tar.gz -C $tools_path/");
    system ("rm $tools_path/sratoolkit.2.10.5-ubuntu64.tar.gz");
    if (! check_exists_command "$tools_path/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump"){
        print "ERROR : Unable to setup sratoolkit.2.10.5-ubuntu64. Please provide a valid path for fastq-dump in the config file: $out_conf\n";
        $fastqdump_path = '<Please provide a valid path for the executable binary>';
    }
    else {
        $fastqdump_path = "$tools_path/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump";
    }

}

# # setup bowtie2
my $bowtie2_path = `which bowtie2`;
my $bowtie2_build_path = `which bowtie2-build`;
chomp $bowtie2_path;
chomp $bowtie2_build_path;
if ($bowtie2_path eq ""){
    system ("wget $bowtie2URL -P $tools_path/");
    system ("unzip -o $tools_path/bowtie2-2.4.5-linux-x86_64.zip -d $tools_path/");
    system ("rm $tools_path/bowtie2-2.4.5-linux-x86_64.zip");
    if (! check_exists_command "$tools_path/bowtie2-2.4.5-linux-x86_64/bowtie2"){
        print "ERROR : Unable to setup bowtie2. Please provide a valid path for bowtie2 in the config file: $out_conf\n";
    }
    else {
        $bowtie2_path = "$tools_path/bowtie2-2.4.5-linux-x86_64/bowtie2";
        $bowtie2_build_path = "$tools_path/bowtie2-2.4.5-linux-x86_64/bowtie2-build"; 
    }
}

# # setup fastqc
my $fastqc_path = `which fastqc`;
chomp $fastqc_path;
if ($fastqc_path eq ""){
    system ("wget $fastqcURL -P $tools_path/");
    system ("unzip -o $tools_path/fastqc_v0.11.9.zip -d $tools_path/");
    system ("rm $tools_path/fastqc_v0.11.9.zip");
    system ("chmod a+x $tools_path/FastQC/fastqc");
    if (! check_exists_command "$tools_path/FastQC/fastqc"){
        print "ERROR : Unable to setup fastqc. Please provide a valid path for fastqc in the config file: $out_conf\n";
    }
    else {
        $fastqc_path = "$tools_path/FastQC/fastqc";
    }
}

# # setup kraken2 database
my $kraken2_path = `which kraken2`;
my $kraken2_db_path = "$base_path/kraken-db";
chomp $kraken2_path;

if ($kraken2_path ne ''){
    system ("mkdir -p $kraken2_db_path");
    if (! -s "$kraken2_db_path/hash.k2d"){
        system ("wget $kraken_viral_db -P $kraken2_db_path/");
        system ("tar -xvf $kraken2_db_path/k2_viral_20220908.tar.gz -C $kraken2_db_path/");
        system ("rm $kraken2_db_path/k2_viral_20220908.tar.gz");
    }
}
else {
    print "ERROR : Unable to find kraken2 in PATH. Please provide a valid path for kraken2 in the config file: $out_conf\n";
    $kraken2_path = '<Please provide a valid path for the executable binary>';
}

foreach my $binary(keys %tools){
    my $tool_path = `which $binary`;
    chomp $tool_path;

    if ($tool_path eq ''){
         print "ERROR : Unable to find $binary in the path. Please provide a valid path for the executable binary in the config file: $out_conf\n";
         $tools{$binary} = '<Please provide a valid path for the executable binary>';
    }
    else {
        # print "$tool_path\n";
        $tools{$binary} = $tool_path;
    }
}

if ($tools{'esearch'}=~ /esearch$/){
    my @esearch_path_tmp = split(/\//,$tools{'esearch'});
    pop @esearch_path_tmp;
    my $esearch_path = join("/",@esearch_path_tmp);
    $tools{'esearch'} = $esearch_path;
}
if ($tools{'xtract'}=~ /xtract$/){
    my @xtract_path_tmp = split(/\//,$tools{'xtract'});
    pop @xtract_path_tmp;
    my $xtract_path = join("/",@xtract_path_tmp);
    $tools{'xtract'} = $xtract_path;
}
if ($tools{'efetch'}=~ /efetch$/){
    my @efetch_path_tmp = split(/\//,$tools{'efetch'});
    pop @efetch_path_tmp;
    my $efetch_path = join("/",@efetch_path_tmp);
    $tools{'efetch'} = $efetch_path;
}

my $config = <<"CONFIG";

## Empty values are not allowed.

[pipeline]

threads = 1
version = 1.6.0
########################################## end of pipeline config #####################################################

[fastq_dump]

execute = 1
## percentage of data to be downloaded using fastq dump:
data_perc = 5
########################################## end fastq dump config ####################################################

[fastqc]

execute_raw_data_qc = 1
execute_trimmed_data_qc = 1
threads = 2
########################################## end fastqc config ########################################################

[trimmomatic]

execute = 1
threads = 2
sliding_window = 4:20
min_len = 18
read_drop_perc_cutoff = 70
########################################## end trimmomatic config ###################################################

[blastn]

execute = 1
perc_identity = 90
qcov_hsp_perc = 80
num_threads = 2
alignment_perc_cutoff = 1
########################################## end blastn config ########################################################

[bowtie2]

execute = 1
threads = 2
end_to_end_preset = sensitive
alignment_perc_cutoff = 1
########################################## end bowtie2 config ########################################################

[kraken2]

execute = 1
threads = 2
kraken2_db_path = $kraken2_db_path
########################################## end bowtie2 config ########################################################

[toolPath]

fastq_dump = $fastqdump_path
fastqc = $fastqc_path
trimmomatic = $trimmomatic_path
adapters_PE = $trimmomatic_adapter_PE
adapters_SE = $trimmomatic_adapter_SE
blastn = $blastn_path
makeblastdb = $makeblastdb_path
bowtie2 = $bowtie2_path
bowtie2-build = $bowtie2_build_path
samtools = $tools{'samtools'}
ncbi_edirect = $tools{'efetch'}
fastx_collapser = $tools{'fastx_collapser'}
kraken2 = $tools{'kraken2'}
########################################## end of tools path config #######################################################
CONFIG


print "Configuration file generated : $out_conf\n";

open (CONF, ">$out_conf");
print CONF "$config";
close CONF;
