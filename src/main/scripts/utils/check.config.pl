#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Config::Simple;
use Log::Log4perl qw(get_logger);

my $input_config = $ARGV[0];
my $tool_config_obj = new Config::Simple($input_config);

my $ncbi_edirect_path   = $tool_config_obj->param('toolPath.ncbi_edirect');
my $fastq_dump          = $tool_config_obj->param('toolPath.fastq_dump');
my $fastqc              = $tool_config_obj->param('toolPath.fastqc');
my $trimmomatic         = $tool_config_obj->param('toolPath.trimmomatic');
my $adapters_PE         = $tool_config_obj->param('toolPath.adapters_PE');
my $adapters_SE         = $tool_config_obj->param('toolPath.adapters_SE');
my $fastx_collapser     = $tool_config_obj->param('toolPath.fastx_collapser');
my $blastn              = $tool_config_obj->param('toolPath.blastn');
my $bowtie2             = $tool_config_obj->param('toolPath.bowtie2');
my $samtools            = $tool_config_obj->param('toolPath.samtools');

my $pipeline_threads	= $tool_config_obj->param('pipeline.threads');

my $fastq_dump_perc		= $tool_config_obj->param('fastq_dump.data_perc');

my $fastqc_threads      = $tool_config_obj->param('fastqc.threads');

my $trimmomatic_threads          = $tool_config_obj->param('trimmomatic.threads');
my $trimmomatic_sliding_win      = $tool_config_obj->param('trimmomatic.sliding_window');
my $trimmomatic_min_len          = $tool_config_obj->param('trimmomatic.min_len');
my $trimmomatic_read_drop_cutoff = $tool_config_obj->param('trimmomatic.read_drop_perc_cutoff');

my $blastn_perc_identity    = $tool_config_obj->param('blastn.perc_identity');
my $blastn_qcov_hsp_perc    = $tool_config_obj->param('blastn.qcov_hsp_perc');
my $blastn_num_threads      = $tool_config_obj->param('blastn.num_threads');
my $blastn_align_cutoff		= $tool_config_obj->param('blastn.alignment_perc_cutoff');
my $blastn_db_path          = $tool_config_obj->param('blastn.db_path');

my $bowtie2_threads         = $tool_config_obj->param('bowtie2.threads');
my $bowtie2_e2e_preset		= $tool_config_obj->param('bowtie2.end_to_end_preset');
my $bowtie2_align_cutoff	= $tool_config_obj->param('bowtie2.alignment_perc_cutoff');
my $bowtie2_index_path      = $tool_config_obj->param('bowtie2.bowtie_index_path');




my $conf=q{
log4perl.rootLogger = DEBUG, Screen
log4perl.appender.Screen = Log::Log4perl::Appender::Screen
log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
log4perl.appender.Screen.layout.ConversionPattern =  %p %d %F{1} %L > %m %n
log4perl.appender.Screen.color.TRACE=cyan
log4perl.appender.Screen.color.DEBUG=blue
log4perl.appender.Screen.color.ERROR=bold red
};

## Initialize Log4perl, if it is not initialized already by the calling script
unless (Log::Log4perl::initialized()) {
	Log::Log4perl::init( \$conf );
}
my $logger = Log::Log4perl->get_logger();

sub check_exists_command { 
    my $check = `sh -c 'command -v $_[0]'`; 
    return $check;
}

my $flag = 0;


if (! check_exists_command "$ncbi_edirect_path/esearch"){
    $logger->error("Incorrect Configuration : toolPath > ncbi_edirect : Please provide the installation directory of edirect utils.");
    $flag = 1;
}

if (! check_exists_command "$fastq_dump"){
    $logger->error("Incorrect Configuration : toolPath > fastq_dump : Please provide the correct executable for fastq-dump.");
    $flag = 1;
}

if (! check_exists_command "$fastqc"){
    $logger->error("Incorrect Configuration : toolPath > fastqc : Please provide the correct executable for fastqc.");
    $flag = 1;
}

if (! check_exists_command "$trimmomatic"){
    $logger->error("Incorrect Configuration : toolPath > trimmomatic : Please provide the correct executable for trimmomatic.");
    $flag = 1;
}

if (! check_exists_command "$fastx_collapser"){
    $logger->error("Incorrect Configuration : toolPath > fastx_collapser : Please provide the correct executable for fastx_collapser.");
    $flag = 1;
}

if (! check_exists_command "$blastn"){
    $logger->error("Incorrect Configuration : toolPath > blastn : Please provide the correct executable for blastn.");
    $flag = 1;
}

if (! check_exists_command "$bowtie2"){
    $logger->error("Incorrect Configuration : toolPath > bowtie2 : Please provide the correct executable for bowtie2.");
    $flag = 1;
}

if (! check_exists_command "$samtools"){
    $logger->error("Incorrect Configuration : toolPath > samtools : Please provide the correct executable for samtools.");
    $flag = 1;
}

# if (! glob ("$blastn_db_path*")){
#     $logger->error("Incorrect Configuration : blastn > db_path : Please provide the path for the blast database.");
#     $flag = 1;
# }

# if (! glob ("$bowtie2_index_path*")){
#     $logger->error("Incorrect Configuration : bowtie2 > bowtie_index_path : Please provide the path for the bowtie2 index.");
#     $flag = 1;
# }

if (! -s "$adapters_PE"){
    $logger->error("Incorrect Configuration : toolPath > adapters_PE : Please check the paired end adapter fasta.");
    $flag = 1;
}

if (! -s "$adapters_SE"){
    $logger->error("Incorrect Configuration : toolPath > adapters_SE : Please check the single end adapter fasta.");
    $flag = 1;
}

if ($pipeline_threads !~ /^\d+$/){
    $logger->error("Incorrect Configuration : pipeline > threads : Numerical value expected.");
    $flag = 1;
}

if (($fastq_dump_perc !~ /^\d+$/) or ($fastq_dump_perc < 1 ) ){
    $logger->error("Incorrect Configuration : fastq_dump > data_perc : Numerical value expected (>=1).");
    $flag = 1;
}

if ($fastqc_threads !~ /^\d+$/){
    $logger->error("Incorrect Configuration : fastqc > threads : Numerical value expected.");
    $flag = 1;
}

if ($trimmomatic_threads !~ /^\d+$/){
    $logger->error("Incorrect Configuration : trimmomatic > threads : Numerical value expected.");
    $flag = 1;
}

if ($trimmomatic_sliding_win !~ /\d+:\d+$/){
    $logger->error("Incorrect Configuration : trimmomatic > sliding_window : expected format: 4:20.");
    $flag = 1;
}

if ($trimmomatic_min_len !~ /^\d+$/){
    $logger->error("Incorrect Configuration : trimmomatic > min_len : Numerical value expected.");
    $flag = 1;
}

if ($trimmomatic_read_drop_cutoff !~ /^\d+$/){
    $logger->error("Incorrect Configuration : trimmomatic > read_drop_perc_cutoff : Numerical value expected.");
    $flag = 1;
}

if ($blastn_perc_identity !~ /^\d+$/){
    $logger->error("Incorrect Configuration : blastn > perc_identity : Numerical value expected.");
    $flag = 1;
}

if ($blastn_qcov_hsp_perc !~ /^\d+$/){
    $logger->error("Incorrect Configuration : blastn > qcov_hsp_perc : Numerical value expected.");
    $flag = 1;
}

if ($blastn_num_threads !~ /^\d+$/){
    $logger->error("Incorrect Configuration : blastn > num_threads : Numerical value expected.");
    $flag = 1;
}

if ($bowtie2_threads !~ /^\d+$/){
    $logger->error("Incorrect Configuration : bowtie2 > threads : Numerical value expected.");
    $flag = 1;
}

if ($bowtie2_e2e_preset !~ /^sensitive$|^very-sensitive$|^fast$|^very-fast$/){
    $logger->error("Incorrect Configuration : bowtie2 > end_to_end_preset : Please provide one of the following bowtie2 presets: 'sensitive' or 'very-sensitive' or 'fast' or 'very-fast'");
    $flag = 1;
}

if ($flag == 1){
    die $logger->error("Configuration error");
}

# check_exists_command "$ncbi_edirect_path/esearch" or $logger->error("Incorrect Configuration : toolPath > ncbi_edirect : Please provide the installation directory of edirect utils.");
# check_exists_command "$fastq_dump" or $logger->error("Incorrect Configuration : toolPath > fastq_dump : Please provide the correct executable for fastq-dump.");
# check_exists_command "$fastqc" or $logger->error("Incorrect Configuration : toolPath > fastqc : Please provide the correct executable for fastqc.");
# check_exists_command "$trimmomatic" or $logger->error("Incorrect Configuration : toolPath > trimmomatic : Please provide the correct executable for trimmomatic.");
# check_exists_command "$fastx_collapser" or $logger->error("Incorrect Configuration : toolPath > fastx_collapser : Please provide the correct executable for fastx_collapser.");
# check_exists_command "$blastn" or $logger->error("Incorrect Configuration : toolPath > blastn : Please provide the correct executable for blastn.");
# check_exists_command "$bowtie2" or $logger->error("Incorrect Configuration : toolPath > bowtie2 : Please provide the correct executable for bowtie2.");
# check_exists_command "$samtools" or $logger->error("Incorrect Configuration : toolPath > samtools : Please provide the correct executable for samtools.");