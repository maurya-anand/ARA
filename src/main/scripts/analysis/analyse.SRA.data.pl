#!/usr/bin/env perl


################ Pipeline to retrieve raw fastq file for every SRA run
################ and execute QC, trimmomatic and blastn.


use strict;
use warnings;
use Data::Dumper;
use Config::Simple;
use Log::Log4perl qw(get_logger);
use FindBin qw($Bin);
use Sys::Hostname;
use File::Basename;

my $host = hostname;
my $base_path = $Bin;
my $tool_config_obj;


my $input_config = $ARGV[0];
my $sra_sample_id = $ARGV[1];
my $sample_type = $ARGV[2];
my $num_spots = $ARGV[3];
my $analysis_dir = $ARGV[4];
my $inp_sequences = $ARGV[5];
my $references_dir = $ARGV[6];
my $fastq_dump_perc = $ARGV[7];


my $conf=q{
log4perl.rootLogger = DEBUG, Logfile
log4perl.appender.Screen = Log::Log4perl::Appender::Screen
log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
log4perl.appender.Screen.layout.ConversionPattern =  %p %d %F{1} > %m %n
log4perl.appender.Screen.color.TRACE=cyan
log4perl.appender.Screen.color.DEBUG=blue
log4perl.appender.Screen.color.ERROR=bold red
log4perl.appender.Logfile = Log::Log4perl::Appender::File
log4perl.appender.Logfile.filename = sub { logfile(); };
log4perl.appender.Logfile.layout = Log::Log4perl::Layout::PatternLayout
log4perl.appender.Logfile.layout.ConversionPattern = %p %d %F{1} > %m %n
};

## Initialize Log4perl, if it is not initialized already by the calling script
unless (Log::Log4perl::initialized()) {
	Log::Log4perl::init( \$conf );
}
my $logger = Log::Log4perl->get_logger();

#For dynamic logfile name
sub logfile {
    return "$analysis_dir/runlog.$sra_sample_id.txt";
}

if ((!$input_config)||(!$sra_sample_id)||(!$sample_type)||(!$num_spots)||(!$analysis_dir)){

	$logger->error("Check input parameters");
	$logger->error("perl sra_rna_seq_wrapper.pl <SRA id> <sample type: PE or SE> <output directory path> <optional config file>");
	$logger->error("perl sra_rna_seq_wrapper.pl SRR7266357 SE /home/anand/Documents/pipeline.SRA.RNA.seq/src/main/test/pipeline");
	exit;
}
## Loading the pipeline configuration
else {

	$logger->info("Pipeline : Config: $input_config");
    
    if ($fastq_dump_perc){
        $logger->info("Pipeline : Params: Run_ID=$sra_sample_id,Run_type=$sample_type,Fastq_subset_percentage=$fastq_dump_perc,Output_dir=$analysis_dir"); 
    }
    else {
        $logger->info("Pipeline : Params: Run_ID=$sra_sample_id,Run_type=$sample_type,Output_dir=$analysis_dir"); 
    }
   $tool_config_obj = new Config::Simple($input_config);
}


 
my $ncbi_edirect_path	= $tool_config_obj->param('toolPath.ncbi_edirect');
my $fastq_dump      	= $tool_config_obj->param('toolPath.fastq_dump');
my $fastqc          	= $tool_config_obj->param('toolPath.fastqc');
my $trimmomatic 		= $tool_config_obj->param('toolPath.trimmomatic');
my $adapters_PE     	= $tool_config_obj->param('toolPath.adapters_PE');
my $adapters_SE 		= $tool_config_obj->param('toolPath.adapters_SE');
my $fastx_collapser 	= $tool_config_obj->param('toolPath.fastx_collapser');
my $blastn          	= $tool_config_obj->param('toolPath.blastn');
my $makeblastdb			= $tool_config_obj->param('toolPath.makeblastdb');
my $bowtie2         	= $tool_config_obj->param('toolPath.bowtie2');
my $bowtie2_build		= $tool_config_obj->param('toolPath.bowtie2-build');
my $kraken2         	= $tool_config_obj->param('toolPath.kraken2');
my $samtools        	= $tool_config_obj->param('toolPath.samtools');
my $pipeline_ver    	= $tool_config_obj->param('pipeline.version');

## Loading tool parameters
my $fastqc_threads             = $tool_config_obj->param('fastqc.threads');

my $trimmomatic_threads          = $tool_config_obj->param('trimmomatic.threads');
my $trimmomatic_sliding_win      = $tool_config_obj->param('trimmomatic.sliding_window');
my $trimmomatic_min_len          = $tool_config_obj->param('trimmomatic.min_len');
my $trimmomatic_read_drop_cutoff = $tool_config_obj->param('trimmomatic.read_drop_perc_cutoff');

my $blastn_perc_identity    = $tool_config_obj->param('blastn.perc_identity');
my $blastn_qcov_hsp_perc    = $tool_config_obj->param('blastn.qcov_hsp_perc');
my $blastn_num_threads      = $tool_config_obj->param('blastn.num_threads');
my $blastn_align_cutoff		= $tool_config_obj->param('blastn.alignment_perc_cutoff');
# my $blastn_db_path          = $tool_config_obj->param('blastn.db_path');

my $bowtie2_threads         = $tool_config_obj->param('bowtie2.threads');
my $bowtie2_e2e_preset		= $tool_config_obj->param('bowtie2.end_to_end_preset');
my $bowtie2_align_cutoff	= $tool_config_obj->param('bowtie2.alignment_perc_cutoff');
my $kraken2_threads         = $tool_config_obj->param('kraken2.threads');
my $kraken2_db         		= $tool_config_obj->param('kraken2.kraken2_db_path');
# my $bowtie2_index_path      = $tool_config_obj->param('bowtie2.bowtie_index_path');

# my $fastq_dump_exec_flag = 1; my $fastqc_exec_flag = 1; my $trimmomatic_exec_flag = 1; my $blastn_exec_flag = 1; 

## Loading the tool execution flags
my $fastq_dump_exec_flag            = $tool_config_obj->param('fastq_dump.execute');
my $fastqc_exec_raw_data_flag       = $tool_config_obj->param('fastqc.execute_raw_data_qc');
my $trimmomatic_exec_flag           = $tool_config_obj->param('trimmomatic.execute');
my $fastqc_exec_trimmed_data_flag   = $tool_config_obj->param('fastqc.execute_trimmed_data_qc'); 
my $blastn_exec_flag                = $tool_config_obj->param('blastn.execute');
my $bowtie2_exec_flag				= $tool_config_obj->param('bowtie2.execute');
my $kraken2_exec_flag				= $tool_config_obj->param('kraken2.execute');

# my $exec_SRA_analysis_flag	= $tool_config_obj->param('sra_metadata.raw_data_complete_analysis');

my $exec_cov_analysis_flag = 0;

if ($fastq_dump_perc) {
    $exec_cov_analysis_flag = 0;
}

if ($fastqc_exec_raw_data_flag ==1){
	$fastq_dump_exec_flag = 1;
}

if ($trimmomatic_exec_flag == 1){
	$fastq_dump_exec_flag = 1;
	$fastqc_exec_raw_data_flag = 1;
	$fastqc_exec_trimmed_data_flag = 1;
}

if (($blastn_exec_flag == 1) or ($bowtie2_exec_flag == 1)){
	$fastq_dump_exec_flag = 1;
	$fastqc_exec_raw_data_flag = 1;
	$fastqc_exec_trimmed_data_flag = 1;
	$trimmomatic_exec_flag = 1;
}

## if decimal is present, then extract the number.. else.. substitue all no digit characters with nothing
if ($blastn_align_cutoff =~ /(^\d+(?:\.\d+)?)/){
	$blastn_align_cutoff = $1;
}
else {
	$blastn_align_cutoff =~ s/\D//g;
}
if ($bowtie2_align_cutoff =~ /(^\d+(?:\.\d+)?)/){
	$bowtie2_align_cutoff = $1;
}
else {
	$bowtie2_align_cutoff =~ s/\D//g;
}

$logger->info("Pipeline : v.$pipeline_ver executed on $host");

# my $references_dir = "$analysis_dir/reference";
system("mkdir -p $references_dir");

my $output_dir = "$analysis_dir/$sra_sample_id";
if (-d "$output_dir"){
	$logger->info("The target working directory $output_dir already exists. Please delete the $sra_sample_id directory or provide a different output destination");
	$logger->info("Analysis skipped for $sra_sample_id");
	exit;
}
else {
	system("mkdir -p $output_dir");
}


# fallback method
my $sra_url;

my $sra_run_info_str = `$ncbi_edirect_path/efetch -db sra -id $sra_sample_id -format runinfo |tail -n1`;
chomp $sra_run_info_str;
# Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
# SRR12548243,2021-04-07 17:04:56,2020-08-30 19:21:08,1051,79876,1051,76,0,,https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12548243/SRR12548243,SRX9036864,,RNA-Seq,cDNA,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,NextSeq 500,SRP279357,PRJNA660263,3,660263,SRS7287368,SAMN15943290,simple,3702,Arabidopsis thaliana,GSM4755568,,,,,,,no,,,,,GEO,SRA1119003,,public,D21CE888E8EF2BA2AE78A4FE7577C51F,9EE437AA19CFF2DAD503DF6FDB4EC570

if ($sra_run_info_str){
	my @tmp=split(",",$sra_run_info_str);
	
	my $samType_re_check = $tmp[15];
	my $sra_url_check = $tmp[9];
	my $run_spots_re_check = $tmp[3];

	if (($samType_re_check eq "PAIRED") or ($samType_re_check eq "SINGLE")){
		if($sample_type ne $samType_re_check){
			$logger->info("Using corrected sample type ($samType_re_check instead of $sample_type) for $sra_sample_id");
			$sample_type = $samType_re_check;
		}
	}
	if ($run_spots_re_check){
		if($run_spots_re_check ne $num_spots){
			$logger->info("Using corrected number of spots ($run_spots_re_check instead of $num_spots) for $sra_sample_id");
			$num_spots = $run_spots_re_check;
		}
	}
	if (($sra_url_check =~ /^http/) && ($sra_url_check =~ m/$sra_sample_id/)){
		$sra_url = $sra_url_check;
	}
	else {
		$sra_url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$sra_sample_id/$sra_sample_id";
		# sratoolkit.2.10.5-ubuntu64/bin/srapath SRR000001 returns https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001
	}
}

## Execute fastq dump ###############################################################
if ($fastq_dump_exec_flag == 1 ) {
    
    my $output_dir = "$analysis_dir/$sra_sample_id/raw_fastq";
    system("mkdir -p $output_dir");

	my $dump_process_flag = download_sample($sra_sample_id, $sample_type, $output_dir);

	if ($dump_process_flag == 0){
		$logger->error("Fastq dump failed for $sra_sample_id");
		my $full_sra_download_flag = download_sample_full_sra ($sra_sample_id, $sample_type, $output_dir);
		if ($full_sra_download_flag ==0){
			$logger->error("Unable to retrive raw data for $sra_sample_id");
			exit;
		}
	}
}

## Execute raw fastQC ###############################################################
if ($fastqc_exec_raw_data_flag == 1 ) {

    my $qc_output_dir = "$analysis_dir/$sra_sample_id/fastQC";
    system("mkdir -p $qc_output_dir");
    
	my @fastFilesList;
    my $qc_output_summary;

	if ($sample_type eq "PAIRED"){

		my $fastqn1="$output_dir/raw_fastq/$sra_sample_id\_1.fastq";
		my $fastqn2="$output_dir/raw_fastq/$sra_sample_id\_2.fastq";

        $qc_output_summary = "$analysis_dir/$sra_sample_id/fastQC/$sra_sample_id\_1_fastqc/summary.txt";
		
		if ((-s "$fastqn1") && (-s "$fastqn2")){

			push @fastFilesList, $fastqn1;

			push @fastFilesList, $fastqn2;
		}
	}

	if ($sample_type eq "SINGLE") {

		my $fastqn1="$output_dir/raw_fastq/$sra_sample_id.fastq";
        
        $qc_output_summary = "$analysis_dir/$sra_sample_id/fastQC/$sra_sample_id\_fastqc/summary.txt";
		
		if(-s "$fastqn1"){

			push @fastFilesList, $fastqn1;
		}
	}
    if (! -s "$qc_output_summary"){

        $logger->info("Executing: FastQC on raw fastq files");
        
        run_qc($sra_sample_id, $qc_output_dir, @fastFilesList);
        
        $logger->info("FastQC on raw fastq files completed");
    }
    else{
        
        $logger->info("FastQC on raw fastq files completed");
    }
	
}

## Execute trimmomatic ##############################################################
if ($trimmomatic_exec_flag == 1 ) {

    my $trimmomatic_output_dir = "$analysis_dir/$sra_sample_id/trimmed_data";
    system("mkdir -p $trimmomatic_output_dir");

	if ($sample_type eq "PAIRED"){

		my $fastqn1="$output_dir/raw_fastq/$sra_sample_id\_1.fastq";
		my $fastqn2="$output_dir/raw_fastq/$sra_sample_id\_2.fastq";
	
		if ((-s "$fastqn1") && (-s "$fastqn2")){

			adapter_trim_PE($sra_sample_id, $fastqn1, $fastqn2, $trimmomatic_output_dir);
			# system("rm $output_dir/raw_fastq/*.fastq");
		}
	}

	if ($sample_type eq "SINGLE") {

		my $fastqn1="$output_dir/raw_fastq/$sra_sample_id.fastq";

		if(-s "$fastqn1"){

			adapter_trim_SE($sra_sample_id, $fastqn1, $trimmomatic_output_dir);
			# system("rm $output_dir/raw_fastq/*.fastq");
		}
	}
}

## Execute blastn ###################################################################
if ($blastn_exec_flag == 1 ) {
	
	my $sampleTrimmedFastq;
	my $trimmedCollapsedReHeadFasta;

	my $trimmomatic_output_dir = "$analysis_dir/$sra_sample_id/trimmed_data";
    
    my $blast_analysis_dir="$analysis_dir/$sra_sample_id/blastn";
    system("mkdir -p $blast_analysis_dir");

    my $blast_results = "$blast_analysis_dir/$sra_sample_id.blast.results.txt";
    my $blast_filtered_results = "$blast_analysis_dir/$sra_sample_id\_blast_best_hits.txt";
    my $blast_filtered_bed = "$blast_analysis_dir/$sra_sample_id\_blast_best_hits_sorted.bed";

	my $blast_hits_status_logfile = "$blast_analysis_dir/blast.stats.txt";

	my $blastn_db_path =  create_blastn_db($inp_sequences,$references_dir);


	if ($sample_type eq "PAIRED") {

        $sampleTrimmedFastq="$analysis_dir/$sra_sample_id/trimmed_data/$sra_sample_id\_COM\_AT.fastq";
    }

	if ($sample_type eq "SINGLE") {

        $sampleTrimmedFastq="$analysis_dir/$sra_sample_id/trimmed_data/$sra_sample_id\_AT.fastq";
    }

	if (-s "$sampleTrimmedFastq"){

        $trimmedCollapsedReHeadFasta = fastq_to_fasta($sra_sample_id, $sampleTrimmedFastq, $trimmomatic_output_dir);
	}
    
    my $blast_outfmt= '"6 qseqid qstart qend length nident qlen qseq sseqid sstart send slen sstrand sseq pident bitscore evalue"';

    if ((! -s "$blast_results") && (-s "$trimmedCollapsedReHeadFasta")){

        $logger->info("Executing: Blastn");
        my $execBlastComm = "$blastn -db $blastn_db_path -query $trimmedCollapsedReHeadFasta -out $blast_results -perc_identity $blastn_perc_identity -qcov_hsp_perc $blastn_qcov_hsp_perc -num_threads $blastn_num_threads -outfmt $blast_outfmt";
        $logger->info("Blastn : Command: $execBlastComm");
        system ("$execBlastComm");

        if (! -s "$blast_results"){
            $logger->error("Blastn : Result : No hits found");
        }
		else {
			$logger->info("Blastn for $sra_sample_id completed");
		}
    }

	if (! -s "$blast_results"){
		$logger->error("Blastn for $sra_sample_id failed");
	}
	
	if ((! -s "$blast_hits_status_logfile") && (-s "$blast_results")){
		my ($percblastHits,$percMappedReads,$retStr)=check_blast_hits($sra_sample_id,$trimmedCollapsedReHeadFasta,$blast_results);
		open(STATUS, ">$blast_hits_status_logfile") or $logger->logdie("Cannot write to file  - $blast_hits_status_logfile: $!");
		print STATUS "Unique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tStatus\n";
		
		if ($percMappedReads >= $blastn_align_cutoff){
			# $logger->info("Blastn : Result : ID\tUnique_seqs\tUnique_Query_Ids\tUnique_tRNA_IDs\tHits_Perc");
			$logger->info("Blastn : Result : Unique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast");
			$logger->info("Blastn : Result : PASS\t$retStr");
			print STATUS "$retStr\tPASS\n";
		}
		else{
			print STATUS "$retStr\tFAIL\n";
			# $logger->error("Blastn : Result : Hits less than 0.01%");
			$logger->error("Blastn : Result : Hits less than $blastn_align_cutoff %");
		}
		close STATUS;
	}
}

if ($bowtie2_exec_flag == 1 ) {
    
    my $trimmomatic_output_dir = "$analysis_dir/$sra_sample_id/trimmed_data";

    my $bowtie_analysis_dir="$analysis_dir/$sra_sample_id/bowtie2";
    my $aligned_bam = "$bowtie_analysis_dir/$sra_sample_id.bam";
    my $alignment_stats = "$bowtie_analysis_dir/alignment.stats.txt";
	my $alignment_stats_out = "$bowtie_analysis_dir/alignment.txt";
    
    my $bowtie2_comm;
	my $blast_analysis_dir="$analysis_dir/$sra_sample_id/blastn";
	my $trimmedCollapsedReHeadFasta;
	
	my $sampleTrimmedFastq;

	my $bowtie2_index_path =  create_bowtie2_index($inp_sequences,$references_dir);

	if ($sample_type eq "PAIRED") {

		$sampleTrimmedFastq="$analysis_dir/$sra_sample_id/trimmed_data/$sra_sample_id\_COM\_AT.fastq";
	}

	if ($sample_type eq "SINGLE") {

		$sampleTrimmedFastq="$analysis_dir/$sra_sample_id/trimmed_data/$sra_sample_id\_AT.fastq";
	}

	if (-s "$sampleTrimmedFastq"){

		$trimmedCollapsedReHeadFasta = fastq_to_fasta($sra_sample_id, $sampleTrimmedFastq, $trimmomatic_output_dir);
		# system("rm $analysis_dir/$sra_sample_id/trimmed_data/*.fastq");
	}

    if ((! -s "$alignment_stats") && (-s "$trimmedCollapsedReHeadFasta")){
        $logger->info("Executing : Bowtie2");
        system("mkdir -p $bowtie_analysis_dir");
		# running bowtie on collapsed unique fasta
		$bowtie2_comm = "$bowtie2 -x $bowtie2_index_path -f $trimmedCollapsedReHeadFasta --no-unal --end-to-end --$bowtie2_e2e_preset -p $bowtie2_threads 2> $alignment_stats | $samtools sort -@ $bowtie2_threads -T $analysis_dir/ -o $aligned_bam";
		$logger->info("Bowtie2 : Command : $bowtie2_comm");
		my $exec_bowtie2_stdout = `$bowtie2_comm`;
    }
	else {
		$logger->error("Bowtie2 for $sra_sample_id failed");
	}
	
	if (-s "$alignment_stats"){
		my ($percbowtieHits,$percTotalMappedReads,$retStr)= check_bowtie2_hits($bowtie_analysis_dir,$aligned_bam,$trimmedCollapsedReHeadFasta);
		
		open(BOWTIESTATS, ">$alignment_stats_out") or $logger->logdie("Cannot write to file  - $alignment_stats_out: $!");

		print BOWTIESTATS "Unique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tStatus\n";
		if ($percTotalMappedReads >= $bowtie2_align_cutoff){
			$logger->info("Bowtie2 : Result : Unique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2");
			$logger->info("Bowtie2 : Result : PASS\t$retStr");

			print BOWTIESTATS "$retStr\tPASS\n";
		}
		else{
			print BOWTIESTATS "$retStr\tFAIL\n";
			$logger->error("Bowtie2 : Result : Hits less than $bowtie2_align_cutoff %");
		}
		close BOWTIESTATS;
	}
}

if ($kraken2_exec_flag == 1 ) {
    my $trimmomatic_output_dir = "$analysis_dir/$sra_sample_id/trimmed_data";

    my $kraken2_output_dir = "$analysis_dir/$sra_sample_id/kraken2";

    my $kraken2_output_file="$kraken2_output_dir/$sra_sample_id.kraken";
	my $kraken2_report_file="$kraken2_output_dir/$sra_sample_id.report";
	my $kraken2_stdout = "$kraken2_output_dir/$sra_sample_id.stdout.txt";

	my $trimmedCollapsedReHeadFasta;
	my $sampleTrimmedFastq;

	if ($sample_type eq "PAIRED") {

		$sampleTrimmedFastq="$analysis_dir/$sra_sample_id/trimmed_data/$sra_sample_id\_COM\_AT.fastq";
	}

	if ($sample_type eq "SINGLE") {

		$sampleTrimmedFastq="$analysis_dir/$sra_sample_id/trimmed_data/$sra_sample_id\_AT.fastq";
	}

	if (-s "$sampleTrimmedFastq"){

		$trimmedCollapsedReHeadFasta = fastq_to_fasta($sra_sample_id, $sampleTrimmedFastq, $trimmomatic_output_dir);
		# system("rm $analysis_dir/$sra_sample_id/trimmed_data/*.fastq");
	}

    if ((! -s "$kraken2_report_file") && (-s "$trimmedCollapsedReHeadFasta")){
        $logger->info("Executing : Kraken2");
        system("mkdir -p $kraken2_output_dir");
		# running kraken on collapsed unique fasta
		my $kraken2_comm = "$kraken2 --db $kraken2_db --threads $kraken2_threads $trimmedCollapsedReHeadFasta --output $kraken2_output_file --report $kraken2_report_file >> $kraken2_stdout 2>&1";
		$logger->info("Kraken2 : Command : $kraken2_comm");
		system("$kraken2_comm");
    }
	
	if (-s "$kraken2_report_file"){
		$logger->info("Kraken2 : Result : Classification file generated: $kraken2_report_file");
	} else {
		$logger->error("Kraken2 : Failed");
	}
}


## Functions ########################################################################
sub adapter_trim_PE {

	my ($id,$fastq1,$fastq2,$outPath)= @_;
	
	my $trimmedFastqPE1="$outPath/$id\_1_AT.fastq";
	my $trimmedFastqPE2="$outPath/$id\_2_AT.fastq";
	my $trimmedFastqUnPaired1="$outPath/$id\_1_unpaired.fastq";
	my $trimmedFastqUnPaired2="$outPath/$id\_2_unpaired.fastq";
	my $trimmedFastqLogPE="$outPath/$id\_trim_log.txt";
	my $trimmedFastqLogStdoutPE="$outPath/$id\_trim_stdout_log.txt";

    my @fastFilesList= ($trimmedFastqPE1,$trimmedFastqPE2,$trimmedFastqUnPaired1,$trimmedFastqUnPaired2);
    my $concatenatedPEfastqAT= "$outPath/$id\_COM_AT.fastq";
	
	my $trimComm="java -jar $trimmomatic PE -threads $trimmomatic_threads -phred33 -trimlog $trimmedFastqLogPE $fastq1 $fastq2 $trimmedFastqPE1 $trimmedFastqUnPaired1 $trimmedFastqPE2 $trimmedFastqUnPaired2 ILLUMINACLIP:$adapters_PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:$trimmomatic_sliding_win MINLEN:$trimmomatic_min_len >> $trimmedFastqLogStdoutPE 2>&1";
    my @exec_trimmomaticPE_stdout;

    if (! -s "$trimmedFastqLogStdoutPE"){

        $logger->info("Executing: Trimmomatic");
        $logger->info("Param: Trimmomatic : adapters = $adapters_PE");
        $logger->info("Param: Trimmomatic : SLIDINGWINDOW = $trimmomatic_sliding_win");
        $logger->info("Param: Trimmomatic : MINLEN = $trimmomatic_min_len");
        $logger->info("Command: $trimComm");

        @exec_trimmomaticPE_stdout=`java -jar $trimmomatic PE -threads $trimmomatic_threads -phred33 -trimlog $trimmedFastqLogPE $fastq1 $fastq2 $trimmedFastqPE1 $trimmedFastqUnPaired1 $trimmedFastqPE2 $trimmedFastqUnPaired2 ILLUMINACLIP:$adapters_PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:$trimmomatic_sliding_win MINLEN:$trimmomatic_min_len 2>&1`;
        open (TRIMLOG, ">$trimmedFastqLogStdoutPE");
        print TRIMLOG @exec_trimmomaticPE_stdout;
        close TRIMLOG;

        if (! -s "$concatenatedPEfastqAT"){
            my $concaPEfastComm= "cat ".join(" ", @fastFilesList)." > $concatenatedPEfastqAT";
            system("$concaPEfastComm");
        }
        $logger->info("Trimmomatic for $id completed");
		system("rm $trimmedFastqLogPE");
    }
    else{

        $logger->info("Trimmomatic for $id  completed");
    }

	if ((-s "$trimmedFastqPE1") && (-s "$trimmedFastqPE2")){

		if ($fastqc_exec_trimmed_data_flag == 1 ) {

            my $qc_output_dir = "$analysis_dir/$id/fastQC";
            system("mkdir -p $qc_output_dir");

            my $qc_output_summary = "$analysis_dir/$id/fastQC/$id\_1_AT_fastqc/summary.txt";
            
            $logger->info("Checking sample quality of $id ");

            if (! -s "$qc_output_summary"){

                $logger->info("Executing: FastQC on adapter trimmed fastq files");
                run_qc($id, $qc_output_dir, @fastFilesList);
                $logger->info("FastQC on adapter trimmed fastq files completed");
            }
            else{

                $logger->info("FastQC on adapter trimmed fastq files completed");
            }

            my @trimLogArr  = `cat $trimmedFastqLogStdoutPE`;
            chomp $trimLogArr[-2];
            my @trimLogStat = split(/\s/,$trimLogArr[-2]); #Input Read Pairs: 1467341 Both Surviving: 772840 (52.67%) Forward Only Surviving: 92309 (6.29%) Reverse Only Surviving: 105347 (7.18%) Dropped: 496845 (33.86%)
			$trimLogStat[-1]=~ s/[()%]//g;

			if ($trimLogStat[-1] >= $trimmomatic_read_drop_cutoff){ # if the trimmomatic read drop percentage criteria is not satisfied 

				$logger->error("Trimmomatic QC failed for $id");
				$logger->error("Sample: $id $trimLogArr[-2]");
				exit;
			}
			else { # open each fastqc summary file and look for "PASS" flag for the given qc criterias

                $logger->info("Trimmomatic : Result: Dropped reads $trimLogStat[-1]%.");
                $logger->info("$id passed trimmomatic minimum read threshold"); 
				checkSampleQuality($id,$outPath,$qc_output_dir);
			}
		}
	}

	else { # if adapter trimmed R1 or R2 files are empty

		$logger->error("Trimmomatic failed for $sra_sample_id");

		exit;
	}
}

sub adapter_trim_SE {

	my ($id,$fastq1,$outPath)= @_;

	my $trimmedFastqSE="$outPath/$id\_AT.fastq";
	my $trimmedFastqLogSE="$outPath/$id\_trim_log.txt";
	my $trimmedFastqLogStdoutSE="$outPath/$id\_trim_stdout_log.txt";
	
	my $trimComm="java -jar $trimmomatic SE -threads $trimmomatic_threads -phred33 -trimlog $trimmedFastqLogSE $fastq1 $trimmedFastqSE ILLUMINACLIP:$adapters_SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:$trimmomatic_sliding_win MINLEN:$trimmomatic_min_len >> $trimmedFastqLogStdoutSE 2>&1";
    my @exec_trimmomaticSE_stdout;

    if (! -s "$trimmedFastqLogStdoutSE"){
        $logger->info("Executing: Trimmomatic");
        $logger->info("Param: Trimmomatic : adapters = $adapters_SE");
        $logger->info("Param: Trimmomatic : SLIDINGWINDOW = $trimmomatic_sliding_win");
        $logger->info("Param: Trimmomatic : MINLEN = $trimmomatic_min_len");
        $logger->info("Command: $trimComm");

        @exec_trimmomaticSE_stdout=`java -jar $trimmomatic SE -threads $trimmomatic_threads -phred33 -trimlog $trimmedFastqLogSE $fastq1 $trimmedFastqSE ILLUMINACLIP:$adapters_SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:$trimmomatic_sliding_win MINLEN:$trimmomatic_min_len 2>&1`;

        open (TRIMLOG, ">$trimmedFastqLogStdoutSE");
        print TRIMLOG @exec_trimmomaticSE_stdout;
        close TRIMLOG;

        $logger->info("Trimmomatic completed");
		system("rm $trimmedFastqLogSE");
    }
    else{

        $logger->info("Trimmomatic completed");
    }

	if (-s "$trimmedFastqSE"){

		if ($fastqc_exec_trimmed_data_flag == 1 ) {
            
            my $qc_output_dir = "$analysis_dir/$id/fastQC";
            system("mkdir -p $qc_output_dir");
            my $qc_output_summary = "$analysis_dir/$id/fastQC/$id\_AT_fastqc/summary.txt";
    		my @fastFilesList= ($trimmedFastqSE);

            if (! -s "$qc_output_summary"){

                $logger->info("Executing: FastQC on adapter trimmed fastq file");
                run_qc($id, $qc_output_dir, @fastFilesList);
                $logger->info("FastQC on adapter trimmed fastq file completed");
            }
            else{

                $logger->info("FastQC on adapter trimmed fastq file completed");
            }

            my @trimLogArr  = `cat $trimmedFastqLogStdoutSE`;
            chomp $trimLogArr[-2];
		    my @trimLogStat = split(/\s/,$trimLogArr[-2]); #Input Read Pairs: 1467341 Both Surviving: 772840 (52.67%) Forward Only Surviving: 92309 (6.29%) Reverse Only Surviving: 105347 (7.18%) Dropped: 496845 (33.86%)
			$trimLogStat[-1]=~ s/[()%]//g;

			if ($trimLogStat[-1] >= $trimmomatic_read_drop_cutoff){  # if the trimmomatic read drop percentage criteria is not satisfied 

				$logger->error("Trimmomatic QC failed for $id");
				$logger->error("Sample: $id $trimLogArr[-2]");
				exit;
			}
			else { # open each fastqc summary file and look for "PASS" flag for the given qc criterias

                $logger->info("Trimmomatic : Result: Dropped reads $trimLogStat[-1]%.");
                $logger->info("$id passed trimmomatic minimum read threshold");
				$logger->info("Checking sample quality of $id ");
				checkSampleQuality($id,$outPath,$qc_output_dir);
			}
		}
	}

	else { # if adapter trimmed fastq file is empty

		$logger->error("Trimmomatic failed for $sra_sample_id");
		exit;
	}
}

sub run_qc {

	my ($id,$outPath,@fastqArr)= @_;
	
	my $fastqFiles = join(" ", @fastqArr);
    $logger->info("Param: FastQC : threads = $fastqc_threads");
    my $qcComm="$fastqc $fastqFiles --threads $fastqc_threads --quiet --extract --outdir $outPath";
    $logger->info("Command: $qcComm");
	system("$qcComm");
}

sub download_sample {

	my ($id,$samType,$outPath)= @_;
	my $dumpComm;
	my $sumStat=0;

	my $dumpCommPELogStdout = "$outPath/fastq_dump.stdout.txt";
	my $dumpCommSELogStdout = "$outPath/fastq_dump.stdout.txt";

    my $dumpCommPE="$fastq_dump $id --split-files -O $outPath/ >> $dumpCommPELogStdout 2>&1";
    my $dumpCommSE="$fastq_dump $id -O $outPath/ >> $dumpCommSELogStdout 2>&1";



    if ($fastq_dump_perc){
        my $fraction=int (($fastq_dump_perc/100)*$num_spots);
        $logger->info("Fastq dump : Number of spots to be downloaded : $fraction");
        $dumpCommPE="$fastq_dump -X $fraction $id --split-files -O $outPath/ >> $dumpCommPELogStdout 2>&1";
        $dumpCommSE="$fastq_dump -X $fraction $id -O $outPath/ >> $dumpCommSELogStdout 2>&1";
    }

	if ($samType eq "PAIRED"){

        my $fastqn1="$outPath/$id\_1.fastq";
		my $fastqn2="$outPath/$id\_2.fastq";
		my $fastqn3="$outPath/$id\_3.fastq";

		my $flagR1=0; my $flagR2=0; my $flagR3=0;	
			
		for (my $attempt = 1; $attempt<=3; $attempt++){
			
			if(-s "$fastqn1"){$flagR1=1;}
			if(-s "$fastqn2"){$flagR2=1;}
			if(-s "$fastqn3"){$flagR3=1;}
			
			$sumStat=$flagR1 + $flagR2 + $flagR3;

			if ( $sumStat >=2 ){

                $logger->info("Fastq Dump: Completed");
				last;
			}

			else {
				
				$sumStat=0;
                $logger->info("Executing: Fastq dump");
				$logger->info("Command: $dumpCommPE : attempt: $attempt");
				# system("$dumpCommPE");
			}
		}
	}
	
	if ($samType eq "SINGLE") { # for SE samples

		my $flagR1=0;
		my $fastqn1="$outPath/$id.fastq";

		for (my $attempt = 1; $attempt<=3; $attempt++){

			if(-s "$fastqn1"){$flagR1=1;}

			$sumStat=$flagR1;

			if ( $sumStat == 1 ){

                $logger->info("Fastq Dump: Completed");
				last;
			}

			else {

                $logger->info("Executing: Fastq dump");
				$logger->info("Command: $dumpCommSE : attempt: $attempt");
				# system("$dumpCommSE");
			}
		}
	}
    # $logger->info("Fastq Dump: Completed");
	return ($sumStat);
}



sub download_sample_full_sra {

	my ($id,$samType,$outPath)= @_;
	
	my $dumpComm;
	my $sumStat=0;
	my $sraWgetFlag=0;
	
	my $wgetLogStdout = "$outPath/wget.full.sra.stdout.txt";
	my $sra_file="-";
	
	for (my $attempt = 1; $attempt<=3; $attempt++){
		
		if (my @files = glob("$outPath/sra/$id*")) {
			system ("cp $files[0] $outPath/$id");
			system("rm $outPath/sra/*");
			$sra_file = "$outPath/$id";
			if(-s "$sra_file"){
				last;
			}
			
		}	
		
		else{
			$logger->info("Retrieving the complete SRA file: $sra_url : attempt: $attempt");
			my $get_sra_comm= `wget -o $wgetLogStdout -c $sra_url -P $outPath/sra/`;
		}
	}
	
	if(-s "$sra_file"){
		$logger->info("Successfully retieved sra file from $sra_url");
	}
	
	else {
		$logger->error("Failed to retieve sra file from $sra_url");
		return 0;
	}
	
	my $dumpCommPELogStdout = "$outPath/fastq_dump.stdout.txt";
	my $dumpCommSELogStdout = "$outPath/fastq_dump.stdout.txt";
	
	my $dumpCommPE="$fastq_dump $sra_file --split-files -O $outPath/ >> $dumpCommPELogStdout 2>&1";
    my $dumpCommSE="$fastq_dump $sra_file -O $outPath/ >> $dumpCommSELogStdout 2>&1";
	
	if ($fastq_dump_perc){
		my $fraction=int (($fastq_dump_perc/100)*$num_spots);
        $logger->info("Fastq dump : Number of spots to be extracted : $fraction");
        $dumpCommPE="$fastq_dump -X $fraction $sra_file --split-files -O $outPath/ >> $dumpCommPELogStdout 2>&1";
        $dumpCommSE="$fastq_dump -X $fraction $sra_file -O $outPath/ >> $dumpCommSELogStdout 2>&1";
    }
	
	if ($samType eq "PAIRED"){

        my $fastqn1="$outPath/$id\_1.fastq";
		my $fastqn2="$outPath/$id\_2.fastq";
		my $fastqn3="$outPath/$id\_3.fastq";
		
		if((! -s "$fastqn1") or (! -s "$fastqn2")){
			$logger->info("Executing: Fastq dump");
			$logger->info("Command: $dumpCommPE");
			system("$dumpCommPE");	
		}
		
		my $flagR1=0; my $flagR2=0; my $flagR3=0;	
		
		if(-s "$fastqn1"){$flagR1=1;}
		if(-s "$fastqn2"){$flagR2=1;}
		if(-s "$fastqn3"){$flagR3=1;}
		
		$sumStat=$flagR1 + $flagR2 + $flagR3;
		
		if ( $sumStat >=2 ){
			$logger->info("Fastq Dump: Completed");
		}
	}
	
	if ($samType eq "SINGLE") { # for SE samples

		my $fastqn1="$outPath/$id.fastq";

		if(! -s "$fastqn1"){
			$logger->info("Executing: Fastq dump");
			$logger->info("Command: $dumpCommSE");
			system("$dumpCommSE");
		}
		
		if (-s "$fastqn1") {
			$logger->info("Fastq Dump: Completed");
			$sumStat=1;
		}
	}
	
	if (-s "$sra_file"){
		system("rm $sra_file");
	}

	return ($sumStat);
}

sub checkSampleQuality {

	my ($id,$concatenated_fastq_dir,$outPath)= @_;
	
	if ($sample_type eq "PAIRED"){

		my $trimmedFastqPE1="$outPath/$id\_1_AT_fastqc/summary.txt";
		my $trimmedFastqPE2="$outPath/$id\_2_AT_fastqc/summary.txt";
		my $trimmedFastqUnPaired1="$outPath/$id\_1_unpaired_fastqc/summary.txt";
        my $trimmedFastqUnPaired2="$outPath/$id\_2_unpaired_fastqc/summary.txt";
        my $concatenatedPEfastqAT= "$concatenated_fastq_dir/$id\_COM_AT.fastq";
		
        my @summaryFileSet = ($trimmedFastqPE1,$trimmedFastqPE2,$trimmedFastqUnPaired1,$trimmedFastqUnPaired2);
        my $sampleQCFlag=0;

		foreach my $file (@summaryFileSet){

			$sampleQCFlag+=parseQCsummary($trimmedFastqPE1);
		}

		if ($sampleQCFlag == 4){

			$logger->info("$id passed all quality criteria");
		}
        else {

			$logger->error("$id failed quality check");
			exit;
		}
	}

	if ($sample_type eq "SINGLE") {

		my $trimmedFastqSE1="$outPath/$id\_AT_fastqc/summary.txt";
		my $sampleQCFlag=0;

		# $logger->info("Checking fastqc summary $trimmedFastqSE1");
		$sampleQCFlag=parseQCsummary($trimmedFastqSE1);

		if ($sampleQCFlag == 1){

			$logger->info("$id passed all quality criteria");
		}
		else {

			$logger->error("$id failed quality check");
			exit;
		}
	}
}

sub parseQCsummary {

	my ($summaryFile) = (@_);

	my $qual=0; my $adapter= 0;

	open(TXT,"$summaryFile");

	while(my $rec=<TXT>){

		chomp $rec;
		my ($status,$category,$fileName)= split(/\t/,$rec);
		
		if (($category eq 'Per base sequence quality') && ($status eq "PASS")){

			$qual=1;
		}

		if (($category eq 'Adapter Content') && ($status eq "PASS")){

			$adapter=1;
		}
	}

	close TXT;

	if ($qual+$adapter == 2){

		return 1;
	}

	else {

		return 0;
	}
}

sub fastq_to_fasta {

    my ($id,$inpFastq,$outPath)= @_;

    my $trimmedCollapsedTemp="$outPath/$id\_AT_COL_TEMP.fa";
	my $trimmedCollapsedReHead="$outPath/$id\_AT_COL.fa";
	my $trimmedCollapsedReHeadTable="$outPath/$id\_AT_COL_tab.txt";

	if (! -s "$trimmedCollapsedReHead"){

		my $inpFastqCount=`cat $inpFastq|wc -l`;
		chomp $inpFastqCount;
		my $inpReads= $inpFastqCount/4;

		my $fastqCollapseWithCount="$fastx_collapser -i $inpFastq -o $trimmedCollapsedTemp";

		$logger->info("Executing: fastx_collapser on $inpFastq");
		$logger->info("Command: $fastqCollapseWithCount");
        system("$fastqCollapseWithCount");

		open (FA,"$trimmedCollapsedTemp") or die $logger->error("Cant open $trimmedCollapsedTemp");
		open (FAOUT,">$trimmedCollapsedReHead") or die $logger->error("Cant write to $trimmedCollapsedReHead");

		my $uniqID="0000001";

		while (my $rec=<FA>) {

			chomp $rec;
			if ($rec=~ /^\>/){

				my ($sn,$count)= split (/\-/,$rec);
				$sn=~ s/\>//g;
				my $normalizedCount=sprintf "%.2f",(($count/$inpReads)*10000000);
				print FAOUT ">$id\_$uniqID\_$count\_$normalizedCount\n";
				$uniqID++;
			}
			else {

				print FAOUT "$rec\n";
			}
		}
		close FA;
		close FAOUT;

		system("rm $trimmedCollapsedTemp");

		$logger->info("Collapsed fasta: $trimmedCollapsedReHead");
		return "$trimmedCollapsedReHead";
	}
	else {

		$logger->info("Collapsed fasta: $trimmedCollapsedReHead");
		return "$trimmedCollapsedReHead";
	}
}

sub check_blast_hits {

    my ($id,$blastInpQuery,$blastFilteredResult)= @_;
    
    my $uniqueReads=0; my $totalReads=0; my $blastHitscount=0; my $countMappedReads=0; my $percBlastHits=0; my $countMappedReadsPerc=0;
    
    # get number of unique sequences from query.fasta
    my $inpSeqCount=`cat $blastInpQuery|wc -l`;
    chomp $inpSeqCount;
    $uniqueReads= $inpSeqCount/2;

    # get the expanded reads using the serial number from query.fasta
    my $allReadHeaders=`cat $blastInpQuery|grep ">"`;
    chomp $allReadHeaders;
    my @allReadIDs = split (/\n/,$allReadHeaders);
    foreach my $read (@allReadIDs){
        my ($run,$sNo,$count,$normCount) = split(/\_/,$read);
        $totalReads+=$count;
    }

    # get number of unique query IDs in the blast filtered table
    my $queryHitIDs=`cut -f1 $blastFilteredResult|sort|uniq`;
    chomp $queryHitIDs;
    my @readIDs = split (/\n/,$queryHitIDs);
    $blastHitscount = scalar @readIDs;

    foreach my $read (@readIDs){
        my ($run,$sNo,$count,$normCount) = split(/\_/,$read);
        $countMappedReads+=$count;
    }

    $percBlastHits=sprintf "%.4f",(($blastHitscount/$uniqueReads)*100);
    $countMappedReadsPerc=sprintf "%.4f",(($countMappedReads/$totalReads)*100);
    my $retStr="$uniqueReads\t$blastHitscount\t$percBlastHits\t$totalReads\t$countMappedReads\t$countMappedReadsPerc";
    
    return ($percBlastHits,$countMappedReadsPerc,$retStr);
}



sub check_bowtie2_hits {
    my ($bowtie_analysis_dir,$aligned_bam,$trimmedCollapsedReHeadFasta)= @_;
    
    my $uniqueReads=0; my $totalReads=0; my $bowtieHitscount=0; my $countMappedReads=0; my $percBowtieHits=0; my $countMappedReadsPerc=0;

    # get the expanded reads using the serial number from query.fasta
    my $allReadHeaders=`cat $trimmedCollapsedReHeadFasta|grep ">"`;
    chomp $allReadHeaders;
    my @allReadIDs = split (/\n/,$allReadHeaders);
    foreach my $read (@allReadIDs){
        my ($run,$sNo,$count,$normCount) = split(/\_/,$read);
        $totalReads+=$count;
    }

    $uniqueReads= scalar @allReadIDs;
    
    # get the number of mapped IDs in the bam file
    my $mappedIDsList = `$samtools view $aligned_bam |cut -f1|uniq`;
    chomp $mappedIDsList;
    my @mappedReadIDs = split (/\n/,$mappedIDsList);
    
    $bowtieHitscount = scalar @mappedReadIDs;

    foreach my $read (@mappedReadIDs){
        my ($run,$sNo,$count,$normCount) = split(/\_/,$read);
        $countMappedReads+=$count;
    }

    $percBowtieHits=sprintf "%.4f",(($bowtieHitscount/$uniqueReads)*100);
    $countMappedReadsPerc=sprintf "%.4f",(($countMappedReads/$totalReads)*100);
    my $retStr="$uniqueReads\t$bowtieHitscount\t$percBowtieHits\t$totalReads\t$countMappedReads\t$countMappedReadsPerc";
    
    return ($percBowtieHits,$countMappedReadsPerc,$retStr);
}


sub create_bowtie2_index {
	my ($inp_sequence,$outDir)= @_;
	
	my @fasta_label = split (/\./, basename($inp_sequence));
	pop @fasta_label;
	my $index_label = join(".",@fasta_label);

	system("mkdir -p $outDir/bowtie2_index");

	my $create_index_comm = "$bowtie2_build $inp_sequence $outDir/bowtie2_index/$index_label";
	if (! glob ("$outDir/bowtie2_index/$index_label*")){
		system ("$create_index_comm");
	}

	return "$outDir/bowtie2_index/$index_label";
    
}


sub create_blastn_db {
	my ($inp_sequence,$outDir)= @_;
	
	my @fasta_label = split (/\./, basename($inp_sequence));
	pop @fasta_label;
	my $index_label = join(".",@fasta_label);

	system("mkdir -p $outDir/blastn_db");

	# my $ = "$bowtie2_build $inp_sequence $outDir/blastn_db/$index_label";
	my $create_db_comm = "$makeblastdb -in $inp_sequence -dbtype 'nucl' -hash_index -out $outDir/blastn_db/$index_label";
	
	if (! glob ("$outDir/blastn_db/$index_label*")){
		system ("$create_db_comm");
	}

	return "$outDir/blastn_db/$index_label";
    
}