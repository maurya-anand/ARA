#!/usr/bin/env perl


################ Pipeline to retrieve sample attributes for every SRA run
################ to categorize the samples into tissue, genotype and treatment group.
################ to align the reads to the provided reference.

use strict;
use warnings;
use Data::Dumper;
use Config::Simple;
use Parallel::ForkManager;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case);
use FindBin qw($Bin);
use Sys::Hostname;
use File::Basename;

my $host = hostname;
my $base_path = $Bin;

my $input_run_info = '';
my $input_seq_fasta = '';
my $output_directory = '';
my $output_file_label = '';
my $run_mode = 'screen';
my $input_config = '';
my $help = 0;

my $default_out_dir = "$base_path/results";
my $default_config_file = "$base_path/conf.txt";

# my $input_config = $ARGV[0];

my ($maxproc) = 1;
my $pm = Parallel::ForkManager->new($maxproc);

my $tool_config_obj;
my $exec_SRA_screen_flag = 1;
my $exec_SRA_analysis_flag = 1;

## Initialize Log4perl, if it is not initialized already by the calling script
unless (Log::Log4perl::initialized()) {
	Log::Log4perl::init("$base_path/src/main/resources/log4perl.config");
}
my $logger = Log::Log4perl->get_logger();

GetOptions ("i|input=s" => \$input_run_info,
			"s|sequences=s" => \$input_seq_fasta,
			"o|output=s"   => \$output_directory,
			"m|mode=s"  => \$run_mode,
			"c|config=s"  => \$input_config,
			"h|help"  => \$help
			) 
			or
			$logger->error("Invalid arguments");

if ($help) {
	$logger->info("\n\nUsage info :\n\nperl ama.pl\n--input or -i  <runInfo_csv>\n--sequences or -s  <sequences.fasta>\n--output or -o  <output_directory>\n--mode or -m  <run_mode> (eg: 'screen' or 'full' or 'both' or 'summary'). Default: screen\n--config or -c  <pipeline_config_file>\n");
	$logger->info("Example usage:\nperl ama.pl --input example/SraRunInfo.csv --sequences example/Arabidopsis_thaliana.TAIR10.ncrna.fa --output src/main/test/ --mode screen --config conf.txt");
	exit;
}
else {
	my $inp_flag = 0;

	if ($input_run_info eq ''){
		$inp_flag = 1; $logger->error("Invalid input. Please provide a valid 'SRA run accession' or a file containing the list of SRA run accessions.");
	}
	
	if ((! -s "$input_seq_fasta") or ($input_seq_fasta eq '')){
		$inp_flag = 1; $logger->error("Invalid fasta file. Please provide the full path of the file.");
	}
	
	if (! -d "$output_directory"){
		$logger->info("Output directory not specified. Using defaults.");
		# $inp_flag = 1; $logger->error("Invalid Output directory.");}
		$output_directory = $default_out_dir;
	}

	if ($run_mode !~ /^screen$|^full$|^summary$|^both$/i){
		$logger->info("Run mode not specified. Using default 'screen' mode. (Choices: 'screen' or 'full' or 'both' or 'summary')");
	}

	if (! -s "$input_config"){
		if (-s "$default_config_file"){
			$logger->info("Configuration file not specified: $input_config. Using default $default_config_file.");
			$input_config = $default_config_file;
		}
		else {
			$inp_flag = 1; $logger->error("Configuration file not found.");
			$logger->info("Please Provide a valid configuration file or run setup.pl to generate one.");
			exit;
		}
	}

	if ($inp_flag == 1){
		$logger->info("Type this command for usage info: perl ama.pl -help");
		exit;
	}
	else {
		## Loading the pipeline configuration
		# $logger->info("Pipeline : Config: $input_config");
		$tool_config_obj = new Config::Simple($input_config);
	}
}

# test the dependencies
my $exit_code = system("perl $base_path/src/main/scripts/utils/check.config.pl $input_config");
if ($exit_code != 0){
	exit;
}

my $fastq_dump_perc		= $tool_config_obj->param('fastq_dump.data_perc');
my $blastn_exec_flag	= $tool_config_obj->param('blastn.execute');
my $bowtie2_exec_flag	= $tool_config_obj->param('bowtie2.execute');
my $ncbi_edirect_path	= $tool_config_obj->param('toolPath.ncbi_edirect');
my $pipeline_threads	= $tool_config_obj->param('pipeline.threads');
my $pipeline_ver		= $tool_config_obj->param('pipeline.version');

my $bowtie2_build	= $tool_config_obj->param('toolPath.bowtie2-build');
my $makeblastdb		= $tool_config_obj->param('toolPath.makeblastdb');

my $single_run_mode = 0;
my $multi_run_mode = 0;

if(! -s "$input_run_info"){
	if (validate_input($input_run_info) == 0){
		$logger->error("Invalid input. Please provide a valid SRA run accession.");
		exit;
	}
	else {
		$single_run_mode = 1;
	}
}
else {
	$multi_run_mode = 1;
}

if (int($pipeline_threads) > 1) {
	$pm = Parallel::ForkManager->new($pipeline_threads);
}
else {
	$pm = Parallel::ForkManager->new($maxproc);
}


if (lc $run_mode eq "screen"){
	$exec_SRA_analysis_flag = 0;
}
if (lc $run_mode eq "full"){
	$exec_SRA_screen_flag = 0;
}
if (lc $run_mode eq "both"){
	$exec_SRA_screen_flag = 1;
	$exec_SRA_analysis_flag = 1;
}
if (lc $run_mode eq "summary"){
	$exec_SRA_screen_flag = 0;
	$exec_SRA_analysis_flag = 0;
}

my @fasta_label = split (/\./, basename($input_seq_fasta));
pop @fasta_label;
$output_file_label = join(".",@fasta_label);
# print "$output_file_label\n";

my $metadata_dir = "$output_directory/$output_file_label/metadata";
system("mkdir -p $metadata_dir");

my $references_dir = "$output_directory/$output_file_label/reference";

# retrieving run and sample info from NCBI SRA
$logger->info("Pipeline : AMA v.$pipeline_ver executed on $host");
$logger->info("Using the following parameters:\n\n--input: $input_run_info\n--sequence: $input_seq_fasta\n--output: $output_directory\n--mode: $run_mode\n--config: $input_config\n\n");

sub validate_input {
	my ($id)=(@_);
	my $runAccn = `$ncbi_edirect_path/esearch -db sra -query $id | $ncbi_edirect_path/efetch -format xml | $ncbi_edirect_path/xtract -pattern EXPERIMENT_PACKAGE -block RUN_SET -element RUN\@accession`;
	chomp $runAccn;
	if ($id eq $runAccn){
		return 1;
	}
	else {
		return 0;
	}
}

if ($exec_SRA_screen_flag == 1){
	my $data_analysis_dir = "$output_directory/$output_file_label/screening_results";
	my $outFileParsed="$metadata_dir/$output_file_label.metadata.txt";
	my $sraIDsBlastStats="$output_directory/$output_file_label/$output_file_label.screening.analysis.stats.txt";
	my $outRunInfoParsedScreened="$metadata_dir/$output_file_label.metadata.screened.txt";
	
	if ($single_run_mode == 1){
		my $metadata_exec_comm = "perl $base_path/src/main/scripts/metadata/fetch.SRA.metadata.pl $input_run_info $metadata_dir $output_file_label $base_path $ncbi_edirect_path single";
		system("$metadata_exec_comm");
	}
	else {
		my $metadata_exec_comm = "perl $base_path/src/main/scripts/metadata/fetch.SRA.metadata.pl $input_run_info $metadata_dir $output_file_label $base_path $ncbi_edirect_path list";
		system("$metadata_exec_comm");
	}
	
	system("mkdir -p $references_dir");
	if ($blastn_exec_flag == 1){
		my $blastn_db_path =  create_blastn_db($input_seq_fasta,$references_dir);
	}
	if ($bowtie2_exec_flag == 1){
		my $bowtie2_index_path =  create_bowtie2_index($input_seq_fasta,$references_dir);
	}

	system("mkdir -p $data_analysis_dir");
	screen_raw_data($outFileParsed,$data_analysis_dir,$sraIDsBlastStats,$outRunInfoParsedScreened);
}

if ($exec_SRA_analysis_flag == 1){
	my $data_analysis_dir = "$output_directory/$output_file_label/full_analyis_results";
	my $runInfoParsed="$metadata_dir/$output_file_label.metadata.txt";
	my $sraIDsAlnStats="$output_directory/$output_file_label/$output_file_label.full.analysis.stats.txt";

	if ($single_run_mode == 1){
		my $metadata_exec_comm = "perl $base_path/src/main/scripts/metadata/fetch.SRA.metadata.pl $input_run_info $metadata_dir $output_file_label $base_path $ncbi_edirect_path single";
		system("$metadata_exec_comm");
	}
	else {
		my $metadata_exec_comm = "perl $base_path/src/main/scripts/metadata/fetch.SRA.metadata.pl $input_run_info $metadata_dir $output_file_label $base_path $ncbi_edirect_path list";
		system("$metadata_exec_comm");
	}
	
	if (lc $run_mode eq "both"){
		$runInfoParsed="$metadata_dir/$output_file_label.metadata.screened.txt";
	}
	
	system("mkdir -p $references_dir");
	if ($blastn_exec_flag == 1){
		my $blastn_db_path =  create_blastn_db($input_seq_fasta,$references_dir);
	}
	if ($bowtie2_exec_flag == 1){
		my $bowtie2_index_path =  create_bowtie2_index($input_seq_fasta,$references_dir);
	}

	system("mkdir -p $data_analysis_dir");
	analyse_raw_data($runInfoParsed,$data_analysis_dir,$sraIDsAlnStats);
}


if (lc $run_mode eq "summary"){
	my $screening_analysis_dir = "$output_directory/$output_file_label/screening_results";
	my $full_analysis_dir = "$output_directory/$output_file_label/full_analyis_results";
	
	my $runInfoParsed="$metadata_dir/$output_file_label.metadata.txt";

	if (-d "$screening_analysis_dir"){
		my $screenAlnStats="$output_directory/$output_file_label/$output_file_label.screening.analysis.stats.txt";
		summarize($runInfoParsed,$screening_analysis_dir,$screenAlnStats);
		$logger->info("Pipeline : Combined screening alignment stats : $screenAlnStats");
	}
	if (-d "$full_analysis_dir"){
		my $fullAlnStats="$output_directory/$output_file_label/$output_file_label.full.analysis.stats.txt";
		summarize($runInfoParsed,$full_analysis_dir,$fullAlnStats);
		$logger->info("Pipeline : Combined alignment stats : $fullAlnStats");
	}
}

$logger->info("Analysis completed");
# $logger->info("Please check the runlog.*.txt for details.");

## Takes the metadata table (parsed uncategorised) and executes blast pipeline on complete dataset from each run
## Generates a combined blast hits percentage table after analysis
sub analyse_raw_data {

	my($metadata,$work_dir,$statsFile)=(@_);

	my %runIDs;
	
	open (TAB,"$metadata") or $logger->logdie("Cannot read metadata - $metadata: $!");
	my $header = <TAB>; chomp $header;
	my @runInfoHeaders = split("\t",$header);
	my $runIDcolIdx; my $runtypeColIdx; my $runspotsColIdx;
	for(my $i=0;$i<=$#runInfoHeaders;$i++){
		if ($runInfoHeaders[$i] eq 'Run'){$runIDcolIdx = $i;}
		if ($runInfoHeaders[$i] eq 'RunSpots'){$runspotsColIdx = $i;}
		if ($runInfoHeaders[$i] eq 'RunType'){$runtypeColIdx = $i;}
	}
	while(my $rec=<TAB>){
		chomp $rec;
		# if ($rec=~ /^Run/) {@runInfoHeaders = split("\t",$rec);}
		if (($rec=~ /^Run/)|| ($rec eq "")) {next;}
		my @cols = split("\t",$rec);
		$runIDs{$cols[0]}{'id'}=$cols[$runIDcolIdx];
		$runIDs{$cols[0]}{'type'}=$cols[$runtypeColIdx];
		$runIDs{$cols[0]}{'spots'}=$cols[$runspotsColIdx];
		$runIDs{$cols[0]}{'row'}=$rec;
	}
	close TAB;

	if (lc $run_mode ne 'summary'){
		## parallel execution of the pipeline
		foreach my $run (keys %runIDs){
			my $pid = $pm->start and next;
			$logger->info("Analysing $runIDs{$run}{'id'} $runIDs{$run}{'type'}");
			analyse_sra($runIDs{$run}{'id'},$runIDs{$run}{'type'},$runIDs{$run}{'spots'},$work_dir);
			$pm->finish;
		}
		$pm->wait_all_children;
	}

	if (-d "$work_dir"){
		
		open (STATS,">$statsFile") or $logger->logdie("Cannot write blast hits stats to - $statsFile: $!");
		
		my $header;

		if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		elsif ($blastn_exec_flag == 1){
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		else{
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		
		print STATS "$header\n";

		## collect stats for each run from the metadata
		foreach my $runAccn (keys %runIDs){
			my $blast_stats = ''; my $bowtie_stats = '';
			
			$blast_stats = collect_blast_stats($runAccn,$work_dir);
			$bowtie_stats = collect_bowtie2_stats($runAccn,$work_dir);
			
			if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
				print STATS "$runAccn\t$blast_stats\t$bowtie_stats\t$runIDs{$runAccn}{'row'}\n";
			}
			elsif ($blastn_exec_flag == 1){
				print STATS "$runAccn\t$blast_stats\t$runIDs{$runAccn}{'row'}\n";
			}
			else{
				print STATS "$runAccn\t$bowtie_stats\t$runIDs{$runAccn}{'row'}\n";
			}
		}
		close STATS;

		reorder($statsFile,$header);

		$logger->info("Pipeline : Combined alignment stats : $statsFile");
		$logger->info("Pipeline : Summary : More details can be found in sample-wise log files at $work_dir/runlog.*.txt");
	}
}

## Takes the metadata table (parsed uncategorised) and executes blast pipeline for 10% reads from each run
## Generates a combined blast hits percentage table after screening
sub screen_raw_data {

	my($metadata,$work_dir,$statsFile,$screenedRunInfo)=(@_);
	my %runIDs;

	open (TAB,"$metadata") or $logger->logdie("Cannot read metadata - $metadata: $!");
	my $header = <TAB>; chomp $header;
	my @runInfoHeaders = split("\t",$header);
	my $runIDcolIdx; my $runtypeColIdx; my $runspotsColIdx;
	for(my $i=0;$i<=$#runInfoHeaders;$i++){
		if ($runInfoHeaders[$i] eq 'Run'){$runIDcolIdx = $i;}
		if ($runInfoHeaders[$i] eq 'RunSpots'){$runspotsColIdx = $i;}
		if ($runInfoHeaders[$i] eq 'RunType'){$runtypeColIdx = $i;}
	}
	
	while(my $rec=<TAB>){
		chomp $rec;
		if (($rec=~ /^Run/)|| ($rec eq "")) {next;}
		my @cols = split("\t",$rec);
		$runIDs{$cols[0]}{'id'}=$cols[$runIDcolIdx];
		$runIDs{$cols[0]}{'type'}=$cols[$runtypeColIdx];
		$runIDs{$cols[0]}{'spots'}=$cols[$runspotsColIdx];
		$runIDs{$cols[0]}{'row'}=$rec;
		# print "$rec\t$cols[$runIDcolIdx]--$cols[$runtypeColIdx]--$cols[$runspotsColIdx]\n";
	}
	close TAB;

	if (lc $run_mode ne 'summary'){
		## parallel execution of the pipeline
		foreach my $run (keys %runIDs){
			my $pid = $pm->start and next;
			# print "$runIDs{$run}{'id'},$runIDs{$run}{'type'},$runIDs{$run}{'spots'},$work_dir,10\n";
			analyse_sra($runIDs{$run}{'id'},$runIDs{$run}{'type'},$runIDs{$run}{'spots'},$work_dir,$fastq_dump_perc);
			$pm->finish;
		}
		$pm->wait_all_children;
	}

	if (-d "$work_dir"){
		open (SCREENED,">$screenedRunInfo") or $logger->logdie("Cannot write screened run info to - $screenedRunInfo: $!");
		open (STATS,">$statsFile") or $logger->logdie("Cannot write blast screeing stats to - $statsFile: $!");

		my $header;

		if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders);
		}
		elsif ($blastn_exec_flag == 1){
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\t".join("\t",@runInfoHeaders);
		}
		else{
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders);
		}
		
		print STATS "$header\n";
		print SCREENED join("\t",@runInfoHeaders)."\n";

		my $passed_runs=0;
		my $total_runs = scalar %runIDs;

		## collect stats for each run from the metadata
		foreach my $run (keys %runIDs){
			my $blast_stats; my $bowtie_stats;
			
			if (($blastn_exec_flag == 1) and ($bowtie2_exec_flag == 1)){
				$blast_stats = collect_blast_stats($run,$work_dir);
				$bowtie_stats = collect_bowtie2_stats($run,$work_dir);
				
				my @blStatsRec = split("\t",$blast_stats);
				my @bowStatsRec = split("\t",$bowtie_stats);
				
				print STATS "$run\t$blast_stats\t$bowtie_stats\t$runIDs{$run}{'row'}\n";
				
				if (($blStatsRec[-1] eq "PASS") or ($bowStatsRec[-1] eq "PASS")){
					print SCREENED "$runIDs{$run}{'row'}\n";
					$passed_runs+=1;
				}
			}
			elsif ($blastn_exec_flag == 1){
				$blast_stats = collect_blast_stats($run,$work_dir);
				
				my @blStatsRec = split("\t",$blast_stats);

				print STATS "$run\t$blast_stats\t$runIDs{$run}{'row'}\n";
				
				if (($blStatsRec[-1] eq "PASS")){
					print SCREENED "$runIDs{$run}{'row'}\n";
					$passed_runs+=1;
				}
			}
			else{
				$bowtie_stats = collect_bowtie2_stats($run,$work_dir);
				
				my @bowStatsRec = split("\t",$bowtie_stats);
				
				print STATS "$run\t$bowtie_stats\t$runIDs{$run}{'row'}\n";
				
				if (($bowStatsRec[-1] eq "PASS")){
					print SCREENED "$runIDs{$run}{'row'}\n";
					$passed_runs+=1;
				}
			}
		}
		close STATS;
		
		reorder($statsFile,$header);

		$logger->info("Pipeline : RAW data screening stats : $statsFile");
		$logger->info("Pipeline : Summary : Total runs processed: $total_runs. Total number of passed runs: $passed_runs");
		$logger->info("Pipeline : Summary : More details can be found in sample-wise log files at $work_dir/runlog.*.txt");

	}
}

sub analyse_sra {

	my($runID,$type,$spots,$work_dir,$fraction)=(@_);

	if ($fraction){
		$logger->info("Pipeline : Processing : $runID");
		my $analysis_comm="perl $base_path/src/main/scripts/analysis/analyse.SRA.data.pl $input_config $runID $type $spots $work_dir $input_seq_fasta $references_dir $fraction";
		system("$analysis_comm");
	}
	else {
		$logger->info("Pipeline : Processing : $runID");
		my $analysis_comm="perl $base_path/src/main/scripts/analysis/analyse.SRA.data.pl $input_config $runID $type $spots $work_dir $input_seq_fasta $references_dir";
		system("$analysis_comm");
	}
}

sub collect_blast_stats {
	my($runID,$work_dir)=(@_);
	my $stats="-\t-\t-\t-\t-\t-\t-";
	my $blast_stats_file="$work_dir/$runID/blastn/blast.stats.txt";
	if (-s "$blast_stats_file"){
		$stats= `tail -n1 $blast_stats_file`;
		chomp $stats;
	}
	return $stats;
}

sub collect_bowtie2_stats {
	my($runID,$work_dir)=(@_);
	my $overall_aln="-\t-\t-\t-\t-\t-\t-";
	my $bowtie2_stats_file="$work_dir/$runID/bowtie2/alignment.txt";
	if (-s "$bowtie2_stats_file"){
		$overall_aln= `tail -n1 $bowtie2_stats_file`;
		chomp $overall_aln; # 2315245   1068    0.0461  3137916 1269    0.0404
	}
	return $overall_aln;
}

sub summarize {
	my($metadata,$work_dir,$statsFile)=(@_);
	
	my %runIDs;
	
	open (TAB,"$metadata") or $logger->logdie("Cannot read metadata - $metadata: $!");
	my $header = <TAB>; chomp $header;
	my @runInfoHeaders = split("\t",$header);
	my $runIDcolIdx; my $runtypeColIdx; my $runspotsColIdx;
	for(my $i=0;$i<=$#runInfoHeaders;$i++){
		if ($runInfoHeaders[$i] eq 'Run'){$runIDcolIdx = $i;}
		if ($runInfoHeaders[$i] eq 'RunSpots'){$runspotsColIdx = $i;}
		if ($runInfoHeaders[$i] eq 'RunType'){$runtypeColIdx = $i;}
	}
	while(my $rec=<TAB>){
		chomp $rec;
		# if ($rec=~ /^Run/) {@runInfoHeaders = split("\t",$rec);}
		if (($rec=~ /^Run/)|| ($rec eq "")) {next;}
		my @cols = split("\t",$rec);
		$runIDs{$cols[0]}{'id'}=$cols[$runIDcolIdx];
		$runIDs{$cols[0]}{'type'}=$cols[$runtypeColIdx];
		$runIDs{$cols[0]}{'spots'}=$cols[$runspotsColIdx];
		$runIDs{$cols[0]}{'row'}=$rec;
	}
	close TAB;

	if (-d "$work_dir"){
		
		open (STATS,">$statsFile") or $logger->logdie("Cannot write blast hits stats to - $statsFile: $!");

		my $header;

		if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders);
		}
		elsif ($blastn_exec_flag == 1){
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\t".join("\t",@runInfoHeaders);
		}
		else{
			$header = "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders);
		}

		print STATS "$header\n";

		## collect stats for each run from the metadata
		foreach my $runAccn (keys %runIDs){
			my $blast_stats = ''; my $bowtie_stats = '';
			
			$blast_stats = collect_blast_stats($runAccn,$work_dir);
			$bowtie_stats = collect_bowtie2_stats($runAccn,$work_dir);
			
			if (($blastn_exec_flag == 1) and ($bowtie2_exec_flag == 1)){
				print STATS "$runAccn\t$blast_stats\t$bowtie_stats\t$runIDs{$runAccn}{'row'}\n";
			}
			elsif ($blastn_exec_flag == 1){
				print STATS "$runAccn\t$blast_stats\t$runIDs{$runAccn}{'row'}\n";
			}
			else{
				print STATS "$runAccn\t$bowtie_stats\t$runIDs{$runAccn}{'row'}\n";
			}
		}
		close STATS;
		reorder($statsFile,$header);
	}
}

sub create_bowtie2_index {
	my ($inp_sequence,$outDir)= @_;
	
	my @fasta_label = split (/\./, basename($inp_sequence));
	pop @fasta_label;
	my $index_label = join(".",@fasta_label);

	system("mkdir -p $outDir/bowtie2_index");

	my $bowtie2_indexLogStdout = "$outDir/bowtie2_index.stdout.txt";

	my $create_index_comm = "$bowtie2_build $inp_sequence $outDir/bowtie2_index/$index_label >> $bowtie2_indexLogStdout 2>&1";
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

	my $blastn_dbLogStdout = "$outDir/makeblastdb.stdout.txt";

	# my $ = "$bowtie2_build $inp_sequence $outDir/blastn_db/$index_label";
	my $create_db_comm = "$makeblastdb -in $inp_sequence -dbtype 'nucl' -hash_index -out $outDir/blastn_db/$index_label >> $blastn_dbLogStdout 2>&1";
	
	if (! glob ("$outDir/blastn_db/$index_label*")){
		system ("$create_db_comm");
	}

	return "$outDir/blastn_db/$index_label";
    
}

sub reorder{
	my($statsFile,$origHeaders)=(@_);
	my $comm= `head -n 1 $statsFile && tail -n +2 $statsFile | sort -k7,7nr`;
	# print "$origHeaders\n";
	# print "---comm\n$comm\n---comm\n";
	open (STATS,">$statsFile") or $logger->logdie("Cannot write blast hits stats to - $statsFile: $!");
	print STATS $comm;
	close STATS;
}