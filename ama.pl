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

my $input_run_info = '';
my $output_directory = '';
my $output_file_label = '';
my $run_mode = '';
my $input_config = '';
my $help = 0;

# my $input_config = $ARGV[0];

my $host = hostname;
my $base_path = $Bin;

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
			"o|output=s"   => \$output_directory,
			"l|label=s"  => \$output_file_label,
			"m|mode=s"  => \$run_mode,
			"c|config=s"  => \$input_config,
			"h|help"  => \$help
			) 
			or
			$logger->error("Invalid arguments");

if ($help) {
	$logger->info("\n\nUsage info :\nperl ama.pl\n--input or -i\t<runInfo_csv>\n--output or -o\t<output_directory>\n--label or -l\t<prefix_label>\n--mode or -m\t<run_mode> (eg: 'screen' or 'full' or 'both' or 'summary')\n--config or -c\t<pipeline_config_file>\n");
	$logger->info("Example usage:\nperl ama.pl --input SraRunInfo.csv --output src/main/test/ --label test --mode screen --config conf.txt");
	exit;
}
else {
	my $inp_flag = 0;
	# if ($input_run_info eq ''){$inp_flag = 1; $logger->error("Run info table is missing.");}
	# if ((! -s "$input_run_info") or ($input_run_info eq '')){$inp_flag = 1; $logger->error("Invalid input.");}
	if ($input_run_info eq ''){
		$inp_flag = 1; $logger->error("Invalid input.");
	}
	# if ($output_directory eq ''){$inp_flag = 1; $logger->error("Output directory is missing.");}
	if (! -d "$output_directory"){$inp_flag = 1; $logger->error("Invalid Output directory.");}
	
	if ($output_file_label eq ''){$inp_flag = 1; $logger->error("Analysis lable/prefix is missing.");}
	
	# if ($run_mode eq ''){$inp_flag = 1; $logger->error("Pipeline run mode is missing. (use: 'screen' or 'full' or 'both')");}
	if ($run_mode !~ /^screen$|^full$|^summary$|^both$/i){$inp_flag = 1; $logger->error("Invalid run mode. (use: 'screen' or 'full' or 'both' or 'summary')");}
	
	# if ($input_config eq ''){$inp_flag = 1; $logger->error("Pipeline configuration file is missing.");}
	if (! -s "$input_config"){$inp_flag = 1; $logger->error("Invalid configuration file");}

	if ($inp_flag == 1){
		$logger->info("Try: perl ama.pl -help");
		exit;
	}
	else {
		## Loading the pipeline configuration
		$logger->info("Pipeline : Config: $input_config");
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

my $blastn_db_path		= $tool_config_obj->param('blastn.db_path');
my $bowtie2_index_path	= $tool_config_obj->param('bowtie2.bowtie_index_path');

my $single_run_mode = 0;
my $multi_run_mode = 0;

if ($blastn_exec_flag == 1){
	if (! glob ("$blastn_db_path*")){
		$logger->error("Incorrect Configuration : blastn > db_path : Please provide the path for the blast database.");
		exit;
	}
}

if ($bowtie2_exec_flag == 1){
	if (! glob ("$bowtie2_index_path*")){
		$logger->error("Incorrect Configuration : bowtie2 > bowtie_index_path : Please provide the path for the bowtie2 index.");
		exit;
	}
}

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

$logger->info("Pipeline : v.$pipeline_ver executed on $host");
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

my $metadata_dir = "$output_directory/$output_file_label/metadata";
system("mkdir -p $metadata_dir");

# retrieving run and sample info from NCBI SRA
$logger->info("USING: $input_run_info");


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

my $exec_comm1 = "perl $base_path/src/main/scripts/metadata/fetch.SRA.metadata.pl $input_run_info $metadata_dir $output_file_label $base_path $ncbi_edirect_path list";
if ($single_run_mode == 1){
	$exec_comm1 = "perl $base_path/src/main/scripts/metadata/fetch.SRA.metadata.pl $input_run_info $metadata_dir $output_file_label $base_path $ncbi_edirect_path single";
}
# $logger->info("COMMAND: $exec_comm1");
system("$exec_comm1");

if ($exec_SRA_screen_flag == 1){
	
	my $data_analysis_dir = "$output_directory/$output_file_label/screening_results";
	
	my $outFileParsed="$metadata_dir/$output_file_label.metadata.txt";
	my $sraIDsBlastStats="$output_directory/$output_file_label/$output_file_label.screening.analysis.stats.txt";
	my $outRunInfoParsedScreened="$metadata_dir/$output_file_label.metadata.screened.txt";
	
	if (-s "$outFileParsed"){
		system("mkdir -p $data_analysis_dir");
		screen_raw_data($outFileParsed,$data_analysis_dir,$sraIDsBlastStats,$outRunInfoParsedScreened);
	}
}

if ($exec_SRA_analysis_flag == 1){

	my $data_analysis_dir = "$output_directory/$output_file_label/full_analyis_results";

	my $runInfoParsed="$metadata_dir/$output_file_label.metadata.txt";
	
	if (lc $run_mode eq "both"){
		$runInfoParsed="$metadata_dir/$output_file_label.metadata.screened.txt";
	}
	
	my $sraIDsAlnStats="$output_directory/$output_file_label/$output_file_label.full.analysis.stats.txt";
	
	if (-s "$runInfoParsed"){
		system("mkdir -p $data_analysis_dir");
		analyse_raw_data($runInfoParsed,$data_analysis_dir,$sraIDsAlnStats);
		# $logger->info("Pipeline : Combined alignment stats : $sraIDsAlnStats");
	}
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

## Takes the metadata table (parsed uncategorised) and executes blast pipeline on complete dataset from each run
## Generates a combined blast hits percentage table after analysis
sub analyse_raw_data {

	my($metadata,$work_dir,$statsFile)=(@_);

	my %runIDs;
	# my @runInfoHeaders;
	
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

		if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		elsif ($blastn_exec_flag == 1){
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		else{
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}

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
		$logger->info("Pipeline : Combined alignment stats : $statsFile");
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
	
		if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		elsif ($blastn_exec_flag == 1){
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		else{
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		
		print SCREENED join("\t",@runInfoHeaders)."\n";

		## collect stats for each run from the metadata
		foreach my $run (keys %runIDs){
			my $blast_stats; my $bowtie_stats;
			
			$blast_stats = collect_blast_stats($run,$work_dir);
			$bowtie_stats = collect_bowtie2_stats($run,$work_dir);
			
			my @blStatsRec = split("\t",$blast_stats);
			my @bowStatsRec = split("\t",$bowtie_stats);
			
			if (($blastn_exec_flag == 1) or ($bowtie2_exec_flag == 1)){
				print STATS "$run\t$blast_stats\t$bowtie_stats\t$runIDs{$run}{'row'}\n";
				
				if (($blStatsRec[-1] eq "PASS") or ($bowStatsRec[-1] eq "PASS")){
					print SCREENED "$runIDs{$run}{'row'}\n";
				}
			}
			elsif ($blastn_exec_flag == 1){
				print STATS "$run\t$blast_stats\t$runIDs{$run}{'row'}\n";
				
				if (($blStatsRec[-1] eq "PASS")){
					print SCREENED "$runIDs{$run}{'row'}\n";
				}
			}
			else{
				print STATS "$run\t$bowtie_stats\t$runIDs{$run}{'row'}\n";
				
				if (($bowStatsRec[-1] eq "PASS")){
					print SCREENED "$runIDs{$run}{'row'}\n";
				}
			}
		}
		close STATS;
		$logger->info("Pipeline : RAW data screening stats : $statsFile");
	}
}

sub analyse_sra {

	my($runID,$type,$spots,$work_dir,$fraction)=(@_);
	
	if ($fraction){
		$logger->info("Pipeline : Processing : $runID");
		my $analysis_comm="perl $base_path/src/main/scripts/analysis/analyse.SRA.data.pl $input_config $runID $type $spots $work_dir $fraction";
		system("$analysis_comm");
	}
	else {
		$logger->info("Pipeline : Processing : $runID");
		my $analysis_comm="perl $base_path/src/main/scripts/analysis/analyse.SRA.data.pl $input_config $runID $type $spots $work_dir";
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
		my $aln_stats= `tail -n1 $bowtie2_stats_file`;
		chomp $aln_stats; # 2315245   1068    0.0461  3137916 1269    0.0404
		$overall_aln = $aln_stats;
	}
	return $overall_aln;
}

sub summarize {
	my($metadata,$work_dir,$statsFile)=(@_);
	
	my %runIDs;
	# my @runInfoHeaders;
	
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

		if (($blastn_exec_flag == 1) && ($bowtie2_exec_flag == 1)){
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		elsif ($blastn_exec_flag == 1){
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Blast\tUnique_hits_perc_in_Blast\tTotal_reads_in_Sample\tTotal_hits_in_Blast\tTotal_hits_perc_in_Blast\tBlast_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}
		else{
			print STATS "SRA_Run_ID\tUnique_reads_in_Sample\tUnique_hits_in_Bowtie2\tUnique_hits_perc_in_Bowtie2\tTotal_reads_in_Sample\tTotal_hits_in_Bowtie2\tTotal_hits_perc_in_Bowtie2\tBowtie2_mapping_cutoff\t".join("\t",@runInfoHeaders)."\n";
		}

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
	}
}
