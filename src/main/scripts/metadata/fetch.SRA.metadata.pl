#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Text::CSV;
use Parallel::ForkManager;
use Cwd;
use Config::Simple;
use Log::Log4perl qw(get_logger);


my $inpSRARunInfo = $ARGV[0];	# qc_PASS_sample_total_data_sorted.txt
my $outDir = $ARGV[1];   		# /mnt/Storage1/hpz/
my $outfileLabel = $ARGV[2];	# test-list
my $script_dir = $ARGV[3];		# parent directory of wrapper script
my $edirectPath = $ARGV[4]; 	# /home/anand/edirect
my $input_type = $ARGV[5];		# single or list

## Initialize Log4perl, if it is not initialized already by the calling script
unless (Log::Log4perl::initialized()) {
	Log::Log4perl::init("$script_dir/src/main/resources/log4perl.config");
}
my $logger = Log::Log4perl->get_logger();
my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1 });


my $outFile="$outDir/$outfileLabel.tmp.txt";
my $outFileParsed="$outDir/$outfileLabel.metadata.txt";


fetch_metadata($inpSRARunInfo,$outFile);
parse_table($outFile,$outFileParsed);
system("rm $outFile");

$logger->info("Metadata table generated:\t$outFileParsed");



sub fetch_metadata {
	my ($inp,$out)=(@_);

	open (OUT,">$out") or $logger->logdie("Cannot write to file - $out: $!");
	print OUT "Run\tRunType\tRunSpots\tSampleAlias\tSampleAccession\tStudyAlias\tStudyAccession\tStudyTitle\n";
	
	if 	($input_type eq 'list'){
		my @fname = split(/\./,$inp);
		open(LIST,"$inp") or die "can't open the file - - - -\n";
		while(my $rec=<LIST>){
			chomp $rec;
			if (($rec=~ /^Run/) or ($rec eq "")){next;}
			
			my $sraRun; my @cols; my $str='-';

			if ($fname[-1] eq 'csv'){
				my $res = $csv->parse($rec);
				@cols = $csv->fields ();
				$sraRun = shift(@cols);
				# $str = join ("\t", @cols);
			}
			else {
				@cols = split(/\s+/,$rec);
				$sraRun = shift(@cols);
				# $str = join ("\t", @cols);
			}
			$logger->info("Fetching metadata for $sraRun");
			my $outString= getRunTable($sraRun);
			print OUT "$outString\n";
		}
	}
	else {
		$logger->info("Fetching metadata for $inp");
		my $outString= getRunTable($inp);
		print OUT "$outString\n";
	}

}

sub getRunTable {
	my ($id)=(@_);
	my ($samHeaders,$samDetails) = getSampleInfo($id);
	my ($runtype,$runspots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,$studyTitle) = getRunInfo($id);
	# my $retStr="$str\t$samHeaders\t$samDetails";
	my $retStr="$id\t$runtype\t$runspots\t$sampleAlias\t$sampleAccn\t$studyAlias\t$studyAccn\t$studyTitle\t$samHeaders\t$samDetails";
	return $retStr;
}

sub getSampleInfo {
	my ($id)=(@_);
	my $saminfoTemp = `$edirectPath/esearch -db sra -query $id | $edirectPath/efetch -format xml | $edirectPath/xtract -pattern EXPERIMENT_PACKAGE  -block SAMPLE -sep '|' -element SAMPLE_ATTRIBUTE/TAG SAMPLE_ATTRIBUTE/VALUE`;
	chomp $saminfoTemp;
	my ($headers,$info)=split(/\t/,$saminfoTemp);
	return ($headers,$info);
}


sub getRunInfo {
	my ($id)=(@_);
	
	my $studyTitle = "-";
	my $runinfoTmp = `$edirectPath/esearch -db sra -query $id | $edirectPath/efetch -format xml | $edirectPath/xtract -pattern EXPERIMENT_PACKAGE -block RUN_SET -element Statistics\@nreads RUN\@total_spots -block SAMPLE -element SAMPLE\@alias SAMPLE\@accession -block STUDY -element STUDY\@alias STUDY\@accession -block STUDY -element STUDY_TITLE`;
	
	chomp $runinfoTmp;
	
	my ($type,$spots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,@title)=split(/\s+/,$runinfoTmp);
	if($type > 1){
		$type= "PAIRED";
	}
	else {
		$type= "SINGLE";
	}
	
	if (scalar @title == 0){
		$studyTitle = "-";
	}
	else{
		$studyTitle = join(" ",@title);
	}
	return ($type,$spots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,$studyTitle);
}


sub parse_table {
	my ($rawMetadataTable,$parsedMetadataTable)=(@_);

	my $getTagsCol= q(awk '{FS="\t"}{print $(NF-1)}' ).$rawMetadataTable."|sort|uniq";
	my $uniqTags=`$getTagsCol`;
	chomp $uniqTags;
	# print "--$uniqTags--\n";

	my @tagArrTemp = split(/\n/,$uniqTags);
	my %tagNames;

	foreach my $tagSet (@tagArrTemp){
		#print "$tagSet\n";
		my @tagTemp=split(/\|/,$tagSet);
		foreach my $tName (@tagTemp){
			if ($tName ne 'RunType'){
				$tagNames{$tName}=$tName;
			}
		}
	}

	# print Dumper(\%tagNames);

	my @tagArrFull;
	foreach my $tag(sort keys %tagNames){
		push @tagArrFull,$tag;
	}

	# print join ("\n",@tagArrFull);

	open (PARSEDTAB,">$parsedMetadataTable") or $logger->logdie("Cannot write to file - $parsedMetadataTable: $!");

	open(META,"$rawMetadataTable") or $logger->logdie("Cannot open raw metadata file - $rawMetadataTable: $!");

	my $header = <META>; 
	chomp $header;

	print PARSEDTAB "$header\t".join ("\t",@tagArrFull)."\n";


	while(my $rec=<META>){
		chomp $rec;
		my @recTemp=split(/\t/,$rec);
		my @tagList; my $tagAttrList;
		my %rowTagAttrMap;
		my @rowAttrsOut;
		if ($recTemp[-2] && $recTemp[-1]){
			my @tagList=split(/\|/,$recTemp[-2]);   
			my @tagAttrList=split(/\|/,$recTemp[-1]);
			
			for(my $i=0;$i<=$#tagList;$i++){
				$rowTagAttrMap{$tagList[$i]}=$tagAttrList[$i];
			}
			# print Dumper(\%rowTagAttrMap);
		}
		foreach my $tag(sort keys %tagNames){
			if (defined $rowTagAttrMap{$tag}){
				push @rowAttrsOut,$rowTagAttrMap{$tag};
			}
			else{
				push @rowAttrsOut, "-";
			}
		}
		pop @recTemp; pop @recTemp;
		
		print PARSEDTAB join ("\t",@recTemp)."\t";
		print PARSEDTAB join ("\t",@rowAttrsOut)."\n";
	}

	close PARSEDTAB;
}
