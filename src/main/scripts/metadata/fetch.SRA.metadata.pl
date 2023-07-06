#!/usr/bin/env perl

# ARA (Automated Record Analysis) : An automatic pipeline for exploration of SRA datasets with sequences as a query
# Copyright (C) 2023 Anand Maurya;Maciej Szymanski;Wojciech Karlowski

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

use strict;
use warnings;
use Data::Dumper;
use Text::CSV;
use Parallel::ForkManager;
use Cwd;
use Config::Simple;
use Text::Unidecode;
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

my @metadata_recs;

if 	($input_type eq 'list'){
    my @fname = split(/\./,$inpSRARunInfo);
    open(LIST,"$inpSRARunInfo") or die "can't open the file - - - -\n";
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
        my $outFile="$outDir/$sraRun.run.metadata.txt";
        if (!-s "$outFile"){
            open (OUT,">$outFile") or $logger->logdie("Cannot write to file - $outFile: $!");
            print OUT "Run\tRunType\tRunSpots\tSampleAlias\tSampleAccession\tStudyAlias\tStudyAccession\tStudyTitle\tStudyAbstract\t\t\n";
            $logger->info("Fetching metadata for $sraRun");
            my $outString= getRunTable($sraRun);
            print OUT "$outString\n";
            close OUT;
        }
    }
}
else {
    my $outFile="$outDir/$inpSRARunInfo.run.metadata.txt";
    if (!-s "$outFile"){
        open (OUT,">$outFile") or $logger->logdie("Cannot write to file - $outFile: $!");
        print OUT "Run\tRunType\tRunSpots\tSampleAlias\tSampleAccession\tStudyAlias\tStudyAccession\tStudyTitle\tStudyAbstract\t\t\n";
        $logger->info("Fetching metadata for $inpSRARunInfo");
        my $outString= getRunTable($inpSRARunInfo);
        print OUT "$outString\n";
        close OUT;
    }
}

my @files = glob("$outDir/*.run.metadata.txt");

foreach my $file (@files) {
    if (-s "$file"){
        my $rec = `tail -n1 $file`;
        chomp $rec;
        push @metadata_recs, $rec;
    }
}

my $outTmp="$outDir/$outfileLabel.tmp.txt";
my $outFileParsed="$outDir/$outfileLabel.metadata.txt";

open (TMP,">$outTmp") or $logger->logdie("Cannot write to file - $outTmp: $!");
print TMP "Run\tRunType\tRunSpots\tSampleAlias\tSampleAccession\tStudyAlias\tStudyAccession\tStudyTitle\tStudyAbstract\n";
print TMP join("\n",@metadata_recs);
close TMP;

parse_table($outTmp,$outFileParsed);
system("rm $outTmp");


sub getRunTable {
	my ($id)=(@_);
	my ($samHeaders,$samDetails) = getSampleInfo($id);
	my ($runtype,$runspots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,$studyTitle,$study_abstract) = getRunInfo($id);
	my $retStr="$id\t$runtype\t$runspots\t$sampleAlias\t$sampleAccn\t$studyAlias\t$studyAccn\t$studyTitle\t$study_abstract\t$samHeaders\t$samDetails";
	return $retStr;
}

sub getSampleInfo {
	my ($id)=(@_);
	my $saminfoTemp = `$edirectPath/esearch -db sra -query $id | $edirectPath/efetch -format xml | $edirectPath/xtract -pattern EXPERIMENT_PACKAGE  -block SAMPLE -sep '|' -element SAMPLE_ATTRIBUTE/TAG SAMPLE_ATTRIBUTE/VALUE`;
	chomp $saminfoTemp;
	my ($headers,$info)=split(/\t/,$saminfoTemp);
	$headers = unidecode($headers);
	$info = unidecode($info);
	return ($headers,$info);
}


sub getRunInfo {
	my ($id)=(@_);
	
	my $runinfoTmp = `$edirectPath/esearch -db sra -query $id | $edirectPath/efetch -format xml | $edirectPath/xtract -pattern EXPERIMENT_PACKAGE -block RUN_SET -element Statistics\@nreads RUN\@total_spots -block SAMPLE -element SAMPLE\@alias SAMPLE\@accession -block STUDY -element STUDY\@alias STUDY\@accession -block STUDY -element STUDY_TITLE -element STUDY_ABSTRACT`;
	
	chomp $runinfoTmp;

	my ($type,$spots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,$studyTitle,$study_abstract)=split(/\t/,$runinfoTmp);
	if($type > 1){
		$type= "PAIRED";
	}
	else {
		$type= "SINGLE";
	}
	
	if ($studyTitle){
		$studyTitle = unidecode($studyTitle);
	}
	else {
		$studyTitle = "-";
	}

	if ($study_abstract){
		$study_abstract = unidecode($study_abstract);
	}
	else {
		$study_abstract = "-";
	}
	# print($type,$spots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,$studyTitle,$study_abstract);
	return ($type,$spots,$sampleAlias,$sampleAccn,$studyAlias,$studyAccn,$studyTitle,$study_abstract);
}



sub parse_table {
	my ($rawMetadataTable,$parsedMetadataTable)=(@_);
	# print "$rawMetadataTable\n";
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
			if (($tName ne 'RunType') and ($tName ne 'StudyTitle')){
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
