#!/usr/bin/perl -w
# Copyright (c) yuhy 2014/4/30
# Writer:         yuhy <yuhy@xxx.com.cn>
# Program Date:   2014/4/30.
# Modifier:       yuhy <yuhy@xxx.com.cn>
# Last Modified:  2014/4/30.

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);

my $programe_dir=basename($0);
my $path=dirname($0);

my $ver    = "1.0";
my $Writer = "yuhy <yuhy\@yuhyxxx@163.com>";
my $Data   = "2014/12/09";
my $BEGIN=time();
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in,$out,$snp,$deg,$indel,$sin);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"i:s"=>\$in,
			) || &help;
&help unless ($in && $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: ......
    Usage:
        -i          infile     must be given

        -o          outfile    must be given

        -h          Help document
	Usage End.
	exit;
}
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
################
$/ = "\n";
open (IN,"$in") || die $!;
open (OUT,">$out") || die $!;
my (%hash,%pep);
my ($l,$spe);
$spe=0;
while (<IN>){
	chomp;
	my $file = $_;
	open (FILE,"$file") || die $!;
	$/ = ">";
	<FILE>;
	while(<FILE>){
		chomp;
		my ($id,$seq) = split /\n/,$_,2;
		$seq=~s/\n//g;
		$seq=~s/ //g;
		my ($key) = $id =~ /^([a-zA-Z]*?)\_/;
		#print $key;die;
		if(exists $hash{$key}){
			$hash{$key} .= $seq;
		}else{
			$hash{$key} = $seq;
		}
	}
	$/ = "\n";
	close FILE;
}
my $gene;
my $info;
foreach $gene(keys %hash){
	$l=length($hash{$gene});
	$spe++;
	if($spe==1){
		$info = $gene . "  " . $hash{$gene};
	}else{
		$info = $info . "\n" . $gene . "  " . $hash{$gene};
	}
}
print OUT "$spe  $l\n$info\n";



###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime # &Runtime($BEGIN);
{
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe_dir elapsed time : [",&sub_time($t),"]\n";
}
sub sub_time
{
	my ($T)=@_;chomp $T;
	my $s=0;my $m=0;my $h=0;
	if ($T>=3600) {
		my $h=int ($T/3600);
		my $a=$T%3600;
		if ($a>=60) {
			my $m=int($a/60);
			$s=$a%60;
			$T=$h."h\-".$m."m\-".$s."s";
		}else{
			$T=$h."h-"."0m\-".$a."s";
		}
	}else{
		if ($T>=60) {
			my $m=int($T/60);
			$s=$T%60;
			$T=$m."m\-".$s."s";
		}else{
			$T=$T."s";
		}
	}
	return ($T);
}

sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub getSeq	{
	my $file=shift;
	my %chr;
	$/=">";
	open FA,$file;
	<FA>;
	while(<FA>)
	{
		chomp;
		my($head,$seq)=split/\n+/,$_,2;
		my $id =(split/\s+/,$head)[0];
		$seq=~s/\n+//g;
		$seq=~s/>$//;
		$chr{$id}=$seq;
	}
	close FA;
	$/="\n";
	return (%chr);
}
