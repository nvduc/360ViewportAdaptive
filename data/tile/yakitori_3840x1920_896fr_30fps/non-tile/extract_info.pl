#!/bin/perl -w
use strict;
use warnings;

my @QP_ar = qw (48 44 40 36 32 28 24);
my $FPS = 30;
my $fnum;
my $vernum = 7;
my @avgBR;
my @avgPSNR;
my $fname;
my $LOG;
my $fid;
my $line;
my @arr;
my $i;
#
$fname = "DASH_frame.txt";
open $LOG, $fname or die "Could not open file '$fname'";
$fid = 0;
while(<$LOG>){
    chomp;
    $line = $_;
#    print "$line\n";
    @arr = $line =~ /(\S+)/g;
    for($i=0; $i < $vernum; $i++){
	$avgBR[$i] += $arr[2*$i];
	$avgPSNR[$i] += $arr[2*$i + 1];
    }
    $fid ++;
}
close $LOG;
#
$fnum = $fid;
print "$fnum frames loaded \n";
$fname = "stat.txt";
open $LOG, '>', $fname or die "Could not open file '$fname'";
for($i=0; $i < $vernum; $i++){
    $avgBR[$i] = $avgBR[$i] / ($fnum / $FPS) / 1000.0;
    $avgPSNR[$i] = $avgPSNR[$i] / $fnum;
    printf $LOG "%.2f\t%.2f\n", $avgBR[$i], $avgPSNR[$i];
}
close $LOG;