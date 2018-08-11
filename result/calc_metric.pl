#! /bin/perl -w
use strict;
use warnings;
use List::Util qw(sum);
sub mean {
    return sum(@_)/@_;
}
#=============================================#
my $vname = "timelapseNY";
my $W = 3840;
my $H = 1920;
my @method_arr = qw{1 2 3 4 5 6 13 11 12};
my @DASH_VER_arr = qw{1 2 3 4 5 6};
my @theta_arr = qw{90};
my $INTER = 24;
my $BUFF = 24;
my $EST = 0;
my $METHOD;
my $phi;
my $phi_tmp;
my $theta;
my $DASH_VER;
my $fileName;
my $outFile;
my $logFile;
my $dir = "video_${vname}/${W}x${H}/fixedDASH/fixedVP/theta=0";
my $line;
my $frame_id;
my $cnt;
my @word;
my @tmp;
my @avg_vpsnr;
my $phi_id;
my $dash_id;
my $theta_id;
my $method_id;
# 
$dash_id = 0;
foreach $DASH_VER (@DASH_VER_arr) {	
	$theta_id = 0;
	foreach $theta (@theta_arr){
		$dir = "video_${vname}/${W}x${H}/fixedDASH/fixedVP/theta=${theta}";
		$fileName = "${dir}/avg_fixed_vp_DASH_VER_${DASH_VER}_theta_${theta}.txt";
		open $outFile, '>', $fileName;
		$phi_id = 0;
		for ($phi= -165;  $phi <= 180; $phi += 15) {
			if($phi < 0){
				$phi_tmp = $phi + 360;
			}else{
				$phi_tmp = $phi;
			}
			$method_id = 0;
			printf $outFile "$phi\t";
			foreach $METHOD (@method_arr){
				$fileName = "${dir}/log_frame_DASH_VER_${DASH_VER}_phi_${phi_tmp}_theta_${theta}_METHOD_${METHOD}_INTER_${INTER}_BUFF_${BUFF}_EST_${EST}.txt";
				print "$fileName\n";
				open $logFile, $fileName or die "Cannot open file '${fileName}'\n";
				$frame_id = 0;
				$cnt = 0;
				while(<$logFile>){
					chomp;
					$line = $_;
					if($line =~ /fid/){
						# print "$line\n";
						# exit;
					}else{
						@word = $line =~ /(\d+)/g;
						$frame_id = "$word[0]";
						if($frame_id >= 2 * $INTER){
							$tmp[$cnt] = "$word[2].$word[3]";
							# print "frame=${frame_id} tmp=$tmp[$cnt]\n";
							$cnt ++;
							# exit;
						}
					}
				}
				$avg_vpsnr[$dash_id][$theta_id][$phi_id][$method_id] = mean(@tmp);
				print "$avg_vpsnr[$dash_id][$theta_id][$phi_id][$method_id]\n";
				printf $outFile "%.2f\t",$avg_vpsnr[$dash_id][$theta_id][$phi_id][$method_id];
			$method_id ++;
			}
			print $outFile "\n";
		$phi_id ++;
		}
	$theta_id++;
	} 
$dash_id ++;
}
