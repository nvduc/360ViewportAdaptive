#! /bin/perl -w
use strict;
use warnings;
# ================================================
my $copied_file_log = "file_list_copy.txt";
my $FILE;
my @array_log;
my $len;
my $fname;
my $cmd;
my $ORI_FOLDER = "yakitori_3840x1920_896fr_30fps";
my @TILING_W = qw(1 8);
my @TILING_H = qw(1 8);
my $NO_TILING = @TILING_W;
my @INTER_LIST = qw(4 16 32 64);
my $TILE_ID;
my $NO_TILE;
my $INTER;
my $TID;
my $old_fname;
my $new_fname;
my $cmd;
for ($TILE_ID = 0; $TILE_ID < $NO_TILING; $TILE_ID++) {
	$NO_TILE = $TILING_W[$TILE_ID] * $TILING_H[$TILE_ID];
	foreach $INTER (@INTER_LIST){
		for ($TID = 0; $TID < $NO_TILE; $TID++) {
			$old_fname = "$ORI_FOLDER/tile_${TILING_W[$TILE_ID]}x${TILING_H[$TILE_ID]}/${INTER}frame/tile_${TID}.txt";
			$new_fname = "$ORI_FOLDER/tile_${TILING_W[$TILE_ID]}x${TILING_H[$TILE_ID]}/${INTER}frame/f0_t${TID}.txt";
			$cmd = "cp $old_fname $new_fname";
			print "$cmd\n";
			system "$cmd";
			$old_fname = "$ORI_FOLDER/tile_${TILING_W[$TILE_ID]}x${TILING_H[$TILE_ID]}/${INTER}frame/per_frame/tile_${TID}.txt";
			$new_fname = "$ORI_FOLDER/tile_${TILING_W[$TILE_ID]}x${TILING_H[$TILE_ID]}/${INTER}frame/per_frame/f0_t${TID}.txt";
			$cmd = "cp $old_fname $new_fname";
			print "$cmd\n";
			# exit;
			system "$cmd";
			# exit;
		}
	}
}
