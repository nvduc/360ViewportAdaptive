#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Projection.h"
#include "common.h"
using namespace std;

#ifdef DEBUG
int main(int argc, char const *argv[])
{
	// if(argc != 2){
	// 	printf("usage: ./projection tile_id\n");
	// 	exit(1);
	// }
	int NO_FACE = 1;
	int NO_TILE_H = 8;
	int NO_TILE_W = 8;
	// int tid = atoi(argv[1]);
	// int N;
	// int *adjTile;
	// adjTile = getAdjacentTileERP2(NO_TILE_H, NO_TILE_W, tid, &N);
	// showArray(adjTile, N);
	int vmask[] = {0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,0,0,
				  0,0,0,1,1,0,0,0,
				  0,0,0,1,1,0,0,0,
				  0,0,0,1,1,0,0,0,
				  0,0,0,1,1,0,0,0,
				  0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,0,0};
	int* vmask_est = extVmask(vmask, NO_FACE, NO_TILE_W, NO_TILE_H, 2);
	showTileVersion(vmask_est, NO_TILE_W, NO_TILE_H);
	int* vmask_ext2 = extVmask2(vmask, NO_FACE, NO_TILE_W, NO_TILE_H, 2);
	printf("# ext2:\n");
	showTileVersion(vmask_ext2, NO_TILE_W, NO_TILE_H);
}
#endif

int* getAdjacentTileERP(int NO_TILE_H, int NO_TILE_W, int tid, int* N){
	int tid_w = tid / NO_TILE_H; /* vertical position */
	int tid_h = tid - tid_w * NO_TILE_H; /* Horizontal position */
	int MAX = 16;
	int cnt = 0, i;
	Tile* adjTiles = new Tile[MAX];
	int* ret;
	// top
	if(tid_w - 1 >= 0){
		// (w-1, h)
		adjTiles[cnt].tid_h = tid_h;
		adjTiles[cnt++].tid_w = tid_w - 1;
		// (w-1, h-1)
		adjTiles[cnt].tid_h = ((tid_h - 1) >= 0)?(tid_h - 1):(NO_TILE_H - 1);
		adjTiles[cnt++].tid_w = tid_w - 1;
		// (w-1, h+1)
		adjTiles[cnt].tid_h = ((tid_h + 1) < NO_TILE_H)?(tid_h + 1):0;
		adjTiles[cnt++].tid_w = tid_w - 1;
	}
	// bottom
	if(tid_w + 1 < NO_TILE_W){
		// (w+1, h)
		adjTiles[cnt].tid_h = tid_h;
		adjTiles[cnt++].tid_w = tid_w + 1;
		// (w+1, h-1)
		adjTiles[cnt].tid_h = ((tid_h - 1) >= 0)?(tid_h - 1):(NO_TILE_H - 1);
		adjTiles[cnt++].tid_w = tid_w + 1;
		// (w+1, h+1)
		adjTiles[cnt].tid_h = ((tid_h + 1) < NO_TILE_H)?(tid_h + 1):0;
		adjTiles[cnt++].tid_w = tid_w + 1;
	}
	// middle
	// (w, h-1)
	adjTiles[cnt].tid_h = ((tid_h - 1) >= 0)?(tid_h - 1):(NO_TILE_H - 1);
	adjTiles[cnt++].tid_w = tid_w;
	// (w, h+1)
	adjTiles[cnt].tid_h = ((tid_h + 1) < NO_TILE_H)?(tid_h + 1):0;
	adjTiles[cnt++].tid_w = tid_w;
	// top and bottom tiles
	if(tid_w == 0 || tid_w == NO_TILE_W - 1){
		for(i=0; i < NO_TILE_H; i++){
			if(i != tid_h && i != ((tid_h - 1 >=0)?(tid_h - 1):(NO_TILE_H-1)) && i != ((tid_h + 1 < NO_TILE_H)?(tid_h+1):0)){
				// printf("i=%d\n", i);
				adjTiles[cnt].tid_h = i;
				adjTiles[cnt].tid_w = tid_w;
				cnt++;
			}
		}
	}
	*N = cnt;
	ret = new int[*N];
	/* convert to tile's id*/
	for(i=0; i < *N; i++){
		ret[i] = adjTiles[i].tid_w * NO_TILE_H + adjTiles[i].tid_h;
	}
	return ret;
}
int* getAdjacentTileERP2(int NO_TILE_H, int NO_TILE_W, int tid, int* N){
	int tid_w; /* vertical position */
	int tid_h; /* Horizontal position */
	int cnt = 0, i, K, MAX = 8;
	int* ret = new int[MAX];
	int* tmp;
	int NO_TILE_W2 = NO_TILE_W + 2;
	bool TOP_TILE = false, BOT_TILE = false;
	//
	if(tid == 0 || tid == 1){ // top or bottom tiles
		tid_w = (tid == 0)?0:(NO_TILE_W-1);
		//printf("tid_w=%d\n", tid_w);
		for(tid_h =0; tid_h < NO_TILE_H; tid_h ++){
			ret[tid_h] = 2 + (tid_w) * NO_TILE_H + tid_h;
		}
		*N = NO_TILE_H;
	}else{
		tid_w = (tid - 2 + NO_TILE_H)/NO_TILE_H;
		tid_h = (tid - 2 + NO_TILE_H)%NO_TILE_H;
		ret = getAdjacentTileERP(NO_TILE_H, NO_TILE_W + 2, tid_w*NO_TILE_H + tid_h, &K);
		//printf("tid_w=%d tid_h=%d\n", tid_w, tid_h);
		showArray(ret, K);
		//
		for(i=0; i < K; i++){
			// indices on MXN grid
			tid_w = ret[i] / NO_TILE_H;
			tid_h = ret[i] % NO_TILE_H;
			// if the tile belongs to poles tiles
			if(tid_w == 0){
				if(!TOP_TILE){
					ret[cnt++] = 0;
					TOP_TILE = true;
				}
			}else{
				if(tid_w == NO_TILE_W2 - 1){
					if(!BOT_TILE){
						ret[cnt++] = 1;
						BOT_TILE = true;
					}
				}else{
					ret[cnt++] = 2 + (tid_w - 1) * NO_TILE_H + tid_h;
				}
			}
			//printf("i=%d tid_w=%d tid_h=%d cnt=%d last=%d\n", i, tid_w, tid_h, cnt, ret[cnt-1]);
		}
		*N = cnt;
	}
	return ret;
}
int* getAdjacentTileCMP(int NO_TILE_H, int NO_TILE_W, int tid, int* N){
	int i, MAX = 8, M, cnt;
	int* ret;
	Tile tile;
	tile.fid = tid / (NO_TILE_H * NO_TILE_W);
	tile.tid_w = (tid - (tile.fid * (NO_TILE_H * NO_TILE_W)))/NO_TILE_H;
	tile.tid_h = tid - tile.fid * (NO_TILE_H * NO_TILE_W) - tile.tid_w * NO_TILE_H;	
	Tile* tiles = new Tile[MAX + 1];
	Tile* adjTile;
	if(NO_TILE_H >= 2 && NO_TILE_W >= 2){
		if(tile.fid == 0){
			/*  type 1: in-face tiles */
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				// tiles[0]  = new Tile();
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w;
				//5
				tiles[5].fid = tile.fid;
				tiles[5].tid_h = tile.tid_h - 1;
				tiles[5].tid_w = tile.tid_w + 1;
				//6
				tiles[6].fid = tile.fid;
				tiles[6].tid_h = tile.tid_h;
				tiles[6].tid_w = tile.tid_w + 1;
				//7
				tiles[7].fid = tile.fid;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			/* type 2: vertex tiles */
			if(tile.tid_h == 0 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h +1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w + 1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == 0 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h +1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 4;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 4;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h -1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 3;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 3;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = NO_TILE_W - 2;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			/* type 3: other tiles */
			if(tile.tid_w == 0 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h-1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w +1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 2;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = tile.tid_h;
				//6
				tiles[6].fid = 2;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = tile.tid_h - 1;
				//7
				tiles[7].fid = 2;
				tiles[7].tid_h = NO_TILE_H - 1;
				tiles[7].tid_w = tile.tid_h + 1;
				*N = 8;
			}
			if(tile.tid_w == NO_TILE_W - 1 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h -1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w -1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w - 1;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = tile.tid_h;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = tile.tid_h - 1;
				//7
				tiles[7].fid = 3;
				tiles[7].tid_h = NO_TILE_H - 1;
				tiles[7].tid_w = tile.tid_h + 1;
				*N = 8;
			}
			if(tile.tid_h == 0 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h + 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 4;
				tiles[7].tid_h = NO_TILE_H - 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h - 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 5;
				tiles[7].tid_h = 0;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}		
		}
		if(tile.fid == 5){
			/*  type 1: in-face tiles */
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				// tiles[0]  = new Tile();
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w;
				//5
				tiles[5].fid = tile.fid;
				tiles[5].tid_h = tile.tid_h - 1;
				tiles[5].tid_w = tile.tid_w + 1;
				//6
				tiles[6].fid = tile.fid;
				tiles[6].tid_h = tile.tid_h;
				tiles[6].tid_w = tile.tid_w + 1;
				//7
				tiles[7].fid = tile.fid;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			/* type 2: vertex tiles */
			if(tile.tid_h == 0 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h + 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w + 1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = NO_TILE_H - 2;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == 0 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h +1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 0;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 0;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h -1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 3;
				tiles[3].tid_h = 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 3;
				tiles[4].tid_h = 0;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = NO_TILE_W - 2;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			/* type 3: other tiles */
			if(tile.tid_w == 0 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h-1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w +1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 2;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 2;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = 0;
				//7
				tiles[7].fid = 2;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = 0;
				*N = 8;
			}
			if(tile.tid_w == NO_TILE_W - 1 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h -1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w -1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w - 1;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = NO_TILE_W - 1;
				//7
				tiles[7].fid = 3;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = NO_TILE_W - 1;
				*N = 8;
			}
			if(tile.tid_h == 0 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h + 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 0;
				tiles[7].tid_h = NO_TILE_H - 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h - 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 1;
				tiles[7].tid_h = 0;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}		
		}
		if(tile.fid == 1){
			/*  type 1: in-face tiles */
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				// tiles[0]  = new Tile();
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w;
				//5
				tiles[5].fid = tile.fid;
				tiles[5].tid_h = tile.tid_h - 1;
				tiles[5].tid_w = tile.tid_w + 1;
				//6
				tiles[6].fid = tile.fid;
				tiles[6].tid_h = tile.tid_h;
				tiles[6].tid_w = tile.tid_w + 1;
				//7
				tiles[7].fid = tile.fid;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			/* type 2: vertex tiles */
			if(tile.tid_h == 0 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h + 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w + 1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = 0;
				tiles[4].tid_w = 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = 0;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == 0 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h +1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 5;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 5;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = NO_TILE_W - 2;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h -1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 3;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 3;
				tiles[4].tid_h = 0;
				tiles[4].tid_w = 1;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = NO_TILE_W - 2;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			/* type 3: other tiles */
			if(tile.tid_w == 0 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h-1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w +1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 2;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = tile.tid_h;
				//6
				tiles[6].fid = 2;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = tile.tid_h - 1;
				//7
				tiles[7].fid = 2;
				tiles[7].tid_h = 0;
				tiles[7].tid_w = tile.tid_h + 1;
				*N = 8;
			}
			if(tile.tid_w == NO_TILE_W - 1 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w -1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w - 1;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = tile.tid_h;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = tile.tid_h - 1;
				//7
				tiles[7].fid = 3;
				tiles[7].tid_h = 0;
				tiles[7].tid_w = tile.tid_h + 1;
				*N = 8;
			}
			if(tile.tid_h == 0 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h + 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 5;
				tiles[7].tid_h = NO_TILE_H - 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h - 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 4;
				tiles[7].tid_h = 0;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}		
		}
		if(tile.fid == 4){
			/*  type 1: in-face tiles */
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				// tiles[0]  = new Tile();
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w;
				//5
				tiles[5].fid = tile.fid;
				tiles[5].tid_h = tile.tid_h - 1;
				tiles[5].tid_w = tile.tid_w + 1;
				//6
				tiles[6].fid = tile.fid;
				tiles[6].tid_h = tile.tid_h;
				tiles[6].tid_w = tile.tid_w + 1;
				//7
				tiles[7].fid = tile.fid;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			/* type 2: vertex tiles */
			if(tile.tid_h == 0 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h + 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w + 1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = 2;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 2;
				tiles[4].tid_h = NO_TILE_H - 2;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = 1;
				*N = 7;
			}
			if(tile.tid_h == 0 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h +1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 1;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 1;
				tiles[4].tid_h = NO_TILE_H - 1;
				tiles[4].tid_w = NO_TILE_W - 2;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = 1;
				tiles[6].tid_w = 0;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h -1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 3;
				tiles[3].tid_h = NO_TILE_H - 1;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 3;
				tiles[4].tid_h = NO_TILE_H - 2;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = NO_TILE_W - 2;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			/* type 3: other tiles */
			if(tile.tid_w == 0 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h-1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w +1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 2;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 2;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = NO_TILE_W - 1;
				//7
				tiles[7].fid = 2;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = NO_TILE_W - 1;
				*N = 8;
			}
			if(tile.tid_w == NO_TILE_W - 1 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w -1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w - 1;
				//5
				tiles[5].fid = 3;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 3;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = 0;
				//7
				tiles[7].fid = 3;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = 0;
				*N = 8;
			}
			if(tile.tid_h == 0 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h + 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = NO_TILE_H - 1;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 1;
				tiles[7].tid_h = NO_TILE_H - 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h - 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = 0;
				tiles[5].tid_w = tile.tid_w;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = 0;
				tiles[6].tid_w = tile.tid_w - 1;
				//7
				tiles[7].fid = 0;
				tiles[7].tid_h = 0;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}		
		}
		if(tile.fid == 3){
			/*  type 1: in-face tiles */
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				// tiles[0]  = new Tile();
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w;
				//5
				tiles[5].fid = tile.fid;
				tiles[5].tid_h = tile.tid_h - 1;
				tiles[5].tid_w = tile.tid_w + 1;
				//6
				tiles[6].fid = tile.fid;
				tiles[6].tid_h = tile.tid_h;
				tiles[6].tid_w = tile.tid_w + 1;
				//7
				tiles[7].fid = tile.fid;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			/* type 2: vertex tiles */
			if(tile.tid_h == 0 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h + 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w + 1;
				//3
				tiles[3].fid = 4;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 4;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = 0;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 0;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			if(tile.tid_h == 0 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h +1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 1;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 1;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h -1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 5;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = NO_TILE_W - 1;
				//4
				tiles[4].fid = 5;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = NO_TILE_W - 1;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = NO_TILE_W - 1;
				*N = 7;
			}
			/* type 3: other tiles */
			if(tile.tid_w == 0 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h-1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w +1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = NO_TILE_W - 1;
				//7
				tiles[7].fid = 4;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = NO_TILE_W - 1;
				*N = 8;
			}
			if(tile.tid_w == NO_TILE_W - 1 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w -1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w - 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = NO_TILE_W - 1;
				//7
				tiles[7].fid = 5;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = NO_TILE_W - 1;
				*N = 8;
			}
			if(tile.tid_h == 0 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h + 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = tile.tid_w;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = tile.tid_w - 1;
				tiles[6].tid_w = NO_TILE_W - 1;
				//7
				tiles[7].fid = 1;
				tiles[7].tid_h = tile.tid_w + 1;
				tiles[7].tid_w = NO_TILE_W - 1;
				*N = 8;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h - 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = tile.tid_w;
				tiles[5].tid_w = NO_TILE_W - 1;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = tile.tid_w - 1;
				tiles[6].tid_w = NO_TILE_W - 1;
				//7
				tiles[7].fid = 0;
				tiles[7].tid_h = tile.tid_w + 1;
				tiles[7].tid_w = NO_TILE_W - 1;
				*N = 8;
			}		
		}
		if(tile.fid == 2){
			/*  type 1: in-face tiles */
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				// tiles[0]  = new Tile();
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w;
				//5
				tiles[5].fid = tile.fid;
				tiles[5].tid_h = tile.tid_h - 1;
				tiles[5].tid_w = tile.tid_w + 1;
				//6
				tiles[6].fid = tile.fid;
				tiles[6].tid_h = tile.tid_h;
				tiles[6].tid_w = tile.tid_w + 1;
				//7
				tiles[7].fid = tile.fid;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = tile.tid_w + 1;
				*N = 8;
			}
			/* type 2: vertex tiles */
			if(tile.tid_h == 0 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h + 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w + 1;
				//3
				tiles[3].fid = 1;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 1;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = 0;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == 0){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w + 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = 5;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 5;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = 0;
				*N = 7;
			}
			if(tile.tid_h == 0 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h +1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 4;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 4;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = 0;
				*N = 7;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w == NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w - 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h -1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = 0;
				tiles[3].tid_h = 0;
				tiles[3].tid_w = 0;
				//4
				tiles[4].fid = 0;
				tiles[4].tid_h = 1;
				tiles[4].tid_w = 0;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = NO_TILE_H - 1;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = NO_TILE_H - 2;
				tiles[6].tid_w = 0;
				*N = 7;
			}
			/* type 3: other tiles */
			if(tile.tid_w == 0 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h-1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h -1;
				tiles[2].tid_w = tile.tid_w +1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w +1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 5;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 5;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = 0;
				//7
				tiles[7].fid = 5;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = 0;
				*N = 8;
			}
			if(tile.tid_w == NO_TILE_W - 1 && tile.tid_h > 0 && tile.tid_h < NO_TILE_H - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h - 1;
				tiles[0].tid_w = tile.tid_w;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h + 1;
				tiles[1].tid_w = tile.tid_w;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h;
				tiles[3].tid_w = tile.tid_w -1;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w - 1;
				//5
				tiles[5].fid = 4;
				tiles[5].tid_h = tile.tid_h;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 4;
				tiles[6].tid_h = tile.tid_h - 1;
				tiles[6].tid_w = 0;
				//7
				tiles[7].fid = 4;
				tiles[7].tid_h = tile.tid_h + 1;
				tiles[7].tid_w = 0;
				*N = 8;
			}
			if(tile.tid_h == 0 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h + 1;
				tiles[2].tid_w = tile.tid_w -1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h + 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h + 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 1;
				tiles[5].tid_h = tile.tid_w;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 1;
				tiles[6].tid_h = tile.tid_w - 1;
				tiles[6].tid_w = 0;
				//7
				tiles[7].fid = 1;
				tiles[7].tid_h = tile.tid_w + 1;
				tiles[7].tid_w = 0;
				*N = 8;
			}
			if(tile.tid_h == NO_TILE_H - 1 && tile.tid_w > 0 && tile.tid_w < NO_TILE_W - 1){
				//0
				tiles[0].fid = tile.fid;
				tiles[0].tid_h = tile.tid_h;
				tiles[0].tid_w = tile.tid_w + 1;
				//1
				tiles[1].fid = tile.fid;
				tiles[1].tid_h = tile.tid_h;
				tiles[1].tid_w = tile.tid_w - 1;
				//2
				tiles[2].fid = tile.fid;
				tiles[2].tid_h = tile.tid_h - 1;
				tiles[2].tid_w = tile.tid_w - 1;
				//3
				tiles[3].fid = tile.fid;
				tiles[3].tid_h = tile.tid_h - 1;
				tiles[3].tid_w = tile.tid_w;
				//4
				tiles[4].fid = tile.fid;
				tiles[4].tid_h = tile.tid_h - 1;
				tiles[4].tid_w = tile.tid_w + 1;
				//5
				tiles[5].fid = 0;
				tiles[5].tid_h = tile.tid_w;
				tiles[5].tid_w = 0;
				//6
				tiles[6].fid = 0;
				tiles[6].tid_h = tile.tid_w - 1;
				tiles[6].tid_w = 0;
				//7
				tiles[7].fid = 0;
				tiles[7].tid_h = tile.tid_w + 1;
				tiles[7].tid_w = 0;
				*N = 8;
			}		
		}
	}
	if(NO_TILE_W == 1){
		/* get in-face adjcent tiles */
		adjTile = get_adj_tile_in_face(NO_TILE_H, NO_TILE_W, tile.tid_w, tile.tid_h, &M);
		for(i=0; i < M; i++){
			tiles[i] = adjTile[i];
			tiles[i].fid = tile.fid;
			// printf("%d %d %d \n", tiles[i].fid,tiles[i].tid_w, tiles[i].tid_h);
		}
		cnt = M;
		if(tile.fid == 0){
			/* boundary tiles */
			if(tile.tid_h == 0){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
			}
			if(tile.tid_h == NO_TILE_H - 1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;	
			}
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H-1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
			}
			*N = cnt;
		}
		if(tile.fid == 5){
			/* boundary tiles */
			if(tile.tid_h == 0){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 2;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 2;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 0;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
			}
			if(tile.tid_h == NO_TILE_H - 1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 1;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;	
			}
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H-1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
			}
			*N = cnt;
		}
		if(tile.fid == 1){
			/* boundary tiles */
			if(tile.tid_h == 0){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
			}
			if(tile.tid_h == NO_TILE_H - 1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;	
			}
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H-1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
			}
			*N = cnt;
		}
		if(tile.fid == 4){
			/* boundary tiles */
			if(tile.tid_h == 0){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 1;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
			}
			if(tile.tid_h == NO_TILE_H - 1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = NO_TILE_H - 2;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = NO_TILE_H - 2;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 0;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;	
			}
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H-1){
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//2
				tiles[cnt].fid = 2;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//3
				tiles[cnt].fid = 3;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
			}
			*N = cnt;
		}
		if(tile.fid == 2){
			if(tile.tid_h == 0){
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = NO_TILE_H-1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = NO_TILE_H-2;
				tiles[cnt++].tid_w = 0;
				//1
				for(i=0; i < NO_TILE_H; i++){
					tiles[cnt].fid = 1;
					tiles[cnt].tid_h = i;
					tiles[cnt++].tid_w = 0;
				}
			}
			if(tile.tid_h == NO_TILE_H-1){
				//4
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = NO_TILE_H-1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = NO_TILE_H-2;
				tiles[cnt++].tid_w = 0;
				//1
				for(i=0; i < NO_TILE_H; i++){
					tiles[cnt].fid = 0;
					tiles[cnt].tid_h = i;
					tiles[cnt++].tid_w = 0;
				}
			}
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H-1){
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
			}
			*N = cnt;
		}
		if(tile.fid == 3){
			if(tile.tid_h == 0){
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = NO_TILE_H-1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = NO_TILE_H-2;
				tiles[cnt++].tid_w = 0;
				//1
				for(i=0; i < NO_TILE_H; i++){
					tiles[cnt].fid = 1;
					tiles[cnt].tid_h = i;
					tiles[cnt++].tid_w = 0;
				}
			}
			if(tile.tid_h == NO_TILE_H-1){
				//4
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = 0;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = NO_TILE_H-1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = NO_TILE_H-2;
				tiles[cnt++].tid_w = 0;
				//1
				for(i=0; i < NO_TILE_H; i++){
					tiles[cnt].fid = 0;
					tiles[cnt].tid_h = i;
					tiles[cnt++].tid_w = 0;
				}
			}
			if(tile.tid_h > 0 && tile.tid_h < NO_TILE_H-1){
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//4
				tiles[cnt].fid = 4;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = tile.tid_h;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = tile.tid_h - 1;
				tiles[cnt++].tid_w = 0;
				//5
				tiles[cnt].fid = 5;
				tiles[cnt].tid_h = tile.tid_h + 1;
				tiles[cnt++].tid_w = 0;
			}
			*N = cnt;
		}
	}
	ret = new int[*N];
	for(i=0; i < *N; i++){
		// printf("%d %d %d \n", tiles[i].fid,tiles[i].tid_w, tiles[i].tid_h);
		ret[i] = NO_TILE_H * NO_TILE_W * tiles[i].fid + NO_TILE_H * tiles[i].tid_w + tiles[i].tid_h;	
	}
	return ret;
}
int* getAdjacentTile(int NO_FACE, int NO_TILE_H, int NO_TILE_W, int tid, int *N){
	switch(NO_FACE){
		case 1:
		return getAdjacentTileERP(NO_TILE_H, NO_TILE_W, tid, N);
		case 2:
		return getAdjacentTileERP2(NO_TILE_H, NO_TILE_W, tid, N);
		break;
		case 6:
		return getAdjacentTileCMP(NO_TILE_H, NO_TILE_W, tid, N);
		break;		
	}
}
void showArray(int *arr, int N){
	int i;
	for(i=0; i < N; i++){
		printf("%d, ", arr[i]);
	}
	printf("\n");
}
void showVmask(int* vmask, int NO_FACE, int NO_TILE_H, int NO_TILE_W){
	int i,j,k;
	for(i=0; i < NO_FACE; i++){
		for(j=0; j < NO_TILE_W; j++){
			for(k=0; k < NO_TILE_H; k++){
				printf("%d ", vmask[i*NO_TILE_H*NO_TILE_W + j * NO_TILE_H + k]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
int* extVmask(int* vmask,int NO_FACE, int NO_TILE_H, int NO_TILE_W, int ext_width){
	int N = NO_TILE_W * NO_TILE_H; /* Number of tiles per face */
	int NO_TILE = N * NO_FACE;
	int* vmask_ext = new int[NO_TILE];
	int tid, tid2;
	int* adj;
	int M,i,a;
	/* init */
	for(tid = 0; tid < NO_TILE; tid ++){
		vmask_ext[tid] = vmask[tid];
	}
	//
	for(a = 1; a <= ext_width; a++){
		for(tid = 0; tid < NO_TILE; tid ++){
			/* For each visible tile */
			if(vmask_ext[tid] == a){
				adj = getAdjacentTile(NO_FACE, NO_TILE_H, NO_TILE_W, tid, &M);
				for(i=0; i < M; i++){
					if(vmask_ext[adj[i]] == 0)
						vmask_ext[adj[i]] = a + 1;
				}
			}
		}
	}
	return vmask_ext;
}
int* extVmask2(int* vmask,int NO_FACE, int NO_TILE_H, int NO_TILE_W, int ext_width){
  int No_tile = NO_FACE * NO_TILE_H * NO_TILE_W;
  int* ret = new int[No_tile];
  int M = NO_TILE_H;
  int N = NO_TILE_W;
  int i,j;
  struct Pos pos;
  if(ext_width == 0) return vmask;
  if(ext_width >= 3){
    for(i=0; i < No_tile; i++)
      ret[i] = 1;
    return ret;
  }
  // vmask after extension
  int** vmask_ext = init2dArrayInt(M, N);
  // convert 'vmask' to rectangle form
  // printf("#[extVmask2]: vmask\n");
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      vmask_ext[i][j] = vmask[i*M + j];
      // printf("%d ", vmask[i*M + j]);
    }
    // printf("\n");
  }
  //
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        // printf("#[extVmask2]: (i,j)=(%d,%d)\n", i,j);
        //TOP
        pos = TOP(i,j,M,N);
        // printf("#[extVmask2]: TOP(i,j)=(%d,%d)\n", pos.i,pos.j);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 2;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //DOWN
        pos = DOWN(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 3;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //LEFT
        pos = LEFT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 4;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //RIGHT
        pos = RIGHT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 5;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //TOP-LEFT
        pos = TOP_LEFT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 4;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //TOP-RIGHT
        pos = TOP_RIGHT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 5;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //DOWN-LEFT
        pos = DOWN_LEFT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 4;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] == 1){
        //DOWN-RIGHT
        pos = DOWN_RIGHT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 5;
      }
    }
  }
  // 2->1
  // printf("#[extVmask2]: vmask_ext\n");
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      // if(vmask_ext[i][j] == 2)
      // 	vmask_ext[i][j] = 1;
      ret[i*N + j]= vmask_ext[i][j];
      // printf("%d ", vmask_ext[i][j]);
    }
    // printf("\n"); 
  }
  if(ext_width == 1) return ret;
  // extension width = 2
  //
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        // printf("#[extVmask2]: (i,j)=(%d,%d)\n", i,j);
        //TOP
        pos = TOP(i,j,M,N);
        // printf("#[extVmask2]: TOP(i,j)=(%d,%d)\n", pos.i,pos.j);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 6;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //DOWN
        pos = DOWN(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 7;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //LEFT
        pos = LEFT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 8;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //RIGHT
        pos = RIGHT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 9;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //TOP-LEFT
        pos = TOP_LEFT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 8;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //TOP-RIGHT
        pos = TOP_RIGHT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 9;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //DOWN-LEFT
        pos = DOWN_LEFT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 8;
      }
    }
  }
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      if(vmask_ext[i][j] >= 1 && vmask_ext[i][j] <= 5){
        //DOWN-RIGHT
        pos = DOWN_RIGHT(i,j,M,N);
        if(pos.i != -1 && vmask_ext[pos.i][pos.j] == 0)
          vmask_ext[pos.i][pos.j] = 9;
      }
    }
  }
  // 2->1
  // printf("#[extVmask2]: vmask_ext_2\n");
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      // if(vmask_ext[i][j] == 2)
      // 	vmask_ext[i][j] = 1;
      ret[i*N + j]= vmask_ext[i][j];
      // printf("%d ", vmask_ext[i][j]);
    }
    // printf("\n"); 
  }
  return ret;
}
Tile* get_adj_tile_in_face(int NO_TILE_H, int NO_TILE_W, int tid_w, int tid_h, int *N){
	int MAX = 8, cnt = 0;
	Tile* adjTiles = new Tile[MAX];
	// top
	if(tid_w - 1 >= 0){
		// (w-1, h)
		adjTiles[cnt].tid_h = tid_h;
		adjTiles[cnt++].tid_w = tid_w - 1;
		// (w-1, h-1)
		if(tid_h - 1 >= 0){
			adjTiles[cnt].tid_h = (tid_h - 1);
			adjTiles[cnt++].tid_w = tid_w - 1;
		}
		// (w-1, h+1)
		if(tid_h + 1 < NO_TILE_H){
			adjTiles[cnt].tid_h = tid_h + 1;
			adjTiles[cnt++].tid_w = tid_w - 1;
		}
	}
	// bottom
	if(tid_w + 1 < NO_TILE_W){
		// (w+1, h)
		adjTiles[cnt].tid_h = tid_h;
		adjTiles[cnt++].tid_w = tid_w + 1;
		// (w+1, h-1)
		if(tid_h - 1 >= 0){
			adjTiles[cnt].tid_h = (tid_h - 1);
			adjTiles[cnt++].tid_w = tid_w + 1;
		}
			// (w+1, h+1)
		if(tid_h + 1 < NO_TILE_H){
			adjTiles[cnt].tid_h = tid_h + 1;
			adjTiles[cnt++].tid_w = tid_w + 1;
		}
	}
	// middle
	// (w, h-1)
	if(tid_h - 1 >= 0){
		adjTiles[cnt].tid_h = (tid_h - 1);
		adjTiles[cnt++].tid_w = tid_w;
	}
	// (w, h+1)
	if(tid_h + 1 < NO_TILE_H){
		adjTiles[cnt].tid_h = tid_h + 1;
		adjTiles[cnt++].tid_w = tid_w;
	}
	*N= cnt;
	return adjTiles;
}
struct Pos TOP(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = -1;
  pos.j = -1;
  if(i==0) return pos;
  pos.i = i-1;  
  pos.j = j;
  // printf("#[TOP] (pos.i, pos.j) = (%d, %d)\n", pos.i, pos.j);
  return pos;
}
struct Pos DOWN(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = -1;
  pos.j = -1;
  if(i==M-1) return pos;
  pos.j = j;
  pos.i = i+1;
  return pos;
}
struct Pos LEFT(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = i;
  if(j==0)
    pos.j = N-1;
  else
    pos.j = j-1;
  return pos;
}
struct Pos RIGHT(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = i;
  if(j==N-1)
    pos.j = 0;
  else
    pos.j = j+1;
  return pos;
}
struct Pos TOP_LEFT(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = -1;
  pos.j = -1;
  if(i==0)
    return pos;
  pos.i = i-1;
  pos.j = (j==0)?(N-1):(j-1);
  return pos;
}
struct Pos TOP_RIGHT(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = -1;
  pos.j = -1;
  if(i==0)
    return pos;
  pos.i = i-1;
  pos.j = (j==N-1)?0:(j+1);
  return pos;
}
struct Pos DOWN_LEFT(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = -1;
  pos.j = -1;
  if(i==M-1)
    return pos;
  pos.i = i+1;
  pos.j = (j==0)?(N-1):(j-1);
  return pos;
}
struct Pos DOWN_RIGHT(int i, int j, int M, int N){
  struct Pos pos;
  pos.i = -1;
  pos.j = -1;
  if(i==M-1)
    return pos;
  pos.i = i+1;
  pos.j = (j==N-1)?0:(j+1);
  return pos;
}