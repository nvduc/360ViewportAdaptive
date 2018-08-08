struct Tile
{
	int fid; /* face id */
	int tid_h; /* Tile id in horizontal axis */
	int tid_w; /* Tile id in vertical axis */
};
int* getAdjacentTileCMP(int NO_TILE_H, int NO_TILE_W, int tid, int* N);
int* getAdjacentTileERP2(int NO_TILE_H, int NO_TILE_W, int tid, int* N);
int* getAdjacentTileERP(int NO_TILE_H, int NO_TILE_W, int tid, int* N);
int* getAdjacentTile(int NO_FACE, int NO_TILE_H, int NO_TILE_W, int tid, int *N);
int* extVmask(int* vmask,int NO_FACE, int NO_TILE_H, int NO_TILE_W, int ext_width);
int* extVmask2(int* vmask,int NO_FACE, int NO_TILE_H, int NO_TILE_W, int ext_width);
void showVmask(int* vmask, int NO_FACE, int NO_TILE_H, int NO_TILE_W);
void showArray(int *arr, int N);
Tile* get_adj_tile_in_face(int NO_TILE_H, int NO_TILE_W, int tid_w, int tid_h, int* N);
// tile position structure
struct Pos {
  int i;
  int j;
};
// list of adjecent tiles of tile (i,j)
struct Pos TOP(int i, int j, int M, int N);
struct Pos DOWN(int i, int j, int M, int N);
struct Pos LEFT(int i, int j, int M, int N);
struct Pos RIGHT(int i, int j, int M, int N);
struct Pos TOP_LEFT(int i, int j, int M, int N);
struct Pos TOP_RIGHT(int i, int j, int M, int N);
struct Pos DOWN_LEFT(int i, int j, int M, int N);
struct Pos DOWN_RIGHT(int i, int j, int M, int N);