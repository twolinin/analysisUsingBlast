#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h> 
#include <algorithm> 


using namespace std;

struct alignInfo;

typedef pair<int, int> Range;
typedef vector<Range> RangeVec;
typedef vector<alignInfo> AlignVec;

struct alignInfo
{
    string contig;
    vector<int> alignLenVec;
    RangeVec contigVec;
    RangeVec refVec;
};

void caculateNA50(int argc, char** argv);
void chimeraDetected(int argc, char** argv);

int Cmp(const void *lhs, const void *rhs);

bool getBlastResult(std::string blastAlignFile, AlignVec &result);

AlignVec filterDuplicate(AlignVec rawAlignData);
	
float alignRegionRatio(int positionArray[][4], int arraySize, int start, int &end);