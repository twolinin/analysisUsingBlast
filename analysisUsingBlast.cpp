#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h> 
#include <algorithm> 

using namespace std;

typedef pair<int, int> Range;
typedef vector<Range> RangeVec;

struct alignInfo
{
    string contig;
    vector<int> alignLenVec;
    RangeVec PBVec;
    RangeVec refVec;
};

typedef vector<alignInfo> AlignVec;

int Cmp(const void *lhs, const void *rhs) {
    return ((const int *)lhs)[1]-((const int *)rhs)[1];
}

int main(int argc, char** argv) 
{
    std::string command(argv[1]);
    
    string queryName;
    string refName;
    float bitScore;
    int alignLen;
    string unknown1;
    string unknown2;
    int queryStart;
    int queryEnd;
    int refStart;
    int refEnd;
    float pValue;
    string unknown3;
    
    AlignVec alignData;
    AlignVec resultData;
    
    ifstream inputFile( command.c_str() );
    
    if(!inputFile.is_open())
    {
        std::cout << "can't open "<< command << "\n";
        return 0;
    }
    
    while(inputFile && inputFile.is_open())
    {
         inputFile >> queryName;
         inputFile >> refName;
         inputFile >> bitScore;
         inputFile >> alignLen;
         inputFile >> unknown1;
         inputFile >> unknown2;
         inputFile >> queryStart;
         inputFile >> queryEnd;
         inputFile >> refStart;
         inputFile >> refEnd;
         inputFile >> pValue;
         inputFile >> unknown3;
        /*
         std::cout<< queryName  << "\t" 
                  << alignLen   << "\t"
                  << queryStart << "\t" 
                  << queryEnd   << "\t" 
                  << refStart   << "\t"
                  << refEnd     << "\n";
        */    
        /*
        std::cout<< queryName  << "\t"
                 << refName    << "\t"
                 << bitScore   << "\t"
                 << alignLen   << "\t"
                 << unknown1   << "\t"
                 << unknown2   << "\t"
                 << queryStart << "\t"
                 << queryEnd   << "\t"
                 << refStart   << "\t"
                 << refEnd     << "\t"
                 << pValue     << "\t"
                 << unknown3   << "\t"
                 << "\n";
        */
        if( pValue != 0 )continue;
        
        if( alignData.size() == 0 )
        {
            alignInfo tmp;
            tmp.contig = queryName;
            tmp.alignLenVec.push_back(alignLen);
            tmp.PBVec.push_back(make_pair(queryStart,queryEnd));
            tmp.refVec.push_back(make_pair(refStart,refEnd));
            alignData.push_back(tmp);
            /*
            std::cout<< queryName  << "\t"
                 << refName    << "\t"
                 << bitScore   << "\t"
                 << alignLen   << "\t"
                 << unknown1   << "\t"
                 << unknown2   << "\t"
                 << queryStart << "\t"
                 << queryEnd   << "\t"
                 << refStart   << "\t"
                 << refEnd     << "\t"
                 << pValue     << "\t"
                 << unknown3   << "\t"
                 << "\n";
                 */
        }
        else
        {
            bool exist = false;
            
            for(AlignVec::iterator iter = alignData.begin() ; iter != alignData.end() ; ++iter)
            {
                if( queryName == (*iter).contig )
                {
                    exist = true;
                    
                    bool contain = false;
                    
                    for(RangeVec::iterator posIter = (*iter).PBVec.begin() ; posIter != (*iter).PBVec.end() ; ++posIter)
                    {
                        // --------------------------------------
                        //        -------------------------
                        if( (*posIter).first <= queryStart && (*posIter).second >= queryEnd )
                        {
                            contain = true; 
                            break;
                        }
                        //        -------------------------------
                        //  ---------------------
                        if( (*posIter).first >= queryStart && (*posIter).second >= queryEnd  && queryEnd >= (*posIter).first )
                        {
                            if( (float)(queryEnd - (*posIter).first) / (float)(queryEnd - queryStart + 1) > 0.7 ) contain = true;
                            if( (*posIter).first - queryStart <= 1000 )
                            {
                                contain = true;
                                break;
                            }
                        }
                        //        -------------------------------
                        //                            ---------------------
                        if( (*posIter).first <= queryStart && (*posIter).second <= queryEnd && (*posIter).second >= queryStart)
                        {
                            if( (float)((*posIter).second - queryStart) / (float)(queryEnd - queryStart + 1) > 0.7 ) contain = true;
                            if( queryEnd - (*posIter).second <= 1000 )
                            {
                                contain = true;
                                break;
                            }
                        }
                    }
                    
                    if( !contain && alignLen >= 500 )
                    {
                        (*iter).alignLenVec.push_back(alignLen);
                        (*iter).PBVec.push_back(make_pair(queryStart,queryEnd));
                        (*iter).refVec.push_back(make_pair(refStart,refEnd));
                        /*
                        std::cout<< queryName  << "\t"
                                 << refName    << "\t"
                                 << bitScore   << "\t"
                                 << alignLen   << "\t"
                                 << unknown1   << "\t"
                                 << unknown2   << "\t"
                                 << queryStart << "\t"
                                 << queryEnd   << "\t"
                                 << refStart   << "\t"
                                 << refEnd     << "\t"
                                 << pValue     << "\t"
                                 << unknown3   << "\t"
                                 << "\n";
                                 */
                    }
                }
            }
            
            if( !exist )
            {
                alignInfo tmp;
                tmp.contig = queryName;
                tmp.alignLenVec.push_back(alignLen);
                tmp.PBVec.push_back(make_pair(queryStart,queryEnd));
                tmp.refVec.push_back(make_pair(refStart,refEnd));
                alignData.push_back(tmp);
                /*
                std::cout<< queryName  << "\t"
                 << refName    << "\t"
                 << bitScore   << "\t"
                 << alignLen   << "\t"
                 << unknown1   << "\t"
                 << unknown2   << "\t"
                 << queryStart << "\t"
                 << queryEnd   << "\t"
                 << refStart   << "\t"
                 << refEnd     << "\t"
                 << pValue     << "\t"
                 << unknown3   << "\t"
                 << "\n";
                 */
            }
        }
    }
    
    int misassembled = 0;
	int continueMisassembled = 0;
    
    // filter duplicate
    for(AlignVec::iterator iter = alignData.begin() ; iter != alignData.end() ; ++iter)
    {
        int max = 0;
        bool* bitArray;
        bool insert = false;
        
        RangeVec::iterator refIter = (*iter).refVec.begin();
        vector<int>::iterator alignIter = (*iter).alignLenVec.begin();
        
        if((*iter).PBVec.size()==1)
        {
            resultData.push_back((*iter));
            continue;
        }
        
        for(RangeVec::iterator posIter = (*iter).PBVec.begin() ; posIter != (*iter).PBVec.end() ; ++posIter)
            if( (*posIter).second > max ) max = (*posIter).second;
        
        bitArray = new bool[max];
        
        for( int i = 0 ; i < max ; i++ ) 
            bitArray[i] = false;
        
        for(RangeVec::iterator posIter = (*iter).PBVec.begin() ; posIter != (*iter).PBVec.end() ; ++posIter)
        {
            int overlap = 0;
            int noOverlap = 0;
            bool filter = false;
            
            for( int i = (*posIter).first ; i < (*posIter).second ; i++ )
            {
                if( !bitArray[i] )
                {
                    bitArray[i] = true;
                    noOverlap++;
                }
                else overlap++;
            }
            
            if( overlap == 0 ) 
            {
                if(!insert)
                {
                    alignInfo tmp;
                    tmp.contig = (*iter).contig;
                    tmp.alignLenVec.push_back((*alignIter));
                    tmp.PBVec.push_back(make_pair((*posIter).first,(*posIter).second));
                    tmp.refVec.push_back(make_pair((*refIter).first,(*refIter).second));
                    resultData.push_back(tmp);
                    ++refIter;
                    insert = true;
                }
                else
                {
                    AlignVec::iterator tmp = resultData.end();
                    --tmp;
                    (*tmp).alignLenVec.push_back((*alignIter));
                    (*tmp).PBVec.push_back(make_pair((*posIter).first,(*posIter).second));
                    (*tmp).refVec.push_back(make_pair((*refIter).first,(*refIter).second));
                    ++refIter;
                    
                }
            }
            else if( (float)noOverlap / (float)(overlap + noOverlap ) > 0.7 )
            {

                if(!insert)
                {
                    alignInfo tmp;
                    tmp.contig = (*iter).contig;
                    tmp.alignLenVec.push_back((*alignIter));
                    tmp.PBVec.push_back(make_pair((*posIter).first,(*posIter).second));
                    tmp.refVec.push_back(make_pair((*refIter).first,(*refIter).second));
                    resultData.push_back(tmp);
                    ++refIter;
                    insert = true;

                }
                else
                {
                    AlignVec::iterator tmp = resultData.end();
                    --tmp;
                    (*tmp).alignLenVec.push_back((*alignIter));
                    (*tmp).PBVec.push_back(make_pair((*posIter).first,(*posIter).second));
                    (*tmp).refVec.push_back(make_pair((*refIter).first,(*refIter).second));
                    ++refIter;
                }
            }
        }
    }
    
    int deviation = 1000;
    int alignLength[100000];
    int alignPoint = 0;
    size_t totalLength = 0;

    for(AlignVec::iterator iter = resultData.begin() ; iter != resultData.end() ; ++iter)
    {
        vector<int>::iterator alnter  = (*iter).alignLenVec.begin(); // record align length
        RangeVec::iterator firPBIter  = (*iter).PBVec.begin();  //record pb align start and end
        RangeVec::iterator firRefIter = (*iter).refVec.begin(); //record ref align start and end
        int positionArray[(*iter).PBVec.size()][4];
        
        if((*iter).PBVec.size()==1)
        {
            /*
            cout<< (*iter).contig       << "\t"
                << (*firPBIter).first   << "\t" 
                << (*firPBIter).second  << "\t" 
                << (*firRefIter).first  << "\t" 
                << (*firRefIter).second << "\t" ;
            */
            alignLength[alignPoint] = (*alnter);
            totalLength += alignLength[alignPoint];
            //cout<< alignLength[alignPoint] << "\n";
            alignPoint++;
            continue;
        }
        
        for(int i = 0 ; i < (*iter).PBVec.size() ; i++ )
        {
            positionArray[i][0] = (*firPBIter).first;
            positionArray[i][1] = (*firPBIter).second;
            positionArray[i][2] = (*firRefIter).first;
            positionArray[i][3] = (*firRefIter).second;
            
            firPBIter++;
            firRefIter++;
        }

        qsort(positionArray,(*iter).PBVec.size(),sizeof(positionArray[0]),Cmp);
        
        //1       25949   739574  765566
        //29479   67738   768163  806592
        //67736   114370  800893  847667
        //
        int startPacbioPosition = 0;
        int referenceAccumulationLen = 0;
		int continueMisassembled = 0;
		
        for(int i = 0 ; i < (*iter).PBVec.size() ; i++ )
        {
            /*
            std::cout << (*iter).contig      << "\t"
                      << positionArray[i][0] << "\t"
                      << positionArray[i][1] << "\t"
                      << positionArray[i][2] << "\t"
                      << positionArray[i][3] << "\t"
                      << "\t";
            */
            if( i == 0 )
            {
                bool  nowDir     = positionArray[i][2]   > positionArray[i][3];
                //cout << ( nowDir ? "   <--" : "-->" ) << "\n";
                
                startPacbioPosition      = positionArray[i][0];
                referenceAccumulationLen = abs( positionArray[i][2] - positionArray[i][3] );
                
				continue;
            }
            
            
            int   pbLen      = abs( startPacbioPosition - positionArray[i][1] );
            int   refLen     = abs( positionArray[i][2] - positionArray[i][3] ) + referenceAccumulationLen;
            float dissimilar = (float)((int)( abs( 1 - (float)pbLen/(float)refLen )*1000 ))/1000;
            bool  beforeDir  = positionArray[i-1][2] > positionArray[i-1][3];
            bool  nowDir     = positionArray[i][2]   > positionArray[i][3];
            bool  connect    = ( positionArray[i-1][3] < 10 || positionArray[i][2] < 10 ) || 
                               ( abs( 1 - (float)positionArray[i-1][3]/(float)positionArray[i][2]) < 0.2 );

            //cout << ( nowDir ? "   <--" : "-->" ) << "\t" << dissimilar  << "\t" ;
                               
            if ( beforeDir == nowDir && dissimilar < 0.2 && connect)
            {
                //cout << "\n";
                
				continueMisassembled--;
				if( continueMisassembled < 0 ) continueMisassembled = 0;
				
                referenceAccumulationLen = refLen;
                continue;
            }
            else
            {
                //cout << referenceAccumulationLen << "\t";
                
                totalLength += referenceAccumulationLen;
                
                alignLength[alignPoint]  = referenceAccumulationLen;
                startPacbioPosition      = positionArray[i][0];
                referenceAccumulationLen = abs( positionArray[i][2] - positionArray[i][3] );
                
                alignPoint++;
                
				if( continueMisassembled == 0 ) misassembled++;
				
				continueMisassembled = 3;
				
				//cout << misassembled << "\t" << "\t" << "\n";
				
            }

        }
        
        //prevent miniasm contig blast result, it exist many local align 

        
        //cout << referenceAccumulationLen << "\t" << misassembled << "\t" << "\n\n";
        alignLength[alignPoint] = referenceAccumulationLen;
        totalLength += referenceAccumulationLen;
        alignPoint++;
    }
    
    sort(alignLength,alignLength+alignPoint);
    totalLength /= 2;
    
    //cout<< alignPoint << "\n";
    /*
    cout << alignLength[alignPoint/10] << "\t"
         << alignLength[alignPoint*2/10] << "\t"
         << alignLength[alignPoint*3/10] << "\t"
         << alignLength[alignPoint*4/10] << "\t"
         << alignLength[alignPoint*5/10] << "\t"
         << alignLength[alignPoint*6/10] << "\t"
         << alignLength[alignPoint*7/10] << "\t"
         << alignLength[alignPoint*8/10] << "\t"
         << alignLength[alignPoint*9/10] << "\t"
         << alignLength[alignPoint-1] << "\n";
    */
    
    for(int i = alignPoint - 1 ; i >= 0; i-- )
    {
        if( alignLength[i] > totalLength )
        {
            cout << alignLength[i] << "\t";
            break;
        }
        totalLength -= alignLength[i];
    }
    
    cout << misassembled << "\t" << "\n";
     
    return 0;
}