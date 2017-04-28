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
    int totalAlignContig = 0;
	
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
        if( queryName == "tig00000142")
		
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
        if( alignLen < 1000 )continue;
        
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
                            // unique length
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
                            // unique length
							if( queryEnd - (*posIter).second <= 1000 )
                            {
                                contain = true;
                                break;
                            }
                        }
                    }
                    
                    if( !contain && alignLen >= 1000 )
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
    int misassembledContig = 0;
	
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
        
        for( int i = 0 ; i < max ; i++ ) bitArray[i] = false;
        
        for(RangeVec::iterator posIter = (*iter).PBVec.begin() ; posIter != (*iter).PBVec.end() ; ++posIter)
        {
			int overlap   = 0;
            int noOverlap = 0;
            bool filter   = true;

            for( int i = (*posIter).first ; i < (*posIter).second ; i++ )
            {
                if( !bitArray[i] )
                {
                    bitArray[i] = true;
                    noOverlap++;
                }
                else overlap++;
            }
			
			if( overlap == 0 ) filter = false;
			else if( (float)noOverlap / (float)(overlap + noOverlap ) > 0.5  || noOverlap > 5000 ) filter = false;
			
			if( !filter )
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
    
    int alignLength[100000];
    int alignPoint = 0;
    size_t totalLength = 0;
	
	totalAlignContig = resultData.size();
	
	bool showDebug = true;
	
    for(AlignVec::iterator iter = resultData.begin() ; iter != resultData.end() ; ++iter)
    {
        vector<int>::iterator alnter  = (*iter).alignLenVec.begin(); // record align length
        RangeVec::iterator firPBIter  = (*iter).PBVec.begin();  //record pb align start and end
        RangeVec::iterator firRefIter = (*iter).refVec.begin(); //record ref align start and end
        int positionArray[(*iter).PBVec.size()][4];
		
        if((*iter).PBVec.size()==1)
        {
            /*
			if(showDebug)
				cout<< (*iter).contig       << "\t"
					<< (*firPBIter).first   << "\t"  << (*firPBIter).second  << "\t" 
					<< (*firRefIter).first  << "\t"  << (*firRefIter).second << "\t" ;
            */
            alignLength[alignPoint] = (*alnter);
            totalLength += alignLength[alignPoint];
			
            //if(showDebug)cout<< alignLength[alignPoint] << "\n";
			
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
        
        int startPacbioPosition = 0;
        int referenceAccumulationLen = 0;
		int continueMisassembled = 0;
		bool misassembledContigStatus = false;
		
        for(int i = 0 ; i < (*iter).PBVec.size() ; i++ )
        {
            if(showDebug)
				std::cout << (*iter).contig      << "\t"
						  << positionArray[i][0] << "\t" << positionArray[i][1] << "\t"
						  << positionArray[i][2] << "\t" << positionArray[i][3] << "\t"
						  << "\t";
            
            if( i == 0 )
            {
                bool  nowDir     = positionArray[i][2]   > positionArray[i][3];
				
                if(showDebug)cout << ( nowDir ? "   <--" : "-->" ) << "\n";
                
                startPacbioPosition      = positionArray[i][0];
                referenceAccumulationLen = abs( positionArray[i][2] - positionArray[i][3] );
                
				continue;
            }
            
            int   pbLen      = abs( startPacbioPosition - positionArray[i][1] );
			int   pbGapLen   = abs( positionArray[i-1][1] - positionArray[i][0] );
            int   refLen     = abs( positionArray[i][2] - positionArray[i][3] ) + referenceAccumulationLen;
            float dissimilar = (float)((int)( abs( 1 - (float)pbLen/(float)refLen )*1000 ))/1000;
            bool  beforeDir  = positionArray[i-1][2] > positionArray[i-1][3];
            bool  nowDir     = positionArray[i][2]   > positionArray[i][3];
            bool  connect    = ( positionArray[i-1][3] < 3000 || positionArray[i][2] < 3000 ) || 
							   ( abs( positionArray[i-1][3] - positionArray[i][2]) < 3000  + pbGapLen ) && ( abs( positionArray[i-1][1] - positionArray[i][0]) < 3000  + pbGapLen );
                               //( abs( 1 - (float)positionArray[i-1][3]/(float)positionArray[i][2]) < 0.2 );
			
            if( positionArray[i-1][2] >= positionArray[i][2] && positionArray[i-1][3] <= positionArray[i][3] ) connect = false;
			if( positionArray[i][2] >= positionArray[i-1][2] && positionArray[i][3] <= positionArray[i-1][3] ) connect = false;
			
			
			if ( !(beforeDir == nowDir && dissimilar < 0.2 && connect) )
			{
				if( i+1 < (*iter).PBVec.size() )
				{
					int   local_pbLen      = abs( startPacbioPosition - positionArray[i+1][1] );
					int   local_pbGapLen   = abs( positionArray[i-1][1] - positionArray[i+1][0] );
					int   local_refLen     = abs( positionArray[i+1][2] - positionArray[i+1][3] ) + referenceAccumulationLen;
					float local_dissimilar = (float)((int)( abs( 1 - (float)pbLen/(float)refLen )*1000 ))/1000;
					bool  local_beforeDir  = positionArray[i-1][2] > positionArray[i-1][3];
					bool  local_nowDir     = positionArray[i+1][2]   > positionArray[i+1][3];
					bool  local_connect    = ( positionArray[i-1][3] < 3000 || positionArray[i+1][2] < 3000 ) || 
									         ( abs( positionArray[i-1][3] - positionArray[i+1][2]) < 3000  + pbGapLen ) && ( abs( positionArray[i-1][1] - positionArray[i+1][0]) < 3000  + pbGapLen );
					
					if( positionArray[i-1][2] >= positionArray[i+1][2] && positionArray[i-1][3] <= positionArray[i+1][3] ) local_connect = false;
			        if( positionArray[i+1][2] >= positionArray[i-1][2] && positionArray[i+1][3] <= positionArray[i-1][3] ) local_connect = false;

					if ( local_beforeDir == local_nowDir && local_dissimilar < 0.2 && local_connect )
					{
						cout << "second chance\n";
						positionArray[i][0] = positionArray[i-1][0];
						positionArray[i][1] = positionArray[i-1][1];
						positionArray[i][2] = positionArray[i-1][2];
						positionArray[i][3] = positionArray[i-1][3];
						
						continue;
					}
					
				}
			}
			
			
			
			
            if(showDebug)
				cout << ( nowDir ? "   <--" : "-->" ) << "\t" 
				     << positionArray[i-1][1] - positionArray[i][0] << "\t"
					 << positionArray[i-1][3] - positionArray[i][2] << "\t"
					 << dissimilar  << "\t" ;
                               
            if ( beforeDir == nowDir && dissimilar < 0.2 && connect)
            {
                if(showDebug)cout << "\n";
                
				continueMisassembled--;
				if( continueMisassembled < 0 ) continueMisassembled = 0;
				
                referenceAccumulationLen = refLen;
                continue;
            }
            else
            {
                if(showDebug)cout << referenceAccumulationLen << "\t";
                
                totalLength += referenceAccumulationLen;
                
                alignLength[alignPoint]  = referenceAccumulationLen;
                startPacbioPosition      = positionArray[i][0];
                referenceAccumulationLen = abs( positionArray[i][2] - positionArray[i][3] );
                
                alignPoint++;
				
                //prevent miniasm contig blast result, it exist many local align 
				if( continueMisassembled == 0 ) misassembled++;
				
				continueMisassembled = 3;
				
				misassembledContigStatus = true;
				
				if(showDebug)cout << misassembled << "\t" << "\t" << "\n";
				
            }

        }
        
		if( misassembledContigStatus ) 
		{
			//cout << (*iter).contig << "\n";
			misassembledContig++;
        }
        if(showDebug)cout << referenceAccumulationLen << "\t" << misassembled << "\t" << "\n\n";
        alignLength[alignPoint] = referenceAccumulationLen;
        totalLength += referenceAccumulationLen;
        alignPoint++;
    }
    
    sort(alignLength,alignLength+alignPoint);
    totalLength /= 2;
    
    //cout<< alignPoint << "\n";
    
    for(int i = alignPoint - 1 ; i >= 0; i-- )
    {
        if( alignLength[i] > totalLength )
        {
            cout << alignLength[i] << "\t";
            break;
        }
        totalLength -= alignLength[i];
    }
    
    cout << totalAlignContig << "\t" << misassembled << "\t" << misassembledContig << "\n";
     
    return 0;
}