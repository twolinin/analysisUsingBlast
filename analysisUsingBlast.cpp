#include "analysisUsingBlast.h"

bool showDebug = false;

int main(int argc, char** argv) 
{

    AlignVec alignData;
    AlignVec resultData;
	
    int totalAlignContig = 0;
	int misassembled = 0;
	
    int alignLength[100000];
    int alignPoint = 0;
    size_t totalLength = 0;
	
	// read raw blast result
	if(!getBlastResult(argv[1],alignData))return 0;
	
    // filter duplicate
	resultData = filterDuplicate(alignData);

	totalAlignContig = resultData.size();
	
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
			{
				cout.width(13);
				cout << (*iter).contig;
				cout.width(13);
				cout << (*firPBIter).first;
				cout.width(13);
				cout <<(*firPBIter).second; 
				cout.width(13);
				cout << (*firRefIter).first; 
				cout.width(13);					
				cout << (*firRefIter).second ;
				cout.width(13);
            }
			*/
            alignLength[alignPoint] = (*alnter);
            totalLength += alignLength[alignPoint];
			
            if(showDebug)cout<< alignLength[alignPoint] << "\n";
			
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
		int continueMisassembled = 0;

        for(int i = 0 ; i < (*iter).PBVec.size() ; i++ )
        {
            if(showDebug)
			{
				cout.width(13);	
				cout << (*iter).contig;
				cout.width(13);	
				cout << positionArray[i][0];
				cout.width(13);	
				cout << positionArray[i][1];
				cout.width(13);			 
				cout << positionArray[i][2];
				cout.width(13);	
				cout << positionArray[i][3];
				cout.width(13);	
            }
            if( i == 0 )
            {
			    if(showDebug)cout << ( ( positionArray[i][2] > positionArray[i][3] ) ? "   <--" : "-->" ) << "\n";
                startPacbioPosition = positionArray[i][0];
				continue;
            }
			
			bool nowDir = positionArray[i][2]   > positionArray[i][3];
            int tmpArray[(*iter).PBVec.size()][4];
			int endBlastResult = i;
			
			if(showDebug) cout << ( nowDir ? "   <--" : "-->" ) << "\t" ;

			std::copy( positionArray, positionArray + (*iter).PBVec.size() , tmpArray);

		    float alignRatio = alignRegionRatio( tmpArray, (*iter).PBVec.size(), i, endBlastResult );
			
			if(showDebug) cout << "\n" <<  alignRatio << "\t" << abs( positionArray[i-1][2] - positionArray[endBlastResult-1][3] ) << "\t";
			
			if( alignRatio < 0.2 )
			{
				alignLength[alignPoint]  = abs( positionArray[i-1][2] - positionArray[endBlastResult-1][3] );
                alignPoint++;
			}
			else
			{
				alignLength[alignPoint]  = abs( positionArray[i-1][2] - positionArray[endBlastResult-1][3] );
                alignPoint++;
				misassembled++;
			}
			i = endBlastResult;
			
			totalLength += alignLength[alignPoint-1];
			
			if(showDebug) cout << totalLength << "\n";
        }
    }
    
	//cout << totalLength << "\n";
	
	sort(alignLength,alignLength+alignPoint);
    
	totalLength /= 2;
    
    //cout<< alignPoint << "\n";
    
	// NA50
    for(int i = alignPoint - 1 ; i >= 0; i-- )
    {
        if( alignLength[i] > totalLength )
        {
            cout << alignLength[i] << "\t";
            break;
        }
        totalLength -= alignLength[i];
    }
    
    cout 
	     //<< totalAlignContig << "\t" 
	     << misassembled     << "\n";
     
    return 0;
}


int Cmp(const void *lhs, const void *rhs) {
    return ((const int *)lhs)[1]-((const int *)rhs)[1];
}

bool getBlastResult(std::string blastAlignFile, AlignVec &result)
{
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
	
	ifstream inputAlignFile ( blastAlignFile.c_str() );
	
	if(!inputAlignFile.is_open())
    {
        std::cout << "can't open "<< blastAlignFile << "\n";
        return false;
    }
	
	while(inputAlignFile && inputAlignFile.is_open())
    {
         inputAlignFile >> queryName;
         inputAlignFile >> refName;
         inputAlignFile >> bitScore;
         inputAlignFile >> alignLen;
         inputAlignFile >> unknown1;
         inputAlignFile >> unknown2;
         inputAlignFile >> queryStart;
         inputAlignFile >> queryEnd;
         inputAlignFile >> refStart;
         inputAlignFile >> refEnd;
         inputAlignFile >> pValue;
         inputAlignFile >> unknown3;
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
        
        if( result.size() == 0 )
        {
            alignInfo tmp;
            tmp.contig = queryName;
            tmp.alignLenVec.push_back(alignLen);
            tmp.PBVec.push_back(make_pair(queryStart,queryEnd));
            tmp.refVec.push_back(make_pair(refStart,refEnd));
            result.push_back(tmp);  
        }
        else
        {
            bool exist = false;
            
            for(AlignVec::iterator iter = result.begin() ; iter != result.end() ; ++iter)
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
                result.push_back(tmp);
            }
        }
    }
	
	return true;
}

AlignVec filterDuplicate(AlignVec rawAlignData)
{
	AlignVec result;
	
	for(AlignVec::iterator iter = rawAlignData.begin() ; iter != rawAlignData.end() ; ++iter)
    {
		int max = 0;
        bool* bitArray;
        bool insert = false;
        
        RangeVec::iterator refIter = (*iter).refVec.begin();
        vector<int>::iterator alignIter = (*iter).alignLenVec.begin();
        
        if((*iter).PBVec.size()==1)
        {
			result.push_back((*iter));
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
                    result.push_back(tmp);
                    ++refIter;
                    insert = true;
                }
                else
                {
                    AlignVec::iterator tmp = result.end();
                    --tmp;
                    (*tmp).alignLenVec.push_back((*alignIter));
                    (*tmp).PBVec.push_back(make_pair((*posIter).first,(*posIter).second));
                    (*tmp).refVec.push_back(make_pair((*refIter).first,(*refIter).second));
                    ++refIter;
                    
                }
			}
        }
    }
	return result;
}

float alignRegionRatio(int positionArray[][4], int arraySize, int start, int &end)
{
	size_t loaclUnAlignLength = 1;
	size_t contigLength =  abs( positionArray[start-1][0] - positionArray[arraySize-1][1]) ;
	
	if(showDebug)
	{
		cout << "\n";
		cout.width(26);	
		cout << positionArray[start-1][0];
		cout.width(13);	
		cout << positionArray[start-1][1];
		cout.width(13);			 
		cout << positionArray[start-1][2];
		cout.width(13);	
		cout << positionArray[start-1][3];
		cout.width(13);	
		cout << ( (positionArray[start][2] > positionArray[start][3]) ? "   <--" : "-->" );
	}
	
	for(int i = start ; i < arraySize ; i++ )
	{
		int dis = 0;
		
		for(int j = i ; j < arraySize && dis < 5; j++, dis++ )
		{
			// check same region by calculate dissimilar ratio
			int   local_pbLen      = abs( positionArray[start-1][0] - positionArray[j][1] );
			int   local_refLen     = abs( positionArray[start-1][2] - positionArray[j][3] );
			float local_dissimilar = (float)((int)( abs( 1 - (float)local_pbLen/(float)local_refLen )*1000 ))/1000;
			// check same region by two fragment direction
			bool  local_beforeDir  = positionArray[start-1][2] > positionArray[start-1][3];
			bool  local_nowDir     = positionArray[j][2]       > positionArray[j][3];
			
			int   local_pbGapLen   = abs( positionArray[i-1][1] - positionArray[j][0] );
			bool  local_connect    = ( positionArray[i-1][3] < 3000 || positionArray[j][2] < 3000 ) || 
									 ( abs( positionArray[i-1][3] - positionArray[j][2]) < 7000  + local_pbGapLen ) && ( abs( positionArray[i-1][1] - positionArray[j][0]) < 7000  + local_pbGapLen );
			
			if( (float)local_pbGapLen/(float)max(positionArray[i-1][1],positionArray[j][0]) > 0.2 ) local_connect = false;
			
			// check contain
			if( positionArray[i-1][2] >= positionArray[j][2] && positionArray[i-1][3] <= positionArray[j][3] ) local_connect = false;
			if( positionArray[j][2] >= positionArray[i-1][2] && positionArray[j][3] <= positionArray[i-1][3] ) local_connect = false;
			// regard two alignment result as one
			
			if(showDebug)
			{
				cout << "\n";
				cout.width(26);	
				cout << positionArray[j][0];
				cout.width(13);	
				cout << positionArray[j][1];
				cout.width(13);			 
				cout << positionArray[j][2];
				cout.width(13);	
				cout << positionArray[j][3];
				cout.width(13);	
				cout << ( local_nowDir ? "   <--" : "-->" );
				cout.width(13);	
				cout << local_dissimilar;
				cout.width(10);	
				cout << positionArray[i-1][1];
				cout.width(10);	
				cout << positionArray[i-1][3];
			}
			
			if ( local_beforeDir == local_nowDir && local_dissimilar < 0.2 && local_connect )
			{	
				i += dis;
				end = i;
				
				if(showDebug)
				{
					cout.width(5);	
					cout << dis;
				}
				
				positionArray[j][0] = positionArray[j-1][0];
				//positionArray[j][1] = positionArray[j-1][1];
				positionArray[j][2] = positionArray[j-1][2];
				//positionArray[j][3] = positionArray[j-1][3];

				break;
			}
			
			loaclUnAlignLength += abs( positionArray[j][0] - positionArray[j][1] );
			
			positionArray[j][0] = positionArray[j-1][0];
			positionArray[j][1] = positionArray[j-1][1];
			positionArray[j][2] = positionArray[j-1][2];
			positionArray[j][3] = positionArray[j-1][3];
			
			
			if( j == arraySize-1 ) return (float)loaclUnAlignLength/(float)contigLength;
		}
		if( dis == 5 ) break;
	}

	return (float)loaclUnAlignLength/(float)contigLength;
};