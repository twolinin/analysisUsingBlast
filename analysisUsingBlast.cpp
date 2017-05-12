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
        RangeVec::iterator firPBIter  = (*iter).contigVec.begin();  //record contig align start and end
        RangeVec::iterator firRefIter = (*iter).refVec.begin(); //record ref align start and end
        int positionArray[(*iter).contigVec.size()][4];
        
        if((*iter).contigVec.size()==1)
        { 
            /*if(showDebug)
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
            }*/
            
            alignLength[alignPoint] = (*alnter);
            totalLength += alignLength[alignPoint];
            
            //if(showDebug)cout<< alignLength[alignPoint] << "\n";
            
            alignPoint++;
            continue;
        }

        for(int i = 0 ; i < (*iter).contigVec.size() ; i++ )
        {
            positionArray[i][0] = (*firPBIter).first;
            positionArray[i][1] = (*firPBIter).second;
            positionArray[i][2] = (*firRefIter).first;
            positionArray[i][3] = (*firRefIter).second;
            
            firPBIter++;
            firRefIter++;    
        }
        
        qsort(positionArray,(*iter).contigVec.size(),sizeof(positionArray[0]),Cmp);
        
        totalLength += abs( positionArray[0][0] - positionArray[(*iter).contigVec.size()-1][1]) + 1;
        
        if(showDebug)
        {
            cout<< "\n";
            for(int i = 0 ; i < (*iter).contigVec.size() ; i++ )
            {
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
                cout << ( (positionArray[i][2] > positionArray[i][3]) ? "   <--" : "-->" );
                cout << "\n";
                
            }
            cout<< "\n";
        }
        
        if(showDebug)cout << (*iter).contig << "\n";
        
        int continueMisassembled = 0;
        
        
        for(int i = 0 ; i < (*iter).contigVec.size() ; i++ )
        {
            int tmpArray[(*iter).contigVec.size()][4];
            int endBlastResult = i;
            
            std::copy( positionArray, positionArray + (*iter).contigVec.size() , tmpArray);
 
            float alignRatio = alignRegionRatio( tmpArray, (*iter).contigVec.size(), i, endBlastResult );

            alignLength[alignPoint] = abs( positionArray[i][0] - positionArray[endBlastResult][1]) + 1;
            alignPoint++;
            
            i = endBlastResult;
    
            if( alignRatio > 0.2 && i != (*iter).contigVec.size() - 1 ) 
            {
                misassembled++;
            }
            

            if(showDebug) cout  << ((float)((int)(alignRatio*1000 )))/1000 << "\t" << alignLength[alignPoint-1] << "\t" << misassembled << "\n";
            
        }
        if(showDebug)cout << "\n";
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


int Cmp(const void *lhs, const void *rhs)
{
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
            tmp.contigVec.push_back(make_pair(queryStart,queryEnd));
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
                    
                    for(RangeVec::iterator posIter = (*iter).contigVec.begin() ; posIter != (*iter).contigVec.end() ; ++posIter)
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
                        (*iter).contigVec.push_back(make_pair(queryStart,queryEnd));
                        (*iter).refVec.push_back(make_pair(refStart,refEnd));      
                    }
                }
            }
            
            if( !exist )
            {
                alignInfo tmp;
                tmp.contig = queryName;
                tmp.alignLenVec.push_back(alignLen);
                tmp.contigVec.push_back(make_pair(queryStart,queryEnd));
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
        
        if((*iter).contigVec.size()==1)
        {
            result.push_back((*iter));
            continue;
        }

        for(RangeVec::iterator posIter = (*iter).contigVec.begin() ; posIter != (*iter).contigVec.end() ; ++posIter)
            if( (*posIter).second > max ) max = (*posIter).second;
        
        bitArray = new bool[max];
        
        for( int i = 0 ; i < max ; i++ ) bitArray[i] = false;
        
        for(RangeVec::iterator posIter = (*iter).contigVec.begin() ; posIter != (*iter).contigVec.end() ; ++posIter)
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
                    tmp.contigVec.push_back(make_pair((*posIter).first,(*posIter).second));
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
                    (*tmp).contigVec.push_back(make_pair((*posIter).first,(*posIter).second));
                    (*tmp).refVec.push_back(make_pair((*refIter).first,(*refIter).second));
                    ++refIter;
                    
                }
            }
        }
    }
    return result;
}

float alignRegionRatio(int positionArray[][4], int arraySize, int start, int &end )
{
    size_t loaclUnAlignLength = 1;
    size_t contigLength =  abs( positionArray[start][0] - positionArray[arraySize-1][1]);
    
    bool  local_connect = true;

    if(showDebug)
    {
        cout.width(26);    
        cout << positionArray[start][0];
        cout.width(13);    
        cout << positionArray[start][1];
        cout.width(13);             
        cout << positionArray[start][2];
        cout.width(13);    
        cout << positionArray[start][3];
        cout.width(13);    
        cout << ( (positionArray[start][2] > positionArray[start][3]) ? "   <--" : "-->" );
        cout << "\n";
    }
    
    // use the parameter of distance for preventing from checking all results
    for(int target = start + 1 , distance = 1 ; target < arraySize && distance <= 10 ; target++, distance++ )
    {
        // check same region by two fragment direction
        bool  local_beforeDir  = positionArray[start][2]  > positionArray[start][3];
        bool  local_nowDir     = positionArray[target][2] > positionArray[target][3];
		
		// check same region by calculate dissimilar ratio
        int   local_contitLen  = abs( positionArray[start][0] - positionArray[target][1] );
        int   local_refLen     = abs( positionArray[start][2] - positionArray[target][3] );
        float local_disSimilar = (float)((int)( abs( 1 - (float)local_contitLen/(float)local_refLen )*1000 ))/1000;
        
		// contig gap size
        int   local_contigGapLen = abs( positionArray[start][1] - positionArray[target][0] );
        int   local_refGapLen    = abs( positionArray[start][3] - positionArray[target][2] );
        
        // check gap size by contig position
        if( local_refGapLen    > 1000  ) local_connect = false;
        // check gap size by reference position
        if( local_contigGapLen > 1000  ) local_connect = false;
        // prevent circular align result
        if(( positionArray[start][3] < 1000 || positionArray[target][2] < 1000 ) ) local_connect = true;
        // prevent large gap size
        if( (float)local_contigGapLen/(float)max(positionArray[start][1],positionArray[target][0]) > 0.2 ) local_connect = false;
        // check contain
        if( positionArray[start][2] >= positionArray[target][2] && positionArray[start][3] <= positionArray[target][3] ) local_connect = false;
        if( positionArray[target][2] >= positionArray[start][2] && positionArray[target][3] <= positionArray[start][3] ) local_connect = false;
        
        if(showDebug)
        {
            cout.width(26);    
            cout << positionArray[target][0];
            cout.width(13);    
            cout << positionArray[target][1];
            cout.width(13);             
            cout << positionArray[target][2];
            cout.width(13);    
            cout << positionArray[target][3];
            cout.width(13);    
            cout << ( local_nowDir ? "   <--" : "-->" );
            cout.width(13);    
            cout << local_disSimilar;
            cout.width(13);    
            cout << (local_connect ? "O":"X") ;
        }
        
        // regard two alignment result as one
        if ( local_beforeDir == local_nowDir && ( local_disSimilar < 0.2 || local_connect ) )
        {
			if(showDebug)cout << "\t" << distance << "\n";
            // refresh distance counter
            distance = 1;
            // refresh number of blast result
            end = target;
            // the new round will begin after replace start with target
            start = target;
        }
        else
        {
            if(showDebug)cout << "\n";
            loaclUnAlignLength += abs( positionArray[target][0] - positionArray[target][1] );
        }
    }
    
    return (float)loaclUnAlignLength/(float)contigLength;
}
