//
//  GraphUtils.cpp
//  
//
//  Created by Yufeng Wu on 3/7/15.
//  Basic graph utilities e.g. BFS
//

#include "GraphUtils.h"
#include <queue>
#include <map>
#include <iostream>
#include "Utils-basic.h"
#include <cstdlib>
#include <fstream>

// ******************************************************************
// Iterator through ngbrs

class AbstractGraphNode;


AbstractGraphNodeNgbrIterator :: AbstractGraphNodeNgbrIterator( AbstractGraphNode *pNodeToVisitIn) : pNodeToVisit(pNodeToVisitIn), indexNgbr(0)
{
    //
}

void AbstractGraphNodeNgbrIterator :: First()
{
    //
    indexNgbr = 0;
}

void AbstractGraphNodeNgbrIterator :: Next()
{
    //
    ++indexNgbr;
}

bool AbstractGraphNodeNgbrIterator :: IsDone() const
{
    //
    return indexNgbr < pNodeToVisit->GetNumAdjNodes();
}

AbstractGraphNode *AbstractGraphNodeNgbrIterator ::  GetCurNgbr() const
{
    //
    vector<AbstractGraphNode *> listNgbrs;
    pNodeToVisit->GetNgbrNodes( listNgbrs );
    return listNgbrs[ indexNgbr ];
}

void AbstractGraphNode :: Dump() const
{
    //
    cout << "[" << GetName() << "]: ";
    cout << "number of neighbors: " << GetNumAdjNodes() << " : ";
    for(int i=0; i<(int)listEdges.size(); ++i)
    {
        AbstractGraphNode *pdest = listEdges[i]->GetDest();
        cout << pdest->GetName() << " ";
    }
    cout << endl;
}


// ******************************************************************
// Iterator through a list of nodes (w/ prioorty)


AbstractGraphNodeListPrioriIterator :: AbstractGraphNodeListPrioriIterator( vector<AbstractGraphNode *> &listNodes, map<AbstractGraphNode *, int> &mapNodeOrder ) : listNodesIn(listNodes), mapNodeOrderIn(mapNodeOrder), indexPri(0)
{
    //First();
    // create the list of ordered indices
    priority_queue< pair<int,int> > queuePriori;
    // queue item format: <-1*order, index of the node>; -1: queue is max priority by default
    for(int i=0; i<(int)listNodesIn.size(); ++i)
    {
        int priVal = 0;
        if( mapNodeOrder.find( listNodesIn[i] ) != mapNodeOrder.end() )
        {
            priVal = -1*mapNodeOrder[ listNodesIn[i] ];
        }
        pair<int,int> pp(priVal, i);
        queuePriori.push(pp);
    }
    // now pop it to get the ordered list
    listOrderIndices.clear();
    while( queuePriori.empty() == false )
    {
        pair<int,int> pp = queuePriori.top();
        queuePriori.pop();
        listOrderIndices.push_back( pp.second );
    }
//cout << "Ordered list of nodes: ";
//for(int i=0; i<(int)listOrderIndices.size(); ++i)
//{
//cout << listOrderIndices[i] << "  ";
//}
//cout << endl;
    
    // move to first
    First();
}

void AbstractGraphNodeListPrioriIterator :: First()
{
    indexPri = 0;
}

void AbstractGraphNodeListPrioriIterator :: Next()
{
    //
    ++indexPri;
}

bool AbstractGraphNodeListPrioriIterator :: IsDone() const
{
    //
    return indexPri >=(int)listOrderIndices.size();
}

AbstractGraphNode * AbstractGraphNodeListPrioriIterator :: GetCurNode() const
{
    //
    if( indexPri >= (int) listOrderIndices.size() )
    {
        cout << "FATAL ERROR: overlfow.\n";
        exit(1);
    }
    int indOrig = listOrderIndices[ indexPri ];
//cout << "indOrig: " << indOrig << ", size of nodes given: " << listNodesIn.size() << endl;
    return listNodesIn[ indOrig ];
}


// ******************************************************************
// abstract class for graph node

AbstractGraphNode :: ~AbstractGraphNode()
{
    //
    for(int i=0; i<(int)listEdges.size(); ++i)
    {
        delete listEdges[i];
    }
}

void AbstractGraphNode :: AddNgbrNodeBasic(AbstractGraphNode *ngbr)
{
    //
    GraphEdge *pedge = new GraphEdge(this, ngbr);
    listEdges.push_back(pedge);
}

void AbstractGraphNode :: GetNgbrNodes( vector<AbstractGraphNode *> &listNgbrs )
{
    //
    listNgbrs.clear();
    for(int i=0; i<(int)listEdges.size(); ++i)
    {
        listNgbrs.push_back( listEdges[i]->GetDest() );
    }
}

bool AbstractGraphNode :: IsNeighbor(AbstractGraphNode *pnode)
{
    for(int i=0; i<(int)listEdges.size(); ++i)
    {
        if( listEdges[i]->GetDest() == pnode )
        {
            return true;
        }
    }
    return false;
}

GraphEdge * AbstractGraphNode :: GetEdgeTo( AbstractGraphNode *pDest )
{
    for(int i=0; i<(int)listEdges.size(); ++i)
    {
        if( listEdges[i]->GetDest() == pDest )
        {
            return listEdges[i];
        }
    }

    return NULL;
}

// ******************************************************************
// Nodes w/ reference node

void GraphNodeRefExt :: AddNgbrRef( GraphNodeRefExt *pNgbr, void *pref, const string &extInfo )
{
    //
    GraphEdgeREfExt *pedge = new GraphEdgeREfExt( this, pNgbr );
    pedge->SetRef(pref);
    pedge->SetExt(extInfo);
    AddEdge( pedge );
}

void GraphNodeRefExt :: AddNgbrRef( GraphNodeRefExt *pNgbr, void *pref, const string &extInfo, double len )
{
    //
    GraphEdgeREfExt *pedge = new GraphEdgeREfExt( this, pNgbr );
    pedge->SetRef(pref);
    pedge->SetExt(extInfo);
    pedge->SetLength(len);
    AddEdge( pedge );
}

// ******************************************************************
// abstract class for graph

AbstractGraph :: ~AbstractGraph()
{
    Reset();
}

void AbstractGraph :: Reset()
{
    //
    for(int i=0; i<(int)listGraphNodes.size(); ++i )
    {
        delete listGraphNodes[i];
    }
    listGraphNodes.clear();
}

void AbstractGraph :: AddNode( AbstractGraphNode *pn )
{
    //
    listGraphNodes.push_back( pn );
}

void AbstractGraph :: InitBFS()
{
    // prepare for BFS by setting visited to be false for all nodes
    for( int i=0; i<(int)listGraphNodes.size(); ++i )
    {
        listGraphNodes[i]->SetVisited(false);
    }
}

void AbstractGraph :: BFS()
{
    // visit each node
    for(int i=0; i<(int)listGraphNodes.size(); ++i)
    {
        BFSFrom( listGraphNodes[i] );
    }
}

void AbstractGraph :: BFSFrom( AbstractGraphNode *pnodeStart )
{
    //
    InitBFS();
    
    queue< AbstractGraphNode * > queueBFS;
    queueBFS.push( pnodeStart );
    pnodeStart->SetVisited(true);
    while( queueBFS.empty() == false )
    {
        AbstractGraphNode *pnodecur = queueBFS.front();
        queueBFS.pop();
        pnodecur->PostVisit( );
        
        // now visit its ngbrs
        AbstractGraphNodeNgbrIterator itor( pnodecur );
        itor.First();
        while( itor.IsDone() == false )
        {
            //
            AbstractGraphNode *pnngbr = itor.GetCurNgbr();
            if( pnngbr->IsVisited() == false )
            {
                pnngbr->PreVisit( pnodecur );
                queueBFS.push( pnngbr );
                pnngbr->SetVisited(true);
            }
        }
        
    }
}

int AbstractGraph :: GetNumEdges() const
{
    //
    int res = 0;
    for(int i=0; i<(int)listGraphNodes.size(); ++i)
    {
        res += listGraphNodes[i]->GetNumEdges();
    }
    return res;
}

void AbstractGraph :: FindSimplePaths( set<vector<AbstractGraphNode *> > &setPaths, int maxPathLen, int maxNodeOccurInPath )
{
    //
    setPaths.clear();
    for(int i=0; i<(int)listGraphNodes.size(); ++i)
    {
        //
        vector<AbstractGraphNode *> pathInit;
        pathInit.push_back( listGraphNodes[i] );
        setPaths.insert( pathInit );
    }
    ExtendSimplePathAt(setPaths, maxPathLen, maxNodeOccurInPath);
}

void AbstractGraph :: ExtendSimplePathAt( set<vector<AbstractGraphNode *> > &setPaths,  int maxPathLen, int maxNodeOccurInPath )
{
    // loop until we are done
//GraphNodeRefExt *pnodeCurRE = dynamic_cast<GraphNodeRefExt *>(pnodeCur);
//cout << "ExtendSimplePathAt from: " << pnodeCurRE->GetExt() << endl;
    
    //const int MAX_NUM_OCCUR_NODE_IN_PATH = 100;
    
    // this is the paths that cannot be further extended
    set< vector<AbstractGraphNode *> > setPathsDone;
    //set< AbstractGraphNode *> setNodesEncountered;
    
    // keep track how many times nodes have been in the constructed paths
    map<AbstractGraphNode *, int> mapCountsNodesOccurInPath;
    // populate the node counts
    for( set<vector<AbstractGraphNode *> > :: iterator it = setPaths.begin(); it != setPaths.end(); ++it  )
    {
        for(int i=0; i<(int)it->size(); ++i)
        {
            //
            AbstractGraphNode *pn = (*it)[i];
            if( mapCountsNodesOccurInPath.find(pn) == mapCountsNodesOccurInPath.end() )
            {
                mapCountsNodesOccurInPath.insert( map<AbstractGraphNode *, int> :: value_type(pn, 0)  );
            }
            ++mapCountsNodesOccurInPath[pn];
        }
    }
    
    while(setPaths.size() > 0)
    {
//cout << "The number of paths being considered: " << setPaths.size() << endl;
//cout << "The number of paths finished: " << setPathsDone.size() << endl;
        
        //bool fCont = false;
    
        // extend all paths from this graph from current node
        // if this node is not added, then just add the node itself as start
        set<vector<AbstractGraphNode *> > setPathsNew;
        
        //bool fAdded = false;
        for( set<vector<AbstractGraphNode *> > :: iterator it = setPaths.begin(); it != setPaths.end(); ++it  )
        {
            bool fPathExt = false;
            
            if( (int)it->size() < maxPathLen )
            {
            
                //
                AbstractGraphNode *pnodePathEnd = (*it)[ it->size()-1 ];
                
                vector<AbstractGraphNode *> listNgbrs;
                pnodePathEnd->GetNgbrNodes( listNgbrs );
                
                
                for(int jj=0; jj<(int)listNgbrs.size(); ++jj)
                {
                    AbstractGraphNode *pnodeStep = listNgbrs[jj];
                    
                    // don't allow to extend to this node if it has been used too often
                    if( mapCountsNodesOccurInPath.find(pnodeStep) == mapCountsNodesOccurInPath.end() )
                    {
                        cout <<  "Fail to find\n";
                        exit(1);
                    }
                    if( mapCountsNodesOccurInPath[pnodeStep] > maxNodeOccurInPath )
                    {
                        continue;
                    }
                    
                    //if( pnodeStep == pnodeCur)
                    //{
                    //    fAdded = true;
                    //}
                
                    // if this node is not yet in the path and there is an edge from the end of the path to this node, add it
                    bool fExist = false;
                    for(int i=0; i<(int)it->size(); ++i)
                    {
                        if( pnodeStep == (*it)[i] )
                        {
                            fExist = true;
                            break;
                        }
                    }
                    if( fExist == false  )
                    {
                        // add to the end
                        vector<AbstractGraphNode *> pathNew=*it;
                        pathNew.push_back( pnodeStep );
    //GraphNodeRefExt *pnodeStepRE = dynamic_cast<GraphNodeRefExt *>(pnodeStep);
    //cout << "  extend path to: " << pnodeStepRE->GetExt() << endl;
                        setPathsNew.insert(pathNew);
                        //fCont = true;
                        fPathExt = true;
                        
                        // add counts by one
                        ++mapCountsNodesOccurInPath[pnodeStep];
                    }
                }
            }
            
            if( fPathExt == false )
            {
                // add this original path since this is done
                setPathsDone.insert( *it );
            }
        }
        
        // if not added yet, add it
        //if( fAdded == false && fInit == false )
        //{
        //    vector<AbstractGraphNode *> path;
        //    //path.push_back( pnodeCur );
        //    setPathsNew.insert( path );
        //    fCont = true;
        //    fInit = true;
        //}
        setPaths = setPathsNew;
        
        // trim duplicate path
        //RemoveSubPathsFrom(setPaths);
        
        
        // stop if nothing new is done
        //if(  )
        //{
        //    break;
        //}
    }
    
    // add back those paths done
    for( set<vector<AbstractGraphNode *> > :: iterator it = setPathsDone.begin(); it != setPathsDone.end(); ++it)
    {
        setPaths.insert(*it);
    }
}


void AbstractGraph:: FindSimplePathsBoundedLength( set<vector<AbstractGraphNode *> > &setPaths, int maxPathLen )
{
    // find all simple paths up to the given length
    // approach: incrementally grow the graph by adding edges between two nodes (each annotated by the set of original
    // list of nodes; so that we can get the path)
    
    setPaths.clear();
    
    // format: <src, dest>, <list of nodes along the path from src to dest>
    map< pair<AbstractGraphNode *, AbstractGraphNode *>, vector<AbstractGraphNode *> > mapEdgeFromToWithPath;
    map<AbstractGraphNode *, set<AbstractGraphNode *> > mapSetNgbrsExpanded;
    map<AbstractGraphNode *, set<AbstractGraphNode *> > mapSetNgbrsToProcess;
    map<AbstractGraphNode *, set<AbstractGraphNode *> > mapSetNgbrsToOutput;    // these are not belonging to subpath of some other
    map<AbstractGraphNode *, set<AbstractGraphNode *> > mapSetNgbrsExcluded;    // should not be output
    
    // init: with the original graph edges
    for(int i=0; i<(int)listGraphNodes.size(); ++i)
    {
        //
        AbstractGraphNode *psrc = listGraphNodes[i];
        vector<AbstractGraphNode *> listNgbrs;
        psrc->GetNgbrNodes( listNgbrs );
        
        for(int i=0; i<(int)listNgbrs.size(); ++i)
        {
            AbstractGraphNode *pDest = listNgbrs[i];
        
            vector<AbstractGraphNode *> pathInit;
            pathInit.push_back( psrc );
            pathInit.push_back( pDest );
            
            pair<AbstractGraphNode *, AbstractGraphNode *> pp( psrc, pDest );
            mapEdgeFromToWithPath.insert( map< pair<AbstractGraphNode *, AbstractGraphNode *>, vector<AbstractGraphNode *> > :: value_type( pp, pathInit ) );
        }
        
        // also update its ngbrs
        set<AbstractGraphNode *> setNgbrsInit;
        for(int i=0; i<(int)listNgbrs.size(); ++i)
        {
            setNgbrsInit.insert( listNgbrs[i] );
        }
        mapSetNgbrsExpanded.insert( map<AbstractGraphNode *, set<AbstractGraphNode *> > :: value_type( psrc, setNgbrsInit ) );
        mapSetNgbrsToOutput.insert(map<AbstractGraphNode *, set<AbstractGraphNode *> > :: value_type( psrc, setNgbrsInit ));
        set<AbstractGraphNode *> setEmpty;
        mapSetNgbrsExcluded.insert( map<AbstractGraphNode *, set<AbstractGraphNode *> > :: value_type( psrc, setEmpty ) );
    }
    
    // upon start, init the ngbrs to process
    mapSetNgbrsToProcess = mapSetNgbrsExpanded;
    
    
    // loop up to that many times
    int loopIndex = 1;
    while( loopIndex < maxPathLen )
    {
        // consider each node to see if
        bool fNewEdgeAdded = false;
        for( int i=0; i<(int)listGraphNodes.size(); ++i )
        {
            //
            AbstractGraphNode *pSrc = listGraphNodes[i];
            
            set<AbstractGraphNode *> setNgbrNext;
            
            set<AbstractGraphNode *> setNgbrsNow = mapSetNgbrsToProcess[ pSrc ];
            for( set<AbstractGraphNode *> :: iterator it = setNgbrsNow.begin(); it != setNgbrsNow.end(); ++it )
            {
                //
                AbstractGraphNode *pDest = *it;
                
                // get the current path from pSrc to pDest
                pair<AbstractGraphNode *, AbstractGraphNode *> pp( pSrc, pDest );
                if( mapEdgeFromToWithPath.find(pp) == mapEdgeFromToWithPath.end() )
                {
                    cout << "Fatal error: quite\n";
                    exit(1);
                }
                vector<AbstractGraphNode *> currPathStep = mapEdgeFromToWithPath[pp];
                
                // get its ngbr of this pDest
                vector<AbstractGraphNode *> listNgbrsStep;
                pDest->GetNgbrNodes( listNgbrsStep );
                
                bool fExtended = false;
                
                for(int jj=0; jj<(int)listNgbrsStep.size(); ++jj)
                {
                    AbstractGraphNode *pnodeStep = listNgbrsStep[jj];
                    
                    // don't allow to extend to this node has already been (indirect) ngbr of pSrc
                    if( pnodeStep == pSrc ||  mapSetNgbrsExpanded[pSrc].find(pnodeStep) != mapSetNgbrsExpanded[pSrc].end() )
                    {
                        continue;
                    }
                    
                    // add to ngbrhood
                    mapSetNgbrsExpanded[pSrc].insert( pnodeStep );
                    
                    // also remember the path
                    pair<AbstractGraphNode *, AbstractGraphNode *> pp2( pSrc, pnodeStep );
                    vector<AbstractGraphNode *> pathNew =  currPathStep;
                    pathNew.push_back( pnodeStep );
                    if( mapEdgeFromToWithPath.find(pp2) == mapEdgeFromToWithPath.end() )
                    {
                        //cout << "FATAL ERROR: the path should not exist.\n";
                        //exit(1);
                        //}
                        mapEdgeFromToWithPath.insert( map< pair<AbstractGraphNode *, AbstractGraphNode *>, vector<AbstractGraphNode *> > :: value_type( pp2, pathNew ) );
                    
                        //GraphNodeRefExt *pnodeStepRE = dynamic_cast<GraphNodeRefExt *>(pnodeStep);
                        //cout << "  extend path to: " << pnodeStepRE->GetExt() << endl;
                        //setPathsNew.insert(pathNew);
                        //fCont = true;
                        setNgbrNext.insert(pnodeStep);
                        
                        if( mapSetNgbrsExcluded[pSrc].find( pnodeStep) == mapSetNgbrsExcluded[pSrc].end() )
                        {
                            mapSetNgbrsToOutput[pSrc].insert( pnodeStep );
                        }
                        // exclude the other sub-paths
                        for(int i=1; i<(int)currPathStep.size(); ++i )
                        {
                            mapSetNgbrsExcluded[currPathStep[i] ].insert( pnodeStep );;
                        }
                        
                        fNewEdgeAdded = true;
                        fExtended = true;
                    }
                }
                
                // if this is expanded, remove it from output list
                if( fExtended == true )
                {
                    mapSetNgbrsToOutput[pSrc].erase( pDest );
                    mapSetNgbrsExcluded[pSrc].insert(pDest);
                }

            }
//cout << "i: " << i << "  size of ngbrs to output: " << mapSetNgbrsToOutput[pSrc].size() << endl;
//vector<AbstractGraphNode *> listNgbrs;
//pSrc->GetNgbrNodes( listNgbrs );
//cout << "List of ngbrs: " << listNgbrs.size() << endl;
            
            mapSetNgbrsToProcess[pSrc] = setNgbrNext;
        }
        
        if( fNewEdgeAdded == false )
        {
            // that is,
            break;
        }
        
        //
        ++loopIndex;
    }
    
    // collect all the paths
    for( map< pair<AbstractGraphNode *, AbstractGraphNode *>, vector<AbstractGraphNode *> > :: iterator it = mapEdgeFromToWithPath.begin(); it != mapEdgeFromToWithPath.end(); ++it)
    {
        // only output a path if it is not part of another (longer) path
        if( mapSetNgbrsToOutput[ it->first.first ].find( it->first.second ) != mapSetNgbrsToOutput[it->first.first].end() &&
           mapSetNgbrsExcluded[it->first.first].find(it->first.second) == mapSetNgbrsExcluded[it->first.first].end()  )
        {
            setPaths.insert(it->second);
        }
    }
    
    // remove duplicate paths
//cout << "Before removing duplicate paths: " << setPaths.size() << endl;
//    RemoveSubPathsFromFast(setPaths);
//cout << "AFTER removing duplicate paths: " << setPaths.size() << endl;
}

void AbstractGraph :: FindSimplePathsTopSort( set<vector<AbstractGraphNode *> > &setPaths, int maxNodeOccurInPath )
{
    // find a set of simple paths (of arbitary length)
    // approach: perform a topological sort using the SCC function
    // which sorts the graph if it is acyclic; if it is not, we may have a component w/ more than 1 node approach: 
    // this should still work since we are using DFS
    vector< set<AbstractGraphNode *> > setSCCs;
    SCC( setSCCs );
//cout << "The number of strongly connected components: " << setSCCs.size() << endl;
//#if 0
//for( int i=0; i<(int)setSCCs.size(); ++i )
//{
//cout << setSCCs[i].size() << "  ";
//}
//cout << endl;
    
     //
    //set<AbstractGraphNode *> setNodesExclude;
    //map<AbstractGraphNode *, int> mapNodesOccurCounts;       // should not find path from these nodes since they are already done
    map<AbstractGraphNode *, int> mapNodesOrderIndex;          // at which order, this node occurs in the topologlical sort
    
    // use the sorted order
    vector<AbstractGraphNode *> listNodesToExp;
    int indexOrder = 1;
    for( int i=0; i<(int)setSCCs.size(); ++i )
    {
        for( set<AbstractGraphNode *> :: iterator it = setSCCs[i].begin(); it != setSCCs[i].end(); ++it)
        {
            //(*it)->SetVisited(false);
            listNodesToExp.push_back( *it );
            
            //mapNodesOccurCounts.insert( map<AbstractGraphNode *,int> :: value_type(*it, 0) );
            
            mapNodesOrderIndex.insert( map<AbstractGraphNode *,int> :: value_type(*it, indexOrder++)  );
        }
    }
    
    // find
    vector<AbstractGraphNode *> listPathRoots;
    FindSimplePathsTopSortStart( true, setSCCs, listPathRoots );
//cout << "Total number of candidate roots: " << listPathRoots.size() << endl;
//for(int i=0; i<(int)listPathRoots.size();++i)
//{
//cout << "  " << listPathRoots[i]->GetName() << endl;
//}
    vector<AbstractGraphNode *> listPathEnds;
    FindSimplePathsTopSortStart( false, setSCCs, listPathEnds );
//cout << "Total number of candidate path ends: " << listPathEnds.size() << endl;
//for(int i=0; i<(int)listPathEnds.size();++i)
//{
//cout << "  " << listPathEnds[i]->GetName() << endl;
//}
    
// find all the components
//set< set<AbstractGraphNode *> > setComponents;
//FindComponents( listNodesToExp, setComponents );
//cout << "+++ Number of weakly connected component: " << setComponents.size() << endl ;
//for( set< set<AbstractGraphNode *> > :: iterator it = setComponents.begin(); it != setComponents.end(); ++it  )
//{
//cout << "Find one component: \n";
//for( set<AbstractGraphNode *> :: iterator it2 = it->begin(); it2 != it->end(); ++it2)
//{
//cout << "  " << (*it2)->GetName() << endl;
//}
//}
    
    

    
    // iterate over
    for(int i=0; i<(int)listPathRoots.size(); ++i)
    {
//cout << "Iterate: " << i << endl;
        AbstractGraphNode *pnodeCurr = listPathRoots[i];
//cout << "Processing candidate root: ";
//pnodeCurr->Dump();
        //if( setNodesExclude.find( pnodeCurr) == setNodesExclude.end() )
        //{
        // set all nodes to be unvisited so we can do it once again
        //for(int j=0; j<(int)listNodesToExp.size(); ++j)
        //{
        //    listNodesToExp[j]->SetVisited(false);
        //}
        
        // using this node as path root to find path
        vector<AbstractGraphNode *> pathSoFar;
        pathSoFar.push_back(pnodeCurr);
        
        set<vector<AbstractGraphNode *> > setPathsCurr;
        FindSimplePathsTopSortFrom( pnodeCurr, listNodesToExp, listPathEnds, setPathsCurr );
        
        
        // take the longest k path of it
//cout << "Path found at this step....\n";
        map<int, set< const vector<AbstractGraphNode *> *> > mapLenToPathPtrs;
        for( set<vector<AbstractGraphNode *> > :: iterator itg = setPathsCurr.begin(); itg != setPathsCurr.end(); ++itg )
        {
//for(int i=0; i<(int)itg->size(); ++i)
//{
//cout << "  " << (*itg)[i]->GetName();
//}
//cout << endl;
            
            int neglen = -1.0* (*itg).size();
            if( mapLenToPathPtrs.find( neglen) == mapLenToPathPtrs.end() )
            {
                set< const vector<AbstractGraphNode *>* > ss;
                mapLenToPathPtrs.insert(  map<int, set< const vector<AbstractGraphNode *> *> > :: value_type(neglen,  ss) );
            }
            mapLenToPathPtrs[neglen].insert( &(*itg) );
        }
        int numOut = 0;
        for( map<int, set< const vector<AbstractGraphNode *> *> > :: iterator it1 = mapLenToPathPtrs.begin(); it1 != mapLenToPathPtrs.end() ; ++it1 )
        {
            for( set< const vector<AbstractGraphNode *> *> :: iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2  )
            {
                setPaths.insert( *(*it2) );
                ++numOut;
                
                if( numOut > maxNodeOccurInPath)
                {
                    break;
                }
            }
            if( numOut > maxNodeOccurInPath)
            {
                break;
            }
        }
        
        //}
    }
   
#if 0
cout << "^^^^^^ The list of found paths: \n";
int pathIndex = 1;
for( set<vector<AbstractGraphNode *> > :: iterator it = setPaths.begin(); it != setPaths.end(); ++it )
{
cout << "PATH" << pathIndex++ << ": ";
for(int i=0; i<(int)it->size(); ++i )
{
cout << "  " << (*it)[i]->GetName();
}
cout << endl;
}
#endif
}

void AbstractGraph :: FindSimplePathsTopSortFrom(AbstractGraphNode *pnodePathRoot, const vector<AbstractGraphNode *> &listNodesToExp, const vector<AbstractGraphNode *> &listPathEnds, set<vector<AbstractGraphNode *> > &setPaths)
{
    // using DP to find, from node pnodePathRoot, set of paths from the spcpfied root to the particular node
    // when there are multiple paths, use the shortest one (if you want longest one, use negitive)
    // here, assume ordered, so ignore edges going the other way
    // create a map to keep track node position in listNodesToExp, which is assumed to at least partially fit topological order
    const double MAX_PATH_LEN = 1.0e100;
    map<AbstractGraphNode *, int> mapNodeRanks;
    for(int i=0; i<(int)listNodesToExp.size(); ++i)
    {
        mapNodeRanks.insert( map<AbstractGraphNode *, int> :: value_type(listNodesToExp[i], i) );
    }
    
    //
    vector< pair< double, vector<AbstractGraphNode *> > > DPTblPathLenPath;
    // init
    for(int i=0;i<(int)listNodesToExp.size(); ++i)
    {
        vector<AbstractGraphNode *> elist;
        pair<double, vector<AbstractGraphNode *> > pp( MAX_PATH_LEN, elist );
        DPTblPathLenPath.push_back(pp);
    }
    int posRoot = mapNodeRanks[pnodePathRoot];
    DPTblPathLenPath[posRoot].first = 0.0;
    DPTblPathLenPath[posRoot].second.push_back( pnodePathRoot );
    
    // now for each node after the root, check their neighbors; update if a better path is found; ignore back-edges
    for(int i=posRoot; i<(int)listNodesToExp.size(); ++i)
    {
//cout << "******i = " << i  << ", IN DP: proc: " << listNodesToExp[i]->GetName() << endl;
        if( DPTblPathLenPath[i].first >= MAX_PATH_LEN )
        {
            // skip if there is no path to this node            
//cout << "No path info. \n";
            continue;
        }
        AbstractGraphNode *pnodeCurr = listNodesToExp[i];
        double pathCur = DPTblPathLenPath[i].first;
//cout << "pathcur: " << pathCur << endl;
        vector<AbstractGraphNode *> listNgbrs;
        pnodeCurr->GetNgbrNodes( listNgbrs );
        for(int kk=0; kk<(int)listNgbrs.size(); ++kk )
        {
            AbstractGraphNode *pnodeNgbr = listNgbrs[kk];
            int posNgbr = mapNodeRanks[pnodeNgbr];
//cout << "posNgbr: " << posNgbr  << endl;
            
            // skip if order is not right
            if( posNgbr < i)
            {
//cout << "Out of order\n";
                continue;
            }
            double pathLen = pnodeCurr->GetEdgeTo(pnodeNgbr)->GetLength();
//cout << "pathLen: " << pathLen << endl;
            
            // is this a better path?
            if( pathCur + pathLen < DPTblPathLenPath[posNgbr].first )
            {
                // update path length
                DPTblPathLenPath[ posNgbr].first = pathCur + pathLen;
                // update path itself
                DPTblPathLenPath[ posNgbr ].second = DPTblPathLenPath[i].second;
                DPTblPathLenPath[ posNgbr ].second.push_back( pnodeNgbr );
//cout << pnodeNgbr->GetName() << "  BETTER: pathlen: " << DPTblPathLenPath[ posNgbr].first << ", list of nodes: ";
//for(int jj=0; jj<(int) DPTblPathLenPath[ posNgbr ].second.size(); ++jj )
//{
//cout << "  " << DPTblPathLenPath[ posNgbr ].second[jj]->GetName();
//}
//cout << endl;
            }
        }
    }
    
    // find paths by checking those candidate roots
    setPaths.clear();
    for(int i=0; i<(int)listPathEnds.size(); ++i )
    {
        //
        int posEnd = mapNodeRanks[ listPathEnds[i] ];
        if( DPTblPathLenPath[posEnd].first < MAX_PATH_LEN )
        {
            // output it
            setPaths.insert( DPTblPathLenPath[posEnd].second );
        }
    }
}

#if 0
void AbstractGraph :: FindSimplePathsTopSortFrom(AbstractGraphNode *pnodeCurr, vector<AbstractGraphNode *> &pathSoFar, map<AbstractGraphNode *, int> &mapNodesOrderIndex, set<vector<AbstractGraphNode *> > &setPaths)
{
    // pathSofar: what nodes are on the path to this node; 
    // mapNodesOrderIndex: for priority of a node
//cout << "--current path length: " << pathSoFar.size() << ", number of found path: " << setPaths.size() << endl;
    // find paths from this current node
    pnodeCurr->SetVisited(true);
    //++mapNodesOccurCounts[pnodeCurr];
    //setNodesExclude.insert( pnodeCurr );
    
    vector<AbstractGraphNode *> listNgbrs;
    pnodeCurr->GetNgbrNodes( listNgbrs );
//cout << "Number of ngbr nodes: " << listNgbrs.size() << endl;
    bool fExt = false;
    //while (itor.IsDone() == false)
    AbstractGraphNodeListPrioriIterator itor( listNgbrs, mapNodesOrderIndex );
    itor.First();
    //for(int i=0; i<(int)listNgbrs.size(); ++i)
    while( itor.IsDone() == false )
    {
        AbstractGraphNode *pnodeStep = itor.GetCurNode();
        //if( mapNodesOccurCounts.find( pnodeStep ) == mapNodesOccurCounts.end() )
        //{
        //    cout << "FATAL ERROR: node not found.\n";
        //    exit(1);
        //}
        //if( pnodeStep->IsVisited() == false && mapNodesOccurCounts[ pnodeStep ] < maxNodeOccurInPath )
        if( pnodeStep->IsVisited() == false  )
        {
            //
            vector<AbstractGraphNode *> pathNext = pathSoFar;
            pathNext.push_back( pnodeStep );
            FindSimplePathsTopSortFrom(pnodeStep, pathNext, mapNodesOrderIndex, setPaths);
            fExt = true;
        }
        
        itor.Next();
    }
    if( fExt == false )
    {
        // output a path if a path exists
        if( pathSoFar.size() > 1 )
        {
            //
            setPaths.insert( pathSoFar );
            
            // update counts for all the nodes along the path
            //for(int jj=0; jj<(int)pathSoFar.size(); ++jj)
            //{
            //    ++mapNodesOccurCounts[ pathSoFar[jj] ];
            //}
            
        }
    }
}
#endif

#if 0
void AbstractGraph :: RemoveSubPathsFrom( set<vector<AbstractGraphNode *> > &setPaths )
{
    // remove a path if it is a subpath of another (longer) path
    // approach: a little heuristic: go from longer path to shorter ones
    // if all nodes of a path belong to a longer path (in any order), don't include it
    set<vector<AbstractGraphNode *> > setPathsKeep;
    
    map<int, map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > > mapPathNodesSetLen;
    
    for( set<vector<AbstractGraphNode *> > :: iterator it = setPaths.begin(); it != setPaths.end(); ++it)
    {
        int szPath = -1*it->size();
        if( mapPathNodesSetLen.find(szPath) == mapPathNodesSetLen.end() )
        {
            map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > mapEmpty;
            mapPathNodesSetLen.insert( map<int, map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > > :: value_type(szPath, mapEmpty) );
        }
        set<AbstractGraphNode *> setNodes;
        for(int i=0; i<(int)it->size(); ++i)
        {
            setNodes.insert( (*it)[i] );
        }
        
        mapPathNodesSetLen[szPath].insert( map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > :: value_type( &(*it), setNodes )  );
    }
    for( map<int, map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > > :: iterator it = mapPathNodesSetLen.begin(); it != mapPathNodesSetLen.end(); ++it  )
    {
        for( map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > :: iterator it0 = it->second.begin(); it0 != it->second.end(); ++it0 )
        {
            bool fAdd = true;
            for( map<int, map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > > :: iterator it2 = mapPathNodesSetLen.begin(); it2 != it; ++it2  )
            {
                //
                for( map<const vector<AbstractGraphNode *> *, set<AbstractGraphNode *> > :: iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3 )
                {
                    // is the set contained in this one?
                    if( IsSetContainerGen1(it3->second, it0->second ) == true )
                    {
                        fAdd = false;
                        break;
                    }
                }
                if( fAdd == false )
                {
                    break;
                }
            }
            if( fAdd == true )
            {
                setPathsKeep.insert( *(it0->first) );
            }
        }
    }
    setPaths = setPathsKeep;
}

void AbstractGraph :: RemoveSubPathsFromFast( set<vector<AbstractGraphNode *> > &setPaths )
{
    // remove a path if it is a subpath of another (longer) path
    // approach: a little heuristic: go from longer path to shorter ones
    // if all nodes of a path belong to a longer path (in any order), don't include it
    set<vector<AbstractGraphNode *> > setPathsKeep;
    
    map<int, set<const vector<AbstractGraphNode *> *  > > mapPathNodesSetLen;
    
    for( set<vector<AbstractGraphNode *> > :: iterator it = setPaths.begin(); it != setPaths.end(); ++it)
    {
        int szPath = -1*it->size();
        if( mapPathNodesSetLen.find(szPath) == mapPathNodesSetLen.end() )
        {
            set<const vector<AbstractGraphNode *> * > mapEmpty;
            mapPathNodesSetLen.insert( map<int, set<const vector<AbstractGraphNode *>* > > :: value_type(szPath, mapEmpty) );
        }
        
        mapPathNodesSetLen[szPath].insert( &(*it) );
    }
    // keep track all pairs of the nodes along the same path that have been reached
    set< pair<  AbstractGraphNode*,  AbstractGraphNode * > > setPathNodesPairsDone;
    set<AbstractGraphNode *> setPathNodesDone;
    for( map<int, set<const vector<AbstractGraphNode *>* > > :: iterator it = mapPathNodesSetLen.begin(); it != mapPathNodesSetLen.end(); ++it  )
    {
        for( set<const vector<AbstractGraphNode *>* > :: iterator it0 = it->second.begin(); it0 != it->second.end(); ++it0 )
        {
            // check only the begining and end
            //pair< AbstractGraphNode *,  AbstractGraphNode *> pp( (*(*it0))[0], (*(*it0))[ (*it0)->size()-1 ] );
            //if( setPathNodesPairsDone.find(pp ) == setPathNodesPairsDone.end() )
            if( setPathNodesDone.find( (*(*it0))[0]) == setPathNodesDone.end()  && setPathNodesDone.find( (*(*it0))[ (*it0)->size()-1 ]  ) == setPathNodesDone.end() )
            {
                // this is not done yet,
                setPathsKeep.insert( *(*it0) );
                
                // also keep track of all pair of paths
                for(int i=0; i<(int)(*it0)->size(); ++i)
                {
                    setPathNodesDone.insert( (*(*it0))[i] );
                    //for(int j=i+1; j<(int)(*it0)->size(); ++j)
                    //{
                    //    pair< AbstractGraphNode *,  AbstractGraphNode *> pp0( (*(*it0))[i], (*(*it0))[j]  );
                    //    setPathNodesPairsDone.insert(pp0);
                    //}
                }
            }
        }
    }
    setPaths = setPathsKeep;
}
#endif

void AbstractGraph :: SCC( vector< set<AbstractGraphNode *> > &setSCCs )
{
    // find all strongly connected components (SCC)
    // note: must needs reference/extension to make it work
    int indexNext = 1;
    stack<GraphNodeSCC *> stackNodesToExp;
    
    vector<GraphNodeSCC *> listGraphNodesSCC;
    map< AbstractGraphNode *, GraphNodeSCC *> mapOrigNodeToSCCNode;
    
    // create (ammended) list of nodes
    for(int i=0; i<(int)listGraphNodes.size(); ++i)
    {
        //
        GraphNodeSCC *psrc = new GraphNodeSCC( listGraphNodes[i] );
        listGraphNodesSCC.push_back(psrc);
        mapOrigNodeToSCCNode.insert( map< AbstractGraphNode *, GraphNodeSCC *> :: value_type( listGraphNodes[i], psrc ) );
    }

    // start finding SCC 
    vector< set<AbstractGraphNode *> > setSCCsRev;
    for(int i=0; i<(int)listGraphNodesSCC.size(); ++i)
    {
//cout << "Processing SCC node: " << i << endl;
        //
        GraphNodeSCC *psrc = listGraphNodesSCC[i];
        if(psrc->GetIndex() < 0 )
        {
            SCCFrom( psrc, indexNext, stackNodesToExp, mapOrigNodeToSCCNode, setSCCsRev );
        }
    }
    
    // reverse the order so that we can get topological order
    setSCCs.clear();
    for(int i=setSCCsRev.size()-1; i>=0; --i)
    {
        setSCCs.push_back( setSCCsRev[i] );
    }
    
    //
    for(int i=0; i<(int)listGraphNodesSCC.size(); ++i)
    {
        delete listGraphNodesSCC[i];
    }
    listGraphNodesSCC.clear();
}

static int GetMinInt(int v1, int v2)
{
    if(v1 <= v2)
    {
        return v1;
    }
    else
    {
        return v2;
    }
}

void AbstractGraph :: SCCFrom( GraphNodeSCC *pnodeCurr, int &indexNext, stack<GraphNodeSCC *> &stackNodesToExp, map< AbstractGraphNode *, GraphNodeSCC *> & mapOrigNodeToSCCNode, vector< set<AbstractGraphNode *> > &setSCCs )
{
    //
    if( pnodeCurr == NULL )
    {
        cout << "pnodeCurr: null\n";
        exit(1);
    }
//cout << "indexNext: " << indexNext << endl;
//cout << "Before SCCFrom: currnode: ";
//pnodeCurr->Dump();
    pnodeCurr->SetIndex( indexNext );
    pnodeCurr->SetLowLink(indexNext);
    ++indexNext;
    stackNodesToExp.push( pnodeCurr );
    pnodeCurr->SetOnStack( true );
    
    // Consider successors of v
    //AbstractGraphNodeNgbrIterator itor( pnodeCurr );
    //itor.First();
    
    vector<AbstractGraphNode *> listNgbrs;
    pnodeCurr->GetNgbrNodes( listNgbrs );
    //while (itor.IsDone() == false)
    for(int i=0; i<(int)listNgbrs.size(); ++i)
    {
        AbstractGraphNode *pNodeNgbr = listNgbrs[i];
        //AbstractGraphNode *pNodeNgbr = itor.GetCurNgbr();
        if( mapOrigNodeToSCCNode.find( pNodeNgbr) == mapOrigNodeToSCCNode.end() )
        {
            cout << "Map: not found\n";
            exit(1);
        }
        GraphNodeSCC *pNodeNgbrSCC = mapOrigNodeToSCCNode[ pNodeNgbr ];
        if( pNodeNgbrSCC == NULL )
        {
            cout << "Node: cannot be NULL\n";
            exit(1);
        }
        
    
        if (pNodeNgbrSCC->GetIndex() < 0 )
        {
            // Successor w has not yet been visited; recurse on it
            SCCFrom( pNodeNgbrSCC, indexNext, stackNodesToExp, mapOrigNodeToSCCNode, setSCCs );
            pnodeCurr->SetLowLink(  GetMinInt(pnodeCurr->GetLowLink(), pNodeNgbrSCC->GetLowLink() ) );
        }
        else if ( pNodeNgbrSCC->IsOnStack()  == true )
        {
            // Successor w is in stack S and hence in the current SCC
            pnodeCurr->SetLowLink(  GetMinInt( pnodeCurr->GetLowLink() , pNodeNgbrSCC->GetIndex() ) );
        }
        
        //itor.Next();
    }
//cout << "After checking ngbrs: currnode: ";
//pnodeCurr->Dump();
    // If v is a root node, pop the stack and generate an SCC
    if( pnodeCurr->GetLowLink() == pnodeCurr->GetIndex() )
    {
        set<AbstractGraphNode *> sccNew;
        while( true  )
        {
            if( stackNodesToExp.empty() == true)
            {
                cout << "FATAL ERROR: stack\n";
                exit(1);
            }
            GraphNodeSCC *pnodew = stackNodesToExp.top();
            stackNodesToExp.pop();
            if( pnodew == NULL )
            {
                cout << "pnodew: cannot be NULL\n";
                exit(1);
            }
            pnodew->SetOnStack(false);
            sccNew.insert( pnodew->GetNodeOrig() );
//cout << "Find a node in the new SCC: ";
//pnodew->Dump();
            if( pnodew == pnodeCurr )
            {
                break;
            }
        }
        // add this SCC
        setSCCs.push_back( sccNew );
//cout << "Find one new SCC: size: " << sccNew.size() << endl;
    }
//cout << "AFTER SCCFrom: currnode: ";
//pnodeCurr->Dump();

}

static void OutputQuotedString(ofstream &outFile, const char *buf)
{
	outFile << '"';
	outFile << buf;
	outFile << '"';
}

void AbstractGraph :: OutputGML(const char *filename)
{
    // output graph in GML format
    // Now output a file in GML format
	// First create a new name
	string name = filename;
    //cout << "num edges = " << listEdges.size() << endl;
    
	//DEBUG("FileName=");
	//DEBUG(name);
	//DEBUG("\n");
	// Now open file to write out
	ofstream outFile( name.c_str() );
    
	// First output some header info
	outFile << "graph [\n";
	outFile << "comment ";
	OutputQuotedString(outFile, "Automatically generated by Graphing tool");
	outFile << "\ndirected  1\n";
	outFile << "id  1\n";
	outFile << "label ";
	OutputQuotedString ( outFile, "To be more meaningful later....\n");
    
	// Now output all the vertices
	int i;
	int numVerts = (int)listGraphNodes.size();
	for( i= 0;  i< numVerts; ++i)
	{
		outFile << "node [\n";
		//char name[100];
		//name[0] = 'v';
		//sprintf(&name[1], "%d", i+1);
		outFile << "id " <<  i+1  << endl;
		outFile << "label ";
		OutputQuotedString (outFile,  listGraphNodes[i]->GetName().c_str()  );
		outFile << endl;
		outFile << "defaultAtrribute   1\n";
		outFile << "]\n";
	}
    
	// Now output all the edges
	for(i=0; i< numVerts; ++i )
	{
		for(int j=0; j<numVerts; ++j)
		{
            if( i == j )
            {
                continue;
            }
            
			if( listGraphNodes[i]->IsNeighbor (  listGraphNodes[j] ) == true )
			{
                
                //cout << "Output an edge \n";
				outFile << "edge [\n";
				outFile << "source " << i+1 << endl;
				outFile << "target  " << j+1 << endl;
				outFile << "label " ;
				OutputQuotedString( outFile,  ""  );
				outFile << "\n";
				outFile << "]\n";
			}
		}
	}
    
    
	// Finally quite after closing file
	outFile << "\n]\n";
	outFile.close();
}

void AbstractGraph :: FindSimplePathsTopSortStart(bool fStart, const vector< set<AbstractGraphNode *> > &setSCCs, vector<AbstractGraphNode *> &listPathRoots)
{
    // fStart: = true if finding the begining, =false: find the ending candidate
    // find list of possible maximal path root: either they don't have incoming edges or all incoming edges come from the same SCC
    map< AbstractGraphNode *, const set<AbstractGraphNode *> * > mapNodeToSCCPtr;
    set<AbstractGraphNode *> setCandidateRoots;
    vector<AbstractGraphNode *> vecNodesAll;
    for(int i=0;i<(int)setSCCs.size(); ++i)
    {
        for( set<AbstractGraphNode *> :: iterator it = setSCCs[i].begin(); it != setSCCs[i].end(); ++it)
        {
            setCandidateRoots.insert( *it );
            vecNodesAll.push_back( *it );
            
            mapNodeToSCCPtr.insert( map< AbstractGraphNode *, const set<AbstractGraphNode *> * >  :: value_type( *it, &setSCCs[i] ) );
        }
    }
    // check all edges of each node
    for(int i=0; i<(int)vecNodesAll.size(); ++i)
    {
        vector<AbstractGraphNode *> listNgbrs;
        vecNodesAll[i]->GetNgbrNodes( listNgbrs );
        for(int jj=0; jj<(int)listNgbrs.size(); ++jj)
        {
            // remove the sink if it comes from a different SCC
            if( mapNodeToSCCPtr[ listNgbrs[jj] ] != mapNodeToSCCPtr[ vecNodesAll[i] ] )
            {
                if( fStart == true )
                {
                    setCandidateRoots.erase( listNgbrs[jj] );
                }
                else
                {
                    // otherwise it is the ending
                    setCandidateRoots.erase( vecNodesAll[i] );
                    break;
                }
            }
        }
    }
    // for each SCC w/ size > 1, if not all the nodes are still in the list, then remove all; otherwise, keep only the first node
    for(int i=0;i<(int)setSCCs.size(); ++i)
    {
        if( setSCCs[i].size() > 1 )
        {
            //
            bool fAllIn = true;
            for( set<AbstractGraphNode *> :: iterator it = setSCCs[i].begin(); it != setSCCs[i].end(); ++it )
            {
                if( setCandidateRoots.find(*it) == setCandidateRoots.end()  )
                {
                    fAllIn = false;
                    break;
                }
            }
            // remove all remaining
            for( set<AbstractGraphNode *> :: iterator it = setSCCs[i].begin(); it != setSCCs[i].end(); ++it  )
            {
                if( (fStart == true && it != setSCCs[i].begin())  || ( fStart == false && *it != *setSCCs[i].rbegin() ) )
                {
                    setCandidateRoots.erase( *it );
                }
            }

            if( fAllIn == false )
            {
                // also remove the first node
                if( fStart == true )
                {
                    setCandidateRoots.erase( *setSCCs[i].begin() );
                }
                else
                {
                    setCandidateRoots.erase( *setSCCs[i].rbegin() );
                }
            }

        }
    }
    
    //
    listPathRoots.clear();
    for( set<AbstractGraphNode *> :: iterator it = setCandidateRoots.begin(); it != setCandidateRoots.end(); ++it)
    {
        listPathRoots.push_back( *it );
    }
}



#if 0
void AbstractGraph :: FindComponents( const vector<AbstractGraphNode *> &listNodesTopOrder, set< set<AbstractGraphNode *> > &setComponents )
{
    // given a list of nodes (topologically sorted if DAG; otherwise from SCC), find components of the graph (weakly connected)
    // approach: assign an id for the component; during traversal, if it meets a known id, then assign all to the old id
    // end: all nodes w/ same id will be in the same component
    
    // set visited to false to everything in the list
    //for(int i=0; i<(int)listNodesTopOrder.size(); ++i)
    //{
    //    listNodesTopOrder[i]->SetVisited(false);
    //}
    
    setComponents.clear();
    map<AbstractGraphNode *, int> mapNodesToCompId;
    int idComponentNext = 1;
    for(int i=0; i<(int)listNodesTopOrder.size(); ++i)
    {
cout << "^^PRocessing node: " << listNodesTopOrder[i]->GetName() << ", idComponentNext = " << idComponentNext << endl;
        if( mapNodesToCompId.find( listNodesTopOrder[i] ) != mapNodesToCompId.end() )   // || listNodesTopOrder[i]->IsVisited() == true )
        {
            // this node is done
            continue;
        }
        //
        set<AbstractGraphNode *> setCompCurr;
        stack<AbstractGraphNode *> stackExp;
        stackExp.push( listNodesTopOrder[i] );
        int compIdExist = -1;
        while( stackExp.empty() == false )
        {
cout << "Top of stack: " << stackExp.top()->GetName() << endl;
            //
            AbstractGraphNode *pcurStackTop = stackExp.top();
            setCompCurr.insert( pcurStackTop );
            pcurStackTop->SetVisited(true);
            stackExp.pop();
            
            vector<AbstractGraphNode *> listNgbrs;
            pcurStackTop->GetNgbrNodes( listNgbrs );
            //while (itor.IsDone() == false)
            for(int jj=0; jj<(int)listNgbrs.size(); ++jj)
            {
                //
                if( mapNodesToCompId.find( listNgbrs[jj] ) != mapNodesToCompId.end()  )
                {
                    // reach a node that has already been discovered; don't exp from here
                    if( compIdExist > 0 &&  compIdExist != mapNodesToCompId[ listNgbrs[jj] ])
                    {
                        cout << "compIdExist: " << compIdExist << ", the neighbor " << listNgbrs[jj]->GetName() << " finds a component " << mapNodesToCompId[ listNgbrs[jj] ] << endl;
                        cout << "This should not happen: since node list is from SCC\n";
                        exit(1);
                    }
                    compIdExist = mapNodesToCompId[ listNgbrs[jj] ];
cout << "Find an existing component: " << compIdExist << endl;
                }
                else if( setCompCurr.find( listNgbrs[jj] ) == setCompCurr.end() )   // && listNgbrs[jj]->IsVisited() == true )
                {
                    // cointue explore it
                    stackExp.push( listNgbrs[jj] );
cout << "Pushing neighbor to stack: " << listNgbrs[jj]->GetName() << endl;
                }
            }
        }
        // save all disvoered node to be the component
        int compIdToUse = compIdExist;
        if( compIdToUse <  0)
        {
            compIdToUse = idComponentNext++;
        }
        for( set<AbstractGraphNode *> :: iterator it = setCompCurr.begin(); it != setCompCurr.end(); ++it )
        {
            //
            mapNodesToCompId.insert( map<AbstractGraphNode *,int> :: value_type( *it, compIdToUse ) );
        }
    }
    
    // all nodes w/ same id belong to the same component
    vector< set<AbstractGraphNode *> > listComponentId( idComponentNext+1 );
    for( map<AbstractGraphNode *,int> :: iterator it = mapNodesToCompId.begin(); it != mapNodesToCompId.end(); ++it )
    {
        listComponentId[ it->second].insert( it->first );
    }
    for( int i=0; i<(int)listComponentId.size(); ++i)
    {
        if( listComponentId[i].size() > 0 )
        {
            //
            setComponents.insert( listComponentId[i] );
        }
    }
}
#endif
