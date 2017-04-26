//
//  GraphUtils.h
//  
//
//  Created by Yufeng Wu on 3/7/15.
//
//

#ifndef ____GraphUtils__
#define ____GraphUtils__

#include <vector>
#include <set>
#include <string>
#include <stack>
#include <map>
#include <iostream>

using namespace std;

// ******************************************************************
// Iterator through ngbrs

class AbstractGraphNode;

class AbstractGraphNodeNgbrIterator
{
public:
    AbstractGraphNodeNgbrIterator( AbstractGraphNode *pNodeToVisit);
    void First();
    void Next();
    bool IsDone() const;
    AbstractGraphNode *GetCurNgbr() const;
    
private:
    AbstractGraphNode *pNodeToVisit;
    int indexNgbr;
};

// ******************************************************************
// Iterator through a list of nodes (w/ prioorty)

class AbstractGraphNodeListPrioriIterator
{
public:
    AbstractGraphNodeListPrioriIterator( vector<AbstractGraphNode *> &listNodes, map<AbstractGraphNode *, int> &mapNodeOrder );
    void First();
    void Next();
    bool IsDone() const;
    AbstractGraphNode *GetCurNode() const;
    
private:
    vector<AbstractGraphNode *> &listNodesIn;
    map<AbstractGraphNode *, int> &mapNodeOrderIn;
    int indexPri;
    vector<int> listOrderIndices;
};



// ******************************************************************
// abstract class for graph node

class GraphEdge;

class AbstractGraphNode
{
    friend class AbstractGraphNodeNgbrIterator;
    
public:
    AbstractGraphNode() : fFlagVisit(false) {  }
    AbstractGraphNode(const string &name) : fFlagVisit(false), nodeName(name) {  }
    virtual ~AbstractGraphNode();
    //virtual int GetId() = 0;
    virtual string GetName() const {return nodeName;}
    virtual void SetName(const string &name) { nodeName = name; }
    virtual void AddNgbrNodeBasic(AbstractGraphNode *ngbr);
    virtual void AddEdge(GraphEdge *pe) { listEdges.push_back(pe); }
    virtual void GetNgbrNodes( vector<AbstractGraphNode *> &listNgbrs );
    virtual bool IsNeighbor(AbstractGraphNode *pnode);
    virtual bool IsVisited() const { return fFlagVisit; }
    virtual void SetVisited(bool f) { fFlagVisit = f; }
    virtual int GetNumAdjNodes() const { return listEdges.size(); }
    virtual void PreVisit( AbstractGraphNode *pSrcNode ) {}
    virtual void PostVisit( ) {}
    virtual int GetNumEdges() const { return listEdges.size(); }
    virtual GraphEdge *GetEdgeTo( AbstractGraphNode *pDest );
    virtual double GetTravPriority() const { return 1.0; }
    virtual void Dump() const;
    
private:
    vector<GraphEdge *> listEdges;
    bool fFlagVisit;
    string nodeName;      // often, a node will have an assigned name
};


// ******************************************************************
// graph node w/ reference to something

class GraphNodeRef : public AbstractGraphNode
{
public:
    GraphNodeRef() : refToSomething(NULL) {}
    void SetRef( void *refptr ) { refToSomething = refptr; }
    void *GetRef() const { return refToSomething; }
    
private:
    void *refToSomething;
};

// ******************************************************************
// graph node w/ info to extend the current node

class GraphNodeRefExt;

class GraphNodeRefExt : public GraphNodeRef
{
public:
    GraphNodeRefExt() {}
    void SetExt( const string &strExt ) { strExtInfo = strExt; }
    string GetExt() const { return strExtInfo; }
    void AddNgbrRef( GraphNodeRefExt *pNgbr, void *pref, const string &extInfo );
    void AddNgbrRef( GraphNodeRefExt *pNgbr, void *pref, const string &extInfo, double len );
    
private:
    string strExtInfo;
};

// ******************************************************************
// graph edge
class GraphEdge
{
public:
    GraphEdge( AbstractGraphNode *pSourceIn, AbstractGraphNode *pDestIn ) : pSource(pSourceIn), pDest(pDestIn), pathLen(0.0) {}
    GraphEdge( const GraphEdge &rhs) : pSource(rhs.pSource), pDest(rhs.pDest), pathLen(rhs.pathLen) {}
    AbstractGraphNode *GetSource() const { return pSource; }
    AbstractGraphNode *GetDest() const { return pDest; }
    double GetLength() const { return pathLen; }
    void SetLength(double len) { pathLen = len; }

private:
    AbstractGraphNode *pSource;
    AbstractGraphNode *pDest;
    double pathLen;
};


// ******************************************************************
// graph edge w/ reference and string

class GraphEdgeREfExt : public GraphEdge
{
public:
    GraphEdgeREfExt(GraphNodeRefExt *pSource, GraphNodeRefExt *pDest) : GraphEdge(pSource, pDest), pref(NULL) {}
    GraphEdgeREfExt(const GraphEdgeREfExt &rhs) : GraphEdge(rhs.GetSource(), rhs.GetDest() ), pref(rhs.pref), ext(rhs.ext) {}
    void SetRef(void *prefIn) { pref = prefIn; }
    void *GetRef() const { return pref; }
    void SetExt(const string &extIn) { ext = extIn; }
    string GetExt() const { return ext; }
    //void AddNgbrRef( GraphNodeRefExt *pNgbr, void *pref, const string &extInfo );
    
private:
    void *pref;
    string ext;
};


// ******************************************************************
// abstract class for graph 

class GraphNodeSCC;

class AbstractGraph
{
public:
    AbstractGraph() {}
    virtual ~AbstractGraph();
    int GetNumNodes() const { return listGraphNodes.size(); }
    int GetNumEdges() const;
    void Reset();
    void AddNode( AbstractGraphNode *pn );
    void InitBFS();
    void BFS();
    void FindSimplePaths( set<vector<AbstractGraphNode *> > &setPaths,  int maxPathLen, int maxNodeOccurInPath  );
    void FindSimplePathsBoundedLength( set<vector<AbstractGraphNode *> > &setPaths, int maxPathLen );
    void FindSimplePathsTopSort( set<vector<AbstractGraphNode *> > &setPaths, int maxNodeOccurInPath );
    void SCC( vector< set<AbstractGraphNode *> > &setSCCs );
    void OutputGML(const char *filename);
    //void FindComponents( const vector<AbstractGraphNode *> &listNodesTopOrder, set< set<AbstractGraphNode *> > &setComponents );
    
private:
    void SCCFrom( GraphNodeSCC *pnodeCurr, int &indexNext, stack<GraphNodeSCC *> &stackNodesToExp, map< AbstractGraphNode *, GraphNodeSCC *> &mapOrigNodeToSCCNode, vector< set<AbstractGraphNode *> > &setSCCs );
    void BFSFrom( AbstractGraphNode *pnodeStart );
    void ExtendSimplePathAt( set<vector<AbstractGraphNode *> > &setPaths, int maxPathLen, int maxNodeOccurInPath  );
    //void RemoveSubPathsFrom( set<vector<AbstractGraphNode *> > &setPaths );
    //void RemoveSubPathsFromFast( set<vector<AbstractGraphNode *> > &setPaths );
    void FindSimplePathsTopSortFrom(AbstractGraphNode *pnodeCurr, const vector<AbstractGraphNode *> &listNodesExp, const vector<AbstractGraphNode *> &listPathEnds, set<vector<AbstractGraphNode *> > &setPaths);
    void FindSimplePathsTopSortStart(bool fStart, const vector< set<AbstractGraphNode *> > &setSCCs, vector<AbstractGraphNode *> &listPathRoots);
    
    vector<AbstractGraphNode *> listGraphNodes;
};

// ******************************************************************
// Ammended node class for SCC finding node

class GraphNodeSCC : public AbstractGraphNode
{
public:
    GraphNodeSCC( AbstractGraphNode *pNodeOrigIn) : pNodeOrig(pNodeOrigIn), index(-1), lowlink(-1), onStack(false)  {}
    int GetIndex() const { return index; }
    void SetIndex(int idn) { index = idn; }
    int GetLowLink() const { return lowlink; }
    void SetLowLink(int idn) { lowlink = idn; }
    bool IsOnStack() const { return onStack; }
    void SetOnStack(bool f) { onStack = f; }
    AbstractGraphNode *GetNodeOrig() const { return pNodeOrig; }
    virtual void Dump() const { cout << "index: " << index << ", lowlink: " << lowlink<< " "; pNodeOrig->Dump(); cout << " "; }
    
    // now create delegation
    virtual string GetName() const {return pNodeOrig->GetName();}
    virtual void SetName(const string &name) { pNodeOrig->SetName( name ); }
    virtual void AddNgbrNodeBasic(AbstractGraphNode *ngbr) { pNodeOrig->AddNgbrNodeBasic(ngbr); }
    virtual void AddEdge(GraphEdge *pe) { pNodeOrig->AddEdge(pe); }
    virtual void GetNgbrNodes( vector<AbstractGraphNode *> &listNgbrs ) { pNodeOrig->GetNgbrNodes(listNgbrs); }
    virtual bool IsNeighbor(AbstractGraphNode *pnode) { return pNodeOrig->IsNeighbor(pnode); }
    virtual bool IsVisited() const { return pNodeOrig->IsVisited(); }
    virtual void SetVisited(bool f) { pNodeOrig->SetVisited(f); }
    virtual int GetNumAdjNodes() const { return pNodeOrig->GetNumAdjNodes(); }
    virtual void PreVisit( AbstractGraphNode *pSrcNode ) {}
    virtual void PostVisit( ) {}
    virtual int GetNumEdges() const { return pNodeOrig->GetNumEdges(); }
    virtual GraphEdge *GetEdgeTo( AbstractGraphNode *pDest ) { return pNodeOrig->GetEdgeTo(pDest); }
    virtual double GetTravPriority() const { return pNodeOrig->GetTravPriority();  }
    
private:
    AbstractGraphNode *pNodeOrig;
    int index;
    int lowlink;
    bool onStack;
};




#endif /* defined(____GraphUtils__) */
