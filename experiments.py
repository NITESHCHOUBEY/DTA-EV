#!/usr/bin/python

# TODO: Use getopt or argparse to read command line arguments
import sys, time, os, numpy
from networkloading import *
from fixedPointAlgorithm import *

# For reading graphs
import networkx as nx

## EVEexample3 Network ##
def genEVExample3Network():
    G = Network()
    G.addNode("s")
    G.addNode("u")
    G.addNode("v")
    G.addNode("t")
    G.addEdge("s", "u", ExtendedRational(2), ExtendedRational(1), ExtendedRational(2))
    G.addEdge("s", "u", ExtendedRational(2), ExtendedRational(2), ExtendedRational(1))
    G.addEdge("u", "v", ExtendedRational(1), ExtendedRational(1), ExtendedRational(0))
    G.addEdge("v", "t", ExtendedRational(2), ExtendedRational(1), ExtendedRational(2))
    G.addEdge("v", "t", ExtendedRational(2), ExtendedRational(2), ExtendedRational(1))
    return G

def getEVExample3Paths(G: Network, s: Node, t: Node) -> List[Path]:
    # Paths for EVExample3 network
    p1 = Path([G.edges[0], G.edges[2], G.edges[4]])
    p2 = Path([G.edges[1], G.edges[2], G.edges[3]])
    p3 = Path([G.edges[1], G.edges[2], G.edges[4]])
    return [p1, p2, p3]


## Leon's example network ##
def genLeonNetwork():
    G = Network()
    G.addNode("s")
    G.addNode("u")
    G.addNode("v")
    G.addNode("t")
    G.addNode("b")
    G.addNode("a")
    G.addEdge("s", "u", ExtendedRational(3,2), ExtendedRational(2))
    G.addEdge("u", "v", ExtendedRational(1), ExtendedRational(2))
    G.addEdge("v", "t", ExtendedRational(1,2), ExtendedRational(2))
    G.addEdge("s", "a", ExtendedRational(3), ExtendedRational(3))
    G.addEdge("a", "v", ExtendedRational(3), ExtendedRational(3))
    G.addEdge("u", "b", ExtendedRational(3), ExtendedRational(3))
    G.addEdge("b", "t", ExtendedRational(3), ExtendedRational(3))
    return G

def getLeonsPaths(G: Network, s: Node, t: Node) -> List[Path]:
    # Paths for the example in Leon's thesis
    p1 = Path([G.edges[0], G.edges[1], G.edges[2]])
    p2 = Path([G.edges[3], G.edges[4], G.edges[2]])
    p3 = Path([G.edges[0], G.edges[5], G.edges[6]])
    return [p1, p2, p3]


## Nguyen network ##
def genNguyenNetwork():
    G = Network()
    # Nodes
    G.addNode("1")
    G.addNode("2")
    G.addNode("3")
    G.addNode("4")
    G.addNode("5")
    G.addNode("6")
    G.addNode("7")
    G.addNode("8")
    G.addNode("9")
    G.addNode("10")
    G.addNode("11")
    G.addNode("12")
    G.addNode("13")
    # Edges
    G.addEdge("1", "12", ExtendedRational(5,6), ExtendedRational(150))
    G.addEdge("1", "5", ExtendedRational(5,6),  ExtendedRational(75 ))
    G.addEdge("12", "6", ExtendedRational(5,6), ExtendedRational(75 ))
    G.addEdge("12", "8", ExtendedRational(5,6), ExtendedRational(150))
    G.addEdge("4", "5", ExtendedRational(5,6),  ExtendedRational(150))
    G.addEdge("5", "6", ExtendedRational(5,6),  ExtendedRational(150))
    G.addEdge("6", "7", ExtendedRational(5,6),  ExtendedRational(150))
    G.addEdge("7", "8", ExtendedRational(5,6),  ExtendedRational(150))
    G.addEdge("4", "9", ExtendedRational(5,6),  ExtendedRational(225))
    G.addEdge("5", "9", ExtendedRational(5,6,), ExtendedRational(75 ))
    G.addEdge("6", "10", ExtendedRational(5,6), ExtendedRational(75 ))
    G.addEdge("7", "11", ExtendedRational(5,6), ExtendedRational(75 ))
    G.addEdge("8", "2", ExtendedRational(5,6),  ExtendedRational(75 ))
    G.addEdge("9", "10", ExtendedRational(5,6), ExtendedRational(150))
    G.addEdge("10", "11", ExtendedRational(5,6),ExtendedRational(150))
    G.addEdge("11", "2", ExtendedRational(5,6), ExtendedRational(150))
    G.addEdge("9", "13", ExtendedRational(5,6), ExtendedRational(225))
    G.addEdge("11", "3", ExtendedRational(5,6), ExtendedRational(75 ))
    G.addEdge("13", "3", ExtendedRational(5,6), ExtendedRational(150))
    return G

def getNguyenPaths(G: Network, s: Node, t: Node) -> List[Path]:
    # Paths for Nguyen network
    pathlist = []
    if int(s.name)==1:
        if int(t.name)==2:
            pathlist.append(Path([G.edges[1-1], G.edges[4-1], G.edges[13-1]]))
            pathlist.append(Path([G.edges[1-1], G.edges[3-1], G.edges[7-1], G.edges[8-1], G.edges[13-1]]))
            pathlist.append(Path([G.edges[1-1], G.edges[3-1], G.edges[7-1], G.edges[12-1], G.edges[16-1]]))
            pathlist.append(Path([G.edges[1-1], G.edges[3-1], G.edges[11-1], G.edges[15-1], G.edges[16-1]]))

            pathlist.append(Path([G.edges[2-1], G.edges[6-1], G.edges[7-1], G.edges[8-1], G.edges[13-1]]))
            pathlist.append(Path([G.edges[2-1], G.edges[6-1], G.edges[7-1], G.edges[12-1], G.edges[16-1]]))
            pathlist.append(Path([G.edges[2-1], G.edges[6-1], G.edges[11-1], G.edges[15-1], G.edges[16-1]]))
            pathlist.append(Path([G.edges[2-1], G.edges[10-1], G.edges[14-1], G.edges[15-1], G.edges[16-1]]))
        elif int(t.name)==3:
            pathlist.append(Path([G.edges[1-1], G.edges[3-1], G.edges[7-1], G.edges[12-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[1-1], G.edges[3-1], G.edges[11-1], G.edges[15-1], G.edges[18-1]]))

            pathlist.append(Path([G.edges[2-1], G.edges[6-1], G.edges[7-1], G.edges[12-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[2-1], G.edges[6-1], G.edges[11-1], G.edges[15-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[2-1], G.edges[10-1], G.edges[14-1], G.edges[15-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[2-1], G.edges[10-1], G.edges[17-1], G.edges[19-1]]))
    elif int(s.name)==4:
        if int(t.name)==2:
            pathlist.append(Path([G.edges[5-1], G.edges[6-1], G.edges[7-1], G.edges[8-1], G.edges[13-1]]))
            pathlist.append(Path([G.edges[5-1], G.edges[6-1], G.edges[7-1], G.edges[12-1], G.edges[16-1]]))
            pathlist.append(Path([G.edges[5-1], G.edges[6-1], G.edges[11-1], G.edges[15-1], G.edges[16-1]]))
            pathlist.append(Path([G.edges[5-1], G.edges[10-1], G.edges[14-1], G.edges[15-1], G.edges[16-1]]))

            pathlist.append(Path([G.edges[9-1], G.edges[14-1], G.edges[15-1], G.edges[16-1]]))
        elif int(t.name)==3:
            pathlist.append(Path([G.edges[5-1], G.edges[6-1], G.edges[7-1], G.edges[12-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[5-1], G.edges[6-1], G.edges[11-1], G.edges[15-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[5-1], G.edges[10-1], G.edges[14-1], G.edges[15-1], G.edges[18-1]]))

            pathlist.append(Path([G.edges[9-1], G.edges[14-1], G.edges[15-1], G.edges[18-1]]))
            pathlist.append(Path([G.edges[9-1], G.edges[17-1], G.edges[19-1]]))

    return pathlist



def readArgs(argv):
    # Arguments passed
    print("\nName of script:", sys.argv[0])
    n = len(sys.argv)
    print("Total arguments passed:", n)
    # print("\nArguments passed:", end = " ")
    # for i in range(1, n):
        # print(sys.argv[i], end = " ")
    # print()
    print(sys.argv)

    # Read the instance name
    insName = argv[0]

    # arglist.append(argv[i] for i in range(1,n))
    # exclude = {","," ","[","]"}
    argList = []
    for i in range(1,n):
        argList.append(argv[i])
        # print(i,arglist[i-1])
    # for i in range(len(arglist)):
        # if arglist[i] not in exclude:
            # args.append(arglist[i])
            # print(i,arglist[i])
    return argList


def readNetwork(edgeList, verbose: bool=False) -> Network:
    #TODO: put checks for a valid network
    # First read as a MultiDiGraph
    Gn = nx.MultiDiGraph()
    Gn = nx.read_edgelist(edgeList, comments='#', nodetype=str,\
            create_using=nx.MultiDiGraph, data=(("nu", ExtendedRational),\
            ("tau", ExtendedRational), ("ec", ExtendedRational),))
    if verbose: print('edges: ', list(Gn.edges(data=True)))
    if verbose: print('nodes: ', list(Gn.nodes()))

    # Convert to a Network object for our purposes
    G = Network()
    for node in Gn.nodes():
        G.addNode(node)
    for u,v,data in Gn.edges(data=True):
        G.addEdge(u,v,data['nu'],data['tau'],data['ec'])

    #TODO: Plot the graph using nx
    return G


def readCommodities(commList) -> List[Tuple[Node, Node, PWConst]]:
    commodities = []
    with open(commList, 'r') as fobj:
        for line in fobj:
            # print(line)
            data = [entry for entry in line.split()]
            # print(data, len(data)-2, data[2:len(data)-1])
            # Create the PWConst function for this commodity
            times = [ExtendedRational(i) for i in data[2:2+math.ceil(len(data)/2 -1)]]
            vals = [ExtendedRational(i) for i in data[2+len(times):len(data)]]
            # pwcList = [ExtendedRational()]
            # TODO: Autosimplify=0 has to be passed otherwise error occurs in
            # networkloaading()
            pwcf = PWConst(times, vals, 0)
            print(pwcf)
            commodities.append((G.getNode(data[0]), G.getNode(data[1]), pwcf))
    print('comm: ', commodities)
    return commodities


if __name__ == "__main__":

    # TODO: Read as part of arglist
    energyBudget = 3

    argList = readArgs(sys.argv)

    G = readNetwork(argList[0])

    commodities = readCommodities(argList[1])

    fname = ""
    for i in range(len(argList)-1):
        fname += argList[i] + "_"
    fname += argList[-1]

    # Read arguments into required variables
    [insName,timeHorizon,maxIter,precision,alpha,timeStep] = argList[2:len(argList)]
    # print("read as: insName,timeHorizon,maxIter,precision,alpha,timeStep")
    [insName,timeHorizon,maxIter,precision,alpha,timeStep] = [str(insName),\
            ExtendedRational(timeHorizon),int(maxIter),float(precision),\
            ExtendedRational(alpha),ExtendedRational(timeStep)]
    print("check args: ",insName,timeHorizon,maxIter,precision,alpha,timeStep)

    # if insName == "leon":
        # G = genLeonNetwork()
        # # pathList = getLeonsPaths(G,G.getNode("s"),G.getNode("t"))
        # pathList = G.findPaths(G.getNode("s"),G.getNode("t"))
    # elif insName == "evExample3":
        # # G = genEVExample3Network()
        # # pathList = getEVExample3Paths(G,G.getNode("s"),G.getNode("t"))
        # pathList = G.findPaths(G.getNode("s"),G.getNode("t"), energyBudget)
    # elif insName == "nguyen":
        # G = genNguyenNetwork()
        # pathList = []
    # else:
        # print("Unknown network name. Exiting.")
        # exit(0)

    # Find list of paths for each commodity
    # TODO: put data checks
    pathList = []
    for i,(s,t,u) in enumerate(commodities):
        print("i ", i,s,t,u)
        pathList.append(G.findPaths(s, t, energyBudget))
        # for p in pathList[i]:
            # print(p)

    # Start
    tStart = time.time()
    # filename = 'alpha=%.2f'%alpha + '_maxIter=%d'%maxIter +\
    # '_timeStep=%.2f'%timeStep + '_precision=%.2f'%precision + '.npz'
    # print("------------------------------------------------\n",filename,\
            # "\n------------------------------------------------")
    # if insName == "nguyen":
        # f, alphaIter, absDiffBwFlowsIter, relDiffBwFlowsIter, travelTime, stopStr,\
        # alphaStr, qopiIter = fixedPointAlgo(G, pathList, precision, [\
                # (G.getNode("1"), G.getNode("2"),PWConst([0,10,20],[3,0],0)),\
                # (G.getNode("1"), G.getNode("3"),PWConst([0,10,20],[3,0],0)),\
                # (G.getNode("4"), G.getNode("2"),PWConst([0,10,20],[3,0],0)),\
                # (G.getNode("4"), G.getNode("3"),PWConst([0,10,20],[3,0],0))],\
                # timeHorizon,maxIter,timeStep,alpha,True)
    # else:
        # comm = [(G.getNode("s"),G.getNode("t"), PWConst([0,10],[3],0))]
        # f, alphaIter, absDiffBwFlowsIter, relDiffBwFlowsIter, travelTime, stopStr,\
        # alphaStr, qopiIter = fixedPointAlgo(G, pathList,\
                # # precision,[(G.getNode("s"),G.getNode("t"), PWConst([0,10,50],[3,0],0))],\
                # # precision,[(G.getNode("s"),G.getNode("t"), PWConst([0,10,20],[3,0],0))],\
                # precision,comm,\
                # # precision,[(G.getNode("s"),G.getNode("t"), PWConst([0,10],[3],0))],\
                # timeHorizon, maxIter, timeStep, alpha, True)
    f, alphaIter, absDiffBwFlowsIter, relDiffBwFlowsIter, travelTime, stopStr,\
            alphaStr, qopiIter = fixedPointAlgo(G, pathList, precision, commodities,\
                    timeHorizon, maxIter, timeStep, alpha, True)

    tEnd = time.time()
    print("travelTimes: ", travelTime)
    eventualFlow = networkLoading(f)
    print("eventualFlow: ", eventualFlow)
    print("f: ", f)
    print("queue at: ")
    for id, e in enumerate(eventualFlow.network.edges):
        print("edge %d: "%id, e, eventualFlow.queues[e])

    # alpha and flowDiff
    ralphaIter = [round(float(b),3) for b in alphaIter]
    rAbsDiffBwFlowsIter = [round(float(b),3) for b in absDiffBwFlowsIter]
    rRelDiffBwFlowsIter = [round(float(b),3) for b in relDiffBwFlowsIter]
    rqopiIter = [round(float(b),3) for b in qopiIter]
    print("alpha ", ralphaIter)
    print("absDiffBwFlowsIter ", rAbsDiffBwFlowsIter)
    print("relDiffBwFlowsIter ", rRelDiffBwFlowsIter)
    print("qopiIter ", rqopiIter)

    print("Termination message: ", stopStr)
    print("\nElasped wall time: ", round(tEnd-tStart,4))

    # Save the results to files
    # dirname = os.path.expanduser('./npzfiles')
    dirname = os.path.expanduser('./miscfiles')
    fname += '_' + alphaStr.replace('/','By')
    numpy.savez(os.path.join(dirname, fname),G=G,f=f,eventualFlow=eventualFlow,time=tEnd-tStart,\
            alphaIter=alphaIter,absDiffBwFlowsIter=absDiffBwFlowsIter,\
            relDiffBwFlowsIter=relDiffBwFlowsIter,travelTime=travelTime,\
            stopStr=stopStr,alphaStr=alphaStr,qopiIter=qopiIter)
    print("output saved to file: %s.npz"%fname)
