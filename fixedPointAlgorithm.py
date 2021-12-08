from typing import List, Dict, Tuple

from networkloading import *
from dynamic_dijkstra import dynamic_dijkstra

# For the root finding problem
from scipy import optimize

def findShortestSTpath(s: Node, t: Node, flow: PartialFlow, time: ExtendedRational) -> Path:
    (arrivalTimes, realizedTimes) = dynamic_dijkstra(time, s, t, flow)
    p = Path([], t)
    while p.getStart() != s:
        v = p.getStart()
        for e in v.incoming_edges:
            if e in realizedTimes and arrivalTimes[e.node_from] + realizedTimes[e] == arrivalTimes[v]:
                p.add_edge_at_start(e)
                break

    return p

def setInitialPathFlows(commodityId: int, G: Network, s: Node, t: Node,	u: PWConst, zeroflow: PartialFlow, pathInflows: PartialFlowPathBased) -> PartialFlowPathBased:
    print("To be implemented! Passing hardcoded path inflows.")
    # print("Setting up the shortest path and other paths.")
    p1 = Path([G.edges[0], G.edges[2], G.edges[4]])
    p2 = Path([G.edges[1], G.edges[2], G.edges[3]])
    # printPathInNetwork(p2,G)
    p3 = Path([G.edges[1], G.edges[2], G.edges[4]])
    # printPathInNetwork(p3,G)
    # exit(0)
    # pathInflows.setPaths(commodityId, [findShortestSTpath(s, t, zeroflow,
	    # ExtendedRational(0)),p2,p3], [u,PWConst([0,50],[0],0),PWConst([0,50],[0],0)])
    pathInflows.setPaths(commodityId, [p1,p2,p3], [u,PWConst([0,50],[0],0),PWConst([0,50],[0],0)])
    # print("Setting up path ", p2)
    # pathInflows.setPaths(commodityId, p2, 0)
    # print("Setting up path ", p3)
    # pathInflows.setPaths(commodityId, p3, 0)
    return pathInflows


def getAllSTpaths(G: Network, s: Node, t: Node, flow: PartialFlow) -> List[Path]:
    print("To be implemented: passing hardcoded paths for now.")
    p1 = Path([G.edges[0], G.edges[2], G.edges[4]])
    p2 = Path([G.edges[1], G.edges[2], G.edges[3]])
    p3 = Path([G.edges[1], G.edges[2], G.edges[4]])
    pathList = [p1,p2,p3]
    # print("Path list :", pathList)
    # print("Freeflow travel times :", (sum(t) for t in p1.edges.tau)
    return pathList


def fixedPointUpdate(oldPathInflows: PartialFlowPathBased, timeHorizon:
        ExtendedRational, alpha: float, timestepSize, u: PWConst, verbose: bool) -> PartialFlowPathBased:
    currentFlow = networkLoading(oldPathInflows, timeHorizon)

    # TODO: Adjust during algorithm?
    # timestepSize = ExtendedRational(1,1)
    threshold = ExtendedRational(1,100)

    newPathInflows = PartialFlowPathBased(oldPathInflows.network, oldPathInflows.getNoOfCommodities())

    for i in range(oldPathInflows.getNoOfCommodities()):
        flowValue = [None]*len(oldPathInflows.fPlus[i])
        travelTime = [None]*len(oldPathInflows.fPlus[i])
        if verbose: print("Considering commodity ", i)
        newPathInflows.setPaths(i,[P for P in oldPathInflows.fPlus[i]],[PWConst([ExtendedRational(0,1)],[],ExtendedRational(0,1)) for P in oldPathInflows.fPlus[i]])
        s = newPathInflows.sources[i]
        t = newPathInflows.sinks[i]
        theta = ExtendedRational(0,1)
        # We subdivide the time into intervals of length timestepSize
        while theta < oldPathInflows.getEndOfInflow(i):
        # while theta < oldPathInflows.getEndOfInflow(i) and theta < 48:
            # For each subinterval i we determine the dual variable v_i
            # (and assume that it will stay the same for the whole interval)
            # if verbose: print("timeinterval [", theta, ",", theta+timestepSize,"]")

	    # Set up the update problem for each subinterval
            print("\nSetting up the update problem for subinterval [", theta,
                    ",", theta+timestepSize,"]\n")
	    # Get path travel times for this subinterval
            for j,P in enumerate(oldPathInflows.fPlus[i]):
                 fP = oldPathInflows.fPlus[i][P]
                 # converting to float (optimize.root does not work with fractions)
                 print("theta arg ", theta)
                 travelTime[j] = float(currentFlow.pathArrivalTime(P, theta) - theta)
                 flowValue[j] = float(fP.getValueAt(theta))
                 print("Path: P",j, "flowValue: ", flowValue[j], "travelTime: ",\
                         travelTime[j], "at theta =", theta, "fp: ", fP)

            # Find integral value, ubar, of (piecewise constant) function u in this
            # subinterval
            uval1 = u.getValueAt(theta)
            uval2 = u.getValueAt(theta + timestepSize)
            if uval1==uval2:
                ubar = uval1*timestepSize
            else:
                if not uval2 == 0:
                    print("Implement me: find integral of u when it has different\
                            positive values within a subinterval")
                    exit(0)

            # TODO: Find a good starting point
            # A trivial guess: assume all terms to be positive and solve for the dual variable
            x0 = ((-sum(flowValue) + alpha*sum(travelTime))*timestepSize +
                    ubar)/(len(flowValue)*timestepSize)
            print("x0 ", round(x0,2))
            # exit(0)
            sol = optimize.root(findDualVar, x0, (alpha, flowValue, travelTime,
                    timestepSize, ubar))
                    # timestepSize, ubar), method='broyden1')
            if not sol.success:
                print(sol.message)
                exit(0)
            else:
                print(sol)
            # print("currentFlow ", currentFlow)
            for j,P in enumerate(oldPathInflows.fPlus[i]):
                newFlowVal = max(flowValue[j] - alpha*travelTime[j] + sol.x[0], 0)
                # print("newFlowVal ", newFlowVal)
                newPathInflows.fPlus[i][P].addSegment(ExtendedRational(theta +
                    timestepSize), ExtendedRational(newFlowVal))
            print("newPathInflows: ", newPathInflows)
            theta = theta + timestepSize
    exit(0)
    return newPathInflows

def findDualVar(x, alpha, flowValue, travelTime, timestepSize, ubar):
    # print("printing args ", x,alpha,flowValue,travelTime,timestepSize,ubar)
    termSum = 0
    # print("termSum ", termSum)
    for j,fv in enumerate(flowValue):
        # print("print terms for j", j, " : ", flowValue[j],alpha,travelTime[j],x,timestepSize)
        # termSum += max(flowValue[j] - alpha*travelTime[j] + x[0], 0)*timestepSize
        termSum += max(flowValue[j] - alpha*travelTime[j] + x, 0)*timestepSize
        # print("termSum in loop ", termSum)
    # print("result ", termSum - ubar)
    # exit(0)
    return termSum - ubar


def differenceBetweenPathInflows(oldPathInflows : PartialFlowPathBased, newPathInflows : PartialFlowPathBased) -> ExtendedRational:
    assert (oldPathInflows.getNoOfCommodities() == newPathInflows.getNoOfCommodities())
    #TODO: Also check if the time horizon for both the pathInflows is same or not

    difference = ExtendedRational(0)

    for i in range(oldPathInflows.getNoOfCommodities()):
        for path in oldPathInflows.fPlus[i]:
            if path in newPathInflows.fPlus[i]:
                difference += (oldPathInflows.fPlus[i][path] + newPathInflows.fPlus[i][path].smul(ExtendedRational(-1,1))).norm()
            else:
                difference += oldPathInflows.fPlus[i][path].norm()
        for path in newPathInflows.fPlus[i]:
            if path not in oldPathInflows.fPlus[i]:
                difference += newPathInflows.fPlus[i][path].norm()

    return difference

# Function arguments: (network, precision, List[source node, sink node, ?], time
# horizon, maximum allowed number of iterations, verbosity on/off)
def fixedPointAlgo(N : Network, precision : float, commodities :
        List[Tuple[Node, Node, PWConst]], timeHorizon:
        ExtendedRational=math.inf, maxSteps: int = None, timeStep: int = None,
        alpha : float = None, verbose : bool = False) -> PartialFlowPathBased:
    step = 0

    ## Init:
    # Create zero-flow (PP: why?)
    pathInflows = PartialFlowPathBased(N,0)
    zeroflow = networkLoading(pathInflows,timeHorizon)

    i = 0
    pathInflows = PartialFlowPathBased(N, len(commodities))
    # print(pathInflows)
    # Initial flow: For every commodity, select the shortest s-t path and send
    # all flow along this path (and 0 flow along all other paths)
    for (s,t,u) in commodities:
        # pathInflows.setPaths(i, [findShortestSTpath(s, t, zeroflow, ExtendedRational(0))], [u])
        # print("u ", u)
        setInitialPathFlows(i, N, s, t, u, zeroflow, pathInflows)
        print(pathInflows)
        i += 1

    if verbose: print("Starting with flow: \n", pathInflows)

    ## Iteration:
    while maxSteps is None or step < maxSteps:
        if verbose: print("Starting iteration #", step)
        # newpathInflows = networkLoading(pathInflows,timeHorizon)
        # print("newpathInflows ", newpathInflows)
        newpathInflows = fixedPointUpdate(pathInflows, timeHorizon, alpha,
                timeStep, u, verbose)
        # if differenceBetweenPathInflows(pathInflows,newpathInflows) < precision:
            # print("Attained required precision!")
            # return newpathInflows
        if verbose: print("Changed amount is ", differenceBetweenPathInflows(pathInflows,newpathInflows))
        if verbose: print("Current flow is\n", newpathInflows)
        pathInflows = newpathInflows
        step += 1

    print("Maximum number of steps reached without attaining required precision!")
    return pathInflows
