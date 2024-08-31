from typing import List, Dict, Tuple

from networkloading import *
from dynamic_dijkstra import dynamic_dijkstra

# For the root finding problem
from scipy import optimize
import sys, time
import numpy as np
import copy

# def findShortestSTpath(s: Node, t: Node, flow: PartialFlow, time: number) -> Path:
#     (arrivalTimes, realizedTimes) = dynamic_dijkstra(time, s, t, flow)
#     p = Path([], t)
#     while p.getStart() != s:
#         v = p.getStart()
#         for e in v.incoming_edges:
#             if e in realizedTimes and arrivalTimes[e.node_from] + realizedTimes[e] == arrivalTimes[v]:
#                 p.add_edge_at_start(e)
#                 break
#     return p


# def findShortestFeasibleSTpath(time: number, s: Node, t: Node, flow:
#         PartialFlow, budget: number) -> Path:
#     (arrivalTimes, realizedTimes) = dynamicFeasDijkstra(time, s, t, flow, budget)
#     for i in enumerate(arrivalTimes):
#         print("times ", i, i[0], i[1])
#     for i in enumerate(realizedTimes):
#         print("times ", i, i[0], i[1])
#     p = Path([], t)
#     while p.getStart() != s:
#         v = p.getStart()
#         for e in v.incoming_edges:
#             if e in realizedTimes and arrivalTimes[e.node_from] + realizedTimes[e] == arrivalTimes[v]:
#                 p.add_edge_at_start(e)
#                 break

#     print("shortest ST path ", printPathInNetwork(p, flow.network))
#     return p

#countIterations=0

class PrioritizedItem:
    # constructing the class to create the first element for priority and store all other elements in items 
    def __init__(self, priority, item):
        self.priority = priority
        self.item = item

    # defining behaviour of the less than operator so that it only uses the priority to compare the objects of this class
    def __lt__(self, other):
        return self.priority < other.priority

    # defining behaviour of equal to regarding the same 
    def __eq__(self, other):
        return self.priority == other.priority


def calculate_cost(energy: float,price: float,time: float,alpha:float,priceToTime:float)->float:
    beta=0
    #alpha=1
    #priceToTime=0
    return alpha*(time+priceToTime*price)+beta*energy


# def dijkstra(nt: dict, src: Node, dest: Node, excluded_paths: set, EB: float, PB: float,alpha:float,priceToTime:float) -> Path:
#     # Priority queue containing (cost, time , node, path as a list, total_energy spent till that node, total_price spent till that node)
#     pq = PriorityQueue()
#     pq.put(PrioritizedItem(0, (0, src, [], 0, 0, {src: 1}))) 

#     # Dictionary to store the shortest path to each node
#     shortest_paths = {src: (0,0, [], 0, 0)}  # cost(priority), time , path_edges, energy spent, price spent

#     while not pq.empty():
#         current_item = pq.get()
#         combined_cost = current_item.priority
#         current_time, current_node, current_path, current_energy, current_price, visited_nodes = current_item.item

#         if current_node in shortest_paths and combined_cost>shortest_paths[current_node][0]:
#             continue
#         # arrived at the destination 
#         if current_node == dest:
#             potential_path = Path(current_path, src)                            # creating the path
#             if potential_path not in excluded_paths:                            # making sure we have a new path
#                 if current_energy <= EB and current_price <= PB:       
#                     return potential_path
#                 else:
#                     continue

#         for edge in current_node.outgoing_edges:
#             neighbor = edge.node_to
           
#             if neighbor in visited_nodes:
#                 if visited_nodes[neighbor] >= 2:
#                     continue  # Skippint the node if its neighbor has already been visited twice
#                 visited_nodes[neighbor] += 1
#             else:
#                 visited_nodes[neighbor] = 1

#             updated_time = current_time + nt[edge]
#             updated_energy = current_energy + edge.ec
#             updated_price = current_price + edge.price
#             updated_path = current_path + [edge]
#             updated_cost=calculate_cost(updated_energy,updated_price,updated_time,alpha,priceToTime)

#             # if we have a new node or if we are getting a faster route we update the dictionary shortest path and pq
#             if (neighbor not in shortest_paths or updated_cost < shortest_paths[neighbor][0]) and updated_energy <= EB and updated_price <= PB:
#                 shortest_paths[neighbor] = (updated_cost , updated_time, updated_path, updated_energy, updated_price)
#                 pq.put(PrioritizedItem(updated_cost, (updated_time, neighbor, updated_path, updated_energy, updated_price,visited_nodes)))

#             visited_nodes[neighbor]-=1
#             if visited_nodes[neighbor]==0:
#                 del visited_nodes[neighbor]

#     # If no path is found we return none
#     return None

def dijkstra(nt: dict, src: Node, dest: Node, excluded_paths: set, EB: float, PB: float, alpha: float, priceToTime: float) -> Path:
    # Priority queue containing (cost, time , node, path as a list, total_energy spent till that node, total_price spent till that node)
    pq = PriorityQueue()
    pq.put(PrioritizedItem(0, (0, src, [], 0, 0, {src: 1}))) 

    # Dictionary to store all feasible paths to each node
    all_paths = {src: [(0, 0, [], 0, 0)]}  # cost(priority), time , path_edges, energy spent, price spent

    while not pq.empty():
        current_item = pq.get()
        combined_cost = current_item.priority
        current_time, current_node, current_path, current_energy, current_price, visited_nodes = current_item.item

        # Skip if a better path has already been found
        # if current_node in all_paths and any(combined_cost > path[0] for path in all_paths[current_node]):
        #     continue

        # Arrived at the destination 
        if current_node == dest:
            potential_path = Path(current_path, src)  # creating the path
            if potential_path not in excluded_paths:  # making sure we have a new path
                if current_energy <= EB and current_price <= PB:
                    return potential_path
                else:
                    continue

        for edge in current_node.outgoing_edges:
            neighbor = edge.node_to

            if neighbor in visited_nodes:
                if visited_nodes[neighbor] >= 2:
                    continue  # Skipping the node if its neighbor has already been visited twice
                visited_nodes[neighbor] += 1
            else:
                visited_nodes[neighbor] = 1

            updated_time = current_time + nt[edge]
            updated_energy = current_energy + edge.ec
            updated_price = current_price + edge.price
            updated_path = current_path + [edge]
            updated_cost = calculate_cost(updated_energy, updated_price, updated_time, alpha, priceToTime)

            # adding new paths if they are feasible
            if updated_energy <= EB and updated_price <= PB:
                if neighbor not in all_paths:
                    all_paths[neighbor] = []
                all_paths[neighbor].append((updated_cost, updated_time, updated_path, updated_energy, updated_price))
                pq.put(PrioritizedItem(updated_cost, (updated_time, neighbor, updated_path, updated_energy, updated_price, visited_nodes)))

            visited_nodes[neighbor] -= 1
            if visited_nodes[neighbor] == 0:
                del visited_nodes[neighbor]

    # If no path is found we return none
    return None


class PreserveReprWrapper:
    def __init__(self, obj):
        self.obj = obj
        self._repr = repr(obj) 
    
    def __repr__(self):
        return self._repr
    
    def unwrap(self):
        return self.obj
    
    def __getattr__(self, attr):
        return getattr(self.obj, attr)
    
    def __len__(self):
        return len(self.obj)
    
    def __iter__(self):
        return iter(self.obj)
    
    def __eq__(self, other):
        if isinstance(other, PreserveReprWrapper):
            return self.obj == other.obj
        return self.obj == other

    def __hash__(self):
        return hash(self.obj)
    
    def __deepcopy__(self, memo):
        return PreserveReprWrapper(copy.deepcopy(self.obj, memo))


def custom_copy_partialflowpathbased(oldPathInflows):
    
    new_flow = PartialFlowPathBased(oldPathInflows.network, oldPathInflows.noOfCommodities)
    
    # Copy fPlus with the keys wrapped to preserve repr
    for i in range(oldPathInflows.noOfCommodities):
        new_flow.fPlus[i] = {PreserveReprWrapper(path): copy.deepcopy(pwconst)
                             for path, pwconst in oldPathInflows.fPlus[i].items()}
    
    # Deep copy sources and sinks
    new_flow.sources = copy.deepcopy(oldPathInflows.sources)
    new_flow.sinks = copy.deepcopy(oldPathInflows.sinks)
    
    # Deep copy u (network inflow rates)
    new_flow.u = [copy.deepcopy(u) for u in oldPathInflows.u]
    
    return new_flow




def nextShortestPath(nt:dict, oldPathInFlowsCommodity: dict, src: Node, dest: Node, EB: float = infinity, PB: float = infinity, alpha: float=1 ,priceToTime:float=0) -> Path:
    # creating a set of old path to make sure we return a new path only
    # global countIterations
    existing_paths = set(oldPathInFlowsCommodity.keys())

    # running Dijkstra's algorithm to find the next shortest path not in existing_paths and meeting EB and PB constraints
    new_path = dijkstra(nt, src, dest, existing_paths, EB, PB,alpha,priceToTime)

    # If no path is found, return None
    if new_path is None:
        return None
    
    
    return new_path


def fixedPointUpdate(N:Network,currentFlow: PartialFlow, oldPathInflows: PartialFlowPathBased, timeHorizon:
        number, alpha: float, timestepSize, priceToTime: float, commodities, verbose:
        bool,generatedPath:List[List[Path]],EB: float=infinity,PB: float=infinity) -> PartialFlowPathBased:

    # for i,_ in enumerate(commodities):
    #     original_paths = list(oldPathInflows.fPlus[i].keys())
    #     print("initial ", original_paths)
    newPathInflows = PartialFlowPathBased(oldPathInflows.network, oldPathInflows.getNoOfCommodities())

    nt={e:e.tau for e in N.edges}


    for i,comd in enumerate(commodities):

        original_oldPathInflows = custom_copy_partialflowpathbased(oldPathInflows)

        flowValue = [None]*len(oldPathInflows.fPlus[i])
        travelTime = [None]*len(oldPathInflows.fPlus[i])
        price = [None]*len(oldPathInflows.fPlus[i])
        if False: print("Considering commodity ", i)
        newPathInflows.setPaths(i,[P for P in oldPathInflows.fPlus[i]],[PWConst([zero],[],zero) for P in oldPathInflows.fPlus[i]])
        s = newPathInflows.sources[i]
        t = newPathInflows.sinks[i]
        theta = zero
        meanIter = 0

        k = -1


        fl=newPathInflows.fPlus[i]
        npaths=[]
        nff=[]

        for p,pc in fl.items():
            npaths.append(p)
            nff.append(pc)
        
        paths = list(oldPathInflows.fPlus[i].keys())
        flow_function = list(oldPathInflows.fPlus[i].values())

        
        itera=1


        while theta < oldPathInflows.getEndOfInflow(i):
            k += 1
            # For each subinterval i we determine the dual variable v_i
            # (and assume that it will stay constant for the whole interval)
            # if verbose: print("timeinterval [", theta, ",", theta+timestepSize,"]")
            # print(i," ",itera)

            # for keyed, plst in currentFlow.queues.items():
            #     print(f"Key: {keyed} (Address: {repr(keyed)}) -> Value: {plst}")
            
            # print("Key in dictionary:", repr(list(currentFlow.queues.keys())[0]))

            for p in oldPathInflows.fPlus[i]:
                for ed in p.edges:
                    nt[ed]=currentFlow.c(ed,theta)
            itera+=1
          #  if(countIterations>30):
           #     print("check4")
            nsp=nextShortestPath(nt,oldPathInflows.fPlus[i],s,t,EB,PB,alpha,priceToTime)
            #if countIterations>30:
             #   print("check5") 
            if nsp!=None:

                # print(i," ",nsp," check ",theta)

                updated_oldpathInflow=PartialFlowPathBased(N,len(commodities))
                generatedPath[i].append(nsp)
                paths.append(nsp)
                # flow_function.append(PWConst([0, comd[4].segmentBorders[-1]], [0], 0))
                flow_function.append(PWConst([zero],[],zero))
                updated_oldpathInflow.setPaths(i,paths,flow_function)
                oldPathInflows=updated_oldpathInflow

                travelTime.append(None)
                flowValue.append(None)
                price.append(None)
                
                npaths.append(nsp)
                nff.append(PWConst([zero],[],zero))

                
                newPathInflows.addPath(i,nsp,PWConst([zero],[],zero))

            # else:
            #     print("No new feasible path")

	    # Set up the update problem for each subinterval
            # Get path travel times for this subinterval
            for j,P in enumerate(oldPathInflows.fPlus[i]):
                 fP = oldPathInflows.fPlus[i][P]
                 # converting to float (optimize.root does not work with fractions)
                 travelTime[j] = float(currentFlow.pathArrivalTime(P,\
                     theta + timestepSize/2) - (theta + timestepSize/2))
                 flowValue[j] = float(fP.getValueAt(theta))
                 price[j] = P.getPrice()

            # Find integral value, ubar, of (piecewise constant) function u in this
            # subinterval
            ubar = comd[4].integrate(theta, theta + timestepSize)
            
            # for j,P in enumerate(newPathInflows.fPlus[i]):
            #      fP = newPathInflows.fPlus[i][P]
            #      print("for path ",N.printPathInNetwork(P)," flow is ",fP)

            
            # print("Now redis")
            # For adjusting alpha
            # uval = ubar/timestepSize
            # alpha = uval/min(travelTime)
            # alpha = 0.5*uval*(1/min(travelTime) + 1/max(travelTime))
            # alpha = 0.5*uval*(1/max(travelTime))
            # alpha = uval/max(travelTime)
            # adjAlpha = [uval/min(travelTime) for j,_ in enumerate(flowValue)]
            # alphaList.append(adjAlpha[0])
            # print('adjAlpha: ', theta + timestepSize/2, uval, [round(i, 4) for
                # i in travelTime], [round(i,4) for i in adjAlpha])
            # adjAlpha = [alpha*flowValue[j]/(2*travelTime[j]) for j,_ in enumerate(flowValue)]
            # adjAlpha = [flowValue[j]/(2*travelTime[j]) for j,_ in enumerate(flowValue)]

            # adjmin = min([i for i in adjAlpha if i > 0])
            # adjmax = max([i for i in adjAlpha if i > 0])
            # adjAlpha = [adjmin if j == 0 else j for j in adjAlpha]
            # adjAlpha = [adjmax if j == 0 else j for j in adjAlpha]

            # TODO: Find a good starting point
            # A trivial guess: assume all terms to be positive and solve for the dual variable
            # TODO: adjust for price
            x0 = ((-sum(flowValue) + alpha*sum(travelTime))*timestepSize +
                    ubar)/(len(flowValue)*timestepSize)
            # optimize.show_options(solver='root', method='broyden1', disp=True)
            # TODO: Find a way to run optimize.root quietly
            bracketLeft = 0
            # TODO: adjust for price
            # bracketRight = abs(max(list(map(float.__sub__, list(map(lambda x: alpha*x,
                # travelTime)), flowValue)))) + ubar + 1
            bracketRight = 0
            for j,_ in enumerate(travelTime):
                bracketRight += max(bracketRight, -flowValue[j] +
                        alpha*(travelTime[j] + priceToTime*price[j] )) + ubar + 1

            # Newton's method using a routine that return value and derivative
            # TODO: pass an aggregation function value based on travel time and price
            # of paths
            # sol = optimize.root_scalar(dualVarRootFuncComb, (adjAlpha, flowValue, travelTime,
            if priceToTime == 0:
                sol = optimize.root_scalar(dualVarRootFuncComb, (alpha, flowValue, travelTime,
                    timestepSize, ubar), x0=x0, bracket=[bracketLeft, bracketRight],
                    fprime=True, method='newton')
            else:
                sol = optimize.root_scalar(dualVarRootFuncCombPrice, (alpha, flowValue,
                    travelTime, priceToTime, price, timestepSize, ubar), x0=x0,
                    bracket=[bracketLeft, bracketRight], fprime=True,
                    method='newton')

            if not sol.converged:
                print("The optimize.root_scalar() method has failed to converge due to the following reason:")
                print("\"", sol.flag, "\"")
                exit(0)
                # alpha = alpha/2 + uval/(2*max(travelTime))
                sol = optimize.root_scalar(dualVarRootFuncComb, (alpha, flowValue, travelTime, price, priceToTime,
                    timestepSize, ubar), x0=x0, bracket=[bracketLeft, bracketRight],
                    fprime=True, method='newton')
                if not sol.converged:
                    print("Adjusted alpha! The optimize.root_scalar() method has still failed to converge because:")
                    print("\"", sol.flag, "\"")
                    exit(0)
                else:
                    meanIter += sol.iterations
            else:
                meanIter += sol.iterations

            
            for j,P in enumerate(oldPathInflows.fPlus[i]):
                # CAUTION: Price term to be included here
                newFlowVal = max(flowValue[j] - alpha*(travelTime[j] + priceToTime*price[j]) + sol.root, 0)
                # print(f" for path {N.printPathInNetwork(P)} flow is {flowValue[j]:.2f} travel time {travelTime[j]:.2f} ans sol ,{sol.root:.2f} finally newflowval {newFlowVal:.2f}")

                newPathInflows.fPlus[i][P].addSegment(makeNumber(theta+timestepSize), makeNumber(newFlowVal))

            theta = theta + timestepSize

        tmpVar = max(timestepSize,1/timestepSize)
        if False: print("Mean # of root.scalar() iterations ",\
                float(round(meanIter/(tmpVar*oldPathInflows.getEndOfInflow(i)),2)),\
                " for ", tmpVar*oldPathInflows.getEndOfInflow(i), " subintervals")
            
        oldPathInflows = custom_copy_partialflowpathbased(original_oldPathInflows)
    
    return newPathInflows, alpha


def dualVarRootFunc(x, alpha, flowValue, travelTime, timestepSize, ubar):
    termSum = 0
    for j,fv in enumerate(flowValue):
        termSum += max(flowValue[j] - alpha*travelTime[j] + x, 0)*timestepSize
    return float(termSum - ubar)


def dualVarRootFuncGrad(x, alpha, flowValue, travelTime, timestepSize, ubar):
    termSum = 0
    for j,fv in enumerate(flowValue):
        if (flowValue[j] - alpha*travelTime[j] + x) > 0:
            termSum += timestepSize
    return float(termSum)


def dualVarRootFuncComb(x, alpha, flowValue, travelTime, timestepSize, ubar):
    # TODO: Read as input with each commodity
    termSum = 0
    gradTermSum = 0
    for j,fv in enumerate(flowValue):
        tmp = flowValue[j] - alpha*travelTime[j] + x
        if tmp > 0:
            termSum += tmp*timestepSize
            gradTermSum += timestepSize
    return float(termSum - ubar), float(gradTermSum)


def dualVarRootFuncCombPrice(x, alpha, flowValue, travelTime, priceToTime, price, timestepSize, ubar):
    termSum = 0
    gradTermSum = 0
    for j,fv in enumerate(flowValue):
        tmp = flowValue[j] - alpha*(travelTime[j] + priceToTime*price[j]) + x
        if tmp > 0:
            termSum += tmp*timestepSize
            gradTermSum += timestepSize
    return float(termSum - ubar), float(gradTermSum)


def differenceBetweenPathInflows(oldPathInflows : PartialFlowPathBased, newPathInflows : PartialFlowPathBased) -> number:
    assert (oldPathInflows.getNoOfCommodities() == newPathInflows.getNoOfCommodities())
    #TODO: Also check if the time horizon for both the pathInflows is same or not

    difference = zero

    for i in range(oldPathInflows.getNoOfCommodities()):
        for path in oldPathInflows.fPlus[i]:
            if path in newPathInflows.fPlus[i]:
                difference += (oldPathInflows.fPlus[i][path] + newPathInflows.fPlus[i][path].smul(-1)).norm()
            else:
                difference += oldPathInflows.fPlus[i][path].norm()
        for path in newPathInflows.fPlus[i]:
            if path not in oldPathInflows.fPlus[i]:
                difference += newPathInflows.fPlus[i][path].norm()

    return difference

# Find the sum of norm (integration) of path inflow fPlus functions
def sumNormOfPathInflows(pathInflows : PartialFlowPathBased) -> number:
    sumNorm = zero
    for i in range(pathInflows.getNoOfCommodities()):
        for path in pathInflows.fPlus[i]:
            sumNorm += (pathInflows.fPlus[i][path]).norm()
    # TODO: This should be equal to the integration of u (if required, put an assert)
    return sumNorm

# Function arguments: (network, precision, List[source node, sink node, ?], time
# horizon, maximum allowed number of iterations, verbosity on/off)
# TODO: warm-start using an available path flow?
def fixedPointAlgo(N : Network, pathList : List[Path], precision : float, commodities :
        List[Tuple[Node, Node, PWConst]], timeHorizon:
        number=infinity, maxSteps: int = None, timeLimit: int = infinity, timeStep: int = None,
        alpha : float = None, priceToTime : float = None,EB: float=infinity,PB: float=infinity, verbose : bool = False) -> PartialFlowPathBased:
    tStartAbs = time.time()
    step = 0

    ## Initialize:
    # Create zero-flow (PP: why?)
    pathInflows = PartialFlowPathBased(N,0)
    # TODO: Conform with LG if this can be removed
    # zeroflow = networkLoading(pathInflows)

    pathInflows = PartialFlowPathBased(N, len(commodities))
    # Initial flow: For every commodity, select the shortest s-t path and send
    # all flow along this path (and 0 flow along all other paths)
    for i,(s,t,_,_,u) in enumerate(commodities):
        flowlist = [PWConst([0,u.segmentBorders[-1]],[0],0)]*(len(pathList[i])-1)
        flowlist.insert(0,u)
        pathInflows.setPaths(i, pathList[i], flowlist)

    if False: print("Starting with flow: \n", pathInflows)

    oldAbsDiffBwFlows = infinity
    oldRelDiffBwFlows = infinity
    gamma = makeNumber(1)
    alphaIter = []
    absDiffBwFlowsIter = []
    relDiffBwFlowsIter = []
    travelTime = []
    qopiIter = []  # qualityOfPathInflows
    qopiMeanIter = []  # mean of qopi
    qopiFlowIter = []  # qopii per unit flow per unit time
    qopiPathComm = []  # mean of qopi
    shouldStop = not (maxSteps is None or step < maxSteps)


    # alphaStr = ''
    # alphaStr = 'uByTmin'
    # alphaStr = 'uBy2Tmin'
    # alphaStr = 'meanUByTminTmax'
    # alphaStr = 'uByTmax'
    # alphaStr = r'($\gamma$)'
    # alphaStr = r'($\gamma\alpha$)'
    alphaStr = r'expoSmooth($\gamma$)'
    # alphaStr = r'expoSmooth($\gamma/2$)'
    # alphaStr = r'relExpoSmooth($\gamma/2$)'
    # alphaStr = r'min2ExpoSmooth($\gamma/2$)'

    totDNLTime = 0
    totFPUTime = 0
    tStart = time.time()
    iterFlow = networkLoading(pathInflows)
    totDNLTime += time.time()-tStart
    print("\nTime taken in networkLoading(): ", round(time.time()-tStart,4))

    ## Iteration:
    while not shouldStop:
        if verbose: print("STARTING ITERATION #", step)
        tStart = time.time()
        # TODO: Read priceToTime as input per commodity
        genPaths = [[] for _ in range(len(commodities))]
        # print(f"i: {i}, len(generatedPath): {len(commodities)}")
        # for i,_ in enumerate(commodities):
        #     original_paths = list(pathInflows.fPlus[i].keys())
        #     print("verify lengths:",original_paths)
        newpathInflows, alpha = fixedPointUpdate(N,iterFlow, pathInflows, timeHorizon, alpha,
                timeStep, priceToTime, commodities, verbose,genPaths,EB,PB)
        
        # original_paths = list(newpathInflows.fPlus[i].keys())
        # print("Original path lengths:", [len(P) for P in original_paths])
        for i in range(newpathInflows.getNoOfCommodities()):
            unwrapped_fPlus = {}
            for path, flow in newpathInflows.fPlus[i].items():
                if isinstance(path, PreserveReprWrapper):
                    unwrapped_fPlus[path.unwrap()] = flow
                else:
                    unwrapped_fPlus[path] = flow
            newpathInflows.fPlus[i] = unwrapped_fPlus


        # print("new path is")
        # for i in range(len(commodities)):
        #     genPaths[i]=list(set(genPaths[i]))
        #     for p in genPaths[i]:
        #         print(N.printPathInNetwork(p))
        #     print(" ends here \n")

        
     #   print(" new inflows \n ",newpathInflows,"\n new inflows end")

        # updated_newpathInflow=PartialFlowPathBased(N,len(commodities))

        # for i in range(len(commodities)):
        #     curr_flow=newpathInflows.fPlus[i]
        #     paths=[]
        #     flow_function=[]
        #     for path,pwc in curr_flow.items():
        #         paths.append(path)
        #         flow_function.append(pwc)
            
        #     for p in genPaths[i]:
        #         paths.append(p)
        #         flow_function.append(PWConst([0, u.segmentBorders[-1]], [0], 0))

        #     updated_newpathInflow.setPaths(i,paths,flow_function)

        # newpathInflows=updated_newpathInflow

        # print("updated \t",updated_newpathInflow)

        totFPUTime += time.time() - tStart
        print("\nTime taken in fixedPointUpdate(): ", round(time.time()-tStart,4))

        newAbsDiffBwFlows = differenceBetweenPathInflows(pathInflows,newpathInflows)
        newRelDiffBwFlows = newAbsDiffBwFlows/sumNormOfPathInflows(pathInflows)

        # Check Stopping Conditions
        if newAbsDiffBwFlows < precision:
            shouldStop = True
            stopStr = "Attained required (absolute) precision!"
        # elif newRelDiffBwFlows < precision/10:
            # shouldStop = True
            # stopStr = "Attained required (relative) precision!"
        elif not (maxSteps is None or step < maxSteps):
            shouldStop = True
            stopStr = "Maximum number of steps reached!"

        elif (time.time() - tStartAbs > timeLimit):
            shouldStop = True
            stopStr = "Maximum time limit reached!"

        qopi = infinity
        qopiMean = infinity
        qopiFlow = infinity
        if not shouldStop:
            # Update Alpha
            if newAbsDiffBwFlows == 0:
                gamma = 0
            else:
                if step > 0: gamma = 1 - abs(newAbsDiffBwFlows - oldAbsDiffBwFlows)/(newAbsDiffBwFlows +
                        oldAbsDiffBwFlows)

            # Alpha Update Rule
            # oldalpha = alpha
            # alpha = gamma # equal to factor
            # alpha = gamma*alpha # multiplied by factor
            alpha = (gamma)*(gamma*alpha) + (1-gamma)*alpha # expo smooth using gamma
            # alpha = (0.5*gamma)*(0.5*gamma*alpha) + (1-0.5*gamma)*alpha # expo smooth using gamma/2
            # alpha = max(0.2, (0.5*gamma)*(0.5*gamma*alpha) + (1-0.5*gamma)*alpha) # expo smooth using gamma/2
            # if step > 1 and oldalpha == alpha:
                # print('Changing alpha')
                # alpha = alpha*(0.5) #step/maxSteps

            # Measure quality of the path inflows
            tStart = time.time()
            iterFlow = networkLoading(newpathInflows)
            tEnd = time.time()
            totDNLTime += tEnd - tStart
            print("\nTime taken in networkLoading(): ", round(tEnd-tStart,4))

            qopi = 0
            qopiFlow = 0
            for i,comd in enumerate(commodities):
                qopiInt = 0
                if False: print('comm ', i)
                fP = newpathInflows.fPlus[i]
                theta = zero
                oldqopi = np.zeros(len(newpathInflows.fPlus[i]))
                while theta < newpathInflows.getEndOfInflow(i):
                    tt = np.empty(len(newpathInflows.fPlus[i]))
                    for j,P in enumerate(newpathInflows.fPlus[i]):
                        tt[j] = iterFlow.pathArrivalTime(P,theta + timeStep/2) - (theta + timeStep/2)
                    tmin = min(tt)
                    fval = []
                    for j,P in enumerate(newpathInflows.fPlus[i]):
                        val = fP[P].getValueAt(theta + timeStep/2)
                        fval.append(val)
                        newqopi = (tt[j] - tmin)*val
                        qopi += newqopi
                        qopiInt += ((newqopi-oldqopi[j])/2 + min(newqopi, oldqopi[j]))*timeStep/tmin
                        oldqopi[j] = newqopi
                    theta = theta + timeStep
                # Integrate also over the interval [finalTheta, T]
                for j,_ in enumerate(newpathInflows.fPlus[i]):
                    qopiInt += ((oldqopi[j]-0)/2)*timeStep/tmin

                commFlow = comd[4].integrate(comd[4].segmentBorders[0], comd[4].segmentBorders[-1])
                qopiFlow += qopiInt/commFlow

            if verbose: print("Norm of change in flow (abs.) ", round(float(newAbsDiffBwFlows),4),\
                    " previous change ", round(float(oldAbsDiffBwFlows),4), " alpha ",\
                    round(float(alpha),4), ' qopi ', round(qopi,4), ' qopiFlow ', round(qopiFlow,4))
            if verbose: print("Norm of change in flow (rel.) ", round(float(newRelDiffBwFlows),4),\
                    " previous change ", round(float(oldRelDiffBwFlows),4))

            # update iteration variables
            # for i,_ in enumerate(commodities):
            #     original_paths = list(pathInflows.fPlus[i].keys())
            #     print("pathinflow initially:", original_paths)

            #     original_paths = list(newpathInflows.fPlus[i].keys())
            #     print("newpath inflow:", original_paths)

            pathInflows = newpathInflows

            # for i,_ in enumerate(commodities):
            #     original_paths = list(pathInflows.fPlus[i].keys())
            #     print("newpathinfowcopy :", original_paths)

            oldAbsDiffBwFlows = newAbsDiffBwFlows
            oldRelDiffBwFlows = newRelDiffBwFlows

        qopiIter.append(qopi)
        qopiMeanIter.append(qopiMean)
        qopiFlowIter.append(qopiFlow)
        # alphaIter.append(alphaList)
        alphaIter.append(alpha)
        absDiffBwFlowsIter.append(newAbsDiffBwFlows)
        relDiffBwFlowsIter.append(newRelDiffBwFlows)

        step += 1

    print(stopStr)
    # Find path travel times for the final flow
    finalFlow = networkLoading(pathInflows, verbose=False)
    for i,comd in enumerate(commodities):
        fP = newpathInflows.fPlus[i]
        ttravelTime = np.empty([len(pathInflows.fPlus[i]),\
                math.ceil(pathInflows.getEndOfInflow(i)/timeStep)])
        qopiPath = np.empty([len(pathInflows.fPlus[i]),\
                math.ceil(pathInflows.getEndOfInflow(i)/timeStep)])
        theta = zero
        k = -1
        commFlow = comd[4].integrate(comd[4].segmentBorders[0], comd[4].segmentBorders[-1])
        while theta < pathInflows.getEndOfInflow(i):
            k += 1
            for j,P in enumerate(pathInflows.fPlus[i]):
                 ttravelTime[j][k] = finalFlow.pathArrivalTime(P,\
                     theta + timeStep/2) - (theta + timeStep/2)
            tmin = np.min(ttravelTime[:,k])
            for j,P in enumerate(pathInflows.fPlus[i]):
                
                # if(i>0):
                #     P=repr(P)
                # print(i," ", j," c ",P," ",theta)
                val = fP[P.unwrap() if isinstance(P, PreserveReprWrapper) else P].getValueAt(theta + timeStep/2)
                # val = fP[P].getValueAt(theta + timeStep/2)
                # print(f"Path {j}, flow value: {val}")
                qopiPath[j][k] = (ttravelTime[j][k] - tmin)*val/(tmin*commFlow)
            theta = theta + timeStep
        travelTime.append(ttravelTime)
        qopiPathComm.append(qopiPath)
    return pathInflows, alphaIter, absDiffBwFlowsIter, relDiffBwFlowsIter,\
            travelTime, stopStr, alphaStr, qopiIter, qopiFlowIter,\
            qopiPathComm, totDNLTime, totFPUTime
