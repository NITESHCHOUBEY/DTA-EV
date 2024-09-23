from __future__ import annotations

import random
from typing import Union, List

from utilities import *
from collections import deque
from queue import PriorityQueue
import itertools

# for parallel computations
from joblib import Parallel, delayed
import sys

globPathList = []
# A directed edge (from https://github.com/Schedulaar/predicted-dynamic-flows/blob/main/predictor/src/core/graph.py)
# with capacity nu, travel time tau, and energy consumption ec
class Edge:
    node_from: Node
    node_to: Node
    tau: number
    nu: number
    ec: number
    price: number

    def __init__(self, node_from: Node, node_to: Node, capacity: number=1,
            traveltime: number=1, energyCons: number=0, price: number=0):
        # Creating an edge from node_from to node_to
        self.node_from = node_from
        self.node_to = node_to

        # Free flow travel time over this edge
        assert(traveltime >= 0)
        self.tau = traveltime

        # Capacity of this edge
        assert(capacity >= 0)
        self.nu = capacity

        # Energy consumption over this edge (can be negative)
        self.ec = energyCons

        # Price of travelling over this edge
        assert(price >= 0)
        self.price = price

    def __str__(self):
        return "("+str(self.node_from)+","+str(self.node_to)+")"

    # Get edge (id) in network
    def getIdInNetwork(self, G: Network):
        return "Implement me."

# A node (from https://github.com/Schedulaar/predicted-dynamic-flows/blob/main/predictor/src/core/graph.py )
class Node:
    name: str
    # TODO: rename the parameter name id
    id: int
    incoming_edges: List[Edge]
    outgoing_edges: List[Edge]

    def __init__(self, name: str, id:int):
        # Create a node with name name and without incoming or outgoing edges
        self.name = name
        self.id = id
        self.incoming_edges = []
        self.outgoing_edges = []

    def __str__(self):
        # Print the nodes name
        if self.name == "":
            return str(self.id)
        else:
            return str(self.name)

class PrioritizedItem:
    # constructing the class to create the first element for priority and store all other elements in items 
    def __init__(self, priority, item):
        self.priority = priority
        self.item = item

    # defining behaviour of the less than operator so that it only uses the priority to compare the objects of this class
    def __lt__(self, other):
        return self.priority < other.priority

    # defining behaviour of equal
    def __eq__(self, other):
        return self.priority == other.priority
    
class Network:
    edges: List[Edge]
    nodes: List[Node]
    idCounter: int

    def __init__(self):
        # Create an empty network
        self.edges = []
        self.nodes = []
        self.idCounter = 0

    def getNode(self,node: Union[str,int,Node]) -> Union[Node,None]:
        if isinstance(node,Node):
            return node
        elif isinstance(node, str):
            for v in self.nodes:
                if v.name == node:
                    return v
        elif isinstance(node,int):
            # TODO: Schnelleres Verfahren fur Suche nach IDs?
            for v in self.nodes:
                if v.id == node:
                    return v

        return None

    def addNode(self,name: str="") -> Node:
        self.nodes.append(Node(name,self.idCounter))
        self.idCounter += 1
        return self.nodes[-1]

    def addEdge(self,node_from: Union[str,int,Node], node_to: Union[str,int,Node],
            nu: number, tau: number, ec: number=zero, price: number=zero):
        v = self.getNode(node_from)
        w = self.getNode(node_to)
        e = Edge(v, w, nu, tau, ec, price)
        v.outgoing_edges.append(e)
        w.incoming_edges.append(e)
        self.edges.append(e)

    def removeEdge(self,edge:Edge):
        edge.node_to.incoming_edges.remove(edge)
        edge.node_from.outgoing_edges.remove(edge)
        self.edges.remove(edge)

    def subdivide(self,edge:Edge, nu: number, tau: number,
            ec: number):
        self.edges.remove(edge)
        v = self.addNode()
        self.addEdge(edge.node_from, v, nu, tau, ec)
        self.addEdge(v, edge.node_to, nu, tau, ec)

    def duplicate(self, edge:Edge, nu: number, tau: number, ec: number):
        self.addEdge(edge.node_from, edge.node_to, nu, tau, ec)

    def __str__(self) -> str:
        s = "Graph with " + str(len(self.nodes)) + " nodes and " + str(len(self.edges)) + " edges:\n"
        for v in self.nodes:
            s += str(v) + ": {"
            for e in v.outgoing_edges:
                s += str(e) + ", "
            s = s[:-2]
            if len(v.outgoing_edges) > 0:
                s += "}\n"
            else:
                s += "\n"
        return s

    # Print path p based on edge ids in the network G
    def printPathInNetwork(self, p: Path):
        s = str()
        for e in p.edges:
            if len([i for i in e.node_from.outgoing_edges\
                    if i.node_to == e.node_to]) > 1:
                s += str(self.edges.index(e))
            s += str(e)
        return s

    # Join two (sub)paths
    def joinPaths(self, p1: Path, p2: Path) -> Path:
        if p1.getEnd() == p2.getStart():
            path = Path(p1.edges)
            for e in p2.edges:
                path.add_edge_at_end(e)
                # print('edge: ', e, end=' ')
                # print('path: ', self.printPathInNetwork(path))
            return path
        else:
            print('Unable to join paths: ', self.printPathInNetwork(p1),
                    self.printPathInNetwork(p2))
            exit(0)
            return None

    # Concatenate two sets of (sub)paths while checking feasibility
    def joinFeasiblePaths(self, P1: List[Path], P2: List[Path], EB: number=infinity,\
            PB: number=infinity) -> List[Path]:
        pathList = []
        for p1 in P1:
            # print('for p1: ', self.printPathInNetwork(p1))
            for p2 in P2:
                # if first edge is a self loop, join (p1,p2)
                if p2.edges[0].node_from == p2.edges[0].node_to:
                    pathList.append(self.joinPaths(p1, p2))
                # else check feasibility
                else:
                    if (p1.getEnergyConsump() + p2.getEnergyConsump() <= EB) and \
                            (p1.getPrice() + p2.getPrice() <= PB):
                        pathList.append(self.joinPaths(p1, p2))
        return pathList

    # Find acyclic energy-feasible paths in the network from source src to
    # destination dest with/without excluding nodes with self loop
    # args: source (src), destination (dest), energy budget (EB), price budget (PB)
    # TODO: Finding paths with cycles in the network
    # This function can find only those paths in which a node is visited just once
    def findPaths(self, src, dest, EB: number=infinity, PB: number=infinity,
            excludeSelfLoopNodes: bool=False, verbose: bool=False) -> List[Path]:
        if excludeSelfLoopNodes:
            # Find nodes with self loops
            selfLoopNodes = [e.node_from for e in self.edges if e.node_from == e.node_to]
        else:
            selfLoopNodes = []
        # print('Nodes with self loop: ', *(n for n in selfLoopNodes))

        # Queue to store (partial) paths
        q = deque()

        # Add edges going out of the source to q
        for e in self.edges:
            if e.node_from == src:
                q.append(Path([e]))

        # List to store the final paths
        pathList = []
        count = 0

        while q:
            count += 1
            if verbose: print("\ncount:%d"%count)

            # Get the (earliest generated) partial path
            if verbose: print("q before pop")
            for p in q:
                if verbose: print(printPathInNetwork(p, self))
            path = q.popleft()
            # print("after pop ", printPathInNetwork(path, self), q)
            for p in q:
                if verbose: printPathInNetwork(p, self)

            # Get the last node in the partial path
            last = path.edges[-1].node_to
            if verbose: print("last ", last)

            # If the last node is the destination node then store the path
            if last == dest:
                if verbose: print("Found s-t Path:", printPathInNetwork(path, self))
                pathList.append(path)

            # Traverse all the nodes connected to the current node and push new partial
            # path to queue
            edgeListCurrNode = [e for e in self.edges if (e.node_from == last and
                e.node_to not in selfLoopNodes)]
            if verbose: print("edgeListCurrNode ", len(edgeListCurrNode), edgeListCurrNode)
            for e in edgeListCurrNode:
                if verbose: print("edge %d" %self.edges.index(e), e,
                        printPathInNetwork(path, self), path.isNodeInPath(e.node_to),
                        path.getEnergyConsump(), e.ec, (path.getEnergyConsump() + e.ec <= EB))
                if (not path.isNodeInPath(e.node_to)) and (path.getEnergyConsump() + e.ec <= EB)\
                        and (path.getPrice() + e.price <= PB):
                    newpath = Path(path.edges)
                    if verbose: print("newpath before append ", printPathInNetwork(newpath,self))
                    newpath.add_edge_at_end(e)
                    q.append(newpath)
                    if verbose: print("newpath after append ", printPathInNetwork(newpath, self))

        # Print pathList
        if verbose:
            print("\nTotal %d paths found from node %s to node %s:"%(len(pathList),src,dest))
            for i,p in enumerate(pathList):
                print(i, len(p), printPathInNetwork(p, self))
        return pathList

    # TODO: Avoid creating duplicate paths (those without recharging edges)
    # Removes duplicate paths from a list of paths.
    def removeDuplicatePaths(self, pathList: List[Path]) -> List[Path]:
        # print('Removing duplicate paths.')
        noDupSet = set()
        # dupes = []
        dupCount = 0
        for x in pathList:
            if x not in noDupSet:
                noDupSet.add(x)
            else:
                dupCount += 1
                # dupes.append(x)
        # print('duplicate: ', *((self.printPathInNetwork(p), p.getPrice(),
            # p.getNetEnergyConsump(), p.getFreeFlowTravelTime()) for p in dupes),
            # sep='\n')
        print('Removed %d duplicate paths.'%dupCount)
        return list(dict.fromkeys(noDupSet))


    # Check if a path has a subpath in a list
    def hasSubpath(self, s: Path) -> bool:
        for p in globPathList:
            if frozenset(s.edges) > frozenset(p.edges):
                return True
        return False


    # Removes the dominated paths (ones that will not have positive flows due to
    # presence of proper subpaths with a lower latency
    # For full recharging, proper superset paths are dominated
    def removeDominatedPathsPar(self) -> List[Path]:
        print('Removing dominated paths in parallel.')
        results = Parallel(n_jobs=int(sys.argv[10]), require='sharedmem', verbose=1)\
                (delayed(self.hasSubpath)(s) for s in globPathList)
        nonDomPathList = [globPathList[i] for i,j in enumerate(results) if not j]
        print('Number of nondominated paths found: %d'%len(nonDomPathList))
        return nonDomPathList


    # Removes the dominated paths (ones that will not have positive flows due to
    # presence of proper subpaths with a lower latency
    # For full recharging, proper superset paths are dominated
    def removeDominatedPaths(self) -> List[Path]:
        print('Removing dominated paths.')
        # pathSet = {frozenset(p.edges) for p in pathList}
        # pathSet = dict.fromkeys(pathList).keys()
        nonDomPathList = []
        # while len(pathList) > 0:
            # print('len(pathList) ', len(pathList))
        for p in globPathList:
            # if any(set(p.edges) > s for s in pathSet):
                # continue
            # else:
                # nonDomPathList.append(p)
            addToPaths = True
            for s in globPathList:
                # if p.edges[0].node_from == s.edges[0].node_from:
                if frozenset(p.edges) > frozenset(s.edges):
                    # print('dominated ', self.printPathInNetwork(p), ' by ', self.printPathInNetwork(s))
                    # pathSet.remove(p)
                    addToPaths = False
                    break
                # else:
                    # nonDomPathList.append(p)
                    # print('not dom. ', self.printPathInNetwork(p), ' by ', self.printPathInNetwork(s))
            if addToPaths:
                # print('not dom. ', self.printPathInNetwork(p))
                nonDomPathList.append(p)
            # pathList.remove(p)
            # break
        print('Number of nondominated paths found: %d'%len(nonDomPathList))
        # exit(0)
        return nonDomPathList


    # Find (cyclic) energy-feasible paths in the network from source src to
    # destination dest including ones generated using concatenation of subpaths (i)
    # from src to nodes with self-loop (ii) from nodes with self-loop to dest (iii)
    # between nodes with self-loop
    # args: source (src), destination (dest), energy budget (EB), price budget (PB)
    def findPathsWithLoops(self, src, dest, EB: number=infinity, PB: number=infinity,
            numThreads: number=1, verbose: bool=False) -> List[Path]:
        global globPathList
        # (i) All energy feasible paths from src to dest
        globPathList = self.findPaths(src, dest, EB, PB, excludeSelfLoopNodes=True)
        if False: print('Number of %s-%s paths with no recharging: '%(src,dest), len(globPathList))
        # print('paths: ',*(printPathInNetwork(p, self) for p in globPathList), sep='\n')

        # (ii) All energy feasible paths involving nodes with self loops
        # Each such path will contain a node with self loop atmost once
        # Find nodes with self loops
        selfLoopNodes = [e.node_from for e in self.edges if e.node_from == e.node_to]
        pathCombs = list(combo for r in range(len(selfLoopNodes)) for combo in
                itertools.combinations(selfLoopNodes, r+1))
        allCombs = []
        for c in pathCombs:
            allCombs = [*allCombs, *list(itertools.permutations(c))]
        # for c in allCombs:
            # print('c: ', c)
            # print('i', *(i for i in (c)))
        # print('combs: ',*(list(c) for c in allCombs), sep='\n')
        # for c in allCombs:
            # print('combs:',*(i for i in list(c)), sep='\n')
            # print('combs:',list(c), sep='\n')
        for comb in allCombs:
            if False: print('\nComb: %s'%src, *(i for i in list(comb)), dest, sep='-')
            pathsComb = self.findPaths(src, comb[0], EB, PB)
            if False: print('number of subpaths from %s-%s: '%(src,comb[0]), len(pathsComb))
            k = 1
            while k < len(comb):
                # pathsComb = self.joinFeasiblePaths(pathsComb,
                        # self.findPaths(comb[k-1], comb[k], EB, PB), EB, PB)
                pathsTmp = self.findPaths(comb[k-1], comb[k], EB, PB)
                if False: print('number of subpaths from %s-%s: '%(comb[k-1],comb[k]), len(pathsTmp))
                pathsComb = self.joinFeasiblePaths(pathsComb, pathsTmp, EB, PB)

                if False: print('number of subpaths from %s--%s: '%(src,comb[k]), len(pathsComb))
                k += 1
            # globPathList = [*globPathList, *self.joinFeasiblePaths(pathsComb,\
                    # self.findPaths(comb[-1], dest, EB, PB), EB, PB)]
            pathsTmp = self.findPaths(comb[-1], dest, EB, PB)
            if False: print('number of subpaths from %s-%s: '%(comb[-1], dest), len(pathsTmp))
            pathsComb = self.joinFeasiblePaths(pathsComb, pathsTmp, EB, PB) 
            globPathList = [*globPathList, *pathsComb]

            # print('pathsComb: ', *((self.printPathInNetwork(p),p.getEnergyConsump())
                # for p in self.joinFeasiblePaths(pathsComb, self.findPaths(comb[-1],
                    # dest, EB), EB)), sep='\n')
            if False: print('paths added: ', len(globPathList))

            globPathList = self.removeDuplicatePaths(globPathList)
            print('paths found: %d'%len(globPathList))

            if numThreads == 1:
                # globPathList = self.removeDominatedPaths(globPathList)
                globPathList = self.removeDominatedPaths()
            elif numThreads == 0:
                print('Not removing dominated paths.')
            else:
                # globPathList = self.removeDominatedPathsPar(globPathList)
                globPathList = self.removeDominatedPathsPar()

        print('Number of paths found: %d'%len(globPathList))
        # print('pathsList: ', *((self.printPathInNetwork(p), p.getPrice(),
            # p.getNetEnergyConsump(), p.getFreeFlowTravelTime()) for p in globPathList),
            # sep='\n')

        # globPathList = self.removeDuplicatePaths(globPathList)

        # TODO: Ideally, we should add a path to the globPathList only if a subpath of it
        # does not exist in the globPathList. The way we are constructing paths right now
        # is compact code-wise but removing redundant paths later is more expensive
        # globPathList = self.removeDominatedPaths(globPathList)
        # exit(0)
        return globPathList
    
    # def dsPath(self, src, dest, EB: float = float('inf'), PB: float = float('inf')) -> List[Path]:
    #     globPathList=[]
    #     pq = PriorityQueue()
    #     pq.put(PrioritizedItem(0, (src, [], 0, 0)))
    #     # pq.put((0, src, [], 0, 0))  # (current_time, current_node, current_path, current_energy, current_price)

    #     shortest_paths = {src: (0, [], 0, 0)}  # cost, path_edges, energy spent, price spent

    #     while not pq.empty():
    #         store = pq.get()
    #         # print(store.priority)
    #         currTime = store.priority
    #         currNode, currPath, currEnergy, currPrice = store.item
    #         # print(currNode," ",currTime," ",pq.qsize())

    #         # If we reached the destination
    #         if currNode == dest:
    #             potential_path = Path(currPath, src)  # Creating the path
    #             # printPathInNetwork("check")
    #             if currEnergy <= EB and currPrice <= PB:
    #                 globPathList.append(potential_path)
    #                 return globPathList

    #         # Visiting each neighbor of the current node
    #         for outedge in currNode.outgoing_edges:
    #             # print(currNode," ",outedge)
    #             neighbor = outedge.node_to
    #             updatedTime = currTime + outedge.tau
    #             updatedEnergy = currEnergy + outedge.ec
    #             updatedPrice = currPrice + outedge.price
    #             updatedPath = currPath + [outedge]

    #             # print(currNode," ",neighbor," ",updatedTime," ",updatedEnergy," ",updatedPrice)
    #             # Checking if we have a new node or if we are getting a faster route
    #             # print(neighbor not in shortest_paths," ",updatedEnergy," ",EB," ",updatedPrice <= PB)
                    
    #             if ((neighbor not in shortest_paths) or (updatedTime < shortest_paths[neighbor][0])) and updatedEnergy <= EB and updatedPrice <= PB:
    #                 # print(currNode," ",neighbor," ",updatedTime," ",updatedEnergy," ",updatedPrice)
    #                 shortest_paths[neighbor] = (updatedTime, updatedPath, updatedEnergy, updatedPrice)
    #                 pq.put(PrioritizedItem(updatedTime, (neighbor, updatedPath, updatedEnergy, updatedPrice)))

        
    #     return globPathList

    def calculate_cost(self,energy: float,price: float,time: float,alpha:float,priceToTime:float)->float:
        beta=0
        #alpha=1
        #priceToTime=0
        return alpha*(time+priceToTime*price)+beta*energy

    def is_charging_station(self,node: Node) -> bool:
        return any(edge.node_to == node and edge.ec < 0 for edge in node.outgoing_edges)

    def get_charging_edge(self,node: Node) -> Edge:
        return next((edge for edge in node.outgoing_edges if edge.node_to == node and edge.ec < 0), None)

    def dsPath(self, src: Node, dest: Node, EB: float = float('inf'), PB: float = float('inf'), 
               alpha: float = 1, priceToTime: float = 0) -> List[Path]:
        globPathList = []
        pq = PriorityQueue()
        pq.put(PrioritizedItem(0, (0, src, [], 0, 0, set())))  
        
        all_paths = {src: [(0, [], 0, 0, set())]} 
        
        while not pq.empty():
            store = pq.get()
            currCost = store.priority
            currTime, currNode, currPath, currEnergy, currPrice, charged_stations = store.item
            
            if currNode == dest:
                potential_path = Path(currPath, src)
                if currEnergy <= EB and currPrice <= PB:
                    globPathList.append(potential_path)
                    return globPathList
            
            for outedge in currNode.outgoing_edges:
                neighbor = outedge.node_to
                updatedTime = currTime + outedge.tau
                updatedEnergy = currEnergy + outedge.ec
                updatedPrice = currPrice + outedge.price
                updatedPath = currPath + [outedge]
                updatedChargedStations = charged_stations.copy()
                
                if updatedEnergy <= EB and updatedPrice <= PB:
                    updatedCost = self.calculate_cost(updatedEnergy, updatedPrice, updatedTime, alpha, priceToTime)
                    pq.put(PrioritizedItem(updatedCost, (updatedTime, neighbor, updatedPath, updatedEnergy, updatedPrice, updatedChargedStations)))
                    
                    if neighbor not in all_paths:
                        all_paths[neighbor] = []
                    all_paths[neighbor].append((updatedTime, updatedPath, updatedEnergy, updatedPrice, updatedChargedStations))

                if updatedEnergy <= EB and updatedPrice <= PB and self.is_charging_station(neighbor) and neighbor not in charged_stations:
                    charging_edge = self.get_charging_edge(neighbor)
                    if charging_edge:
                        recharged_energy = 0
                        recharged_price = updatedPrice + charging_edge.price
                        recharged_time = updatedTime + charging_edge.tau
                        recharged_path = updatedPath + [charging_edge]
                        recharged_charged_stations = updatedChargedStations.copy()
                        recharged_charged_stations.add(neighbor)

                        recharged_cost = self.calculate_cost(recharged_energy, recharged_price, recharged_time, alpha, priceToTime)
                        pq.put(PrioritizedItem(recharged_cost, (recharged_time, neighbor, recharged_path, recharged_energy, recharged_price, recharged_charged_stations)))
                        
                        all_paths[neighbor].append((recharged_time, recharged_path, recharged_energy, recharged_price, recharged_charged_stations))
        
        return globPathList

    # def removeSlowEquivalentPaths(self, paths):
    #     # This method removes paths that take more energy and price but take longer time
    #     filtered_paths = []
    #     for path in paths:
    #         time, path_edges, energy, price = path
    #         if not any(other_time <= time and other_energy <= energy and other_price <= price 
    #                 for other_time, _, other_energy, other_price in paths):
    #             filtered_paths.append(path)
    #     return filtered_paths



# A directed path
class Path:
    edges: List[Edge]
    firstNode: Node

    def __init__(self, edges: List[Edge], start: Node=None):
        # For the empty path we explicitely have to specify the start node
        assert(len(edges) > 0 or start is not None)
        if start is None:
            start = edges[0].node_from
        self.firstNode = start

        currentNode = start
        self.edges = []
        for e in edges:
            # Check whether edges do form a path:
            assert(currentNode == e.node_from)
            self.edges.append(e)
            currentNode = e.node_to

    def add_edge_at_start(self,e : Edge):
        # Adds an edge at the beginning of a path
        assert(e.node_to == self.getStart())
        self.edges.insert(0,e)
        self.firstNode = e.node_from

    def add_edge_at_end(self,e : Edge):
        # Adds an edge at the end of a path
        assert (e.node_from == self.getEnd())
        self.edges.append(e)

    def getStart(self) -> Node:
        return self.firstNode

    def getEnd(self) -> Node:
        if len(self.edges) == 0:
            return self.firstNode
        else:
            return self.edges[-1].node_to

    def getFreeFlowTravelTime(self):
        # TODO: put in checks
        fftt = 0
        for e in self.edges:
            fftt += e.tau
        return fftt

    def getEnergyConsump(self):
        # TODO: put in checks
        ec = 0
        for e in self.edges:
            ec += e.ec
        return ec

    def getNetEnergyConsump(self):
        # TODO: put in checks
        ec = 0
        for e in self.edges:
            if not (e.node_from == e.node_to): ec += abs(e.ec)
        return ec

    def getPrice(self):
        totprice = 0
        for e in self.edges:
            totprice += e.price
        return totprice

    def getNodesInPath(self) -> List[Node]:
        nodeList = [self.firstNode]
        for e in self.edges:
            nodeList.append(e.node_to)
        return nodeList

    def __len__(self):
        return len(self.edges)

    def __str__(self) -> str:
        s = str(self.firstNode) + "-" + str(self.getEnd()) + " path: "
        for e in self.edges:
            s += str(e)
        return s

    def __eq__(self, other : Path) -> bool:
        if not isinstance(other,Path):
            return False
        if len(self.edges) != len(other.edges):
            return False
        for i in range(len(self.edges)):
            if self.edges[i] != other.edges[i]:
                return False
        return True

    def __hash__(self):
        # TODO: This is a very bad way of hashing
        # If possible, this should be avoided at all as Path-objects are not immutable
        h = hash(self.firstNode)
        for e in self.edges:
            h = 2*h + hash(e)
        return h

    # Function to check if a node is part of the path
    def isNodeInPath(self, x: Node) -> bool:
        # print("nodes in path ", str(self.getNodesInPath()))
        if x in self.getNodesInPath():
            # print("False")
            return True
        else:
            # print("True")
            return False


# Creates a random series parallel network with m edges
# The source node will be named s, the sink node t. All other nodes are nameless
def createRandomSPnetwork(m:int) -> Network:
    # TODO: Also create random capacities and travel times
    # In some structured way? Maybe let the user specify some point between the following two extremes:
    # - completely deterministic: Subdividing an edge subdivides the lengths, duplicating splits the capacity
    # - completely random: All edges get random capacities and travel times
    assert(m>0)

    network = Network()
    network.addNode("s")
    network.addNode("t")
    network.addEdge("s","t",1,1)
    m -= 1

    random.seed()
    for _ in range(m):
        e = random.choice(network.edges)
        if random.choice([True,False]):
            network.subdivide(e,1,1,1)
        else:
            network.duplicate(e,1,1,1)
    return network

# Print path p based on edge ids in the network G
def printPathInNetwork(p: Path, G: Network):
    s = str()
    for e in p.edges:
        if len([i for i in e.node_from.outgoing_edges\
                if i.node_to == e.node_to]) > 1:
            s += str(G.edges.index(e))
        s += str(e)
    return s

# TODO: a function to get all energy feasbile s-t paths in a network N, given
# nodes s and t
def getPathList(self, network: Network, sourceNode: Node, sinkNode: Node):
    print("To be implemented")


