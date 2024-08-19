import heapq
from collections import namedtuple, defaultdict

Edge = namedtuple('Edge', ['node_to', 'cost'])
Node = namedtuple('Node', ['name'])

class Path:
    def __init__(self, edges, source):
        self.edges = edges
        self.source = source

    def __repr__(self):
        return f"Path({[edge.node_to.name for edge in self.edges]})"

# Graph definition with a negative self-loop at E (no direct connection from E to D)
graph = {
    Node('A'): [Edge(Node('B'), 3)],
    Node('B'): [Edge(Node('C'), 2)],
    Node('C'): [Edge(Node('D'), 6), Edge(Node('E'), 1)],
    Node('E'): [Edge(Node('B'), 2), Edge(Node('E'), -20)],  # Self-loop with negative cost at E
    Node('D'): []
}

# Normal Dijkstra's Algorithm To observe infinite loop we can uncomment the following function observe how the normal dijkstra would get stuck in a loop of negative edge weight
# def normal_dijkstra(graph, src, dest, initial_fuel):
#     pq = []
#     heapq.heappush(pq, (0, src, [], initial_fuel))
#     shortest_paths = {src: (0, initial_fuel)}

#     while pq:
#         current_cost, current_node, current_path, current_fuel = heapq.heappop(pq)

#         if current_node == dest:
#             return current_path + [current_node], current_cost

#         for edge in graph[current_node]:
#             neighbor = edge.node_to
#             fuel_cost = edge.cost
#             new_cost = current_cost + fuel_cost
#             new_fuel = current_fuel - fuel_cost

#             if new_fuel < 0:
#                 continue  # If fuel goes below 0, skip this path

#             if neighbor not in shortest_paths or new_cost < shortest_paths[neighbor][0]:
#                 shortest_paths[neighbor] = (new_cost, new_fuel)
#                 heapq.heappush(pq, (new_cost, neighbor, current_path + [current_node], new_fuel))

#     return None  # If no path is found

# Modified Dijkstra's Algorithm allowing one revisit to a node
def modified_dijkstra(graph, src, dest, initial_fuel):
    pq = []
    heapq.heappush(pq, (0, src, [], initial_fuel, defaultdict(int)))
    shortest_paths = {src: (0, initial_fuel)}

    while pq:
        current_cost, current_node, current_path, current_fuel, visited_nodes = heapq.heappop(pq)

        if current_node == dest:
            return current_path + [current_node], current_cost

        for edge in graph[current_node]:
            neighbor = edge.node_to
            fuel_cost = edge.cost
            new_cost = current_cost + fuel_cost
            new_fuel = current_fuel - fuel_cost

            if new_fuel < 0:
                continue  # If fuel goes below 0, we would skip this path

            # Allow visiting a node twice at most
            if neighbor in visited_nodes:
                if visited_nodes[neighbor] >= 2:
                    continue
                visited_nodes[neighbor] += 1
            else:
                visited_nodes[neighbor] = 1

            if neighbor not in shortest_paths or new_cost < shortest_paths[neighbor][0]:
                shortest_paths[neighbor] = (new_cost, new_fuel)
                visited_nodes[neighbor] += 1
                heapq.heappush(pq, (new_cost, neighbor, current_path + [current_node], new_fuel, visited_nodes))

            visited_nodes[neighbor] -= 1
            if visited_nodes[neighbor] == 0:
                del visited_nodes[neighbor]

    return None  # If no path is found

# Testing the algorithms
src = Node('A')
dest = Node('D')
initial_fuel = 6

print("Normal Dijkstra's Algorithm:")#for normal dijkstra 
# normal_result = normal_dijkstra(graph, src, dest, initial_fuel)
# if normal_result is None:
#     print("No valid path found.")
# else:
#     normal_path, normal_cost = normal_result
#     print(f"Path: {normal_path} with cost: {normal_cost}")

print("\nModified Dijkstra's Algorithm:")
modified_result = modified_dijkstra(graph, src, dest, initial_fuel)
if modified_result is None:
    print("No valid path found.")
else:
    modified_path, modified_cost = modified_result
    print(f"Path: {modified_path} with cost: {modified_cost}")
