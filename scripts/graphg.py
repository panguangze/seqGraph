# Python3 program for the above approach

class Graph:

    # Function to initialize the graph
    def __init__(self, V):

        self.V = V
        self.adj = [[] for i in range(self.V)]

    # Function that adds directed edges
    # between node v with node w
    def addEdge(self, v, w):
        self.adj[v].append(w);

    # Function to perform DFS Traversal
    # and return True if current node v
    # formes cycle
    def printNodesNotInCycleUtil(self, v, visited,recStack, cyclePart):

        # If node v is unvisited
        if (visited[v] == False):

            # Mark the current node as
            # visited and part of
            # recursion stack
            visited[v] = True;
            recStack[v] = True;

            # Traverse the Adjacency
            # List of current node v
            for child in self.adj[v]:

                # If child node is unvisited
                if (not visited[child] and self.printNodesNotInCycleUtil(child, visited,recStack, cyclePart)):

                    # If child node is a part
                    # of cycle node
                    cyclePart[child] = 1;
                    return True;

                # If child node is visited
                elif (recStack[child]):
                    cyclePart[child] = 1;
                    return True;

        # Remove vertex from recursion stack
        recStack[v] = False;
        return False;

    # Function that print the nodes for
    # the given directed graph that are
    # not present in any cycle
    def printNodesNotInCycle(self):

        # Stores the visited node
        visited = [False for i in range(self.V)];

        # Stores nodes in recursion stack
        recStack = [False for i in range(self.V)];

        # Stores the nodes that are
        # part of any cycle
        cyclePart = [False for i in range(self.V)]

        # Traverse each node
        for i in range(self.V):

            # If current node is unvisited
            if (not visited[i]):

                # Perform DFS Traversal
                if(self.printNodesNotInCycleUtil(
                        i, visited, recStack,
                        cyclePart)):

                    # Mark as cycle node
                    # if it return True
                    cyclePart[i] = 1;

        # Traverse the cyclePart[]
        for i in range(self.V):

            # If node i is not a part
            # of any cycle
            if (cyclePart[i] == 0) :
                print(i,)

# Function that print the nodes for
# the given directed graph that are
# not present in any cycle
def solve( N, E, Edges):

    # Initialize the graph g
    g = Graph(N);

    # Create a directed Graph
    for i in range(E):
        # print(i, Edges[i][1])
        # print(Edges[i][0])
        # print(node_idxes_map[Edges[i][0]])
        # g.addEdge(node_idxes_map[Edges[i][0]],
        #           node_idxes_map[Edges[i][1]]);
        g.addEdge(Edges[i][0],Edges[i][1]);

    # Function Call
    # print('gggg')
    g.printNodesNotInCycle();

# Driver Code
if __name__=='__main__':

    # Given Number of nodes
    N = 6;

    # Given Edges
    E = 7;

    Edges = [ [ 0, 1 ], [ 0, 2 ],
              [ 1, 3 ], [ 2, 1 ],
              [ 2, 5 ], [ 3, 0 ],
              [ 4, 5 ] ];

    # Function Call
    solve(N, E, Edges);

# This code is contributed by rutvik_56
