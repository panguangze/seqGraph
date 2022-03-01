#  Python 3 Program
#  Find nodes which are not part of any cycle in a directed graph

#  Define edge of graph node
class AjlistNode :
    #  Vertices id
    def __init__(self, id) :
        self.id = id
        self.next = None


#  Define node of graph
class GraphNode :
    #  node key value
    def __init__(self, data) :
        self.data = data
        self.edge = None


#  Define structure of graph
class MyGraph :
    #  Number of graph nodes
    def __init__(self, size) :
        if (size < 0) :
            print("\n Invalid size of nodes ", size ," ")
        else :
            self.size = size
            self.node = [None] * (self.size)
            i = 0
            #  Set graph node level and e
            while (i < self.size) :
                self.node[i] = GraphNode(i)
                i += 1



    # This function are connect nodes with one edge
    def connectEdge(self, a, b) :
        if (self.size <= 0) :
            print("\n Graph have no nodes ")
            return

        #  Create Adjacency node
        newEdge = AjlistNode(b)
        if (newEdge != None) :
            #  Get first edge of [i] node
            edge = self.node[a].edge
            if (edge == None) :
                #  Whe no edge of graph node [i] then
                #  Add a first edge of a [i] node
                self.node[a].edge = newEdge
            else :
                # Find last location to add new edge
                while (edge.next != None) :
                    edge = edge.next

                #  Add edge at last position
                edge.next = newEdge

        else :
            print("\nMemory overflow, when connect a new edge from ( ", a ," to ", b ," ) nodes. ")


    #  This function is add new edge
    def addEdge(self, a, b) :
        if (self.size > 0 and a < self.size and b < self.size) :
            #  connect edge
            #  a---->b
            self.connectEdge(a, b)
        else :
            # not valid Vertices
            print("Invalid Node Vertices ", a ," ", b, end = "")


    # Display Adjacency list of vertex
    def printGraph(self) :
        if (self.size > 0) :
            temp = None
            #  Loop controlling variable
            i = 0
            i = 0
            while (i < self.size) :
                print("\n Adjacency list of vertex ", i ," :", end = "")
                #  Get first edge of [i] node
                temp = self.node[i].edge
                while (temp != None) :
                    # temp->id is graph node vertices
                    # in this case temp->id is same as
                    # node[temp->id].data
                    print("  ", self.node[temp.id].data, end = "")
                    temp = temp.next

                i += 1

        else :
            print("Empty Graph Node", end = "")


    def resetDefault(self, visit, n) :
        i = 0
        while (i < n) :
            visit[i] = False
            i += 1


    #  This function are detect loop in given source code
    def isCyclic(self, visit, source, current) :
        if (source == current and visit[current] == True) :
            #  When loop exist
            return True
        else :
            if (visit[current] == True) :
                #  When intermediate node contains loop
                #  back to previous process
                return False


        temp = self.node[current].edge
        #  Active the visiting node status
        visit[current] = True
        #  iterating the nodes edges
        while (temp != None) :
            if (self.isCyclic(visit, source, temp.id) == True) :
                #  When source node contains cycle
                return True

            #  Visit to next edge
            temp = temp.next

        #  Deactivate the current visiting node status
        visit[current] = False
        return False

    def nonCyclicNode(self) :
        if (self.size > 0) :
            visit = [False] * (self.size)
            cycle = [False] * (self.size)
            #  Set that initial have no cyclic node
            self.resetDefault(cycle, self.size)
            i = 0
            while (i < self.size) :
                if (cycle[i] == False) :
                    self.resetDefault(visit, self.size)
                    #  Check that cycle
                    if (self.isCyclic(visit, i, i) == True) :
                        #  When i contain cycle
                        #  Then collect the cycle node
                        j = 0
                        while (j < self.size) :
                            if (visit[j] == True) :
                                #  contain cycle
                                cycle[j] = True

                            j += 1



                i += 1

            #  result counter
            count = 0
            print("\n Result : ", end = "")
            #  Print resultant node
            j = 0
            while (j < self.size) :
                if (cycle[j] == False) :
                    #  Count the number of resultant nodes
                    count += 1
                    print(" ", j, end = "")

                j += 1

            if (count == 0) :
                #  In case no result
                print("\n None ")




def main() :
    #  8 is number of nodes
    # graph = MyGraph(8)
    # # Connected two node with Edges
    # graph.addEdge(0, 1)
    # graph.addEdge(1, 2)
    # graph.addEdge(2, 4)
    # graph.addEdge(2, 7)
    # graph.addEdge(3, 2)
    # graph.addEdge(3, 4)
    # graph.addEdge(3, 5)
    # graph.addEdge(4, 2)
    # graph.addEdge(6, 2)
    # graph.addEdge(6, 4)
    # graph.addEdge(7, 0)
    g = MyGraph(5)
    g.addEdge(0,1)
    g.addEdge(1,2)
    g.addEdge(2,3)
    g.addEdge(3,1)
    g.addEdge(2,4)
    g.addEdge(2,0)
    g.printGraph()
    g.nonCyclicNode()

if __name__ == "__main__": main()
