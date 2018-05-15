'''

    node = patroon met 1 vouws opties
    parent = vorige vouw state
    children = alle mogelijke nodes (alle mogelijke 1 staps vouw opties)
    queue = meest promesing patronen staan bovenaan
    start value = stabiliteit van 0
    eind value = stabiliteit van -9
    value = huidige vouwpatroon
    cost = kosten van een vouw (bijv 0.1)
    distance = stabiliteit - eind value
    path = alle vouwpatronen in 1 vouws stapjes die die is afgelopen

'''


from queue import PriorityQueue


class State(object):
    ''' Base template class to be inherited from for subclasses '''

    def __init__(self, value, parent,
                 start=0, goal=0):
        ''' Constructor that stores parent, value, start and goal '''
        self.children = []
        self.parent = parent
        self.value = value
        self.stability = 0
        self.dist = 0
        self.orientations = ['+X', '-X', '+Y', '-Y']
        if parent:
            self.path = parent.path[:]
            self.path.append(value[:])
            self.start = parent.start
            self.goal = parent.goal
        else:
            self.path = [value[:]]
            self.start = start
            self.goal = goal

    def GetDist(self):
        pass

    def CreateChildren(self):
        pass


class StateProtein(State):
    def __init__(self, value, parent, start=0, goal=0):

        super(StateProtein, self).__init__(value, parent, start, goal)
        self.dist = self.GetDist()

    def GetDist(self):
        ''' Get the current distance to goal, which will be state's value '''

        self.stability = len(self.value) / 2

        # if goal is reached
        if self.stability == self.goal:
            return 0

        # calculate distance
        dist = self.stability - self.goal
        return dist

    def CreateChildren(self):
        ''' Generate all possible children,
        by getting all 1-step folding options of parent '''

        if not self.children:


            for i in range(len(self.value) - 2):
                for orientation in self.orientations:
                    val = self.value
                    val[i + 2] = orientation
                    child = StateProtein(val, self)
                    self.children.append(child)


class AStar_Solver:
    def __init__(self, start, goal):
        ''' Store the start and goal of program, and set up vars '''
        self.path = []
        self.visitedQueue = []
        self.priorityQueue = PriorityQueue()
        self.start = start
        self.goal = goal

    def Solve(self):
        ''' Create start state, then organize children
        based on value and check children of highest value child '''
        startState = StateProtein(['0', '+Y', '+Y', '+Y'],
                                  0,
                                  self.start,
                                  self.goal)
        count = 0
        self.priorityQueue.put((0, count, startState))
        while (not self.path and self.priorityQueue.qsize()):
            closestChild = self.priorityQueue.get()[2]
            closestChild.CreateChildren()
            self.visitedQueue.append(closestChild.value)
            for child in closestChild.children:

                #TODO: blijft hier hangen: maakt geen self.path aan
                if child.value not in self.visitedQueue:
                    count += 1
                    if not child.dist:
                        self.path = child.path
                        break
                    self.priorityQueue.put((child.dist, count, child))

        if not self.path:
            print("Goal of " + str(self.goal) + " is not possible!")
        return self.path

##=======================================
## MAIN

if __name__ == "__main__":

    start1 = 0
    goal1 = 0

    print('starting...')

    a = AStar_Solver(start1, goal1)
    a.Solve()
    for i in range(len(a.path)):
        print("%d) " % i + a.path[i])