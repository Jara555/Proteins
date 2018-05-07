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
        self.dist = 0
        if parent:
            self.path = parent.path[:]
            self.path.append(value)
            self.start = parent.start
            self.goal = parent.goal
        else:
            self.path = [value]
            self.start = start
            self.goal = goal

    def getStability(self):

    def GetDist(self):
        pass

    def CreateChildren(self):
        pass


class StateProtein(State):
    def __init__(self, value, parent, start=0, goal=0):

        super(StateProtein, self).__init__(value, parent, start, goal)
        self.dist = self.GetDist()

    def getStability(self):

        stability = len(self.value) / 2
        return stability

    def GetDist(self):
        ''' Get the current distance to goal, which will be state's value '''

        # get stability
        stability = self.getStability()

        # if goal is reached
        if self.value == self.goal:
            return 0

        # calculate distance
        dist = stability - self.goal
        return dist

    def CreateChildren(self):
        ''' Generate all possible children, by switching two adjacent letters around '''

        if not self.children:



            for i in range(len(self.goal) - 1):
                val = self.value
                val = val[:i] + val[i + 1] + val[i] + val[i + 2:]
                child = State_String(val, self)
                self.children.append(child)