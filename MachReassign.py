import sys
import sys
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math


class Assignment:

    def __init__(self):
        return


    def readAssignmentFile(self, filename):
        ## oldAssignment gives the initial assignment
        ## currentAssignment gives the current assignment
        ## Both are arrays in which the jth entry gives the machine to which process j is assigned
        
        with open(filename, "r") as f:
            line = f.readline()
            x = line.rstrip().split(" ")
            self.oldAssignment = list(map(int, x))
        self.currentAssignment = []
        for i in range(len(self.oldAssignment)):
            self.currentAssignment.append(self.oldAssignment[i])

    def readModelFile(self, filename):
        with open(filename, "r") as f:
            data = []
            for line in f.readlines():
                x = line.rstrip().split(" ")
                y = list(map(int, x))
                data.extend(y)
        
        ### Get Resource data
        ## numResources = number of resources
        ## isTransient = array with True if rth resource is transient, false otherwise
        ## loadCostWeight = array with weight_loadCost for each resource
        ## transientResources = list of the transient resources
        self.numResources = data.pop(0)
        self.isTransient = [False for i in range(self.numResources)]
        self.loadCostWeight = [0 for i in range(self.numResources)]
        self.transientResources = []
        for r in range(self.numResources):
            if data.pop(0) == 1:
                self.isTransient[r] = True
                self.transientResources.append(r)
            self.loadCostWeight[r] = data.pop(0)

        ### Get Machine data
        ## numMachines = number of machines
        ## location = array with ith entry the location of machine i
        ## locations = list of possible locations
        ## neighborhood = array with ith entry the neighborhood of machine i
        ## capacities = doubly-indexed array with entry [i][r] the capacity of machine i for resource r
        ## safetyCapacities = doubly-indexed array with entry [i][r] the safety capacity of machine i for resource r
        ## moveCost = doubly-indexed array with entry [i1][i2] the machine move cost from machine i1 to machine i2
        ## numNeighborhoods = highest indexed neighborhood
        ## machinesInNeighborhood = array with list of machines in neighborhood n in the nth entry

        self.numMachines = data.pop(0)
        self.location = [i for i in range(self.numMachines)]
        self.locations = []
        self.neighborhood = [0 for i in range(self.numMachines)]
        self.capacities = [[0 for i in range(self.numResources)] for j in range(self.numMachines)]
        self.safetyCapacities = [[0 for i in range(self.numResources)] for j in range(self.numMachines)]
        self.moveCost = [[0 for i in range(self.numMachines)] for j in range(self.numMachines)]
        for i in range(self.numMachines):
            self.neighborhood[i] = data.pop(0)
            self.location[i] = data.pop(0)
            if self.location[i] not in self.locations:
                self.locations.append(self.location[i])
            self.capacities[i] = data[:self.numResources]
            del data[:self.numResources]
            self.safetyCapacities[i] = data[:self.numResources]
            del data[:self.numResources]
            self.moveCost[i] = data[:self.numMachines]
            del data[:self.numMachines]

        self.numNeighborhoods = max(self.neighborhood)
        self.machinesInNeighborhood = [[i for i in range(self.numMachines) if self.neighborhood[i] == n] for n in range(self.numNeighborhoods)]


        ### Get Service data
        ## numServices = number of services
        ## spreadMin = array with SpreadMin[s] in the sth entry
        ## serviceDependencies = array with list of services in the sth entry that service s depends on

        self.numServices = data.pop(0)
        self.spreadMin = [0 for i in range(self.numServices)]
        self.serviceDependencies = [[] for i in range(self.numServices)]
        for i in range(self.numServices):
            self.spreadMin[i] = data.pop(0)
            numDependencies = data.pop(0)
            self.serviceDependencies[i] = data[:numDependencies]
            del data[:numDependencies]

        ### Get Process data
        ## numProcesses = number of processes
        ## processService = array with the service of process j in the jth entry
        ## requirements = doubly-indexed array with [j][r] entry giving the requirement of job j for resource r
        ## processMoveCost = array with cost of moving process j in the jth entry
        ## serviceProcesses = array with a list of the processes for service s in the sth entry

        self.numProcesses = data.pop(0)
        self.processService = [0 for i in range(self.numProcesses)]
        self.requirements = [[0 for i in range(self.numResources)] for j in range(self.numProcesses)]
        self.processMoveCost = [0 for i in range(self.numProcesses)]
        self.serviceProcesses = [[] for s in range(self.numServices)]
        for i in range(self.numProcesses):
            self.processService[i] = data.pop(0)
            self.serviceProcesses[self.processService[i]].append(i)
            self.requirements[i] = data[:self.numResources]
            del data[:self.numResources]
            self.processMoveCost[i] = data.pop(0)

        ### Get Balance data
        ## numBalance = number of balance constraints
        ## balanceWeight = array with the weight balance cost of the bth constraint in entry b
        ## balanceTriple = array of 3 element lists [r1, r2, target] in the bth entry for the bth balance constraint

        self.numBalance = data.pop(0)
        self.balanceWeight = [0 for i in range(self.numBalance)]
        self.balanceTriple = [[] for i in range(self.numBalance)]
        for i in range(self.numBalance):
            self.balanceTriple[i] = data[:3]
            del data[:3]
            self.balanceWeight[i] = data.pop(0)

        ### Get weights
        self.weightProcessMoveCost = data.pop(0)
        self.weightServiceMoveCost = data.pop(0)
        self.weightMachineMoveCost = data.pop(0)

        return

    def jobsAssignedToMachine(self, i, a):
        ### Returns all the processes assigned to machine i in assignment a
        return [j for j, e in enumerate(a) if e == i]

    def machineJobAssignedTo(self, j):
        ### Returns the machine job j is assigned to
        return self.currentAssignment[j]

    def resourceUsedOnMachine(self, i, r, a):
        ### Returns the amount of resource r being used on machine i by processes assigned to i in assignment a
        resourceUsage = [self.requirements[j][r] for j in self.jobsAssignedToMachine(i, a)]
        return sum(resourceUsage)
    #the below two methods are used in isFeasible() and are also given by TA
    def resourceUsedOnMachineTR(self, i, r, old, new):
        jobs_union = set(self.jobsAssignedToMachine(i, old) + self.jobsAssignedToMachine(i, new))
        resourceUsage= [self.requirements[j][r] for j in jobs_union]
        #print(self.jobsAssignedToMachine(i,a))
        return sum(resourceUsage)
        
    ### Check transient resource assignment
    def CheckTransience(self):
      for i in range(self.numMachines):
        for r in self.transientResources:
          TR_cap = self.resourceUsedOnMachineTR(i, r, self.oldAssignment, self.currentAssignment)
          if TR_cap > self.capacities[i][r]:
            return False
      return True
      
######################################PART1START####################################
######################################PART1START####################################

    #Max Wulff Barbara Sudol Sep 25 2018
    ### Check capacity constraints
    def CheckCapac(self):
      for r in range(self.numResources):
        for i in range(self.numMachines):
          if self.capacities[i][r] < self.resourceUsedOnMachine(i,r,self.currentAssignment):
            return False
      return True
        
    ### Check conflict constraints
    def CheckConflict(self):
      for s in range(self.numServices):
        for i in self.serviceProcesses[s]:
          for j in self.serviceProcesses[s]:
            if not i == j:
              if self.machineJobAssignedTo(i) == self.machineJobAssignedTo(j):
                return False
      return True
        
    ### Check spread constraints
    def CheckSpread(self):    
      for s in range(self.numServices):
        locs = []
        for i in self.serviceProcesses[s]:
          if self.location[self.machineJobAssignedTo(i)] not in locs:
            locs.append(self.location[self.machineJobAssignedTo(i)]) 
        if self.spreadMin[s] > len(locs):
          return False
      return True

    ### Check dependencies
    def CheckDepend(self):
      for s in range(self.numServices):
        for t in self.serviceDependencies[s]:
          neigh = [];
          for j in self.serviceProcesses[t]:
            a = self.machineJobAssignedTo(j)  
            neigh.append(self.neighborhood[a])
          for i in self.serviceProcesses[s]:
            a = self.machineJobAssignedTo(i)
            if self.neighborhood[a] not in neigh:
              return False
      return True

    def isFeasible(self):
      return (self.CheckCapac() and self.CheckConflict() and self.CheckSpread() and self.CheckDepend() and self.CheckTransience())
    
######################################PART1END####################################
######################################PART1END####################################

##############################PART2START############################################
##############################PART2START############################################
### Project Part 2
### Max Wulff mcw232 Barbara Sudol bas334
### Oct 25, 2018


### Calculate Load Cost

    def loadCost(self,r):
      cost = 0
      for m in range(self.numMachines):
          cost = cost + max(0 , self.resourceUsedOnMachine(m,r,self.currentAssignment) - self.safetyCapacities[m][r])
      return cost

    ### Calculate Balance Cost

    def balanceCost(self , b):
      a = self.balanceTriple[b]
      cost = 0;
      for m in range(self.numMachines):
        cost = cost + max(0, a[2]*self.A(m,a[0]) - self.A(m,a[1]))   
      return cost

    ### Calculate Process Move Cost

    def processMoveCosts(self):
      cost = 0;
      for p in range(self.numProcesses): 
        if not self.currentAssignment[p] == self.oldAssignment[p]:
          cost = cost + self.processMoveCost[p]
      return cost

    ### Calculate Service Move Cost

    def serviceMoveCost(self):
      cost = 0
      for s in range(self.numServices):
        count = 0
        for p in self.serviceProcesses[s]:
          if not self.currentAssignment[p] == self.oldAssignment[p]:
            count = count + 1
        if count > cost:
          cost = count

      return cost    

    ### Calculate Machine Move Cost

    def machineMoveCost(self):
      cost = 0
      for p in range(self.numProcesses):
        cost = cost + self.moveCost[self.oldAssignment[p]][self.currentAssignment[p]]
      return cost

    ### Calculate Total Cost

    def assignmentCost(self):
      cost = 0     
      for r in range(self.numResources):
        cost = cost + self.loadCostWeight[r]*self.loadCost(r)

      for b in range(self.numBalance):  
        cost = cost + self.balanceWeight[b]*self.balanceCost(b)
      cost = cost + self.weightProcessMoveCost*self.processMoveCosts()+ self.weightServiceMoveCost*self.serviceMoveCost()+ self.weightMachineMoveCost*self.machineMoveCost()
      return cost

    ### Swap Processes Method

    def swapProcesses(self,p1,p2):
      i1 = self.currentAssignment[p1]
      i2 = self.currentAssignment[p2]
      self.currentAssignment[p1] = i2
      self.currentAssignment[p2] = i1

    ### Move Process Method

    def moveProcess(self,p,m):
      self.currentAssignment[p] = m
      
    def A(self,m,r):
      return self.capacities[m][r] - self.resourceUsedOnMachine(m,r,self.currentAssignment)
    
    def outputAssignment(self, assignment, filename):
      F = open(filename, "a")
      for a in range(len(self.currentAssignment)):
        F.write(str(assignment[a]))
      F.close() 

###################################PART2END#########################
###################################PART2END#########################

###################################BEGINNINGPART3###################
###################################BEGINNINGPART3###################
#November 12,2018

    def localSearch(self):
      for i in range(self.numProcesses):
        for j in range(i, self.numProcesses):
          if not self.machineJobAssignedTo(i) == self.machineJobAssignedTo(j):
            first = self.assignmentCost()
            self.swapProcesses(i,j)
            second = self.assignmentCost()
            if first < second: #if this didn't improve the cost, switch back
              self.swapProcesses(i,j)
            if not self.isFeasible():
              self.swapProcesses(i,j)
      return self.currentAssignment
    
    def uphillLocalSearch(self):
      #simulated annealing using our above localSearch method:
      T = 1
      To = 0.1
      alpha = 0.95
      k = self.numProcesses
      c = 0
      cprime = 0
      while (T > To):
        for i in range(k):
          a, b = self.findRandomFeasibleSwap()
          c = self.assignmentCost()
          self.swapProcesses(a,b)
          cprime = self.assignmentCost()
          delta = cprime - c
          if(delta > 0):
            r = np.random.rand()
            if(r >= math.exp(-delta/T)):
              self.swapProcesses(a,b)
        T = T*alpha
      return self.currentAssignment
    
    def findRandomFeasibleSwap(self):
      #finds a random process a and b where swapping the two processes creates an assignment that is still feasible
      isFeas = False
      count = 0
      while (isFeas == False):
        a = np.random.randint(0,self.numProcesses)
        b = np.random.randint(0,self.numProcesses)
        self.swapProcesses(a,b)
        isFeas = self.isFeasible()
        self.swapProcesses(a,b)
        if(count >= self.numProcesses*self.numProcesses):
          return(a,a)
        count += 1
      return(a,b)
    
###################################PART3END#########################
###################################PART3END#########################

def main(argv):
    #a = Assignment()
    #a.readModelFile("model_a1_1.txt")
    #a.readAssignmentFile("assignment_a1_1.txt")
    #print(a.isFeasible())
    
    #b = Assignment()
    #b.readModelFile("model_a1_2.txt")
    #b.readAssignmentFile("assignment_a1_2.txt")
    #print(b.isFeasible())
    
    #c = Assignment()
    #c.readModelFile("model_a1_3.txt")
    #c.readAssignmentFile("assignment_a1_3.txt")
    #print(c.isFeasible())
    
    #d = Assignment()
    #d.readModelFile("model_a1_4.txt")
    #d.readAssignmentFile("assignment_a1_4.txt")
    #print(d.isFeasible())
    
    #e = Assignment()
    #e.readModelFile("model1.txt")
    #e.readAssignmentFile("solution1.txt")
    #print(e.isFeasible())
    
    #f = Assignment()
    #f.readModelFile("model2.txt")
    #f.readAssignmentFile("solution2.txt")
    #print(f.isFeasible())
    
    g = Assignment()
    g.readModelFile("model_a1_1.txt")
    g.readAssignmentFile("assignment_a1_1.txt")
    print(g.assignmentCost())
    g.localSearch()
    print(g.assignmentCost())
    g.uphillLocalSearch()
    print(g.assignmentCost())
    
    h = Assignment()
    h.readModelFile("model.txt")
    h.readAssignmentFile("solution1.txt")
    print(h.assignmentCost())
    h.localSearch()
    print(h.assignmentCost())
    h.uphillLocalSearch()
    print(h.assignmentCost())
    h.outputAssignment(h.currentAssignment,"poop.txt")
    
    i = Assignment()
    i.readModelFile("model.txt")
    i.readAssignmentFile("solution2.txt")
    print(i.assignmentCost())
    i.localSearch()
    print(i.assignmentCost())
    i.uphillLocalSearch()
    print(i.assignmentCost())

if __name__ == "__main__":
    main(sys.argv[1:])
    