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
      best = math.inf
      newbest = best - 1
      while(newbest < best):
        best = self.assignmentCost()
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
        newbest = self.assignmentCost()
                
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
#################################PART4BEGIN#########################
#################################PART4BEGIN#########################


    def infeasibleSearch(self):
      self.localSearch() #Start with a Local Search to Find the Locally Optimal Solution
      tried = [] #Will Keep a List of Tried Swaps
      best = self.assignmentCost() #Current Cost of Best Feasible Solution
      real = self.currentAssignment #Current Best Feasible Solution
      for i in range(15):
        
        a = self.findBestSwap(tried) #Find Best Infeasible Swap of What We Haven't Tried
        
        if a == [-1,-1]: #If No Swap Gives Better Solution: Break
          break
          
        tried.append(a) # Add this Swap to the List of Tried Swaps 
        
        self.swapProcesses(a[0],a[1]) # Swap the Processes for the Test
        
        count = 0 #Initialize a Counter
        b = self.showInfeasible() # Show Where the Infeasibility Is
        reall = real
        bestt = best
        
        while b[0] > 0 and count < 20: # While This is Still Infeasible and We Haven't Tried More than n Times
          
          count = count + 1
          something = False
          if b[0] == 1: #If Infeasible Based on Capacity
            for j in range(self.numProcesses): 
              if self.machineJobAssignedTo(j) == b[1]: #If the Process is from the Infeasible Machine
                for m in range(self.numMachines): # Try Moving Process to Every other Machine
                  
                  self.moveProcess(j,m)
                  
                  l = True
                  for r in range(self.numResources): #Loop to check if we fixed the problem
                    if self.capacities[b[1]][r] < self.resourceUsedOnMachine(i,r,self.currentAssignment):
                      l = False
                      
                  boo = True
                  d = self.CheckSpread()
                  e = self.CheckDepend()
                  f = self.CheckTransience()
                  if not d[0]:
                    boo = False
                  if not e[0]:
                    boo = False
                  if not f[0]:
                    boo = False
                    
                  if self.singleFeasibility(j,m) and self.assignmentCost() < best and l and boo: #If the single swap is feasible for capacity and conflict and we lower the assignment cost and our swap fixes the original problem and the whole assignment is still feasible for the rest of the conditions, store the swap
                    reall = self.currentAssignment #Update Best Feasible Solution and cost
                    bestt = self.assignmentCost()
                    something = True   #mark that something happened
                  self.moveProcess(j,b[1])
                  
          else: #Infeasibility based on conflict: Try moving each of the conflicting process to any of the other machines
            for m in [x for x in range(self.numMachines) if x != self.machineJobAssignedTo(b[1])]: #Loop over every machine except one currently assigned to
              self.moveProcess(b[1],m)
                
              boo = True
              d = self.CheckSpread()
              e = self.CheckDepend()
              f = self.CheckTransience()
              if not d[0]:
                boo = False
              if not e[0]:
                boo = False
              if not f[0]:
                boo = False  
            
              if self.singleFeasibility(b[1],m) and self.assignmentCost() < best and boo: #If the single swap is feasible for capacity and conflict and we lower the assignment cost and our swap fixes the original problem and the whole assignment is still feasible for the rest of the conditions, store the swap
                reall = self.currentAssignment #Update Best Solution and cost
                bestt = self.assignmentCost()
                something = True   #mark that something happened
              self.moveProcess(b[1],self.machineJobAssignedTo(b[2]))
                               
              ##
              
              self.moveProcess(b[2],m)
                
              boo = True
              d = self.CheckSpread()
              e = self.CheckDepend()
              f = self.CheckTransience()
              if not d[0]:
                boo = False
              if not e[0]:
                boo = False
              if not f[0]:
                boo = False  
            
              if self.singleFeasibility(b[2],m) and self.assignmentCost() < best and boo: #If the single swap is feasible for capacity and conflict and we lower the assignment cost and our swap fixes the original problem and the whole assignment is still feasible for the rest of the conditions, store the swap
                reall = self.currentAssignment #Update Best Solution and cost
                bestt = self.assignmentCost()
                something = True   #mark that something happened
              self.moveProcess(b[2],self.machineJobAssignedTo(b[1]))
                               
                          
              
          if not something: #If nothing happened that mean that the current problem being worked on can't be fixed and therefore we shouldn't continue
            count = 20
          
          b = self.showInfeasible() #update with the next infeasibility
        ##
                               
        self.currentAssignment = reall
                                                     
        if self.isFeasible():
          real = reall
          best = bestt 
        else:
          self.currentAssignment = real                     
                               
                               
                               
                               
    def singleFeasibility(self,j,m): #Check to see if a single process and single machine is clear for capacity and conflict
      
      for r in range(self.numResources): #Check to see if machine moved to still follows capacity constraints
        if self.capacities[m][r] < self.resourceUsedOnMachine(m,r,self.currentAssignment):
          return False
      
      for s in range(self.numServices): #Check to make sure the process doesn't violate conflict 
        if j in self.serviceProcesses[s]:
          for i in self.serviceProcesses[s]:
            if not i == j:
              if self.machineJobAssignedTo(i) == self.machineJobAssignedTo(j):
                return [False,i,j]
      
      
      
      
      
    def findBestSwap(self,dont): #find best swap period. If none exists return [-1,-1]
      best = math.inf
      a = self.currentAssignment
      ii = -1
      jj = -1 
      for i in range(self.numProcesses):
        for j in range(i,self.numProcesses):
          if [i,j] not in dont:
            self.swapProcesses(i,j)
            if (self.assignmentCost() < best and self.showInfeasible()[0] <= 2):
                a = self.currentAssignment
                ii = i
                jj = j
            self.swapProcesses(i,j)    
      return [ii,jj]      

    def computeLowerBound(self):
 
            m = Model("lowerBound")
 
            x = m.addVars(self.numMachines, self.numProcesses, lb=0.0, ub=1.0, obj=0.0, vtype=GRB.CONTINUOUS, name = "x")

            m.addConstrs((x.sum('*', j) == 1.0 for j in range(self.numProcesses)), "assign")

            m.addConstrs((quicksum(self.requirements[j][r] * x[i,j] 
                         for j in range(self.numProcesses)) 
                         <= self.capacities[i][r]
                             for i in range(self.numMachines) 
                             for r in range(self.numResources)), "capacity")
          
            m.addConstrs((quicksum(x[i,j]
                        for j in self.serviceProcesses[s])
                        <= 1
                          for s in range(self.numServices)
                          for i in range(self.numMachines)), "conflict")
          
            #balance constraint was written in collaboration with Zoe, balance constraint gets more of our sample instances to run 
            #ASK IF THIS IS OK to credit Zoe
            balanceVars = m.addVars(self.numBalance, lb=0.0, obj=0.0, vtype=GRB.CONTINUOUS, name="b")           
            #add into objective
            for i in range(self.numBalance): balanceVars[i].obj = self.balanceWeight[i]           
            #add balance constraints
            m.addConstrs(quicksum(-(self.capacities[i][self.balanceTriple[balance][1]] - quicksum(self.requirements[j][self.balanceTriple[balance][1]] * x[i, j]
                    for j in range(self.numProcesses))) + self.balanceTriple[balance][2]*(self.capacities[i][self.balanceTriple[balance][0]]- quicksum(self.requirements[k][self.balanceTriple[balance][0]] * x[i, k]
				    for k in range(self.numProcesses)) <= balanceVars[balance])
                        for balance in range(self.numBalance))
				            for i in range(self.numMachines))
            
            
          
          #y = m.addVars(self.numServices, len(self.Locations), lb=0.0, ub =1.0, obj = 0, vtype=GRB.CONTINUOUS, name = "y")
          
          #m.addConstrs((quicksum(y[s,l] 
          #             for l in range(len(self.locations)))
          #             >= self.spreadMin[s]
          #                 for s in range(self.numServices)), "spread constr")
          
          #m.addConstrs(quicksum(quicksum
          
          #make loadcost var and constraints for it
            loadCost = m.addVars(self.numMachines, self.numResources, lb = 0.0, obj = 0.0, vtype = GRB.CONTINUOUS, name = "a")

            for i in range(self.numMachines):
                for r in range(self.numResources):
                    loadCost[i,r].obj = self.loadCostWeight[r]

            m.addConstrs(loadCost[i,r] >= -self.safetyCapacities[i][r] + quicksum(self.requirements[j][r] * x[i,j] 
                         for j in range(self.numProcesses))
                             for i in range(self.numMachines) 
                             for r in range(self.numResources))    

            m.optimize()
            return m.objVal
   
          
###############################PART4END################################### 
###############################PART4END################################### 

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
    
