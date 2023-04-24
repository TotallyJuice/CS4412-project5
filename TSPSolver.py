#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
import copy
from queue import PriorityQueue
import random


class TSPSolver:
    def __init__(self, gui_view):
        self.scenario = None
        self.startTime = 0;
        self.endTime = 0;
        self.count = 0;  # number solutions
        self.max = 0;  # max queue size
        self.total = 0;  # total states
        self.pruned = 0;  # pruned states

    def setupWithScenario(self, scenario):
        self.scenario = scenario

    def startTimer(self):
        self.startTime = time.time()

    def stopTimer(self):
        self.endTime = time.time()

    # generates a solution as well as all the statistical data
    def createSolution(self, route):

        solution = TSPSolution(route)

        results = {}
        results['time'] = self.endTime - self.startTime
        results['cost'] = solution.cost
        results['count'] = self.count
        results['soln'] = solution
        results['max'] = self.max
        results['total'] = self.total
        results['pruned'] = self.pruned
        return results

    # returns the cost of a route
    def checkSolutionCost(self, route):
        solution = TSPSolution(route)
        return solution.cost

    # resets the class variables at the start of each solve
    def reset(self):
        self.startTime = 0;
        self.endTime = 0;
        self.count = 0;  # number solutions
        self.max = 0;  # max queue size
        self.total = 0;  # total states
        self.pruned = 0;  # pruned states

    ''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self.scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

    def greedy(self, time_allowance=60.0):
        self.reset()

        cities = self.scenario.getCities()

        startCity = 0

        self.startTimer()

        bestRoute = None
        bestRouteLen = math.inf

        # check every start city
        while startCity < len(cities):
            # self.count += 1

            # create lists of cities in and out of tour
            unvisitedCities = copy.deepcopy(cities)
            visitedCities = []

            # move the starting city over
            visitedCities.append(unvisitedCities[startCity])
            unvisitedCities.pop(startCity)

            validPath = True

            # continue until all cities are visited or a dead end if found
            while len(unvisitedCities) > 0 and validPath:
                minlen = math.inf
                bestCity = None

                lastVisited = visitedCities[-1]

                # find closest unvisited city
                for city in unvisitedCities:
                    l = lastVisited.costTo(city)
                    if l < minlen:
                        minlen = l
                        bestCity = city

                if (bestCity == None):
                    validPath = False
                else:
                    visitedCities.append(bestCity)
                    unvisitedCities.remove(bestCity)

            startCity += 1

            # remember which path is shortest
            routeLen = self.checkSolutionCost(visitedCities)

            if routeLen < bestRouteLen and validPath:
                bestRouteLen = routeLen
                bestRoute = visitedCities

            # check for timeout
            if (time.time() - self.startTime > time_allowance):
                self.stopTimer()
                results = self.createSolution(bestRoute)
                return results

        self.stopTimer()

        results = self.createSolution(bestRoute)
        return results

    ''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''

    def branchAndBound(self, time_allowance=60.0):
        self.reset()

        queue = PriorityQueue()
        cities = self.scenario.getCities()

        self.startTimer()

        bssf = self.greedy(time_allowance / 10)

        # create state with all cities unvisited
        initState = PartialSolution(cities)

        queue.put(initState)

        # continue until queue is empty
        while (not queue.empty()):
            state = queue.get()
            # check if branch should be continued
            if (state.bound < bssf['cost']):
                # check if complete route
                if (len(state.route) == len(cities)):
                    # create new bssf
                    route = []
                    for cityNum in state.route:
                        route.append(cities[cityNum])
                    bssf = self.createSolution(route)
                    self.count += 1
                else:
                    # expand tree
                    children = state.createSubProblems()

                    for child in children:

                        self.total += 1
                        # check if branch is dead end
                        if (child.bound < bssf['cost']):
                            queue.put(child)
                        else:
                            self.pruned += 1
                    # record longest queue
                    if (queue.qsize() > self.max):
                        self.max = queue.qsize()

            else:
                self.pruned += 1

            # check for timeout
            if (time.time() - self.startTime > time_allowance):
                self.stopTimer()

                return self.createSolution(bssf['soln'].route)

        self.stopTimer()

        return self.createSolution(bssf['soln'].route)

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''

	#return the city that is not in the tour that is the shortest distance from a city in the tour
	def getClosest(tour, unvisitedCities):
		pass

	#return a random city from the list
	def getRandom(unvisitedCities):
		pass
    # return a random city from the list
    def getRandom(cities):

        return random.choice(cities.getRandom())

	#return the city that is the furthest distance (but not inf) from a city in the tour
	def getFurthest(tour, unvisitedCities):
		pass
		
	def fancy(self,time_allowance=60.0 ):
		self.reset()

		self.startTimer()

		#select two start cities
		#default first two cities in list
		cities = self.scenario.getCities()
		city1 = cities[0]
		city2 = cities[1]

		tour = []
		unvisitedCities = copy.deepcopy(cities)

		tour.append(city1)
		unvisitedCities.remove(city1)
		tour.append(city2)
		unvisitedCities.remove(city2)

		#while tour < cities
		while(len(tour) < len(cities)):

			#select next city
			newCity = getRandom(unvisitedCities)
			#newCity = getClosest(tour, unvisitedCities)
			#newCity = getFurthest(tour, unvisitedCities)

			bestCost = math.inf
			insertBefore = 1

			for i in range(1,len(tour)):
			#check if position is new min tour route
				testTour = tour[:i] + [newCity] + tour[i:]
				cost = checkSolutionCost(self, testTour)
				if(cost < bestCost):
					bestCost = cost
					insertBefore = i


			#insert chosen city 
			tour = tour[:insertBefore] + [newCity] + tour[insertBefore:]

		self.stopTimer()
		
		return self.createSolution(tour)
	
		



