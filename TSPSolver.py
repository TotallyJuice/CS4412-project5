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



class TSPSolver:
	def __init__( self, gui_view ):
		self.scenario = None
		self.startTime = 0;
		self.endTime = 0;
		self.count = 0; # number solutions
		self.max = 0; # max queue size
		self.total = 0; # total states
		self.pruned = 0; # pruned states 

	def setupWithScenario( self, scenario ):
		self.scenario = scenario


	def startTimer(self):
		self.startTime = time.time()

	def stopTimer(self):
		self.endTime = time.time()

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

	def checkSolutionCost(self, route):
		solution = TSPSolution(route)
		return solution.cost

	def reset(self):
		self.startTime = 0;
		self.endTime = 0;
		self.count = 0; # number solutions
		self.max = 0; # max queue size
		self.total = 0; # total states
		self.pruned = 0; # pruned states 


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
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self.scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
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

	def greedy( self,time_allowance=60.0 ):
		self.reset()

		cities = self.scenario.getCities()

		startCity = 0
		
		self.startTimer()

		bestRoute = None
		bestRouteLen = math.inf

		while startCity < len(cities):
			unvisitedCities = copy.deepcopy(cities)
			visitedCities = []

			visitedCities.append(unvisitedCities[startCity])
			unvisitedCities.pop(startCity)

			validPath = True

			while len(unvisitedCities) > 0 and validPath:
				minlen = math.inf
				bestCity = None

				lastVisited = visitedCities[-1]

				for city in unvisitedCities:
					l = lastVisited.costTo(city)
					if l < minlen:
						minlen = l
						bestCity = city

				if(bestCity == None):
					validPath = False
				else:
					visitedCities.append(bestCity)
					unvisitedCities.remove(bestCity)
				
			startCity += 1

			routeLen = self.checkSolutionCost(visitedCities)

			if routeLen < bestRouteLen and validPath:
				bestRouteLen = routeLen
				bestRoute = visitedCities

			if(time.time() - self.startTime > time_allowance):
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
		
	def branchAndBound( self, time_allowance=60.0 ):
		self.reset()

		queue = PriorityQueue()
		cities = self.scenario.getCities()

		self.startTimer()

		bssf = self.greedy(time_allowance / 10)

		initState = PartialSolution(cities)

		queue.put(initState)

		while(not queue.empty()):
			state = queue.get()
			if(state.bound < bssf['cost']):
				if(len(state.route) == len(cities)):
					route = []
					for cityNum in state.route:
						route.append(cities[cityNum])
					bssf = self.createSolution(route)
				else:
					children = state.createSubProblems()
					
					
					for child in children:
						if(child.bound < math.inf):
							queue.put(child)
							
					


		self.stopTimer()

		return bssf





	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		self.reset()
		



