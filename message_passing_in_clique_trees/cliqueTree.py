from numpy import *
from factor import *

class CliqueTree:
    cliques = []
    edges = []
    neighbors = {}
    factors = []

    potentials = []
    messages = {}
    beliefs = []

    def __init__(self, clique_list, edge_list, factors):
        self.cliques = clique_list
        self.edges = edge_list
        for (i,j) in edge_list:
            if i in self.neighbors:
                self.neighbors[i].append(j)
            else:
                self.neighbors[i] = [j]
            if j in self.neighbors:
                self.neighbors[j].append(i)
            else:
                self.neighbors[j] = [i]
        # TODO check the runnning intersection property
        if len(factors) != len(clique_list):
            raise Exception("Factors must be a list of a same length as cliques.")
        self.factors = factors
        # TODO check the family presevation property

    def beliefPropagation(self):
        message_queue = self.edges
        message_queue.extend([[j,i] for (i,j) in message_queue])
        self.calculatePotentials()
        while len(message_queue) > 0:
            p = self.getNextTransmitPair(message_queue)
            m = self.calculateMessage(p)
            print "PASSING MESSAGE: " + str(p[0]) + " -> " + str(p[1])
            print m
            print
            self.messages[tuple(p)] = m
        self.calculateBeliefs()
        

    def getNextTransmitPair(self, message_queue):
        for i in range(len(message_queue)):
            if self.canTransmit(message_queue[i]):
                return message_queue.pop(i)

    def canTransmit(self, from_to):
        (i, j) = from_to
        neigh = self.neighbors[i]
        if len(neigh) == 1:
            return True
        return all([((x,i) in self.messages) 
            for x in neigh if x != j])

    def calculateMessage(self, from_to):
        i, j = from_to
        psi = self.potentials[i]
        messages_to_i = [self.messages[(x,i)] 
            for x in self.neighbors[i] if x != j ]
        #print (i,j)
        #print "MESS: " + str(messages_to_i)
        #print "PSI: " + str(psi)
        for m in messages_to_i:
            psi *= m
        elimVars = setdiff1d(self.cliques[i], self.cliques[j])
        tau = psi.margin(elimVars)
        return tau

    def calculatePotentials(self):
        self.potentials = [prod(i) for i in self.factors]

    def calculateBeliefs(self):
        self.beliefs = [self.calculateBelief(i) 
            for i in range(len(self.cliques))]

    def calculateBelief(self, i):
        messages = [self.messages[(x,i)] 
            for x in self.neighbors[i]]
        beta = self.potentials[i] * prod(messages)
        return beta
