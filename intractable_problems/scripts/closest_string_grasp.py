import argparse
import math
import numpy as np

class ClosestStringGRASP:
    '''
    GRASP approach to solve closest string problem.
    
    Args:
        strings (list of str): List of n strings to search for closest string
        max_iterations (int): Max iterations for running a constructive heuristics
            and a local search
        alpha (float): ratio of elements that will be candidates at generating
            a randomized constructive heuristics
    '''
    def __init__(self, strings, max_iterations=100, alpha=0.1, r=2):
        self.strings = strings
        self.n = len(self.strings)
        self.m = len(self.strings[0])
        for string in self.strings:
            assert self.m == len(string)
        self.max_iterations = max_iterations
        self.alpha = alpha
        self.r = r
        
        # Initial solution
        self.solution = 'a'*self.m
        
        self.candidates()
        self.GRASP()
        self.report()
        
    def candidates(self):
        '''Generate the candidates for each position.'''
        self.S = []
        for i in range(self.m):
            S_i = []
            for j in range(self.n):
                S_i.append(self.strings[j][i])
            S_i = sorted(set(S_i), key=S_i.count, reverse=True)
            self.S.append(S_i)
        
    def GRASP(self):
        '''GRASP main function.'''
        
        for k in range(self.max_iterations):
            solution = self.random_constructive()
            solution = self.local_search(solution)
            self.update_solution(solution)
            
    def random_constructive(self):
        '''Randomized constructive heuristics.
        
        Considering the most common characters for position i,
        it selects the alpha-best candidates to be sorted.'''
        solution = ''
        for i in range(self.m):
            candidates = self.S[i]
            # candidates = sorted(set(candidates), key=candidates.count, reverse=True)
            n = math.ceil(len(candidates)*self.alpha)
            # Restricted candidates list
            RCL = candidates[:n]
            lucky = np.random.choice(RCL)
            solution += lucky
        return solution
    
    def local_search(self, solution):
        '''Local search to solution.
        
        Args:
            solution (str): Current solution
        '''
        indices = np.random.choice(list(range(self.m)), size=self.r, replace=False)
        solution, cost = self.best_solution(solution, indices)
        return solution
    
    def best_solution(self, s, indices):
        if len(indices) == 0:
            return s, self.cost(s)
        index = indices[0]
        indices = indices[1:]
        cost = np.inf
        for char in self.S[index]:
            temp_solution, temp_cost = self.best_solution(s[:index] + char + s[index+1:], indices)
            if temp_cost < cost:
                solution = temp_solution
                cost = temp_cost
        return solution, cost 
            
    def cost(self, s):
        '''Calculate cost of string s, the maximum Hamming distance to one
        of the initial strings.

        Args:
            s (str): String to have cost calculated

        Returns:
            cost (int): Cost of the string
        '''
        ds = [self.hamming(s, s_i) for s_i in self.strings]
        return max(ds)
    
    def hamming(self, s_1, s_2):
        '''Hamming distance between two strings.'''
        return sum([c_1 != c_2 for c_1, c_2 in zip(s_1, s_2)])
    
    def update_solution(self, solution):
        '''Update solution if it is a better one.'''
        if self.cost(solution) < self.cost(self.solution):
            self.solution = solution
    
    def report(self):
        print('Solution is {}.'.format(self.solution))
        print('The cost is {}.'.format(self.cost(self.solution)))

def main():
    parser = argparse.ArgumentParser(description='Closest string problem using GRASP.')
    parser.add_argument('--strings', nargs='+', help='Strings for the problem.')
    args = parser.parse_args()
    ClosestStringGRASP(args.strings)

if __name__ == '__main__':
    main()
