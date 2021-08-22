import argparse
import math

class ClosestStringTree:
    '''
    ClosestStringTree solves the closest string problem
    in a tree approach, using branch and bound.

    Args:
        strings (list of str): List of strings for the problem.

    Attributes:
        S (list of set of str): Contain the possible characters for each position
        costs (list of int): The cost for each position in relation to the
            most common character
    '''
    def __init__(self, strings):
        self.strings = strings.copy()
        self.n = len(self.strings)
        self.m = len(self.strings[0])
        for string in self.strings:
            assert len(string) == self.m

        self.characters()
        
        # Initial solution
        self.s = 'a'*self.m
        self.d = self.cost(self.s)      
        
        self.iterations = 0
        self.branch()
        self.report()

    def characters(self):
        '''Extract the characters and also calculate the cost of each position,
        in relation to the most common character for that position.'''
        self.S = []
        self.costs = []
        for i in range(self.m):
            S_i = []
            for j in range(self.n):
                S_i.append(self.strings[j][i])
            most_common = max(set(S_i), key=S_i.count)
            cost = sum([most_common != c for c in S_i])
            self.costs.append(cost)
            self.S.append(set(S_i))      
    
    def branch(self, s=''):
        '''Branch the tree.
        This is the main function of the algorithm.

        Args:
            s (str): The string currently being handled
        '''
        self.iterations += 1
        i = len(s)
        if i == self.m:
            d = self.cost(s)
            if d < self.d:
                self.d = d
                self.s = s
            return
        if self.best_cost(s) > self.d:
            # Bound
            return
        S_i = self.S[i]
        for char in S_i:
            self.branch(s + char)
            
    def best_cost(self, s):
        '''Calculate a lower bound for the optimal cost of a subtree.
        The subtree is well represented by the current substring.

        Args:
            s (str): Current substring representing subtree

        Returns:
            cost (int): A lower bound for the optimal cost of the subtree
        '''
        i = len(s)
        actual_cost = self.cost(s, [string[:i] for string in self.strings])
        best_after_cost = math.ceil(sum(self.costs[i:])/(self.m - i))
        return actual_cost + best_after_cost
            
    def cost(self, s, strings=None):
        '''Calculate cost of string s, the maximum Hamming distance to one
        of the initial strings.

        Args:
            s (str): String to have cost calculated
            strings (list of str): Optional list of strings to base the calculation

        Returns:
            cost (int): Cost of the string
        '''
        if strings is None:
            strings = self.strings
        ds = [self.hamming(s, s_i) for s_i in strings]
        return max(ds)
    
    def hamming(self, s_1, s_2):
        '''Hamming distance between two strings.'''
        return sum([c_1 != c_2 for c_1, c_2 in zip(s_1, s_2)])
    
    def report(self):
        '''Report the results.'''
        print('The solution is {}.'.format(self.s))
        print('The cost is {}.'.format(self.d))
        print('The number of iterations was {}.'.format(self.iterations))

def main():
    parser = argparse.ArgumentParser(description='Closest string problem using branch and bound.')
    parser.add_argument('--strings', nargs='+', help='Strings for the problem.')
    args = parser.parse_args()
    ClosestStringTree(args.strings)

if __name__ == '__main__':
    main()
