import math

class ClosestStringTree:
    def __init__(self, strings):
        self.strings = strings.copy()
        self.n = len(self.strings)
        self.m = len(self.strings[0])
        for string in self.strings:
            assert len(string) == self.m

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
        
        # Initial solution
        self.s = 'a'*self.m
        self.d = self.cost(self.s)      
        
        self.iterations = 0
        self.branch()
        self.report()
    
    def branch(self, s=''):
        self.iterations += 1
        i = len(s)
        if i == self.m:
            d = self.cost(s)
            if d < self.d:
                self.d = d
                self.s = s
            return
        if self.best_cost(s) > self.d:
#             print('here')
            # Bound
            return
        S_i = self.S[i]
        for char in S_i:
            self.branch(s + char)
            
    def best_cost(self, s):
        start = len(s)
        actual_cost = self.cost(s, [string[:start] for string in self.strings])
        best_after_cost = math.ceil(sum(self.costs[start:])/(self.m - start))
        return actual_cost + best_after_cost
            
    def cost(self, s, strings=None):
        if strings is None:
            ds = [self.hamming(s, s_i) for s_i in self.strings]
        else:
            ds = [self.hamming(s, s_i) for s_i in strings]
        return max(ds)
    
    def hamming(self, s_1, s_2):
        return sum([c_1 != c_2 for c_1, c_2 in zip(s_1, s_2)])
    
    def report(self):
        print('The solution is {}.'.format(self.s))
        print('The cost is {}.'.format(self.d))
        print('The number of iterations was {}.'.format(self.iterations))
        