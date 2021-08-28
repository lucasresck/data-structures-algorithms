from queue import PriorityQueue
import numpy as np

class BoundedPriorityQueue(PriorityQueue):
    '''Bounded priority queue (BPQ) implementation.
    
    Bounded priority queue inherited from a priority queue implementation. The
    differences are:
    - when the queue is full and the new element can be has enough priority to
    be inserted, we remove the less priority element;
    - the priority is reverted, that is, we remove the element with the 
    greatest priority, instead of the one with the least priority. This is
    because, in our application, small number means high priority.

    Args:
        k: size of the priority queue.

    Attributes:
        max_priority: maximum priority.
    '''
    def __init__(self, k):
        super().__init__(k)
        self.max_priority = 0

    def put_item(self, node,  priority):
        '''Put the item into the bounded priority queue.
        
        Args:
            node: the item to be added.
            priority: priority of the node.
        '''
        if not self.full():
            self.put((-priority, node))
            if priority > self.max_priority:
                self.max_priority = priority
        else:
            if priority < self.max_priority:
                self.get()
                self.put((-priority, node))

class KDTree:
    '''k-d tree implementation.'''
    def __init__(self, k):
        self.tree = []
        self.k = k

    def add(self, node, dim=0, subtree=None):
        '''Add a node to the tree.
        
        Args:
            node: node to be added.
            dim: which dim to be compared. By default, it starts with dim 0 at
                root node and increases by one for each level of the tree.
            subtree: optional; in which subtree to add the node. By standard,
                it is the tree itself.'''
        if subtree is None:
            subtree = self.tree
        if not len(subtree):
            subtree.extend(self.create_node(node))
        else:
            if node[dim] < subtree[1][dim]:
                self.add(node, (dim+1) % self.k, subtree[0])
            else:
                self.add(node, (dim+1) % self.k, subtree[2])

    def create_node(self, node):
        return [[], node, []]

if __name__ == '__main__':
    tree = KDTree(2)
    tree.add([3, 6])
    tree.add([17, 15])
    tree.add([13, 15])
    tree.add([6, 12])
    tree.add([9, 1])
    tree.add([2, 7])
    tree.add([10, 19])
    print(tree.tree)
