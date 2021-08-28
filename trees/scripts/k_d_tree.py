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
                self.max_priority = priority

class KDTree:
    '''k-d tree implementation.'''
    def __init__(self, k):
        self.tree = []
        self.k = k

    def add(self, node, dim=0, subtree=None):
        '''Add a node to the tree.
        
        Args:
            node: node to be added.
            dim: in which dim to be compared. By default, it starts with dim 0
                at root node and increases by one for each level of the tree.
            subtree: optional; in which subtree to add the node. By standard,
                it is the tree itself.
        '''
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

    def knn(self, node, k):
        '''Search for the k nearest neighbors.
        
        Args:
            node: which node to search for the neighbors of.
            k: how many neighbors.

        Returns:
            bpq: bounded priority queue with the k-NN.
        '''
        bpq = BoundedPriorityQueue(k)
        self.knn_search(node, 0, self.tree, bpq)
        return bpq

    def knn_search(self, node, dim, subtree, bpq):
        '''k-NN recursion search.
        
        It searchs for k-NN in both left and right subtrees, but it won't when
        the entire subtree is not worth it, that is, it is not possible to have
        any nearest node.
        
        Args:
            node: which node to search for the neighbors of.
            dim: in which dim to be compared. By default, it starts with dim 0
                at root node and increases by one for each level of the tree.
            subtree: in which subtree to add the node.
            bpq: bounded priority queue.
        '''
        # If the subtree is empty
        if not len(subtree):
            return

        bpq.put_item(subtree[1], self.distance(node, subtree[1]))

        if node[dim] < subtree[1][dim]:
            # Recursively search the left subtree
            self.knn_search(node, (dim+1) % self.k, subtree[0], bpq)
            # If it is worth it to search the other subtree
            if not bpq.full() or \
                np.abs(node[dim] - subtree[1][dim]) < bpq.max_priority:
                self.knn_search(node, (dim+1) % self.k, subtree[2], bpq)
        else:
            self.knn_search(node, (dim+1) % self.k, subtree[2], bpq)
            if not bpq.full() or \
                np.abs(node[dim] - subtree[1][dim]) < bpq.max_priority:
                self.knn_search(node, (dim+1) % self.k, subtree[0], bpq)

    def distance(self, node_1, node_2):
        return np.linalg.norm(np.array(node_1) - np.array(node_2))

if __name__ == '__main__':
    tree = KDTree(2)
    tree.add([51, 75])
    tree.add([25, 40])
    tree.add([70, 70])
    tree.add([10, 30])
    tree.add([35, 90])
    tree.add([55, 1])
    tree.add([60, 80])
    tree.add([1, 10])
    tree.add([50, 50])
    print(tree.tree)
    print(tree.knn([52, 52], 3).queue)
