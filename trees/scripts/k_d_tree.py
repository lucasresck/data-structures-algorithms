class KDTree:
    def __init__(self, k):
        self.tree = []
        self.k = k

    def add(self, node, dim=0, subtree=None):
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
