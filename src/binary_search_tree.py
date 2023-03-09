class BinarySearchTreeNode:
    def __init__(self, node, node_pos):
        self.node = node
        self.node_pos = node_pos # set node position
        self.left = None
        self.right = None

    def add_node(self, node, node_pos):
        if node == self.node:
            return # node already exists

        if node < self.node:
            if self.left:
                self.left.add_node(node, node_pos)
            else:
                self.left = BinarySearchTreeNode(node,node_pos)
        else:
            if self.right:
                self.right.add_node(node, node_pos)
            else:
                self.right = BinarySearchTreeNode(node,node_pos)


    def search(self, val):
        # a match found means there are no gap alignments with the k-mers in the database that 
        # fall within the HSSP threshold score 
        if self.node == val:
            return self.node_pos

        if val < self.node:
            if self.left:
                return self.left.search(val)
            else:
                return False

        if val > self.node:
            if self.right:
                return self.right.search(val)
            else:
                return False

    def in_order_traversal(self):
        elements = []
        if self.left:
            elements += self.left.in_order_traversal()

        elements.append(self.node)

        if self.right:
            elements += self.right.in_order_traversal()

        return elements
