from algorithm import mapper
from algorithm.parameters import params
from representation import individual


def check_ind(ind):
    """
    Checks all aspects of an individual to ensure everything is correct.
    
    :param ind: An instance of the representation.individaul.Individual class.
    :return: Nothing.
    """
    
    # Re-map individual using fast genome mapper to check everything is ok
    new_ind = individual.Individual(ind.genome, None)

    # Get attributes of both individuals.
    attributes_0 = vars(ind)
    attributes_1 = vars(new_ind)
    
    if params['GENOME_OPERATIONS']:
        # If this parameter is set then the new individual will have no tree.
        attributes_0['tree'] = None
    
    else:
        if attributes_0['tree'] != attributes_1['tree']:
            print("Error: utilities.check_methods.check_ind."
                  "Individual trees do not match.")

    # Must remove the tree variable as instances of classes cannot be
    # directly compared at the moment.
    attributes_0['tree'], attributes_1['tree'] = None, None
    
    if attributes_0 != attributes_1:
        print("Error: utilities.check_methods.check_ind."
              "Individual attributes do not match correct attributes.")
        print("Original attributes:\n", attributes_0)
        print("Correct attributes:\n", attributes_1)
        quit()


def print_dual_info(p1, g1, n1, i1, d1, c1,
                    p2, n2, i2, d2, c2):
    """Print information about two sets of data."""
    
    print("Phenotypes\n  ", p1, "\n  ", p2)
    print("Nodes\n  ", n1, "\n  ", n2)
    print("Depth\n  ", d1, "\n  ", d2)
    print("Invalid\n  ", i1, "\n  ", i2)
    print("Used Codons\n  ", c1, "\n  ", c2)
    print("Genome\n  ", g1)


def check_mapping(genome):
    """
    Checks both mapping methods to ensure a single genome gives the same
    results. Checks mapper.map_tree_from_genome and mapper.map_ind_from_genome.
    
    :param genome: the genome of an individual
    :return: Nothing.
    """

    p0, g0, t0, n0, i0, d0, c0 = mapper.map_tree_from_genome(genome)
    p1, g1, t1, n1, i1, d1, c1 = mapper.map_ind_from_genome(genome)

    if not i0:
        if i0 == i1:
            # Both individuals are valid.
            
            if p1 != p0:
                print("Error: utilities.check_methods.check_mapping."
                      "Phenotypes don't match")
                print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
                quit()
            
            elif n1 != n0:
                print("Error: utilities.check_methods.check_mapping."
                      "Nodes don't match")
                print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
                quit()
            
            elif d1 != d0:
                print("Error: utilities.check_methods.check_mapping."
                      "Tree depth doesn't match")
                print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
                quit()
            
            elif c1 != c0:
                print("Error: utilities.check_methods.check_mapping."
                      "Used codons don't match")
                print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
                quit()
            
            elif d1 > params['MAX_TREE_DEPTH']:
                print("Error: utilities.check_methods.check_mapping."
                      "Max tree depth limit exceeded")
                print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
                quit()
        
        else:
            print("Error: utilities.check_methods.check_mapping."
                  "Invalid doesn't match.")
            print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
            quit()

    else:
        print("Error: utilities.check_methods.check_mapping."
              "Individual is not valid.")
        print_dual_info(p0, g0, n0, i0, d0, c0, p1, n1, i1, d1, c1)
        quit()


def check_genome_from_tree(ind_tree):
    """
    Goes through a tree and checks each codon to ensure production choice is
    correct.
    
    :param ind_tree: The representation.tree.Tree class derivation tree of
    an individual.
    :return: Nothing.
    """

    if ind_tree.children:
        # This node has children and thus must have an associated codon.
        
        if not ind_tree.codon:
            print("Error: "
                  "utilities.check_methods.check_genome_from_tree. "
                  "Node with children has no codon.")
            print(ind_tree.children)
            quit()
        
        # Check production choices for node root.
        productions = params['BNF_GRAMMAR'].rules[ind_tree.root]
        
        # Select choice based on node codon.
        selection = ind_tree.codon % len(productions)
        chosen_prod = productions[selection]
        
        # Build list of roots of the chosen production.
        prods = [prod[0] for prod in chosen_prod]
        roots = []
        
        # Build list of the roots of all node children.
        for kid in ind_tree.children:
            roots.append(kid.root)
        
        # Match production roots with children roots.
        if roots != prods:
            print("Error: "
                  "utilities.check_methods.check_genome_from_tree. "
                  "Codons are incorrect for given tree.")
            print("Codon productions:\t", prods)
            print("Actual children:  \t", roots)
            quit()
    
    for kid in ind_tree.children:
        # Recurse over all children.
        check_genome_from_tree(kid)


def check_expansion(tree, nt_keys):
    """
    Check if a given tree is completely expanded or not. Return boolean
    True if the tree IS NOT completely expanded, i.e. if tree is invalid.
    
    :param tree: An individual's derivation tree.
    :param nt_keys: The list of all non-terminals.
    :return: True if tree is not fully expanded, else False.
    """
    
    check = False
    if tree.root in nt_keys:
        # Current node is a NT and should have children
        if tree.children:
            # Everything is as expected
            for child in tree.children:
                # Recurse over all children.
                check = child.check_expansion(nt_keys)
                
                if check:
                    # End recursion.
                    break
        
        else:
            # Current node is not completely expanded.
            check = True
    
    return check


def build_genome(tree, genome):
    """
    Goes through a tree and builds a genome from all codons in the subtree.

    :param tree: An individual's derivation tree.
    :param genome: The list of all codons in a subtree.
    :return: The fully built genome of a subtree.
    """
    
    if tree.codon:
        # If the current node has a codon, append it to the genome.
        genome.append(tree.codon)
    
    for child in tree.children:
        # Recurse on all children.
        genome = child.build_genome(genome)
    
    return genome


def get_nodes_and_depth(tree, nodes=0, max_depth=0):
    """
    Get the number of nodes and the max depth of the tree.
    
    :param tree: An individual's derivation tree.
    :param nodes: The number of nodes in a tree.
    :param max_depth: The maximum depth of any node in the tree.
    :return: number, max_depth.
    """
    
    # Increment number of nodes in the tree.
    nodes += 1

    # Set the depth of the current node.
    if tree.parent:
        tree.depth = tree.parent.depth + 1
    else:
        tree.depth = 1
        
    # Check the recorded max_depth.
    if tree.depth > max_depth:
        max_depth = tree.depth
        
    # Create list of all non-terminal children of current node.
    NT_kids = [kid for kid in tree.children if kid.root in
               params['BNF_GRAMMAR'].non_terminals]
    
    if not NT_kids:
        # Current node has only terminal children.
        nodes += 1
    
    else:
        for child in NT_kids:
            # Recurse over all children.
            nodes, max_depth = get_nodes_and_depth(child, nodes, max_depth)
    
    return nodes, max_depth


def get_max_tree_depth(tree, max_depth=1):
    """
    Returns the maximum depth of the tree from the current node.

    :param tree: The tree we wish to find the maximum depth of.
    :param max_depth: The maximum depth of the tree.
    :return: The maximum depth of the tree.
    """
    
    curr_depth = get_current_depth(tree)
    if curr_depth > max_depth:
        max_depth = curr_depth
    for child in tree.children:
        max_depth = get_max_tree_depth(child, max_depth)
    return max_depth


def get_current_depth(tree):
    """
    Get the depth of the current node by climbing back up the tree until no
    parents remain (i.e. the root node has been reached).

    :param tree: An individual's derivation tree.
    :return: The depth of the current node.
    """
    
    # Set the initial depth at 1.
    depth = 1
    
    # Set the current parent.
    current_parent = tree.parent
    
    while current_parent is not None:
        # Recurse until the root node of the tree has been reached.
        
        # Increment depth.
        depth += 1
        
        # Set new parent.
        current_parent = current_parent.parent
    
    return depth


def get_output(ind_tree):
    """
    Calls the recursive build_output(self) which returns a list of all
    node roots. Joins this list to create the full phenotype of an
    individual. This two-step process speeds things up as it only joins
    the phenotype together once rather than at every node.

    :param ind_tree: a full tree for which the phenotype string is to be built.
    :return: The complete built phenotype string of an individual.
    """
    
    def build_output(tree):
        """
        Recursively adds all node roots to a list which can be joined to
        create the phenotype.

        :return: The list of all node roots.
        """
        
        output = []
        for child in tree.children:
            if not child.children:
                # If the current child has no children it is a terminal.
                # Append it to the output.
                output.append(child.root)
            
            else:
                # Otherwise it is a non-terminal. Recurse on all
                # non-terminals.
                output += build_output(child)
        
        return output
    
    return "".join(build_output(ind_tree))


def ret_true(obj):
    """
    Returns "True" if an object is there. E.g. if given a list, will return
    True if the list contains some data, but False if the list is empty.
    
    :param obj: Some object (e.g. list)
    :return: True if something is there, else False.
    """

    if obj:
        return True
    else:
        return False
