from random import choice, randrange, randint

from algorithm.parameters import params
from representation.tree import Tree


def generate_tree(tree, genome, output, method, nodes, depth, max_depth,
                  depth_limit):
    """
    Recursive function to derive a tree using a given method.
    
    :param tree: An instance of the Tree class.
    :param genome: The list of all codons in a tree.
    :param output: The list of all terminal nodes in a subtree. This is
    joined to become the phenotype.
    :param method: A string of the desired tree derivation method,
    e.g. "full" or "random".
    :param nodes: The total number of nodes in the tree.
    :param depth: The depth of the current node.
    :param max_depth: The maximum depth of any node in the tree.
    :param depth_limit: The maximum depth the tree can expand to.
    :return: genome, output, nodes, depth, max_depth.
    """
        
    # Increment nodes and depth, set depth of current node.
    nodes += 1
    depth += 1
    tree.depth = depth

    # Find the productions possible from the current root.
    productions = params['BNF_GRAMMAR'].rules[tree.root]

    if depth_limit:
        # Set remaining depth.
        remaining_depth = depth_limit - depth
    
    else:
        remaining_depth = depth_limit
    
    # Find which productions can be used based on the derivation method.
    available = legal_productions(method, remaining_depth, tree.root,
                                  productions['choices'])
    
    # Randomly pick a production choice.
    chosen_prod = choice(available)

    # Find the index of the chosen production and set a matching codon based
    # on that index.
    prod_index = productions['choices'].index(chosen_prod)
    codon = randrange(productions['no_choices'],
                      params['BNF_GRAMMAR'].codon_size,
                      productions['no_choices']) + prod_index
    
    # Set the codon for the current node and append codon to the genome.
    tree.codon = codon
    genome.append(codon)
    
    # Initialise empty list of children for current node.
    tree.children = []

    for symbol in chosen_prod['choice']:
        # Iterate over all symbols in the chosen production.
        if symbol["type"] == "T":
            # The symbol is a terminal. Append new node to children.
            tree.children.append(Tree(symbol["symbol"], tree))
            
            # Append the terminal to the output list.
            output.append(symbol["symbol"])
        
        elif symbol["type"] == "NT":
            # The symbol is a non-terminal. Append new node to children.
            tree.children.append(Tree(symbol["symbol"], tree))
            
            # recurse on the new node.
            genome, output, nodes, d, max_depth = \
                generate_tree(tree.children[-1], genome, output, method,
                              nodes, depth, max_depth, depth_limit)

    NT_kids = [kid for kid in tree.children if kid.root in
               params['BNF_GRAMMAR'].non_terminals]

    if not NT_kids:
        # Then the branch terminates here
        depth += 1
        nodes += 1

    if depth > max_depth:
        # Set new maximum depth
        max_depth = depth
    
    return genome, output, nodes, depth, max_depth


def legal_productions(method, depth_limit, root, productions):
    """
    Returns the available production choices for a node given a specific
    depth limit.
    
    :param method: A string specifying the desired tree derivation method.
    Current methods are "random" or "full".
    :param depth_limit: The overall depth limit of the desired tree from the
    current node.
    :param root: The root of the current node.
    :param productions: The full list of production choices from the current
    root node.
    :return: The list of available production choices based on the specified
    derivation method.
    """

    # Get all information about root node
    root_info = params['BNF_GRAMMAR'].non_terminals[root]
    
    if method == "random":
        # Randomly build a tree.
        
        if not depth_limit:
            # There is no depth limit, any production choice can be used.
            available = productions
        
        elif depth_limit > params['BNF_GRAMMAR'].max_arity + 1:
            # If the depth limit is greater than the maximum arity of the
            # grammar, then any production choice can be used.
            available = productions

        elif depth_limit < 0:
            # If we have already surpassed the depth limit, then list the
            # choices with the shortest terminating path.
            available = root_info['min_path']
        
        else:
            # The depth limit is less than or equal to the maximum arity of
            # the grammar + 1. We have to be careful in selecting available
            # production choices lest we generate a tree which violates the
            # depth limit.
            available = [prod for prod in productions if prod['max_path'] <=
                         depth_limit - 1]

            if not available:
                # There are no available choices which do not violate the depth
                # limit. List the choices with the shortest terminating path.
                available = root_info['min_path']
    
    elif method == "full":
        # Build a "full" tree where every branch extends to the depth limit.
        
        if not depth_limit:
            # There is no depth limit specified for building a Full tree.
            # Raise an error as a depth limit HAS to be specified here.
            s = "representation.derivation.legal_productions\n" \
                "Error: Depth limit not specified for `Full` tree derivation."
            raise Exception(s)
        
        elif depth_limit > params['BNF_GRAMMAR'].max_arity + 1:
            # If the depth limit is greater than the maximum arity of the
            # grammar, then only recursive production choices can be used.
            available = root_info['recursive']

            if not available:
                # There are no recursive production choices for the current
                # rule. Pick any production choices.
                available = productions

        else:
            # The depth limit is less than or equal to the maximum arity of
            # the grammar + 1. We have to be careful in selecting available
            # production choices lest we generate a tree which violates the
            # depth limit.
            available = [prod for prod in productions if prod['max_path'] ==
                         depth_limit - 1]
                        
            if not available:
                # There are no available choices which extend exactly to the
                # depth limit. List the NT choices with the longest terminating
                # paths that don't violate the limit.
                available = [prod for prod in productions if prod['max_path']
                             < depth_limit - 1]

    return available

