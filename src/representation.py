from math import floor
from random import choice, randint, randrange
from re import DOTALL, MULTILINE, finditer, match
from sys import maxsize

import numpy as np

from parameters import params
from mapper import mapper, map_ind_from_genome

"""
derivation.py
"""

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

"""
grammar.py
"""

class Grammar(object):
    """
    Parser for Backus-Naur Form (BNF) Context-Free Grammars.
    """

    def __init__(self, file_name):
        """
        Initialises an instance of the grammar class. This instance is used
        to parse a given file_name grammar.

        :param file_name: A specified BNF grammar file.
        """

        if file_name.endswith("pybnf"):
            # Use python filter for parsing grammar output as grammar output
            # contains indented python code.
            self.python_mode = True

        else:
            # No need to filter/interpret grammar output, individual
            # phenotypes can be evaluated as normal.
            self.python_mode = False

        # Initialise empty dict for all production rules in the grammar.
        # Initialise empty dict of permutations of solutions possible at
        # each derivation tree depth.
        self.rules, self.permutations = {}, {}

        # Initialise dicts for terminals and non terminals, set params.
        self.non_terminals, self.terminals = {}, {}
        self.start_rule, self.codon_size = None, params['CODON_SIZE']
        self.min_path, self.max_arity, self.min_ramp = None, None, None

        # Set regular expressions for parsing BNF grammar.
        self.ruleregex = '(?P<rulename><\S+>)\s*::=\s*(?P<production>(?:(?=\#)\#[^\r\n]*|(?!<\S+>\s*::=).+?)+)'
        self.productionregex = '(?=\#)(?:\#.*$)|(?!\#)\s*(?P<production>(?:[^\'\"\|\#]+|\'.*?\'|".*?")+)'
        self.productionpartsregex = '\ *([\r\n]+)\ *|([^\'"<\r\n]+)|\'(.*?)\'|"(.*?)"|(?P<subrule><[^>|\s]+>)|([<]+)'

        # Read in BNF grammar, set production rules, terminals and
        # non-terminals.
        self.read_bnf_file(file_name)

        # Check the minimum depths of all non-terminals in the grammar.
        self.check_depths()

        # Check which non-terminals are recursive.
        self.check_recursion(self.start_rule["symbol"], [])

        # Set the minimum path and maximum arity of the grammar.
        self.set_arity()

        # Generate lists of recursive production choices and shortest
        # terminating path production choices for each NT in the grammar.
        # Enables faster tree operations.
        self.set_grammar_properties()

        # Calculate the total number of derivation tree permutations and
        # combinations that can be created by a grammar at a range of depths.
        self.check_permutations()

        if params['MIN_INIT_TREE_DEPTH']:
            # Set the minimum ramping tree depth from the command line.
            self.min_ramp = params['MIN_INIT_TREE_DEPTH']

        elif hasattr(params['INITIALISATION'], "ramping"):
            # Set the minimum depth at which ramping can start where we can
            # have unique solutions (no duplicates).
            self.get_min_ramp_depth()

        if params['REVERSE_MAPPING_TARGET'] or params['TARGET_SEED_FOLDER']:
            # Initialise dicts for reverse-mapping GE individuals.
            self.concat_NTs, self.climb_NTs = {}, {}

            # Find production choices which can be used to concatenate
            # subtrees.
            self.find_concatenation_NTs()

    def read_bnf_file(self, file_name):
        """
        Read a grammar file in BNF format. Parses the grammar and saves a
        dict of all production rules and their possible choices.

        :param file_name: A specified BNF grammar file.
        :return: Nothing.
        """

        with open(file_name, 'r') as bnf:
            # Read the whole grammar file.
            content = bnf.read()

            for rule in finditer(self.ruleregex, content, DOTALL):
                # Find all rules in the grammar

                if self.start_rule is None:
                    # Set the first rule found as the start rule.
                    self.start_rule = {"symbol": rule.group('rulename'),
                                       "type": "NT"}

                # Create and add a new rule.
                self.non_terminals[rule.group('rulename')] = {
                    'id': rule.group('rulename'),
                    'min_steps': maxsize,
                    'expanded': False,
                    'recursive': True,
                    'b_factor': 0}

                # Initialise empty list of all production choices for this
                # rule.
                tmp_productions = []

                for p in finditer(self.productionregex,
                                  rule.group('production'), MULTILINE):
                    # Iterate over all production choices for this rule.
                    # Split production choices of a rule.

                    if p.group('production') is None or p.group(
                            'production').isspace():
                        # Skip to the next iteration of the loop if the
                        # current "p" production is None or blank space.
                        continue

                    # Initialise empty data structures for production choice
                    tmp_production, terminalparts = [], None

                    # special cases: GE_RANGE:dataset_n_vars will be
                    # transformed to productions 0 | 1 | ... |
                    # n_vars-1, and similar for dataset_n_is,
                    # dataset_n_os
                    GE_RANGE_regex = r'GE_RANGE:(?P<range>\w*)'
                    m = match(GE_RANGE_regex, p.group('production'))
                    if m:
                        try:
                            if m.group('range') == "dataset_n_vars":
                                # number of columns from dataset
                                n = params['FITNESS_FUNCTION'].n_vars
                            elif m.group('range') == "dataset_n_is":
                                # number of input symbols (see
                                # if_else_classifier.py)
                                n = params['FITNESS_FUNCTION'].n_is
                            elif m.group('range') == "dataset_n_os":
                                # number of output symbols
                                n = params['FITNESS_FUNCTION'].n_os
                            else:
                                # assume it's just an int
                                n = int(m.group('range'))
                        except (ValueError, AttributeError):
                            raise ValueError("Bad use of GE_RANGE: "
                                             + m.group())

                        for i in range(n):
                            # add a terminal symbol
                            tmp_production, terminalparts = [], None
                            symbol = {
                                "symbol": str(i),
                                "type": "T",
                                "min_steps": 0,
                                "recursive": False}
                            tmp_production.append(symbol)
                            if str(i) not in self.terminals:
                                self.terminals[str(i)] = \
                                    [rule.group('rulename')]
                            elif rule.group('rulename') not in \
                                self.terminals[str(i)]:
                                self.terminals[str(i)].append(
                                    rule.group('rulename'))
                            tmp_productions.append({"choice": tmp_production,
                                                    "recursive": False,
                                                    "NT_kids": False})
                        # don't try to process this production further
                        # (but later productions in same rule will work)
                        continue

                    for sub_p in finditer(self.productionpartsregex,
                                          p.group('production').strip()):
                        # Split production into terminal and non terminal
                        # symbols.

                        if sub_p.group('subrule'):
                            if terminalparts is not None:
                                # Terminal symbol is to be appended to the
                                # terminals dictionary.
                                symbol = {"symbol": terminalparts,
                                          "type": "T",
                                          "min_steps": 0,
                                          "recursive": False}
                                tmp_production.append(symbol)
                                if terminalparts not in self.terminals:
                                    self.terminals[terminalparts] = \
                                        [rule.group('rulename')]
                                elif rule.group('rulename') not in \
                                    self.terminals[terminalparts]:
                                    self.terminals[terminalparts].append(
                                        rule.group('rulename'))
                                terminalparts = None

                            tmp_production.append(
                                {"symbol": sub_p.group('subrule'),
                                 "type": "NT"})

                        else:
                            # Unescape special characters (\n, \t etc.)
                            if terminalparts is None:
                                terminalparts = ''
                            terminalparts += ''.join(
                                [part.encode().decode('unicode-escape') for
                                 part in sub_p.groups() if part])

                    if terminalparts is not None:
                        # Terminal symbol is to be appended to the terminals
                        # dictionary.
                        symbol = {"symbol": terminalparts,
                                  "type": "T",
                                  "min_steps": 0,
                                  "recursive": False}
                        tmp_production.append(symbol)
                        if terminalparts not in self.terminals:
                            self.terminals[terminalparts] = \
                                [rule.group('rulename')]
                        elif rule.group('rulename') not in \
                            self.terminals[terminalparts]:
                            self.terminals[terminalparts].append(
                                rule.group('rulename'))
                    tmp_productions.append({"choice": tmp_production,
                                            "recursive": False,
                                            "NT_kids": False})

                if not rule.group('rulename') in self.rules:
                    # Add new production rule to the rules dictionary if not
                    # already there.
                    self.rules[rule.group('rulename')] = {
                        "choices": tmp_productions,
                        "no_choices": len(tmp_productions)}

                    if len(tmp_productions) == 1:
                        # Unit productions.
                        print("Warning: Grammar contains unit production "
                              "for production rule", rule.group('rulename'))
                        print("         Unit productions consume GE codons.")
                else:
                    # Conflicting rules with the same name.
                    raise ValueError("lhs should be unique",
                                     rule.group('rulename'))

    def check_depths(self):
        """
        Run through a grammar and find out the minimum distance from each
        NT to the nearest T. Useful for initialisation methods where we
        need to know how far away we are from fully expanding a tree
        relative to where we are in the tree and what the depth limit is.

        :return: Nothing.
        """

        # Initialise graph and counter for checking minimum steps to Ts for
        # each NT.
        counter, graph = 1, []

        for rule in sorted(self.rules.keys()):
            # Iterate over all NTs.
            choices = self.rules[rule]['choices']

            # Set branching factor for each NT.
            self.non_terminals[rule]['b_factor'] = self.rules[rule][
                'no_choices']

            for choice in choices:
                # Add a new edge to our graph list.
                graph.append([rule, choice['choice']])

        while graph:
            removeset = set()
            for edge in graph:
                # Find edges which either connect to terminals or nodes
                # which are fully expanded.
                if all([sy["type"] == "T" or
                        self.non_terminals[sy["symbol"]]['expanded']
                        for sy in edge[1]]):
                    removeset.add(edge[0])

            for s in removeset:
                # These NTs are now expanded and have their correct minimum
                # path set.
                self.non_terminals[s]['expanded'] = True
                self.non_terminals[s]['min_steps'] = counter

            # Create new graph list and increment counter.
            graph = [e for e in graph if e[0] not in removeset]
            counter += 1

    def check_recursion(self, cur_symbol, seen):
        """
        Traverses the grammar recursively and sets the properties of each rule.

        :param cur_symbol: symbol to check.
        :param seen: Contains already checked symbols in the current traversal.
        :return: Boolean stating whether or not cur_symbol is recursive.
        """

        if cur_symbol not in self.non_terminals.keys():
            # Current symbol is a T.
            return False

        if cur_symbol in seen:
            # Current symbol has already been seen, is recursive.
            return True

        # Append current symbol to seen list.
        seen.append(cur_symbol)

        # Get choices of current symbol.
        choices = self.rules[cur_symbol]['choices']
        nt = self.non_terminals[cur_symbol]

        recursive = False
        for choice in choices:
            for sym in choice['choice']:
                # Recurse over choices.
                recursive_symbol = self.check_recursion(sym["symbol"], seen)
                recursive = recursive or recursive_symbol

        # Set recursive properties.
        nt['recursive'] = recursive
        seen.remove(cur_symbol)

        return nt['recursive']

    def set_arity(self):
        """
        Set the minimum path of the grammar, i.e. the smallest legal
        solution that can be generated.

        Set the maximum arity of the grammar, i.e. the longest path to a
        terminal from any non-terminal.

        :return: Nothing
        """

        # Set the minimum path of the grammar as the minimum steps to a
        # terminal from the start rule.
        self.min_path = self.non_terminals[self.start_rule["symbol"]][
            'min_steps']

        # Set the maximum arity of the grammar as the longest path to
        # a T from any NT.
        self.max_arity = max(self.non_terminals[NT]['min_steps']
                             for NT in self.non_terminals)

        # Add the minimum terminal path to each production rule.
        for rule in self.rules:
            for choice in self.rules[rule]['choices']:
                NT_kids = [i for i in choice['choice'] if i["type"] == "NT"]
                if NT_kids:
                    choice['NT_kids'] = True
                    for sym in NT_kids:
                        sym['min_steps'] = self.non_terminals[sym["symbol"]][
                            'min_steps']

        # Add boolean flag indicating recursion to each production rule.
        for rule in self.rules:
            for prod in self.rules[rule]['choices']:
                for sym in [i for i in prod['choice'] if i["type"] == "NT"]:
                    sym['recursive'] = self.non_terminals[sym["symbol"]][
                        'recursive']
                    if sym['recursive']:
                        prod['recursive'] = True

    def set_grammar_properties(self):
        """
        Goes through all non-terminals and finds the production choices with
        the minimum steps to terminals and with recursive steps.

        :return: Nothing
        """

        for nt in self.non_terminals:
            # Loop over all non terminals.
            # Find the production choices for the current NT.
            choices = self.rules[nt]['choices']

            for choice in choices:
                # Set the maximum path to a terminal for each produciton choice
                choice['max_path'] = max([item["min_steps"] for item in
                                      choice['choice']])

            # Find shortest path to a terminal for all production choices for
            # the current NT. The shortest path will be the minimum of the
            # maximum paths to a T for each choice over all chocies.
            min_path = min([choice['max_path'] for choice in choices])

            # Set the minimum path in the self.non_terminals dict.
            self.non_terminals[nt]['min_path'] = [choice for choice in
                                                  choices if choice[
                                                      'max_path'] == min_path]

            # Find recursive production choices for current NT. If any
            # constituent part of a production choice is recursive,
            # it is added to the recursive list.
            self.non_terminals[nt]['recursive'] = [choice for choice in
                                                   choices if choice[
                                                       'recursive']]

    def check_permutations(self):
        """
        Calculates how many possible derivation tree combinations can be
        created from the given grammar at a specified depth. Only returns
        possible combinations at the specific given depth (if there are no
        possible permutations for a given depth, will return 0).

        :param ramps:
        :return: Nothing.
        """

        # Set the number of depths permutations are calculated for
        # (starting from the minimum path of the grammar)
        ramps = params['PERMUTATION_RAMPS']

        perms_list = []
        if self.max_arity > self.min_path:
            for i in range(max((self.max_arity + 1 - self.min_path), ramps)):
                x = self.check_all_permutations(i + self.min_path)
                perms_list.append(x)
                if i > 0:
                    perms_list[i] -= sum(perms_list[:i])
                    self.permutations[i + self.min_path] -= sum(perms_list[:i])
        else:
            for i in range(ramps):
                x = self.check_all_permutations(i + self.min_path)
                perms_list.append(x)
                if i > 0:
                    perms_list[i] -= sum(perms_list[:i])
                    self.permutations[i + self.min_path] -= sum(perms_list[:i])

    def check_all_permutations(self, depth):
        """
        Calculates how many possible derivation tree combinations can be
        created from the given grammar at a specified depth. Returns all
        possible combinations at the specific given depth including those
        depths below the given depth.

        :param depth: A depth for which to calculate the number of
        permutations of solution that can be generated by the grammar.
        :return: The permutations possible at the given depth.
        """

        if depth < self.min_path:
            # There is a bug somewhere that is looking for a tree smaller than
            # any we can create
            s = "representation.grammar.Grammar.check_all_permutations\n" \
                "Error: cannot check permutations for tree smaller than the " \
                "minimum size."
            raise Exception(s)

        if depth in self.permutations.keys():
            # We have already calculated the permutations at the requested
            # depth.
            return self.permutations[depth]

        else:
            # Calculate permutations at the requested depth.
            # Initialise empty data arrays.
            pos, depth_per_symbol_trees, productions = 0, {}, []

            for NT in self.non_terminals:
                # Iterate over all non-terminals to fill out list of
                # productions which contain non-terminal choices.
                a = self.non_terminals[NT]

                for rule in self.rules[a['id']]['choices']:
                    if rule['NT_kids']:
                        productions.append(rule)

            # Get list of all production choices from the start symbol.
            start_symbols = self.rules[self.start_rule["symbol"]]['choices']

            for choice in productions:
                # Generate a list of the symbols of each production choice
                key = str([sym['symbol'] for sym in choice['choice']])

                # Initialise permutations dictionary with the list
                depth_per_symbol_trees[key] = {}

            for i in range(2, depth + 1):
                # Find all the possible permutations from depth of min_path up
                # to a specified depth

                for choice in productions:
                    # Iterate over all production choices
                    sym_pos = 1

                    for j in choice['choice']:
                        # Iterate over all symbols in a production choice.
                        symbol_arity_pos = 0

                        if j["type"] is "NT":
                            # We are only interested in non-terminal symbols
                            for child in self.rules[j["symbol"]]['choices']:
                                # Iterate over all production choices for
                                # each NT symbol in the original choice.

                                if len(child['choice']) == 1 and \
                                   child['choice'][0]["type"] == "T":
                                    # If the child choice leads directly to
                                    # a single terminal, increment the
                                    # permutation count.
                                    symbol_arity_pos += 1

                                else:
                                    # The child choice does not lead
                                    # directly to a single terminal.
                                    # Generate a key for the permutations
                                    # dictionary and increment the
                                    # permutations count there.
                                    key = [sym['symbol'] for sym in child['choice']]
                                    if (i - 1) in depth_per_symbol_trees[str(key)].keys():
                                        symbol_arity_pos += depth_per_symbol_trees[str(key)][i - 1]

                            # Multiply original count by new count.
                            sym_pos *= symbol_arity_pos

                    # Generate new key for the current production choice and
                    # set the new value in the permutations dictionary.
                    key = [sym['symbol'] for sym in choice['choice']]
                    depth_per_symbol_trees[str(key)][i] = sym_pos

            # Calculate permutations for the start symbol.
            for sy in start_symbols:
                key = [sym['symbol'] for sym in sy['choice']]
                if str(key) in depth_per_symbol_trees:
                    pos += depth_per_symbol_trees[str(key)][depth] if depth in depth_per_symbol_trees[str(key)] else 0
                else:
                    pos += 1

            # Set the overall permutations dictionary for the current depth.
            self.permutations[depth] = pos

            return pos

    def get_min_ramp_depth(self):
        """
        Find the minimum depth at which ramping can start where we can have
        unique solutions (no duplicates).

        :param self: An instance of the representation.grammar.grammar class.
        :return: The minimum depth at which unique solutions can be generated
        """

        max_tree_depth = params['MAX_INIT_TREE_DEPTH']
        size = params['POPULATION_SIZE']

        # Specify the range of ramping depths
        depths = range(self.min_path, max_tree_depth + 1)

        if size % 2:
            # Population size is odd
            size += 1

        if size / 2 < len(depths):
            # The population size is too small to fully cover all ramping
            # depths. Only ramp to the number of depths we can reach.
            depths = depths[:int(size / 2)]

        # Find the minimum number of unique solutions required to generate
        # sufficient individuals at each depth.
        unique_start = int(floor(size / len(depths)))
        ramp = None

        for i in sorted(self.permutations.keys()):
            # Examine the number of permutations and combinations of unique
            # solutions capable of being generated by a grammar across each
            # depth i.
            if self.permutations[i] > unique_start:
                # If the number of permutations possible at a given depth i is
                # greater than the required number of unique solutions,
                # set the minimum ramp depth and break out of the loop.
                ramp = i
                break
        self.min_ramp = ramp

    def find_concatenation_NTs(self):
        """
        Scour the grammar class to find non-terminals which can be used to
        combine/reduce_trees derivation trees. Build up a list of such
        non-terminals. A concatenation non-terminal is one in which at least
        one production choice contains multiple non-terminals. For example:

            <e> ::= (<e><o><e>)|<v>

        is a concatenation NT, since the production choice (<e><o><e>) can
        reduce_trees multiple NTs together. Note that this choice also includes
        a combination of terminals and non-terminals.

        :return: Nothing.
        """

        # Iterate over all non-terminals/production rules.
        for rule in sorted(self.rules.keys()):

            # Find rules which have production choices leading to NTs.
            concat = [choice for choice in self.rules[rule]['choices'] if
                      choice['NT_kids']]

            if concat:
                # We can reduce_trees NTs.
                for choice in concat:

                    symbols = [[sym['symbol'], sym['type']] for sym in
                               choice['choice']]

                    NTs = [sym['symbol'] for sym in choice['choice'] if
                           sym['type'] == "NT"]

                    for NT in NTs:
                        # We add to our self.concat_NTs dictionary. The key is
                        # the root node we want to reduce_trees with another
                        # node. This way when we have a node and wish to see
                        # if we can reduce_trees it with anything else, we
                        # simply look up this dictionary.
                        conc = [choice['choice'], rule, symbols]

                        if NT not in self.concat_NTs:
                            self.concat_NTs[NT] = [conc]
                        else:
                            if conc not in self.concat_NTs[NT]:
                                self.concat_NTs[NT].append(conc)

    def __str__(self):
        return "%s %s %s %s" % (self.terminals, self.non_terminals,
                                self.rules, self.start_rule)

"""
tree.py
"""

class Tree:

    def __init__(self, expr, parent):
        """
        Initialise an instance of the tree class.
        
        :param expr: A non-terminal from the params['BNF_GRAMMAR'].
        :param parent: The parent of the current node. None if node is tree
        root.
        """
        
        self.parent = parent
        self.codon = None
        self.depth = 1
        self.root = expr
        self.children = []
        self.snippet = None

    def __str__(self):
        """
        Builds a string of the current tree.
        
        :return: A string of the current tree.
        """
        
        # Initialise the output string.
        result = "("
        
        # Append the root of the current node to the output string.
        result += str(self.root)
        
        for child in self.children:
            # Iterate across all children.
            
            if len(child.children) > 0:
                # Recurse through all children.
                result += " " + str(child)
            
            else:
                # Child is a terminal, append root to string.
                result += " " + str(child.root)
        
        result += ")"
        
        return result

    def __copy__(self):
        """
        Creates a new unique copy of self.
        
        :return: A new unique copy of self.
        """

        # Copy current tree by initialising a new instance of the tree class.
        tree_copy = Tree(self.root, self.parent)
        
        # Set node parameters.
        tree_copy.codon, tree_copy.depth = self.codon, self.depth

        tree_copy.snippet = self.snippet

        for child in self.children:
            # Recurse through all children.
            new_child = child.__copy__()
            
            # Set the parent of the copied child as the copied parent.
            new_child.parent = tree_copy
            
            # Append the copied child to the copied parent.
            tree_copy.children.append(new_child)

        return tree_copy

    def __eq__(self, other, same=True):
        """
        Set the definition for comparison of two instances of the tree
        class by their attributes. Returns True if self == other.

        :param other: Another instance of the tree class with which to compare.
        :return: True if self == other.
        """

        # Get attributes of self and other.
        a_self, a_other = vars(self), vars(other)
                
        # Don't look at the children as they are class instances themselves.
        taboo = ["parent", "children", "snippet", "id"]
        self_no_kids = {k: v for k, v in a_self.items() if k not in taboo}
        other_no_kids = {k: v for k, v in a_other.items() if k not in taboo}
                
        # Compare attributes
        if self_no_kids != other_no_kids:
            # Attributes are not the same.
            return False
            
        else:
            # Attributes are the same
            child_list = [self.children, other.children]
            
            if len(list(filter(lambda x: x is not None, child_list))) % 2 != 0:
                # One contains children, the other doesn't.
                return False

            elif self.children and len(self.children) != len(other.children):
                # Number of children differs between self and other.
                return False

            elif self.children:
                # Compare children recursively.
                for i, child in enumerate(self.children):
                    same = child.__eq__(other.children[i], same)

        return same

    def get_target_nodes(self, array, target=None):
        """
        Returns the all NT nodes which match the target NT list in a
        given tree.
        
        :param array: The array of all nodes that match the target.
        :param target: The target nodes to match.
        :return: The array of all nodes that match the target.
        """
            
        if self.root in target:
            # Check if the current node matches the target.
            
            # Add the current node to the array.
            array.append(self)
        
        # Find all non-terminal children of the current node.
        NT_kids = [kid for kid in self.children if kid.root in
                   params['BNF_GRAMMAR'].non_terminals]
        
        for child in NT_kids:
            if NT_kids:
                # Recursively call function on any non-terminal children.
                array = child.get_target_nodes(array, target=target)
        
        return array

    def get_node_labels(self, labels):
        """
        Recurses through a tree and appends all node roots to a set.
        
        :param labels: The set of roots of all nodes in the tree.
        :return: The set of roots of all nodes in the tree.
        """
        
        # Add the current root to the set of all labels.
        labels.add(self.root)

        for child in self.children:
            # Recurse on all children.
            labels = child.get_node_labels(labels)
        
        return labels

    def get_tree_info(self, nt_keys, genome, output, invalid=False,
                      max_depth=0, nodes=0):
        """
        Recurses through a tree and returns all necessary information on a
        tree required to generate an individual.
        
        :param genome: The list of all codons in a subtree.
        :param output: The list of all terminal nodes in a subtree. This is
        joined to become the phenotype.
        :param invalid: A boolean flag for whether a tree is fully expanded.
        True if invalid (unexpanded).
        :param nt_keys: The list of all non-terminals in the grammar.
        :param nodes: the number of nodes in a tree.
        :param max_depth: The maximum depth of any node in the tree.
        :return: genome, output, invalid, max_depth, nodes.
        """

        # Increment number of nodes in tree and set current node id.
        nodes += 1
        
        if self.parent:
            # If current node has a parent, increment current depth from
            # parent depth.
            self.depth = self.parent.depth + 1
        
        else:
            # Current node is tree root, set depth to 1.
            self.depth = 1
        
        if self.depth > max_depth:
            # Set new max tree depth.
            max_depth = self.depth

        if self.codon:
            # If the current node has a codon, append it to the genome.
            genome.append(self.codon)

        # Find all non-terminal children of current node.
        NT_children = [child for child in self.children if child.root in
                       nt_keys]
        
        if not NT_children:
            # The current node has only terminal children, increment number
            # of tree nodes.
            nodes += 1

            # Terminal children increase the current node depth by one.
            # Check the recorded max_depth.
            if self.depth + 1 > max_depth:
                # Set new max tree depth.
                max_depth = self.depth + 1

        if self.root in nt_keys and not self.children:
            # Current NT has no children. Invalid tree.
            invalid = True

        for child in self.children:
            # Recurse on all children.

            if not child.children:
                # If the current child has no children it is a terminal.
                # Append it to the phenotype output.
                output.append(child.root)
                
                if child.root in nt_keys:
                    # Current non-terminal node has no children; invalid tree.
                    invalid = True
            
            else:
                # The current child has children, recurse.
                genome, output, invalid, max_depth, nodes = \
                    child.get_tree_info(nt_keys, genome, output, invalid,
                                        max_depth, nodes)

        return genome, output, invalid, max_depth, nodes

    def print_tree(self):
        """
        Prints out all nodes in the tree, indented according to node depth.
        
        :return: Nothing.
        """

        print(self.depth, "".join([" " for _ in range(self.depth)]), self.root)

        for child in self.children:
            if not child.children:
                print(self.depth + 1,
                      "".join([" " for _ in range(self.depth + 1)]),
                      child.root)
            else:
                child.print_tree()

"""
individual.py
"""

class Individual(object):
    """
    A GE individual.
    """

    def __init__(self, genome, ind_tree, map_ind=True):
        """
        Initialise an instance of the individual class (i.e. create a new
        individual).

        :param genome: An individual's genome.
        :param ind_tree: An individual's derivation tree, i.e. an instance
        of the representation.tree.Tree class.
        :param map_ind: A boolean flag that indicates whether or not an
        individual needs to be mapped.
        """

        if map_ind:
            # The individual needs to be mapped from the given input
            # parameters.
            self.phenotype, self.genome, self.tree, self.nodes, self.invalid, \
                self.depth, self.used_codons = mapper(genome, ind_tree)

        else:
            # The individual does not need to be mapped.
            self.genome, self.tree = genome, ind_tree

        
        self.fitness = params['FITNESS_FUNCTION'].default_fitness
        self.runtime_error = False
        self.name = None

    def __lt__(self, other):
        """
        Set the definition for comparison of two instances of the individual
        class by their fitness values. Allows for sorting/ordering of a
        population of individuals. Note that numpy NaN is used for invalid
        individuals and is used by some fitness functions as a default fitness.
        We implement a custom catch for these NaN values.

        :param other: Another instance of the individual class (i.e. another
        individual) with which to compare.
        :return: Whether or not the fitness of the current individual is
        greater than the comparison individual.
        """

        if np.isnan(self.fitness): return True
        elif np.isnan(other.fitness): return False
        else: return self.fitness < other.fitness if params['FITNESS_FUNCTION'].maximise else other.fitness < self.fitness

    def __le__(self, other):
        """
        Set the definition for comparison of two instances of the individual
        class by their fitness values. Allows for sorting/ordering of a
        population of individuals. Note that numpy NaN is used for invalid
        individuals and is used by some fitness functions as a default fitness.
        We implement a custom catch for these NaN values.

        :param other: Another instance of the individual class (i.e. another
        individual) with which to compare.
        :return: Whether or not the fitness of the current individual is
        greater than or equal to the comparison individual.
        """

        if np.isnan(self.fitness): return True
        elif np.isnan(other.fitness): return False
        else: return self.fitness <= other.fitness if params['FITNESS_FUNCTION'].maximise else other.fitness <= self.fitness

    def __str__(self):
        """
        Generates a string by which individuals can be identified. Useful
        for printing information about individuals.

        :return: A string describing the individual.
        """
        return ("Individual: " +
                str(self.phenotype) + "; " + str(self.fitness))

    def deep_copy(self):
        """
        Copy an individual and return a unique version of that individual.

        :return: A unique copy of the individual.
        """

        if not params['GENOME_OPERATIONS']:
            # Create a new unique copy of the tree.
            new_tree = self.tree.__copy__()

        else:
            new_tree = None

        # Create a copy of self by initialising a new individual.
        new_ind = Individual(self.genome.copy(), new_tree, map_ind=False)

        # Set new individual parameters (no need to map genome to new
        # individual).
        new_ind.phenotype, new_ind.invalid = self.phenotype, self.invalid
        new_ind.depth, new_ind.nodes = self.depth, self.nodes
        new_ind.used_codons = self.used_codons
        new_ind.runtime_error = self.runtime_error

        return new_ind

    def evaluate(self):
        """
        Evaluates phenotype in using the fitness function set in the params
        dictionary. For regression/classification problems, allows for
        evaluation on either training or test distributions. Sets fitness
        value.

        :return: Nothing unless multicore evaluation is being used. In that
        case, returns self.
        """

        # Evaluate fitness using specified fitness function.
        self.fitness = params['FITNESS_FUNCTION'](self)

        if params['MULTICORE']:
            return self


def map_tree_from_genome(genome):
    """
    Maps a full tree from a given genome.

    :param genome: A genome to be mapped.
    :return: All components necessary for a fully mapped individual.
    """

    # Initialise an instance of the tree class
    tree = Tree(str(params['BNF_GRAMMAR'].start_rule["symbol"]), None)

    # Map tree from the given genome
    output, used_codons, nodes, depth, max_depth, invalid = \
        genome_tree_map(tree, genome, [], 0, 0, 0, 0)

    # Build phenotype.
    phenotype = "".join(output)

    if invalid:
        # Return "None" phenotype if invalid
        return None, genome, tree, nodes, invalid, max_depth, \
           used_codons

    else:
        return phenotype, genome, tree, nodes, invalid, max_depth, \
           used_codons


def genome_tree_map(tree, genome, output, index, depth, max_depth, nodes,
                    invalid=False):
    """
    Recursive function which builds a tree using production choices from a
    given genome. Not guaranteed to terminate.

    :param tree: An instance of the representation.tree.Tree class.
    :param genome: A full genome.
    :param output: The list of all terminal nodes in a subtree. This is
    joined to become the phenotype.
    :param index: The index of the current location on the genome.
    :param depth: The current depth in the tree.
    :param max_depth: The maximum overall depth in the tree so far.
    :param nodes: The total number of nodes in the tree thus far.
    :param invalid: A boolean flag indicating whether or not the individual
    is invalid.
    :return: index, the index of the current location on the genome,
             nodes, the total number of nodes in the tree thus far,
             depth, the current depth in the tree,
             max_depth, the maximum overall depth in the tree,
             invalid, a boolean flag indicating whether or not the
             individual is invalid.
    """

    if not invalid and index < len(genome) * (params['MAX_WRAPS'] + 1):
        # If the solution is not invalid thus far, and if we still have
        # remaining codons in the genome, then we can continue to map the tree.

        if params['MAX_TREE_DEPTH'] and (max_depth > params['MAX_TREE_DEPTH']):
            # We have breached our maximum tree depth limit.
            invalid = True

        # Increment and set number of nodes and current depth.
        nodes += 1
        depth += 1
        tree.id, tree.depth = nodes, depth

        # Find all production choices and the number of those production
        # choices that can be made by the current root non-terminal.
        productions = params['BNF_GRAMMAR'].rules[tree.root]['choices']
        no_choices = params['BNF_GRAMMAR'].rules[tree.root]['no_choices']

        # Set the current codon value from the genome.
        tree.codon = genome[index % len(genome)]

        # Select the index of the correct production from the list.
        selection = tree.codon % no_choices

        # Set the chosen production
        chosen_prod = productions[selection]

        # Increment the index
        index += 1

        # Initialise an empty list of children.
        tree.children = []

        for symbol in chosen_prod['choice']:
            # Add children to the derivation tree by creating a new instance
            # of the representation.tree.Tree class for each child.

            if symbol["type"] == "T":
                # Append the child to the parent node. Child is a terminal, do
                # not recurse.
                tree.children.append(Tree(symbol["symbol"], tree))
                output.append(symbol["symbol"])

            elif symbol["type"] == "NT":
                # Append the child to the parent node.
                tree.children.append(Tree(symbol["symbol"], tree))

                # Recurse by calling the function again to map the next
                # non-terminal from the genome.
                output, index, nodes, d, max_depth, invalid = \
                    genome_tree_map(tree.children[-1], genome, output,
                                    index, depth, max_depth, nodes,
                                    invalid=invalid)

    else:
        # Mapping incomplete, solution is invalid.
        return output, index, nodes, depth, max_depth, True

    # Find all non-terminals in the chosen production choice.
    NT_kids = [kid for kid in tree.children if kid.root in
               params['BNF_GRAMMAR'].non_terminals]

    if not NT_kids:
        # There are no non-terminals in the chosen production choice, the
        # branch terminates here.
        depth += 1
        nodes += 1

    if not invalid:
        # The solution is valid thus far.

        if depth > max_depth:
            # Set the new maximum depth.
            max_depth = depth

        if params['MAX_TREE_DEPTH'] and (max_depth > params['MAX_TREE_DEPTH']):
            # If our maximum depth exceeds the limit, the solution is invalid.
            invalid = True

    return output, index, nodes, depth, max_depth, invalid
