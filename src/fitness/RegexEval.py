import re
import time, timeit
import traceback
import sys
from algorithm.parameters import params
from representation import individual
from multiprocessing import Process, Queue
#from multiprocessing.pool import ThreadPool # could be better performance to keep the thread around (is overhead cost of recreating a thread, less than the cost of exponential regex? no)

# http://stackoverflow.com/questions/24812253/how-can-i-capture-return-value-with-python-timeit-module/
timeit.template = """
def inner(_it, _timer{init}):
    {setup}
    _t0 = _timer()
    for _i in _it:
        retval = {stmt}
    _t1 = _timer()
    return _t1 - _t0, retval
"""

"""
TODO
Apache log file
([0-9a-f.:]+)\s+(-|.+?)\s+(-|.+?)\s+\[([0-9]{2}\/[a-z]{3}\/[0-9]{4}\:[0-9]{2}:[0-9]{2}:[0-9]{2}[^\]]*)\] \"(\S+?)\s(\S*?)\s{0,1}(\S+?)\" ([0-9|\-]+) ([0-9|\-]+)

"""

# pool = ThreadPool(processes=1)

class RegexEval:
    """
    Fitness function for checking regex matching which sums functionality error.
    The regex is presented with a number of strings.
    The resulting matched values are checked against known correct answers
    """

    maximise = False
    
    def __init__(self):
        self.test_cases = list()
        self.generate_tests()
        self.time=True


    def call_fitness(self, individual, q):
        regex_string = individual.phenotype
#        regex_string = "^(?=.{1,254}$)(?=.{1,64}@)[-!#$%&'*+0-9=?A-Z^_`a-z{|}~]+(\.[-!#$%&'*+0-9=?A-Z^_`a-z{|}~]+)*@[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?(\.[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?)*$"
#        regex_string = "^\s*(-|\+)?(\d+|(\d*(\.\d*)))([eE][+-]?\d+)?\s*$"
        #regex_string = "\d.*\w.$|e"
        #        regex_string="!|\w.[6-l]\w\w[^!]\w\w[^!]\w\w[^!]\d\w[^P]\w\w"
        try:
            compiled_regex = re.compile(regex_string)
            eval_results = self.test_regex(compiled_regex)
            result_error,time_sum = self.calculate_fitness(eval_results)
            #            if "5" in regex_string:
            # (we should use multi-objective/pareto front)
            #print(regex_string + ": {}".format(fitness))
            #sys.exit()
            fitness = result_error + time_sum
#            if fitness == 0 : # no error
#                fitness += time_sum
            if 'SEED_GENOME' in params and params['SEED_GENOME']:
                similarity_score = self.calculate_similarity_score(regex_string)
                #if(similarity_score > 1):
                #    fitness += (similarity_score)
                #else:
                #    fitness += time_sum


            # phenotype or genome length?
            # return fitness + (len(individual.phenotype)/10000) #fitness # error is first, length second
            # if fitness >= 1: # If there is a functionality error, we don't really care about time 
            #     return fitness + (len(individual.genome)/10000) #fitness # error is first, length second
            #     #return fitness + (len(regex_string)/10000) #fitness # error is first, length second
            # else:
            #     self.time=True
            #     return fitness
            q.put(fitness)
#            print("\n {}".format(fitness))

        except:
            #print(" ------------ exception")
#            print(traceback.format_exc())
#            sys.exit()
            q.put(100000)

    def calculate_similarity_score(self, regex_string):
        char_same_count = 0
        if 'SEED_GENOME' in params and params['SEED_GENOME']:
            seed_ind = individual.Individual(params['SEED_GENOME'], None)
            seed_string = seed_ind.phenotype
            if regex_string == seed_string:
                return len(seed_string)
            for i in range(len(regex_string)):
                for j in range(i,len(seed_string)):
                    if regex_string[i] == (seed_string)[j]:
                        char_same_count += 1
                        break
        #char_same_count += (abs(len(seed_string) - len(regex_string))/100)
        return char_same_count # / len(seed_string)
         

    def calculate_fitness(self,eval_results):
        result_error=0
        time_sum=0.0
        # print(" ")
        # print("Calculate_fitness ")
        for a_result in eval_results:
            time_sum += a_result[0] # /  a_result[2]
            # print("time_val : {} ".format(a_result[0]))
            # if a_result[1] == None: # no match
            # result_error += 100 * (len(a_result[3].search_string)) #+ len(a_result[3].matched_string))
            # else: # a match which may be the empty string
            result_error += a_result[3].calc_match_errors(list(a_result[1]))
        #if result_error >1:
        #    fitness = result_error
        #else:
        #    fitness = time_sum
        # fitness = result_error + time_sum
        # if fitness == seed_fitness:
        # fitness = 100 * len(a_result) # identical result to seed penalised (plucking the centre from spiderweb)
        return result_error,time_sum

    def test_regex(self,compiled_regex):
        results = list()
        testing_iterations=1
        # do a quick test to time the longest test case (which is also the last in the list)
        quick_test = self.time_regex_test_case(compiled_regex, self.test_cases[len(self.test_cases)-1], testing_iterations)
        #        print("quick_test time: {}".format(quick_test[0]))
        # aim for entire test suite to take less than a second
        if quick_test[3].calc_match_errors(list(quick_test[1])) < 0 : # Ideally we only time a program if it is funtionally correct
            # eval_time = .05 # seconds
            testing_iterations = 10000000 # int(( eval_time / (quick_test[0]/10))/len(self.test_cases)) # change to half second?
        # print("Iterations {}".format(testing_iterations))
        for test_case in self.test_cases:
            results.append(self.time_regex_test_case(compiled_regex, test_case, testing_iterations))
        return results

    def time_regex_test_case(self, compiled_regex, test_case, iterations):
        repeats = 3
        iterations_per_repeat = iterations
        search_string = test_case.get_search_string()
        def wrap():
            # Timing bug, lazy eval defers computation if we don't convert to list here
            # https://swizec.com/blog/python-and-lazy-evaluation/swizec/5148
            return list(compiled_regex.finditer(search_string)) 
        
        t = timeit.Timer(wrap)

        repeat_iterations = t.repeat(repeat=repeats, number = iterations)
        
        best_run = list(repeat_iterations[0])
        for repeated_timeit in repeat_iterations:
            if best_run[0] > list(repeated_timeit)[0]:
                best_run = list(repeated_timeit)
        #return_tuple = t.timeit(number = iterations)
        #return_vals = t.timeit()
        return_vals = list(best_run)
        return_vals.append(iterations)
        return_vals.append(test_case)
        return return_vals

    """ 
    This is a 'booster' for test suite generation. We know a single good match, 
    and we can use that to find search strings which do not contain a match.
    Given a regex, generate/discoverB a test suite of examples which match, and those that don't.
    The test suite is used to define (or outline) the functionality boundaries of the regex.
    When we go to evolve new regexs, we can use the test suite to measure functionality equivalence 
    with the original test regex.
    """
    def generate_equivalence_test_suite_replacement(self, a_match, a_regex):
        # go through the whole known search string, changing letters until you find one which does not match.
        compiled_regex = re.compile(a_regex)
        if len(a_match.matches) > 0 :
            for i in range(0, len(a_match.search_string)):
                for char in [ a for a in range(ord('0'), ord('9'))] + [ord('a'), ord('Z') ]:
                    new_search_string = a_match.search_string[:i] + chr(char) + a_match.search_string[i+1:]
                    a_test_case_string = RegexTestString(new_search_string)
                    vals = self.time_regex_test_case(compiled_regex, a_test_case_string, 1)
                    if len(list(vals[1])) == 0:
                        self.test_cases.append(a_test_case_string)


    def generate_equivalence_test_suite_length(self, a_match, a_regex):
        # improve performance of a bad regex.
        # taken from here: https://www.loggly.com/blog/regexes-the-bad-better-best/
        
        # best (by hand) regex: [12]\d{3}-[01]\d-[0-3]\d ([^ \[]*?)\[([^\]]*?)\]:.*
        """
        Technically this regex is still not the “best,” as there are even
        more optimizations that you can make. Those optimizations, in
        my opinion, greatly complicate the regex and take quite a lot
        of time and expertise to craft, but they do boost the
        performance of the regex by a non-trivial amount. The
        optimizations I describe in this post give by far the greatest
        bang for your buck performance-wise while maintaining the
        readability of the regex. If you’d like to learn more about
        how to optimize even further, check out this
        post:http://www.rexegg.com/regex-quantifiers.html#explicit_greed.
        """
        # worst regex: .* (.*)\[(.*)\]:.*
        # seed genome to go in regex.txt parameters file
        # [10082, 41053, 89472, 13943, 83234, 23739, 51472, 59723, 52910, 4880, 72071, 68514, 53866, 25754, 39758, 51472, 59723, 29655, 22034, 52910, 4880, 72071, 5630, 83669, 54528, 45876, 76718, 53866, 25754, 39758, 51472, 59723, 29655, 22034, 52910, 4880, 72071, 56000, 67442, 49372, 13302, 70567, 37415, 81973, 40820, 18637, 88657, 51472, 59723, 52910, 4880, 72071]
        
        # add and remove characters from the string until we find a regex which fails
        compiled_regex = re.compile(a_regex)
        if len(a_match.matches) > 0 :
            new_search_string = 'a'+ a_match.search_string # check string with one character added at the front
            self.add_test_case_if_fails(new_search_string, compiled_regex)
            new_search_string = a_match.search_string+'a' # check string with one character added at the end
            self.add_test_case_if_fails(new_search_string, compiled_regex)

            for i in range(len(a_match.search_string)-1):
                new_search_string = a_match.search_string[i:] # TODO: refactor this
                self.add_test_case_if_fails(new_search_string, compiled_regex)
            for i in range(len(a_match.search_string)-1):
                new_search_string = a_match.search_string[:i] # TODO: refactor this
                self.add_test_case_if_fails(new_search_string, compiled_regex)


    def add_test_case_if_fails(self,new_search_string, compiled_regex):
        a_test_case_string = RegexTestString(new_search_string)
        vals = self.time_regex_test_case(compiled_regex, a_test_case_string, 1)
        if len(list(vals[1])) == 0:
            self.test_cases.append(a_test_case_string)


    def add_test(self,regex_string, match):
        a_test_string = RegexTestString(regex_string)
        if match and len(match) == 2: # should be e.g. [0,27]
            a_test_string.add_match(match[0],match[1])
        self.test_cases.append(a_test_string)
        return a_test_string
            
    """
    Multiple search_strings should be used to guide toward generality.
    """
    def generate_iso8601_datetime_tests(self):
        # target: ^\d{4,}-[01]\d-[0-3]\dT[0-2]\d:[0-5]\d:[0-5]\d\.\d+(?:[+-][0-2]\d:[0-5]\d|Z)$
        a_test_string = self.add_test("2016-12-09T08:21:15.9+00:00",[0,27])
        self.generate_equivalence_test_suite_replacement(a_test_string,"^\d{4,}-[01]\d-[0-3]\dT[0-2]\d:[0-5]\d:[0-5]\d\.\d+(?:[+-][0-2]\d:[0-5]\d|Z)$")

        # this does not match at all! (what will our fitness function throw?)
        self.add_test("2016-12-09T08:21:15.9+00:0", None)
        self.add_test("2016-22-09T08:21:15.9+00:00000000000",None) # 
        self.add_test("2016-22-09T08:21:15.9+00:00",None) # no match, as 22 is not a valid month
        self.add_test("1911-02-19T22:35:42.3+08:43",[0,27])
        self.add_test("2016-09-05T15:22:26.286Z",[0,24])

    """
    From: https://github.com/angular/angular.js/blob/a24777a2c4ad2ac087d9e3aa278fa2e61e8cc740/src/ng/directive/input.js
    """
    def generate_scientific_number_tests(self):
        self.add_test("230.234E-10", [0,10])
        self.add_test("971.829E+26", [0,10])
        self.add_test("3566", [0,3])
        self.add_test("4", [0,0])
        self.add_test("-7", [0,1])
        self.add_test("+94", [0,2])
        self.add_test("            36", [0,13])
        self.add_test("78      ", [0,8])
        self.add_test("87465.345345", [0,11])
        self.add_test("2346.533", [0,7])
        a_test_string = self.add_test("0.045e-10", [0,8])
        known_regex = "^\s*(-|\+)?(\d+|(\d*(\.\d*)))([eE][+-]?\d+)?\s*$"
        self.generate_equivalence_test_suite_length(a_test_string, known_regex)
        self.generate_equivalence_test_suite_replacement(a_test_string, known_regex)
        # self.add_test("3566.", None) # this is allowed per the original regex 
        # self.add_test(".3456", None) # this is allowed per the original regex 

    """
    From: https://github.com/jmmcd/PonyGE2/blob/c256f9a36331078b9ca298af4d73034b623dd8a0/src/representation/grammar.py
    """
    def generate_PonyGE2_Grammar_File_Rule_tests(self):
        known_regex = "(?P<rulename><\S+>)\s*::=\s*(?P<production>(?:(?=\#)\#[^\r\n]*|(?!<\S+>\s*::=).+?)+)"
        a_test_string = self.add_test("<string> ::= <letter>|<letter><string>", [0,37])
        self.generate_equivalence_test_suite_length(a_test_string, known_regex)
        self.generate_equivalence_test_suite_replacement(a_test_string, known_regex)

    """
    From: https://github.com/jmmcd/PonyGE2/blob/c256f9a36331078b9ca298af4d73034b623dd8a0/src/representation/grammar.py
    """
    def generate_catastrophic_QT3TS_tests(self):
        known_regex = ".X(.+)+XX"
        self.add_test("hryxioXcXXdornct", [5,9])
        a_test_string = self.add_test("bbbbXcyXXaaa", [3,8])
        self.generate_equivalence_test_suite_length(a_test_string, known_regex)
        self.generate_equivalence_test_suite_replacement(a_test_string, known_regex)

    """
    From: https://github.com/d3/d3/commit/4ffe6b40f367c696107b9c01a67812559159633a
    (This is similar to generate_scientific_number_tests, though used for extraction, as opposed to validation)
    """
    def generate_d3_interpolate_number(self):    
        known_regex = "[-+]?(?:\d+\.?\d*|\d*\.?\d+)(?:[eE][-+]?\d+)?"
        self.add_test("230.234E-10", [0,10])
        self.add_test("971.829E+26", [0,10])
        self.add_test("3566", [0,3])
        self.add_test("4", [0,0])
        self.add_test("-7", [0,1])
        self.add_test("+94", [0,2])
        self.add_test("            36", [12,13])
        self.add_test("78      ", [0,1])
        self.add_test("87465.345345", [0,11])
        self.add_test("2346.533", [0,7])
        self.add_test("  3566.   ", [2,6]) # this is allowed per the original regex 
        self.add_test(" .3456  ", [1,5]) # this is allowed per the original regex
        self.add_test("a46b  ", [1,2])
        a_test_string = self.add_test("0.045e-10", [0,8])
        self.generate_equivalence_test_suite_length(a_test_string, known_regex)
        self.generate_equivalence_test_suite_replacement(a_test_string, known_regex)

    """
    From https://github.com/ghiscoding/angular-validation/wiki/Regular-Expression-Pattern
    """        
    def generate_macaddress_validation_tests(self):
        a_test_string = RegexTestString("5C0A5B634A82")
        a_test_string.add_match(0,12)
        self.test_cases.append(a_test_string)

        self.generate_equivalence_test_suite_length(a_test_string, "^[0-9A-F]{12}$")
        self.generate_equivalence_test_suite_replacement(a_test_string, "^[0-9A-F]{12}$")

        
    def generate_regex_mac_search_string_tests(self):
        #a_test_string = RegexTestString("Jan 12 06:26:19: ACCEPT service http from 119.63.193.196 to firewall(pub-nic), prefix: \"none\" (in: eth0 119.63.193.196(5c:0a:5b:63:4a:82):4399 -> 14") 
       #a_test_string.add_match(119,136) # 5c:0a:5b:63:4a:82

        a_test_string = RegexTestString("Jan 12 06:26:19: ACCEPT service http from 119.63.193.196 to firewall(pub-nic), prefix: \"none\" (in: eth0 119.63.193.196(5c:0a:5b:63:4a:82):4399 -> 140.105.63.164(50:06:04:92:53:44):80 TCP flags: ****S* len:60 ttl")
        a_test_string.add_match(119,136) # 5c:0a:5b:63:4a:82
        a_test_string.add_match(161,178) # 50:06:04:92:53:44
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("26:19: ACCEPT service http from 119.63.193.196 to firewall(pub-nic), prefix: \"none\" (in: eth0 119.63.193.196(5c:0a:5b:63:4a:82):4399 -> 140.105.63.1sdkfjhakljwesjdhfglksjhdfgk")
        a_test_string.add_match(109,126) # 5c:0a:5b:63:4a:82
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString(" -> 140.105.63.164(50:j6:04:92:53:44):80 TCP flags: ****S* len:60 ttl:32)sdkfjhaklsjdhfglksjhdfgk")
        # a_test_string.add_match(19,36) # 50:06:04:92:53:44
        # self.test_cases.append(a_test_string)

        # Longer is better for appropriate testing of regex time. (Why is a longer test better than more iterations on a smaller test?) Smaller tests introduce more variability, but why?
        a_test_string = RegexTestString(" -> 140.105.63.16(50:06:04:9r:53:44):80 TCP flags: ****S* len:60 ttl:32)ssjdhfglksjhdfgk")
        #a_test_string.add_match(18,35) # 50:06:04:92:53:44
        # self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("Jan 12 06:26:20: ACCEPT service dns from 140.105.48.16 to firewall(pub-nic-dns), prefix: \"none\" (in: eth0 140.105.48.16(00:21:dd:bc:95:44):4263 -> 140.105.63.158(00:14:31:83:c6:8d):53 UDP len:76 ")
        a_test_string.add_match(120,137)
        a_test_string.add_match(162,179)
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("Jan 12 06:27:09: DROP service 68->67(udp) from 216.34.211.83 to 216.34.253.94, prefix: \"spoof iana-0/8\" (in: eth0 213.92.153.78(00:1f:d6:19:0a:80):68 -> 69.43.177.110(00:30:fe:fd:d6:51):67 UDP le")
        a_test_string.add_match(128,145)
        a_test_string.add_match(167,184)
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("105.63.1650:06:04:92:53:44:80")
        a_test_string.add_match(7,27) # 50:06:04:92:53:44
        # self.test_cases.append(a_test_string)
        self.generate_equivalence_test_suite_replacement(a_test_string, "([0-9A-Fa-f]{2}[:-]){5}([0-9A-Fa-f]{2})")
        
        # a_test_string = RegexTestString(" -> 140.105.63.164(50:06:g4:92:53:44):80 TCP flags: ****S* len:60 ttl:32)") # negative match!
        # self.test_cases.append(a_test_string)
        # a_test_string = RegexTestString(" -> 140.105.63.164(50:06:54:92:r3:44):80 TCP flags: ****S* len:60 ttl:32)") # negative match!
        # self.test_cases.append(a_test_string)


    def generate_catastrophic_csv(self):
        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,12,13777,5P,5,5,6,5P") 
        self.test_cases.append(a_test_string)

        # http://www.regular-expressions.info/catastrophic.html
        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,12,13777,24,5P")
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,12,13777,243,3P") # catastrophic
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,12,13777,P")
        a_test_string.add_match(0,33)
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,12,P")
        a_test_string.add_match(0,27)
        self.test_cases.append(a_test_string)
        
        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,P")
        a_test_string.add_match(0,24)
        self.test_cases.append(a_test_string)

        #        self.generate_equivalence_test_suite_replacement(a_test_string,"^(.*?,){11}P")
        self.generate_equivalence_test_suite_length(a_test_string,"^(.*?,){11}P")

        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,11,3P")
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("1,2,3,4,5,6,7,8,9,10,P")
        self.test_cases.append(a_test_string)

    def generate_email_validation_tests(self):
        a_test_string = RegexTestString("codykenny@gmail.com")
        a_test_string.add_match(0,18)
        self.test_cases.append(a_test_string)

        self.generate_equivalence_test_suite_replacement(a_test_string,"^(?=.{1,254}$)(?=.{1,64}@)[-!#$%&'*+/0-9=?A-Z^_`a-z{|}~]+(\.[-!#$%&'*+/0-9=?A-Z^_`a-z{|}~]+)*@[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?(\.[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?)*$")
        self.generate_equivalence_test_suite_length(a_test_string,"^(?=.{1,254}$)(?=.{1,64}@)[-!#$%&'*+/0-9=?A-Z^_`a-z{|}~]+(\.[-!#$%&'*+/0-9=?A-Z^_`a-z{|}~]+)*@[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?(\.[A-Za-z0-9]([A-Za-z0-9-]{0,61}[A-Za-z0-9])?)*$")

        a_test_string = RegexTestString("codykennygmail.com")
        self.test_cases.append(a_test_string)

        a_test_string = RegexTestString("codykenny@gmailcom")
        a_test_string.add_match(0,17)
        self.test_cases.append(a_test_string)

        
    def generate_tests(self):
        test_gen_method = getattr(self, params['FITNESS_TEST_SUITE'])
        test_gen_method()
        # can be any of these:
        # self.generate_regex_mac_search_string_tests()
        # self.generate_catastrophic_csv()
        # self.generate_iso8601_datetime_tests()
        # self.generate_macaddress_validation_tests()
        # self.generate_email_validation_tests()
        # self.generate_scientific_number_tests()
        # self.generate_PonyGE2_Grammar_File_Rule_tests()
        # self.generate_catastrophic_QT3TS_tests()
        # self.generate_d3_interpolate_number()

        print("Number of test cases: {}".format(len(self.test_cases)))

    def __call__(self, individual):
        q = Queue()
        p = Process(target=self.call_fitness, name="self.call_fitness", args=(individual,q))
        p.start()

        p.join(5)
        
        # If thread is active
        if p.is_alive():
            print("      ------------      ------------      ------------    timeout reached, killing process")
            # Terminate foo
            p.terminate()
            p.join()
            return 100000
        else:
            fitness = q.get()
#            print("\n again: {}".format(fitness))            
            return fitness

            

class RegexTestString:
    def __init__(self,search_string):
        # print("Added regex search string: "+search_string)
        self.search_string = search_string
        self.matches = list()

    def add_match(self, start, end):
        self.matches.append({'start': start,'end': end, 'matched_string': self.search_string[start:end]}) # ew?
        # print("Added the following match: "+self.matches[-1].get("matched_string"))
#         for i in range(start,end-1):
#             for j in range(start+1,end):
#                 self.match.append(start=i,end=j, matched_string=self.search_string[i:j])
        
    def calc_match_errors(self,match_candidates):
        match_ranges=list()
        undesired_range = missing_range = 0
        #        for match in match_candidates:
        #           match_ranges.append(match)

        for a_known_match in self.matches:
            # missing any of the desired extraction costs a lot
            missing_range += self.find_missing_range(a_known_match, match_candidates)
        for match_candidate in match_candidates:
            undesired_range += self.find_undesired_range(match_candidate, self.matches)

        match_error = missing_range + undesired_range
        match_error += (abs(len(match_candidates) - len(self.matches)))

        return (match_error) 

    def find_missing_range(self, a_known_match, match_ranges):
        start = a_known_match.get("start")
        end = a_known_match.get("end")
        missing = end - start
        for i in range(start,end):
            found = False
            for m_range in match_ranges:
                if i >= m_range.start() and i < m_range.end():
                    found = True
            if found:
                missing -= 1
        return missing

    def find_undesired_range(self,match_candidate, known_matches):
        undesired_matched = 0
        for i in range(match_candidate.start(), match_candidate.end()):
            in_range=False
            for a_known_match in known_matches:
                start = a_known_match.get("start")
                end = a_known_match.get("end")
                if i >= start and i <= end:
                    in_range = True
            if not in_range:
                undesired_matched += 1
        return undesired_matched


    #             for a_match_candidate in match_candidates: # how much of our known match has gone unmatched?
#                 temp_match=0
#                 if a_matched_candidate.start < a_known_match.end or a_matched_candidate.end > a_known_match.start:
#                     unmatched_start = a_matched_candidate.start - a_known_match.start # should be positive, if negative then there is none unmatched
#                     unmatched_end = a_matched_candidate.end - a_known_match.start # should be positive
#                     temp_match = max(unmatched_end,0) + max(unmatched_start,0)


        
#     def compare(self,attempt_string):
#         error = distance(self.matched_string, attempt_string.group(0))
#         attempti=0
#         attempti += self.search_string.find(attempt_string.group(0))
#         if attempti < self.starti:
#             error+= self.starti-attempti
#         elif attempti > self.endi:
#             error+= attempti-self.endi
#         return error

    def get_search_string(self):
        return self.search_string

    