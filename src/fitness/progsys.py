from algorithm.parameters import params
from os import listdir, getcwd, path
import subprocess
import json
import sys

class progsys:
    """"""

    INSERTCODE = "<insertCodeHere>"

    INDENTSPACE = "  "
    LOOPBREAK = "loopBreak"
    LOOPBREAKUNNUMBERED = "loopBreak%"
    LOOPBREAK_INITIALISE = "loopBreak% = 0"
    LOOPBREAK_IF = "if loopBreak% >"
    LOOPBREAK_INCREMENT = "loopBreak% += 1"
    FUNCTIONSTART = "def evolve():"
    FORCOUNTER = "forCounter"
    FORCOUNTERUNNUMBERED = "forCounter%"

    maximise = False

    def __init__(self):
        self.training, self.test, self.embed_header, self.embed_footer = self.get_data(params['DATASET'])
        self.create_eval_process()
        if params['MULTICORE']:
            print("Warming: Multicore is not supported with progsys as fitness function.\n"
                  "Fitness function only allows sequential evaluation.")

    def __call__(self, individual, dist="training"):
        program = self.format_program(individual.phenotype, self.embed_header, self.embed_footer)
        data = self.training if dist == "training" else self.test
        program = "{}\n{}\n".format(data, program)
        eval_json = json.dumps({'script': program, 'timeout': 1.0, 'variables': ['cases', 'caseQuality', 'quality']})

        self.eval.stdin.write((eval_json+'\n').encode())
        self.eval.stdin.flush()
        result_json = self.eval.stdout.readline()

        result = json.loads(result_json.decode())

        if 'exception' in result and 'JSONDecodeError' in result['exception']:
            self.eval.stdin.close()
            self.create_eval_process()

        if 'quality' not in result:
            result['quality'] = sys.maxsize
        return result['quality']

    def create_eval_process(self):
        self.eval = subprocess.Popen(['python', 'fitness/python_script_evaluation.py'],
                                stdout=subprocess.PIPE,
                                stdin=subprocess.PIPE)

    def format_program(self, individual, header, footer):
        lastNewLine = header.rindex('\n')
        indent = header[lastNewLine + len('\n'):len(header)]
        return header + self.format_individual(individual, indent) + footer

    def format_individual(self, code, additional_indent=""):
        parts = code.split('\n')
        indent = 0
        stringBuilder = ""
        forCounterNumber = 0
        first = True
        for part in parts:
            line = part.strip()
            # remove indentation if bracket is at the beginning of the line
            while line.startswith(":}"):
                indent -= 1
                line = line[2:len(line) - 2].strip()

            # add indent
            if not first:
                stringBuilder += additional_indent
            else:
                first = False

            for i in range(0, indent):
                stringBuilder += self.INDENTSPACE

            # add indentation
            while line.endswith("{:"):
                indent += 1
                line = line[0:len(line) - 2].strip()
            # remove indentation if bracket is at the end of the line
            while line.endswith(":}"):
                indent -= 1
                line = line[0:len(line) - 2].strip()

            if self.LOOPBREAKUNNUMBERED in line:
                if self.LOOPBREAK_INITIALISE in line:
                    line = ""  # remove line
                elif self.LOOPBREAK_IF in line:
                    line = line.replace(self.LOOPBREAKUNNUMBERED, self.LOOPBREAK)
                elif self.LOOPBREAK_INCREMENT in line:
                    line = line.replace(self.LOOPBREAKUNNUMBERED, self.LOOPBREAK)
                else:
                    raise Exception("Python 'while break' is malformed.")
            elif self.FORCOUNTERUNNUMBERED in line:
                line = line.replace(self.FORCOUNTERUNNUMBERED, self.FORCOUNTER + str(forCounterNumber))
                forCounterNumber += 1

            # add line to code
            stringBuilder += line
            stringBuilder += '\n'

        return stringBuilder

    def get_data(self, experiment):
        """ Return the training and test data for the current experiment.
        """
        train_set = path.join("..", "datasets", "progsys",
                              (experiment + "-Train.txt"))
        test_set = path.join("..", "datasets", "progsys",
                             (experiment + "-Test.txt"))

        embed_file = path.join("..", "grammars", "progsys", (experiment + "-Embed.txt"))
        with open(embed_file, 'r') as embed:
            embed_code = embed.read()
        insert = embed_code.index(self.INSERTCODE)
        embed_header, embed_footer = "", ""
        if insert > 0:
            embed_header = embed_code[:insert]
            embed_footer = embed_code[insert + len(self.INSERTCODE):]
        with open(train_set, 'r') as train_file, open(test_set, 'r') as test_file:
            return train_file.read(), test_file.read(), embed_header, embed_footer
