from utilities.helper_methods import get_Xy_train_test_separate
import subprocess
import json

class ProgSys:
    """"""

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

    def __init__(self, alternate):
        self.training, self.test, self.embed_header, self.embed_footer = get_data(alternate)
        self.eval = subprocess.Popen(['python', 'fitness/python_script_evaluation.py'],
                                stdout=subprocess.PIPE,
                                stdin=subprocess.PIPE)

        self.count = 0

    def __call__(self, individual):  # , dist):
        program = self.format_program(individual, self.embed_header, self.embed_footer)

        program = self.training + '\n' + program + '\n'
        eval_json = json.dumps({'script' : program, 'timeout' : 1.0, 'variables' : ['cases', 'caseQuality', 'quality']})

        self.eval.stdin.write((eval_json+'\n').encode())
        self.eval.stdin.flush()

        result_json = self.eval.stdout.readline()

        result = json.loads(result_json.decode())
        # if dist == "test":
        #     x = self.test_in
        #     y = self.test_exp
        # elif dist == "training":
        #     x = self.training_in
        #     y = self.training_exp

        # self.count += 1
        # if self.count == 10 and result['quality'] == 1559.046069859538:
        #     print(individual)
        #     print(result['quality'])

        # if self.count == 10 and result['quality'] == 1696.7149730636077:
        # print(individual)
        # print(result['quality'])
        return result['quality']

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
                line = line.Substring[0:len(line) - 2].strip()

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


INSERTCODE = "<insertCodeHere>"


def get_data(alternate):
    train_set = "datasets/" + alternate + "-Train.txt"
    test_set = "datasets/" + alternate + "-Test.txt"
    embed_file = "grammars/" + alternate + "-Embed.txt"
    embed_header, embed_footer = "", ""
    with open(embed_file, 'r') as embed:
        embedCode = embed.read()
        insert = embedCode.index(INSERTCODE)
        if insert > 0:
            embed_header = embedCode[:insert]
            embed_footer = embedCode[insert + len(INSERTCODE):]
    with open(train_set, 'r') as train, open(test_set, 'r') as test:
        return train.read(), test.read(), embed_header, embed_footer
