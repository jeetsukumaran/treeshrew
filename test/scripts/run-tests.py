#! /usr/bin/env python

"""
Run tests.
"""

import sys
import os
import argparse
import inspect
import subprocess
import re

import dendropy
from dendropy import treecalc
from dendropy.interop import paup

class AnsiColorMeta(type):

    ##############################################################################
    ## Color infrastructure modified from:
    ##
    ## http://dev.pocoo.org/hg/pygments-main/file/b2deea5b5030/pygments/console.py
    ##
    ## pygments.console
    ## ~~~~~~~~~~~~~~~~
    ## Format colored console output.
    ## :copyright: Copyright 2006-2009 by the Pygments team, see AUTHORS.
    ## :license: BSD, see LICENSE for details.
    ##

    ansiesc = "\x1b["

    @staticmethod
    def get_ansicodes():
        ansicodes = {}
        ansicodes[""]          = ""
        ansicodes["reset"]     = AnsiColorMeta.ansiesc + "39;49;00m"
        ansicodes["bold"]      = AnsiColorMeta.ansiesc + "01m"
        ansicodes["faint"]     = AnsiColorMeta.ansiesc + "02m"
        ansicodes["standout"]  = AnsiColorMeta.ansiesc + "03m"
        ansicodes["underline"] = AnsiColorMeta.ansiesc + "04m"
        ansicodes["blink"]     = AnsiColorMeta.ansiesc + "05m"
        ansicodes["overline"]  = AnsiColorMeta.ansiesc + "06m"
        dark_colors  = ["black", "darkred", "darkgreen", "brown", "darkblue",
                        "purple", "teal", "lightgray"]
        light_colors = ["darkgray", "red", "green", "yellow", "blue",
                        "fuchsia", "turquoise", "white"]
        x = 30
        for d, l in zip(dark_colors, light_colors):
            ansicodes[d] = AnsiColorMeta.ansiesc + "%im" % x
            ansicodes[l] = AnsiColorMeta.ansiesc + "%i;01m" % x
            x += 1
        ansicodes["darkteal"]   = ansicodes["turquoise"]
        ansicodes["darkyellow"] = ansicodes["brown"]
        ansicodes["fuscia"]     = ansicodes["fuchsia"]
        # ansicodes["white"]      = ansicodes["bold"]
        return ansicodes

    def reset_color(cls):
        return cls.ansicodes["reset"]

    def colorize(cls, color_key, text):
        return cls.ansicodes[color_key] + text + cls.ansicodes["reset"]

    def ansiformat(cls, attr, text):
        """
        Format ``text`` with a color and/or some attributes::

            color       normal color
            *color*     bold color
            _color_     underlined color
            +color+     blinking color
        """
        result = []
        if attr[:1] == attr[-1:] == '+':
            result.append(cls.ansicodes['blink'])
            attr = attr[1:-1]
        if attr[:1] == attr[-1:] == '*':
            result.append(cls.ansicodes['bold'])
            attr = attr[1:-1]
        if attr[:1] == attr[-1:] == '_':
            result.append(cls.ansicodes['underline'])
            attr = attr[1:-1]
        result.append(cls.ansicodes[attr])
        result.append(text)
        result.append(cls.ansicodes['reset'])
        return ''.join(result)

    def __new__(cls, name, bases, dct):
        return type.__new__(cls, name, bases, dct)

    def __init__(cls, name, bases, dct):
        super(AnsiColorMeta, cls).__init__(name, bases, dct)
        # setattr(cls, "ansicodes", AnsiColorMeta.get_ansicodes())
        # setattr(cls, "ansiesc", AnsiColorMeta.ansiesc)
        cls.ansicodes = AnsiColorMeta.get_ansicodes()
        cls.ansiesc = AnsiColorMeta.ansiesc

class AnsiColor(object):

    __metaclass__ = AnsiColorMeta

    def __init__(self, stream=sys.stdout, colorize=True):
        self.stream = stream
        self.colorize = colorize
        self.color_pattern = re.compile(r"@(\w+)@<<(.*?)>>")

    def format_color(self, message):
        if not self.colorize or not self.color_pattern.findall(message):
            return message
        else:
            output = []
            cur_pos = 0
            for match in self.color_pattern.finditer(message):
                start, end = match.span()
                output.append(message[cur_pos:start])
                output.append(AnsiColor.ansiformat(match.group(1), match.group(2)))
                cur_pos = end
            output.append(message[cur_pos:])
            output = "".join(output)
            return output

    def write(self, message):
        self.stream.write(self.format_color(message))

    def __call__(self, message):
        self.write(message)

class TestRunner(object):

    PASS = 0
    FAIL = 1
    ERROR = 2

    dna_to_partials =  {
        'A' : [1.0, 0.0, 0.0, 0.0],
        'C' : [0.0, 1.0, 0.0, 0.0],
        'G' : [0.0, 0.0, 1.0, 0.0],
        'T' : [0.0, 0.0, 0.0, 1.0],
        'U' : [0.0, 0.0, 0.0, 1.0],
        'N' : [1.0, 1.0, 1.0, 1.0],
        'X' : [1.0, 1.0, 1.0, 1.0],
        '-' : [1.0, 1.0, 1.0, 1.0],
        '?' : [1.0, 1.0, 1.0, 1.0],
        'R' : [1.0, 0.0, 1.0, 0.0],
        'Y' : [0.0, 1.0, 0.0, 1.0],
        'M' : [1.0, 1.0, 0.0, 0.0],
        'W' : [1.0, 0.0, 0.0, 1.0],
        'S' : [0.0, 1.0, 1.0, 0.0],
        'K' : [0.0, 0.0, 1.0, 1.0],
        'V' : [1.0, 1.0, 1.0, 0.0],
        'H' : [1.0, 1.0, 0.0, 1.0],
        'D' : [1.0, 0.0, 1.0, 1.0],
        'B' : [0.0, 1.0, 1.0, 1.0]
    }

    @staticmethod
    def get_node_tag(nd):
        if nd.taxon:
            return nd.taxon.label
        else:
            return nd.label

    def __init__(self, opts):
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.data_dir = os.path.join(self.script_dir, "data")
        self.verbosity = opts.verbosity
        self.break_on_fail = opts.break_on_fail
        self.test_command = None
        self.test_retcode = None
        self.test_stdout = None
        self.test_stderr = None
        self.test_result = None
        self.test_fail_message = None
        self.test_pass_message = None
        self.cout = AnsiColor(sys.stdout, colorize=True)

    def execute_test(self, test_program, args=None):
        self.test_command = None
        self.test_stdout = None
        self.test_stderr = None
        self.test_retcode = None
        self.test_fail_message = None
        cmd = [os.path.abspath(os.path.join(self.script_dir, test_program))]
        if args:
            cmd.extend(args)
        self.test_command = " ".join(cmd)
        p = subprocess.Popen(cmd,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                cwd=self.script_dir)
        self.test_stdout, self.test_stderr = p.communicate()
        self.test_retcode = p.returncode
        return self.test_retcode

    def fail(self, message):
        self.test_fail_message = message
        return TestRunner.FAIL

    def is_almost_equal(self, v1, v2, prec=1e-4):
        if v1 is None and v2 is None:
            return True
        if v1 is None:
            return self.is_almost_equal(v2, 0.0)
        if v2 is None:
            return self.is_almost_equal(v1, 0.0)
        if abs(v1-v2) > prec:
            return False
        return True

    def compare_trees(self, tree1, tree2):
        status = TestRunner.PASS
        tree1.update_splits()
        tree2.update_splits()

        splits = set(tree1.split_edges.keys() + tree2.split_edges.keys())
        for split in splits:
            if split not in tree2.split_edges:
                return self.fail("Split {} not found on tree 2: {}".format(split, tree1.taxon_set.split_as_newick_string(split)))
            if split not in tree1.split_edges:
                return self.fail("Split {} not found on tree 1: {}".format(split, tree2.taxon_set.split_as_newick_string(split)))
            edge1_len = tree1.split_edges[split].length
            edge2_len = tree2.split_edges[split].length
            if not self.is_almost_equal(edge1_len, edge2_len):
                return self.fail("Unequal edge length for split {} {}: {} vs. {}".format(
                    split,
                    tree1.taxon_set.split_as_newick_string(split),
                    edge1_len,
                    edge2_len))
        if len(tree1.split_edges) != len(tree2.split_edges):
            return self.fail("Different number of splits on trees: {} vs. {}".format(len(tree1.split_edges), len(tree2.split_edges)))
        return status

    def compare_tree_traversal(self, tree1, tree2, traverse_func):
        status = TestRunner.PASS
        tree1_nodes = [nd for nd in getattr(tree1, traverse_func)()]
        tree2_nodes = [nd for nd in getattr(tree2, traverse_func)()]
        if len(tree1_nodes) != len(tree2_nodes):
            return self.fail("Trees have different number of nodes: {} vs. {}".format(len(tree1_nodes), len(tree2_nodes)))
        for nd_idx, node1 in enumerate(tree1_nodes):
            node2 = tree2_nodes[nd_idx]
            if node1.taxon is not node2.taxon:
                return self.fail("Different taxa found during postorder traversal of nodes: {} vs. {}".format(node1.taxon, node2.taxon))
            if node1.label is not node2.label:
                return self.fail("Different labels found during postorder traversal of nodes: {} vs. {}".format(node1.label, node2.label))
            if not self.is_almost_equal(node1.edge.length, node2.edge.length):
                return self.fail("Different edge lengths found during postorder traversal of nodes: {} vs. {}".format(node1.edge.length, node2.edge.length))
        return status

    def test_tree_postorder_iter(self):
        treefile = os.path.join(self.data_dir, "basic", "pythonidae.postorder.newick")
        self.execute_test("tree_postorder_iter",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_visits = []
        test_edge_lens = []
        for item in self.test_stdout.split("\n"):
            if not item:
                continue
            node, edge_len = item.split("\t")
            test_visits.append(node)
            test_edge_lens.append(float(edge_len))
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", taxon_set=taxa)
        check_visits = []
        check_edge_lens = []
        for check_node in check_tree.postorder_node_iter():
            label = self.get_node_tag(check_node)
            edge_len = check_node.edge.length if check_node.edge.length else 0.0
            check_visits.append(label)
            check_edge_lens.append(edge_len)
        for idx in range(len(check_visits)):
            if idx > len(test_visits):
                return self.fail("Insufficient visits: expecting {} but found {}".format(len(check_visits), len(test_visits)))
            n1 = check_visits[idx]
            n2 = test_visits[idx]
            if n1 != n2:
                return self.fail("Incorrect visit {}: '{}' vs. '{}'".format(idx+1, n1, n2))
            e1 = check_edge_lens[idx]
            e2 = test_edge_lens[idx]
            if not self.is_almost_equal(e1, e2):
                return self.fail("Incorrect node edge length on visit {}: {} vs. {}".format(idx+1, e1, e2))
        return TestRunner.PASS

    def test_tree_leaf_iter(self):
        treefile = os.path.join(self.data_dir, "basic", "pythonidae.postorder.newick")
        self.execute_test("tree_leaf_iter",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_leaves = []
        for item in self.test_stdout.split("\n"):
            if not item:
                continue
            idx, label = item.split("\t")
            test_leaves.append(label)
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", taxon_set=taxa)
        check_leaves = [ self.get_node_tag(nd) for nd in check_tree.leaf_iter() ]
        test_leaves.sort()
        check_leaves.sort()
        if set(check_leaves) != set(test_leaves):
            return self.fail("Unequal leaf set: {} vs. {}".format(set(check_leaves), set(test_leaves)))
        if len(check_leaves) != len(test_leaves):
            return self.fail("Duplicate leaves: {}".format(set(test_leaves)))
        return TestRunner.PASS

    def test_tree_child_iter(self):
        treefile = os.path.join(self.data_dir, "basic", "pythonidae.postorder.newick")
        self.execute_test("tree_child_iter",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_node_children = {}
        for item in self.test_stdout.split("\n"):
            if not item:
                continue
            nodes = item.split("\t")
            test_node_children[nodes[0]] = nodes[1:]
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", taxon_set=taxa)
        check_node_children = {}
        for check_node in check_tree.postorder_node_iter():
            check_node_children[self.get_node_tag(check_node)] = [self.get_node_tag(child) for child in check_node.child_nodes()]
        for check_node, check_children in check_node_children.items():
            if check_node not in test_node_children:
                return self.fail("Node not visited: '{}'".format(check_node))
            if test_node_children[check_node] != check_children:
                return self.fail("Incorrect children: '{}' vs. '{}'".format(check_children, test_node_children[check_node]))
        return TestRunner.PASS

    def compare_tree_scores(self, tree_filename, data_filename):
        full_tree_filepath = os.path.join(self.data_dir, "basic", tree_filename)
        full_data_filepath = os.path.join(self.data_dir, "basic", data_filename)
        self.execute_test("score_phylogenetic_tree",
                [full_tree_filepath, full_data_filepath])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_ln_like = float(self.test_stdout)
        taxa = dendropy.TaxonSet()
        tree = dendropy.Tree.get_from_path(full_tree_filepath, "newick", taxon_set=taxa)
        data = dendropy.DnaCharacterMatrix.get_from_path(full_data_filepath, "fasta", taxon_set=taxa)
        paup_tree, results = paup.estimate_model(
                char_matrix=data,
                tree_model=tree,
                num_states=1,
                unequal_base_freqs=False,
                gamma_rates=False,
                prop_invar=False,
                tree_est_criterion="likelihood",
                tree_user_brlens=True,
                paup_path="paup")
        check_ln_like = float(results['likelihood'])
        if not self.is_almost_equal(check_ln_like, test_ln_like):
            return self.fail("Unequal log-likelihoods (tree: '{}', data: '{}'): {} vs. {}".format(tree_filename, data_filename, check_ln_like, test_ln_like))
        return TestRunner.PASS

    def test_tree_score1(self):
        return self.compare_tree_scores("primates.beast.mcct.medianh.newick.tre", "primates.chars.fasta")

    def test_tree_score2(self):
        return self.compare_tree_scores("pythonidae.tree.newick", "pythonidae.chars.fasta")

    def test_tree_read_from_file(self):
        treefile = os.path.join(self.data_dir, "basic", "bird_orders.nex")
        self.execute_test("read_tree",
                [treefile, "nexus"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "nexus", taxon_set=taxa)
        test_tree = dendropy.Tree.get_from_string(self.test_stdout, "newick", taxon_set=taxa)
        # return self.compare_tree_traversal(check_tree, test_tree, "postorder_node_iter")
        return self.compare_trees(check_tree, test_tree)

    def test_read_dna_sequences(self):
        datafile = os.path.join(self.data_dir, "basic", "pythonidae.chars.fasta")
        self.execute_test("read_dna_sequences",
                [datafile, "fasta"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        dna_matrix = dendropy.DnaCharacterMatrix.get_from_path(datafile, "fasta", row_type="STR")
        expected_partials = {}
        for taxon in dna_matrix:
            # label = taxon.label.replace(" ", "_")
            label = taxon.label
            expected_partials[label] = []
            states = dna_matrix[taxon]
            for state in states:
                sub_partials = self.dna_to_partials[state]
                expected_partials[label].extend(list(sub_partials))
        rows = self.test_stdout.split("\n")
        observed_partials = {}
        for row in rows:
            if not row:
                continue
            label, partials = row.split(":")
            partials = [float(v) for v in partials.split(";") if v]
            observed_partials[label] = list(partials)
        for label in expected_partials:
            if label not in observed_partials:
                return self.fail("Sequence '{}' not found: {}".format(label,
                    ",".join(["'{}'".format(t) for t in observed_partials])))
            p1 = expected_partials[label]
            p2 = observed_partials[label]
            if len(p1) != len(p2):
                return self.fail("Sequence '{}': expected {} elements but found {}".format(label,
                    len(p1), len(p2)))
            for idx, i1 in enumerate(p1):
                i2 = p2[idx]
                if not self.is_almost_equal(i1, i2):
                    return self.fail("Sequence '{}': character {}: expected {} but found {}".format(label,
                        idx, i1, i2))
        return TestRunner.PASS

    def run(self):
        tests_to_run = []
        for name, value in inspect.getmembers(self, callable):
            if name.startswith("test"):
                tests_to_run.append((name, value))
        passes = []
        fails = []
        for test_idx, (test_name, test_call) in enumerate(tests_to_run):
            self.cout("@turquoise@<<{: 4d}/{:<4d}>>: {}: ".format(test_idx+1, len(tests_to_run), test_name))
            result = test_call()
            if result == TestRunner.PASS:
                self.cout("@green@<<PASS>>\n")
                passes.append(test_name)
                if self.test_pass_message and self.verbosity > 3:
                    self.cout("         :   - {}\n".format(self.test_pass_message))
            elif result == TestRunner.FAIL:
                self.cout("@fuchsia@<<FAIL>>\n")
                fails.append(test_name)
                if self.test_fail_message:
                    self.cout("         :   - {}\n".format(self.test_fail_message))
                if self.break_on_fail:
                    self.summarize(passes, fails)
                    return
            else:
                self.cout("@red@<<ERROR>>\n")
                self.cout("         @red@<<:        Executed:>> {}\n".format(self.test_command))
                self.cout("         @red@<<:     Return Code:>> {}\n".format(self.test_retcode))
                self.cout("         @red@<<: Standard Output:>> {}\n".format(self.test_stdout))
                self.cout("         @red@<<:  Standard Error:>> {}\n".format(self.test_stderr))
                if self.break_on_fail:
                    self.summarize(passes, fails)
                    return
        self.summarize(passes, fails)

    def summarize(self, passes, fails):
        self.cout("\n--\nTests completed.\n")
        self.cout("@turquoise@<<{}>> tests run with @green@<<{}>> successes and @fuchsia@<<{}>> failures.\n".format(len(passes)+len(fails), len(passes), len(fails)))

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbosity",
            default=1,
            help="control messaging level")

    parser.add_argument("-b", "--break-on-fail",
            action="store_true",
            default=False,
            help="terminate tests after first failure")

    args = parser.parse_args()

    test_runner = TestRunner(args)
    test_runner.run()

if __name__ == '__main__':
    main()


