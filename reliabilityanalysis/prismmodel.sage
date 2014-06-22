import os
import glob
import unittest
import sys
import textwrap
from subprocess import call
from itertools import groupby


def generate_ctmc(G, failure_rates):
    """
    Return the digraph and transition rate matrix for the continuous-time
    Markov chain corresponding to the design space graph G, whose elements have
    the failure rates specified in failure_rates.

    """
    # Use dictionary as a workaround for non-existence of 'nonlocal' statement
    # in python 2.x
    d = {
        # digraph to which we will be adding subgraphs as vertices
        'dg': DiGraph(),
        # the last vertex that has been added to the digraph dg
        'last_vertex': -1,
        # dictionary of the coefficients of the transition rate matrix
        'M': {}
    }

    def generate_ctmc_aux(
            # graph to which we apply recursion
            current_graph):

        for v in d['dg'].vertices():
            if current_graph == d['dg'].get_vertex(v):
                # current_graph has already been added to dg.
                return v

        # add current_graph as a new vertex to dg
        next_vertex = d['last_vertex'] + 1
        d['dg'].add_vertex(next_vertex)
        d['dg'].set_vertex(next_vertex, current_graph)
        d['last_vertex'] = next_vertex
        root_vertex = d['last_vertex']

        # find new subgraphs through single edge deletion
        eL = current_graph.edges()
        if eL != []:
            for e in eL:
                graph_copy = copy(current_graph)
                graph_copy.delete_edge(e)
                subgraph = graph_copy
                child_vertex = generate_ctmc_aux(subgraph)
                d['dg'].add_edge(root_vertex, child_vertex)
                d['M'][(root_vertex, child_vertex)] = failure_rates[e[:2]]

        # find new subgraphs through single vertex deletion
        vL = current_graph.vertices()
        if vL != []:
            for v in vL:
                graph_copy = copy(current_graph)
                graph_copy.delete_vertex(v)
                subgraph = graph_copy
                child_vertex = generate_ctmc_aux(subgraph)
                d['dg'].add_edge(root_vertex, child_vertex)
                d['M'][(root_vertex, child_vertex)] = failure_rates[v]
        return root_vertex

    def insert_diagonals(M):
        for row_number in range(M.nrows()):
            row_sum = sum(M.row(row_number))
            M[row_number, row_number] = -row_sum

    generate_ctmc_aux(G)
    matrix_size = d['last_vertex'] + 1
    trans_matrix = Matrix(QQ, matrix_size, d['M'], sparse=False)
    insert_diagonals(trans_matrix)
    return d['dg'], trans_matrix




def create_prism_files(trans_rate_matrix, output_files_basename):
    """
    Generates the explicit PRISM model files corresponding to the continuous
    time Markov chain specified by the transition rate matrix trans_rate_matrix
    and writes it to files named <output_files_basename>.{tra,sta}.

    See http://www.prismmodelchecker.org/manual/Appendices/ExplicitModelFiles
    for the format of an explicit PRISM model file.
    """
    def create_transitions_file(trans_rate_matrix, output_files_basename):
        assert trans_rate_matrix.nrows() == trans_rate_matrix.ncols()

        transition_count = 0
        content = []
        for row_number in range(trans_rate_matrix.nrows()):
            for column_number in range(trans_rate_matrix.ncols()):
                failure_rate = trans_rate_matrix[row_number, column_number]
                if failure_rate > 0:
                    transition_count += 1
                    content.append("{} {} {}".format(row_number, column_number,
                        n(failure_rate)))
        transitions_file_name = output_files_basename + ".tra"
        transitions_file = open(transitions_file_name, "w")
        state_count = trans_rate_matrix.nrows()
        transitions_file.write("{} {}\n".format(state_count, transition_count))
        transitions_file.write("\n".join(content))
        transitions_file.write("\n")
        transitions_file.close()
        return transitions_file_name

    def create_states_file(trans_rate_matrix, output_files_basename):
        states_file_name = output_files_basename + ".sta"
        states_file = open(states_file_name, "w")
        # we identify each state by means of a single variable called 's'
        states_file.write("(s)\n")
        assert trans_rate_matrix.nrows() == trans_rate_matrix.ncols()
        # each row of the matrix corresponds to a state
        for row_number in range(trans_rate_matrix.nrows()):
            states_file.write("{}:({})\n".format(row_number, row_number))
        states_file.close()
        return states_file_name

    transitions_file_name = create_transitions_file(trans_rate_matrix,
        output_files_basename)
    states_file_name = create_states_file(trans_rate_matrix, output_files_basename)
    # the labels file will be created in run_experiments() because the labels
    # file selects the initial state
    return transitions_file_name, states_file_name



def filter_graphs(digraph, filter_function):
    """
    Returns a list of vertices of the digraph such that each vertex of the list
    has associated a graph that passes the filter function.
    """
    filter_passing_vertices = []
    for v in digraph.vertices():
        graph = digraph.get_vertex(v)
        if filter_function(graph) == True:
            filter_passing_vertices.append(v)
    return filter_passing_vertices


def prism_properties(non_faulty_states, output_file_name):
    """
    Generates a PRISM properties file called output_file_name that gives the
    probability that at an instant of time T the PRISM model is in one of the
    states listed in non_faulty_states.
    """
    properties_file = open(output_file_name, "w")
    file_content_template = """
        // This file has automatically been generated by {this_file}.

        // Instant of time to be used to determine transient probabilities
        const int T;

        label "nonfaulty" = {expression};

        // Returns the probability that, starting from the initial state, in
        // the future, at instant of time T, we are in a non-faulty state. In
        // other words, the following property returns the reliability at time
        // T.
        P=? [ F[T,T] "nonfaulty" ];
    """
    expression_strings = []
    for state in non_faulty_states:
        expression_strings.append("(s={})".format(state))
    file_content = file_content_template.format(
        this_file=sys.argv[0],
        expression=('|').join(expression_strings))
    properties_file.write(textwrap.dedent(file_content))


def run_experiments(
        transitions_file_name,
        states_file_name,
        properties_file_name,
        results_dir_name,
        initial_states,
        time):
    """
    Run experiments with the PRISM model checker.

    INPUT:
        - transitions_file_name: the path to the file that contains the
          explicit transitions file for PRISM.
        - states_file_name: the path to the file that contains the
          explicit states file for PRISM.
        - properties_file_name: the path to the file that contains the
          properties to verify with the PRISM model checker.
        - results_dir_name: name of the directory where to store the results.
        - initial_states: list of states that should be chosen as initial
          states of the markov chain specified in model_file_name.
        - time: instant of time up to which to measure transient reliability in
          the experiments.
    """
    for state in initial_states:

        stdout_file_name = "{}/{}.stdout".format(results_dir_name, state)
        stdout_file = open(stdout_file_name, "w")

        stderr_file_name = "{}/{}.stderr".format(results_dir_name, state)
        stderr_file = open(stderr_file_name, "w")

        results_file_name = "{}/{}.rslts".format(results_dir_name, state)

        base_name, file_extension = os.path.splitext(states_file_name)
        labels_file_name = base_name + ".lab"
        labels_file = open(labels_file_name, "w")
        labels_file.write('0="init"\n{}: 0'.format(state))
        labels_file.close()
        call(["prism", "-importtrans", transitions_file_name, "-ctmc",
            properties_file_name,
            "-importstates", states_file_name,
            "-importlabels", labels_file_name,
            "-const", "T=0:{}".format(time),
            "-exportresults", results_file_name],
            stdout=stdout_file, stderr=stderr_file)
        stdout_file.close()
        stderr_file.close()



def plot_results(results_dir_name, time, max_curves=float("inf"),
        decimal_places=8):

    """
    INPUT

        - time: maximum value for the X axis in the plot.

        - max_curves: the maximum number of curves to plot. Note that the
          number of curves may differ from the number of results because
          several results may be plotted by the same curve. This option is
          useful in case that the number of results in results_dir_name leads
          to so many curves that the plot becomes unreadable.

        - decimal_places: number of decimal places to consider in the results.
          Reducing the number of decimal places can ensure that values that
          should be identical are not different due to the inherent inaccuracy
          of the floating-point number results given by Prism.
    """
    def collect_results(results_dir_name):
        """
        """
        plot_points_dictionary = {}
        result_file_names = glob.glob("{}/*.rslts".format(results_dir_name))
        for results_path in result_file_names:

            results_file = open(results_path, "r")
            # Skip first line because it contains a header
            next(results_file)
            plot_points = []
            for line in results_file:
                xy = line.split("\t")
                x = float(xy[0])
                y = float(xy[1])
                plot_points.append((x, y))
            basename = os.path.basename(results_path)
            label, file_extension = os.path.splitext(basename)
            plot_points_dictionary[label] = plot_points
            results_file.close()
        return plot_points_dictionary

    d = collect_results(results_dir_name)
    # Get a list representation of the dictionary sorted by the reliability
    # of the states corresponding to the dictionary keys
    sorted_data = sorted(d.items(),
        key=lambda (k, lst): lst[int(time/2)][1],
        reverse=True)

    # since the results in the dictionary d are floating-point numbers, and
    # floating-point numbers have inherent inaccuracy, some of the Y
    # coordinates that should be identical in d are actually slightly
    # different. To fix this, we round the results to the number of decimal
    # places given by decimal_places.
    sorted_data_reduced_accuracy = []
    for (k, lst) in sorted_data:
        lst_reduced_accuracy = [
            (x, numerical_approx(y, digits=decimal_places))
            for (x, y) in lst]
        sorted_data_reduced_accuracy.append((k, lst_reduced_accuracy))
    # sorted_data_reduced_accuracy may contain items (x1, y1) and (x2, y2)
    # such that y1 == y2. We create a new list labeled_data with such items
    # merged, i.e., the new list will contain ("x1, x2", y1) instead of
    # both (x1, y1) and (x2, y2).  The new list also contains all items
    # (x3, y3) of sorted_data_reduced_accuracy which have a unique y3.
    labeled_data = []
    for key, group in groupby(sorted_data_reduced_accuracy, lambda x: x[1]):
        new_label = [item[0] for item in group]
        new_label.sort(key=int)
        new_label = ", ".join(new_label)
        labeled_data.append((new_label, key))

    # create the plot
    loop_count = 0
    reliability_plot = line2d([(0, 0)], axes_labels=["$t$", "$R(t)$"])
    colors = rainbow(len(labeled_data))
    label_xpos = time * 0.2
    label_ypos = 1.1
    for label, plot_points in labeled_data:
        # plot the label
        reliability_plot += text(label, (label_xpos, label_ypos),
            horizontal_alignment="left")
        label_arrow_start = (label_xpos - time * 0.01, label_ypos)
        label_arrow_end = \
            plot_points[(int(len(plot_points)/10) + loop_count) %
                len(plot_points)]
        reliability_plot += line2d([label_arrow_start, label_arrow_end])
        label_xpos += time * 0.05
        label_ypos -= 0.05

        # plot the reliability data
        current_color = colors.pop()
        reliability_plot += line2d(plot_points,
            linestyle=":",
            color=current_color,
            marker=".",
            markerfacecolor=current_color)
        loop_count += 1
        if loop_count >= max_curves:
            break
    return reliability_plot

########################################
# unit tests
########################################


def generate_subgraphs(G):
    generated_subgraphs = []
    for V in powerset(G.vertices()):
        for E in powerset(G.edges()):
            accept = True
            for e in E:
                if e[0] not in V or e[1] not in V:
                    accept = False
            if accept:
                H = Graph()
                H.add_vertices(V)
                H.add_edges(E)
                generated_subgraphs.append(H)
    return generated_subgraphs



class AbstractTestGenerateCtmc(object):
    """
    Test that the function generate_ctmc() generates a correct transition rate
    matrix using a given graph and a given set of failure rates for each
    element of the graph as input.
    """
    # G is the graph to use for testing. Subclasses must override this.
    G = None
    # dictionary indicating the failure rate of each edge of the graph
    edge_failure_rates = {}
    # dictionary indicating the failure rate of each vertex of the graph
    vertex_failure_rates = {}
    # The transition rate matrix that the function generate_ctmc()
    # should return if it receives as input the graph G and the failure
    # rate dictionary failure_rates.
    expected_trans_matrix = None

    def setUp(self):
        self.failure_rates = self.edge_failure_rates.copy()
        self.failure_rates.update(self.vertex_failure_rates)
        self.dg, self.trans_matrix = generate_ctmc(self.G, self.failure_rates)

        self.G_subgraphs = generate_subgraphs(self.G)
        self.G_subgraph_count = len(self.G_subgraphs)
        # list of graphs that make up the vertices of the digraph
        self.dg_graphs = [self.dg.get_vertex(v) for v in self.dg.vertices()]

    def test_digraph_order_matches_G_subgraph_count(self):
        """
        Test that the number of vertices of the digraph matches the number of
        subgraphs of G.
        """
        self.assertEqual(self.dg.order(), self.G_subgraph_count)

    def test_digraph_has_subgraphs_of_G_as_vertices(self):
        """
        Test that each of the subgraphs of G is a vertex of the digraph.
        """
        for g in self.G_subgraphs:
            self.assertIn(g, self.dg_graphs)

    def test_digraph_first_vertex_is_G(self):
        """
        Test that the first vertex that has been added to the digraph dg is the
        graph G.
        """
        first_vertex = self.dg.vertices()[0]
        self.assertEqual(self.G, self.dg.get_vertex(first_vertex))

    def test_digraph_first_vertex_has_in_degree_0(self):
        """
        Test that the vertex corresponding to the graph G does not have any
        incoming edges.
        """
        first_vertex = self.dg.vertices()[0]
        self.assertEqual(self.dg.in_degree(first_vertex), 0)

    def test_digraph_has_empty_graph_as_vertex(self):
        """
        Test that the empty graph is a vertex of the digraph.
        """
        empty_graph = Graph()
        self.assertIn(empty_graph, self.dg_graphs)

    def test_digraph_has_G_as_vertex(self):
        """
        Test that the graph G is a vertex of the digraph.
        """
        self.assertIn(self.G, self.dg_graphs)

    def test_each_vertex_has_correct_out_degree(self):
        """
        Test that each vertex of dg has the correct number of outgoing edges,
        which is the number of edges of the graph corresponding to the vertex
        plus the number of vertices of that graph.  This is so because outgoing
        edges represent either the failure of an edge or the failure of a
        vertex.
        """
        for v in self.dg.vertices():
            graph_of_v = self.dg.get_vertex(v)
            self.assertEqual(self.dg.out_degree(v),
                graph_of_v.order() + graph_of_v.size())

    def test_in_neighbors_have_1_more_edge_or_vertex(self):
        """
        Test that all the in-neighbors of a given vertex v are graphs with one
        additional vertex or one additional edge.
        """
        for v in self.dg.vertices():
            graph_of_v = self.dg.get_vertex(v)
            for n in self.dg.neighbors_in(v):
                graph_of_n = self.dg.get_vertex(n)
                self.assertTrue(
                    (graph_of_n.order() == graph_of_v.order() + 1) or
                    (graph_of_n.size() == graph_of_v.size() + 1)
                    )

    def test_out_neighbors_have_1_less_edge_or_vertex(self):
        """
        Test that all the out-neighbors of a given vertex v are graphs with one
        less vertex or one less edge.
        """
        for v in self.dg.vertices():
            graph_of_v = self.dg.get_vertex(v)
            for n in self.dg.neighbors_out(v):
                graph_of_n = self.dg.get_vertex(n)
                self.assertTrue(
                    (graph_of_n.order() == graph_of_v.order() - 1) or
                    (graph_of_n.size() == graph_of_v.size() - 1)
                    )

    def test_empty_graph_vertex_has_out_degree_0(self):
        """
        Test that the vertex corresponding to the empty graph does not have any
        outgoing edges.
        """
        # find the vertex corresponding to the empty graph
        v = self.dg.vertices()[0]
        while self.dg.get_vertex(v) != Graph():
            v += 1

        self.assertEqual(self.dg.out_degree(v), 0)

    def test_is_square_matrix(self):
        """
        Test that the transition rate matrix is a square matrix.
        """
        self.assertTrue(self.trans_matrix.is_square())

    def test_rows_sum_zero(self):
        """
        Test that the elements of each row sum zero.
        """
        for row in self.trans_matrix.rows():
            self.assertEqual(sum(row), 0)

    def test_all_failure_rates_are_elements_of_matrix(self):
        """
        Test that each of the failure rates is a coefficient of the matrix
        returned by generate_ctmc().
        """
        matrix_elements = self.trans_matrix.list()
        for failure_rate in self.failure_rates.values():
            self.assertIn(failure_rate, matrix_elements)

    def test_trans_rate_matrix_is_correct(self):
        """
        Test that the function generate_ctmc() returns the expected matrix by
        checking whether the returned and the expected matrix are equal.
        """
        self.assertEqual(self.trans_matrix, self.expected_trans_matrix)

    def test_trans_rate_matrix_coefficients_are_correct(self):
        """
        Test that the function generate_ctmc() returns the expected matrix by
        comparing each coefficient of the returned and the expected matrix.
        """
        for row, column in zip(range(self.trans_matrix.nrows()),
            range(self.expected_trans_matrix.ncols())):
            self.assertEqual(self.trans_matrix[row, column],
                self.expected_trans_matrix[row, column],
                    msg="row {} differ:\n{}\nvs\n{}\n".format(row,
                    self.trans_matrix[row],
                    self.expected_trans_matrix[row]))





class TestGenerateCtmcWithEmptyGraph(
        AbstractTestGenerateCtmc,
        unittest.TestCase):
    """
    Test the function generate_ctmc() using the empty graph as input.
    """
    G = Graph()
    edge_failure_rates = {}
    vertex_failure_rates = {}
    expected_trans_matrix = Matrix(QQ, [0])

class TestGenerateCtmcWithTrivialGraph(
        AbstractTestGenerateCtmc,
        unittest.TestCase):
    """
    Test the function generate_ctmc() using the trivial graph as input.
    """
    def setUp(self):
        v0 = 'v0'
        self.G = Graph({v0: []})
        efr = self.edge_failure_rates = {}
        vfr = self.vertex_failure_rates = {
            v0: 5/100
        }
        self.expected_trans_matrix = Matrix(QQ, [
            [-vfr[v0], vfr[v0]],
            [0, 0]
        ])

        super(TestGenerateCtmcWithTrivialGraph, self).setUp()

class TestGenerateCtmcWithP2(AbstractTestGenerateCtmc,
        unittest.TestCase):
    """
    Test the function generate_ctmc() using the path graph P2 as input.
    """
    def setUp(self):
        v0 = 'v0'
        v1 = 'v1'
        self.G = Graph({v0: [v1]})
        efr = self.edge_failure_rates = {
            (v0, v1): 1/100
        }
        vfr = self.vertex_failure_rates = {
            v0: 2/100,
            v1: 3/100
        }
        # number of states of the CTMC corresponding to the P2 graph
        # (this equals the number of subgraphs of P2)
        ctmc_state_count = 5
        # initialize the matrix as a null matrix of size 5x5
        M = Matrix(QQ, ctmc_state_count)

        # set up row 0 of the transition rate matrix
        M[0, 1] = efr[(v0, v1)]
        M[0, 2] = vfr[v0]
        M[0, 4] = vfr[v1]
        row_sum = sum(M.row(0))
        M[0, 0] = -row_sum

        # set up row 1 of the transition rate matrix
        M[1, 2] = vfr[v0]
        M[1, 4] = vfr[v1]
        row_sum = sum(M.row(1))
        M[1, 1] = -row_sum

        # set up row 2 of the transition rate matrix
        M[2, 3] = vfr[v1]
        row_sum = sum(M.row(2))
        M[2, 2] = -row_sum

        # set up row 3 of the transition rate matrix
        pass

        # set up row 4 of the transition rate matrix
        M[4, 3] = vfr[v0]
        row_sum = sum(M.row(4))
        M[4, 4] = -row_sum

        self.expected_trans_matrix = M

        super(TestGenerateCtmcWithP2, self).setUp()


class TestGenerateCtmcWithP3(AbstractTestGenerateCtmc,
        unittest.TestCase):
    """
    Test the function generate_ctmc() using the path graph P3 as input.
    """
    def setUp(self):
        v0 = 'v0'
        v1 = 'v1'
        v2 = 'v2'
        self.G = Graph({v0: [v1],  v1: [v2]})
        efr = self.edge_failure_rates = {
            (v0, v1): 1/100,
            (v1, v2): 2/100
        }
        vfr = self.vertex_failure_rates = {
            v0: 3/100,
            v1: 4/100,
            v2: 5/100
        }
        # number of states of the CTMC corresponding to the P3 graph
        # (this equals the number of subgraphs of P3)
        ctmc_state_count = 13
        # initialize the matrix as a null matrix of size 13x13
        M = Matrix(QQ, ctmc_state_count)

        # set up row 0 of the transition rate matrix
        M[0, 1] = efr[(v0, v1)]
        M[0, 7] = vfr[v1]
        M[0, 10] = vfr[v0]
        M[0, 11] = efr[(v1, v2)]
        M[0, 12] = vfr[v2]
        row_sum = sum(M.row(0))
        M[0, 0] = -row_sum

        # set up row 1 of the transition rate matrix
        M[1, 2] = efr[(v1, v2)]
        M[1, 7] = vfr[v1]
        M[1, 9] = vfr[v2]
        M[1, 10] = vfr[v0]
        row_sum = sum(M.row(1))
        M[1, 1] = -row_sum

        # set up row 2 of the transition rate matrix
        M[2, 3] = vfr[v0]
        M[2, 7] = vfr[v1]
        M[2, 9] = vfr[v2]
        row_sum = sum(M.row(2))
        M[2, 2] = -row_sum

        # set up row 3 of the transition rate matrix
        M[3, 4] = vfr[v1]
        M[3, 6] = vfr[v2]
        row_sum = sum(M.row(3))
        M[3, 3] = -row_sum

        # set up row 4 of the transition rate matrix
        M[4, 5] = vfr[v2]
        row_sum = sum(M.row(4))
        M[4, 4] = -row_sum

        # set up row 5 of the transition rate matrix
        pass

        # set up row 6 of the transition rate matrix
        M[6, 5] = vfr[v1]
        row_sum = sum(M.row(6))
        M[6, 6] = -row_sum

        # set up row 7 of the transition rate matrix
        M[7, 4] = vfr[v0]
        M[7, 8] = vfr[v2]
        row_sum = sum(M.row(7))
        M[7, 7] = -row_sum

        # set up row 8 of the transition rate matrix
        M[8, 5] = vfr[v0]
        row_sum = sum(M.row(8))
        M[8, 8] = -row_sum

        # set up row 9 of the transition rate matrix
        M[9, 6] = vfr[v0]
        M[9, 8] = vfr[v1]
        row_sum = sum(M.row(9))
        M[9, 9] = -row_sum

        # set up row 10 of the transition rate matrix
        M[10, 3] = efr[(v1, v2)]
        M[10, 4] = vfr[v1]
        M[10, 6] = vfr[v2]
        row_sum = sum(M.row(10))
        M[10, 10] = -row_sum

        # set up row 11 of the transition rate matrix
        M[11, 2] = efr[(v0, v1)]
        M[11, 3] = vfr[v0]
        M[11, 7] = vfr[v1]
        M[11, 12] = vfr[v2]
        row_sum = sum(M.row(11))
        M[11, 11] = -row_sum

        # set up row 12 of the transition rate matrix
        M[12, 6] = vfr[v0]
        M[12, 8] = vfr[v1]
        M[12, 9] = efr[(v0, v1)]
        row_sum = sum(M.row(12))
        M[12, 12] = -row_sum

        self.expected_trans_matrix = M

        super(TestGenerateCtmcWithP3, self).setUp()

class TestGenerateCtmcWithC3(AbstractTestGenerateCtmc,
        unittest.TestCase):
    """
    Test the function generate_ctmc() using the path graph C3 as input.
    """
    def setUp(self):
        v0 = 'v0'
        v1 = 'v1'
        v2 = 'v2'
        self.G = Graph({v0: [v1, v2],  v1: [v2]})
        efr = self.edge_failure_rates = {
            (v0, v1): 1/100,
            (v0, v2): 2/100,
            (v1, v2): 3/100
        }
        vfr = self.vertex_failure_rates = {
            v0: 4/100,
            v1: 5/100,
            v2: 6/100
        }
        # number of states of the CTMC corresponding to the P3 graph
        # (this equals the number of subgraphs of P3)
        ctmc_state_count = 18
        # initialize the matrix as a null matrix of size 18x18
        M = Matrix(QQ, ctmc_state_count)

        # set up row 0 of the transition rate matrix
        M[0, 1] = efr[(v0, v1)]
        M[0, 11] = vfr[v0]
        M[0, 13] = vfr[v1]
        M[0, 14] = efr[(v0, v2)]
        M[0, 16] = vfr[v2]
        M[0, 17] = efr[(v1, v2)]
        row_sum = sum(M.row(0))
        M[0, 0] = -row_sum

        # set up row 1 of the transition rate matrix
        M[1, 2] = efr[(v0, v2)]
        M[1, 10] = vfr[v2]
        M[1, 11] = vfr[v0]
        M[1, 12] = efr[(v1, v2)]
        M[1, 13] = vfr[v1]
        row_sum = sum(M.row(1))
        M[1, 1] = -row_sum

        # set up row 2 of the transition rate matrix
        M[2, 3] = efr[(v1, v2)]
        M[2, 8] = vfr[v1]
        M[2, 10] = vfr[v2]
        M[2, 11] = vfr[v0]
        row_sum = sum(M.row(2))
        M[2, 2] = -row_sum

        # set up row 3 of the transition rate matrix
        M[3, 4] = vfr[v0]
        M[3, 8] = vfr[v1]
        M[3, 10] = vfr[v2]
        row_sum = sum(M.row(3))
        M[3, 3] = -row_sum

        # set up row 4 of the transition rate matrix
        M[4, 5] = vfr[v1]
        M[4, 7] = vfr[v2]
        row_sum = sum(M.row(4))
        M[4, 4] = -row_sum

        # set up row 5 of the transition rate matrix
        M[5, 6] = vfr[v2]
        row_sum = sum(M.row(5))
        M[5, 5] = -row_sum

        # set up row 6 of the transition rate matrix
        pass

        # set up row 7 of the transition rate matrix
        M[7, 6] = vfr[v1]
        row_sum = sum(M.row(7))
        M[7, 7] = -row_sum

        # set up row 8 of the transition rate matrix
        M[8, 5] = vfr[v0]
        M[8, 9] = vfr[v2]
        row_sum = sum(M.row(8))
        M[8, 8] = -row_sum

        # set up row 9 of the transition rate matrix
        M[9, 6] = vfr[v0]
        row_sum = sum(M.row(9))
        M[9, 9] = -row_sum

        # set up row 10 of the transition rate matrix
        M[10, 7] = vfr[v0]
        M[10, 9] = vfr[v1]
        row_sum = sum(M.row(10))
        M[10, 10] = -row_sum

        # set up row 11 of the transition rate matrix
        M[11, 4] = efr[(v1, v2)]
        M[11, 5] = vfr[v1]
        M[11, 7] = vfr[v2]
        row_sum = sum(M.row(11))
        M[11, 11] = -row_sum

        # set up row 12 of the transition rate matrix
        M[12, 3] = efr[(v0, v2)]
        M[12, 4] = vfr[v0]
        M[12, 10] = vfr[v2]
        M[12, 13] = vfr[v1]
        row_sum = sum(M.row(12))
        M[12, 12] = -row_sum

        # set up row 13 of the transition rate matrix
        M[13, 5] = vfr[v0]
        M[13, 8] = efr[(v0, v2)]
        M[13, 9] = vfr[v2]
        row_sum = sum(M.row(13))
        M[13, 13] = -row_sum

        # set up row 14 of the transition rate matrix
        M[14, 2] = efr[(v0, v1)]
        M[14, 8] = vfr[v1]
        M[14, 11] = vfr[v0]
        M[14, 15] = efr[(v1, v2)]
        M[14, 16] = vfr[v2]
        row_sum = sum(M.row(14))
        M[14, 14] = -row_sum

        # set up row 15 of the transition rate matrix
        M[15, 3] = efr[(v0, v1)]
        M[15, 4] = vfr[v0]
        M[15, 8] = vfr[v1]
        M[15, 16] = vfr[v2]
        row_sum = sum(M.row(15))
        M[15, 15] = -row_sum

        # set up row 16 of the transition rate matrix
        M[16, 7] = vfr[v0]
        M[16, 9] = vfr[v1]
        M[16, 10] = efr[(v0, v1)]
        row_sum = sum(M.row(16))
        M[16, 16] = -row_sum

        # set up row 17 of the transition rate matrix
        M[17, 4] = vfr[v0]
        M[17, 12] = efr[(v0, v1)]
        M[17, 13] = vfr[v1]
        M[17, 15] = efr[(v0, v2)]
        M[17, 16] = vfr[v2]
        row_sum = sum(M.row(17))
        M[17, 17] = -row_sum

        self.expected_trans_matrix = M

        super(TestGenerateCtmcWithC3, self).setUp()


if __name__ == '__main__':
    unittest.main()
