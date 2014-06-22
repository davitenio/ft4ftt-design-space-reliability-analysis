import prismmodel as pm

def make_filter_function(slaves):
    def exists_connected_component_spanning_slaves(graph):
        connected_components = graph.connected_components()
        for cc in connected_components:
            if set(slaves).issubset(cc):
                return True
        return False
    return exists_connected_component_spanning_slaves



# Define the graph corresponding to the replicated switch
slaves = ["$s_1$", "$s_2$"]
switches = ["$w_1$", "$w_2$"]
replicated_switch = Graph()
replicated_switch.add_edges([(s, w) for s in slaves for w in switches])
replicated_switch.add_edges([(wa, wb) for wa in switches for wb in switches])

edge_failure_rates = {edge[:2]: 1/100 for edge in replicated_switch.edges()}
vertex_failure_rates = {vertex: 1/100 for vertex in replicated_switch.vertices()}

failure_rates = dict(edge_failure_rates.items() + vertex_failure_rates.items())

p = replicated_switch.plot()
p = p.save('design_space_graph.png')
digraph, transition_rate_matrix = pm.generate_ctmc(replicated_switch,
    failure_rates)
p = digraph.plot(layout='circular')
p = p.save('ctmc_digraph.png')
pm.prism_model(transition_rate_matrix, 'design_space.sm')
non_faulty_states = pm.filter_graphs(digraph,
    make_filter_function(slaves))
pm.prism_properties(non_faulty_states, 'design_space.pctl')

results_plot = pm.run_experiments('design_space.sm', 'design_space.pctl',
    non_faulty_states, 100)
results_plot.save("res.png")
