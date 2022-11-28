
from qiskit import QuantumCircuit, Aer, execute

def get_size_dist(counts, sizes):
    """ For given measurement outcomes, i.e. combinations of cuts, counts and sizes, return counts corresponding to each cut size.
    Args:
        counts (list): List of integers, denoting number of times each cut was measured 
        sizes (list): List of integers. Cut size corresponding to each measured cut
    Returns:
        unique_counts (list) : Number of times each size was obtained
        unique_sizes (list) : List of all sizes
    """
    unique_sizes = list(set(sizes))
    unique_counts = [0] * len(unique_sizes)
    
    for i, size in enumerate(unique_sizes):
        # corresp_counts = [counts[ind] for ind,s in enumerate(sizes) if s == size]
        corresp_counts = [c for (c,s) in zip(counts, sizes) if s == size]
        unique_counts[i] = sum(corresp_counts)
    
    # Make sure that the scores are in non-decreasing order
    s_and_c_list = [[a,b] for (a,b) in zip(unique_sizes, unique_counts)] #size and cut list
    s_and_c_list = sorted(s_and_c_list, key = lambda x : x[0]) # sort according to sizes
    unique_sizes = [x[0] for x in s_and_c_list]
    unique_counts = [x[1] for x in s_and_c_list]
    return unique_counts, unique_sizes
    
    
class qaoa_anstaz:
    def __init__(self, nodes, edges, betas, gammas, num_shots):
        """
        Initialize a class for a QAOA ansatz circuit

        Args:
            nodes (int): number of nodes in the graph
            edges (list of tuples): list of (undirected) edges of the graph. Each element is a two-tuple of node labels, which should be between 0 and (nodes-1) 
            betas (list of floats): List of β angles. Number of rounds is equal to the length of this list.
            gammas (list of floats): List of γ angles. Number of rounds is equal to the length of this list.
            num_shots (int): Number of times the qaoa ansatz circuit is to be measured
        """
        self.nodes = nodes
        self.edges = edges
        self.betas = betas
        self.gammas = gammas
        self.num_shots = num_shots
        
        assert len(self.betas) == len(self.gammas), "Betas and Gammas lists should be of the same length"
        self.rounds = len(self.betas) 
        
        self.create_ansatz_circuit()
        
    def eval_cut(self, solution):
        """Compute the cut size corresponding to a given solution

        Args:
            solution : List of 0s and 1s, denoting a cut

        Returns:
            int: cut size
        """
        cut_size = 0
        for i,j in self.edges:
            if solution[i] != solution[j]:
                cut_size += 1

        return cut_size
        
    def execute_circuit(self):
        """
        Implement a noiseless simulation of the ansatz circuit.
        Compute the measured cuts, along with the number of times they were measured
        Creates:
            cuts (list): List of bitstrings, denoting measured cuts
            counts (list): List of integers, denoting number of times each cut was measured 
            sizes (list): List of integers. Cut size corresponding to each measured cut
        """
        # Implement a noiseless simulation using the qasm_simulator 
        backend = Aer.get_backend('qasm_simulator')
        cuts_counts_dict = backend.run(self.circuit, nshots=self.num_shots).result().get_counts() # Dictionary of ('bitstring', number of times it was measured) pairs
        
        self.cuts = list(cuts_counts_dict.keys())
        self.counts = list(cuts_counts_dict.values())
        self.sizes = [self.eval_cut(solution[::-1]) for solution in self.cuts] # Reverse each cut passed ot eval_cut, since qiskit uses the little-endian format
        
    
    def create_ansatz_circuit(self):
        """
        Create the ansatz circuit (determined by the β and γ parameters)
        """
        self.circuit = QuantumCircuit(self.nodes)
        for i in range(self.nodes):
            # Apply the Hadamard gate on each qubit
            self.circuit.h(i)

        for roundIndex in range(self.rounds):
            # For each round, apply RZZ gates followed by RX gates
            beta = self.betas[roundIndex]
            gamma = self.gammas[roundIndex]

            for i,j in self.edges:
                # exp(i gamma/2 Z_i Z_j)
                self.circuit.rzz(- gamma, i, j)
                # circuit.rzz( -gamma, i, j)
            for i in range(self.nodes):
                # exp(-i beta X_i)
                self.circuit.rx(2 * beta, i)
                
        self.circuit.measure_all()
        




