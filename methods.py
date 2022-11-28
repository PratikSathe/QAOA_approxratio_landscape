
from qiskit import QuantumCircuit, Aer, execute

reverseStep = -1 # To take into account the little-endian format used in qiskit

def get_metric_for_angles(nodes, edges, betas, gammas, num_shots, metric_type = 'approx_ratio'):
    cuts, counts = get_cut_distribution(nodes, edges, betas, gammas, num_shots)
    
    
    
    
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
        self.sizes = [eval_cut(self.nodes, self.edges, solution[::-1]) for solution in cuts] # Reverse each cut passed ot eval_cut, since qiskit uses the little-endian format
        
    
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
        




