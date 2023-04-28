
import pandas as pd

# read reactions names from file
with open("reactions_list", "r") as f:
    reactions_names = f.read().splitlines()
# read parameter names from file
with open("metabolites_list.txt", "r") as f:
    metabolites_names = f.read().splitlines()

    
# read parameter names from file
with open("parameters_list.txt", "r") as f:
    param_names = f.read().splitlines()

# create params string
params = {
    "Vmax": 1000,
    "UpperLimits": {name: 1000 for name in param_names},
    "LowerLimits": {name: -1000 for name in param_names}
}

S = {}

with open('s_matrix.txt', 'r') as f:
    for line in f:
        data = line.strip().split()
        i = data[0]
        j = data[1]
        value = float(data[2])
        if i not in S:
            S[i] = {}
        S[i][j] = value

n = set(range(1, 9))
# read gene names from file
with open("gene_list.txt", "r") as f:
    gene_names = f.read().splitlines()

rules_reactions = {}
with open('Rules_Reaction_matrix.txt', 'r') as f:
    for line in f:
        data = line.strip().split()
        i = data[0]
        j = data[1]
        value = data[2]
        if i not in rules_reactions:
            rules_reactions[i] = {}
        rules_reactions[i][j] = value

# Load the data from a file into a DataFrame
df = pd.read_csv('Rules_Reaction_matrix.txt', sep=r'\s+', header=None, names=['j', 'n', 'value'])

# Display the DataFrame
print(df)