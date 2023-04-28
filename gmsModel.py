from cobra import Model, Reaction, Metabolite, Gene
from cobra.io import load_model

model = Model('gms_model')

reaction = Reaction('R_3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

#load the reactions from a list of reactions
with open("reactions_list", "r") as f:
    line = f.read().splitlines()
    reactions_names = line[0].split(",")

reactions_dict = {}
# Loop through the list of reaction names and perform desired operations on each reaction
for name in reactions_names:
    reaction = Reaction(name)
    #reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
    #reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = -1000.
    reaction.upper_bound = 1000.
    reactions_dict[name] = reaction

# read parameter names from file
with open("metabolites_list.txt", "r") as f:
    metabolites_names = f.read().splitlines()


# Create an empty dictionary to store metabolites and coefficients
metabolites_dict = {}

# Read the file and parse each line
with open('s_matrix.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        parts = line.strip().split()
        reaction_name = parts[1].lstrip(".")
        reaction = reactions_dict[reaction_name]

        metabolite_coeff = float(parts[2])
        metabolite_obj = Metabolite(parts[0])
        
        reaction.add_metabolites({metabolite_obj: metabolite_coeff})

print(reactions_dict['ZNabcpp'])
print(reactions_dict['TPI'])

# Create an empty dictionary to store genes and reactions
genes_dict = {}
with open('Rules_Reaction_matrix.txt', 'r') as f:
    for line in f:
        if line.startswith('/'):
            # End of the reaction block
            continue
        reaction_name, stoichiometry, gene_id = line.strip().split()
        stoichiometry = float(stoichiometry)

        gene_id = gene_id.lstrip(".")
        gene_obj = Gene(gene_id)
        if reaction_name in genes_dict:
            genes_dict[reaction_name].append(gene_id)
        else:
            genes_dict[reaction_name] = [gene_id]

#for reaction_name in genes_dict.keys():
#    print(genes_dict[reaction_name])
#    reaction = reactions_dict[reaction_name]
#    reaction.genes = genes_dict[reaction_name]

print("Example model:")


model_example = load_model("textbook")
for kei in model_example.reactions._dict.keys():
    reaction = model_example.reactions.get_by_id(kei)
    print(reaction.build_reaction_string())
    kei_pp = kei+"pp"
    if kei in reactions_dict:
        print(reactions_dict[kei])
        print(model_example.reactions._dict[kei])
    elif kei_pp in reactions_dict:
        print("pp: ") 
        print(reactions_dict[kei_pp])
        print(model_example.reactions._dict[kei])
    else:
        print(kei + " missing")
#reaction_sample = model_example.reactions['ZNabcpp']
