import pandas as pd
import cobra
from cobra.io import read_sbml_model

# Load the SBML model
model = read_sbml_model("model.xml")

# Load the ABC values
abc_values = pd.read_csv("abc_values.csv")

# Set the ABC values in the model
for index, row in abc_values.iterrows():
    rxn_id = row["reaction_id"]
    slope = row["slope"]
    intercept1 = row["intercept1"]
    intercept2 = row["intercept2"]
    reaction = model.reactions.get_by_id(rxn_id)
    reaction.bounds = (intercept2, intercept1)

# Load the carbon source uptake fluxes
carbon_uptake = pd.read_csv("carbon_uptake.csv")

# Load the expression data
expression_data = pd.read_csv("geneexpression_rxn.csv")

# Set the expression data in the model
for index, row in expression_data.iterrows():
    gene_id = row["gene_id"]
    reaction_id = row["reaction_id"]
    expression_value = row["expression_value"]
    gene = model.genes.get_by_id(gene_id)
    reaction = model.reactions.get_by_id(reaction_id)
    gene.reactions.append(reaction)
    reaction.genes.append(gene)
    reaction.add_metabolites({gene: expression_value})

# Define the list of measured reactions
measured_rxns = [
    "EX_glc_e",
    "EX_etoh_e",
    "EX_ac_e",
    "EX_lac-D_e",
    "EX_succ_e",
    "EX_pyr_e",
    "EX_for_e",
    "Ec_biomass_iJO1366_core_53p95M",
]

# Define the list of knockouts
knockouts = [
    "b0756",
    "b2388",
    "b0688",
    "b4025",
    "b3916",
    "b1723",
    "b4232",
    "b2097",
    "pseudogene",
    "b0755",
    "b4395",
    "b1854",
    "b1676",
    "b1702",
    "b1852",
    "b0767",
    "b2029",
    "b3386",
    "b2914",
    "b4090",
    "b2935",
    "b2465",
    "b2464",
    "b0008",
]

# Define the list of conditions
conditions = ["c{}".format(i) for i in range(1, 30)]

# Define the test set
test_set = set(conditions) - set(["c{}".format(i) for i in range(1, 25)])

# Set up the LP problem
lbfba = cobra.Model("lbfba")

# Add the reactions from the original model to the LP problem
for reaction in model.reactions:
    lbfba.add_reaction(reaction.copy())

# Add the gene expression variables to the LP problem
for reaction in lbfba.reactions:
    for gene in reaction.genes:
        lbfba.add_variable(
            "{}_{}".format(gene.id, reaction.id),
            lb=0.0,
            ub=1.0,
            type="binary",
        )

# Add the mass balance constraints
for metabolite in lbfba.metabolites:
    for condition in conditions:
        lbfba.add_constraint(
            cobra.Model.Constraint(
                0,
                lb=0,
                ub=0,
                name="mass_balance_{}_{}".format(metabolite.id, condition),
            )
        )
        for reaction in metabolite.reactions:
            stoichiometry = reaction.get_coefficient(metabolite)
            if stoichiometry != 0:
                lbfba.constraints[
                    "mass_balance_{}_{}".format(metabolite.id, condition)
                ].set_linear_coefficients(
                    {
                        "{}_{}".format(condition, reaction.id): stoichiometry,
                    }
                )

# Add the P/O ratio constraints
for condition in conditions:
    lbfba.add_constraint(
        cobra.Model.Constraint(
            0,
            lb=0,
            ub=0,
            name="PO_ratio_1_{}".format(condition),
        )
    )
    lbfba.add_constraint(
        cobra.Model.Constraint(
            0,
            lb=0,
            ub=0,
            name="PO_ratio_2_{}".format(condition),
        )
    )
    lbfba.add_constraint(
        cobra.Model.Constraint(
            0,
            lb=0,
            ub=0,
            name="PO_ratio_3_{}".format(condition),
        )
    )
    lbfba.constraints["PO_ratio_1_{}".format(condition)].set_linear_coefficients(
        {
            "{}_{}".format(condition, "NADH16pp"): 1,
            "{}_{}".format(condition, "NADH5"): -1,
        }
    )
    lbfba.constraints["PO_ratio_2_{}".format(condition)].set_linear_coefficients(
        {
            "{}_{}".format(condition, "NADH17pp"): 1,
            "{}_{}".format(condition, "NADH10"): -1,
        }
    )
    lbfba.constraints["PO_ratio_3_{}".format(condition)].set_linear_coefficients(
        {
            "{}_{}".format(condition, "NADH18pp"): 1,
            "{}_{}".format(condition, "NADH9"): -1,
        }
    )

# Add the forward and backward flux variables
for condition in conditions:
    for reaction in lbfba.reactions:
        lbfba.add_variable(
            "vfor_{}_{}".format(condition, reaction.id),
            lb=0.0,
            ub=1000.0,
            type="continuous",
        )
        lbfba.add_variable(
            "vback_{}_{}".format(condition, reaction.id),
            lb=0.0,
            ub=1000.0,
            type="continuous",
        )

# Add the objective function and the lower and upper bound constraints
minflux = cobra.Model.Variable("minflux", lb=-1000.0, ub=1000.0, type="continuous")
lbfba.add_variable(minflux)

# Set lower and upper bound factors
lower_bound_factor = 0.01
upper_bound_factor = 100.0

# Load gene expression data
gene_expression_data = pd.read_csv("geneexpression_rxn.csv", index_col=0)

# Add the intracellular flux variables to the LP problem
for condition in conditions:
    for reaction in lbfba.reactions:
        if reaction.id in intracellular_fluxes:
            variable_name = "{}_{}".format(condition, reaction.id)
            lbfba.add_variable(
                cobra.Model.Variable(
                    variable_name, lb=-1000.0, ub=1000.0, type="continuous"
                )
            )
            if variable_name in gene_expression_data.index:
                gene_expression_value = gene_expression_data.loc[
                    variable_name, condition
                ]
                if not pd.isna(gene_expression_value):
                    gene_expression_value = float(gene_expression_value)
                    lbfba.constraints[variable_name].set_linear_coefficients(
                        {
                            variable_name: -1.0,
                            "vfor_{}_{}".format(condition, reaction.id): 1.0,
                            "vback_{}_{}".format(condition, reaction.id): -1.0,
                        }
                    )
                    lbfba.constraints[variable_name].lb = (
                        lower_bound_factor * gene_expression_value
                    )
                    lbfba.constraints[variable_name].ub = (
                        upper_bound_factor * gene_expression_value
                    )
                    # Add the mass balance constraints
for metabolite in lbfba.metabolites:
    for condition in conditions:
        lbfba.add_constraint(
            cobra.Model.Constraint(
                0,
                lb=0,
                ub=0,
                name="mass_balance_{}_{}".format(metabolite.id, condition),
            )
        )
        for reaction in metabolite.reactions:
            stoichiometry = reaction.get_coefficient(metabolite)
            if stoichiometry != 0:
                lbfba.constraints[
                    "mass_balance_{}_{}".format(metabolite.id, condition)
                ].set_linear_coefficients(
                    {
                        "{}_{}".format(condition, reaction.id): stoichiometry,
                    }
                )

# Add the P/O ratio constraints
for condition in conditions:
    lbfba.add_constraint(
        cobra.Model.Constraint(
            0,
            lb=0,
            ub=0,
            name="PO_ratio_1_{}".format(condition),
        )
    )
    lbfba.add_constraint(
        cobra.Model.Constraint(
            0,
            lb=0,
            ub=0,
            name="PO_ratio_2_{}".format(condition),
        )
    )
    lbfba.add_constraint(
        cobra.Model.Constraint(
            0,
            lb=0,
            ub=0,
            name="PO_ratio_3_{}".format(condition),
        )
    )
    lbfba.constraints["PO_ratio_1_{}".format(condition)].set_linear_coefficients(
        {
            "{}_{}".format(condition, "NADH16pp"): 1,
            "{}_{}".format(condition, "NADH5"): -1,
        }
    )
    lbfba.constraints["PO_ratio_2_{}".format(condition)].set_linear_coefficients(
        {
            "{}_{}".format(condition, "NADH17pp"): 1,
            "{}_{}".format(condition, "NADH10"): -1,
        }
    )
    lbfba.constraints["PO_ratio_3_{}".format(condition)].set_linear_coefficients(
        {
            "{}_{}".format(condition, "NADH18pp"): 1,
            "{}_{}".format(condition, "NADH9"): -1,
        }
    )

# Add the forward and backward flux variables
for condition in conditions:
    for reaction in lbfba.reactions:
        lbfba.add_variable(
            "vfor_{}_{}".format(condition, reaction.id),
            lb=0.0,
            ub=1000.0,
            type="continuous",
        )
        lbfba.add_variable(
            "vback_{}_{}".format(condition, reaction.id),
            lb=0.0,
            ub=1000.0,
            type="continuous",
        )

# Add the objective function and the lower and upper bound constraints
minflux = cobra.Model.Variable("minflux", lb=-1000.0, ub=1000.0, type="continuous")
lbfba.add_variable(minflux)

for condition in conditions:
    for reaction in lbfba.reactions:
        if reaction.id in knockouts:
            lbfba.constraints[
                "{}_{}".format(condition, reaction.id)
            ].set_linear_coefficients({"vfor_{}_{}".format(condition, reaction.id): 1})
        else:
            lbfba.constraints[
                "{}_{}".format(condition, reaction.id)
            ].set_linear_coefficients(
                {
                    "vfor_{}_{}".format(condition, reaction.id): 1,
                    "vback_{}_{}".format(condition, reaction.id): -1,
                }
            )