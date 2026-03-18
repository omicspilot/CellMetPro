"""Development script to explore COBRA model structure for pathway analysis.

Run this script to see what the input data looks like before implementing
get_pathway_mapping().
"""

import cobra

# Load the E. coli core model (built into COBRApy)
# This is a small model with ~95 reactions and well-defined subsystems
print("Loading E. coli core model...")
model = cobra.io.load_model("RECON1")

print(f"\nModel: {model.id}")
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes: {len(model.genes)}")


# Manual counting to show what you need to implement
subsystem_counts = {}
for rxn in model.reactions:
    subsystem = rxn.subsystem
    if subsystem:  # Skip empty/None
        if subsystem not in subsystem_counts:
            subsystem_counts[subsystem] = 0
        subsystem_counts[subsystem] += 1

# ============================================================================
# Explore the subsystem attribute
# ============================================================================

print("\n" + "=" * 60)
print("EXPLORING REACTION SUBSYSTEMS")
print("=" * 60)

# ============================================================================
# Count reactions per subsystem
# ============================================================================

print("\n" + "=" * 60)
print("SUBSYSTEM SUMMARY")
print("=" * 60)

# Manual counting to show what you need to implement
subsystem_counts = {}
for rxn in model.reactions:
    subsystem = rxn.subsystem
    if subsystem:  # Skip empty/None
        if subsystem not in subsystem_counts:
            subsystem_counts[subsystem] = 0
        subsystem_counts[subsystem] += 1

print(f"\nFound {len(subsystem_counts)} unique subsystems:\n")
for subsystem, count in sorted(subsystem_counts.items(), key=lambda x: -x[1]):
  if rxn.subsystem is not None or rxn.subsystem != "":
    print(f"  {subsystem:40s} : {count} reactions")

# # ============================================================================
# # Show expected output format
# # ============================================================================

# print("\n" + "=" * 60)
# print("EXPECTED OUTPUT FORMAT FOR get_pathway_mapping()")
# print("=" * 60)

# print("""
# Your function should return a dict like this:

# {
#     "Glycolysis/Gluconeogenesis": ["PFK", "FBA", "TPI", "GAPD", ...],
#     "Citric Acid Cycle": ["CS", "ACONT", "ICDHyr", "AKGDH", ...],
#     "Oxidative Phosphorylation": ["NADH16", "CYTBD", "ATPS4r", ...],
#     ...
# }
# """)

# ============================================================================
# Show edge cases
# ============================================================================

print("=" * 60)
print("EDGE CASES TO HANDLE")
print("=" * 60)

# Check for reactions without subsystems
no_subsystem = [rxn.id for rxn in model.reactions if not rxn.subsystem]
print(f"\nReactions without subsystem: {len(no_subsystem)}")
if no_subsystem:
    print(f"  Examples: {no_subsystem[:5]}")

# Check for None vs empty string
none_subsystem = [rxn.id for rxn in model.reactions if rxn.subsystem is None]
empty_subsystem = [rxn.id for rxn in model.reactions if rxn.subsystem == ""]
print(f"\nReactions with subsystem=None: {len(none_subsystem)}")
print(f"Reactions with subsystem='': {len(empty_subsystem)}")


for rxn in model.reactions:
    print(f"flux_expression: {rxn.flux_expression}")
    print(f"forward_variable: {rxn.forward_variable}")
    print(f"reverse_variable: {rxn.reverse_variable}")
    print(f"objective_coefficient: {rxn.objective_coefficient}")
    print(f"bounds: {rxn.bounds}")
    print(f"upper_bound: {rxn.upper_bound}")
    # print(f"flux: {rxn.flux}")
    # print(f"reduced_cost: {rxn.reduced_cost}")
    print(f"metabolites: {rxn.metabolites}")
    print(f"genes: {rxn.genes}")
    print(f"gene_reaction_rule: {rxn.gene_reaction_rule}")
    print(f"gene_name_reaction_rule: {rxn.gene_name_reaction_rule}")
    print(f"gpr: {rxn.gpr}")
    print(f"functional: {rxn.functional}")
    print(f"reversibility: {rxn.reversibility}")
    print(f"reactants: {rxn.reactants}")
    print(f"products: {rxn.products}")
    print(f"reaction: {rxn.reaction}")
    print(f"compartments: {rxn.compartments}")
    print("=" * 60)
# # ============================================================================
# # Your task
# # ============================================================================

# print("\n" + "=" * 60)
# print("YOUR TASK")
# print("=" * 60)
# print("""
# Implement get_pathway_mapping(model) in:
#   cellmetpro/analysis/pathway.py

# It should:
# 1. Take a cobra.Model as input
# 2. Return a dict[str, list[str]] mapping subsystem -> reaction IDs
# 3. Skip reactions with empty/None subsystems

# Test it with:
#   from cellmetpro.analysis.pathway import get_pathway_mapping
#   mapping = get_pathway_mapping(model)
#   print(mapping)
# """)
