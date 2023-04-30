import msprime, tskit, pyslim

# Based on tutorial at: https://tskit.dev/pyslim/docs/latest/tutorial.html#adding-neutral-mutations-to-a-slim-simulation

# 1. Recapitation:
print("\n Loading 10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg.trees with tskit.")
orig_ts = tskit.load("10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg.trees")
print(f" The original tree sequence now has {orig_ts.num_mutations} mutations,\n"
      f"  and mean pairwise nucleotide diversity is {orig_ts.diversity():0.3e}.")
      
print("\n Recapitating 10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg.trees with neutral mutations at 7*10^-8 mutation rate:\n")
rts = pyslim.recapitate(orig_ts,
            recombination_rate=7e-8,
            ancestral_Ne=10000, random_seed=1)
print(" DONE: 10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg.trees recapitated.\n")

# Adding neutral mutations:
print("\n Adding neutral mutations to 10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg.trees at a 1.4*10^-8 mutation rate:\n")

next_id = pyslim.next_slim_mutation_id(rts)
ts = msprime.sim_mutations(
           rts,
           rate=1.4e-8,
           model=msprime.SLiMMutationModel(type=0, next_id=next_id),
           keep=True,
)

print(f" DONE: The tree sequence now has {ts.num_mutations} mutations,\n"
      f"  and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")

ts.dump("10Ki_50Ms_5o6m8_4r8_simp100_dfeG_Ne__60Kg_recap_7m8_overlaid_1o4m8.trees")
print(" Saved trees file 10Ki_50Ms_5o6m8_4r8_simp100_dfeG_Ne__60Kg_recap_7m8_overlaid_1o4m8.trees")
#_____________________________________
# ts = tskit.load("10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg.trees")

# # Checks for coalescence
# print('Test coalescence before simplification:')
# for t in ts.trees():
    # assert t.num_roots == 1, ("not coalesced! on segment {} to {}".
        # format(t.interval[0], t.interval[1]))

# ts = ts.simplify()
# print('Test coalescence after simplification:')
# for t in ts.trees():
    # assert t.num_roots == 1, ("not coalesced! on segment {} to {}".
        # format(t.interval[0], t.interval[1]))
        
# mutated = msprime.sim_mutations(ts, rate=1.4e-8, random_seed=1, keep=True)
# mutated.dump("10Ki_50Ms_56m8_4r8_simp100_dfeG_Ne__60Kg_overlaid_7m8.trees")