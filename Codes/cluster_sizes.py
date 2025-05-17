import json

# Generate standard cluster sizes (N)
def generate_standard_N(limit=20):
    N_set = set()
    for i in range(limit):
        for k in range(limit):
            if i == 0 and k == 0:
                continue
            N = i**2 + i*k + k**2
            N_set.add(N)
    return sorted(N_set)

# Save to JSON file for future use
N_std = generate_standard_N()
with open("cluster_sizes.json", "w") as f:
    json.dump(N_std, f)