import glob
import re
import numpy as np

def get_delta(name):
    m = re.search(r'([+-]\d+\.\d+)', name)
    return float(m.group(1))

value_files = sorted(glob.glob("eigenvalues*.dat"), key=get_delta)
vector_files = sorted(glob.glob("eigenvector_*.dat"), key=get_delta)

deltas = [get_delta(f) for f in value_files]

energies = {}
vectors = {}

for f in value_files:
    d = get_delta(f)
    data = np.loadtxt(f, skiprows=2)
    energies[d] = data[:, 1]

for f in vector_files:
    d = get_delta(f)
    data = np.loadtxt(f, comments="#")
    nstate = int(data[:,1].max())
    nbasis = int(data[:,2].max())
    C = np.zeros((nbasis, nstate))
    for row in data:
        state = int(row[1]) - 1
        basis = int(row[2]) - 1
        coeff = row[7]
        C[basis, state] = coeff
    vectors[d] = C

deltas = sorted(set(energies.keys()) & set(vectors.keys()))

nstate = len(energies[deltas[0]])
track_state = list(range(nstate))

tracks = {t: [] for t in range(nstate)}

for t, s in enumerate(track_state):
    tracks[t].append((deltas[0], energies[deltas[0]][s]))

for a, b in zip(deltas[:-1], deltas[1:]):
    C1 = vectors[a]
    C2 = vectors[b]

    overlap = np.abs(C1.T @ C2)

    new_track_state = [-1] * nstate
    used = set()

    for t, old_state in enumerate(track_state):
        candidates = np.argsort(overlap[old_state, :])[::-1]
        for s2 in candidates:
            if s2 not in used:
                new_track_state[t] = s2
                used.add(s2)
                break

    track_state = new_track_state

    for t, s in enumerate(track_state):
        tracks[t].append((b, energies[b][s]))

with open("tracks.dat", "w") as f:
    for t in range(nstate):
        for d, e in tracks[t]:
            f.write(f"{d: .5f} {e: .10f} {t+1}\n")
        f.write("\n")