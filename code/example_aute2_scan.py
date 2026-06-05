"""
Example: Planning an IXS scan for AuTe2
"""

from sixcircle_minimal import init_sixcircle
from aute2_structure import AuTe2

print("=" * 60)
print("AuTe2 IXS Measurement Planning")
print("=" * 60)

# Create AuTe2 crystal
crystal = AuTe2()
crystal.info()

# Initialize for BL43
print("\n" + "=" * 60)
print("Initialize for BL43LXU (Si 11,11,11)")
print("=" * 60)
sc = init_sixcircle(21.747)
sc.setup_crystal(crystal.a, crystal.b, crystal.c,
                 crystal.alpha, crystal.beta, crystal.gamma,
                 crystal.name)

# Calculate key reflections
print("\n" + "=" * 60)
print("Key Reflections")
print("=" * 60)

reflections = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
]

results = []
for h, k, l in reflections:
    result = sc.ca(h, k, l)
    if result:
        results.append(result)

print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Sample: {crystal.name}")
print(f"Energy: {sc.energy_kev:.3f} keV")
print(f"Accessible reflections: {len(results)}")
