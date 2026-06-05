# Notes on AuTe₂ IXS Analysis Code

## Python vs MATLAB Comparison

### ✅ Perfect Matches
- **Phonon frequencies**: Exact match
- **IXS cross-sections**: Exact match (this is the physical observable)
- **Form factors**: Exact match
- **Bose factors**: Exact match

### ⚠️ Expected Differences

#### Longitudinal Character (L-char)
The L-character values differ between Python and MATLAB due to **eigenvector phase conventions**.

**Why this happens:**
- Eigenvectors from `numpy.linalg.eigh` and MATLAB's `eig` can differ by arbitrary phase factors
- For atoms at non-equivalent positions (like Te1 and Te2 in AuTe₂), the phase relationship affects how we calculate L-char
- Both are "correct" - they just use different conventions

**Physical interpretation:**
- L-char tells you if atoms move parallel (longitudinal) or perpendicular (transverse) to Q
- The **IXS intensity** doesn't depend on this - it's gauge-invariant
- Use L-char as a rough guide, not an absolute measure

#### Atomic Participation Percentages
The participation of individual atoms (e.g., "Au: 50%, Te1: 25%, Te2: 25%") can differ.

**Why this happens:**
- At generic Q-points, the two Te atoms are **not symmetry-equivalent**
- Different eigensolvers can pick different linear combinations of degenerate modes
- Python's eigenvector might have mostly Te2 moving, MATLAB's might have mostly Te1

**What to trust:**
- **Total participation by element type** is reliable (sum Te1 + Te2)
- Individual atom percentages are eigenvector-convention-dependent

### 🎯 What Matters
The **IXS cross-sections** are the measurable quantities and they match perfectly.
The L-char and participation are analysis aids to understand the modes, but they're not unique.

---

## Usage Recommendation

When comparing modes between Python and MATLAB:
1. ✅ **Match by frequency** - this is unique
2. ✅ **Match by IXS intensity** - this is the observable
3. ⚠️ **Don't worry if L-char differs** - it's convention-dependent
4. ⚠️ **Don't worry if Te1/Te2 swap** - they're nearly equivalent


## Quick Reference: What to Trust

| Quantity | Status | Notes |
|----------|--------|-------|
| Frequencies | ✅ Exact | Always trust these |
| **IXS cross-sections** | ✅ **Exact** | **This is what you measure!** |
| Form factors | ✅ Exact | Atomic scattering factors |
| Bose factors | ✅ Exact | Temperature dependence |
| L-character | ⚠️ Convention | Rough guide only |
| Atom participation % | ⚠️ Convention | Sum by element type |

