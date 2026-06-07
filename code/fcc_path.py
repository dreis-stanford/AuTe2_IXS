"""
High-symmetry path for FCC Brillouin zone
Define in conventional cubic reciprocal coordinates, then convert to primitive
"""

import numpy as np
from .fcc_structure import fcc_conv2prim_k


def get_fcc_high_symmetry_points_conventional():
    """
    Return high-symmetry points for FCC in CONVENTIONAL cubic reciprocal coordinates
    These are the standard textbook definitions
    
    Returns:
    --------
    points : dict
        Dictionary of high-symmetry points in conventional coordinates
    """
    points = {
        'G':     np.array([0.0, 0.0, 0.0]),      # Gamma - zone center
        'X':     np.array([0.0, 0.5, 0.5]),      # Face center
        'L':     np.array([0.5, 0.5, 0.5]),      # Corner
        'W':     np.array([0.25, 0.75, 0.5]),    # Edge
        'K':     np.array([0.375, 0.75, 0.375]), # Between W and X
        'U':     np.array([0.25, 0.625, 0.625]), # Edge (less common)
    }
    return points


def generate_path(points_list, n_points_per_segment=100):
    """
    Generate q-point path through high-symmetry points
    
    Parameters:
    -----------
    points_list : list of (name, coordinates) tuples
        Coordinates should be in CONVENTIONAL reciprocal coordinates
    n_points_per_segment : int or list of int
        Number of points between consecutive high-symmetry points
    
    Returns:
    --------
    q_path : ndarray (N, 3)
        Q-points along path in PRIMITIVE reciprocal coordinates
    q_distances : ndarray (N,)
        Distance along path for plotting x-axis
    labels : list of (position, name)
        High-symmetry point labels and positions for plotting
    """
    
    if isinstance(n_points_per_segment, int):
        n_points_per_segment = [n_points_per_segment] * (len(points_list) - 1)
    
    q_path_conv = []
    
    # Generate path in conventional coordinates
    for i in range(len(points_list) - 1):
        start_name, start_pt = points_list[i]
        end_name, end_pt = points_list[i + 1]
        n_pts = n_points_per_segment[i]
        
        # Generate points along segment (excluding end to avoid duplicates)
        for j in range(n_pts):
            t = j / n_pts
            q_pt = (1 - t) * start_pt + t * end_pt
            q_path_conv.append(q_pt)
    
    # Add final point
    q_path_conv.append(points_list[-1][1])
    
    # Convert entire path to primitive coordinates using existing function
    q_path_conv = np.array(q_path_conv)
    q_path_prim = fcc_conv2prim_k(q_path_conv)
    
    # Calculate distances in primitive reciprocal space
    q_distances = [0.0]
    current_distance = 0.0
    
    for i in range(1, len(q_path_prim)):
        delta = np.linalg.norm(q_path_prim[i] - q_path_prim[i-1])
        current_distance += delta
        q_distances.append(current_distance)
    
    q_distances = np.array(q_distances)
    
    # Create labels at high-symmetry points
    labels = [(0.0, points_list[0][0])]
    
    cumulative_pts = 0
    for i in range(len(points_list) - 1):
        cumulative_pts += n_points_per_segment[i]
        labels.append((q_distances[cumulative_pts], points_list[i + 1][0]))
    
    return q_path_prim, q_distances, labels


def get_qe_fcc_path(n_points=100):
    """
    Generate QE/MatdynX standard FCC path: G - K - X - G - L - X - W - L
    This matches the path shown in your reference plot
    
    Parameters:
    -----------
    n_points : int
        Number of points per segment
    
    Returns:
    --------
    q_path : ndarray
        Q-points in PRIMITIVE reciprocal coordinates (for dynamical matrix calculation)
    q_distances : ndarray
        Distances for plotting
    labels : list
        High-symmetry point labels
    """
    
    hs_points = get_fcc_high_symmetry_points_conventional()
    
    # QE standard path: G - K - X - G - L - X - W - L
    path = [
        ('G', hs_points['G']),
        ('K', hs_points['K']),
        ('X', hs_points['X']),
        ('G', hs_points['G']),
        ('L', hs_points['L']),
        ('X', hs_points['X']),
        ('W', hs_points['W']),
        ('L', hs_points['L']),
    ]
    
    return generate_path(path, n_points)


if __name__ == "__main__":
    # Test the conversion
    print("FCC Reciprocal Lattice Coordinate Conversion")
    print("=" * 60)
    
    hs_conv = get_fcc_high_symmetry_points_conventional()
    
    print("\nHigh-symmetry points:")
    print(f"{'Point':<8} {'Conventional':<25} {'Primitive':<25}")
    print("-" * 60)
    
    for name, pt_conv in hs_conv.items():
        pt_prim = fcc_conv2prim_k(pt_conv)
        conv_str = f"[{pt_conv[0]:.3f}, {pt_conv[1]:.3f}, {pt_conv[2]:.3f}]"
        prim_str = f"[{pt_prim[0]:.3f}, {pt_prim[1]:.3f}, {pt_prim[2]:.3f}]"
        print(f"{name:<8} {conv_str:<25} {prim_str:<25}")
    
    # Test path generation
    print("\n" + "=" * 60)
    print("QE Standard FCC Path: G - K - X - G - L - X - W - L")
    print("=" * 60)
    
    q_path, q_dist, labels = get_qe_fcc_path(50)
    
    print(f"\nTotal q-points: {len(q_path)}")
    print(f"Path length: {q_dist[-1]:.3f}")
    print("\nHigh-symmetry points:")
    for dist, name in labels:
        print(f"  {name:<8} at distance {dist:.3f}")
