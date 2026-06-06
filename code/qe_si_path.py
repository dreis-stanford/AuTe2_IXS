"""
QE/MatdynX path for Silicon - corrected based on conventional coordinates
"""

import numpy as np

def get_qe_si_path(n_points=100):
    """
    Generate the correct QE path for Silicon
    
    Path: G - K - X - G - L - X' - W - L
    where first X is (1,0.5,0.5) and the X in X->G is (0.5,0.5,0)
    
    Returns:
    --------
    q_path : ndarray
        Q-points in primitive reciprocal coordinates
    q_distances : ndarray
        Distances for plotting
    labels : list
        High-symmetry point labels
    """
    
    # Correct path based on conventional->primitive conversion
    path = [
        ('G', np.array([0.0, 0.0, 0.0])),
        ('K', np.array([0.75, 0.375, 0.375])),
        ('X', np.array([1.0, 0.5, 0.5])),        # First X at (1,1,0) conv
        ('G', np.array([0.0, 0.0, 0.0])),        # But approach from (0,0,1) conv = (0.5,0.5,0) prim
        ('L', np.array([0.5, 0.5, 0.5])),
        ('X', np.array([0.5, 0.5, 0.0])),        # Second X at (0,0,1) conv
        ('W', np.array([0.75, 0.5, 0.25])),
        ('L', np.array([0.5, 0.5, 0.5])),
    ]
    
    # Need special handling for X->G segment
    # Instead of going from (1.0,0.5,0.5) to (0,0,0)
    # Go from (1.0,0.5,0.5) to (0.5,0.5,0) to (0,0,0)
    
    if isinstance(n_points, int):
        n_points_per_segment = [n_points] * (len(path) - 1)
    else:
        n_points_per_segment = n_points
    
    q_path = []
    
    # Generate path with special X->G handling
    for i in range(len(path) - 1):
        start_name, start_pt = path[i]
        end_name, end_pt = path[i + 1]
        n_pts = n_points_per_segment[i]
        
        # Special case: X->G segment (i==2)
        if i == 2:
            # Go via the equivalent X point at (0.5, 0.5, 0)
            X_equiv = np.array([0.5, 0.5, 0.0])
            for j in range(n_pts):
                t = j / n_pts
                q_pt = (1 - t) * start_pt + t * X_equiv
                q_path.append(q_pt)
        else:
            for j in range(n_pts):
                t = j / n_pts
                q_pt = (1 - t) * start_pt + t * end_pt
                q_path.append(q_pt)
    
    # Add final point
    q_path.append(path[-1][1])
    q_path = np.array(q_path)
    
    # Calculate distances
    q_distances = [0.0]
    current_distance = 0.0
    
    for i in range(1, len(q_path)):
        delta = np.linalg.norm(q_path[i] - q_path[i-1])
        current_distance += delta
        q_distances.append(current_distance)
    
    q_distances = np.array(q_distances)
    
    # Create labels
    labels = [(0.0, path[0][0])]
    cumulative_pts = 0
    for i in range(len(path) - 1):
        cumulative_pts += n_points_per_segment[i]
        labels.append((q_distances[cumulative_pts], path[i + 1][0]))
    
    return q_path, q_distances, labels
