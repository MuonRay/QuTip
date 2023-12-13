# -*- coding: utf-8 -*-
def trace_distance_from_angle_list(angle_list: Sequence[float]) -> float:
    """Given a list of arguments of the eigenvalues of a unitary matrix,
    calculates the trace distance bound of the unitary effect.

    The maximum provided angle should not exceed the minimum provided angle
    by more than 2π.
    """
    angles = np.sort(angle_list)
    maxim = 2 * np.pi + angles[0] - angles[-1]
    for i in range(1, len(angles)):
        maxim = max(maxim, angles[i] - angles[i - 1])
    if maxim <= np.pi:
        return 1.0
    return max(0.0, np.sin(0.5 * maxim))