# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 21:53:50 2022

@author: cosmi
"""

def _term(rho, j, dims, axis=0):
    # (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I) for all j's
    # in the system we want to trace out.
    # This function returns the jth term in the sum, namely
    # (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I).
    a = sp.coo_matrix(([1.0], ([0], [0])))
    b = sp.coo_matrix(([1.0], ([0], [0])))
    for (i_axis, dim) in enumerate(dims):
        if i_axis == axis:
            v = sp.coo_matrix(([1], ([j], [0])), shape=(dim, 1))
            a = sp.kron(a, v.T)
            b = sp.kron(b, v)
        else:
            eye_mat = sp.eye(dim)
            a = sp.kron(a, eye_mat)
            b = sp.kron(b, eye_mat)
    return a @ rho @ b

def partial_trace(rho, dims, axis=0):
    rho = Expression.cast_to_const(rho)
    if rho.ndim < 2 or rho.shape[0] != rho.shape[1]:
        raise ValueError("Only supports square matrices.")
    if axis < 0 or axis >= len(dims):
        raise ValueError(f"Invalid axis argument, should be between 0 and {len(dims)}, got {axis}.")
    if rho.shape[0] != np.prod(dims):
        raise ValueError("Dimension of system doesn't correspond to dimension of subsystems.")
    return sum([_term(rho, j, dims, axis) for j in range(dims[axis])]