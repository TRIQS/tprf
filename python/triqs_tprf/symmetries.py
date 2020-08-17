import numpy as np

def enforce_symmetry(gf, variables, symmetries):
    """Symmetrize Green's function in the given variables

    Parameters
    ----------
    gf : Gf,
         One-particle fermionic Green's function with a MeshProduct containing
         a MeshImFreq on first and a MeshBrillouinZone in second position.
    variables : str or iterator of str,
                Tells what variable(s) shall be symmetrized, e.g. "momentum"
                or ["frequency", "momentum"]
    symmetries : str or iterator of str,
                 Gives the symmetry for the respective variable, e.g. "even"
                 or ["odd", "even"]

    Returns
    -------
    gf : Gf
    """
    if type(variables) == str:
        variables = [variables]
    if type(symmetries) == str:
        symmetries = [symmetries]
    if len(variables) != len(symmetries):
        raise ValueError("Variables and symmetries must be of equal length.")
    
    variable_symmetrize_fct = {"frequency" : _symmetrize_frequency,
                               "momentum" : _symmetrize_momentum,
                               "orbital" : _symmetrize_orbital,}

    for variable in variables:
        if variable not in list(variable_symmetrize_fct.keys()):
            raise ValueError("No symmetrize function for this variable exists.") 

    for symmetry in symmetries:
        if symmetry not in ['even', 'odd']:
            raise ValueError("Symmetry can only be 'even' or 'odd'.") 

    gf_symmetrized = gf.copy()
    
    for variable, symmetry in zip(variables, symmetries):
        symmetrize_fct = variable_symmetrize_fct[variable]
        symmetrize_fct(gf_symmetrized, symmetry)

    return gf_symmetrized

def check_symmetry(gf, atol=1e-08):
    """Check the symmetry of a Green's function for various variables

    Parameters
    ----------
    gf : Gf,
         One-particle fermionic Green's function with a MeshProduct containing
         a MeshImFreq on first and a MeshBrillouinZone in second position.
    atol : float,
           Absolute tolerance used as parameter `atol` in `np.allclose`.

    Returns
    -------
    variable_symmetry : dict,
                        Keys give the variable and values the symmetry.
                        If even, +1, if odd, -1, else None.
    """
    variable_check_symmetry = {"frequency" : _check_frequency_symmetry,
                               "momentum" : _check_momentum_symmetry,
                               "orbital" : _check_orbital_symmetry,}
    
    variable_symmetry = {}
    for variable, check_symmetry_fct in list(variable_check_symmetry.items()):
        variable_symmetry[variable] = check_symmetry_fct(gf, atol)

    return variable_symmetry

def _average_halfs(half_1, half_2):
    """Stub that can be used as an averager
    """
    return half_2

def _overall_sign(signs):
    """Return +/- 1 if all elements of signs are +/- 1, None else

    Parameters
    ----------
    signs : list of +/- 1
    """
    signs = set(signs)
    if len(signs) == 1:
        return signs.pop()
    return None

# -- Frequency
# ============================================================================
def _split_frequency(gf):
    """Split Green's function data in positive and negative frequencies

    Parameters
    ---------
    gf : Gf,
         Fermionic Green's function with a MeshProduct containing
         a MeshImFreq on first position. 
    
    Returns
    -------
    negative_half : np.array,
    positive_half : np.array,
    """
    if not gf.mesh[0].statistic == 'Fermion':
        raise ValueError("The Green's function must be a fermionic one")

    nw_half = gf.data.shape[0]//2

    negative_half = gf.data[:nw_half]
    positive_half = gf.data[nw_half:]

    return negative_half, positive_half

def _check_frequency_symmetry(gf, atol=1e-08):
    """Check if frequency symmetry of Green's function is even or odd

    Parameters
    ----------
    gf : Gf,
         Fermionic Green's function with a MeshProduct containing
         a MeshImFreq on first position. 
    atol : float,
           Absolute tolerance used as parameter `atol` in `np.allclose`.

    Returns
    -------
    +1 if the Green's function is even in frequency space, -1 if odd,
    and None if undefined.
    """
    negative_half, positive_half = _split_frequency(gf)

    if np.allclose(negative_half[::-1], positive_half, atol=atol):
        return +1
    elif np.allclose(negative_half[::-1], -1*positive_half, atol=atol):
        return -1
    return None

def _symmetrize_frequency(gf, symmetry='even'):
    r"""Symmetrize the data of a Green's function in frequency space

    Parameters
    ---------
    gf : Gf,
         Fermionic Green's function with a MeshProduct containing
         a MeshImFreq on first position. 
    symmetry : str, ['even', 'odd'], optional
               What frequency symmetry shall be enforced:
                   
                   'even' : no sign change :math:`\nu_n\rightarrow\nu_{-n}`
                   'odd'  : sign change :math:`\nu_n\rightarrow\nu_{-n}`
    """
    negative_half, positive_half = _split_frequency(gf)
    avg = _average_halfs(negative_half[::-1], positive_half)

    # Use slice access so that the data in the Green's function get changed
    positive_half[:] = avg
    if symmetry == "even":
        negative_half[:] = avg[::-1]
    else:
        negative_half[:] = -1* avg[::-1]

# -- Momentum
# ============================================================================
def _invert_momentum(momentum, momentum_mesh):
    """Returns the indices corresponding to the inverted momentum

    Parameters
    ----------
    momentum : array_like,
               The indices that correspond to a momentum, e.g. (0, 3).
    momentum_mesh : array_like,
                    The number of points used in the corresponding dimension,
                    e.g. [4, 4].

    Returns
    -------
    inv_k : tuple,
            The indices that correspond to the inverted momentum, 
            e.g. (0, 1).
    """
    momentum = np.array(momentum)
    momentum_mesh = np.array(momentum_mesh)
    inv_k = momentum_mesh - momentum
    inv_k = inv_k % momentum_mesh
    inv_k = tuple(inv_k)
    return inv_k

def _split_momentum(gf):
    """Split Green's function data in momentum and inversed momentum

    Parameters
    ----------
    gf : Gf,
         Green's function with a MeshProduct containing and a 
         MeshBrillouinZone in second position.

    Yields
    ------
    positive_half : np.array
    negative_half : np.array
    """
    nk = gf.data.shape[1]
    momentum_mesh = gf.mesh[1].linear_dims
    # Drop dimensions which are not meshed over, i.e. value of 1
    momentum_mesh = [mesh for mesh in momentum_mesh if mesh != 1]

    for idx in range(nk):
        unraveled_idx = np.unravel_index(idx, momentum_mesh) 
        unraveled_inv_idx = _invert_momentum(unraveled_idx, momentum_mesh)
        inv_idx = np.ravel_multi_index(unraveled_inv_idx,  momentum_mesh)

        positive_half = gf.data[slice(None), idx]
        negative_half = gf.data[slice(None), inv_idx]

        if idx == inv_idx: 
            # Give same id to both halfs for identification of k = -k
            yield positive_half, positive_half
        else:       
            yield positive_half, negative_half

def _check_momentum_symmetry(gf, atol=1e-08):
    """Check if momentum symmetry of Green's function is even or odd

    Parameters
    ----------
    gf : Gf,
         Green's function with a MeshProduct containing and a 
         MeshBrillouinZone in second position.
    atol : float,
           Absolute tolerance used as parameter `atol` in `np.allclose`.

    Returns
    -------
    +1 if the Green's function is even in momentum space, -1 if odd,
    and None if undefined.
    """
    signs= []
    for positive_half, negative_half in _split_momentum(gf):
        
        # If the k-point and its inverse are numercial zero the symmetry does
        # not matter
        if np.allclose(positive_half, 0.0, atol=atol) and \
           np.allclose(negative_half, 0.0, atol=atol):
            continue

        # Check if k = -k, if not equal to 0.0 the gf must be even 
        if id(positive_half) == id(negative_half):
            if not np.allclose(0.0, positive_half):
                signs.append(+1)
            continue

        if np.allclose(positive_half, negative_half, atol=atol):
            signs.append(+1)
        elif np.allclose(-1*positive_half, negative_half, atol=atol):
            signs.append(-1)
        else:
            return None

    return _overall_sign(signs)

def _symmetrize_momentum(gf, symmetry='even'):
    r"""Symmetrize the data of a Green's function in momentum space

    Parameters
    ----------
    gf : Gf,
         Green's function with a MeshProduct containing and a 
         MeshBrillouinZone in second position.
    symmetry : str, ['even', 'odd'], optional
               What momentum symmetry shall be enforced:
                   
                   'even' : no sign change :math:`\mathbf{k}\rightarrow\mathbf{k}`
                   'odd'  : sign change :math:`\mathbf{k}\rightarrow\mathbf{k}`
    """
    for positive_half, negative_half in _split_momentum(gf):
        # Check if k = -k, i.e. needs different symmetry treatment
        if id(positive_half) == id(negative_half):
            if symmetry == "odd":
                positive_half *= 0.0
            continue

        avg = _average_halfs(positive_half, negative_half)

        # Use slice access so that the data in the Green's function get changed
        positive_half[:] = avg
        if symmetry == "even":
            negative_half[:] = avg
        else:
            negative_half[:] = -1 * avg

# -- Orbitals
# ============================================================================
def _split_orbital_triangle(gf):
    """Split Green's function data in upper and lower triangle without diagonal

    Parameters
    ----------
    gf : Gf,
         One-particle Green's function with a MeshProduct with two meshes.
    
    Yields
    ------
    upper_triangle : np.array
    lower_triangle : np.array
    """
    target_shape = gf.target_shape
    nparticle = len(target_shape)
    if nparticle != 2:
        raise ValueError("The Green's function must be a one-particle one.")
    
    norb = target_shape[0]
    mesh_slice = (slice(None),) * gf.mesh.rank

    upper_triangle_slice = np.triu_indices(norb, k=1)
    lower_triangle_slice = np.tril_indices(norb, k=-1)
    combined_slice = upper_triangle_slice + lower_triangle_slice

    for upper_1, upper_2, lower_1, lower_2 in zip(*combined_slice):
        upper_triangle = gf.data[mesh_slice + (upper_1, upper_2)]
        lower_triangle = gf.data[mesh_slice + (lower_1, lower_2)]

        yield upper_triangle, lower_triangle

def _split_orbital_diagonal(gf):
    """Extract the diagonal part of the Green's function data 

    Parameters
    ----------
    gf : Gf,
         One-particle Green's function with a MeshProduct with two meshes.
    
    Yields
    ------
    diagonal : np.array
    """
    target_shape = gf.target_shape
    nparticle = len(target_shape)
    if nparticle != 2:
        raise ValueError("The Green's function must be a one-particle one.")

    norb = target_shape[0]
    mesh_slice = (slice(None),) * gf.mesh.rank

    for orb in range(norb):
        yield gf.data[mesh_slice + (orb,)*nparticle]

def _check_orbital_symmetry(gf, atol=1e-08):
    """Check if orbital symmetry of Green's function is even or odd

    Parameters
    ----------
    gf : Gf,
         One-particle Green's function with a MeshProduct with two meshes.
    atol : float,
           Absolute tolerance used as parameter `atol` in `np.allclose`.

    Returns
    -------
    +1 if the Green's function is even in orbital space, -1 if odd,
    and None if undefined.
    """
    signs = []
    for upper_triangle, lower_triangle in _split_orbital_triangle(gf):
        if np.allclose(upper_triangle, lower_triangle, atol=atol):
            signs.append(+1)
        elif np.allclose(upper_triangle, -1*lower_triangle, atol=atol):
            signs.append(-1)
        else:
            return None

    for diagonal in _split_orbital_diagonal(gf):
        if not np.allclose(diagonal, 0.0, atol=atol):
            signs.append(+1)

    return _overall_sign(signs)
    
def _symmetrize_orbital(gf, symmetry="even"):
    r"""Symmetrize the data of a Green's function in orbital space

    Parameters
    ---------
    gf : Gf,
         One-particle Green's function with a MeshProduct with two meshes.
    """
    for upper_triangle, lower_triangle in _split_orbital_triangle(gf):
        avg = _average_halfs(upper_triangle, lower_triangle)

        # Use slice access so that the data in the Green's function get changed
        upper_triangle[:] = avg
        if symmetry == "even":
            lower_triangle[:] = avg
        else:
            lower_triangle[:] = -1 * avg

    if symmetry == "odd":
        for diagonal in _split_orbital_diagonal(gf):
            diagonal *= 0.0

