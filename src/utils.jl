"""
    all_bound(atom)

Returns `true` if all orbitals in `atom` are bound orbitals.
"""
all_bound(atom::Atom) = all(isbound.(atom.orbitals))
