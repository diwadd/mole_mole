struct BoundAtom
    id::Int32
    symbol::String
    z_an::Int32 # atomic number
    v::Int32
    cr::Float32
    xyz::Array{Float32}
end

function get_z(symbol::String)

    d = Dict{String,Float32}([("H", 1),
                              ("C", 6),
                              ("N", 7),
                              ("O", 8),
                              ("F", 9)])

    return d[symbol]
end


function get_v(symbol::String)

    d = Dict{String,Int32}([("H", 1),
                            ("C", 4),
                            ("N", 5),
                            ("O", 6),
                            ("F", 7)])

    return d[symbol]
end


function get_cr(symbol::String)

    h_tol = 1.135  # Hydrogen covalent radii tolerance
    c_tol = 1.135
    n_tol = 1.135
    o_tol = 1.135
    f_tol = 1.135

    h_cr = 0.31  # Hydrogen covalent radii
    c_cr = 0.76
    n_cr = 0.71
    o_cr = 0.66
    f_cr = 0.57

    d = Dict{String,Float32}([("H", h_tol*h_cr),
                              ("C", c_tol*c_cr),
                              ("N", n_tol*n_cr),
                              ("O", o_tol*o_cr),
                              ("F", f_tol*f_cr)])

    return d[symbol]
end

function get_atom_string_rep(atom::BoundAtom)

    id = atom.id
    symbol = atom.symbol
    z_an = atom.z_an
    v = atom.v
    cr = atom.cr
    x = atom.xyz[1]
    y = atom.xyz[2]
    z = atom.xyz[3]

    return "id: $id symbol: $symbol z_an: $z_an v: $v cr: $cr x: $x y: $y z: $z"
end
