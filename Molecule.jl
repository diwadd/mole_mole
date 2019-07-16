struct Molecule
    molecule_name::String
    atoms::Array{BoundAtom}
    molecule_graph::Array{Array{UInt32}}
end


function print_atoms(atoms::Array{BoundAtom})

    n = length(atoms)
    println("  --- Printing atoms --- ")
    println(string("Number of atoms: ",n))
    for i in 1:1:n
        s = get_atom_string_rep(atoms[i])
        println(s)
    end
end


function make_molecule_graphy(atoms::Array{BoundAtom})

    sorted_atoms = sort(atoms, by = a -> a.z_an, rev=true)
    # println(" --- After sort, inside make_molecule_graphy --- ")
    # print_atoms(sorted_atoms)

    n = length(sorted_atoms)
    # molecule_graph = Array{Array{UInt32,1},1}([[] for i=1:n])
    molecule_graph = Array{Array{UInt32,1},1}([Array{UInt32,1}(undef,0) for i=1:n])

    for i = 1:n
        a = sorted_atoms[i]

        for j = i+1:n
            b = sorted_atoms[j]
            r = LinearAlgebra.norm(a.xyz - b.xyz)
            max_bond = a.cr + b.cr

            if r <= max_bond

                if length(molecule_graph[a.id]) < a.v && length(molecule_graph[b.id]) < b.v
                    push!(molecule_graph[a.id], b.id)
                    push!(molecule_graph[b.id], a.id)
                end

                # if length(molecule_graph[b.id]) < b.v
                #     push!(molecule_graph[b.id], a.id)
                # end
            end
        end
    end

    is_graph_valid = true
    for i = 1:n
        if length(molecule_graph[i]) == 0
            s = get_atom_string_rep(atoms[i])
            msg = string("Unbount atom ", s)
            #error(error_msg)
            println(msg)
            is_graph_valid = false
        end
    end

    return is_graph_valid, molecule_graph
end


function get_number_of_edges(molecule_graph::Array{Array{UInt32}})

    n_edges = 0
    n_nodes = length(molecule_graph)
    for i = 1:n_nodes
        n_edges = n_edges + length(molecule_graph[i])
    end

    return n_edges
end



function print_molecular_graph(molecule)

    println(" --- Printing molecule graph --- ")
    n = length(molecule.atoms)
    for i = 1:n
        println(molecule.atoms[i])
        m = length(molecule.molecule_graph[i])
        for j = 1:m
            neighbour_id = molecule.molecule_graph[i][j]
            neighbour_atom = molecule.atoms[neighbour_id]

            s = string(" --> ", get_atom_string_rep(neighbour_atom))
            println(s)
        end
    end

end
