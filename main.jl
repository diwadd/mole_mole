include("BoundAtom.jl")
include("Molecule.jl")
include("ScalarCouplingConstantDataPoint.jl")

using Pkg
Pkg.add("Plots")

using Plots
using LinearAlgebra

TRAIN_CVS_FILE_NAME = "../champs-scalar-coupling/train.csv"
STRUCTURE_XYZ_FILES_PATH = "../champs-scalar-coupling/structures/"


function read_xyz(file_name::String)

    open(file_name) do file
        lines = readlines(file)
        # println(lines[1])

        n = length(lines)
        atoms = Array{BoundAtom}(undef, n-2)
        for i = 3:n
            # println(lines[i])
            l = split(lines[i], " ")

            symbol = String(l[1])
            z_an = get_z(symbol)
            v = get_v(symbol)
            cr = get_cr(symbol)

            x = parse(Float32, l[2])
            y = parse(Float32, l[3])
            z = parse(Float32, l[4])
            xyz = [x, y, z]

            atoms[i-2] = BoundAtom(i-2, symbol, z_an, v, cr, xyz)

        end

        return atoms
    end
end

function check_if_dict_has_key(d, k)
    try
        d[k]
        catch error
           if isa(error, KeyError)
               return false
           end
    end
    return true
end

function read_data(file_name::String,
                   structures::String=STRUCTURE_XYZ_FILES_PATH)
    @time begin

        println(" --- Reading data --- ")

        molecules = Dict{String, Molecule}()
        not_valid_molecules = Dict{String, Molecule}()

        open(file_name) do file
            lines = readlines(file)

            n = length(lines)
            data = Array{ScalarCouplingConstantDataPoint, 1}(undef, n)
            for i = 2:n

                if i % 10000 == 0
                    p = i/n
                    println("We are at $i out of $n which is $p")
                end

                l = split(lines[i], ",")
                # print(string(l) * "\n")

                id = parse(Int32, l[1])
                molecule_name = String(l[2])
                atom_one = parse(Int32, l[3])
                atom_two = parse(Int32, l[4])
                type = String(l[5])
                scalar_coupling_constant = 0.0

                if size(l) == 6
                    scalar_coupling_constant = parse(Float32, l[6])
                end

                if check_if_dict_has_key(molecules, molecule_name) == false &&
                   check_if_dict_has_key(not_valid_molecules, molecule_name) == false
                    atoms = read_xyz(structures * molecule_name * ".xyz")

                    is_graph_valid, molecule_graph = make_molecule_graphy(atoms)
                    molecule = Molecule(molecule_name, atoms, molecule_graph)

                    # print_atoms(atoms)
                    # print_molecular_graph(molecule)
                    if is_graph_valid == true
                        molecules[molecule_name] = molecule
                    else
                        not_valid_molecules[molecule_name] = molecule
                    end
                end

                data[i-1] = ScalarCouplingConstantDataPoint(id, molecule_name, atom_one, atom_two, type, scalar_coupling_constant)

            end
            n_molecules = length(molecules)
            n_not_valid_molecules = length(not_valid_molecules)

            println("Number of valid molecules $n_molecules")
            println("Number of not valid molecules $n_not_valid_molecules")
            return data, molecules
        end

    end
end

data, molecules = read_data(TRAIN_CVS_FILE_NAME)
#read_xyz(STRUCTURE_XYZ_FILES_PATH * "dsgdb9nsd_109615.xyz")
