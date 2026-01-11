using HDF5

h5file = "data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"

println("Examining file: ", basename(h5file))
println("="^80)

function explore_h5(file, prefix="")
    for name in keys(file)
        obj = file[name]
        full_path = isempty(prefix) ? name : "$prefix/$name"

        if obj isa HDF5.Group
            println("\n$(prefix)📁 Group: $name")
            explore_h5(obj, "  $prefix")
        elseif obj isa HDF5.Dataset
            ds_size = size(obj)
            ds_type = eltype(obj)
            println("$(prefix)📄 Dataset: $name")
            println("$(prefix)   Size: $ds_size")
            println("$(prefix)   Type: $ds_type")

            # Show attributes if any
            attrs = HDF5.attributes(obj)
            if length(keys(attrs)) > 0
                println("$(prefix)   Attributes:")
                for attr_name in keys(attrs)
                    attr_val = read(attrs[attr_name])
                    println("$(prefix)     - $attr_name: $attr_val")
                end
            end

            # Show a sample of the data if it's small enough
            if prod(ds_size) <= 10
                println("$(prefix)   Data: ", read(obj))
            elseif length(ds_size) == 1 && ds_size[1] <= 20
                println("$(prefix)   Data: ", read(obj))
            end
        end
    end
end

HDF5.h5open(h5file, "r") do file
    # Show root attributes
    root_attrs = HDF5.attributes(file)
    if length(keys(root_attrs)) > 0
        println("\n🏷️  Root Attributes:")
        for attr_name in keys(root_attrs)
            attr_val = read(root_attrs[attr_name])
            println("  - $attr_name: $attr_val")
        end
    end

    println("\n📂 File Structure:")
    explore_h5(file)
end

println("\n" * "="^80)
