using HDF5

h5file = "data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"

println("Examining camera data in: ", basename(h5file))
println("="^80)

HDF5.h5open(h5file, "r") do file
    camera_group = file["Main"]["camera"]

    println("\n📷 Camera Group Contents:")
    for name in keys(camera_group)
        obj = camera_group[name]
        if obj isa HDF5.Dataset
            println("\n  Dataset: $name")
            println("    Size: ", size(obj))
            println("    Type: ", eltype(obj))

            # Show attributes
            attrs = HDF5.attributes(obj)
            if length(keys(attrs)) > 0
                println("    Attributes:")
                for attr_name in keys(attrs)
                    attr_val = read(attrs[attr_name])
                    println("      - $attr_name: $attr_val")
                end
            end

            # Show sample data for smaller datasets
            if prod(size(obj)) <= 100
                println("    Data: ", read(obj))
            elseif name == "data"
                # For the main data, just show some stats
                println("    Shape: (width=$(size(obj,1)), height=$(size(obj,2)), frames=$(size(obj,3)))")
                println("    Total frames: ", size(obj,3))
                println("    Sample frame (first 5x5 of frame 1):")
                sample = read(obj[1:5, 1:5, 1])
                display(sample)
            end
        elseif obj isa HDF5.Group
            println("\n  Group: $name")
            for subname in keys(obj)
                subobj = obj[subname]
                if subobj isa HDF5.Dataset
                    println("    - $subname: ", size(subobj), " (", eltype(subobj), ")")
                    # Show small datasets
                    if prod(size(subobj)) <= 20
                        println("      Data: ", read(subobj))
                    end
                end
            end
        end
    end

    # Check for other useful groups
    println("\n\n📋 Other Main Groups:")
    for name in keys(file["Main"])
        if name != "camera"
            obj = file["Main"][name]
            if obj isa HDF5.Group
                println("  Group: $name")
                for subname in keys(obj)
                    println("    - $subname")
                end
            end
        end
    end
end

println("\n" * "="^80)
