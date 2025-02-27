# Include packages

using NetworkInference
using LightGraphs

algorithm = PIDCNetworkInference()

dataset_name = string(ARGS[1])

@time genes = get_nodes(dataset_name);

@time network = InferredNetwork(algorithm, genes);

# Define output file path
output_file = string(ARGS[2])

# Write network to a temporary file
temp_file = output_file * ".tmp"
write_network_file(temp_file, network)

# Add header and combine with the network data
open(output_file, "w") do f
    write(f, "TF\ttarget\timportance\n")
    write(f, read(temp_file, String))
end

# Remove the temporary file
rm(temp_file)