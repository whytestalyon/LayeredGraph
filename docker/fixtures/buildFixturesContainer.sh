#!/usr/bin/env bash
echo "Preparing to download graph data, please enter Morgan username:"
read USERNAME
echo "Downloading graph data from Morgan, scp will prompt for password..."
scp -rp ${USERNAME}@login.morgan.haib.org:/gpfs/gpfs1/home/jholt/HPO_graph_data .

docker build -t docker-registry.haib.org/sdi/layered-graph-fixtures .
if [ $? -ne 0 ]; then
    echo "Failed to build the fixtures container..."
else
    echo "Fixtures image built!"
fi

echo "Cleaning up..."
rm -rf HPO_graph_data

