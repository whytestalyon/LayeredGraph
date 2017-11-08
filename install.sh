#!/usr/bin/env bash
pip3 install flask flask-restful requests numpy pronto six scipy xmltodict

if [ ! -d ./HPO_graph_data ]; then
    echo "Preparing to download graph data, please enter Morgan username:"
    read USERNAME
    echo "Downloading graph data from Morgan, scp will prompt for password..."
    scp -rp ${USERNAME}@login.morgan.haib.org:/gpfs/gpfs1/home/jholt/HPO_graph_data .
fi

