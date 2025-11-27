### TKNA script
python ../reconstruction/intake_data.py --data-dir input_folder/ --out-file output/all_data_and_metadata.zip
python ../reconstruction/run.py --data-source output/all_data_and_metadata.zip --config-file config.json --out-file output/network_output.zip
python ../reconstruction/to_csv.py --data-file output/network_output.zip --config-file config.json --out-dir output/network_output
python ../analysis/assess_network.py --file output/network_output/correlations_bw_signif_measurables.csv --out-dir output/network_output/
python ../analysis/infomap_assignment.py --network output/network_output/network_output_comp.csv --network-format csv --map type_map.csv --out-dir output/network_output/
python ../analysis/louvain_partition.py --network output/network_output/network_output_comp.csv --network-format csv --map type_map.csv --out-dir output/network_output/
python ../analysis/find_all_shortest_paths_bw_subnets.py --network output/network_output/network_output_comp.csv --network-format csv --map type_map.csv --node-groups virus bacteria --out-dir output/network_output/
python ../analysis/calc_network_properties.py --network output/network_output/network_output_comp.csv --bibc --bibc-groups node_types --bibc-calc-type rbc --map type_map.csv --node-groups virus bacteria --out-dir output/network_output/
python ../random_networks/create_random_networks.py --template-network output/network_output/network.pickle --networks-file output/network_output/all_random_nws.zip 
python ../random_networks/compute_network_stats.py --networks-file output/network_output/all_random_nws.zip --bibc-groups node_types --bibc-calc-type rbc --stats-file output/network_output/random_network_analysis.zip --node-map type_map.csv --node-groups virus bacteria
python ../random_networks/synthesize_network_stats.py --network-stats-file output/network_output/random_network_analysis.zip --synthesized-stats-file output/network_output/random_networks_synthesized.csv
python ../visualization/dot_plots.py --pickle output/network_output/network.pickle --node-props output/network_output/node_properties.txt --network-file output/network_output/network_output_comp.csv --propx BiBC_virus_bacteria --propy Node_degrees --top-num 5 --top-num-per-type 5 --plot-dir output/network_output/plots/ --file-dir output/network_output/
python ../visualization/plot_density.py --rand-net output/network_output/random_networks_synthesized.csv --pickle output/network_output/inputs_for_downstream_plots.pickle --bibc-name BiBC_virus_bacteria
python ../visualization/plot_abundance.py --pickle output/network_output/inputs_for_downstream_plots.pickle --abund-data data.csv data.csv --metadata data_graph_map.csv data_graph_map.csv --x-axis Experiment --group-names NO_HAP HAP --group-colors blue pink

