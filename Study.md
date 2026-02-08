====================================================================================================
UNDERSTANDING CIS vs TRANS DIHEDRAL ANGLES
====================================================================================================

1. WHAT IS A DIHEDRAL ANGLE?
----------------------------------------------------------------------------------------------------
The dihedral angle describes the orientation of the H-bond relative to the base pair.
It's measured in degrees, typically ranging from -180° to +180°.

The angle tells us how the H-bond is oriented:
  - Around 0°: H-bond is in 'CIS' orientation (same side)
  - Around ±180°: H-bond is in 'TRANS' orientation (opposite side)

2. HOW DO WE CLASSIFY CIS vs TRANS?
----------------------------------------------------------------------------------------------------
We classify based on the ANGLE VALUE itself, not the edge type notation!

Classification rules:
  CIS range: -50° to +50°
    → H-bonds with dihedral angles in this range are 'CIS'
    → Example: -30°, 0°, +25°, +45° are all CIS

  TRANS range: >= 140° OR <= -140°
    → H-bonds with dihedral angles in this range are 'TRANS'
    → Example: 150°, 170°, -150°, -170° are all TRANS

  Forbidden zone: 50° to 140° and -50° to -140°
    → These are strained/unusual orientations
    → Example: 60°, 100°, -80° are in forbidden zone


    Commands to run for all RNA's
    -python3 cache_metadata_parallel.py to cache all metadata like nucleotide count and validation metrics
    -python3 run_all_rnas_fast.py for all rnas in data folder or python run_all_unique_rna_fast.py for unique rnas

    sbatch run_all_rnas_cluster.sh

    Commands to run for all motifs
    python3 cache_all_unique_rnas.py
    sbatch process_motifs_cluster.sh
    python3 merge_motif_csvs.py --csv-dir motif_csvs --output scores_motifs_summary.csv

    To score basepairs and write results for all motifs base-pairs in csv
      sbatch run_export_motif_basepairs_cluster.sh
      python3 export_motif_basepairs.py --merge-shards --shard-dir motif_base_pair --output motif_basepairs.csv

    To do the above locally
    python3 cache_all_unique_rnas.py
    python3 export_motif_basepairs.py --motifs-dir unique_motifs --shard-by-pdb --shard-dir motif_base_pair
    python3 export_motif_basepairs.py --merge-shards --shard-dir motif_base_pair --output motif_basepairs.csv




# Updated Use default 'unique_motifs' folder =  python app.py --pdb_id 6V3A --motif-name HAIRPIN-2-CGAG-7O7Y-1

# Use custom folder = python app.py --pdb_id 6V3A --motif-name HAIRPIN-2-CGAG-7O7Y-1 --motif-dir /path/to/my/motifs

# Use relative custom folder = python app.py --pdb_id 6V3A --motif-name HAIRPIN-2-CGAG-7O7Y-1 --motif-dir custom_motifs


Usage:                                                                                                                                       
                                                                                                                                               
  # Default: 500 EM + 500 X-ray                                                                                                               
  python torsion_scores_analysis.py                                                                                                            
                                                                                                                                               
  #Custom sampling                                                                                                                    python torsion_scores_analysis.py --em 1000 --xray 1000                                                                                      
                                                                                                   
  # Process all base-pairs                                                                                                                  python torsion_scores_analysis.py --all                                                                                                      
                                                                                                                                               
  # Specify output file                                                                                                                      
  python torsion_scores_analysis.py -o my_analysis.csv


//for config geom and hbon thresholds generation
  python analyze_by_edge_type.py
  python generate_threshold_config.py 
      