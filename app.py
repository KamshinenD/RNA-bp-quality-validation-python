"""Main entry point for RNA quality scorer."""

import sys
import json
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import Config
from utils.data_loader import DataLoader
from utils.report_generator import ReportGenerator
from analyzers.base_pair_analyzer_bc import BasePairAnalyzer
from analyzers.hbond_analyzer_bc import HBondAnalyzer
from scorer2 import Scorer
import pandas as pd

# Standard amino acid 3-letter codes (normalized to uppercase for comparison)
AMINO_ACIDS_3LETTER = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}


def _is_amino_acid_residue(residue_id: str) -> bool:
    """
    Check if a residue ID corresponds to an amino acid.

    Residue format is like "Q-ARG-52-" where the second part is the residue name.

    Args:
        residue_id: Residue identifier string (e.g., "Q-ARG-52-")

    Returns:
        True if the residue is an amino acid, False otherwise
    """
    try:
        parts = residue_id.split('-')
        if len(parts) >= 2:
            residue_name = parts[1].upper()  # Normalize to uppercase
            return residue_name in AMINO_ACIDS_3LETTER
    except (IndexError, AttributeError):
        pass
    return False


def _get_ligand_name(residue_id: str) -> str:
    """
    Extract the ligand name from a residue ID.

    Args:
        residue_id: Residue identifier string (e.g., "Q-MG-101-")

    Returns:
        The ligand name (e.g., "MG")
    """
    try:
        parts = residue_id.split('-')
        if len(parts) >= 2:
            return parts[1].upper()
    except (IndexError, AttributeError):
        pass
    return "UNKNOWN"


def analyze_protein_bindings(basepair_scores, all_hbond_data, baseline_threshold=75):
    """
    Analyze protein and ligand bindings for problematic base pairs.

    Protein binding: H-bond to at least one amino acid residue.
    Ligand binding: H-bond to non-RNA, non-amino acid molecules (e.g., ions, small molecules).

    Args:
        basepair_scores: List of base pair score dictionaries
        all_hbond_data: DataFrame with ALL H-bonds (RNA-RNA, RNA-PROTEIN, RNA-LIGAND)
        baseline_threshold: Score threshold for problematic pairs

    Returns:
        Tuple of (protein_explanations, ligand_explanations):
            - protein_explanations: Dict mapping bp_id to list of protein binding descriptions
            - ligand_explanations: Dict mapping bp_id to list of ligand binding descriptions
    """
    if all_hbond_data is None or all_hbond_data.empty:
        return {}, {}

    # Filter to RNA-PROTEIN and RNA-LIGAND interactions
    external_bindings = all_hbond_data[
        ((all_hbond_data['res_type_1'] == 'RNA') & (all_hbond_data['res_type_2'].isin(['PROTEIN', 'LIGAND']))) |
        ((all_hbond_data['res_type_1'].isin(['PROTEIN', 'LIGAND'])) & (all_hbond_data['res_type_2'] == 'RNA'))
    ]

    if external_bindings.empty:
        return {}, {}

    # Build separate maps for protein and ligand bindings
    # RNA residue -> list of (external_res, atom_1, atom_2, distance, is_protein)
    rna_protein_map = {}  # H-bonds to amino acids
    rna_ligand_map = {}   # H-bonds to non-amino acid molecules

    for _, row in external_bindings.iterrows():
        rna_res = row['res_1'] if row['res_type_1'] == 'RNA' else row['res_2']
        external_res = row['res_2'] if row['res_type_1'] == 'RNA' else row['res_1']
        atom_1 = row['atom_1']
        atom_2 = row['atom_2']
        distance = row['distance']

        # Determine if external residue is an amino acid (protein) or ligand
        is_amino_acid = _is_amino_acid_residue(external_res)

        # Format binding description
        binding_desc = f"{external_res}({atom_1}-{atom_2}, {distance:.2f}Å)"

        if is_amino_acid:
            # Protein binding (amino acid)
            if rna_res not in rna_protein_map:
                rna_protein_map[rna_res] = []
            rna_protein_map[rna_res].append(binding_desc)
        else:
            # Ligand binding (non-amino acid)
            ligand_name = _get_ligand_name(external_res)
            if rna_res not in rna_ligand_map:
                rna_ligand_map[rna_res] = []
            rna_ligand_map[rna_res].append(f"{ligand_name}:{binding_desc}")

    # Analyze problematic base pairs for both protein and ligand bindings
    protein_explanations = {}
    ligand_explanations = {}

    for bp_score in basepair_scores:
        if bp_score.get('score', 100) < baseline_threshold:
            bp_info = bp_score.get('bp_info', {})
            res_1 = bp_info.get('res_1', '')
            res_2 = bp_info.get('res_2', '')
            bp_id = f"{res_1}-{res_2}"

            # Check for protein bindings
            protein_bindings = []
            if res_1 in rna_protein_map:
                protein_bindings.extend([f"{res_1}:{b}" for b in rna_protein_map[res_1]])
            if res_2 in rna_protein_map:
                protein_bindings.extend([f"{res_2}:{b}" for b in rna_protein_map[res_2]])

            if protein_bindings:
                protein_explanations[bp_id] = protein_bindings

            # Check for ligand bindings
            ligand_bindings = []
            if res_1 in rna_ligand_map:
                ligand_bindings.extend([f"{res_1}:{b}" for b in rna_ligand_map[res_1]])
            if res_2 in rna_ligand_map:
                ligand_bindings.extend([f"{res_2}:{b}" for b in rna_ligand_map[res_2]])

            if ligand_bindings:
                ligand_explanations[bp_id] = ligand_bindings

    return protein_explanations, ligand_explanations


def filter_motif_data(basepair_data, hbond_data, motif_residues=None, start_res=None, end_res=None, chain=None):
    """
    Filter base pairs and H-bonds for a specific motif.
    
    Args:
        basepair_data: List of all base pairs
        hbond_data: DataFrame of all H-bonds
        motif_residues: Set of residue IDs (e.g., {"B-A-74-", "B-C-75-", ...}) - PREFERRED
        start_res: Starting residue number (int) - fallback if motif_residues not provided
        end_res: Ending residue number (int) - fallback if motif_residues not provided
        chain: Optional chain ID filter (str)
        
    Returns:
        Tuple of (filtered_basepairs, filtered_hbonds)
    """
    
    # Filter base pairs
    motif_bps = []
    for bp in basepair_data:
        # Extract residue numbers from res_1 and res_2
        # Format: "X-A-2104-" -> extract chain and residue number
        res1_parts = bp['res_1'].split('-')
        res2_parts = bp['res_2'].split('-')
        
        # Skip malformed residue IDs lacking numeric position
        try:
            res1_chain = res1_parts[0]
            res1_num = int(res1_parts[2])
            res2_num = int(res2_parts[2])
        except (IndexError, ValueError):
            continue
        
        # Check chain filter if specified
        if chain and res1_chain != chain:
            continue
        
        # PREFERRED: Filter by exact residue list from CIF file
        if motif_residues is not None:
            # Both residues must be in the motif residue list
            if bp['res_1'] in motif_residues and bp['res_2'] in motif_residues:
                motif_bps.append(bp)
        # FALLBACK: Filter by range (for backward compatibility)
        elif start_res is not None and end_res is not None:
            if (start_res <= res1_num <= end_res) and (start_res <= res2_num <= end_res):
                motif_bps.append(bp)
    
    # Filter H-bonds
    if motif_residues is not None:
        # Filter by exact residue list
        def residue_in_motif(res_id):
            return res_id in motif_residues
        
        motif_hbonds = hbond_data[
            hbond_data['res_1'].apply(residue_in_motif) &
            hbond_data['res_2'].apply(residue_in_motif)
        ]
    else:
        # Fallback: Filter by range
        def residue_in_range(res_id):
            try:
                parts = res_id.split('-')
                if chain and parts[0] != chain:
                    return False
                res_num = int(parts[2])
                return start_res <= res_num <= end_res
            except (IndexError, ValueError):
                return False
        
        motif_hbonds = hbond_data[
            hbond_data['res_1'].apply(residue_in_range) &
            hbond_data['res_2'].apply(residue_in_range)
        ]
    
    return motif_bps, motif_hbonds


def app():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="RNA Structure Quality Scorer - Baseline Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Score entire structure
  python app.py --pdb_id 1A9N

  # Score specific motif by name (RECOMMENDED - uses exact residues from CIF)
  python app.py --motif-name HAIRPIN-2-CGAG-7O7Y-1

  # Score specific motif (residues 2104-2169, chain AN1)
  python app.py --pdb_id 6V3A --motif 2104 2169 --chain AN1

  # Score motif without chain filter
  python app.py --pdb_id 6V3A --motif 2104 2169

  # Score all base pairs involving a specific residue
  python app.py --pdb_id 1A9N --residue 52

  # Score base pairs for residue with chain filter
  python app.py --pdb_id 6V3A --residue 2104 --chain AN1
  
Output Files:
  - report.json: Detailed baseline quality assessment (entire structure)
  - motif_report.json: Detailed motif quality assessment (when --motif is used)
  - basepair_report.json: Detailed base pair report (when --residue is used)
  - baseline_summary.csv: Summary table (appends/updates)

Data Requirements:
  - Base pair JSON: data/basepairs/{PDB_ID}.json
  - H-bond CSV: data/hbonds/{PDB_ID}.csv

Scoring:
  Each base pair scored individually (0-100).
  Overall score = average of all base pair scores."""
    )
    
    parser.add_argument(
        '--pdb_id',
        required=False,
        help='PDB ID of the RNA structure to analyze'
    )
    
    parser.add_argument(
        '--motif-name',
        type=str,
        help='Name of motif (e.g., HAIRPIN-2-CGAG-7O7Y-1). Automatically extracts CIF file and uses exact residues from CIF for filtering.'
    )
    
    parser.add_argument(
        '--motif',
        nargs=2,
        type=int,
        metavar=('START', 'END'),
        help='Score specific motif by residue range (e.g., --motif 2104 2169)'
    )

    parser.add_argument(
        '--residue',
        type=int,
        metavar='RES_NUM',
        help='Score all base pairs involving a specific residue number (e.g., --residue 52)'
    )
    
    parser.add_argument(
        '--chain',
        type=str,
        help='Chain ID for motif filtering (optional, e.g., --chain AN1)'
    )
    
    parser.add_argument(
        '--motif-residues',
        type=str,
        help='Comma-separated list of residue IDs from CIF file (e.g., "B-A-74-,B-C-75-,B-A-246-")'
    )
    
    parser.add_argument(
        '--csv',
        default='scores_summary.csv',
        help='CSV file for summary output (default: scores_summary.csv)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Directory to save individual motif reports (JSON files). If provided, each motif gets its own report file.'
    )
    
    parser.add_argument(
        '--csv-dir',
        type=str,
        help='Directory to save individual motif CSV files. If provided, each motif gets its own CSV file (avoids race conditions).'
    )
    
    parser.add_argument(
        '--motif-dir',
        type=str,
        default='unique_motifs',
        help='Directory containing motif CIF files (default: unique_motifs)'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.pdb_id and not args.motif_name:
        parser.error("Either --pdb_id or --motif-name must be provided")
    
    # If motif-name is provided, parse CIF file and extract information
    if args.motif_name:
        motif_cif_file = Path(args.motif_dir) / f"{args.motif_name}.cif"
        if not motif_cif_file.exists():
            parser.error(f"Motif CIF file not found: {motif_cif_file}")
        
        # Extract PDB ID from motif name (e.g., HAIRPIN-2-CGAG-7O7Y-1 -> 7O7Y)
        import re
        pdb_match = re.search(r'([0-9][A-Z0-9]{3})', args.motif_name)
        if not pdb_match:
            parser.error(f"Could not extract PDB ID from motif name: {args.motif_name}")
        args.pdb_id = pdb_match.group(1)
        
        # Parse CIF file to extract chain and exact residues
        chain = None
        residues = []
        with open(motif_cif_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    parts = line.split()
                    if len(parts) >= 6:
                        if chain is None:
                            chain = parts[5]  # Chain is in column 6
                        res_num = parts[4]   # Residue number is in column 5
                        base = parts[3]      # Base type is in column 4
                        residue_id = f"{chain}-{base}-{res_num}-"
                        if residue_id not in residues:
                            residues.append(residue_id)
        
        if not chain or not residues:
            parser.error(f"Could not parse chain or residues from CIF file: {motif_cif_file}")
        
        # Set chain and motif-residues automatically
        args.chain = chain
        args.motif_residues = ','.join(sorted(set(residues)))
        
        # Also set motif range for display purposes (min to max residue numbers)
        res_nums = sorted(set([int(r.split('-')[2]) for r in residues]))
        args.motif = [res_nums[0], res_nums[-1]]
        
        print(f"Parsed motif from CIF file:")
        print(f"  PDB ID: {args.pdb_id}")
        print(f"  Chain: {args.chain}")
        print(f"  Residue range: {args.motif[0]}-{args.motif[1]}")
        print(f"  Exact residues: {len(residues)} unique residues from CIF")
    
    # Validate arguments
    if args.chain and not args.motif and not args.motif_name and not args.residue:
        parser.error("--chain can only be used with --motif, --motif-name, or --residue")
    
    # Initialize components
    config = Config()
    data_loader = DataLoader(config)
    bp_analyzer = BasePairAnalyzer(config)
    hb_analyzer = HBondAnalyzer(config)
    report_gen = ReportGenerator(config)
    
    scorer = Scorer(config, bp_analyzer, hb_analyzer)
    
    # Run analysis
    try:
        # Load data
        print(f"Loading data for {args.pdb_id}...")
        basepair_data = data_loader.load_basepairs(args.pdb_id)
        hbond_data = data_loader.load_hbonds(args.pdb_id)  # RNA-RNA only for scoring
        all_hbond_data = data_loader.load_all_hbonds(args.pdb_id)  # All H-bonds for protein binding analysis
        torsion_data = data_loader.load_torsions(args.pdb_id)  # Backbone torsion angles
        
        if basepair_data is None or hbond_data is None:
            print(f"Error: Could not load data for {args.pdb_id}")
            sys.exit(2)
        
        print(f"Loaded {len(basepair_data)} base pairs")
        print(f"Loaded {len(hbond_data)} RNA-RNA H-bonds")
        if all_hbond_data is not None:
            external_count = len(all_hbond_data) - len(hbond_data)
            if external_count > 0:
                print(f"  (Also loaded {external_count} RNA-protein/ligand interactions for binding analysis)")
        
        # Get nucleotide count and validation metrics
        print(f"\nCounting nucleotides for {args.pdb_id}...")
        
        # Get nucleotide count and validation metrics
        # Skip downloads in motif mode - we calculate from range
        if args.motif:
            num_nucleotides = 0  # Don't download - we calculate from range
            validation_metrics = None  # Don't fetch - not needed for motifs, can get from full RNA CSV
            print(f"Skipped (motif mode)")
        else:
            num_nucleotides = data_loader.get_nucleotide_count(args.pdb_id)
            validation_metrics = data_loader.get_validation_metrics(args.pdb_id)
            print(f"Found {num_nucleotides} nucleotides")
            if validation_metrics:
                print(f"Extracted validation metrics: {list(validation_metrics.keys())}")
        
        # MOTIF MODE: Score both entire structure and motif
        if args.motif:
            start_res, end_res = args.motif
            
            # Parse motif residues from CIF file if provided
            motif_residues = None
            if args.motif_residues:
                motif_residues = set(args.motif_residues.split(','))
                print(f"Using {len(motif_residues)} exact residues from CIF file for filtering")
            
            # ========================================
            # STEP 1: Get full structure score (from cache if available)
            # ========================================
            cache_file = Path('full_structure_cache') / f"{args.pdb_id}.json"
            full_score = None
            full_grade = None
            
            if cache_file.exists():
                try:
                    with open(cache_file, 'r') as f:
                        cache_data = json.load(f)
                        full_score = cache_data.get('full_structure_score')
                        full_grade = cache_data.get('full_structure_grade', 'N/A')
                    print(f"\n{'='*60}")
                    print("STEP 1: Using CACHED full structure score")
                    print(f"{'='*60}")
                    print(f"→ Full structure score: {full_score}/100 ({full_grade}) [from cache]")
                except Exception as e:
                    print(f"⚠ Error reading cache: {e}, computing...")
            
            if full_score is None:
                print(f"\n{'='*60}")
                print("STEP 1: Scoring ENTIRE structure for comparison...")
                print(f"{'='*60}")
                
                full_result = scorer.score_structure(basepair_data, hbond_data, torsion_data=torsion_data)
                full_score = full_result.overall_score
                full_grade = full_result.grade
                
                print(f"\n→ Full structure score: {full_score}/100 ({full_grade})")
                
                # Save to cache for future use
                cache_dir = Path('full_structure_cache')
                cache_dir.mkdir(exist_ok=True)
                cache_data = {
                    'pdb_id': args.pdb_id,
                    'full_structure_score': full_score,
                    'full_structure_grade': full_grade,
                    'total_base_pairs': full_result.total_base_pairs,
                    'num_nucleotides': num_nucleotides
                }
                with open(cache_file, 'w') as f:
                    json.dump(cache_data, f, indent=2)
                
                # Convert full result to dictionary
                full_result_dict = scorer.export_to_dict(full_result)
                full_result_dict['pdb_id'] = args.pdb_id
                full_result_dict['analysis_type'] = 'baseline'
                full_result_dict['num_nucleotides'] = num_nucleotides
                if validation_metrics:
                    full_result_dict.update(validation_metrics)
                
                # Save full structure report
                full_output_file = "report.json"
                with open(full_output_file, 'w') as f:
                    json.dump(full_result_dict, f, indent=2)
                print(f"→ Full structure report saved to: {full_output_file}")
            
            # ========================================
            # STEP 2: Score the motif
            # ========================================
            print(f"\n{'='*60}")
            if motif_residues:
                print(f"STEP 2: Scoring MOTIF ({len(motif_residues)} residues from CIF file)")
            else:
                print(f"STEP 2: Scoring MOTIF (Residues {start_res}-{end_res})")
            if args.chain:
                print(f"Chain: {args.chain}")
            print(f"{'='*60}")
            
            motif_basepairs, motif_hbonds = filter_motif_data(
                basepair_data, hbond_data, 
                motif_residues=motif_residues,
                start_res=start_res, 
                end_res=end_res, 
                chain=args.chain
            )
            
            print(f"Filtered to {len(motif_basepairs)} base pairs in motif")
            print(f"Filtered to {len(motif_hbonds)} H-bonds in motif")
            
            if len(motif_basepairs) == 0:
                print("Warning: No base pairs found in specified motif range!")
                sys.exit(1)
            
            # Score the motif
            motif_result = scorer.score_structure(motif_basepairs, motif_hbonds, torsion_data=torsion_data)
            motif_score = motif_result.overall_score
            motif_grade = motif_result.grade
            
            # Convert motif result to dictionary
            temp_motif_dict = scorer.export_to_dict(motif_result)
            
            # Calculate num_problematic_bps from basepair_scores
            # Use BASELINE threshold (75) to match Detailed_Issues column
            num_problematic_bps = sum(
                1 for bp in temp_motif_dict.get('basepair_scores', [])
                if bp['score'] < config.BASELINE
            )
            
            # Count actual unique residues in motif (from base pairs and H-bonds)
            # This handles non-contiguous motifs (e.g., multi-way junctions)
            # NOTE: This is different from the filtering criteria - this counts what was actually found
            actual_motif_residues = set()
            for bp in motif_basepairs:
                actual_motif_residues.add(bp['res_1'])
                actual_motif_residues.add(bp['res_2'])
            for _, hbond in motif_hbonds.iterrows():
                actual_motif_residues.add(hbond['res_1'])
                actual_motif_residues.add(hbond['res_2'])
            
            # Count unique paired nucleotides (only those in base pairs)
            paired_nucleotides = set()
            for bp in motif_basepairs:
                paired_nucleotides.add(bp['res_1'])
                paired_nucleotides.add(bp['res_2'])
            
            # ========================================
            # REORGANIZE: Put important info at TOP
            # ========================================
            motif_result_dict = {
                # CRITICAL INFORMATION FIRST
                'pdb_id': args.pdb_id,
                'analysis_type': 'motif',
                'motif_range': f"{start_res}-{end_res}",
                'motif_chain': args.chain if args.chain else "all",
                
                # COMPARISON METRICS
                'motif_score': motif_score,
                #'motif_grade': motif_grade,
                'full_structure_score': full_score,
                #'full_structure_grade': full_grade,
                'full_structure_num_nucleotides': num_nucleotides,
                # Motif length: actual number of unique residues in motif (handles non-contiguous)
                'motif_num_nucleotides': len(actual_motif_residues),
                # Count unique nucleotides that are paired (appear in at least one base pair)
                'num_paired_nucleotides': len(paired_nucleotides),
                'score_difference': round(motif_score - full_score, 1),
                
                # MOTIF STATISTICS
                'total_base_pairs': motif_result.total_base_pairs,
                'num_problematic_bps': num_problematic_bps,
                
                # DETAILED ANALYSIS BELOW
                'overall_score': motif_result.overall_score,
                #'grade': motif_result.grade,
                'avg_basepair_score': motif_result.avg_basepair_score,
                
                # Issue counts and fractions
                'geometry_issues': temp_motif_dict['geometry_issues'],
                'geometry_fractions': temp_motif_dict['geometry_fractions'],
                'hbond_issues': temp_motif_dict['hbond_issues'],
                'hbond_fractions': temp_motif_dict['hbond_fractions'],
                'summary': temp_motif_dict['summary'],
                
                # Individual base pair details
                'basepair_scores': temp_motif_dict['basepair_scores'],
                
                # Structure-level metadata
            }
            
            # Analyze protein/ligand bindings for problematic base pairs
            protein_bindings, ligand_bindings = analyze_protein_bindings(
                temp_motif_dict.get('basepair_scores', []),
                all_hbond_data,
                baseline_threshold=config.BASELINE
            )
            motif_result_dict['protein_binding_explanations'] = protein_bindings
            motif_result_dict['ligand_binding_explanations'] = ligand_bindings
            
            # Add validation metrics at the end
            if validation_metrics:
                motif_result_dict.update(validation_metrics)
            
            # Get motif name (needed for file naming)
            if args.motif_name:
                motif_name = args.motif_name
            else:
                # Generate motif name from range
                chain_str = f"{args.chain}_" if args.chain else ""
                motif_name = f"{args.pdb_id}_{chain_str}{start_res}-{end_res}"
            
            # Save motif report
            if args.output_dir:
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                motif_output_file = output_dir / f"{motif_name}.json"
            else:
                motif_output_file = Path("motif_report.json")
            
            with open(motif_output_file, 'w') as f:
                json.dump(motif_result_dict, f, indent=2)
            
            # ========================================
            # STEP 3: Print comparison summary
            # ========================================
            print(f"\n{'='*60}")
            print("COMPARISON SUMMARY")
            print(f"{'='*60}")
            print(f"Full Structure: {full_score}/100 ({full_grade})")
            print(f"Motif Score:    {motif_score}/100 ({motif_grade})")
            print(f"Difference:     {motif_score - full_score:+.1f} points")
            
            if motif_score < full_score:
                print(f"→ Motif is {full_score - motif_score:.1f} points WORSE than full structure")
            elif motif_score > full_score:
                print(f"→ Motif is {motif_score - full_score:.1f} points BETTER than full structure")
            else:
                print(f"→ Motif score matches full structure")
            
            print(f"\n{'='*60}")
            #print(f"Full structure report: {full_output_file}")
            print(f"Motif report:          {motif_output_file}")
            print(f"{'='*60}\n")
            
            # Save/update CSV summary for motif (separate CSV file)
            
            report_gen.save_motifs_summary_csv(
                motif_result_dict, 
                motif_name=motif_name,
                csv_dir=args.csv_dir
            )
            
            # Exit code based on motif quality
            result = motif_result
            
        # SINGLE RESIDUE MODE: Score all base pairs involving a specific residue
        elif args.residue:
            residue_num = args.residue

            print(f"\n{'='*60}")
            print(f"SINGLE RESIDUE MODE: Finding base pairs for residue {residue_num}")
            if args.chain:
                print(f"Chain filter: {args.chain}")
            print(f"{'='*60}")

            # Find all base pairs that involve this residue
            matching_bps = []
            for bp in basepair_data:
                res1_parts = bp['res_1'].split('-')
                res2_parts = bp['res_2'].split('-')

                try:
                    res1_chain = res1_parts[0]
                    res1_num = int(res1_parts[2])
                    res2_chain = res2_parts[0]
                    res2_num = int(res2_parts[2])
                except (IndexError, ValueError):
                    continue

                # Check if this base pair involves our residue
                residue_match = (res1_num == residue_num or res2_num == residue_num)

                # Apply chain filter if specified
                if args.chain:
                    chain_match = (
                        (res1_num == residue_num and res1_chain == args.chain) or
                        (res2_num == residue_num and res2_chain == args.chain)
                    )
                    if residue_match and chain_match:
                        matching_bps.append(bp)
                elif residue_match:
                    matching_bps.append(bp)

            if len(matching_bps) == 0:
                print(f"Error: No base pairs found involving residue {residue_num}")
                if args.chain:
                    print(f"  (with chain filter: {args.chain})")
                sys.exit(1)

            print(f"Found {len(matching_bps)} base pair(s) involving residue {residue_num}")

            # Collect residue IDs for H-bond filtering
            residue_ids = set()
            for bp in matching_bps:
                residue_ids.add(bp['res_1'])
                residue_ids.add(bp['res_2'])

            # Filter H-bonds to those between residues in our base pairs
            filtered_hbonds = hbond_data[
                hbond_data['res_1'].isin(residue_ids) &
                hbond_data['res_2'].isin(residue_ids)
            ]

            print(f"Found {len(filtered_hbonds)} H-bonds for these base pairs")

            # Score each base pair individually
            bp_results = []
            for bp in matching_bps:
                bp_score_dict = scorer._score_base_pair(bp, filtered_hbonds)

                # Add full geometry parameters for detailed report
                bp_score_dict['geometry_params'] = {
                    'shear': bp.get('shear', 0),
                    'stretch': bp.get('stretch', 0),
                    'stagger': bp.get('stagger', 0),
                    'buckle': bp.get('buckle', 0),
                    'propeller': bp.get('propeller', 0),
                    'opening': bp.get('opening', 0),
                }

                # Get H-bonds for this specific base pair
                nt1_id = bp.get('res_1', '')
                nt2_id = bp.get('res_2', '')
                bp_hbonds = scorer._get_basepair_hbonds(nt1_id, nt2_id, filtered_hbonds)

                # Add detailed H-bond info
                hbond_details = []
                for _, hb in bp_hbonds.iterrows():
                    hbond_details.append({
                        'atom_1': hb.get('atom_1', ''),
                        'atom_2': hb.get('atom_2', ''),
                        'distance': round(hb.get('distance', 0), 3),
                        'angle_1': round(hb.get('angle_1', 0), 1),
                        'angle_2': round(hb.get('angle_2', 0), 1),
                        'dihedral_angle': round(hb.get('dihedral_angle', 0), 1),
                        'quality_score': round(hb.get('score', 0), 3),
                    })
                bp_score_dict['hbond_details'] = hbond_details

                bp_results.append(bp_score_dict)

            # Calculate summary statistics
            total_score = sum(bp['score'] for bp in bp_results)
            avg_score = total_score / len(bp_results) if bp_results else 0

            # Determine grade
            if avg_score >= config.GRADE_EXCELLENT:
                grade = "EXCELLENT"
            elif avg_score >= config.GRADE_GOOD:
                grade = "GOOD"
            elif avg_score >= config.GRADE_FAIR:
                grade = "FAIR"
            else:
                grade = "POOR"

            # Build the report
            report = {
                'pdb_id': args.pdb_id,
                'analysis_type': 'single_residue',
                'query_residue': residue_num,
                'query_chain': args.chain if args.chain else 'all',

                # Summary
                'num_base_pairs': len(bp_results),
                'grade': grade,

                # Individual base pair scores
                'base_pairs': bp_results,
            }

            # Analyze protein/ligand bindings
            protein_bindings, ligand_bindings = analyze_protein_bindings(
                bp_results,
                all_hbond_data,
                baseline_threshold=config.BASELINE
            )
            report['protein_binding_explanations'] = protein_bindings
            report['ligand_binding_explanations'] = ligand_bindings

            # Save report
            output_file = Path("basepair_report.json")
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)

            # Print summary
            print(f"\n{'='*60}")
            print("BASE PAIR REPORT")
            print(f"{'='*60}")
            print(f"Query: Residue {residue_num}" + (f" (Chain {args.chain})" if args.chain else ""))
            print(f"Base pairs found: {len(bp_results)}")
            print(f"Average score: {avg_score:.1f}/100 ({grade})")
            print(f"\nIndividual base pairs:")
            for bp in bp_results:
                info = bp['bp_info']
                score = bp['score']
                bp_type = info.get('bp_type', 'unknown')
                edge = info.get('edge_type', 'unknown')
                print(f"  {info['res_1']} <-> {info['res_2']}")
                print(f"    Type: {bp_type} ({edge}), Score: {score}/100")
                if bp['geometry_issues']:
                    print(f"    Geometry issues: {list(bp['geometry_issues'].keys())}")
                if bp['hbond_issues']:
                    print(f"    H-bond issues: {list(bp['hbond_issues'].keys())}")

            print(f"\n{'='*60}")
            print(f"Detailed report saved to: {output_file}")
            print(f"{'='*60}\n")

            sys.exit(0)

        # FULL STRUCTURE MODE: Score only entire structure
        else:
            print(f"\n{'='*60}")
            print("Scoring ENTIRE structure...")
            print(f"{'='*60}")
            
            result = scorer.score_structure(basepair_data, hbond_data, torsion_data=torsion_data)
            
            # Convert result to dictionary for JSON export
            result_dict = scorer.export_to_dict(result)
            result_dict['pdb_id'] = args.pdb_id
            result_dict['analysis_type'] = 'baseline'
            result_dict['num_nucleotides'] = num_nucleotides
            
            # If no base pairs, mark scores as N/A instead of treating as failure
            if result.total_base_pairs == 0:
                result_dict['overall_score'] = 'N/A'
                result_dict['grade'] = 'N/A'
                result_dict['avg_basepair_score'] = 'N/A'
            
            # Analyze protein/ligand bindings and add directly to base pair objects in JSON
            protein_bindings, ligand_bindings = analyze_protein_bindings(
                result_dict.get('basepair_scores', []),
                all_hbond_data,
                baseline_threshold=config.BASELINE
            )

            # Add protein binding info directly to each base pair object
            for bp_score in result_dict.get('basepair_scores', []):
                bp_info = bp_score.get('bp_info', {})
                res_1 = bp_info.get('res_1', '')
                res_2 = bp_info.get('res_2', '')
                bp_id = f"{res_1}-{res_2}"

                # Add bindings to this base pair if it has any
                if bp_id in protein_bindings:
                    bp_score['protein_bindings'] = protein_bindings[bp_id]
                else:
                    bp_score['protein_bindings'] = []

                # Add ligand bindings to this base pair if it has any
                if bp_id in ligand_bindings:
                    bp_score['ligand_bindings'] = ligand_bindings[bp_id]
                else:
                    bp_score['ligand_bindings'] = []

            # Keep top-level for backward compatibility (can be removed later if not needed)
            result_dict['protein_binding_explanations'] = protein_bindings
            result_dict['ligand_binding_explanations'] = ligand_bindings
            
            if validation_metrics:
                result_dict.update(validation_metrics)
            
            # Save detailed JSON report
            output_file = "report.json"
            with open(output_file, 'w') as f:
                json.dump(result_dict, f, indent=2)
            
            print(f"\nDetailed report saved to: {output_file}")
            
            # Save/update CSV summary
            report_gen.save_score_summary_csv(
                result_dict, args.csv, 
                hbond_data=hbond_data, 
                validation_metrics=validation_metrics
            )

        # Always exit success unless an exception occurred
        sys.exit(0)
            
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(2)

if __name__ == "__main__":
    app()
    
    
    
    # python app.py --pdb_id 8IFE
    # to run the app on full structure

    # python main.py --pdb_id 6V3A --motif 2104 2169 --chain AN1
    # to run the app on motif with chain filter

    # python app.py --pdb_id 6V3A --motif-name HAIRPIN-2-CGAG-7O7Y-1
    # to run the app on motif with name

    # python app.py --pdb_id 1A9N --residue 52
    # to score all base pairs involving residue 52

    # python app.py --pdb_id 6V3A --residue 2104 --chain AN1
    # to score base pairs for residue 2104 in chain AN1
    
    # Score all base pairs involving residue 52                                                                
    #python app.py --pdb_id 1A9N --residue 52                                                                   
                                                                                                             
    # Score base pairs for residue 2104 in chain AN1                                                           
    #python app.py --pdb_id 6V3A --residue 2104 --chain AN1   
    