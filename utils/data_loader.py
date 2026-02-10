"""Utility functions for loading base pair and H-bond data."""

import json
import pandas as pd
from pathlib import Path
from typing import Optional
import requests


class DataLoader:
    """Loads precomputed RNA structural data."""
    
    def __init__(self, config):
        self.config = config
        self.basepair_dir = Path(config.BASEPAIR_DIR)
        self.hbond_dir = Path(config.HBOND_DIR)
    
    # def load_basepairs(self, pdb_id: str) -> Optional[list]:
    #     """Load base pair data from JSON file."""
    #     file_path = self.basepair_dir / f"{pdb_id.upper()}.json"
        
    #     if not file_path.exists():
    #         print(f"Warning: Base pair file not found: {file_path}")
    #         return None
        
    #     try:
    #         with open(file_path, 'r') as f:
    #             data = json.load(f)
    #         print(f"✓ Loaded {len(data)} base pairs from {file_path.name}")
    #         return data
    #     except Exception as e:
    #         print(f"Error loading base pairs: {e}")
    #         return None
    
    
    def load_basepairs(self, pdb_id: str, quiet: bool = False) -> list:
        """
        Load base pair data from JSON file.
        Filters out stacking interactions (adjacent residues).
        """
        bp_file = self.basepair_dir / f"{pdb_id}.json"
        
        if not bp_file.exists():
            if not quiet:
                print(f"Error: Base pair file not found: {bp_file}")
            return None
        
        try:
            with open(bp_file, 'r') as f:
                data = json.load(f)
            
            # data is already a list, not a dict
            # Handle both formats: list or dict with 'base_pairs' key
            if isinstance(data, list):
                base_pairs = data
            elif isinstance(data, dict):
                base_pairs = data.get('base_pairs', [])
            else:
                if not quiet:
                    print(f"Error: Unexpected data format in {bp_file}")
                return None
            
            # Filter out stacking interactions (adjacent residues)
            filtered_bps = []
            stacking_count = 0
            
            for bp in base_pairs:
                res1_parts = bp['res_1'].split('-')
                res2_parts = bp['res_2'].split('-')
                
                try:
                    res1_num = int(res1_parts[2])
                    res2_num = int(res2_parts[2])
                    
                    # Skip if residues are adjacent (stacking)
                    if abs(res1_num - res2_num) <= 1:
                        stacking_count += 1
                        continue
                    
                    filtered_bps.append(bp)
                    
                except (IndexError, ValueError):
                    # If parsing fails, keep the base pair
                    filtered_bps.append(bp)
            
            if stacking_count > 0 and not quiet:
                print(f"  Filtered out {stacking_count} stacking interactions")
            
            return filtered_bps
            
        except Exception as e:
            if not quiet:
                print(f"Error loading base pairs: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def load_hbonds(self, pdb_id: str, quiet: bool = False) -> Optional[pd.DataFrame]:
        """Load H-bond data from CSV file (RNA-RNA only)."""
        file_path = self.hbond_dir / f"{pdb_id.upper()}.csv"
        # json_files = list(self.hbond_dir.glob("*.csv"))
        # print(f"@@@@@@@@@@@@@@@@@@@@@@@@@Found {len(json_files)} CSV files in {self.hbond_dir}")
        
        if not file_path.exists():
            if not quiet:
                print(f"Warning: H-bond file not found: {file_path}")
            return None
        
        try:
            df = pd.read_csv(file_path)
            
            # Filter to only RNA-RNA interactions
            rna_rna = df[(df['res_type_1'] == 'RNA') & (df['res_type_2'] == 'RNA')]
            
            if not quiet:
                print(f"✓ Loaded {len(rna_rna)} RNA-RNA H-bonds from {file_path.name}")
                print(f"  (Filtered out {len(df) - len(rna_rna)} RNA-protein interactions)")
            
            return rna_rna
        except Exception as e:
            if not quiet:
                print(f"Error loading H-bonds from {file_path}: {e}")
            return None
    
    def load_all_hbonds(self, pdb_id: str, quiet: bool = False) -> Optional[pd.DataFrame]:
        """Load all H-bond data from CSV file (including RNA-PROTEIN and RNA-LIGAND)."""
        file_path = self.hbond_dir / f"{pdb_id.upper()}.csv"
        
        if not file_path.exists():
            if not quiet:
                print(f"Warning: H-bond file not found: {file_path}")
            return None
        
        try:
            df = pd.read_csv(file_path)
            if not quiet:
                print(f"✓ Loaded {len(df)} total H-bonds from {file_path.name}")
            return df
        except Exception as e:
            if not quiet:
                print(f"Error loading H-bonds: {e}")
            return None
        
    def load_torsions(self, pdb_id: str, quiet: bool = False) -> dict:
        """Load per-residue torsion angles from JSON file.

        Returns:
            Dict keyed by residue ID (e.g. 'A-C-1-') with torsion angle values,
            or None if file not found.
        """
        file_path = Path('data/torsions') / f"{pdb_id.upper()}.json"

        if not file_path.exists():
            if not quiet:
                print(f"Warning: Torsion file not found: {file_path}")
            return None

        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            if not quiet:
                print(f"✓ Loaded torsion data for {len(data)} residues from {file_path.name}")
            return data
        except Exception as e:
            if not quiet:
                print(f"Error loading torsion data: {e}")
            return None

    def download_cif(self, pdb_id: str, output_path: str = "temp_structure.cif") -> bool:
        """
        Download CIF file from RCSB PDB.
        
        Args:
            pdb_id: PDB ID (e.g., '3QG9')
            output_path: Where to save the CIF file
            
        Returns:
            True if successful, False otherwise
        """
        pdb_id = pdb_id.upper().strip()
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                with open(output_path, 'w') as f:
                    f.write(response.text)
                return True
            else:
                print(f"Warning: Could not download {pdb_id}.cif (HTTP {response.status_code})")
                return False
                
        except Exception as e:
            print(f"Warning: Error downloading {pdb_id}.cif: {e}")
            return False
    
    """
CORRECTED CIF PARSER for utils/data_loader.py

Replace the count_nucleotides_from_cif method with this one.
"""

    def count_nucleotides_from_cif(self, cif_file: str) -> int:
        """
        Count unique RNA nucleotides from CIF file.
        
        CIF format uses columnar data after headers like:
        _atom_site.group_PDB
        _atom_site.label_comp_id (residue name)
        _atom_site.auth_asym_id (chain)
        _atom_site.auth_seq_id (residue number)
        
        Args:
            cif_file: Path to CIF file
            
        Returns:
            Number of unique nucleotides
        """
        rna_residues = {
            # Standard RNA
            'A', 'C', 'G', 'U'
            # Common modified
            'PSU', 'H2U', '5MU', '4SU', '1MA', 'M2G', 'OMC', 'OMG',
        }
        
        nucleotides = set()
        in_atom_site = False
        
        try:
            with open(cif_file, 'r') as f:
                for line in f:
                    # Check if we're in the atom_site loop
                    if line.startswith('_atom_site.'):
                        in_atom_site = True
                        continue
                    
                    # End of atom_site loop
                    if in_atom_site and (line.startswith('#') or line.startswith('loop_') or 
                                        (line.startswith('_') and not line.startswith('_atom_site.'))):
                        in_atom_site = False
                        continue
                    
                    # Parse atom lines in the atom_site loop
                    if in_atom_site and line.startswith(('ATOM', 'HETATM')):
                        parts = line.split()
                        
                        # CIF format (columns may vary, but typically):
                        # 0: ATOM/HETATM
                        # 1: atom_id
                        # 2: type_symbol
                        # 3: label_atom_id
                        # 4: alt_id (. or A/B/C)
                        # 5: label_comp_id (residue name) ← WE NEED THIS
                        # 6: label_asym_id
                        # 7: label_entity_id
                        # 8: label_seq_id
                        # ... (more fields)
                        # Near end: auth_seq_id, auth_comp_id, auth_asym_id
                        
                        if len(parts) >= 20:  # CIF has many columns
                            try:
                                # Get residue name (label_comp_id - column 5)
                                res_name = parts[5]
                                
                                # Get chain (auth_asym_id - usually near end)
                                # In the sample: position varies, but auth_asym_id is typically ~column 18-19
                                # Let's use label_asym_id (column 6) as it's more reliable
                                chain = parts[6]
                                
                                # Get residue number (auth_seq_id - near end)
                                # Let's use label_seq_id (column 8) 
                                res_num = parts[8]
                                
                                # Check if it's RNA/DNA
                                if res_name in rna_residues:
                                    nt_id = f"{chain}-{res_name}-{res_num}-"
                                    nucleotides.add(nt_id)
                                    
                            except (IndexError, ValueError):
                                continue
            
            return len(nucleotides)
            
        except Exception as e:
            print(f"Warning: Error reading CIF file: {e}")
            return 0
    
    
    def get_nucleotide_count(self, pdb_id: str) -> int:
        """
        Get nucleotide count by downloading CIF from RCSB PDB.
        
        Downloads CIF, counts nucleotides, deletes CIF.
        
        Args:
            pdb_id: PDB ID
            
        Returns:
            Number of unique RNA nucleotides (0 if download fails)
        """
        temp_cif = "temp_structure.cif"
        
        # Download CIF
        if not self.download_cif(pdb_id, temp_cif):
            return 0
        
        # Count nucleotides
        count = self.count_nucleotides_from_cif(temp_cif)
        
        # Clean up temp file
        try:
            Path(temp_cif).unlink()
        except:
            pass
        
        return count
    
    
    def get_validation_metrics(self, pdb_id: str) -> Optional[dict]:
        """
        Extract validation metrics and structure metadata from RCSB API for a given PDB ID.
        
        Args:
            pdb_id: PDB ID
            
        Returns:
            Dictionary containing validation metrics and metadata, or None if failed
        """
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code != 200:
                print(f"Warning: Could not fetch data for {pdb_id} (HTTP {response.status_code})")
                return None
            
            data = response.json()
            metrics = {}
            
            # Primary validation metrics from pdbx_vrpt_summary_geometry
            if 'pdbx_vrpt_summary_geometry' in data and len(data['pdbx_vrpt_summary_geometry']) > 0:
                geom_data = data['pdbx_vrpt_summary_geometry'][0]
                metrics.update({
                    'clashscore': geom_data.get('clashscore'),
                    'angles_rmsz': geom_data.get('angles_rmsz'),
                    'bonds_rmsz': geom_data.get('bonds_rmsz'),
                    'percent_ramachandran_outliers': geom_data.get('percent_ramachandran_outliers'),
                    'percent_rotamer_outliers': geom_data.get('percent_rotamer_outliers')
                })
                
            # RNA backbone quality from pdbx_vrpt_summary
            if 'pdbx_vrpt_summary' in data:
                vrpt_summary = data['pdbx_vrpt_summary']
                metrics['rnasuiteness'] = vrpt_summary.get('rnasuiteness', None)
            
            # Experimental method and resolution from rcsb_entry_info
            if 'rcsb_entry_info' in data:
                entry_info = data['rcsb_entry_info']
                metrics['Experimental_method'] = entry_info.get('experimental_method', None)
                
                # Resolution - try multiple sources
                resolution = None
                if 'resolution_combined' in entry_info and entry_info['resolution_combined']:
                    resolution = entry_info['resolution_combined'][0]  # Get first/best resolution
                elif 'diffrn_resolution_high' in entry_info and entry_info['diffrn_resolution_high']:
                    resolution = entry_info['diffrn_resolution_high'].get('value', None)
                
                # Set resolution based on experimental method
                if metrics.get('Experimental_method') == 'EM':
                    metrics['EM Resolution (Å)'] = resolution
                    metrics['EM Diffraction Resolution (Å)'] = None
                else:
                    metrics['EM Resolution (Å)'] = None
                    metrics['EM Diffraction Resolution (Å)'] = resolution
            
            # Deposition date from rcsb_accession_info
            if 'rcsb_accession_info' in data:
                accession_info = data['rcsb_accession_info']
                deposit_date = accession_info.get('deposit_date', None)
                if deposit_date:
                    # Extract just the date part (YYYY-MM-DD) from ISO format
                    metrics['Deposition_Date'] = deposit_date.split('T')[0] if 'T' in deposit_date else deposit_date
                else:
                    metrics['Deposition_Date'] = None
            
            # Refinement statistics from refine[0]
            if 'refine' in data and len(data['refine']) > 0:
                refine_data = data['refine'][0]
                metrics.update({
                    'R_free': refine_data.get('ls_rfactor_rfree', None),
                    'R_work': refine_data.get('ls_rfactor_rwork', None),
                    'Refinement_resolution': refine_data.get('ls_dres_high', None),
                    'Average_B_factor': refine_data.get('biso_mean', None)
                })
            
            # Structure determination method from exptl[0]
            if 'exptl' in data and len(data['exptl']) > 0:
                exptl_data = data['exptl'][0]
                metrics['Structure Determination Method'] = exptl_data.get('method', None)
            
            # Data collection quality from pdbx_vrpt_summary_diffraction[0]
            if 'pdbx_vrpt_summary_diffraction' in data and len(data['pdbx_vrpt_summary_diffraction']) > 0:
                diff_data = data['pdbx_vrpt_summary_diffraction'][0]
                metrics.update({
                    'data_completeness': diff_data.get('data_completeness', None),
                    'iover_sigma': diff_data.get('iover_sigma', None),
                    'data_anisotropy': diff_data.get('data_anisotropy', None),
                    'fo_fc_correlation': diff_data.get('fo_fc_correlation', None)
                })
            
            # RNA structural features from rcsb_entry_info
            if 'rcsb_entry_info' in data:
                entry_info = data['rcsb_entry_info']
                if 'ndb_struct_conf_na_feature_combined' in entry_info:
                    features = entry_info['ndb_struct_conf_na_feature_combined']
                    metrics['RNA_structural_features'] = ', '.join(features) if features else None
            
            # Publication year from rcsb_primary_citation
            if 'rcsb_primary_citation' in data:
                citation = data['rcsb_primary_citation']
                metrics['Publication_Year'] = citation.get('year', None)
            
            # Software used
            if 'software' in data:
                software_list = []
                for sw in data['software']:
                    name = sw.get('name', '')
                    version = sw.get('version', '')
                    if version:
                        software_list.append(f"{name} {version}")
                    else:
                        software_list.append(name)
                metrics['Software_used'] = ', '.join(software_list) if software_list else None
            
            # Revision information from rcsb_accession_info
            if 'rcsb_accession_info' in data:
                accession_info = data['rcsb_accession_info']
                metrics.update({
                    'Revision_Date': accession_info.get('revision_date', '').split('T')[0] if accession_info.get('revision_date') else None,
                    'Major_Revision': accession_info.get('major_revision', None),
                    'Minor_Revision': accession_info.get('minor_revision', None)
                })
            
            print(f"✓ Extracted {len(metrics)} validation metrics and metadata for {pdb_id}")
            return metrics
        
        except Exception as e:
            print(f"Error fetching validation metrics: {e}")
            return None
