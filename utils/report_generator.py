"""Generate quality assessment reports with hotspot analysis."""

import json
from typing import Dict
import csv
import os
from pathlib import Path
from collections import Counter

class ReportGenerator:
    """Generates human-readable and JSON reports."""
    
    def __init__(self, config):
        self.config = config
    
    def print_report(self, results: Dict):
        """Print formatted console report - CRITICAL ISSUES ONLY."""
        if 'error' in results:
            print(f"\n{'='*60}")
            print(f"ERROR: {results['error']}")
            print(f"{'='*60}\n")
            return
        
        print(f"\n{'='*60}")
        print("RNA STRUCTURE QUALITY REPORT")
        print(f"{'='*60}")
        print(f"PDB ID: {results['pdb_id']}")
        print(f"Base Score: {results['base_score']}")
        print(f"Total Penalty: -{results['total_penalty']}")
        print(f"FINAL SCORE: {results['final_score']}/100")
        print(f"Quality Grade: {results['grade']}")
        
        # Base pair analysis - Summary only
        if 'base_pairs' in results['analyses']:
            bp = results['analyses']['base_pairs']
            print(f"\n{'='*60}")
            print(f"BASE PAIR ANALYSIS")
            print(f"{'='*60}")
            print(f"Total Pairs: {bp['total_pairs']}")
            print(f"  Canonical (cWW): {bp['stats']['canonical_wc']}")
            print(f"  Non-canonical: {bp['stats']['non_canonical']}")
            print(f"Penalty: -{bp['penalty']} points")
        
        # H-bond analysis - Summary only
        if 'hbonds' in results['analyses']:
            hb = results['analyses']['hbonds']
            print(f"\n{'='*60}")
            print(f"HYDROGEN BOND ANALYSIS")
            print(f"{'='*60}")
            print(f"Total H-bonds: {hb['total_hbonds']}")
            print(f"Penalty: -{hb['penalty']} points")
        
        # Hotspot analysis - CRITICAL and SEVERE in console
        if 'hotspots' in results and results['hotspots']:
            # Filter to CRITICAL and SEVERE issues for console display
            critical_severe_hotspots = [h for h in results['hotspots'] 
                                       if h['severity'] in ['CRITICAL', 'SEVERE']]
            
            print(f"\n{'='*60}")
            print(f"DAMAGED REGIONS SUMMARY (HOTSPOTS)")
            print(f"{'='*60}")
            
            summary = results.get('hotspot_summary', {})
            if summary:
                total_all = summary['total_hotspots']
                critical_count = summary.get('critical', 0)
                severe_count = summary.get('severe', 0)
                moderate_count = summary.get('moderate', 0)
                minor_count = summary.get('minor', 0)
                
                print(f"Total Issues Found: {total_all}")
                print(f"  Critical: {critical_count}")
                print(f"  Severe: {severe_count}")
                print(f"  Moderate: {moderate_count}")
                print(f"  Minor: {minor_count}")
            print(f"\n{'='*60}\n")
    
    def save_json_report(self, results: Dict, output_file: str):
        """Save results as JSON file - ALL ISSUES INCLUDED."""
        try:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            # Count issues by severity for confirmation message
            if 'hotspots' in results and results['hotspots']:
                severity_counts = {}
                for hotspot in results['hotspots']:
                    sev = hotspot['severity']
                    severity_counts[sev] = severity_counts.get(sev, 0) + 1
                
                print(f"Full report saved to: {output_file}")
                print(f"  All {len(results['hotspots'])} issues saved:", end=" ")
                print(f"{severity_counts.get('CRITICAL', 0)} critical, "
                      f"{severity_counts.get('SEVERE', 0)} severe, "
                      f"{severity_counts.get('MODERATE', 0)} moderate, "
                      f"{severity_counts.get('MINOR', 0)} minor\n")
            else:
                print(f"Report saved to: {output_file}\n")
        except Exception as e:
            print(f"Error saving report: {e}")
            
    def save_hotspot_summary(self, results: Dict, output_file: str):
        """Save a simplified hotspot-only summary - ALL ISSUES INCLUDED."""
        if 'hotspots' not in results or not results['hotspots']:
            print("No hotspots to save.")
            return
    
        # Count severity levels
        critical_count = sum(1 for h in results['hotspots'] if h['severity'] == 'CRITICAL')
        severe_count = sum(1 for h in results['hotspots'] if h['severity'] == 'SEVERE')
        moderate_count = sum(1 for h in results['hotspots'] if h['severity'] == 'MODERATE')
        minor_count = sum(1 for h in results['hotspots'] if h['severity'] == 'MINOR')
    
        # Calculate critical+severe count for percentage and recommendation
        critical_severe_count = critical_count + severe_count
    
        # Get total base pairs from results
        total_base_pairs = results['analyses']['base_pairs']['total_pairs']
    
        # Calculate percentage (only CRITICAL + SEVERE)
        hotspot_percentage = (critical_severe_count / total_base_pairs * 100) if total_base_pairs > 0 else 0
    
        # Generate recommendation based on score and critical/severe count
        recommendation = self._generate_recommendation(
            results['final_score'], 
            results['grade'],
            critical_severe_count,
            critical_count
        )
        
        # Convert hotspots to include issues_details
        hotspot_data = []

        for h in results['hotspots']:
            hotspot_dict = {
                'region': h['region'],
                'chain': h['chain'],
                'start_res': h['start_res'],
                'end_res': h['end_res'],
                #'residue_count': h['residue_count'],
                'score': h['score'],
                'severity': h['severity'],
                'issue_density': h['issue_density'],
                #'dominant_issues': h['dominant_issues'],
                #'base_pairs_affected': h['base_pairs_affected'],
                #'hbonds_affected': h['hbonds_affected'],
                'detailed_issues': h.get('detailed_issues', {}),
            }
            hotspot_data.append(hotspot_dict)
    
        summary = {
            'pdb_id': results['pdb_id'],
            'overall_score': results['final_score'],
            'overall_grade': results['grade'],
            'total_base_pairs': total_base_pairs,
            'hotspot_count': len(results['hotspots']),
            #'critical_severe_count': critical_severe_count,
            'hotspot_percentage': f"{hotspot_percentage:.2f}%",
            'recommendation': recommendation,
            'hotspot_breakdown': {
            'critical': critical_count,
            'severe': severe_count,
            #'moderate': moderate_count,
            #'minor': minor_count
        },
            # 'hotspots': results['hotspots']
            'hotspots': hotspot_data
        }
    
        try:
            with open(output_file, 'w') as f:
                json.dump(summary, f, indent=2)
            print(f"Hotspot summary (all severities) saved to: {output_file}")
        except Exception as e:
            print(f"Error saving hotspot summary: {e}")

    def _generate_recommendation(self, score: float, grade: str, critical_severe_count: int, critical_count: int) -> str:
        """Generate recommendation text based on score and hotspot counts."""
    
        # Critical issues present
        if critical_count > 0:
            if critical_count == 1:
                return f"Structure has {critical_severe_count} problematic region{'s' if critical_severe_count > 1 else ''} including 1 critical area requiring immediate attention"
            else:
                return f"Structure has {critical_severe_count} problematic regions including {critical_count} critical areas requiring immediate attention"
    
        # Severe issues only
        elif critical_severe_count > 0:
            if grade == "EXCELLENT QUALITY":
                return f"Excellent structure with {critical_severe_count} minor localized region{'s' if critical_severe_count > 1 else ''} that could benefit from refinement"
            elif grade == "GOOD QUALITY":
                return f"Generally well-refined structure with {critical_severe_count} localized region{'s' if critical_severe_count > 1 else ''} requiring attention"
            elif grade == "FAIR QUALITY":
                return f"Fair quality structure with {critical_severe_count} problematic region{'s' if critical_severe_count > 1 else ''} needing refinement"
            else:
                return f"Poorly refined structure with {critical_severe_count} severe region{'s' if critical_severe_count > 1 else ''} requiring significant work"
    
        # No critical or severe issues
        else:
            if grade == "EXCELLENT QUALITY":
                return "Excellent quality structure with no significant structural issues"
            elif grade == "GOOD QUALITY":
                return "Good quality structure with only minor issues distributed throughout"
            elif grade == "FAIR QUALITY":
                return "Fair quality structure with distributed minor issues"
            else:
                return "Structure quality could be improved with additional refinement"
            
            
    
    def save_issues_json(self, results: Dict, output_file: str = "issues.json"):
        """Save all base pair and hydrogen bond issues to a JSON file."""
        
        if 'error' in results:
            print(f"No issues to save: {results['error']}")
            return
        
        all_issues = []
        
        # Process base pair issues from detailed_issues if available
        if 'base_pairs' in results['analyses']:
            bp_analysis = results['analyses']['base_pairs']
            
            # Check if we have detailed issues stored
            if 'detailed_issues' in bp_analysis:
                all_issues.extend(bp_analysis['detailed_issues'])
            else:
                # Fallback to parsing issue strings
                self._parse_bp_issue_strings(bp_analysis.get('issues', []), all_issues)
        
        # Process hydrogen bond issues
        if 'hbonds' in results['analyses']:
            hb_analysis = results['analyses']['hbonds']
            
            if 'detailed_issues' in hb_analysis:
                all_issues.extend(hb_analysis['detailed_issues'])
            else:
                # Fallback to parsing issue strings
                self._parse_hb_issue_strings(hb_analysis.get('issues', []), all_issues)
        
        # Save to file
        try:
            with open(output_file, 'w') as f:
                json.dump({"issues": all_issues}, f, indent=2)
            print(f"All issues saved to: {output_file}")
            print(f"Total issues found: {len(all_issues)}")
            
            # Print summary by type
            bp_count = sum(1 for issue in all_issues if issue.get('type') == 'base_pair')
            hb_count = sum(1 for issue in all_issues if issue.get('type') == 'hydrogen_bond')
            print(f"  - Base pair issues: {bp_count}")
            print(f"  - Hydrogen bond issues: {hb_count}")
            
        except Exception as e:
            print(f"Error saving issues JSON: {e}")

    def _parse_bp_issue_strings(self, issue_strings: list, all_issues: list):
        """Parse base pair issue strings (fallback method)."""
        for issue in issue_strings:
            if ":" in issue:
                issue_type, rest = issue.split(":", 1)
                residues_part = rest.split("(")[0].strip()
                
                # Extract residue information
                residues = residues_part.split("-")
                if len(residues) >= 4:
                    chain_1 = residues[0]
                    residue_1 = int(residues[2])
                    chain_2 = residues[3] 
                    residue_2 = int(residues[5]) if len(residues) > 5 else 0
                    
                    bp_issue = {
                        "type": "base_pair",
                        "residues": f"{chain_1}-{residue_1} - {chain_2}-{residue_2}",
                        "chain_1": chain_1,
                        "chain_2": chain_2,
                        "residue_1": residue_1,
                        "residue_2": residue_2,
                        "specific_issues": [issue_type.lower().replace(" ", "_")],
                        "description": issue
                    }
                    all_issues.append(bp_issue)

    def _parse_hb_issue_strings(self, issue_strings: list, all_issues: list):
        """Parse hydrogen bond issue strings (fallback method)."""
        for issue in issue_strings:
            if ":" in issue:
                issue_type, rest = issue.split(":", 1)
                residues_part = rest.split("(")[0].strip()
                
                # Extract residue and atom information
                parts = residues_part.split()
                if len(parts) >= 3:
                    residues = parts[0]  # "res1-res2"
                    atoms = parts[1].strip("()")  # "atom1-atom2"
                    
                    # Parse residues
                    res_parts = residues.split("-")
                    if len(res_parts) >= 6:
                        chain_1 = res_parts[0]
                        residue_1 = int(res_parts[2])
                        chain_2 = res_parts[3]
                        residue_2 = int(res_parts[5])
                        
                        hb_issue = {
                            "type": "hydrogen_bond", 
                            "residues": f"{chain_1}-{residue_1} - {chain_2}-{residue_2}",
                            "chain_1": chain_1,
                            "chain_2": chain_2,
                            "residue_1": residue_1,
                            "residue_2": residue_2,
                            "atoms": atoms,
                            "specific_issues": [issue_type.lower().replace(" ", "_")],
                            "description": issue
                        }
                        all_issues.append(hb_issue)     
    
    def save_hotspot_csv(self, results: Dict, csv_file: str = "hotspot.csv"):
        """
        Save hotspot data to CSV file (replace if PDB exists, append if new).
        Preserves manual annotations (Looks_Bad, Notes) when updating.
        
        Args:
            results: Analysis results dictionary
            csv_file: Path to CSV file (default: hotspot.csv in root)
        """

        
        if 'hotspots' not in results or not results['hotspots']:
            return
        
        # Define CSV columns
        fieldnames = [
            'PDB_ID',
            'Hotspot_ID',
            'Region',
            'Chain',
            'Start_Res',
            'End_Res',
            'Span_Length',
            'Score',
            'Severity',
            'Issue_Density',
            'Num_Base_Pairs',
            'Num_Issues',
            'Issues_Per_BP',
            'Dominant_Issue_Types',
            'Has_Geometry_Issues',
            'Has_HBond_Issues',
            'Has_Count_Issues',
            'Avg_DSSR_Score',
            'Overall_Structure_Score',
            'Overall_Grade',
            'Looks_Bad',
            'Notes'
        ]
        
        pdb_id = results['pdb_id']
        
        # Read existing data if file exists
        existing_rows = []
        manual_annotations = {}
        file_exists = Path(csv_file).exists()
        
        if file_exists:
            try:
                with open(csv_file, 'r', newline='') as f:
                    reader = csv.DictReader(f)
                    
                    for row in reader:
                        if row['PDB_ID'] == pdb_id:
                            # Save manual annotations from current PDB
                            key = (row['PDB_ID'], row['Hotspot_ID'])
                            manual_annotations[key] = {
                                'Looks_Bad': row.get('Looks_Bad', ''),
                                'Notes': row.get('Notes', '')
                            }
                        else:
                            # Keep rows from other PDBs
                            existing_rows.append(row)
                
                print(f"\n✓ Updating existing CSV: {csv_file}")
                print(f"  Replacing data for {pdb_id} (preserving manual annotations)")
            except Exception as e:
                print(f"Warning: Could not read existing CSV: {e}")
                existing_rows = []
                manual_annotations = {}
        else:
            print(f"\n✓ Creating new CSV file: {csv_file}")
        
        # Prepare new rows for current PDB
        new_rows = []
        overall_score = results.get('final_score', 0)
        overall_grade = results.get('grade', 'UNKNOWN')
        
        # CRITICAL FIX: Filter out None hotspots and re-index
        #valid_hotspots = [h for h in valid_hotspots['hotspots'] if h is not None]
        
        for idx, hotspot in enumerate(results['hotspots'], start=1):
            # Calculate metrics
            span_length = hotspot['end_res'] - hotspot['start_res'] + 1
            
            # ===== CRITICAL FIX: Get counts from details =====
            details = hotspot.get('details', {})
            
            print(f"\n=== DEBUG Hotspot {idx} ===")
            print(f"Details keys: {details.keys()}")
            print(f"all_base_pairs exists: {'all_base_pairs' in details}")
            print(f"problematic_base_pairs_count: {details.get('problematic_base_pairs_count', 'MISSING')}")
            
            if 'all_base_pairs' in details:
                print(f"all_base_pairs length: {len(details['all_base_pairs'])}")
            else:
                print("all_base_pairs is MISSING!")
            
            # Get ALL base pairs (not just problematic ones)
            all_base_pairs = details.get('all_base_pairs', [])
            num_base_pairs = len(all_base_pairs)
            
            # Get PROBLEMATIC base pairs count
            num_issues = details.get('problematic_base_pairs_count', 0)
            
            # Calculate issues per base pair
            issues_per_bp = round(num_issues / num_base_pairs, 2) if num_base_pairs > 0 else 0
            
            # Extract dominant issue types
            dominant_issues = self._extract_dominant_issue_types(hotspot)
            
            # Check issue type presence
            has_geometry, has_hbond, has_count = self._check_issue_categories(hotspot)
            
            # Calculate average DSSR score
            avg_dssr_score = self._calculate_avg_dssr_score(hotspot)
            
            # Restore manual annotations if they exist
            key = (pdb_id, str(idx))
            annotations = manual_annotations.get(key, {'Looks_Bad': '', 'Notes': ''})
            
            
            row = {
                'PDB_ID': pdb_id,
                'Hotspot_ID': idx,
                'Region': hotspot['region'],
                'Chain': hotspot['chain'],
                'Start_Res': hotspot['start_res'],
                'End_Res': hotspot['end_res'],
                'Span_Length': span_length,
                'Score': hotspot['score'],
                'Severity': hotspot['severity'],
                'Issue_Density': hotspot['issue_density'],
                'Num_Base_Pairs': num_base_pairs,  # ← NOW using all_base_pairs
                'Num_Issues': num_issues,
                'Issues_Per_BP': issues_per_bp,
                'Dominant_Issue_Types': dominant_issues,
                'Has_Geometry_Issues': has_geometry,
                'Has_HBond_Issues': has_hbond,
                'Has_Count_Issues': has_count,
                'Avg_DSSR_Score': avg_dssr_score,
                'Overall_Structure_Score': overall_score,
                'Overall_Grade': overall_grade,
                'Looks_Bad': annotations['Looks_Bad'],
                'Notes': annotations['Notes']
            }
            
            new_rows.append(row)
        
        # Combine: existing rows (without current PDB) + new rows
        all_rows = existing_rows + new_rows
        
        # Write all data back to file
        try:
            with open(csv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_rows)
            
            print(f"✓ Added {len(new_rows)} hotspot(s) for {pdb_id} to CSV")
            if manual_annotations:
                print(f"  Preserved {len(manual_annotations)} manual annotation(s)")
            print(f"  Total hotspots in CSV: {len(all_rows)}")
            
        except Exception as e:
            print(f"Error saving hotspot CSV: {e}")   
    
    def _extract_dominant_issue_types(self, hotspot: Dict) -> str:
        """
        Extract the most common issue types from a hotspot.
        
        Returns a comma-separated string of dominant issues.
        """
        from collections import Counter
        
        # Count all issue types across all detailed_issues
        issue_counter = Counter()
        
        for detail in hotspot.get('detailed_issues', []):
            for issue in detail.get('specific_issues', []):
                issue_counter[issue] += 1
        
        # Get top 3 most common issues
        if issue_counter:
            top_issues = [issue for issue, count in issue_counter.most_common(3)]
            return '; '.join(top_issues)
        else:
            return 'unknown'


    def _check_issue_categories(self, hotspot: Dict) -> tuple:
        """
        Check which issue categories are present in the hotspot.
        
        Returns:
            (has_geometry, has_hbond, has_count): Boolean tuple
        """
        has_geometry = False
        has_hbond = False
        has_count = False
        
        for detail in hotspot.get('detailed_issues', []):
            issue_types = detail.get('issue_types', [])
            
            if 'base_pair_geometry' in issue_types:
                has_geometry = True
            if 'hbond_geometry' in issue_types:
                has_hbond = True
            if 'hbond_count' in issue_types:
                has_count = True
        
        return has_geometry, has_hbond, has_count


    def _calculate_avg_dssr_score(self, hotspot: Dict) -> float:
        """
        Calculate average DSSR H-bond score for all base pairs in hotspot.
        
        Returns:
            Average DSSR score (0.0-1.0), or -1 if no scores available
        """
        dssr_scores = []
        
        for detail in hotspot.get('detailed_issues', []):
            score = detail.get('dssr_hbond_score')
            if score is not None:
                dssr_scores.append(score)
        
        if dssr_scores:
            return round(sum(dssr_scores) / len(dssr_scores), 2)
        else:
            return -1.0  # Indicates no DSSR scores available
        
        
        
        
        
    def save_issues_csv(self, issues_file: str = "issues.json", csv_file: str = "issues_analysis.csv"):
        """
        Generate issues_analysis.csv from issues.json file.
        - Appends new PDB data to existing CSV
        - Replaces data if PDB_ID already exists
        
        Args:
            issues_file: Path to issues.json file
            csv_file: Path to output CSV file
        """
        # Read issues.json
        try:
            with open(issues_file, 'r') as f:
                data = json.load(f)
        except FileNotFoundError:
            print(f"Warning: {issues_file} not found. Skipping issues CSV generation.")
            return
        except Exception as e:
            print(f"Error reading {issues_file}: {e}")
            return
        
        issues = data.get('issues', [])
        
        if not issues:
            print("No issues found in issues.json. Skipping CSV generation.")
            return
        
        # Get PDB_ID from report.json
        pdb_id = "UNKNOWN"
        try:
            with open('report.json', 'r') as f:
                report = json.load(f)
                pdb_id = report.get('pdb_id', 'UNKNOWN')
                overall_structure_score = report.get('final_score', 0.0)
                overall_structure_quality = report.get('grade', 'UNKNOWN')                
        except:
            pass
        
        # Define CSV columns
        fieldnames = [
            'PDB_ID',
            'Chain',
            'Residue_1',
            'Residue_2',
            'Base_Pair',
            'Issues',
            'overall_structure_score',
            'overall_structure_quality',
            'Comments'
        ]
        
        # Read existing data if file exists
        existing_rows = []
        file_exists = Path(csv_file).exists()
        
        if file_exists:
            try:
                with open(csv_file, 'r', newline='') as f:
                    reader = csv.DictReader(f)
                    
                    for row in reader:
                        # Keep rows from other PDBs (exclude current PDB)
                        if row['PDB_ID'] != pdb_id:
                            existing_rows.append(row)
                
                print(f"\n✓ Updating existing CSV: {csv_file}")
                print(f"  Replacing data for {pdb_id}")
            except Exception as e:
                print(f"Warning: Could not read existing CSV: {e}")
                existing_rows = []
        else:
            print(f"\n✓ Creating new CSV file: {csv_file}")
        
        # Prepare new rows for current PDB
        new_rows = []
        
        for issue in issues:
            # Only process base_pair type issues (skip hydrogen_bond for now)
            if issue.get('type') != 'base_pair':
                continue
            
            # Extract data
            chain = issue.get('chain_1', '')  # Use chain_1 as the primary chain
            residue_1 = issue.get('residue_1', '')
            residue_2 = issue.get('residue_2', '')
            base_pair = issue.get('bp_type', '')
            
            # Join specific issues into comma-separated string
            specific_issues = issue.get('specific_issues', [])
            issues_str = ', '.join(specific_issues) if specific_issues else ''
            
            # Comments can be left empty or filled with additional info
            comments = ''
            
            row = {
                'PDB_ID': pdb_id,
                'Chain': chain,
                'Residue_1': residue_1,
                'Residue_2': residue_2,
                'Base_Pair': base_pair,
                'Issues': issues_str,
                'overall_structure_score': overall_structure_score,
                'overall_structure_quality': overall_structure_quality,
                'Comments': comments
            }
            
            new_rows.append(row)
        
        # Combine: existing rows (without current PDB) + new rows
        all_rows = existing_rows + new_rows
        
        # Write all data back to file
        try:
            with open(csv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_rows)
            
            print(f"✓ Added {len(new_rows)} issue(s) for {pdb_id} to CSV")
            print(f"  Total issues in CSV: {len(all_rows)}")
            
        except Exception as e:
            print(f"Error saving issues CSV: {e}")
            

    def save_score_summary_csv(self, result_dict: Dict, csv_file: str = "scores_summary.csv", hbond_data=None, validation_metrics=None):
        """
        Save baseline scoring results to CSV file.
        Appends new PDB data or replaces if PDB_ID already exists.
        
        Args:
            result_dict: Result dictionary from Scorer
            csv_file: Path to output CSV file
        """
        from collections import Counter
        
        pdb_id = result_dict.get('pdb_id', 'UNKNOWN')
        
        # Check if this is a motif analysis (has basepair_scores with detailed info)
        is_motif = result_dict.get('analysis_type') == 'motif'
        has_detailed_scores = 'basepair_scores' in result_dict and result_dict['basepair_scores']
        
        fieldnames = [
            'PDB_ID',
            'Overall_Score',
            #'Grade',
            'Num_Nucleotides',
            'clashscore',
            'rnasuiteness',
            'angles_rmsz',
            'bonds_rmsz',
            'percent_ramachandran_outliers',
            'percent_rotamer_outliers',
            'Total_Base_Pairs',
            'Avg_BasePair_Score',
            'Num_Problematic_BPs',
            'Problematic_Percentage',
            
            # Geometry issues
            'Geom_Misaligned',
            'Geom_Twisted',
            'Geom_NonCoplanar',
            'Geom_PoorHBond',
            'Geom_ZeroHBond',
            
            # H-bond issues
            'HBond_BadDistance',
            'HBond_BadAngles',
            'HBond_BadDihedral',
            'HBond_WeakQuality',
            'HBond_IncorrectCount',
            
            'Dominant_Issues',
            
            # Summary
            'Summary',
            
            # API metadata (matching user's CSV column names)
            'EM Resolution (Å)',
            'EM Diffraction Resolution (Å)',
            'Experimental_method',
            'Deposition_Date',
            'Refinement_resolution',
            'Average_B_factor',
            'R_free',
            'R_work',
            'Structure Determination Method'
        ]
        
        # Add detailed issues column for motifs only
        if is_motif and has_detailed_scores:
            fieldnames.append('Detailed_Issues')
        
        # Note: Protein_Binding_Explanations only for motifs CSV, not full RNA CSV
        # Full RNA has bindings in JSON base pair objects instead
        
        # Read existing data if file exists
        existing_rows = []
        file_exists = Path(csv_file).exists()
        
        if file_exists:
            try:
                with open(csv_file, 'r', newline='') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        if row['PDB_ID'] != pdb_id:
                            # Ensure existing rows have new columns if we're adding them
                            if 'Detailed_Issues' not in row and is_motif and has_detailed_scores:
                                row['Detailed_Issues'] = 'N/A'
                            # Add API metadata columns if missing
                            api_columns = [
                                'EM Resolution (Å)', 'EM Diffraction Resolution (Å)', 'Experimental_method',
                                'Deposition_Date', 'Refinement_resolution', 'Average_B_factor',
                                'R_free', 'R_work', 'Structure Determination Method'
                            ]
                            for col in api_columns:
                                if col not in row:
                                    row[col] = 'N/A'
                            existing_rows.append(row)
                
                print(f"\n✓ Updating existing CSV: {csv_file}")
                print(f"  Replacing data for {pdb_id}")
            except Exception as e:
                print(f"Warning: Could not read existing CSV: {e}")
                existing_rows = []
        else:
            print(f"\n✓ Creating new CSV file: {csv_file}")
        
        # Extract data from result
        total_bps = result_dict.get('total_base_pairs', 0)
        
        clashscore = validation_metrics.get('clashscore', 'N/A') if validation_metrics else 'N/A'
        rnasuiteness = validation_metrics.get('rnasuiteness', 'N/A') if validation_metrics else 'N/A'
        angles_rmsz = validation_metrics.get('angles_rmsz', 'N/A') if validation_metrics else 'N/A'
        bonds_rmsz = validation_metrics.get('bonds_rmsz', 'N/A') if validation_metrics else 'N/A'
        percent_ramachandran_outliers = validation_metrics.get('percent_ramachandran_outliers', 'N/A') if validation_metrics else 'N/A'
        percent_rotamer_outliers = validation_metrics.get('percent_rotamer_outliers', 'N/A') if validation_metrics else 'N/A'   
        
        # Calculate number of unique nucleotides from base pairs
        #num_nucleotides = self._calculate_num_nucleotides(result_dict, hbond_data)
        # Get nucleotide count from result_dict (pre-calculated from PDB CIF)
        num_nucleotides = result_dict.get('num_nucleotides', 0)
        
        # Fallback to old method if CIF download failed
        if num_nucleotides == 0:
            num_nucleotides = self._calculate_num_nucleotides(result_dict, hbond_data)


        
        # Recompute issue counts strictly from below-baseline base pairs when detailed scores are available
        basepair_scores = result_dict.get('basepair_scores', [])
        geom_counts = {
            'misaligned': 0,
            'twisted': 0,
            'non_coplanar': 0,
            'poor_hbond': 0,   # poor DSSR H-bond quality stored under hbond_issues
            'zero_hbond': 0
        }
        hbond_counts = {
            'bad_distance': 0,
            'bad_angles': 0,
            'bad_dihedral': 0,
            'weak_quality': 0,
            'incorrect_count': 0
        }
        num_problematic = result_dict.get('num_problematic_bps', 0)
        
        if basepair_scores:
            num_problematic = 0
            for bp_score in basepair_scores:
                if bp_score.get('score', 100) < self.config.BASELINE:
                    num_problematic += 1
                    
                    geom_bp = bp_score.get('geometry_issues', {})
                    hbond_bp = bp_score.get('hbond_issues', {})
                    
                    for issue in ['misaligned', 'twisted', 'non_coplanar', 'zero_hbond']:
                        if geom_bp.get(issue, False):
                            geom_counts[issue] += 1
                    
                    if hbond_bp.get('poor_hbond', False):
                        geom_counts['poor_hbond'] += 1
                    
                    for issue in ['bad_distance', 'bad_angles', 'bad_dihedral', 'weak_quality', 'incorrect_count']:
                        if hbond_bp.get(issue, False):
                            hbond_counts[issue] += 1
        else:
            # Fallback to pre-aggregated counts if detailed scores are unavailable
            geom_source = result_dict.get('geometry_issues', {})
            hbond_source = result_dict.get('hbond_issues', {})
            for issue in geom_counts:
                geom_counts[issue] = geom_source.get(issue, 0)
            for issue in hbond_counts:
                hbond_counts[issue] = hbond_source.get(issue, 0)
        
        problematic_pct = (num_problematic / total_bps * 100) if total_bps > 0 else 0
        
        # Collect issue counts (below-baseline only when detailed scores exist)
        all_issue_counts = {}
        for issue, count in geom_counts.items():
            if count > 0:
                all_issue_counts[f"geom_{issue}"] = count
        for issue, count in hbond_counts.items():
            if count > 0:
                all_issue_counts[f"hbond_{issue}"] = count
        
        # Sort by count (descending) and create semicolon-separated list
        sorted_issues = sorted(all_issue_counts.items(), key=lambda x: x[1], reverse=True)
        dominant_issues_str = "; ".join([f"{issue}({count})" for issue, count in sorted_issues])
        if not dominant_issues_str:
            dominant_issues_str = "none"
        
        # CHANGED: Clean up summary text - remove formatting artifacts
        summary_text = result_dict.get('summary', '')
        # Remove pipes and newlines
        summary_text = summary_text.replace('|', '').replace('\n', ' ')
        # Remove leading dashes and extra spaces
        import re
        summary_text = re.sub(r'\s*-\s*', ' ', summary_text)  # Remove all dashes
        summary_text = re.sub(r'\s+', ' ', summary_text)  # Collapse multiple spaces
        summary_text = summary_text.strip()
        
        # Create row
        row = {
            'PDB_ID': pdb_id,
            'Overall_Score': result_dict.get('overall_score', 0),
            #'Grade': result_dict.get('grade', 'N/A'),
            'Num_Nucleotides': num_nucleotides,
            'clashscore': clashscore,
            'rnasuiteness': rnasuiteness,
            'angles_rmsz': angles_rmsz,
            'bonds_rmsz': bonds_rmsz,
            'percent_ramachandran_outliers': percent_ramachandran_outliers,
            'percent_rotamer_outliers': percent_rotamer_outliers,
            'Total_Base_Pairs': total_bps,
            'Avg_BasePair_Score': result_dict.get('avg_basepair_score', 0),
            'Num_Problematic_BPs': num_problematic,
            'Problematic_Percentage': round(problematic_pct, 1),
            
            # Geometry issues
            'Geom_Misaligned': geom_counts.get('misaligned', 0),
            'Geom_Twisted': geom_counts.get('twisted', 0),
            'Geom_NonCoplanar': geom_counts.get('non_coplanar', 0),
            'Geom_PoorHBond': geom_counts.get('poor_hbond', 0),
            'Geom_ZeroHBond': geom_counts.get('zero_hbond', 0),
            
            # H-bond issues
            'HBond_BadDistance': hbond_counts.get('bad_distance', 0),
            'HBond_BadAngles': hbond_counts.get('bad_angles', 0),
            'HBond_BadDihedral': hbond_counts.get('bad_dihedral', 0),
            'HBond_WeakQuality': hbond_counts.get('weak_quality', 0),
            'HBond_IncorrectCount': hbond_counts.get('incorrect_count', 0),
            
            # All issues ordered by frequency
            'Dominant_Issues': dominant_issues_str,
            
            # Cleaned summary
            'Summary': summary_text,
            
            # API metadata (using exact column names from user's CSV)
            'EM Resolution (Å)': validation_metrics.get('EM Resolution (Å)', 'N/A') if validation_metrics else 'N/A',
            'EM Diffraction Resolution (Å)': validation_metrics.get('EM Diffraction Resolution (Å)', 'N/A') if validation_metrics else 'N/A',
            'Experimental_method': validation_metrics.get('Experimental_method', 'N/A') if validation_metrics else 'N/A',
            'Deposition_Date': validation_metrics.get('Deposition_Date', 'N/A') if validation_metrics else 'N/A',
            'Refinement_resolution': validation_metrics.get('Refinement_resolution', 'N/A') if validation_metrics else 'N/A',
            'Average_B_factor': validation_metrics.get('Average_B_factor', 'N/A') if validation_metrics else 'N/A',
            'R_free': validation_metrics.get('R_free', 'N/A') if validation_metrics else 'N/A',
            'R_work': validation_metrics.get('R_work', 'N/A') if validation_metrics else 'N/A',
            'Structure Determination Method': validation_metrics.get('Structure Determination Method', 'N/A') if validation_metrics else 'N/A'
        }
        
        # Add detailed issues for motifs (for data analysis)
        # Only include base pairs with score below baseline threshold
        if is_motif and has_detailed_scores:
            problematic_pairs = []
            basepair_scores = result_dict.get('basepair_scores', [])
            
            # Collect base pairs with score below baseline
            for bp_score in basepair_scores:
                bp_info = bp_score.get('bp_info', {})
                res_1 = bp_info.get('res_1', '?')
                res_2 = bp_info.get('res_2', '?')
                bp_type = bp_info.get('bp_type', '?')
                score = bp_score.get('score', 100)
                
                # Only include pairs with score below baseline threshold
                if score < self.config.BASELINE:
                    # Collect all issues for this base pair
                    issues = []
                    geom_issues_bp = bp_score.get('geometry_issues', {})
                    hbond_issues_bp = bp_score.get('hbond_issues', {})
                    
                    for issue, is_present in geom_issues_bp.items():
                        if is_present:
                            issues.append(f"geom_{issue}")
                    for issue, is_present in hbond_issues_bp.items():
                        if is_present:
                            issues.append(f"hbond_{issue}")
                    
                    if issues:
                        bp_id = f"{res_1}-{res_2}"
                        issues_str = ",".join(issues)
                        num_issues = len(issues)
                        # Store: (score, num_issues, formatted_string) for sorting
                        problematic_pairs.append((score, num_issues, f"{bp_id}({bp_type}):{issues_str}"))
            
            # Sort by worst first: lowest score, then most issues
            problematic_pairs.sort(key=lambda x: (x[0], -x[1]))  # Sort by score (ascending), then by num_issues (descending)
            
            # Extract just the formatted strings
            detailed_issues_list = [pair[2] for pair in problematic_pairs]
            
            # Join with pipe separator for CSV (easy to parse)
            row['Detailed_Issues'] = "|".join(detailed_issues_list) if detailed_issues_list else "none"
        elif is_motif:
            # Motif but no detailed scores available
            row['Detailed_Issues'] = "N/A"
        
        # Note: Protein_Binding_Explanations removed from full RNA CSV
        # Bindings are now in JSON base pair objects (bp_score['protein_bindings']) for easier access
        
        # Combine rows
        all_rows = existing_rows + [row]
        
        # Debug: Check row has required fields
        missing_fields = [f for f in fieldnames if f not in row]
        if missing_fields:
            print(f"Warning: Row missing fields: {missing_fields[:5]}...")
        
        # Safety check: Remove any duplicates (keep latest)
        seen_pdbs = set()
        unique_rows = []
        for r in reversed(all_rows):
            if r['PDB_ID'] not in seen_pdbs:
                seen_pdbs.add(r['PDB_ID'])
                unique_rows.append(r)
        all_rows = list(reversed(unique_rows))
        
        # Write to file
        try:
            if len(all_rows) == 0:
                print(f"Warning: No rows to write to CSV!")
                return
            
            with open(csv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                
                # Ensure all rows have all required fields
                for r in all_rows:
                    # Add missing fields with default values
                    for field in fieldnames:
                        if field not in r:
                            r[field] = 'N/A'
                    writer.writerow(r)
            
            print(f"✓ Baseline summary saved to: {csv_file}")
            print(f"  Total structures in CSV: {len(all_rows)}")
            
        except Exception as e:
            import traceback
            print(f"Error saving baseline summary CSV: {e}")
            print(f"Traceback: {traceback.format_exc()}")

    def _calculate_num_nucleotides(self, result_dict, hbond_data=None):
        """
        Calculate number of unique RNA nucleotides.
        
        This counts ALL RNA nucleotides, whether they're in base pairs or not.
        Priority: hbond_data > basepair_scores > 0
        
        Args:
            result_dict: Result dictionary from scorer
            hbond_data: Optional DataFrame of H-bonds (recommended for accurate count)
            
        Returns:
            Number of unique RNA nucleotides
        """
        # Method 1: Count from H-bond data (most accurate)
        if hbond_data is not None and not hbond_data.empty:
            unique_nucleotides = set()
            
            # Get all RNA nucleotides from res_1
            rna_mask_1 = hbond_data['res_type_1'] == 'RNA'
            if rna_mask_1.any():
                unique_nucleotides.update(hbond_data[rna_mask_1]['res_1'].unique())
            
            # Get all RNA nucleotides from res_2
            rna_mask_2 = hbond_data['res_type_2'] == 'RNA'
            if rna_mask_2.any():
                unique_nucleotides.update(hbond_data[rna_mask_2]['res_2'].unique())
            
            if unique_nucleotides:
                return len(unique_nucleotides)
        
        # Method 2: Count from base pairs (fallback)
        basepair_scores = result_dict.get('basepair_scores', [])
        if basepair_scores:
            unique_nucleotides = set()
            for bp in basepair_scores:
                bp_info = bp.get('bp_info', {})
                res_1 = bp_info.get('res_1')
                res_2 = bp_info.get('res_2')
                if res_1:
                    unique_nucleotides.add(res_1)
                if res_2:
                    unique_nucleotides.add(res_2)
            
            if unique_nucleotides:
                return len(unique_nucleotides)
        
        # Method 3: No data available
        return 0

    def save_motifs_summary_csv(self, motif_data: Dict, motif_name: str = None, csv_file: str = "scores_motifs_summary.csv", csv_dir: str = None):
        """
        Save motif scoring results to CSV file.
        If csv_dir is provided, writes to individual CSV file per motif (no race conditions).
        Otherwise, appends to single CSV file.
        
        Args:
            motif_data: Dictionary containing motif scoring data
            motif_name: Name of the motif (e.g., from --motif-name or generated)
            csv_file: Path to output CSV file (used if csv_dir is None)
            csv_dir: Directory to write individual CSV files (one per motif)
        """
        fieldnames = [
            'Motif_Name', 'PDB_ID', 'Chain', 'Residue_Range',
            'Motif_Score', 'Full_Structure_Score', 'Score_Difference',
            'Motif_Length', 'Full_Structure_Length',
            'Total_Base_Pairs', 'Num_Problematic_BPs',
            'Num_Paired_Nucleotides', 'Pairing_Percentage',
            'Geom_Misaligned', 'Geom_Twisted', 'Geom_NonCoplanar',
            'Geom_PoorHBond', 'Geom_ZeroHBond',
            'HBond_BadDistance', 'HBond_BadAngles', 'HBond_BadDihedral',
            'HBond_WeakQuality', 'HBond_IncorrectCount',
            'Dominant_Issues',
            'Detailed_Issues',
            'Protein_Binding_Explanations'
        ]
        
        # Extract and format data from motif_data
        pdb_id = motif_data.get('pdb_id', 'UNKNOWN')
        motif_range = motif_data.get('motif_range', 'N/A')
        chain = motif_data.get('motif_chain', 'all')
        
        # Get motif name (use provided or generate from data)
        if not motif_name:
            motif_name = motif_data.get('motif_name', f"{pdb_id}_{motif_range}")
        
        # Calculate pairing percentage
        num_paired = motif_data.get('num_paired_nucleotides', 0)
        motif_length = motif_data.get('motif_num_nucleotides', 0)
        pairing_pct = (num_paired / motif_length * 100) if motif_length > 0 else 0
        
        # Recompute issue counts strictly from below-baseline base pairs
        basepair_scores = motif_data.get('basepair_scores', [])
        geom_counts = {
            'misaligned': 0,
            'twisted': 0,
            'non_coplanar': 0,
            'poor_hbond': 0,   # stored as geometry_issues in scorer for DSSR quality == poor
            'zero_hbond': 0
        }
        hbond_counts = {
            'bad_distance': 0,
            'bad_angles': 0,
            'bad_dihedral': 0,
            'weak_quality': 0,
            'incorrect_count': 0
        }
        num_problematic_bps = motif_data.get('num_problematic_bps', 0)
        
        detailed_issues = "none"
        if basepair_scores:
            problematic_pairs = []
            num_problematic_bps = 0
            
            for bp_score in basepair_scores:
                if bp_score.get('score', 100) < self.config.BASELINE:
                    num_problematic_bps += 1
                    
                    bp_info = bp_score.get('bp_info', {})
                    res_1 = bp_info.get('res_1', '?')
                    res_2 = bp_info.get('res_2', '?')
                    bp_type = bp_info.get('bp_type', '?')
                    score = bp_score.get('score', 100)
                    
                    issues = []
                    geom_issues_bp = bp_score.get('geometry_issues', {})
                    hbond_issues_bp = bp_score.get('hbond_issues', {})
                    
                    # Track geometry issues
                    for issue in ['misaligned', 'twisted', 'non_coplanar', 'zero_hbond']:
                        if geom_issues_bp.get(issue, False):
                            geom_counts[issue] += 1
                            issues.append(f"geom_{issue}")
                    
                    # DSSR poor_hbond is stored under hbond_issues, but we want it counted in geom_* column per legacy naming
                    if hbond_issues_bp.get('poor_hbond', False):
                        geom_counts['poor_hbond'] += 1
                        issues.append("hbond_poor_hbond")
                    
                    # Track H-bond issues
                    for issue in ['bad_distance', 'bad_angles', 'bad_dihedral', 'weak_quality', 'incorrect_count']:
                        if hbond_issues_bp.get(issue, False):
                            hbond_counts[issue] += 1
                            issues.append(f"hbond_{issue}")
                    
                    if issues:
                        bp_id = f"{res_1}-{res_2}"
                        issues_str = ",".join(issues)
                        num_issues = len(issues)
                        problematic_pairs.append((score, num_issues, f"{bp_id}({bp_type}):{issues_str}"))
            
            # Sort by worst first
            problematic_pairs.sort(key=lambda x: (x[0], -x[1]))
            detailed_issues_list = [pair[2] for pair in problematic_pairs]
            detailed_issues = "|".join(detailed_issues_list) if detailed_issues_list else "none"
        else:
            # Fallback to existing aggregated counts if no per-basepair detail is available
            geom_source = motif_data.get('geometry_issues', {})
            hbond_source = motif_data.get('hbond_issues', {})
            for issue in geom_counts:
                geom_counts[issue] = geom_source.get(issue, 0)
            for issue in hbond_counts:
                hbond_counts[issue] = hbond_source.get(issue, 0)
        
        # Build dominant issues string from the recomputed counts
        all_issue_counts = {}
        for issue, count in geom_counts.items():
            if count > 0:
                all_issue_counts[f"geom_{issue}"] = count
        for issue, count in hbond_counts.items():
            if count > 0:
                all_issue_counts[f"hbond_{issue}"] = count
        
        sorted_issues = sorted(all_issue_counts.items(), key=lambda x: x[1], reverse=True)
        dominant_issues_str = "; ".join([f"{issue}({count})" for issue, count in sorted_issues]) if sorted_issues else "none"
        
        # Get protein binding explanations
        protein_bindings = motif_data.get('protein_binding_explanations', {})
        if protein_bindings:
            binding_strings = []
            for bp_id, bindings in protein_bindings.items():
                bindings_str = ",".join(bindings)
                binding_strings.append(f"{bp_id}:{bindings_str}")
            protein_bindings_str = "|".join(binding_strings)
        else:
            protein_bindings_str = "none"
        
        # Create row dictionary
        row = {
            'Motif_Name': motif_name,
            'PDB_ID': pdb_id,
            'Chain': chain,
            'Residue_Range': motif_range,
            'Motif_Score': motif_data.get('motif_score', 0),
            'Full_Structure_Score': motif_data.get('full_structure_score', 0),
            'Score_Difference': motif_data.get('score_difference', 0),
            'Motif_Length': motif_length,
            'Full_Structure_Length': motif_data.get('full_structure_num_nucleotides', 0),
            'Total_Base_Pairs': motif_data.get('total_base_pairs', 0),
            'Num_Problematic_BPs': num_problematic_bps,
            'Num_Paired_Nucleotides': num_paired,
            'Pairing_Percentage': round(pairing_pct, 1),
            'Geom_Misaligned': geom_counts.get('misaligned', 0),
            'Geom_Twisted': geom_counts.get('twisted', 0),
            'Geom_NonCoplanar': geom_counts.get('non_coplanar', 0),
            'Geom_PoorHBond': geom_counts.get('poor_hbond', 0),
            'Geom_ZeroHBond': geom_counts.get('zero_hbond', 0),
            'HBond_BadDistance': hbond_counts.get('bad_distance', 0),
            'HBond_BadAngles': hbond_counts.get('bad_angles', 0),
            'HBond_BadDihedral': hbond_counts.get('bad_dihedral', 0),
            'HBond_WeakQuality': hbond_counts.get('weak_quality', 0),
            'HBond_IncorrectCount': hbond_counts.get('incorrect_count', 0),
            'Dominant_Issues': dominant_issues_str,
            'Detailed_Issues': detailed_issues,
            'Protein_Binding_Explanations': protein_bindings_str
        }
        
        # If csv_dir is provided, write to individual CSV file (no race conditions)
        if csv_dir:
            csv_dir_path = Path(csv_dir)
            csv_dir_path.mkdir(parents=True, exist_ok=True)
            individual_csv = csv_dir_path / f"{motif_name}.csv"
            
            try:
                # Write individual CSV file with header
                with open(individual_csv, 'w', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerow(row)
                
                print(f"✓ Motif summary saved to: {individual_csv}")
                return True
            except Exception as e:
                print(f"Error saving individual motif CSV: {e}")
                import traceback
                print(f"Traceback: {traceback.format_exc()}")
                return False
        
        # Otherwise, append to single CSV file (original behavior)
        file_exists = Path(csv_file).exists()
        
        try:
            with open(csv_file, 'a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                
                # Write header only if file is new
                if not file_exists:
                    writer.writeheader()
                
                # Append the new row
                writer.writerow(row)
            
            print(f"✓ Motif summary saved to: {csv_file}")
            return True
            
        except Exception as e:
            print(f"Error saving motifs summary CSV: {e}")
            import traceback
            print(f"Traceback: {traceback.format_exc()}")
            return False