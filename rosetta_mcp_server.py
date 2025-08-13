#!/usr/bin/env python3
"""
Rosetta MCP Server - Provides access to Rosetta/PyRosetta functions and properties
"""

import json
import sys
import os
from typing import Dict, Any, List
import subprocess
import importlib.util

# Try to import PyRosetta if available
try:
    import pyrosetta
    PYROSETTA_AVAILABLE = True
except ImportError:
    PYROSETTA_AVAILABLE = False

class RosettaMCPServer:
    def __init__(self):
        self.rosetta_path = self._find_rosetta_path()
        self.available_functions = self._get_rosetta_functions()
        
    def _find_rosetta_path(self) -> str:
        """Find Rosetta installation path"""
        possible_paths = [
            "/Users/arielben-sasson/dev/open_repos/rosetta",
            "/opt/rosetta",
            "/usr/local/rosetta",
            os.path.expanduser("~/rosetta")
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        return "Not found"
    
    def _get_rosetta_functions(self) -> Dict[str, Any]:
        """Get all available Rosetta functions and properties"""
        functions = {
            "rosetta_path": self.rosetta_path,
            "pyrosetta_available": PYROSETTA_AVAILABLE,
            "score_functions": self._get_score_functions(),
            "common_movers": self._get_common_movers(),
            "common_filters": self._get_common_filters(),
            "residue_selectors": self._get_residue_selectors(),
            "task_operations": self._get_task_operations(),
            "common_parameters": self._get_common_parameters(),
            "xml_protocols": self._get_xml_protocols(),
            "command_line_options": self._get_command_line_options()
        }
        return functions
    
    def _get_score_functions(self) -> List[str]:
        """Get available score functions"""
        if PYROSETTA_AVAILABLE:
            try:
                return pyrosetta.get_score_function_names()
            except:
                pass
        
        return [
            "ref2015", "ref15", "talaris2014", "talaris2013",
            "score12", "score13", "soft_rep", "hard_rep"
        ]
    
    def _get_common_movers(self) -> List[str]:
        """Get common Rosetta movers"""
        return [
            "FastDesign", "FastRelax", "PackRotamersMover",
            "MinMover", "ShearMover", "SmallMover",
            "AddCompositionConstraintMover", "ClearConstraintsMover",
            "SavePoseMover", "DumpPdb", "InterfaceAnalyzerMover"
        ]
    
    def _get_common_filters(self) -> List[str]:
        """Get common Rosetta filters"""
        return [
            "ScoreType", "Ddg", "Sasa", "ShapeComplementarity",
            "BuriedUnsatHbonds", "NetCharge", "Time",
            "InterfaceBindingEnergyDensityFilter", "ContactMolecularSurface"
        ]
    
    def _get_residue_selectors(self) -> List[str]:
        """Get common residue selectors"""
        return [
            "Chain", "Neighborhood", "And", "Not", "Or",
            "True", "False", "ResidueName", "ResidueIndex",
            "SecondaryStructure", "ResidueProperty"
        ]
    
    def _get_task_operations(self) -> List[str]:
        """Get common task operations"""
        return [
            "RestrictToRepacking", "ExtraRotamersGeneric",
            "InitializeFromCommandline", "IncludeCurrent",
            "LimitAromaChi2", "DesignRestrictions",
            "OperateOnResidueSubset", "PreventRepackingRLT"
        ]
    
    def _get_common_parameters(self) -> Dict[str, str]:
        """Get common Rosetta parameters"""
        return {
            "nstruct": "Number of structures to generate",
            "ex1": "Extra rotamers level 1",
            "ex2aro": "Extra rotamers level 2 for aromatics",
            "use_input_sc": "Use input side chains",
            "overwrite": "Overwrite existing output files",
            "out:suffix": "Output file suffix",
            "scorefile_format": "Score file format (json, silent, etc.)",
            "holes:dalphaball": "Path to DAlphaBall executable",
            "score:symmetric_gly_tables": "Use symmetric glycine tables",
            "nblist_autoupdate": "Auto-update neighbor list"
        }
    
    def _get_xml_protocols(self) -> List[str]:
        """Get available XML protocol files"""
        xml_dir = os.path.join(self.rosetta_path, "inputs", "xml")
        if os.path.exists(xml_dir):
            try:
                return [f for f in os.listdir(xml_dir) if f.endswith('.xml')]
            except:
                pass
        return ["2chain_interface_design_bias01_v2.xml", "interface_default.xml"]
    
    def _get_command_line_options(self) -> List[str]:
        """Get common command line options"""
        return [
            "-in:file:s", "-out:path:all", "-parser:protocol",
            "-holes:dalphaball", "-score:symmetric_gly_tables",
            "-out:suffix", "-nstruct", "-jd2:ntrials",
            "-out:pdb", "-overwrite", "-ex1", "-ex2aro",
            "-use_input_sc", "-nblist_autoupdate", "-chemical:exclude_patches"
        ]
    
    def get_rosetta_info(self) -> Dict[str, Any]:
        """Get comprehensive Rosetta information"""
        return self.available_functions
    
    def validate_xml(self, xml_content: str) -> Dict[str, Any]:
        """Validate XML protocol file"""
        try:
            # Basic XML validation
            import xml.etree.ElementTree as ET
            ET.fromstring(xml_content)
            return {"valid": True, "message": "XML is valid"}
        except Exception as e:
            return {"valid": False, "error": str(e)}
    
    def get_rosetta_help(self, topic: str = None) -> str:
        """Get help for specific Rosetta topics"""
        help_topics = {
            "score_functions": "Score functions define how Rosetta evaluates protein structures. Common ones include ref2015, ref15, and talaris2014.",
            "movers": "Movers are the core operations in Rosetta protocols. They modify protein structures, apply constraints, or analyze results.",
            "filters": "Filters evaluate structures and can halt protocols if conditions aren't met. Examples include DDG, SASA, and shape complementarity.",
            "xml": "XML files define Rosetta protocols by specifying movers, filters, and their order of execution.",
            "parameters": "Command line parameters control Rosetta behavior, from input/output to algorithm settings."
        }
        
        if topic and topic in help_topics:
            return help_topics[topic]
        return "Available topics: " + ", ".join(help_topics.keys())

def main():
    """Main function for MCP server"""
    server = RosettaMCPServer()
    
    # Simple command line interface
    if len(sys.argv) > 1:
        command = sys.argv[1]
        if command == "info":
            print(json.dumps(server.get_rosetta_info(), indent=2))
        elif command == "help":
            topic = sys.argv[2] if len(sys.argv) > 2 else None
            print(server.get_rosetta_help(topic))
        elif command == "validate":
            if len(sys.argv) > 2:
                xml_file = sys.argv[2]
                try:
                    with open(xml_file, 'r') as f:
                        xml_content = f.read()
                    result = server.validate_xml(xml_content)
                    print(json.dumps(result, indent=2))
                except Exception as e:
                    print(json.dumps({"error": str(e)}, indent=2))
            else:
                print("Usage: python rosetta_mcp_server.py validate <xml_file>")
        else:
            print("Available commands: info, help [topic], validate <xml_file>")
    else:
        print("Rosetta MCP Server")
        print("Available commands: info, help [topic], validate <xml_file>")

if __name__ == "__main__":
    main()
