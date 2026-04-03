#!/usr/bin/env node
/**
 * Rosetta MCP Server Wrapper
 * Makes the Python Rosetta server compatible with MCP protocol
 */

const { spawn } = require('child_process');
const path = require('path');
const readline = require('readline');
const fs = require('fs');
const os = require('os');
const crypto = require('crypto');

// Data-driven mapping table for Rosetta <-> Biotite translation
const BIOTITE_ROSETTA_MAPPINGS = [
    {
        id: 'structure_loading_pdb',
        category: 'Structure I/O',
        rosetta: {
            name: 'pose_from_pdb',
            module: 'pyrosetta',
            signature: 'pose_from_pdb(filename)',
            description: 'Load a PDB file into a Rosetta Pose object',
            example: [
                'import pyrosetta',
                'pyrosetta.init("-mute all")',
                'pose = pyrosetta.pose_from_pdb("input.pdb")'
            ]
        },
        biotite: {
            name: 'PDBFile.read / get_structure',
            module: 'biotite.structure.io.pdb',
            signature: 'PDBFile.read(filepath); get_structure(pdb_file)',
            description: 'Load a PDB file into a Biotite AtomArray',
            example: [
                'from biotite.structure.io.pdb import PDBFile',
                'pdb_file = PDBFile.read("input.pdb")',
                'structure = pdb_file.get_structure(model=1)'
            ]
        },
        notes: 'Biotite returns an AtomArray (numpy-based). Rosetta returns a Pose object with full energy-function support.',
        equivalence: 'full'
    },
    {
        id: 'structure_saving_pdb',
        category: 'Structure I/O',
        rosetta: {
            name: 'pose.dump_pdb',
            module: 'pyrosetta',
            signature: 'pose.dump_pdb(filename)',
            description: 'Save a Pose to a PDB file',
            example: ['pose.dump_pdb("output.pdb")']
        },
        biotite: {
            name: 'PDBFile / set_structure / write',
            module: 'biotite.structure.io.pdb',
            signature: 'PDBFile(); pdb_file.set_structure(atoms); pdb_file.write(filepath)',
            description: 'Write an AtomArray to a PDB file',
            example: [
                'from biotite.structure.io.pdb import PDBFile',
                'pdb_file = PDBFile()',
                'pdb_file.set_structure(structure)',
                'pdb_file.write("output.pdb")'
            ]
        },
        notes: 'Both produce standard PDB files. Biotite also supports mmCIF via CIFFile.',
        equivalence: 'full'
    },
    {
        id: 'structure_loading_cif',
        category: 'Structure I/O',
        rosetta: {
            name: 'pose_from_file (mmCIF)',
            module: 'pyrosetta',
            signature: 'pose_from_file(filename)',
            description: 'Load an mmCIF file into a Rosetta Pose',
            example: [
                'import pyrosetta',
                'pyrosetta.init("-mute all")',
                'pose = pyrosetta.pose_from_file("input.cif")'
            ]
        },
        biotite: {
            name: 'CIFFile.read / get_structure',
            module: 'biotite.structure.io.pdbx',
            signature: 'CIFFile.read(filepath); get_structure(cif_file)',
            description: 'Load an mmCIF file into a Biotite AtomArray',
            example: [
                'from biotite.structure.io.pdbx import CIFFile',
                'cif_file = CIFFile.read("input.cif")',
                'structure = cif_file.get_structure(model=1)'
            ]
        },
        notes: 'Biotite has strong mmCIF/PDBx support including BinaryCIF.',
        equivalence: 'full'
    },
    {
        id: 'sasa',
        category: 'Surface Analysis',
        rosetta: {
            name: 'SasaMetric / calc_per_res_sasa',
            module: 'pyrosetta.rosetta.core.scoring',
            signature: 'calc_per_res_sasa(pose, atom_sasa, residue_sasa, probe_radius)',
            description: 'Calculate solvent-accessible surface area per residue',
            example: [
                'from pyrosetta.rosetta.core.scoring import calc_per_res_sasa',
                'from pyrosetta.rosetta.core.id import AtomID_Map_double_t',
                'from pyrosetta.rosetta.utility import vector1_double',
                'atom_sasa = AtomID_Map_double_t()',
                'residue_sasa = vector1_double()',
                'calc_per_res_sasa(pose, atom_sasa, residue_sasa, 1.4)'
            ]
        },
        biotite: {
            name: 'sasa',
            module: 'biotite.structure',
            signature: 'biotite.structure.sasa(atoms, probe_radius=1.4)',
            description: 'Calculate SASA using Shrake-Rupley algorithm',
            example: [
                'import biotite.structure as struc',
                'atom_sasa = struc.sasa(structure, probe_radius=1.4)',
                'print(f"Total SASA: {atom_sasa.sum():.2f} A^2")'
            ]
        },
        notes: 'Biotite returns a numpy array of per-atom SASA values. Rosetta provides both per-atom and per-residue. Biotite API is simpler.',
        equivalence: 'full'
    },
    {
        id: 'superimpose',
        category: 'Structural Alignment',
        rosetta: {
            name: 'SuperimposeMover',
            module: 'pyrosetta.rosetta.protocols.simple_moves',
            signature: 'SuperimposeMover(reference_pose)',
            description: 'Superimpose a pose onto a reference structure',
            example: [
                'from pyrosetta.rosetta.protocols.simple_moves import SuperimposeMover',
                'mover = SuperimposeMover(reference_pose)',
                'mover.apply(mobile_pose)'
            ]
        },
        biotite: {
            name: 'superimpose',
            module: 'biotite.structure',
            signature: 'biotite.structure.superimpose(fixed, mobile, atom_mask=None)',
            description: 'Superimpose mobile structure onto fixed structure',
            example: [
                'import biotite.structure as struc',
                'fitted, transformation = struc.superimpose(fixed_structure, mobile_structure)',
                'rmsd = struc.rmsd(fixed_structure, fitted)'
            ]
        },
        notes: 'Biotite returns the transformed structure and the transformation matrix. Rosetta modifies the pose in-place.',
        equivalence: 'full'
    },
    {
        id: 'rmsd',
        category: 'Structural Comparison',
        rosetta: {
            name: 'all_atom_rmsd',
            module: 'pyrosetta.rosetta.core.scoring',
            signature: 'all_atom_rmsd(pose1, pose2)',
            description: 'Calculate RMSD between two poses',
            example: [
                'from pyrosetta.rosetta.core.scoring import all_atom_rmsd',
                'rmsd_value = all_atom_rmsd(pose1, pose2)',
                'print(f"RMSD: {rmsd_value:.3f} A")'
            ]
        },
        biotite: {
            name: 'rmsd',
            module: 'biotite.structure',
            signature: 'biotite.structure.rmsd(reference, subject)',
            description: 'Calculate RMSD between two AtomArrays',
            example: [
                'import biotite.structure as struc',
                'rmsd_value = struc.rmsd(reference, subject)',
                'print(f"RMSD: {rmsd_value:.3f} A")'
            ]
        },
        notes: 'Both compute all-atom RMSD. Structures should be pre-aligned for meaningful comparison.',
        equivalence: 'full'
    },
    {
        id: 'secondary_structure',
        category: 'Structure Analysis',
        rosetta: {
            name: 'DsspMover',
            module: 'pyrosetta.rosetta.protocols.moves',
            signature: 'DsspMover().apply(pose)',
            description: 'Assign secondary structure using DSSP',
            example: [
                'from pyrosetta.rosetta.protocols.moves import DsspMover',
                'dssp = DsspMover()',
                'dssp.apply(pose)',
                'ss = pose.secstruct()'
            ]
        },
        biotite: {
            name: 'annotate_sse',
            module: 'biotite.structure',
            signature: 'biotite.structure.annotate_sse(atoms, chain_id)',
            description: 'Annotate secondary structure elements',
            example: [
                'import biotite.structure as struc',
                'sse = struc.annotate_sse(structure, chain_id="A")',
                '# Returns array: "a" = alpha helix, "b" = beta sheet, "c" = coil'
            ]
        },
        notes: 'Biotite uses its own P-SEA algorithm, not DSSP. Results may differ slightly from Rosetta/DSSP.',
        equivalence: 'partial'
    },
    {
        id: 'sequence_extraction',
        category: 'Sequence Analysis',
        rosetta: {
            name: 'pose.sequence()',
            module: 'pyrosetta',
            signature: 'pose.sequence()',
            description: 'Get the amino acid sequence from a pose',
            example: [
                'sequence = pose.sequence()',
                'print(f"Sequence: {sequence}")'
            ]
        },
        biotite: {
            name: 'get_residues',
            module: 'biotite.structure',
            signature: 'biotite.structure.get_residues(atoms)',
            description: 'Get residue identities from structure',
            example: [
                'import biotite.structure as struc',
                'residues = struc.get_residues(structure)',
                '# Returns (residue_ids, residue_names) tuple'
            ]
        },
        notes: 'Rosetta returns a single string. Biotite returns structured arrays of residue IDs and names.',
        equivalence: 'full'
    },
    {
        id: 'distance',
        category: 'Geometry',
        rosetta: {
            name: 'AtomPairConstraint / distance',
            module: 'pyrosetta.rosetta.core.scoring.constraints',
            signature: 'xyz1.distance(xyz2)',
            description: 'Calculate or constrain distances between atoms',
            example: [
                'xyz1 = pose.residue(1).xyz("CA")',
                'xyz2 = pose.residue(10).xyz("CA")',
                'dist = xyz1.distance(xyz2)'
            ]
        },
        biotite: {
            name: 'distance',
            module: 'biotite.structure',
            signature: 'biotite.structure.distance(atoms1, atoms2)',
            description: 'Calculate pairwise distances between atoms',
            example: [
                'import biotite.structure as struc',
                'ca_atoms = structure[structure.atom_name == "CA"]',
                'dist = struc.distance(ca_atoms[0], ca_atoms[9])',
                'print(f"CA-CA distance: {dist:.2f} A")'
            ]
        },
        notes: 'Biotite works directly on numpy coordinate arrays. Rosetta uses its own XYZ vector types.',
        equivalence: 'full'
    },
    {
        id: 'angle',
        category: 'Geometry',
        rosetta: {
            name: 'pose.phi / pose.psi / pose.omega',
            module: 'pyrosetta',
            signature: 'pose.phi(resid), pose.psi(resid), pose.omega(resid)',
            description: 'Measure backbone torsion angles',
            example: [
                'phi = pose.phi(5)',
                'psi = pose.psi(5)',
                'omega = pose.omega(5)',
                'print(f"Phi: {phi:.1f}, Psi: {psi:.1f}, Omega: {omega:.1f}")'
            ]
        },
        biotite: {
            name: 'dihedral',
            module: 'biotite.structure',
            signature: 'biotite.structure.dihedral(atoms1, atoms2, atoms3, atoms4)',
            description: 'Calculate dihedral angles between atoms',
            example: [
                'import biotite.structure as struc',
                'import numpy as np',
                'ca = structure[structure.atom_name == "CA"]',
                'dihedrals = struc.dihedral(ca[:-3], ca[1:-2], ca[2:-1], ca[3:])',
                'print(f"Dihedrals (deg): {np.rad2deg(dihedrals)}")'
            ]
        },
        notes: 'Biotite returns angles in radians. Rosetta provides convenient per-residue phi/psi/omega accessors.',
        equivalence: 'full'
    },
    {
        id: 'interface_analysis',
        category: 'Interface Analysis',
        rosetta: {
            name: 'InterfaceAnalyzerMover',
            module: 'pyrosetta.rosetta.protocols.analysis',
            signature: 'InterfaceAnalyzerMover(interface)',
            description: 'Comprehensive interface analysis (dSASA, dG, shape complementarity)',
            example: [
                'from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover',
                'iam = InterfaceAnalyzerMover("A_B")',
                'iam.apply(pose)',
                'print(f"dSASA: {iam.get_interface_delta_sasa()}")',
                'print(f"dG: {iam.get_interface_dG()}")'
            ]
        },
        biotite: {
            name: 'sasa + residue selection',
            module: 'biotite.structure',
            signature: 'Compute SASA of complex vs isolated chains',
            description: 'Calculate interface buried surface area via SASA difference',
            example: [
                'import biotite.structure as struc',
                'chain_a = structure[structure.chain_id == "A"]',
                'chain_b = structure[structure.chain_id == "B"]',
                'sasa_complex = struc.sasa(structure)',
                'sasa_a = struc.sasa(chain_a)',
                'sasa_b = struc.sasa(chain_b)',
                'delta_sasa = sasa_a.sum() + sasa_b.sum() - sasa_complex.sum()',
                'print(f"Buried surface area: {delta_sasa:.2f} A^2")'
            ]
        },
        notes: 'Biotite can compute dSASA but lacks Rosetta\'s dG, shape complementarity, and other interface metrics.',
        equivalence: 'partial'
    },
    {
        id: 'pdb_fetch',
        category: 'Database Access',
        rosetta: {
            name: 'toolbox.rcsb.pose_from_rcsb',
            module: 'pyrosetta.toolbox',
            signature: 'pyrosetta.toolbox.rcsb.pose_from_rcsb("1UBQ")',
            description: 'Fetch structure from RCSB PDB',
            example: [
                'from pyrosetta.toolbox import rcsb',
                'pose = rcsb.pose_from_rcsb("1UBQ")'
            ]
        },
        biotite: {
            name: 'rcsb.fetch',
            module: 'biotite.database.rcsb',
            signature: 'biotite.database.rcsb.fetch(pdb_ids, format, target_path)',
            description: 'Fetch structures from RCSB PDB',
            example: [
                'import biotite.database.rcsb as rcsb',
                'from biotite.structure.io.pdb import PDBFile',
                'pdb_path = rcsb.fetch("1UBQ", "pdb", target_path=".")',
                'pdb_file = PDBFile.read(pdb_path)',
                'structure = pdb_file.get_structure(model=1)'
            ]
        },
        notes: 'Both can fetch from RCSB. Biotite also supports UniProt and NCBI Entrez databases.',
        equivalence: 'full'
    },
    {
        id: 'residue_selection',
        category: 'Atom/Residue Selection',
        rosetta: {
            name: 'ResidueSelectors (ChainSelector, etc.)',
            module: 'pyrosetta.rosetta.core.select.residue_selector',
            signature: 'ChainSelector("A")',
            description: 'Select residues by chain, property, index, etc.',
            example: [
                'from pyrosetta.rosetta.core.select.residue_selector import ChainSelector',
                'selector = ChainSelector("A")',
                'selection = selector.apply(pose)',
                '# Returns a boolean vector'
            ]
        },
        biotite: {
            name: 'Array indexing / boolean masks',
            module: 'biotite.structure',
            signature: 'structure[structure.chain_id == "A"]',
            description: 'Select atoms using numpy-style boolean indexing',
            example: [
                '# Select by chain',
                'chain_a = structure[structure.chain_id == "A"]',
                '# Select by residue name',
                'ala_residues = structure[structure.res_name == "ALA"]',
                '# Select CA atoms',
                'ca_atoms = structure[structure.atom_name == "CA"]',
                '# Combine conditions',
                'chain_a_ca = structure[(structure.chain_id == "A") & (structure.atom_name == "CA")]'
            ]
        },
        notes: 'Biotite uses numpy-style boolean indexing which is often more concise. Rosetta selectors are composable objects with And/Or/Not combinators.',
        equivalence: 'full'
    },
    {
        id: 'scoring',
        category: 'Energy Scoring',
        rosetta: {
            name: 'ScoreFunction (ref2015)',
            module: 'pyrosetta.rosetta.core.scoring',
            signature: 'get_score_function()(pose)',
            description: 'Full physics-based energy scoring of protein structures',
            example: [
                'from pyrosetta.rosetta.core.scoring import get_score_function',
                'sfxn = get_score_function("ref2015")',
                'score = sfxn(pose)',
                'print(f"Total score: {score:.3f} REU")'
            ]
        },
        biotite: null,
        notes: 'Biotite is an analysis library and does NOT provide energy scoring functions. There is no Biotite equivalent for Rosetta energy calculations.',
        equivalence: 'none_from_biotite'
    },
    {
        id: 'fast_relax',
        category: 'Structure Optimization',
        rosetta: {
            name: 'FastRelax',
            module: 'pyrosetta.rosetta.protocols.relax',
            signature: 'FastRelax(scorefxn, rounds)',
            description: 'All-atom relaxation/energy minimization',
            example: [
                'from pyrosetta.rosetta.protocols.relax import FastRelax',
                'from pyrosetta.rosetta.core.scoring import get_score_function',
                'relax = FastRelax()',
                'relax.set_scorefxn(get_score_function("ref2015"))',
                'relax.apply(pose)'
            ]
        },
        biotite: null,
        notes: 'Biotite does NOT perform structure optimization, relaxation, or energy minimization. These are Rosetta-specific capabilities.',
        equivalence: 'none_from_biotite'
    },
    {
        id: 'fast_design',
        category: 'Protein Design',
        rosetta: {
            name: 'FastDesign',
            module: 'pyrosetta.rosetta.protocols.denovo_design.movers',
            signature: 'FastDesign(scorefxn, rounds)',
            description: 'Combinatorial protein sequence design',
            example: [
                'from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign',
                'from pyrosetta.rosetta.core.scoring import get_score_function',
                'design = FastDesign()',
                'design.set_scorefxn(get_score_function("ref2015"))',
                'design.apply(pose)'
            ]
        },
        biotite: null,
        notes: 'Biotite does NOT perform protein design. Biotite is analysis-only; Rosetta/PyRosetta is needed for design tasks.',
        equivalence: 'none_from_biotite'
    },
    {
        id: 'contact_map',
        category: 'Contact Analysis',
        rosetta: { name: 'ContactMap / distance matrices', module: 'pyrosetta.rosetta.core.scoring', signature: 'Distance-based contact detection', description: 'Compute residue contact maps and contact order', example: ['for i in range(1, pose.total_residue()+1):', '    for j in range(i+1, pose.total_residue()+1):', '        d = pose.residue(i).xyz("CA").distance(pose.residue(j).xyz("CA"))', '        if d < 8.0: print(f"Contact: {i}-{j}")'] },
        biotite: { name: 'CellList / distance', module: 'biotite.structure', signature: 'biotite.structure.CellList(atoms, cell_size)', description: 'Efficient spatial neighbor search and contact detection', example: ['import biotite.structure as struc', 'ca = structure[structure.atom_name == "CA"]', 'cell_list = struc.CellList(ca, cell_size=8.0)', 'contacts = cell_list.get_atoms(ca.coord[0], radius=8.0)'] },
        notes: 'Biotite CellList provides efficient O(n) spatial lookups. Rosetta contact analysis is integrated with its energy framework.',
        equivalence: 'partial'
    },
    {
        id: 'ramachandran',
        category: 'Structure Validation',
        rosetta: { name: 'RamachandranEnergy / pose.phi/psi', module: 'pyrosetta.rosetta.core.scoring', signature: 'pose.phi(resid), pose.psi(resid)', description: 'Ramachandran analysis via backbone torsion angles', example: ['for i in range(1, pose.total_residue()+1):', '    phi = pose.phi(i)', '    psi = pose.psi(i)', '    print(f"Res {i}: phi={phi:.1f} psi={psi:.1f}")'] },
        biotite: { name: 'dihedral_backbone', module: 'biotite.structure', signature: 'biotite.structure.dihedral_backbone(atoms)', description: 'Calculate phi/psi dihedral angles for Ramachandran analysis', example: ['import biotite.structure as struc', 'import numpy as np', 'phi, psi, omega = struc.dihedral_backbone(structure)', 'phi_deg = np.rad2deg(phi)', 'psi_deg = np.rad2deg(psi)'] },
        notes: 'Rosetta integrates Ramachandran scoring into its energy function. Biotite provides raw dihedral calculation.',
        equivalence: 'partial'
    },
    {
        id: 'hydrogen_bonds',
        category: 'Hydrogen Bond Analysis',
        rosetta: { name: 'HBondSet', module: 'pyrosetta.rosetta.core.scoring.hbonds', signature: 'HBondSet(pose)', description: 'Identify and analyze hydrogen bonds', example: ['from pyrosetta.rosetta.core.scoring.hbonds import HBondSet', 'hbond_set = HBondSet()', 'pose.update_residue_neighbors()', 'hbond_set.setup_for_residue_pair_energies(pose, False, False)', 'print(f"Total H-bonds: {hbond_set.nhbonds()}")'] },
        biotite: { name: 'hbond', module: 'biotite.structure', signature: 'biotite.structure.hbond(atoms)', description: 'Identify hydrogen bonds using geometric criteria', example: ['import biotite.structure as struc', 'triplets = struc.hbond(structure)', '# Returns (donor, hydrogen, acceptor) index triplets', 'print(f"Found {len(triplets)} H-bonds")'] },
        notes: 'Biotite uses geometric criteria. Rosetta uses energy-based detection integrated with scoring.',
        equivalence: 'partial'
    },
    {
        id: 'b_factor',
        category: 'Structure Properties',
        rosetta: { name: 'pose.pdb_info().bfactor', module: 'pyrosetta', signature: 'pose.pdb_info().bfactor(res, atom)', description: 'Access B-factors from PDB data', example: ['pdb_info = pose.pdb_info()', 'for i in range(1, pose.total_residue()+1):', '    bf = pdb_info.bfactor(i, 1)', '    print(f"Res {i}: B-factor = {bf:.2f}")'] },
        biotite: { name: 'AtomArray.b_factor', module: 'biotite.structure', signature: 'structure.b_factor', description: 'Access B-factors directly from AtomArray attribute', example: ['b_factors = structure.b_factor', 'print(f"Mean B-factor: {b_factors.mean():.2f}")', 'print(f"Max B-factor: {b_factors.max():.2f}")'] },
        notes: 'Biotite provides direct numpy array access. Rosetta accesses them per residue/atom through PdbInfo.',
        equivalence: 'full'
    },
    {
        id: 'center_of_mass',
        category: 'Geometry',
        rosetta: { name: 'center_of_mass', module: 'pyrosetta.rosetta.core.pose', signature: 'center_of_mass(pose, start, end)', description: 'Calculate center of mass of a pose or residue range', example: ['from pyrosetta.rosetta.core.pose import center_of_mass', 'com = center_of_mass(pose, 1, pose.total_residue())', 'print(f"COM: ({com.x:.2f}, {com.y:.2f}, {com.z:.2f}")'] },
        biotite: { name: 'mass_center', module: 'biotite.structure', signature: 'biotite.structure.mass_center(atoms)', description: 'Calculate the mass-weighted center of an AtomArray', example: ['import biotite.structure as struc', 'com = struc.mass_center(structure)', 'print(f"COM: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f}")'] },
        notes: 'Both compute mass-weighted centers. Biotite returns numpy array, Rosetta returns xyzVector.',
        equivalence: 'full'
    }
];

// Data-driven mapping for XML element -> PyRosetta class translation
const XML_TO_PYROSETTA_MAPPINGS = {
    // Movers - with expanded attribute coverage
    'FastRelax': { category: 'mover', pyclass: 'FastRelax', module: 'pyrosetta.rosetta.protocols.relax', attrs: {
        scorefxn: { setter: 'set_scorefxn', type: 'scorefxn' },
        repeats: { setter: 'set_default_repeats', type: 'int' },
        disable_design: { setter: 'set_enable_design', type: 'bool_invert' },
        cartesian: { setter: 'cartesian', type: 'bool_call' },
        min_type: { setter: 'min_type', type: 'string' },
        ramp_down_constraints: { setter: 'ramp_down_constraints', type: 'bool_call' }
    }},
    'FastDesign': { category: 'mover', pyclass: 'FastDesign', module: 'pyrosetta.rosetta.protocols.denovo_design.movers', attrs: {
        scorefxn: { setter: 'set_scorefxn', type: 'scorefxn' },
        repeats: { setter: 'set_default_repeats', type: 'int' },
        cartesian: { setter: 'cartesian', type: 'bool_call' },
        min_type: { setter: 'min_type', type: 'string' }
    }},
    'MinMover': { category: 'mover', pyclass: 'MinMover', module: 'pyrosetta.rosetta.protocols.minimization_packing', attrs: {
        scorefxn: { setter: 'score_function', type: 'scorefxn' },
        type: { setter: 'min_type', type: 'string' },
        tolerance: { setter: 'tolerance', type: 'float' },
        max_iter: { setter: 'max_iter', type: 'int' },
        cartesian: { setter: 'cartesian', type: 'bool_call' }
    }},
    'PackRotamersMover': { category: 'mover', pyclass: 'PackRotamersMover', module: 'pyrosetta.rosetta.protocols.minimization_packing', attrs: {
        scorefxn: { setter: 'score_function', type: 'scorefxn' },
        nloop: { setter: 'nloop', type: 'int' }
    }},
    'ShearMover': { category: 'mover', pyclass: 'ShearMover', module: 'pyrosetta.rosetta.protocols.simple_moves', attrs: {
        nmoves: { setter: 'nmoves', type: 'int' },
        temperature: { setter: 'temperature', type: 'float' }
    }},
    'SmallMover': { category: 'mover', pyclass: 'SmallMover', module: 'pyrosetta.rosetta.protocols.simple_moves', attrs: {
        nmoves: { setter: 'nmoves', type: 'int' },
        temperature: { setter: 'temperature', type: 'float' }
    }},
    'SavePoseMover': { category: 'mover', pyclass: 'SavePoseMover', module: 'pyrosetta.rosetta.protocols.moves', attrs: {} },
    'DumpPdb': { category: 'mover', pyclass: 'DumpPdb', module: 'pyrosetta.rosetta.protocols.moves', attrs: { fname: { setter: 'fname', type: 'string' } } },
    'InterfaceAnalyzerMover': { category: 'mover', pyclass: 'InterfaceAnalyzerMover', module: 'pyrosetta.rosetta.protocols.analysis', attrs: {
        interface: { setter: 'set_interface', type: 'string' }
    }},
    'AddCompositionConstraintMover': { category: 'mover', pyclass: 'AddCompositionConstraintMover', module: 'pyrosetta.rosetta.protocols.aa_composition', attrs: {} },
    'ClearConstraintsMover': { category: 'mover', pyclass: 'ClearConstraintsMover', module: 'pyrosetta.rosetta.protocols.constraint_movers', attrs: {} },
    // Filters - with expanded attributes
    'ScoreType': { category: 'filter', pyclass: 'ScoreTypeFilter', module: 'pyrosetta.rosetta.protocols.score_filters', attrs: { score_type: { setter: 'score_type', type: 'string' }, threshold: { setter: 'threshold', type: 'float' } } },
    'Ddg': { category: 'filter', pyclass: 'DdgFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: { threshold: { setter: 'threshold', type: 'float' }, repeats: { setter: 'repeats', type: 'int' } } },
    'Sasa': { category: 'filter', pyclass: 'SasaFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: { threshold: { setter: 'threshold', type: 'float' } } },
    'ShapeComplementarity': { category: 'filter', pyclass: 'ShapeComplementarityFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: {} },
    'BuriedUnsatHbonds': { category: 'filter', pyclass: 'BuriedUnsatHbondsFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: {} },
    'NetCharge': { category: 'filter', pyclass: 'NetChargeFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: {} },
    'Time': { category: 'filter', pyclass: 'TimeFilter', module: 'pyrosetta.rosetta.protocols.filters', attrs: {} },
    'ContactMolecularSurface': { category: 'filter', pyclass: 'ContactMolecularSurfaceFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: {} },
    'InterfaceBindingEnergyDensityFilter': { category: 'filter', pyclass: 'InterfaceBindingEnergyDensityFilter', module: 'pyrosetta.rosetta.protocols.simple_filters', attrs: {} },
    // Residue Selectors - with expanded attributes
    'Chain': { category: 'selector', pyclass: 'ChainSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: { chains: { setter: 'set_chain_strings', type: 'string_call' } } },
    'Neighborhood': { category: 'selector', pyclass: 'NeighborhoodResidueSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: { distance: { setter: 'set_distance', type: 'float_call' }, selector: { setter: 'set_focus_selector', type: 'selector_ref' } } },
    'And': { category: 'selector', pyclass: 'AndResidueSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: {} },
    'Not': { category: 'selector', pyclass: 'NotResidueSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: {} },
    'Or': { category: 'selector', pyclass: 'OrResidueSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: {} },
    'ResidueName': { category: 'selector', pyclass: 'ResidueNameSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: { residue_name3: { setter: 'residue_name3', type: 'string' } } },
    'ResidueIndex': { category: 'selector', pyclass: 'ResidueIndexSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: { resnums: { setter: 'set_index', type: 'string_call' } } },
    'SecondaryStructure': { category: 'selector', pyclass: 'SecondaryStructureSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: { ss: { setter: 'set_selected_ss', type: 'string_call' } } },
    'True': { category: 'selector', pyclass: 'TrueResidueSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: {} },
    'False': { category: 'selector', pyclass: 'FalseResidueSelector', module: 'pyrosetta.rosetta.core.select.residue_selector', attrs: {} },
    // Task Operations
    'RestrictToRepacking': { category: 'task_op', pyclass: 'RestrictToRepacking', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: {} },
    'ExtraRotamersGeneric': { category: 'task_op', pyclass: 'ExtraRotamersGeneric', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: { ex1: { setter: 'ex1', type: 'bool_call' }, ex2: { setter: 'ex2', type: 'bool_call' } } },
    'InitializeFromCommandline': { category: 'task_op', pyclass: 'InitializeFromCommandline', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: {} },
    'IncludeCurrent': { category: 'task_op', pyclass: 'IncludeCurrent', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: {} },
    'LimitAromaChi2': { category: 'task_op', pyclass: 'LimitAromaChi2', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: {} },
    'OperateOnResidueSubset': { category: 'task_op', pyclass: 'OperateOnResidueSubset', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: {} },
    'PreventRepackingRLT': { category: 'task_op', pyclass: 'PreventRepackingRLT', module: 'pyrosetta.rosetta.core.pack.task.operation', attrs: {} },
    // Child elements that need special handling (flagged as 'child' category)
    'MoveMap': { category: 'child', pyclass: 'MoveMap', module: 'pyrosetta.rosetta.core.kinematics', attrs: { bb: { setter: 'set_bb', type: 'bool_call' }, chi: { setter: 'set_chi', type: 'bool_call' }, jump: { setter: 'set_jump', type: 'bool_call' } } },
    'Span': { category: 'child', pyclass: null, module: null, attrs: { begin: { type: 'int' }, end: { type: 'int' }, bb: { type: 'bool' }, chi: { type: 'bool' } } },
    'Reweight': { category: 'child', pyclass: null, module: null, attrs: { scoretype: { type: 'string' }, weight: { type: 'float' } } },
    'ScoreFunction': { category: 'child', pyclass: null, module: null, attrs: { weights: { type: 'string' } } },
};

class RosettaMCPServer {
    constructor() {
        this.pythonPath = process.env.PYTHON_BIN || 'python3';
        this.serverPath = path.join(__dirname, 'rosetta_mcp_server.py');
    }

    async pythonEnvInfo() {
        return new Promise((resolve) => {
            const py = this.pythonPath;

            console.error('\x1b[36m[Python] Gathering environment information...\x1b[0m');
            
            const script = `
import json, sys, subprocess
info = {
  'python_executable': sys.executable,
  'python_version': sys.version,
}
try:
    out = subprocess.check_output([sys.executable, '-m', 'pip', 'list', '--format', 'json'], stderr=subprocess.STDOUT, text=True)
    info['pip_list'] = json.loads(out)
except Exception as e:
    info['pip_error'] = str(e)
print(json.dumps(info))
`;
            const proc = spawn(py, ['-c', script]);
            let output = '';
            let error = '';
            proc.stdout.on('data', d => { output += d.toString(); });
            proc.stderr.on('data', d => { error += d.toString(); });
            proc.on('close', () => {
                try { 
                    const result = JSON.parse(output.trim() || '{}');
                    console.error(`\x1b[32m[Python] ✅ Environment info gathered (${result.pip_list ? result.pip_list.length : 0} packages)\x1b[0m`);
                    resolve(result);
                } catch (_) { 
                    console.error(`\x1b[31m[Python] ❌ Failed to gather environment info: ${error}\x1b[0m`);
                    resolve({ stdout: output, stderr: error }); 
                }
            });
        });
    }

    async checkPyRosetta() {
        return new Promise((resolve) => {
            const py = this.pythonPath;

            console.error('\x1b[36m[PyRosetta] Checking if PyRosetta is available...\x1b[0m');
            
            const script = `
import json
resp = {'available': False}
try:
    import pyrosetta
    resp['available'] = True
    resp['version'] = getattr(pyrosetta, '__version__', None)
except Exception as e:
    resp['error'] = str(e)
print(json.dumps(resp))
`;
            const proc = spawn(py, ['-c', script]);
            let output = '';
            let error = '';
            proc.stdout.on('data', d => { output += d.toString(); });
            proc.stderr.on('data', d => { error += d.toString(); });
            proc.on('close', () => {
                try { 
                    const result = JSON.parse(output.trim() || '{}');
                    if (result.available) {
                        console.error(`\x1b[32m[PyRosetta] ✅ Available${result.version ? ` (v${result.version})` : ''}\x1b[0m`);
                    } else {
                        console.error(`\x1b[31m[PyRosetta] ❌ Not available: ${result.error || 'Unknown error'}\x1b[0m`);
                    }
                    resolve(result);
                } catch (_) { 
                    console.error(`\x1b[31m[PyRosetta] ❌ Check failed: ${error}\x1b[0m`);
                    resolve({ stdout: output, stderr: error }); 
                }
            });
        });
    }

    async installPyRosettaViaInstaller({ silent = true } = {}) {
        return new Promise((resolve) => {
            const py = this.pythonPath;

            // Show warning about long installation time
            console.error('\x1b[33m⚠️  WARNING: PyRosetta installation can take 10-30 minutes on first run!\x1b[0m');
            console.error('\x1b[33m   This involves downloading and compiling large scientific libraries.\x1b[0m');
            console.error('\x1b[33m   Please be patient and do not interrupt the process.\x1b[0m\n');
            
            const cmd = `${silent ? 'silent=True' : 'silent=False'}`;
            const script = `
import json, sys, subprocess, time
result = {'ok': False, 'progress': []}

def log_progress(message):
    result['progress'].append({'time': time.time(), 'message': message})
    print(json.dumps({'type': 'progress', 'message': message}))

try:
    log_progress('Upgrading pip, setuptools, and wheel...')
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', 'pip', 'setuptools', 'wheel'])
    
    log_progress('Installing pyrosetta-installer...')
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pyrosetta-installer'])
    
    log_progress('Importing pyrosetta-installer...')
    import pyrosetta_installer as I
    
    log_progress('Starting PyRosetta installation (this will take 10-30 minutes)...')
    I.install_pyrosetta(${cmd}, skip_if_installed=False)
    
    log_progress('PyRosetta installation completed successfully!')
    result['ok'] = True
except Exception as e:
    result['error'] = str(e)
    log_progress(f'Installation failed: {str(e)}')

print(json.dumps({'type': 'final', 'result': result}))
`;
            const proc = spawn(py, ['-c', script]);
            let output = '';
            let error = '';
            
            proc.stdout.on('data', (d) => { 
                const data = d.toString();
                output += data;
                
                // Parse progress messages
                try {
                    const lines = data.split('\n').filter(line => line.trim());
                    for (const line of lines) {
                        try {
                            const parsed = JSON.parse(line);
                            if (parsed.type === 'progress') {
                                console.error(`\x1b[36m[PyRosetta Install] ${parsed.message}\x1b[0m`);
                            } else if (parsed.type === 'final') {
                                // Final result, don't log here
                            }
                        } catch (e) {
                            // Not JSON, ignore
                        }
                    }
                } catch (e) {
                    // Ignore parsing errors
                }
            });
            
            proc.stderr.on('data', (d) => { 
                error += d.toString();
                // Show stderr output in real-time
                console.error(`\x1b[31m[PyRosetta Install Error] ${d.toString()}\x1b[0m`);
            });
            
            proc.on('close', () => {
                try {
                    // Extract the final result from the last line
                    const lines = output.split('\n').filter(line => line.trim());
                    for (let i = lines.length - 1; i >= 0; i--) {
                        try {
                            const parsed = JSON.parse(lines[i]);
                            if (parsed.type === 'final') {
                                resolve(parsed.result);
                                return;
                            }
                        } catch (e) {
                            // Continue searching
                        }
                    }
                    // Fallback to old parsing method
                    resolve(JSON.parse(output.trim() || '{}'));
                } catch (e) {
                    resolve({ stdout: output, stderr: error });
                }
            });
        });
    }

    async findRosettaScripts({ exe_path } = {}) {
        try {
            const resolved = this.resolveRosettaScriptsPath(exe_path);
            return { resolved };
        } catch (e) {
            return { error: e.message };
        }
    }

    resolveRosettaScriptsPath(preferredPath) {

        const resolveFrom = (p) => {
            if (!p || !p.length) return null;
            try {
                const stat = fs.existsSync(p) ? fs.statSync(p) : null;
                if (stat && stat.isFile()) return p;
                if (stat && stat.isDirectory()) {
                    const files = fs.readdirSync(p);
                    const match = files.find(f => f.startsWith('rosetta_scripts'));
                    if (match) return path.join(p, match);
                }
            } catch (_) {}
            return null;
        };

        // 1) Preferred explicit path or directory
        const fromPreferred = resolveFrom(preferredPath);
        if (fromPreferred) return fromPreferred;

        // 2) Environment override: file or directory
        const fromEnv = resolveFrom(process.env.ROSETTA_BIN || '');
        if (fromEnv) return fromEnv;

        // 3) Known candidate directories
        const candidateDirs = [
            // Typical source build locations
            path.join(process.env.HOME || '', 'rosetta', 'main', 'source', 'bin'),
            '/opt/rosetta/main/source/bin',
            '/usr/local/rosetta/main/source/bin'
        ].filter(Boolean);

        for (const dir of candidateDirs) {
            try {
                const files = fs.readdirSync(dir);
                const match = files.find(f => f.startsWith('rosetta_scripts') && fs.existsSync(path.join(dir, f)));
                if (match) {
                    return path.join(dir, match);
                }
            } catch (_) {
                // ignore
            }
        }

        // Fallback to assuming it's on PATH
        return 'rosetta_scripts';
    }

    async getRosettaInfo() {
        return new Promise((resolve, reject) => {
            const proc = spawn(this.pythonPath, [this.serverPath, 'info']);
            let output = '';
            let error = '';

            proc.stdout.on('data', (data) => {
                output += data.toString();
            });

            proc.stderr.on('data', (data) => {
                error += data.toString();
            });

            proc.on('close', (code) => {
                if (code === 0) {
                    try {
                        const result = JSON.parse(output);
                        resolve(result);
                    } catch (e) {
                        reject(new Error('Failed to parse server output'));
                    }
                } else {
                    reject(new Error(`Server exited with code ${code}: ${error}`));
                }
            });
        });
    }

    async getRosettaHelp(topic = null) {
        // Alias map: common concepts -> specific Rosetta names agents can look up
        const TOPIC_ALIASES = {
            'relax': 'FastRelax', 'relaxation': 'FastRelax', 'minimize': 'MinMover', 'minimization': 'MinMover',
            'design': 'FastDesign', 'interface design': 'FastDesign', 'protein design': 'FastDesign',
            'pack': 'PackRotamersMover', 'packing': 'PackRotamersMover', 'repack': 'PackRotamersMover',
            'dock': 'DockingProtocol', 'docking': 'DockingProtocol',
            'constraints': 'constraint_files', 'constraint': 'constraint_files',
            'loop': 'LoopModeler', 'loop modeling': 'LoopModeler', 'loops': 'LoopModeler',
            'mutate': 'MutateResidue', 'mutation': 'MutateResidue',
            'residue selectors': 'ResidueSelectors', 'residue selector': 'ResidueSelectors', 'selectors': 'ResidueSelectors',
            'score function': 'score_functions', 'scoring': 'score_functions', 'ref2015': 'score_functions',
            'pose_from_file': 'input_options', 'input': 'input_options',
        };

        // Extended static help for common topics
        const EXTENDED_HELP = {
            'score_functions': 'Score functions evaluate protein structures using physics-based energy terms. The default is ref2015 (recommended for most tasks). Other options:\n- ref2015: Standard all-atom energy function (default, best for most uses)\n- ref2015_cart: For cartesian-space minimization\n- beta_nov16: Beta energy function with updated parameters\n- soft_rep: Soft repulsive for early-stage design (reduces clashes)\n- score12: Legacy function (not recommended)\n\nUsage in XML: <ScoreFunction name="sfxn" weights="ref2015"/>\nUsage in PyRosetta: sfxn = get_score_function("ref2015")\n\nKey energy terms: fa_atr (attractive), fa_rep (repulsive), fa_sol (solvation), fa_elec (electrostatics), hbond_* (hydrogen bonds), rama_prepro (Ramachandran), p_aa_pp (rotamer probability)',
            'movers': 'Movers modify protein structures. Common movers:\n- FastRelax: All-atom relaxation/minimization (most used)\n- FastDesign: Sequence design with relaxation\n- MinMover: Energy minimization only\n- PackRotamersMover: Side-chain repacking\n- MutateResidue: Change a residue identity\n- DockingProtocol: Protein-protein docking\n- LoopModeler: Loop modeling/refinement\n- InterfaceAnalyzerMover: Analyze protein interfaces\n\nUse pyrosetta_introspect to search for specific movers by name.',
            'filters': 'Filters evaluate structures and can halt protocols. Common filters:\n- Ddg: Binding free energy (interface stability)\n- Sasa: Solvent-accessible surface area\n- ShapeComplementarity: Interface shape match (sc > 0.65 is good)\n- BuriedUnsatHbonds: Buried unsatisfied H-bond donors/acceptors\n- ScoreType: Filter by any score term\n- NetCharge: Total charge of selection\n\nUse pyrosetta_introspect to search for specific filters.',
            'xml': 'RosettaScripts XML defines protocols with these sections:\n<ROSETTASCRIPTS>\n  <SCOREFXNS> - Define score functions\n  <RESIDUE_SELECTORS> - Select residues by chain, property, etc.\n  <TASKOPERATIONS> - Control packing behavior\n  <MOVERS> - Define structural modifications\n  <FILTERS> - Define quality checks\n  <PROTOCOLS> - Order of operations (<Add mover="x"/>)\n  <OUTPUT scorefxn="x"/> - Final scoring\n</ROSETTASCRIPTS>\n\nUse validate_xml to check syntax, xml_to_pyrosetta to convert to Python.',
            'parameters': 'Common command-line flags:\n- -in:file:s input.pdb  (input structure)\n- -out:path:all output/  (output directory)\n- -nstruct 10  (number of structures)\n- -parser:protocol protocol.xml  (XML protocol file)\n- -ex1 -ex2aro  (extra rotamers)\n- -use_input_sc  (keep input side chains)\n- -overwrite  (overwrite existing outputs)\n- -score:symmetric_gly_tables true\n- -nblist_autoupdate true\n- -holes:dalphaball /path/to/DAlphaBall',
        };

        // Check extended static help first
        const lowerTopic = (topic || '').toLowerCase();
        if (topic && EXTENDED_HELP[lowerTopic]) {
            return EXTENDED_HELP[lowerTopic];
        }
        if (topic && EXTENDED_HELP[topic]) {
            return EXTENDED_HELP[topic];
        }

        // Resolve aliases
        const resolvedTopic = TOPIC_ALIASES[lowerTopic] || topic;

        // Check extended help with resolved topic too
        const resolvedLower = resolvedTopic.toLowerCase();
        if (EXTENDED_HELP[resolvedLower]) {
            return EXTENDED_HELP[resolvedLower];
        }

        // If no topic given, return available topics
        if (!topic) {
            return 'Available topics: score_functions, movers, filters, xml, parameters\n\nYou can also ask for specific mover/filter names (e.g., FastRelax, Ddg, ChainSelector) or concepts (relax, docking, constraints, design, loop modeling, mutate).';
        }

        // Try fetching live docs with resolved topic
        const searchResult = await this.searchRosettaWebDocs({ query: resolvedTopic, max_results: 1 });
        if (searchResult.results && searchResult.results.length > 0) {
            const docResult = await this.getRosettaWebDoc({ url: searchResult.results[0].url, max_chars: 6000 });
            if (docResult.text && docResult.text.length > 200) {
                return `# ${docResult.title || resolvedTopic}\nSource: ${docResult.url}\n\n${docResult.text}`;
            }
        }

        // Fallback with helpful guidance
        return `No detailed help found for "${topic}"${resolvedTopic !== topic ? ` (resolved to "${resolvedTopic}")` : ''}. Try:\n- pyrosetta_introspect with query="${resolvedTopic}" for live API details\n- get_rosetta_web_doc with a specific rosettacommons.org URL`;
    }

    async validateXML(xmlContent, validateAgainstSchema = false) {
        return new Promise((resolve, reject) => {
            const schemaPath = path.join(process.cwd(), 'docs_cache', 'rosetta_scripts_schema.xsd');
            const script = `
import json, sys, os
import xml.etree.ElementTree as ET

try:
    xml_content = sys.stdin.read()
    check_schema = os.environ.get('MCP_CHECK_SCHEMA', '') == '1'
    schema_path = os.environ.get('MCP_SCHEMA_PATH', '')

    result = {}
    root = ET.fromstring(xml_content)
    result['valid'] = True
    result['message'] = 'XML is well-formed'
    result['root_tag'] = root.tag

    xml_elements = set()
    for elem in root.iter():
        xml_elements.add(elem.tag)
    result['element_count'] = len(xml_elements)
    result['elements_found'] = sorted(xml_elements)

    if check_schema:
        if schema_path and os.path.exists(schema_path):
            try:
                schema_tree = ET.parse(schema_path)
                schema_root = schema_tree.getroot()
                schema_names = set()
                skip = {'schema','annotation','element','complexType','sequence',
                        'choice','attribute','documentation','simpleType',
                        'restriction','enumeration','group','extension',
                        'simpleContent','complexContent','all','any',
                        'attributeGroup','union','list'}
                for elem in schema_root.iter():
                    name_val = elem.get('name')
                    if name_val and name_val not in skip:
                        schema_names.add(name_val)

                standard_sections = {'ROSETTASCRIPTS','SCOREFXNS','RESIDUE_SELECTORS',
                                     'TASKOPERATIONS','MOVERS','FILTERS','PROTOCOLS',
                                     'OUTPUT','APPLY_TO_POSE','IMPORT','RESOURCES',
                                     'Add','SIMPLE_METRICS','ScoreFunction'}
                unknown = [e for e in xml_elements if e not in standard_sections and e not in schema_names]

                result['schema_validation'] = {
                    'checked': True,
                    'schema_element_count': len(schema_names),
                    'unknown_elements': unknown,
                    'all_known': len(unknown) == 0
                }
                if unknown:
                    result['message'] = f'XML is well-formed but {len(unknown)} element(s) not found in schema: {", ".join(unknown)}'
                else:
                    result['message'] = 'XML is well-formed and all elements match the Rosetta schema'
            except Exception as se:
                result['schema_validation'] = {'checked': False, 'error': f'Schema parsing failed: {str(se)}'}
        else:
            result['schema_validation'] = {'checked': False, 'error': 'No cached schema found. Run rosetta_scripts_schema first.'}

except ET.ParseError as e:
    result = {'valid': False, 'error': str(e)}
except Exception as e:
    result = {'valid': False, 'error': str(e)}

print(json.dumps(result))
`;
            const proc = spawn(this.pythonPath, ['-c', script], {
                env: {
                    ...process.env,
                    MCP_CHECK_SCHEMA: validateAgainstSchema ? '1' : '',
                    MCP_SCHEMA_PATH: schemaPath
                }
            });
            proc.stdin.write(xmlContent);
            proc.stdin.end();
            let output = '';
            let error = '';

            proc.stdout.on('data', (d) => { output += d.toString(); });
            proc.stderr.on('data', (d) => { error += d.toString(); });

            proc.on('close', (code) => {
                try {
                    resolve(JSON.parse(output.trim() || '{}'));
                } catch (e) {
                    reject(new Error('Failed to parse validation result'));
                }
            });
        });
    }

    async runRosettaScripts({ exe_path, xml_path, input_pdb, out_dir, extra_flags }) {
        return new Promise((resolve, reject) => {
            if (!xml_path || !input_pdb || !out_dir) {
                reject(new Error('xml_path, input_pdb, and out_dir are required'));
                return;
            }

            const rosettaExe = this.resolveRosettaScriptsPath(exe_path);
            try {
                if (!fs.existsSync(out_dir)) {
                    fs.mkdirSync(out_dir, { recursive: true });
                }
            } catch (e) {
                reject(new Error(`Failed to ensure out_dir exists: ${e.message}`));
                return;
            }

            const args = [
                '-in:file:s', input_pdb,
                '-out:path:all', out_dir,
                '-parser:protocol', xml_path
            ];

            if (Array.isArray(extra_flags)) {
                for (const flag of extra_flags) {
                    if (typeof flag === 'string' && flag.trim().length > 0) {
                        // Split on whitespace to allow flags with values, e.g., "-nstruct 5"
                        const parts = flag.split(/\s+/).filter(Boolean);
                        args.push(...parts);
                    }
                }
            }

            const proc = spawn(rosettaExe, args, { env: process.env });
            let stdout = '';
            let stderr = '';

            proc.stdout.on('data', d => { stdout += d.toString(); });
            proc.stderr.on('data', d => { stderr += d.toString(); });
            proc.on('error', (e) => {
                reject(new Error(`Failed to start rosetta_scripts: ${e.message}`));
            });
            proc.on('close', (code) => {
                resolve({ exit_code: code, stdout, stderr, out_dir });
            });
        });
    }

    async pyrosettaScore({ pdb_path, scorefxn, per_residue = false }) {
        if (!pdb_path) {
            throw new Error('pdb_path is required');
        }

        // Check if PyRosetta is available
        const pyrosettaStatus = await this.checkPyRosetta();
        if (!pyrosettaStatus.available) {
            console.error('\x1b[33m⚠️  PyRosetta not found. Auto-installing...\x1b[0m');
            console.error('\x1b[33m   This will take 10-30 minutes. Please be patient.\x1b[0m\n');
            const installResult = await this.installPyRosettaViaInstaller({ silent: false });
            if (!installResult.ok) {
                throw new Error(`PyRosetta installation failed: ${installResult.error || 'Unknown error'}`);
            }
            console.error('\x1b[32m✅ PyRosetta installed successfully! Retrying score operation...\x1b[0m\n');
        }

        return new Promise((resolve) => {
            const py = this.pythonPath;
            const script = `
import json, os, sys
# Redirect stdout to stderr during init to prevent banner from polluting JSON output
_orig_stdout = sys.stdout
sys.stdout = sys.stderr
try:
    import pyrosetta
except Exception as e:
    sys.stdout = _orig_stdout
    print(json.dumps({"error": f"PyRosetta not available: {str(e)}"}))
    raise SystemExit(0)

pyrosetta.init('-mute all')
sys.stdout = _orig_stdout

from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory, get_score_function
pdb_path = os.environ['MCP_PDB_PATH']
if not os.path.exists(pdb_path):
    print(json.dumps({"error": f"File not found: {pdb_path}"}))
    raise SystemExit(0)
sfxn_name = os.environ.get('MCP_SCOREFXN', '')
per_residue = os.environ.get('MCP_PER_RESIDUE', '') == '1'
try:
    pose = pose_from_pdb(pdb_path)
    sfxn = ScoreFunctionFactory.create_score_function(sfxn_name) if sfxn_name else get_score_function()
    score = sfxn(pose)
except Exception as e:
    print(json.dumps({"error": f"Scoring failed: {str(e)}"}))
    raise SystemExit(0)
result = {"score": float(score)}
if per_residue:
    energies = []
    for i in range(1, pose.total_residue() + 1):
        energies.append({
            'residue': i,
            'name': pose.residue(i).name3(),
            'score': float(pose.energies().residue_total_energy(i))
        })
    result['per_residue'] = energies
    result['residue_count'] = len(energies)
print(json.dumps(result))
`;
            const proc = spawn(py, ['-c', script], {
                env: { ...process.env, MCP_PDB_PATH: pdb_path, MCP_SCOREFXN: scorefxn || '', MCP_PER_RESIDUE: per_residue ? '1' : '' }
            });
            let output = '';
            let error = '';
            proc.stdout.on('data', d => { output += d.toString(); });
            proc.stderr.on('data', d => { error += d.toString(); });
            proc.on('close', (code) => {
                try {
                    const parsed = JSON.parse(output.trim() || '{}');
                    if (parsed && typeof parsed === 'object') {
                        resolve(parsed);
                    } else {
                        resolve({ exit_code: code, stdout: output, stderr: error });
                    }
                } catch (e) {
                    resolve({ exit_code: code, stdout: output, stderr: error });
                }
            });
        });
    }

    async pyrosettaIntrospect({ query, kind, max_results }) {
        const py = this.pythonPath;

        // Check if PyRosetta is available
        const pyrosettaStatus = await this.checkPyRosetta();
        if (!pyrosettaStatus.available) {
            console.error('\x1b[33m⚠️  PyRosetta not found. Auto-installing...\x1b[0m');
            console.error('\x1b[33m   This will take 10-30 minutes. Please be patient.\x1b[0m\n');
            const installResult = await this.installPyRosettaViaInstaller({ silent: false });
            if (!installResult.ok) {
                return { error: `PyRosetta installation failed: ${installResult.error || 'Unknown error'}` };
            }
            console.error('\x1b[32m✅ PyRosetta installed successfully! Retrying introspect operation...\x1b[0m\n');
        }

        const limit = Number.isInteger(max_results) && max_results > 0 ? max_results : 50;
        return new Promise((resolve) => {
            const script = `
import json, sys, os, importlib, inspect
_orig_stdout = sys.stdout
sys.stdout = sys.stderr
try:
    import pyrosetta
except Exception as e:
    sys.stdout = _orig_stdout
    print(json.dumps({"error": f"PyRosetta not available: {str(e)}"}))
    raise SystemExit(0)
sys.stdout = _orig_stdout

namespaces = [
    'pyrosetta.rosetta.protocols.moves',
    'pyrosetta.rosetta.protocols.simple_moves',
    'pyrosetta.rosetta.protocols.minimization_packing',
    'pyrosetta.rosetta.protocols.relax',
    'pyrosetta.rosetta.protocols.rosetta_scripts',
    'pyrosetta.rosetta.protocols.simple_filters',
    'pyrosetta.rosetta.protocols.filters',
    'pyrosetta.rosetta.core.select.residue_selector',
    'pyrosetta.rosetta.core.pack.task.operation',
]

q = os.environ.get('MCP_QUERY', '').lower()
kind = os.environ.get('MCP_KIND', '').lower()
max_results = int(os.environ.get('MCP_MAX_RESULTS', '50'))

def kind_allows(module_name: str) -> bool:
    if not kind:
        return True
    m = module_name.lower()
    if kind in ('mover','movers'):
        return '.protocols.' in m and ('moves' in m or 'simple_moves' in m or 'relax' in m or 'minimization_packing' in m)
    if kind in ('filter','filters'):
        return '.protocols.' in m and ('filters' in m or 'simple_filters' in m)
    if kind in ('selector','residue_selector','residue_selectors'):
        return 'residue_selector' in m
    if kind in ('task','task_operation','task_operations'):
        return '.task.operation' in m
    return True

results = []
for mod_name in namespaces:
    try:
        m = importlib.import_module(mod_name)
    except Exception:
        continue
    try:
        for name, obj in inspect.getmembers(m, inspect.isclass):
            if q and q not in name.lower():
                continue
            if not kind_allows(getattr(obj, '__module__', '')):
                continue
            entry = {
                'name': name,
                'module': getattr(obj, '__module__', ''),
                'bases': [b.__name__ for b in getattr(obj, '__mro__', [])[1:3] if hasattr(b,'__name__')],
            }
            try:
                doc = inspect.getdoc(obj)
                if doc:
                    entry['doc'] = doc[:2000]
                else:
                    try:
                        if hasattr(obj, '__init__'):
                            sig = inspect.signature(obj.__init__)
                            entry['doc'] = f"Constructor signature: {sig}"
                        else:
                            entry['doc'] = f"Class {name} from {mod_name}"
                    except:
                        entry['doc'] = f"Class {name} from {mod_name}"
            except Exception:
                entry['doc'] = f"Class {name} from {mod_name}"
            try:
                sig = str(inspect.signature(obj.__init__))
                entry['init'] = sig
            except Exception:
                pass
            results.append(entry)
            if len(results) >= max_results:
                raise StopIteration
    except StopIteration:
        break

print(json.dumps({'results': results, 'count': len(results)}))
`;
            const proc = spawn(py, ['-c', script], {
                env: {
                    ...process.env,
                    MCP_QUERY: query || '',
                    MCP_KIND: (kind || '').toString().toLowerCase(),
                    MCP_MAX_RESULTS: String(limit)
                }
            });
            let output = '';
            let error = '';
            proc.stdout.on('data', d => { output += d.toString(); });
            proc.stderr.on('data', d => { error += d.toString(); });
            proc.on('close', () => {
                try {
                    const parsed = JSON.parse(output.trim() || '{}');
                    resolve(parsed);
                } catch (_) {
                    resolve({ stdout: output, stderr: error });
                }
            });
        });
    }

    async rosettaScriptsSchema({ exe_path, cache_dir, extract_elements }) {
        return new Promise((resolve, reject) => {
            const rosettaExe = this.resolveRosettaScriptsPath(exe_path);
            const baseDir = cache_dir && cache_dir.length ? cache_dir : path.join(process.cwd(), 'docs_cache');
            try {
                if (!fs.existsSync(baseDir)) fs.mkdirSync(baseDir, { recursive: true });
            } catch (e) {
                reject(new Error(`Failed to ensure cache_dir exists: ${e.message}`));
                return;
            }
            const schemaPath = path.join(baseDir, 'rosetta_scripts_schema.xsd');
            const proc = spawn(rosettaExe, ['-parser:output_schema', schemaPath]);
            let stdout = '';
            let stderr = '';
            proc.stdout.on('data', d => { stdout += d.toString(); });
            proc.stderr.on('data', d => { stderr += d.toString(); });
            proc.on('error', (e) => reject(new Error(`Failed to start rosetta_scripts: ${e.message}`)));
            proc.on('close', (code) => {
                if (code !== 0) {
                    resolve({ exit_code: code, stdout, stderr });
                    return;
                }
                try {
                    const content = fs.readFileSync(schemaPath, 'utf8');
                    let elements = undefined;
                    if (extract_elements) {
                        const matches = [...content.matchAll(/\bname=\"([^\"]+)\"/g)].map(m => m[1]);
                        // de-duplicate and filter common XML names
                        const skip = new Set(['schema','annotation','element','complexType','sequence','choice','attribute','documentation']);
                        elements = Array.from(new Set(matches)).filter(n => !skip.has(n)).slice(0, 1000);
                    }
                    resolve({ exit_code: 0, schema_path: schemaPath, size: content.length, elements });
                } catch (e) {
                    resolve({ exit_code: 0, error: `Schema generated but could not be read: ${e.message}`, schema_path: schemaPath });
                }
            });
        });
    }

    async cacheCliDocs({ exe_path, cache_dir }) {
        return new Promise((resolve, reject) => {
            const rosettaExe = this.resolveRosettaScriptsPath(exe_path);
            const baseDir = cache_dir && cache_dir.length ? cache_dir : path.join(process.cwd(), 'docs_cache');
            try { if (!fs.existsSync(baseDir)) fs.mkdirSync(baseDir, { recursive: true }); } catch (e) { reject(new Error(`Failed to ensure cache_dir exists: ${e.message}`)); return; }

            const runHelp = (args, outfile) => new Promise((res) => {
                const proc = spawn(rosettaExe, args);
                let output = '';
                let error = '';
                proc.stdout.on('data', d => { output += d.toString(); });
                proc.stderr.on('data', d => { error += d.toString(); });
                proc.on('close', () => {
                    try { fs.writeFileSync(path.join(baseDir, outfile), output + (error ? `\nSTDERR:\n${error}` : '')); } catch (e) { if (process.env.MCP_DEBUG) console.error(`[cache] write failed: ${e.message}`); }
                    res({ stdout: output, stderr: error });
                });
            });

            (async () => {
                const r1 = await runHelp(['-help'], 'help.txt');
                // Best-effort: try additional parser info flags; ignore failures
                const r2 = await runHelp(['-parser:info'], 'parser_info.txt');
                resolve({ saved: [path.join(baseDir, 'help.txt'), path.join(baseDir, 'parser_info.txt')] });
            })();
        });
    }

    async getCachedDocs({ cache_dir, query, max_lines }) {
        const baseDir = cache_dir && cache_dir.length ? cache_dir : path.join(process.cwd(), 'docs_cache');
        const files = ['help.txt', 'parser_info.txt'].map(f => path.join(baseDir, f));

        // Auto-cache: if no cache files exist, run cacheCliDocs first
        const anyExist = files.some(f => fs.existsSync(f));
        if (!anyExist) {
            try {
                await this.cacheCliDocs({ cache_dir: baseDir });
            } catch (e) {
                if (process.env.MCP_DEBUG) console.error(`[get_cached_docs] auto-cache failed: ${e.message}`);
            }
        }

        const q = (query || '').toLowerCase();
        try {
            let results = [];
            for (const file of files) {
                if (!fs.existsSync(file)) continue;
                const content = fs.readFileSync(file, 'utf8');
                const lines = content.split(/\r?\n/);
                for (let i = 0; i < lines.length; i++) {
                    const line = lines[i];
                    if (!q || line.toLowerCase().includes(q)) {
                        results.push({ file, line: i + 1, text: line });
                    }
                }
            }
            const limit = Number.isInteger(max_lines) && max_lines > 0 ? max_lines : 200;
            return { count: results.length, matches: results.slice(0, limit) };
        } catch (e) {
            throw new Error(`Failed to search cache: ${e.message}`);
        }
    }

    // Simple HTTPS fetch utility with redirect support
    async fetchUrl(url, maxRedirects = 5) {
        const https = require('https');
        const http = require('http');
        
        const fetchWithRedirect = (url, redirectCount = 0) => {
            return new Promise((resolve, reject) => {
                if (redirectCount >= maxRedirects) {
                    reject(new Error(`Too many redirects (${redirectCount})`));
                    return;
                }
                
                const isHttps = url.startsWith('https://');
                const client = isHttps ? https : http;
                
                try {
                    const req = client.get(url, { 
                        headers: { 
                            'User-Agent': 'rosetta-mcp/1.0 (+https://npmjs.com/package/rosetta-mcp-server)' 
                        } 
                    }, (res) => {
                        // Handle redirects
                        if (res.statusCode >= 300 && res.statusCode < 400 && res.headers.location) {
                            const newUrl = res.headers.location.startsWith('http') 
                                ? res.headers.location 
                                : `${isHttps ? 'https' : 'http'}://${new URL(url).host}${res.headers.location}`;
                            resolve(fetchWithRedirect(newUrl, redirectCount + 1));
                            return;
                        }
                        
                        let body = '';
                        res.on('data', (d) => { body += d.toString(); });
                        res.on('end', () => {
                            resolve({ status: res.statusCode, headers: res.headers, body });
                        });
                    });
                    req.on('error', (err) => reject(err));
                } catch (e) {
                    reject(e);
                }
            });
        };
        
        return fetchWithRedirect(url);
    }

    stripHtml(html) {
        if (!html) return '';
        // Remove scripts/styles then tags, collapse whitespace
        return html
            .replace(/<script[\s\S]*?<\/script>/gi, '')
            .replace(/<style[\s\S]*?<\/style>/gi, '')
            .replace(/<[^>]+>/g, ' ')
            .replace(/\s+/g, ' ')
            .trim();
    }

    async searchRosettaWebDocs({ query, max_results }) {
        const max = Number.isInteger(max_results) && max_results > 0 ? max_results : 3;
        const q = encodeURIComponent(`site:rosettacommons.org/docs ${query || ''}`.trim());
        // Use DuckDuckGo HTML endpoint (no JS) for simple scraping
        const searchUrl = `https://duckduckgo.com/html/?q=${q}`;
        let results = [];
        try {
            const resp = await this.fetchUrl(searchUrl);
            if (resp.status === 200) {
                const html = resp.body || '';
                const linkRegex = /<a[^>]+href="([^"]+)"[^>]*>([\s\S]*?)<\/a>/gi;
                let m;
                while ((m = linkRegex.exec(html)) !== null && results.length < max * 3) {
                    let href = m[1];
                    const text = this.stripHtml(m[2]);
                    const uddgMatch = href.match(/[?&]uddg=([^&]+)/);
                    if (uddgMatch) {
                        try { href = decodeURIComponent(uddgMatch[1]); } catch (_) {}
                    }
                    if (/rosettacommons\.org\/.*/i.test(href)) {
                        results.push({ title: text.slice(0, 150), url: href });
                    }
                }
            }
        } catch (e) {
            // DuckDuckGo failed, fall through to direct URL fallback
        }

        // Fallback: construct direct Rosetta docs URLs when search fails or returns nothing
        if (results.length === 0 && query) {
            const slug = query.replace(/\s+/g, '-').replace(/[^a-zA-Z0-9_-]/g, '');
            const slugNoHyphen = query.replace(/\s+/g, '').replace(/[^a-zA-Z0-9_]/g, '');
            const base = 'https://www.rosettacommons.org/docs/latest';
            const candidateUrls = [
                // Movers
                { title: `${query} (Mover)`, url: `${base}/scripting_documentation/RosettaScripts/Movers/movers_pages/${slugNoHyphen}Mover` },
                { title: `${query} (Mover)`, url: `${base}/scripting_documentation/RosettaScripts/Movers/movers_pages/${slugNoHyphen}` },
                // Filters
                { title: `${query} (Filter)`, url: `${base}/scripting_documentation/RosettaScripts/Filters/filter_pages/${slugNoHyphen}Filter` },
                { title: `${query} (Filter)`, url: `${base}/scripting_documentation/RosettaScripts/Filters/filter_pages/${slugNoHyphen}` },
                // Selectors
                { title: `${query} (Selector)`, url: `${base}/scripting_documentation/RosettaScripts/ResidueSelectors/residue_selectors/${slugNoHyphen}` },
                // Task Operations
                { title: `${query} (TaskOp)`, url: `${base}/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/${slugNoHyphen}` },
                // Concept/general pages
                { title: `${query} (Concept)`, url: `${base}/rosetta_basics/${slug}` },
                { title: `${query} (Scoring)`, url: `${base}/rosetta_basics/scoring/${slug}` },
                { title: `${query} (File formats)`, url: `${base}/rosetta_basics/file_types/${slug}` },
                { title: `${query} (Application)`, url: `${base}/application_documentation/${slug}` },
                { title: `${query} (RosettaScripts)`, url: `${base}/scripting_documentation/RosettaScripts/${slugNoHyphen}` },
                // With hyphenated slug
                { title: `${query}`, url: `${base}/rosetta_basics/${slug}` },
                { title: `${query}`, url: `${base}/scripting_documentation/RosettaScripts/Movers/movers_pages/${slug}` },
            ];
            // Probe each URL, keep those that return 200
            for (const candidate of candidateUrls) {
                if (results.length >= max) break;
                try {
                    const probe = await this.fetchUrl(candidate.url);
                    if (probe.status === 200 && probe.body && probe.body.length > 500) {
                        results.push(candidate);
                    }
                } catch (_) {}
            }
            if (results.length > 0) {
                return { query, results: results.slice(0, max), count: results.length, source: 'direct_url_probe' };
            }
            return { query, results: [], count: 0, error: 'Search engine unavailable and no direct doc URLs found. Try get_rosetta_web_doc with a specific URL.' };
        }

        // De-duplicate and cap
        const seen = new Set();
        const deduped = [];
        for (const r of results) {
            const key = r.url.split('#')[0];
            if (seen.has(key)) continue;
            seen.add(key);
            deduped.push(r);
            if (deduped.length >= max) break;
        }
        return { query, results: deduped, count: deduped.length };
    }

    async getRosettaWebDoc({ url, max_chars }) {
        const limit = Number.isInteger(max_chars) && max_chars > 0 ? max_chars : 4000;
        if (!url || typeof url !== 'string') return { error: 'url is required' };
        try {
            const resp = await this.fetchUrl(url);
            if (resp.status !== 200) {
                let suggestion = 'Check if the URL is correct and accessible';
                if (resp.status >= 300 && resp.status < 400) {
                    suggestion = 'URL may have redirected - try the search tool first';
                } else if (resp.status === 403) {
                    suggestion = 'Access denied - try using the search tool to find accessible URLs';
                } else if (resp.status === 404) {
                    suggestion = 'Page not found - URL may be outdated, try the search tool';
                }
                
                return { 
                    error: `Fetch failed with status ${resp.status}`, 
                    url,
                    suggestion
                };
            }
            const html = resp.body || '';
            const titleMatch = html.match(/<title>([\s\S]*?)<\/title>/i);
            const title = titleMatch ? this.stripHtml(titleMatch[1]).slice(0, 200) : undefined;
            const text = this.stripHtml(html).slice(0, limit);
            return { url, title, text, length: text.length };
        } catch (e) {
            return { error: `Fetch error: ${e.message}`, url };
        }
    }

    async xmlToPyRosetta({ xml_content, include_comments = true, output_format = 'python' }) {
        return new Promise((resolve) => {
            const py = this.pythonPath;
            const mappingsJson = JSON.stringify(XML_TO_PYROSETTA_MAPPINGS);
            const script = `
import json, sys, os
import xml.etree.ElementTree as ET

try:
    xml_content = sys.stdin.read()
    include_comments = os.environ.get('MCP_INCLUDE_COMMENTS', 'True') == 'True'
    mappings = json.loads(os.environ.get('MCP_XML_MAPPINGS', '{}'))
    root = ET.fromstring(xml_content)

    found = {'mover': [], 'filter': [], 'selector': [], 'task_op': [], 'child': []}
    unrecognized = []
    skip_tags = {'ROSETTASCRIPTS','SCOREFXNS','RESIDUE_SELECTORS','TASKOPERATIONS',
                 'MOVERS','FILTERS','PROTOCOLS','OUTPUT','APPLY_TO_POSE',
                 'IMPORT','RESOURCES','Add','SIMPLE_METRICS'}

    for elem in root.iter():
        tag = elem.tag
        if not isinstance(tag, str):
            continue
        if tag in skip_tags:
            continue
        if tag in mappings:
            m = mappings[tag]
            found[m['category']].append((tag, dict(elem.attrib), m))
        else:
            unrecognized.append(tag)

    lines = []
    if include_comments:
        lines.append('# Generated PyRosetta code from RosettaScripts XML')
        lines.append('# ' + '=' * 50)
        lines.append('')

    modules = set()
    for cat_items in found.values():
        for _, _, m in cat_items:
            if m.get('module'):
                modules.add(m['module'])

    lines.append('import pyrosetta')
    lines.append('from pyrosetta import pose_from_pdb')
    lines.append('from pyrosetta.rosetta.core.scoring import get_score_function')
    for mod in sorted(modules):
        lines.append(f'from {mod} import *')
    lines.append('')
    lines.append('pyrosetta.init("-mute all")')
    lines.append('')
    lines.append('pose = pose_from_pdb("your_protein.pdb")')
    lines.append('')

    category_labels = {'mover': 'Movers', 'filter': 'Filters', 'selector': 'Residue Selectors', 'task_op': 'Task Operations'}
    var_counters = {}

    for cat in ['selector', 'task_op', 'filter', 'mover']:
        items = found[cat]
        if not items:
            continue
        if include_comments:
            lines.append(f'# {category_labels[cat]}')
        for tag, attrs, m in items:
            pyclass = m['pyclass']
            base_var = pyclass[0].lower() + pyclass[1:]
            count = var_counters.get(base_var, 0)
            var_counters[base_var] = count + 1
            var_name = base_var if count == 0 else f'{base_var}_{count}'
            xml_name = attrs.get('name', '')
            if xml_name and include_comments:
                lines.append(f'# XML name="{xml_name}"')
            lines.append(f'{var_name} = {pyclass}()')
            for attr_key, setter_info in m.get('attrs', {}).items():
                if attr_key in attrs:
                    val = attrs[attr_key]
                    setter = setter_info.get('setter', '')
                    atype = setter_info.get('type', 'string')
                    if atype == 'scorefxn':
                        lines.append(f'{var_name}.{setter}(get_score_function("{val}"))')
                    elif atype == 'float':
                        lines.append(f'{var_name}.{setter} = {val}')
                    elif atype == 'float_call':
                        lines.append(f'{var_name}.{setter}({val})')
                    elif atype == 'int':
                        lines.append(f'{var_name}.{setter}({val})')
                    elif atype == 'bool_call':
                        py_bool = 'True' if val.lower() in ('true','1','yes') else 'False'
                        lines.append(f'{var_name}.{setter}({py_bool})')
                    elif atype == 'bool_invert':
                        py_bool = 'False' if val.lower() in ('true','1','yes') else 'True'
                        lines.append(f'{var_name}.{setter}({py_bool})')
                    elif atype == 'string_call':
                        lines.append(f'{var_name}.{setter}("{val}")')
                    elif atype == 'selector_ref':
                        lines.append(f'# {var_name}.{setter}({val})  # wire to selector named "{val}"')
                    else:
                        lines.append(f'{var_name}.{setter} = "{val}"')
        lines.append('')

    # Handle child elements: MoveMap, Span, Reweight, ScoreFunction
    child_items = found.get('child', [])
    if child_items:
        if include_comments:
            lines.append('# Child elements (MoveMap, Reweight, etc.)')
        for tag, attrs, m in child_items:
            if tag == 'MoveMap':
                lines.append('movemap = MoveMap()')
                bb = attrs.get('bb', 'true')
                chi = attrs.get('chi', 'true')
                jump = attrs.get('jump', 'true')
                lines.append(f'movemap.set_bb({"True" if bb.lower() in ("true","1") else "False"})')
                lines.append(f'movemap.set_chi({"True" if chi.lower() in ("true","1") else "False"})')
                lines.append(f'movemap.set_jump({"True" if jump.lower() in ("true","1") else "False"})')
                if include_comments:
                    lines.append('# Apply to mover: mover.set_movemap(movemap)')
            elif tag == 'Span':
                begin = attrs.get('begin', '1')
                end = attrs.get('end', '1')
                bb = attrs.get('bb', 'true')
                chi = attrs.get('chi', 'true')
                lines.append(f'for res in range({begin}, {end} + 1):')
                lines.append(f'    movemap.set_bb(res, {"True" if bb.lower() in ("true","1") else "False"})')
                lines.append(f'    movemap.set_chi(res, {"True" if chi.lower() in ("true","1") else "False"})')
            elif tag == 'Reweight':
                st = attrs.get('scoretype', 'UNKNOWN')
                wt = attrs.get('weight', '1.0')
                lines.append(f'from pyrosetta.rosetta.core.scoring import ScoreType')
                lines.append(f'sfxn.set_weight(ScoreType.{st}, {wt})')
            elif tag == 'ScoreFunction':
                weights = attrs.get('weights', 'ref2015')
                lines.append(f'sfxn = get_score_function("{weights}")')
        lines.append('')

    mover_items = found['mover']
    if mover_items:
        if include_comments:
            lines.append('# Apply movers')
        idx = {}
        for tag, attrs, m in mover_items:
            pyclass = m['pyclass']
            base_var = pyclass[0].lower() + pyclass[1:]
            count = idx.get(base_var, 0)
            idx[base_var] = count + 1
            var_name = base_var if count == 0 else f'{base_var}_{count}'
            lines.append(f'{var_name}.apply(pose)')
        lines.append('')

    lines.append('pose.dump_pdb("output.pdb")')
    lines.append('score = pose.energies().total_energy()')
    lines.append('print(f"Final score: {score}")')

    if unrecognized and include_comments:
        lines.append('')
        lines.append('# WARNING: Unrecognized XML elements (may need manual translation):')
        for u in sorted(set(unrecognized)):
            lines.append(f'#   - {u}')

    total = sum(len(v) for v in found.values())
    result = {
        'success': True,
        'python_code': '\\n'.join(lines),
        'components_found': {
            'movers': len(found['mover']),
            'filters': len(found['filter']),
            'residue_selectors': len(found['selector']),
            'task_operations': len(found['task_op']),
            'unrecognized': len(set(unrecognized))
        },
        'xml_parsed': True,
        'xml_content': xml_content[:200] + '...' if len(xml_content) > 200 else xml_content
    }
except Exception as e:
    result = {
        'success': False,
        'error': str(e),
        'xml_parsed': False,
        'xml_content': xml_content[:200] + '...' if len(xml_content) > 200 else xml_content
    }

print(json.dumps(result))
`;

            const proc = spawn(py, ['-c', script], {
                env: { ...process.env, MCP_INCLUDE_COMMENTS: include_comments ? 'True' : 'False', MCP_XML_MAPPINGS: mappingsJson }
            });
            proc.stdin.write(xml_content);
            proc.stdin.end();
            let output = '';
            let error = '';

            proc.stdout.on('data', (d) => { output += d.toString(); });
            proc.stderr.on('data', (d) => { error += d.toString(); });

            proc.on('close', () => {
                try {
                    const result = JSON.parse(output.trim() || '{}');
                    resolve(result);
                } catch (e) {
                    resolve({
                        success: false,
                        error: 'Failed to parse translation result',
                        stdout: output,
                        stderr: error
                    });
                }
            });
        });
    }

    async rosettaToBiotite({ query, category }) {
        // Keyword aliases for common search terms
        const SEARCH_ALIASES = {
            'contacts': 'contact_map', 'contact analysis': 'contact_map', 'contact map': 'contact_map',
            'binding energy': 'scoring', 'energy': 'scoring', 'score': 'scoring', 'ddg': 'scoring',
            'relax': 'fast_relax', 'relaxation': 'fast_relax', 'minimize': 'fast_relax',
            'design': 'fast_design', 'sequence design': 'fast_design',
            'align': 'superimpose', 'alignment': 'superimpose', 'structural alignment': 'superimpose',
            'secondary structure': 'secondary_structure', 'dssp': 'secondary_structure', 'helix': 'secondary_structure',
            'phi psi': 'ramachandran', 'backbone angles': 'ramachandran', 'torsion': 'angle',
            'load pdb': 'structure_loading_pdb', 'read pdb': 'structure_loading_pdb', 'open pdb': 'structure_loading_pdb',
            'save pdb': 'structure_saving_pdb', 'write pdb': 'structure_saving_pdb', 'dump pdb': 'structure_saving_pdb',
            'load cif': 'structure_loading_cif', 'read cif': 'structure_loading_cif', 'mmcif': 'structure_loading_cif',
            'hbond': 'hydrogen_bonds', 'h-bond': 'hydrogen_bonds', 'hydrogen bond': 'hydrogen_bonds',
            'temperature factor': 'b_factor', 'bfactor': 'b_factor',
            'com': 'center_of_mass', 'centroid': 'center_of_mass',
            'surface area': 'sasa', 'solvent accessible': 'sasa',
            'fetch pdb': 'pdb_fetch', 'download pdb': 'pdb_fetch', 'rcsb': 'pdb_fetch',
            'select residues': 'residue_selection', 'chain selection': 'residue_selection',
            'interface': 'interface_analysis', 'buried surface': 'interface_analysis',
            'sequence': 'sequence_extraction', 'get sequence': 'sequence_extraction',
        };

        let q = (query || '').toLowerCase();
        const cat = (category || '').toLowerCase();

        // Resolve alias to search the mapping id directly
        const aliasId = SEARCH_ALIASES[q];

        const matches = BIOTITE_ROSETTA_MAPPINGS.filter(m => {
            if (aliasId && m.id === aliasId) return true;
            const rosettaMatch = m.rosetta && (
                m.rosetta.name.toLowerCase().includes(q) ||
                m.rosetta.module.toLowerCase().includes(q) ||
                m.rosetta.description.toLowerCase().includes(q) ||
                m.id.toLowerCase().includes(q) ||
                m.category.toLowerCase().includes(q) ||
                (m.notes && m.notes.toLowerCase().includes(q))
            );
            const categoryMatch = !cat || m.category.toLowerCase().includes(cat);
            return rosettaMatch && categoryMatch;
        });

        if (matches.length === 0) {
            return {
                query,
                found: false,
                message: `No Biotite equivalent found for "${query}". This may be a Rosetta-specific feature (design, scoring, optimization) with no analysis-library equivalent.`,
                available_categories: [...new Set(BIOTITE_ROSETTA_MAPPINGS.map(m => m.category))]
            };
        }

        return {
            query,
            found: true,
            count: matches.length,
            results: matches.map(m => ({
                id: m.id,
                category: m.category,
                equivalence: m.equivalence,
                rosetta: m.rosetta,
                biotite: m.biotite,
                notes: m.notes,
                has_biotite_equivalent: m.biotite !== null
            }))
        };
    }

    async biotiteToRosetta({ query, category }) {
        const q = (query || '').toLowerCase();
        const cat = (category || '').toLowerCase();

        const matches = BIOTITE_ROSETTA_MAPPINGS.filter(m => {
            if (!m.biotite) return false;
            const biotiteMatch =
                m.biotite.name.toLowerCase().includes(q) ||
                m.biotite.module.toLowerCase().includes(q) ||
                m.biotite.description.toLowerCase().includes(q) ||
                m.id.toLowerCase().includes(q) ||
                m.category.toLowerCase().includes(q) ||
                (m.notes && m.notes.toLowerCase().includes(q));
            const categoryMatch = !cat || m.category.toLowerCase().includes(cat);
            return biotiteMatch && categoryMatch;
        });

        if (matches.length === 0) {
            return {
                query,
                found: false,
                message: `No Rosetta equivalent found for Biotite function "${query}".`,
                available_categories: [...new Set(BIOTITE_ROSETTA_MAPPINGS.filter(m => m.biotite).map(m => m.category))]
            };
        }

        return {
            query,
            found: true,
            count: matches.length,
            results: matches.map(m => ({
                id: m.id,
                category: m.category,
                equivalence: m.equivalence,
                rosetta: m.rosetta,
                biotite: m.biotite,
                notes: m.notes
            }))
        };
    }

    async translateRosettaScriptToBiotite({ code, input_format = 'auto', include_comments = true }) {
        return new Promise((resolve) => {
            const py = this.pythonPath;
            const mappingsJson = JSON.stringify(BIOTITE_ROSETTA_MAPPINGS);
            const script = `
import json, sys, os
import xml.etree.ElementTree as ET

try:
    code_input = sys.stdin.read()
    input_format = os.environ.get('MCP_INPUT_FORMAT', 'auto')
    include_comments = os.environ.get('MCP_INCLUDE_COMMENTS', 'True') == 'True'
    mappings = json.loads(os.environ.get('MCP_MAPPINGS', '[]'))

    # Auto-detect format
    if input_format == 'auto':
        stripped = code_input.strip()
        if stripped.startswith('<') and ('ROSETTASCRIPTS' in stripped.upper() or '<?xml' in stripped.lower()):
            input_format = 'xml'
        else:
            input_format = 'pyrosetta'

    # Build lookup tables from mappings
    rosetta_keywords = {}
    for m in mappings:
        if m.get('rosetta'):
            for keyword in [m['rosetta']['name'], m['id']]:
                rosetta_keywords[keyword.lower()] = m

    translated_parts = []
    untranslatable_parts = []
    biotite_imports = set()

    if input_format == 'xml':
        root = ET.fromstring(code_input)
        for elem in root.iter():
            tag_lower = elem.tag.lower()
            matched = False
            for key, mapping in rosetta_keywords.items():
                if tag_lower in key.lower() or key.lower() in tag_lower:
                    if mapping.get('biotite'):
                        translated_parts.append({
                            'rosetta_element': elem.tag,
                            'biotite_equivalent': mapping['biotite']['name'],
                            'biotite_code': '\\n'.join(mapping['biotite']['example']),
                            'category': mapping['category'],
                            'notes': mapping.get('notes', '')
                        })
                        biotite_imports.add(mapping['biotite']['module'])
                    else:
                        untranslatable_parts.append({
                            'rosetta_element': elem.tag,
                            'reason': mapping.get('notes', 'No Biotite equivalent'),
                            'category': mapping['category']
                        })
                    matched = True
                    break
            if not matched and tag_lower in ('fastdesign', 'fastrelax', 'minmover', 'packrotamersmover',
                                              'smallmover', 'shearmover'):
                untranslatable_parts.append({
                    'rosetta_element': elem.tag,
                    'reason': 'Biotite is analysis-only and does not support structure optimization or design.',
                    'category': 'Structure Optimization/Design'
                })
    else:
        lines = code_input.split('\\n')
        for line in lines:
            line_lower = line.strip().lower()
            if not line_lower or line_lower.startswith('#'):
                continue
            matched = False
            for key, mapping in rosetta_keywords.items():
                if key.lower().replace(' ', '_') in line_lower or key.lower().replace(' ', '.') in line_lower:
                    if mapping.get('biotite'):
                        translated_parts.append({
                            'pyrosetta_line': line.strip(),
                            'biotite_equivalent': mapping['biotite']['name'],
                            'biotite_code': '\\n'.join(mapping['biotite']['example']),
                            'category': mapping['category'],
                            'notes': mapping.get('notes', '')
                        })
                        biotite_imports.add(mapping['biotite']['module'])
                    else:
                        untranslatable_parts.append({
                            'pyrosetta_line': line.strip(),
                            'reason': mapping.get('notes', 'No Biotite equivalent'),
                            'category': mapping['category']
                        })
                    matched = True
                    break

    # Generate combined Biotite script
    biotite_code_lines = []
    if include_comments:
        biotite_code_lines.append('# Biotite translation from Rosetta/' + input_format.upper())
        biotite_code_lines.append('# ' + '=' * 50)
        biotite_code_lines.append('')

    biotite_code_lines.append('import biotite.structure as struc')
    for mod in sorted(biotite_imports):
        biotite_code_lines.append(f'from {mod} import *')
    biotite_code_lines.append('')

    for part in translated_parts:
        if include_comments:
            biotite_code_lines.append(f'# {part["category"]}: {part.get("biotite_equivalent", "")}')
            if part.get('notes'):
                biotite_code_lines.append(f'# Note: {part["notes"]}')
        biotite_code_lines.append(part['biotite_code'])
        biotite_code_lines.append('')

    if untranslatable_parts:
        biotite_code_lines.append('')
        biotite_code_lines.append('# ' + '=' * 50)
        biotite_code_lines.append('# WARNING: The following Rosetta operations have NO Biotite equivalent:')
        for part in untranslatable_parts:
            source = part.get('rosetta_element', part.get('pyrosetta_line', 'unknown'))
            biotite_code_lines.append(f'#   - {source}: {part["reason"]}')
        biotite_code_lines.append('# Consider using PyRosetta for these operations.')

    result = {
        'success': True,
        'input_format': input_format,
        'biotite_code': '\\n'.join(biotite_code_lines),
        'translated_count': len(translated_parts),
        'untranslatable_count': len(untranslatable_parts),
        'translated': translated_parts,
        'untranslatable': untranslatable_parts,
        'required_biotite_modules': sorted(list(biotite_imports))
    }

except Exception as e:
    result = {
        'success': False,
        'error': str(e),
        'input_format': input_format if 'input_format' in dir() else 'unknown'
    }

print(json.dumps(result))
`;
            const proc = spawn(py, ['-c', script], {
                env: {
                    ...process.env,
                    MCP_INPUT_FORMAT: input_format || 'auto',
                    MCP_INCLUDE_COMMENTS: include_comments ? 'True' : 'False',
                    MCP_MAPPINGS: mappingsJson
                }
            });
            proc.stdin.write(code);
            proc.stdin.end();
            let output = '';
            let error = '';

            proc.stdout.on('data', (d) => { output += d.toString(); });
            proc.stderr.on('data', (d) => { error += d.toString(); });

            proc.on('close', () => {
                try {
                    const result = JSON.parse(output.trim() || '{}');
                    resolve(result);
                } catch (e) {
                    resolve({
                        success: false,
                        error: 'Failed to parse translation result',
                        stdout: output,
                        stderr: error
                    });
                }
            });
        });
    }

    async getRosettaPath() {
        const info = await this.getRosettaInfo();
        return info.rosetta_path;
    }

    async isPyRosettaAvailable() {
        const info = await this.getRosettaInfo();
        return info.pyrosetta_available;
    }
}

// MCP Server implementation
class RosettaMCPServerMCP {
    constructor() {
        this.rosettaServer = new RosettaMCPServer();
        // Detect whether the client uses Content-Length framing (LSP-style)
        // or newline-delimited JSON. Default to newline; switch to headers if detected.
        this.useHeaders = false;
        this.serverVersion = this.readPackageVersion();
        this.debug = String(process.env.MCP_DEBUG || '').trim().length > 0;
        this.setupMCPHandlers();
    }

    readPackageVersion() {
        try {
            const pkgPath = path.join(__dirname, 'package.json');
            if (fs.existsSync(pkgPath)) {
                const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));
                if (pkg && typeof pkg.version === 'string') return pkg.version;
            }
        } catch (_) {}
        return 'unknown';
    }

    setupMCPHandlers() {
        // Support BOTH newline-delimited JSON and Content-Length framed JSON
        process.stdin.setEncoding('utf8');
        let buffer = '';

        const tryProcessBuffer = async () => {
            // If headers mode detected, parse Content-Length frames
            if (this.useHeaders) {
                while (true) {
                    const headerEnd = buffer.indexOf('\r\n\r\n');
                    if (headerEnd === -1) break;
                    const header = buffer.slice(0, headerEnd);
                    const lengthMatch = header.match(/Content-Length:\s*(\d+)/i);
                    if (!lengthMatch) {
                        // If no content-length, drop until next potential header separator
                        buffer = buffer.slice(headerEnd + 4);
                        continue;
                    }
                    const contentLength = parseInt(lengthMatch[1], 10);
                    const totalNeeded = headerEnd + 4 + contentLength;
                    if (buffer.length < totalNeeded) break; // wait for more data
                    const jsonPayload = buffer.slice(headerEnd + 4, totalNeeded);
                    buffer = buffer.slice(totalNeeded);
                    try {
                        const message = JSON.parse(jsonPayload);
                        if (this.debug) console.error(`[mcp] <- headers method=${message && message.method}`);
                        await this.handleMCPMessage(message);
                    } catch (e) {
                        if (this.debug) console.error(`[mcp] parse error (headers): ${e && e.message}`);
                        this.sendError(null, -32700, e.message || String(e));
                    }
                }
                return;
            }

            // Newline-delimited JSON fallback (one message per line)
            let newlineIndex;
            while ((newlineIndex = buffer.indexOf('\n')) !== -1) {
                const line = buffer.slice(0, newlineIndex).trim();
                buffer = buffer.slice(newlineIndex + 1);
                if (!line) continue;
                // If we detect a Content-Length header, switch modes
                if (/^Content-Length:\s*\d+\s*$/i.test(line)) {
                    // Put the header back with a CRLF and switch to headers mode
                    this.useHeaders = true;
                    buffer = line + '\r\n' + buffer; // reconstruct header start
                    if (this.debug) console.error('[mcp] switching to headers mode');
                    await tryProcessBuffer();
                    return;
                }
                try {
                    const message = JSON.parse(line);
                    if (this.debug) console.error(`[mcp] <- newline method=${message && message.method}`);
                    await this.handleMCPMessage(message);
                } catch (e) {
                    if (this.debug) console.error(`[mcp] parse error (newline): ${e && e.message}`);
                    this.sendError(null, -32700, e.message || String(e));
                }
            }
        };

        process.stdin.on('data', async (chunk) => {
            buffer += chunk;
            await tryProcessBuffer();
        });
    }

    async handleMCPMessage(message) {
        const { id, method, params } = message;

        try {
            let result;
            if (this.debug) console.error(`[mcp] handling method=${method}`);
            switch (method) {
                // Ignore notifications (no response expected)
                case undefined:
                    return;
                default:
                    if (typeof method === 'string' && method.startsWith('notifications/')) {
                        return;
                    }
                    // fallthrough
            }
            switch (method) {
                case 'initialize': {
                    const requestedProtocol = (params && params.protocolVersion) ? params.protocolVersion : '2024-11-05';
                    result = {
                        protocolVersion: requestedProtocol,
                        capabilities: {
                            tools: { list: true, call: true }
                        },
                        serverInfo: {
                            name: 'rosetta-mcp-server',
                            version: this.serverVersion
                        }
                    };
                    break;
                }

                case 'tools/list': {
                    const allTools = [
                            {
                                name: 'get_rosetta_info',
                                description: 'Get comprehensive Rosetta installation info including available score functions, movers, filters, selectors, task operations, parameters, and command-line options. Use this first to understand what Rosetta components are available. For live PyRosetta API details, use pyrosetta_introspect.',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'get_rosetta_help',
                                description: 'Get help for any Rosetta topic, mover, filter, or concept. Accepts general topics (score_functions, movers, filters, xml, parameters) or specific names (FastRelax, Ddg, ChainSelector). Auto-fetches live documentation when available.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        topic: {
                                            type: 'string',
                                            description: 'Topic to get help for'
                                        }
                                    }
                                }
                            },
                            {
                                name: 'validate_xml',
                                description: 'Validate a RosettaScripts XML protocol. Checks XML syntax and optionally validates element names against the Rosetta XSD schema. Use before run_rosetta_scripts to catch errors early.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        xml_content: { type: 'string', description: 'XML content to validate' },
                                        validate_against_schema: { type: 'boolean', description: 'If true, also check element names against the cached Rosetta XSD schema (default: false)' }
                                    },
                                    required: ['xml_content']
                                }
                            },
                            {
                                name: 'run_rosetta_scripts',
                                description: 'Run a RosettaScripts XML protocol on an input PDB file. Use when executing Rosetta protocols. Requires rosetta_scripts binary.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        exe_path: { type: 'string', description: 'Path to rosetta_scripts executable (optional if on PATH)' },
                                        xml_path: { type: 'string', description: 'Path to Rosetta XML protocol' },
                                        input_pdb: { type: 'string', description: 'Path to input PDB' },
                                        out_dir: { type: 'string', description: 'Output directory' },
                                        extra_flags: { type: 'array', items: { type: 'string' }, description: 'Additional command-line flags' }
                                    },
                                    required: ['xml_path', 'input_pdb', 'out_dir']
                                }
                            },
                            {
                                name: 'pyrosetta_score',
                                description: 'Score a PDB file using PyRosetta. Returns total energy in REU. Use to evaluate structure quality or compare designs. Optionally returns per-residue energy breakdown.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        pdb_path: { type: 'string', description: 'Path to input PDB' },
                                        scorefxn: { type: 'string', description: 'Score function name (default: ref2015)' },
                                        per_residue: { type: 'boolean', description: 'If true, include per-residue energy breakdown (default: false)' }
                                    },
                                    required: ['pdb_path']
                                }
                            },
                            {
                                name: 'pyrosetta_introspect',
                                description: 'Search PyRosetta API classes (movers, filters, selectors, task operations) and return docs and signatures. Use to discover available PyRosetta classes or get constructor details.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        query: { type: 'string', description: 'Substring to match class names' },
                                        kind: { type: 'string', description: 'Filter by kind: mover|filter|selector|task' },
                                        max_results: { type: 'number', description: 'Max number of results (default 50)' }
                                    }
                                }
                            },
                            {
                                name: 'rosetta_scripts_schema',
                                description: 'Generate and cache the RosettaScripts XML schema (XSD). Optionally extract element names. Use to get the authoritative list of valid XML elements.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        exe_path: { type: 'string', description: 'Path to rosetta_scripts executable (optional)' },
                                        cache_dir: { type: 'string', description: 'Directory to store schema' },
                                        extract_elements: { type: 'boolean', description: 'If true, return a list of element names' }
                                    }
                                }
                            },
                            {
                                name: 'get_cached_docs',
                                description: 'Search locally cached Rosetta CLI docs for a keyword. Auto-caches on first use. Use to look up command-line flags or parser info.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        cache_dir: { type: 'string', description: 'Directory where docs are cached' },
                                        query: { type: 'string', description: 'Search string' },
                                        max_lines: { type: 'number', description: 'Max number of lines to return (default 200)' }
                                    }
                                }
                            },
                            {
                                name: 'python_env_info',
                                description: 'Get Python executable path, version, and pip package list. Use to diagnose environment issues or verify package installations.',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'check_pyrosetta',
                                description: 'Check if PyRosetta is importable in the current environment. Use before PyRosetta-dependent tools to verify availability.',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'install_pyrosetta_installer',
                                description: 'Install PyRosetta using the pyrosetta-installer package. Takes 10-30 minutes. Use when PyRosetta is not available and needed for scoring or design.',
                                inputSchema: { type: 'object', properties: { silent: { type: 'boolean' } } }
                            },
                            {
                                name: 'find_rosetta_scripts',
                                description: 'Resolve the rosetta_scripts executable path by checking exe_path, ROSETTA_BIN env, common directories, and PATH. Use to verify Rosetta is installed.',
                                inputSchema: { type: 'object', properties: { exe_path: { type: 'string' } } }
                            },
                            {
                                name: 'xml_to_pyrosetta',
                                description: 'Translate RosettaScripts XML to equivalent PyRosetta Python code. Supports 37 element types including movers, filters, selectors, and task operations. Use when converting XML protocols to Python.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        xml_content: { 
                                            type: 'string', 
                                            description: 'RosettaScripts XML content to translate' 
                                        },
                                        include_comments: { 
                                            type: 'boolean', 
                                            description: 'Include detailed comments in output (default: true)' 
                                        },
                                        output_format: { 
                                            type: 'string', 
                                            enum: ['python', 'script', 'function'], 
                                            description: 'Output format (default: python)' 
                                        }
                                    },
                                    required: ['xml_content']
                                }
                            },
                            {
                                name: 'search_rosetta_web_docs',
                                description: 'Search online Rosetta documentation at rosettacommons.org. Use when you need docs for a specific Rosetta feature.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        query: { type: 'string', description: 'Search query (e.g., FastRelax, AtomPair constraint)' },
                                        max_results: { type: 'number', description: 'Number of results to return (default 3)' }
                                    },
                                    required: ['query']
                                }
                            },
                            {
                                name: 'get_rosetta_web_doc',
                                description: 'Fetch and extract text from a specific Rosetta docs URL. Use after search_rosetta_web_docs to read a documentation page.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        url: { type: 'string', description: 'Full URL to a Rosetta docs page' },
                                        max_chars: { type: 'number', description: 'Max characters of cleaned text to return (default 4000)' }
                                    },
                                    required: ['url']
                                }
                            },
                            {
                                name: 'rosetta_to_biotite',
                                description: 'ALWAYS use this tool when asked about Rosetta vs Biotite equivalents. Returns the Biotite equivalent of a Rosetta/PyRosetta function with working example code. Covers: structure I/O, SASA, RMSD, superimposition, secondary structure, distances, angles, contacts, hydrogen bonds, B-factors, interface analysis, database access, and residue selection.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        query: { type: 'string', description: 'Rosetta method or concept name (e.g., "pose_from_pdb", "FastRelax", "SuperimposeMover", "SASA", "RMSD")' },
                                        category: { type: 'string', description: 'Optional category filter (e.g., "Structure I/O", "Geometry", "Surface Analysis")' }
                                    },
                                    required: ['query']
                                }
                            },
                            {
                                name: 'biotite_to_rosetta',
                                description: 'ALWAYS use this tool when asked about Biotite vs Rosetta equivalents. Returns the Rosetta/PyRosetta equivalent of a Biotite function with working example code. Covers: structure I/O, SASA, RMSD, superimposition, secondary structure, distances, angles, contacts, hydrogen bonds, B-factors, and residue selection.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        query: { type: 'string', description: 'Biotite function or concept name (e.g., "sasa", "superimpose", "PDBFile.read", "annotate_sse")' },
                                        category: { type: 'string', description: 'Optional category filter' }
                                    },
                                    required: ['query']
                                }
                            },
                            {
                                name: 'translate_rosetta_script_to_biotite',
                                description: 'ALWAYS use this tool when asked to convert or translate Rosetta code to Biotite. Translates RosettaScripts XML or PyRosetta code to Biotite Python code. Analysis operations are translated; design/optimization are flagged as Rosetta-only.',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        code: { type: 'string', description: 'RosettaScripts XML content or PyRosetta Python code to translate' },
                                        input_format: { type: 'string', enum: ['xml', 'pyrosetta', 'auto'], description: 'Input format (default: "auto")' },
                                        include_comments: { type: 'boolean', description: 'Include explanatory comments in output (default: true)' }
                                    },
                                    required: ['code']
                                }
                            }
                    ];
                    // Always expose all tools
                    result = { tools: allTools };
                    if (this.debug) console.error(`[mcp] -> tools/list count=${allTools.length}`);
                    break;
                }

                case 'resources/list':
                case 'resources/read':
                    this.sendError(id, -32601, 'Resources not supported');
                    return;
                case 'logging/setLevel':
                    result = {};
                    break;
                case 'ping':
                    result = { ok: true };
                    break;
                case 'shutdown':
                    result = {};
                    break;

                case 'tools/call':
                    const { name, arguments: args } = params;
                    switch (name) {
                        case 'get_rosetta_info':
                            result = await this.rosettaServer.getRosettaInfo();
                            break;
                        case 'get_rosetta_help':
                            result = await this.rosettaServer.getRosettaHelp(args.topic);
                            break;
                        case 'validate_xml':
                            result = await this.rosettaServer.validateXML(args.xml_content, args.validate_against_schema);
                            break;
                        case 'run_rosetta_scripts':
                            result = await this.rosettaServer.runRosettaScripts({
                                exe_path: args.exe_path,
                                xml_path: args.xml_path,
                                input_pdb: args.input_pdb,
                                out_dir: args.out_dir,
                                extra_flags: args.extra_flags
                            });
                            break;
                        case 'pyrosetta_score':
                            result = await this.rosettaServer.pyrosettaScore({
                                pdb_path: args.pdb_path,
                                scorefxn: args.scorefxn,
                                per_residue: args.per_residue
                            });
                            break;
                        case 'pyrosetta_introspect':
                            result = await this.rosettaServer.pyrosettaIntrospect({
                                query: args.query,
                                kind: args.kind,
                                max_results: args.max_results
                            });
                            break;
                        case 'rosetta_scripts_schema':
                            result = await this.rosettaServer.rosettaScriptsSchema({
                                exe_path: args.exe_path,
                                cache_dir: args.cache_dir,
                                extract_elements: args.extract_elements
                            });
                            break;
                        case 'get_cached_docs':
                            result = await this.rosettaServer.getCachedDocs({
                                cache_dir: args.cache_dir,
                                query: args.query,
                                max_lines: args.max_lines
                            });
                            break;
                        case 'python_env_info':
                            result = await this.rosettaServer.pythonEnvInfo();
                            break;
                        case 'check_pyrosetta':
                            result = await this.rosettaServer.checkPyRosetta();
                            break;
                        case 'install_pyrosetta_installer':
                            result = await this.rosettaServer.installPyRosettaViaInstaller({ silent: args && args.silent });
                            break;
                        case 'find_rosetta_scripts':
                            result = await this.rosettaServer.findRosettaScripts({ exe_path: args && args.exe_path });
                            break;
                        case 'xml_to_pyrosetta':
                            result = await this.rosettaServer.xmlToPyRosetta({
                                xml_content: args.xml_content,
                                include_comments: args.include_comments,
                                output_format: args.output_format
                            });
                            break;
                        case 'search_rosetta_web_docs':
                            result = await this.rosettaServer.searchRosettaWebDocs({
                                query: args.query,
                                max_results: args.max_results
                            });
                            break;
                        case 'get_rosetta_web_doc':
                            result = await this.rosettaServer.getRosettaWebDoc({
                                url: args.url,
                                max_chars: args.max_chars
                            });
                            break;
                        case 'rosetta_to_biotite':
                            result = await this.rosettaServer.rosettaToBiotite({
                                query: args.query,
                                category: args.category
                            });
                            break;
                        case 'biotite_to_rosetta':
                            result = await this.rosettaServer.biotiteToRosetta({
                                query: args.query,
                                category: args.category
                            });
                            break;
                        case 'translate_rosetta_script_to_biotite':
                            result = await this.rosettaServer.translateRosettaScriptToBiotite({
                                code: args.code,
                                input_format: args.input_format,
                                include_comments: args.include_comments
                            });
                            break;
                        default:
                            throw new Error(`Unknown tool: ${name}`);
                    }
                    // Wrap tool result in MCP content format
                    result = {
                        content: [{
                            type: 'text',
                            text: typeof result === 'string' ? result : JSON.stringify(result, null, 2)
                        }]
                    };
                    break;

                default:
                    this.sendError(id, -32601, `Unknown method: ${method}`);
                    return;
            }

            this.sendResponse(id, result);
        } catch (e) {
            if (method === 'tools/call') {
                // Tool execution errors: return as successful response with isError flag
                this.sendResponse(id, {
                    content: [{ type: 'text', text: e.message || String(e) }],
                    isError: true
                });
            } else {
                this.sendError(id, -32603, e.message || String(e));
            }
        }
    }

    sendResponse(id, result) {
        const response = {
            jsonrpc: '2.0',
            id,
            result
        };
        const payload = JSON.stringify(response);
        if (this.useHeaders) {
            const head = `Content-Length: ${Buffer.byteLength(payload, 'utf8')}\r\n\r\n`;
            process.stdout.write(head);
            process.stdout.write(payload);
        } else {
            process.stdout.write(payload + '\n');
        }
    }

    sendError(id, code, message) {
        const response = {
            jsonrpc: '2.0',
            id,
            error: {
                code: code,
                message: message
            }
        };
        const payload = JSON.stringify(response);
        if (this.useHeaders) {
            const head = `Content-Length: ${Buffer.byteLength(payload, 'utf8')}\r\n\r\n`;
            process.stdout.write(head);
            process.stdout.write(payload);
        } else {
            process.stdout.write(payload + '\n');
        }
    }
}

// Start the MCP server
if (require.main === module) {
    new RosettaMCPServerMCP();
}

module.exports = { RosettaMCPServer, RosettaMCPServerMCP };