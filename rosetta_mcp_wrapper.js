#!/usr/bin/env node
/**
 * Rosetta MCP Server Wrapper
 * Makes the Python Rosetta server compatible with MCP protocol
 */

const { spawn } = require('child_process');
const path = require('path');

class RosettaMCPServer {
    constructor() {
        this.pythonPath = 'python3';
        this.serverPath = path.join(__dirname, 'rosetta_mcp_server.py');
    }

    resolveRosettaScriptsPath(preferredPath) {
        if (preferredPath && preferredPath.length > 0) {
            return preferredPath;
        }

        if (process.env.ROSETTA_BIN && process.env.ROSETTA_BIN.includes('rosetta_scripts')) {
            return process.env.ROSETTA_BIN;
        }

        const candidateDirs = [
            // Typical source build locations
            path.join(process.env.HOME || '', 'rosetta', 'main', 'source', 'bin'),
            '/opt/rosetta/main/source/bin',
            '/usr/local/rosetta/main/source/bin',
            '/Users/arielben-sasson/dev/open_repos/rosetta/main/source/bin'
        ].filter(Boolean);

        const fs = require('fs');
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
            const process = spawn(this.pythonPath, [this.serverPath, 'info']);
            let output = '';
            let error = '';

            process.stdout.on('data', (data) => {
                output += data.toString();
            });

            process.stderr.on('data', (data) => {
                error += data.toString();
            });

            process.on('close', (code) => {
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
        return new Promise((resolve, reject) => {
            const args = [this.serverPath, 'help'];
            if (topic) args.push(topic);

            const process = spawn(this.pythonPath, args);
            let output = '';
            let error = '';

            process.stdout.on('data', (data) => {
                output += data.toString();
            });

            process.stderr.on('data', (data) => {
                error += data.toString();
            });

            process.on('close', (code) => {
                if (code === 0) {
                    resolve(output.trim());
                } else {
                    reject(new Error(`Server exited with code ${code}: ${error}`));
                }
            });
        });
    }

    async validateXML(xmlContent) {
        return new Promise((resolve, reject) => {
            // Write XML content to temporary file
            const fs = require('fs');
            const tmpFile = path.join(__dirname, 'temp_validation.xml');
            
            try {
                fs.writeFileSync(tmpFile, xmlContent);
                
                const process = spawn(this.pythonPath, [this.serverPath, 'validate', tmpFile]);
                let output = '';
                let error = '';

                process.stdout.on('data', (data) => {
                    output += data.toString();
                });

                process.stderr.on('data', (data) => {
                    error += data.toString();
                });

                process.on('close', (code) => {
                    // Clean up temp file
                    try {
                        fs.unlinkSync(tmpFile);
                    } catch (e) {
                        // Ignore cleanup errors
                    }

                    if (code === 0) {
                        try {
                            const result = JSON.parse(output);
                            resolve(result);
                        } catch (e) {
                            reject(new Error('Failed to parse validation result'));
                        }
                    } else {
                        reject(new Error(`Validation failed with code ${code}: ${error}`));
                    }
                });
            } catch (e) {
                reject(new Error(`Failed to create temporary file: ${e.message}`));
            }
        });
    }

    async runRosettaScripts({ exe_path, xml_path, input_pdb, out_dir, extra_flags }) {
        return new Promise((resolve, reject) => {
            if (!xml_path || !input_pdb || !out_dir) {
                reject(new Error('xml_path, input_pdb, and out_dir are required'));
                return;
            }

            const rosettaExe = this.resolveRosettaScriptsPath(exe_path);
            const fs = require('fs');
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

    async pyrosettaScore({ pdb_path, scorefxn }) {
        return new Promise((resolve, reject) => {
            if (!pdb_path) {
                reject(new Error('pdb_path is required'));
                return;
            }

            const py = this.pythonPath;
            const script = `
import json
try:
    import pyrosetta
except Exception as e:
    print(json.dumps({"error": f"PyRosetta not available: {str(e)}"}))
    raise SystemExit(0)

pyrosetta.init('-mute all')
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import get_score_function
pose = pose_from_pdb(r'''${pdb_path.replace(/'/g, "'\''")}''')
sfxn = get_score_function()
score = sfxn(pose)
print(json.dumps({"score": float(score)}))
`;

            const proc = spawn(py, ['-c', script]);
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

    async listAvailableFunctions() {
        const info = await this.getRosettaInfo();
        return {
            score_functions: info.score_functions,
            movers: info.common_movers,
            filters: info.common_filters,
            residue_selectors: info.residue_selectors,
            task_operations: info.task_operations,
            parameters: info.common_parameters,
            command_line_options: info.command_line_options
        };
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
        this.setupMCPHandlers();
    }

    setupMCPHandlers() {
        // Handle MCP protocol messages
        process.stdin.setEncoding('utf8');
        process.stdin.on('data', async (data) => {
            try {
                const message = JSON.parse(data);
                await this.handleMCPMessage(message);
            } catch (e) {
                this.sendError('Failed to parse message', e.message);
            }
        });
    }

    async handleMCPMessage(message) {
        const { id, method, params } = message;

        try {
            let result;
            switch (method) {
                case 'initialize':
                    result = {
                        protocolVersion: '2024-11-05',
                        capabilities: {
                            tools: {},
                            resources: {}
                        },
                        serverInfo: {
                            name: 'rosetta-mcp-server',
                            version: '1.0.0'
                        }
                    };
                    break;

                case 'tools/list':
                    result = {
                        tools: [
                            {
                                name: 'get_rosetta_info',
                                description: 'Get comprehensive Rosetta information',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'get_rosetta_help',
                                description: 'Get help for specific Rosetta topics',
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
                                description: 'Validate XML protocol file',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        xml_content: {
                                            type: 'string',
                                            description: 'XML content to validate'
                                        }
                                    },
                                    required: ['xml_content']
                                }
                            },
                            {
                                name: 'list_functions',
                                description: 'List available Rosetta functions',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'run_rosetta_scripts',
                                description: 'Run RosettaScripts with a given XML and input PDB',
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
                                description: 'Score a PDB using PyRosetta (requires PyRosetta installed)',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        pdb_path: { type: 'string', description: 'Path to input PDB' },
                                        scorefxn: { type: 'string', description: 'Score function name (optional)' }
                                    },
                                    required: ['pdb_path']
                                }
                            }
                        ]
                    };
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
                            result = await this.rosettaServer.validateXML(args.xml_content);
                            break;
                        case 'list_functions':
                            result = await this.rosettaServer.listAvailableFunctions();
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
                                scorefxn: args.scorefxn
                            });
                            break;
                        default:
                            throw new Error(`Unknown tool: ${name}`);
                    }
                    break;

                default:
                    throw new Error(`Unknown method: ${method}`);
            }

            this.sendResponse(id, result);
        } catch (e) {
            this.sendError(id, e.message);
        }
    }

    sendResponse(id, result) {
        const response = {
            jsonrpc: '2.0',
            id,
            result
        };
        process.stdout.write(JSON.stringify(response) + '\n');
    }

    sendError(id, error) {
        const response = {
            jsonrpc: '2.0',
            id,
            error: {
                code: -1,
                message: error
            }
        };
        process.stdout.write(JSON.stringify(response) + '\n');
    }
}

// Start the MCP server
if (require.main === module) {
    new RosettaMCPServerMCP();
}

module.exports = { RosettaMCPServer, RosettaMCPServerMCP };
