#!/usr/bin/env node
/**
 * Rosetta MCP Server Wrapper
 * Makes the Python Rosetta server compatible with MCP protocol
 */

const { spawn } = require('child_process');
const path = require('path');
const readline = require('readline');

class RosettaMCPServer {
    constructor() {
        this.pythonPath = process.env.PYTHON_BIN || 'python3';
        this.serverPath = path.join(__dirname, 'rosetta_mcp_server.py');
    }

    async pythonEnvInfo() {
        return new Promise((resolve) => {
            const py = this.pythonPath;
            const { spawn } = require('child_process');
            
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
            const { spawn } = require('child_process');
            
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
            const { spawn } = require('child_process');
            
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

    async searchPyRosettaWheels({ directory }) {
        const fs = require('fs');
        const path = require('path');
        const results = [];
        try {
            const files = fs.readdirSync(directory || '.');
            for (const f of files) {
                if (/pyrosetta.*\.whl$/i.test(f) || /PyRosetta.*\.whl$/i.test(f)) {
                    results.push(path.join(directory, f));
                }
            }
        } catch (e) {
            return { error: e.message };
        }
        return { matches: results };
    }

    resolveRosettaScriptsPath(preferredPath) {
        const fs = require('fs');

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
            '/usr/local/rosetta/main/source/bin',
            '/Users/arielben-sasson/dev/open_repos/rosetta/main/source/bin'
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
        return new Promise(async (resolve, reject) => {
            if (!pdb_path) {
                reject(new Error('pdb_path is required'));
                return;
            }

            // Check if PyRosetta is available
            const pyrosettaStatus = await this.checkPyRosetta();
            if (!pyrosettaStatus.available) {
                // Auto-install PyRosetta if not available
                console.error('\x1b[33m⚠️  PyRosetta not found. Auto-installing...\x1b[0m');
                console.error('\x1b[33m   This will take 10-30 minutes. Please be patient.\x1b[0m\n');
                
                try {
                    const installResult = await this.installPyRosettaViaInstaller({ silent: false });
                    if (!installResult.ok) {
                        reject(new Error(`PyRosetta installation failed: ${installResult.error || 'Unknown error'}`));
                        return;
                    }
                    console.error('\x1b[32m✅ PyRosetta installed successfully! Retrying score operation...\x1b[0m\n');
                } catch (installError) {
                    reject(new Error(`PyRosetta installation failed: ${installError.message}`));
                    return;
                }
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

    async pyrosettaIntrospect({ query, kind, max_results }) {
        return new Promise(async (resolve) => {
            const py = this.pythonPath;
            
            // Check if PyRosetta is available
            const pyrosettaStatus = await this.checkPyRosetta();
            if (!pyrosettaStatus.available) {
                // Auto-install PyRosetta if not available
                console.error('\x1b[33m⚠️  PyRosetta not found. Auto-installing...\x1b[0m');
                console.error('\x1b[33m   This will take 10-30 minutes. Please be patient.\x1b[0m\n');
                
                try {
                    const installResult = await this.installPyRosettaViaInstaller({ silent: false });
                    if (!installResult.ok) {
                        resolve({ error: `PyRosetta installation failed: ${installResult.error || 'Unknown error'}` });
                        return;
                    }
                    console.error('\x1b[32m✅ PyRosetta installed successfully! Retrying introspect operation...\x1b[0m\n');
                } catch (installError) {
                    resolve({ error: `PyRosetta installation failed: ${installError.message}` });
                    return;
                }
            }
            const safeQuery = (typeof query === 'string' ? query : '').replace(/`/g, '');
            const limit = Number.isInteger(max_results) && max_results > 0 ? max_results : 50;
            const script = `
import json, sys, importlib, inspect
try:
    import pyrosetta
except Exception as e:
    print(json.dumps({"error": f"PyRosetta not available: {str(e)}"}))
    raise SystemExit(0)

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

q = '${safeQuery}'.lower()
kind = '${(kind||'').toString().toLowerCase()}'

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
                    # Try to get more info if no doc
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
            if len(results) >= ${limit}:
                raise StopIteration
    except StopIteration:
        break

print(json.dumps({'results': results, 'count': len(results)}))
`;
            const proc = spawn(py, ['-c', script]);
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
            const fs = require('fs');
            const os = require('os');
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
            const fs = require('fs');
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
                    try { fs.writeFileSync(path.join(baseDir, outfile), output + (error ? `\nSTDERR:\n${error}` : '')); } catch (_) {}
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
        return new Promise((resolve, reject) => {
            const fs = require('fs');
            const baseDir = cache_dir && cache_dir.length ? cache_dir : path.join(process.cwd(), 'docs_cache');
            const files = ['help.txt', 'parser_info.txt'].map(f => path.join(baseDir, f));
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
                resolve({ count: results.length, matches: results.slice(0, limit) });
            } catch (e) {
                reject(new Error(`Failed to search cache: ${e.message}`));
            }
        });
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
        try {
            const resp = await this.fetchUrl(searchUrl);
            if (resp.status !== 200) return { error: `Search failed with status ${resp.status}` };
            const html = resp.body || '';
            // Extract anchors; prefer links pointing to rosettacommons.org
            const linkRegex = /<a[^>]+href="([^"]+)"[^>]*>([\s\S]*?)<\/a>/gi;
            const results = [];
            let m;
            while ((m = linkRegex.exec(html)) !== null && results.length < max * 3) {
                let href = m[1];
                const text = this.stripHtml(m[2]);
                // DuckDuckGo may wrap URLs with '/l/?kh=-1&uddg=<ENCODED>'
                const uddgMatch = href.match(/[?&]uddg=([^&]+)/);
                if (uddgMatch) {
                    try { href = decodeURIComponent(uddgMatch[1]); } catch (_) {}
                }
                if (/rosettacommons\.org\/.*/i.test(href)) {
                    results.push({ title: text.slice(0, 150), url: href });
                }
            }
            // De-duplicate and cap to max
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
        } catch (e) {
            return { error: `Search error: ${e.message}` };
        }
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
            const includeCommentsPython = include_comments ? 'True' : 'False';
            const script = `
import json
import xml.etree.ElementTree as ET

try:
    # Parse XML
    xml_content = '''${xml_content.replace(/'/g, "'\\''")}'''
    root = ET.fromstring(xml_content)
    
    # Extract components
    movers = []
    filters = []
    residue_selectors = []
    task_operations = []
    
    for elem in root.iter():
        tag = elem.tag
        if tag == 'FastRelax':
            movers.append(('FastRelax', elem.attrib))
        elif tag == 'MinMover':
            movers.append(('MinMover', elem.attrib))
        elif tag == 'PackRotamersMover':
            movers.append(('PackRotamersMover', elem.attrib))
        elif tag == 'ScoreType':
            filters.append(('ScoreType', elem.attrib))
        elif tag == 'Chain':
            residue_selectors.append(('Chain', elem.attrib))
        elif tag == 'RestrictToRepacking':
            task_operations.append(('RestrictToRepacking', elem.attrib))
    
    # Generate Python code
    python_code = []
    
    if ${includeCommentsPython}:
        python_code.append('# Generated PyRosetta code from RosettaScripts XML')
        python_code.append('# ==================================================')
        python_code.append('')
    
    # Import statements
    python_code.append('import pyrosetta')
    python_code.append('from pyrosetta import pose_from_pdb')
    python_code.append('from pyrosetta.rosetta.protocols.moves import *')
    python_code.append('from pyrosetta.rosetta.protocols.simple_moves import *')
    python_code.append('from pyrosetta.rosetta.protocols.filters import *')
    python_code.append('from pyrosetta.rosetta.core.select.residue_selector import *')
    python_code.append('from pyrosetta.rosetta.core.pack.task.operation import *')
    python_code.append('from pyrosetta.rosetta.core.scoring import get_score_function')
    python_code.append('')
    
    # Initialize PyRosetta
    python_code.append('# Initialize PyRosetta')
    python_code.append('pyrosetta.init("-mute all")')
    python_code.append('')
    
    # Load pose
    python_code.append('# Load your PDB file')
    python_code.append('pose = pose_from_pdb("your_protein.pdb")')
    python_code.append('')
    
    # Create movers
    if movers:
        if ${includeCommentsPython}:
            python_code.append('# Create movers')
        for mover_name, attrs in movers:
            if mover_name == 'FastRelax':
                python_code.append('fast_relax = FastRelax()')
                if 'scorefxn' in attrs:
                    python_code.append(f'fast_relax.score_function = get_score_function("{attrs["scorefxn"]}")')
            elif mover_name == 'MinMover':
                python_code.append('min_mover = MinMover()')
                if 'scorefxn' in attrs:
                    python_code.append(f'min_mover.score_function = get_score_function("{attrs["scorefxn"]}")')
            elif mover_name == 'PackRotamersMover':
                python_code.append('pack_rotamers = PackRotamersMover()')
                if 'scorefxn' in attrs:
                    python_code.append(f'pack_rotamers.score_function = get_score_function("{attrs["scorefxn"]}")')
        python_code.append('')
    
    # Create filters
    if filters:
        if ${includeCommentsPython}:
            python_code.append('# Create filters')
        for filter_name, attrs in filters:
            if filter_name == 'ScoreType':
                python_code.append('score_filter = ScoreType()')
                if 'score_type' in attrs:
                    python_code.append(f'score_filter.score_type = "{attrs["score_type"]}"')
                if 'threshold' in attrs:
                    python_code.append(f'score_filter.threshold = {attrs["threshold"]}')
        python_code.append('')
    
    # Create residue selectors
    if residue_selectors:
        if ${includeCommentsPython}:
            python_code.append('# Create residue selectors')
        for selector_name, attrs in residue_selectors:
            if selector_name == 'Chain':
                python_code.append('chain_selector = Chain()')
                if 'chains' in attrs:
                    python_code.append(f'chain_selector.chain = "{attrs["chains"]}"')
        python_code.append('')
    
    # Create task operations
    if task_operations:
        if ${includeCommentsPython}:
            python_code.append('# Create task operations')
        for op_name, attrs in task_operations:
            if op_name == 'RestrictToRepacking':
                python_code.append('restrict_repack = RestrictToRepacking()')
        python_code.append('')
    
    # Apply movers
    if movers:
        if ${includeCommentsPython}:
            python_code.append('# Apply movers to pose')
        for mover_name, _ in movers:
            if mover_name == 'FastRelax':
                python_code.append('fast_relax.apply(pose)')
            elif mover_name == 'MinMover':
                python_code.append('min_mover.apply(pose)')
            elif mover_name == 'PackRotamersMover':
                python_code.append('pack_rotamers.apply(pose)')
        python_code.append('')
    
    # Save result
    python_code.append('# Save the result')
    python_code.append('pose.dump_pdb("output.pdb")')
    python_code.append('')
    
    # Print score
    python_code.append('# Print final score')
    python_code.append('score = pose.energies().total_energy()')
    python_code.append('print(f"Final score: {score}")')
    
    # If no components found, provide a basic template
    if not movers and not filters and not residue_selectors and not task_operations:
        python_code = [
            '# Generated PyRosetta code from RosettaScripts XML',
            '# ==================================================',
            '#',
            '# No specific components found in XML, providing basic template:',
            '',
            'import pyrosetta',
            'from pyrosetta import pose_from_pdb',
            'from pyrosetta.rosetta.protocols.moves import *',
            'from pyrosetta.rosetta.protocols.simple_moves import *',
            'from pyrosetta.rosetta.protocols.filters import *',
            'from pyrosetta.rosetta.core.select.residue_selector import *',
            'from pyrosetta.rosetta.core.pack.task.operation import *',
            '',
            '# Initialize PyRosetta',
            'pyrosetta.init("-mute all")',
            '',
            '# Load your PDB file',
            'pose = pose_from_pdb("your_protein.pdb")',
            '',
            '# Add your PyRosetta code here based on the XML',
            '# Example:',
            '# fast_relax = FastRelax()',
            '# fast_relax.apply(pose)',
            '',
            '# Save the result',
            'pose.dump_pdb("output.pdb")',
            '',
            '# Print final score',
            'score = pose.energies().total_energy()',
            'print(f"Final score: {score}")'
        ]
    
    result = {
        'success': True,
        'python_code': '\\n'.join(python_code),
        'components_found': {
            'movers': len(movers),
            'filters': len(filters),
            'residue_selectors': len(residue_selectors),
            'task_operations': len(task_operations)
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

            const proc = spawn(py, ['-c', script]);
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
        // Detect whether the client uses Content-Length framing (LSP-style)
        // or newline-delimited JSON. Default to newline; switch to headers if detected.
        this.useHeaders = false;
        this.serverVersion = this.readPackageVersion();
        this.debug = String(process.env.MCP_DEBUG || '').trim().length > 0;
        this.setupMCPHandlers();
    }

    isToolAllowed(name) {
        // Reverted to always allow tools (restore 1.1.1 behavior)
        return true;
    }

    getToolFilterState() {
        const normalize = (s) => String(s || '')
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, '_')
            .replace(/^_+|_+$/g, '');
        const parseList = (s) => (s || '')
            .split(',')
            .map(x => x.trim())
            .filter(Boolean)
            .map(normalize);
        return {
            allowList: parseList(process.env.MCP_TOOLS),
            denyList: parseList(process.env.MCP_TOOLS_DENY)
        };
    }

    readPackageVersion() {
        try {
            const fs = require('fs');
            const path = require('path');
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
                        this.sendError('Failed to parse message', e.message || String(e));
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
                    this.sendError('Failed to parse message', e.message || String(e));
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
                            // Advertise both canonical and alternate capability names
                            tools: { 
                                list: true, 
                                call: true,
                                // Some clients expect these names
                                listTools: true,
                                callTool: true
                            },
                            resources: { 
                                list: true, 
                                read: true,
                                // Some clients expect these names
                                listResources: true,
                                readResource: true
                            }
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
                            },
                            {
                                name: 'pyrosetta_introspect',
                                description: 'Search PyRosetta API classes (movers/filters/selectors) and return docs and signatures',
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
                                description: 'Generate and cache RosettaScripts XML schema (XSD) and optionally extract element names',
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
                                name: 'cache_cli_docs',
                                description: 'Run rosetta_scripts help and save outputs for offline search',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        exe_path: { type: 'string', description: 'Path to rosetta_scripts executable (optional)' },
                                        cache_dir: { type: 'string', description: 'Directory to save docs' }
                                    }
                                }
                            },
                            {
                                name: 'get_cached_docs',
                                description: 'Search locally cached CLI docs for a query',
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
                                description: 'Get Python executable, version, and pip list for the MCP environment',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'check_pyrosetta',
                                description: 'Check if PyRosetta is importable in the MCP Python environment',
                                inputSchema: { type: 'object', properties: {} }
                            },
                            {
                                name: 'install_pyrosetta_installer',
                                description: 'Install PyRosetta using pyrosetta-installer (non-conda path)',
                                inputSchema: { type: 'object', properties: { silent: { type: 'boolean' } } }
                            },
                            {
                                name: 'find_rosetta_scripts',
                                description: 'Resolve the rosetta_scripts executable path (from exe_path/env/known dirs)',
                                inputSchema: { type: 'object', properties: { exe_path: { type: 'string' } } }
                            },
                            {
                                name: 'search_pyrosetta_wheels',
                                description: 'Search a directory for PyRosetta wheel files',
                                inputSchema: {
                                    properties: {
                                        directory: { type: 'string', description: 'Directory to search' }
                                    }
                                }
                            },
                            {
                                name: 'xml_to_pyrosetta',
                                description: 'Translate RosettaScripts XML to PyRosetta Python code',
                                inputSchema: {
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
                                description: 'Search online Rosetta documentation and return top matches',
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
                                description: 'Fetch and summarize a specific Rosetta docs URL',
                                inputSchema: {
                                    type: 'object',
                                    properties: {
                                        url: { type: 'string', description: 'Full URL to a Rosetta docs page' },
                                        max_chars: { type: 'number', description: 'Max characters of cleaned text to return (default 4000)' }
                                    },
                                    required: ['url']
                                }
                            }
                    ];
                    // Always expose all tools
                    result = { tools: allTools };
                    if (this.debug) console.error(`[mcp] -> tools/list count=${allTools.length}`);
                    break;
                }

                case 'resources/list':
                    result = { resources: [] };
                    break;
                case 'resources/read': {
                    const uri = params && params.uri ? params.uri : '';
                    result = { uri, mimeType: 'text/plain', text: '' };
                    break;
                }
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
                        case 'cache_cli_docs':
                            result = await this.rosettaServer.cacheCliDocs({
                                exe_path: args.exe_path,
                                cache_dir: args.cache_dir
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
                        case 'search_pyrosetta_wheels':
                            result = await this.rosettaServer.searchPyRosettaWheels({ directory: args && args.directory });
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
        const payload = JSON.stringify(response);
        if (this.useHeaders) {
            const head = `Content-Length: ${Buffer.byteLength(payload, 'utf8')}\r\n\r\n`;
            process.stdout.write(head);
            process.stdout.write(payload);
        } else {
            process.stdout.write(payload + '\n');
        }
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