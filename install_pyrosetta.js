#!/usr/bin/env node
/**
 * Standalone PyRosetta installer for rosetta-mcp-server
 * Run this after npm install to pre-install PyRosetta
 */

const { RosettaMCPServer } = require('./rosetta_mcp_wrapper.js');
const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const os = require('os');

async function createFreshEnvironment() {
    const homeDir = os.homedir();
    const venvDir = path.join(homeDir, '.venv');
    const envDir = path.join(venvDir, 'pyrosetta_mcp_env');
    
    // Ensure .venv directory exists
    if (!fs.existsSync(venvDir)) {
        fs.mkdirSync(venvDir, { recursive: true });
    }
    
    console.log('\n🆕 Creating fresh Python environment...');
    console.log(`   Location: ${envDir}`);
    
    return new Promise((resolve, reject) => {
        // Use user's home directory and ensure proper permissions
        const proc = spawn('python3', ['-m', 'venv', envDir], {
            cwd: homeDir,
            stdio: 'inherit'
        });
        
        proc.on('close', (code) => {
            if (code === 0) {
                console.log('✅ Fresh environment created successfully!');
                resolve(envDir);
            } else {
                reject(new Error(`Failed to create environment (exit code: ${code})`));
            }
        });
        
        proc.on('error', (err) => {
            reject(new Error(`Failed to create environment: ${err.message}`));
        });
    });
}

async function installInEnvironment(envDir) {
    const pythonPath = path.join(envDir, 'bin', 'python');
    const pipPath = path.join(envDir, 'bin', 'pip');
    
    console.log('\n📦 Installing PyRosetta in fresh environment...');
    
    return new Promise((resolve, reject) => {
        // First upgrade pip in the new environment
        console.log('   Upgrading pip...');
        const upgradeProc = spawn(pipPath, ['install', '--upgrade', 'pip', 'setuptools', 'wheel']);
        
        upgradeProc.on('close', (code) => {
            if (code === 0) {
                console.log('✅ pip upgraded successfully!');
                
                // Now install pyrosetta-installer
                console.log('   Installing pyrosetta-installer...');
                const installProc = spawn(pipPath, ['install', 'pyrosetta-installer'], {
                    stdio: 'inherit'
                });
                
                installProc.on('close', (installCode) => {
                    if (installCode === 0) {
                        console.log('✅ pyrosetta-installer installed successfully!');
                        resolve(pythonPath);
                    } else {
                        reject(new Error(`Failed to install pyrosetta-installer (exit code: ${installCode})`));
                    }
                });
                
                installProc.on('error', (err) => {
                    reject(new Error(`Failed to install pyrosetta-installer: ${err.message}`));
                });
            } else {
                reject(new Error(`Failed to upgrade pip (exit code: ${code})`));
            }
        });
        
        upgradeProc.on('error', (err) => {
            reject(new Error(`Failed to upgrade pip: ${err.message}`));
        });
    });
}

async function checkSystemRequirements() {
    console.log('🔍 Checking system requirements...');
    
    // Check Python version
    try {
        const { spawn } = require('child_process');
        const pythonVersion = await new Promise((resolve, reject) => {
            const proc = spawn('python3', ['--version']);
            let output = '';
            proc.stdout.on('data', (d) => output += d.toString());
            proc.stderr.on('data', (d) => output += d.toString());
            proc.on('close', (code) => {
                if (code === 0) resolve(output.trim());
                else reject(new Error(`Python check failed (exit code: ${code})`));
            });
        });
        
        console.log(`   Python: ${pythonVersion}`);
        
        // Check if Python version is compatible (3.8+)
        const versionMatch = pythonVersion.match(/Python (\d+)\.(\d+)/);
        if (versionMatch) {
            const major = parseInt(versionMatch[1]);
            const minor = parseInt(versionMatch[2]);
            if (major < 3 || (major === 3 && minor < 8)) {
                throw new Error(`Python ${major}.${minor} is not supported. PyRosetta requires Python 3.8+`);
            }
        }
        
        // Check available disk space (need at least 2GB)
        const fs = require('fs');
        const homeDir = os.homedir();
        const stats = fs.statfsSync(homeDir);
        const freeGB = (stats.bavail * stats.bsize) / (1024 * 1024 * 1024);
        
        if (freeGB < 2) {
            console.warn(`   ⚠️  Low disk space: ${freeGB.toFixed(1)}GB available (recommended: 2GB+)`);
        } else {
            console.log(`   Disk space: ${freeGB.toFixed(1)}GB available ✅`);
        }
        
        console.log('✅ System requirements check passed!\n');
        
    } catch (error) {
        console.error('❌ System requirements check failed:');
        console.error('   Error:', error.message);
        console.error('\n🔧 Please ensure:');
        console.error('   - Python 3.8+ is installed and accessible');
        console.error('   - You have at least 2GB free disk space');
        console.error('   - You have proper permissions in your home directory');
        process.exit(1);
    }
}

async function main() {
    console.log('🚀 Rosetta MCP Server - PyRosetta Installer');
    console.log('============================================\n');
    
    // Check system requirements first
    await checkSystemRequirements();
    
    const server = new RosettaMCPServer();
    
    try {
        console.log('📋 Checking current PyRosetta status...');
        const status = await server.checkPyRosetta();
        
        if (status.available) {
            console.log(`✅ PyRosetta is already available${status.version ? ` (v${status.version})` : ''}`);
            console.log('   No installation needed!');
            return;
        }
        
        console.log('❌ PyRosetta not found. Starting installation...\n');
        
        // Check if current Python environment has issues
        console.log('🔍 Checking Python environment...');
        try {
            const envInfo = await server.pythonEnvInfo();
            if (envInfo.pip_error || envInfo.pip_list?.some(pkg => pkg.name === 'conda-repo-cli')) {
                console.log('⚠️  Detected potential environment conflicts (Anaconda/Conda issues)');
                console.log('   Creating fresh Python environment for PyRosetta...\n');
                
                const envDir = await createFreshEnvironment();
                const cleanPythonPath = await installInEnvironment(envDir);
                
                // Update server to use clean environment
                server.pythonPath = cleanPythonPath;
                
                console.log('\n⚠️  WARNING: PyRosetta installation will take 10-30 minutes!');
                console.log('   Please be patient and do not interrupt the process.\n');
                
                try {
                    const result = await server.installPyRosettaViaInstaller({ silent: false });
                    
                    if (result.ok) {
                        console.log('\n🎉 PyRosetta installation completed successfully!');
                        console.log(`   Environment: ${envDir}`);
                        console.log('   You can now use PyRosetta tools in the MCP server.');
                        console.log('\n💡 To use this environment in the future, set:');
                        console.log(`   export PYTHON_BIN="${cleanPythonPath}"`);
                    } else {
                        console.error('\n❌ PyRosetta installation failed:');
                        console.error('   Error:', result.error || 'Unknown error');
                        
                        // Provide troubleshooting tips
                        console.error('\n🔧 Troubleshooting tips:');
                        console.error('   1. Check if you have sufficient disk space');
                        console.error('   2. Ensure you have a stable internet connection');
                        console.error('   3. Try running: source ~/.venv/pyrosetta_mcp_env/bin/activate');
                        console.error('   4. Then manually: pip install pyrosetta-installer');
                        console.error('   5. Contact support if the issue persists');
                        
                        process.exit(1);
                    }
                } catch (installError) {
                    console.error('\n💥 PyRosetta installation crashed:');
                    console.error('   Error:', installError.message);
                    console.error('\n🔧 Try manual installation:');
                    console.error(`   source ${envDir}/bin/activate`);
                    console.error('   pip install pyrosetta-installer');
                    console.error('   python -c "import pyrosetta_installer as I; I.install_pyrosetta()"');
                    console.error('   Or use the standard location:');
                    console.error('   source ~/.venv/pyrosetta_mcp_env/bin/activate');
                    process.exit(1);
                }
            } else {
                // Standard installation in current environment
                console.log('⚠️  WARNING: This will take 10-30 minutes!');
                console.log('   Please be patient and do not interrupt the process.\n');
                
                const result = await server.installPyRosettaViaInstaller({ silent: false });
                
                if (result.ok) {
                    console.log('\n🎉 PyRosetta installation completed successfully!');
                    console.log('   You can now use PyRosetta tools in the MCP server.');
                } else {
                    console.error('\n❌ PyRosetta installation failed:');
                    console.error('   Error:', result.error || 'Unknown error');
                    process.exit(1);
                }
            }
        } catch (envError) {
            console.log('⚠️  Could not check environment, proceeding with standard installation...\n');
            
            console.log('⚠️  WARNING: This will take 10-30 minutes!');
            console.log('   Please be patient and do not interrupt the process.\n');
            
            const result = await server.installPyRosettaViaInstaller({ silent: false });
            
            if (result.ok) {
                console.log('\n🎉 PyRosetta installation completed successfully!');
                console.log('   You can now use PyRosetta tools in the MCP server.');
            } else {
                console.error('\n❌ PyRosetta installation failed:');
                console.error('   Error:', result.error || 'Unknown error');
                process.exit(1);
            }
        }
        
    } catch (error) {
        console.error('\n💥 Installation failed with error:', error.message);
        process.exit(1);
    }
}

if (require.main === module) {
    main();
}

module.exports = { main };
