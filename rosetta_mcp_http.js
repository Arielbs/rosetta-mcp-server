#!/usr/bin/env node
/**
 * Rosetta MCP Server - Streamable HTTP Transport
 *
 * Serves a subset of Rosetta MCP tools over HTTP for use as a
 * remote MCP connector (Claude.ai, Smithery, etc.).
 *
 * Remote tools (no local deps needed): help, validation, translation, docs, Biotite
 * Local-only tools: scoring, introspection, running protocols (redirected to get_local_setup)
 */

const http = require('http');
const { RosettaMCPServer } = require('./rosetta_mcp_wrapper');

const PORT = parseInt(process.env.PORT || '8080', 10);
const RATE_LIMIT = 1000; // requests per hour per IP
const RATE_WINDOW = 60 * 60 * 1000; // 1 hour in ms

// Tools that work without local PyRosetta/Rosetta binary
const REMOTE_TOOLS = new Set([
    'get_rosetta_info', 'get_rosetta_help', 'validate_xml', 'xml_to_pyrosetta',
    'rosetta_to_biotite', 'biotite_to_rosetta', 'translate_rosetta_script_to_biotite',
    'search_rosetta_web_docs', 'get_rosetta_web_doc', 'get_local_setup'
]);

// Rate limiter: per-IP request counts
const rateLimits = new Map();

function checkRateLimit(ip) {
    const now = Date.now();
    let entry = rateLimits.get(ip);
    if (!entry || now - entry.windowStart > RATE_WINDOW) {
        entry = { count: 0, windowStart: now };
        rateLimits.set(ip, entry);
    }
    entry.count++;
    return entry.count <= RATE_LIMIT;
}

// Clean up expired rate limit entries every 10 minutes
setInterval(() => {
    const now = Date.now();
    for (const [ip, entry] of rateLimits) {
        if (now - entry.windowStart > RATE_WINDOW) rateLimits.delete(ip);
    }
}, 10 * 60 * 1000);

// The Rosetta server instance (shared across requests)
const rosettaServer = new RosettaMCPServer();

// Read package version
let serverVersion = 'unknown';
try {
    serverVersion = require('./package.json').version;
} catch (_) {}

// Build the server card for /.well-known/mcp/server-card.json
function getServerCard() {
    return {
        serverInfo: { name: 'rosetta-mcp-server', version: serverVersion },
        authentication: { required: false },
        tools: Array.from(REMOTE_TOOLS).map(name => ({
            name,
            description: `Rosetta MCP tool: ${name}`,
            inputSchema: { type: 'object' }
        })),
        resources: [],
        prompts: []
    };
}

// Handle MCP JSON-RPC messages (filtered for remote tools only)
async function handleMCPRequest(message) {
    const { id, method, params } = message;

    switch (method) {
        case undefined:
            return null; // notification, no response
        case 'initialize': {
            const requestedProtocol = (params && params.protocolVersion) ? params.protocolVersion : '2024-11-05';
            return {
                jsonrpc: '2.0', id,
                result: {
                    protocolVersion: requestedProtocol,
                    capabilities: { tools: { list: true, call: true } },
                    serverInfo: { name: 'rosetta-mcp-server', version: serverVersion }
                }
            };
        }
        case 'tools/list': {
            // Only expose remote-safe tools
            const allTools = require('./rosetta_mcp_wrapper').RosettaMCPServerMCP;
            // Build tool list inline from REMOTE_TOOLS set
            const tools = [];
            // We need to get the tool definitions - instantiate a temporary MCP handler to extract them
            // Instead, define them directly from the wrapper's exported tools
            const tempHandler = { rosettaServer, useHeaders: false, debug: false, serverVersion };
            // Get full tool list by sending tools/list through the wrapper
            const fullResponse = await new Promise((resolve) => {
                const origWrite = process.stdout.write;
                let captured = '';
                process.stdout.write = (data) => { captured += data; return true; };
                const mcpServer = new (require('./rosetta_mcp_wrapper').RosettaMCPServerMCP)();
                // Restore stdout immediately
                process.stdout.write = origWrite;
                // Parse the tools from the captured initialize response...
                // This approach is too hacky. Let's just define the filtered tools directly.
                resolve(null);
            });

            // Direct approach: get tool definitions from the allTools array
            // We import and filter from the tools/list response
            const toolDefs = getRemoteToolDefs();
            return {
                jsonrpc: '2.0', id,
                result: { tools: toolDefs }
            };
        }
        case 'tools/call': {
            const { name, arguments: args } = params || {};
            try {
                let result;
                if (name === 'get_local_setup' || !REMOTE_TOOLS.has(name)) {
                    // Return setup instructions for local-only tools
                    result = {
                        message: `The tool "${name}" requires a local installation of the Rosetta MCP server with PyRosetta.`,
                        install_steps: [
                            '1. npm install -g rosetta-mcp-server',
                            '2. Set up Python env: pip install pyrosetta-installer biotite && python -c "import pyrosetta_installer as I; I.install_pyrosetta()"',
                            '3. Add to your MCP client config (e.g., ~/.cursor/mcp.json):',
                            '   { "mcpServers": { "rosetta": { "command": "rosetta-mcp-server", "env": { "PYTHON_BIN": "/path/to/python" } } } }',
                            '4. Restart your editor'
                        ],
                        local_only_tools: ['pyrosetta_score', 'pyrosetta_introspect', 'run_rosetta_scripts', 'rosetta_scripts_schema', 'check_pyrosetta', 'install_pyrosetta_installer', 'find_rosetta_scripts', 'get_cached_docs', 'python_env_info'],
                        npm_package: 'rosetta-mcp-server',
                        docs: 'https://github.com/Arielbs/rosetta-mcp-server'
                    };
                } else {
                    // Dispatch to the actual tool implementation
                    switch (name) {
                        case 'get_rosetta_info':
                            result = await rosettaServer.getRosettaInfo(); break;
                        case 'get_rosetta_help':
                            result = await rosettaServer.getRosettaHelp(args && args.topic); break;
                        case 'validate_xml':
                            result = await rosettaServer.validateXML(args.xml_content, args.validate_against_schema); break;
                        case 'xml_to_pyrosetta':
                            result = await rosettaServer.xmlToPyRosetta({ xml_content: args.xml_content, include_comments: args.include_comments, output_format: args.output_format }); break;
                        case 'rosetta_to_biotite':
                            result = await rosettaServer.rosettaToBiotite({ query: args.query, category: args.category }); break;
                        case 'biotite_to_rosetta':
                            result = await rosettaServer.biotiteToRosetta({ query: args.query, category: args.category }); break;
                        case 'translate_rosetta_script_to_biotite':
                            result = await rosettaServer.translateRosettaScriptToBiotite({ code: args.code, input_format: args.input_format, include_comments: args.include_comments }); break;
                        case 'search_rosetta_web_docs':
                            result = await rosettaServer.searchRosettaWebDocs({ query: args.query, max_results: args.max_results }); break;
                        case 'get_rosetta_web_doc':
                            result = await rosettaServer.getRosettaWebDoc({ url: args.url, max_chars: args.max_chars }); break;
                        case 'get_local_setup':
                            // Already handled above
                            break;
                        default:
                            throw new Error(`Unknown tool: ${name}`);
                    }
                }
                return {
                    jsonrpc: '2.0', id,
                    result: {
                        content: [{ type: 'text', text: typeof result === 'string' ? result : JSON.stringify(result, null, 2) }]
                    }
                };
            } catch (e) {
                return {
                    jsonrpc: '2.0', id,
                    result: {
                        content: [{ type: 'text', text: e.message || String(e) }],
                        isError: true
                    }
                };
            }
        }
        case 'ping':
            return { jsonrpc: '2.0', id, result: { ok: true } };
        case 'logging/setLevel':
        case 'shutdown':
            return { jsonrpc: '2.0', id, result: {} };
        default:
            if (typeof method === 'string' && method.startsWith('notifications/')) {
                return null; // no response for notifications
            }
            return {
                jsonrpc: '2.0', id,
                error: { code: -32601, message: `Unknown method: ${method}` }
            };
    }
}

// Remote tool definitions (extracted to avoid circular dependency)
function getRemoteToolDefs() {
    return [
        { name: 'get_rosetta_info', description: 'Get comprehensive Rosetta installation info including available score functions, movers, filters, selectors, task operations, parameters, and command-line options. Use this first to understand what Rosetta components are available.', inputSchema: { type: 'object', properties: {} }, annotations: { readOnlyHint: true, openWorldHint: false } },
        { name: 'get_rosetta_help', description: 'Get help for any Rosetta topic, mover, filter, or concept. Accepts general topics (score_functions, movers, filters, xml, parameters) or specific names (FastRelax, Ddg, ChainSelector). Auto-fetches live documentation.', inputSchema: { type: 'object', properties: { topic: { type: 'string', description: 'Topic to get help for' } } }, annotations: { readOnlyHint: true, openWorldHint: true } },
        { name: 'validate_xml', description: 'Validate a RosettaScripts XML protocol. Checks XML syntax and optionally validates element names against the Rosetta XSD schema.', inputSchema: { type: 'object', properties: { xml_content: { type: 'string', description: 'XML content to validate' }, validate_against_schema: { type: 'boolean', description: 'Also check element names against Rosetta XSD schema (default: false)' } }, required: ['xml_content'] }, annotations: { readOnlyHint: true, openWorldHint: false } },
        { name: 'xml_to_pyrosetta', description: 'Translate RosettaScripts XML to equivalent PyRosetta Python code. Supports 37 element types.', inputSchema: { type: 'object', properties: { xml_content: { type: 'string', description: 'RosettaScripts XML content to translate' }, include_comments: { type: 'boolean', description: 'Include comments (default: true)' } }, required: ['xml_content'] }, annotations: { readOnlyHint: true, openWorldHint: false } },
        { name: 'rosetta_to_biotite', description: 'ALWAYS use this tool when asked about Rosetta vs Biotite equivalents. Returns the Biotite equivalent with working example code. Covers: SASA, RMSD, superimposition, secondary structure, distances, angles, contacts, hydrogen bonds, B-factors, and more.', inputSchema: { type: 'object', properties: { query: { type: 'string', description: 'Rosetta method or concept name' }, category: { type: 'string', description: 'Optional category filter' } }, required: ['query'] }, annotations: { readOnlyHint: true, openWorldHint: false } },
        { name: 'biotite_to_rosetta', description: 'ALWAYS use this tool when asked about Biotite vs Rosetta equivalents. Returns the Rosetta/PyRosetta equivalent with working example code.', inputSchema: { type: 'object', properties: { query: { type: 'string', description: 'Biotite function or concept name' }, category: { type: 'string', description: 'Optional category filter' } }, required: ['query'] }, annotations: { readOnlyHint: true, openWorldHint: false } },
        { name: 'translate_rosetta_script_to_biotite', description: 'Translate RosettaScripts XML or PyRosetta code to Biotite Python. Analysis operations translated; design/optimization flagged as Rosetta-only.', inputSchema: { type: 'object', properties: { code: { type: 'string', description: 'RosettaScripts XML or PyRosetta code' }, input_format: { type: 'string', enum: ['xml', 'pyrosetta', 'auto'], description: 'Input format (default: auto)' } }, required: ['code'] }, annotations: { readOnlyHint: true, openWorldHint: false } },
        { name: 'search_rosetta_web_docs', description: 'Search online Rosetta documentation at rosettacommons.org.', inputSchema: { type: 'object', properties: { query: { type: 'string', description: 'Search query' }, max_results: { type: 'number', description: 'Number of results (default 3)' } }, required: ['query'] }, annotations: { readOnlyHint: true, openWorldHint: true } },
        { name: 'get_rosetta_web_doc', description: 'Fetch and extract text from a specific Rosetta docs URL.', inputSchema: { type: 'object', properties: { url: { type: 'string', description: 'Full URL to a Rosetta docs page' }, max_chars: { type: 'number', description: 'Max characters (default 4000)' } }, required: ['url'] }, annotations: { readOnlyHint: true, openWorldHint: true } },
        { name: 'get_local_setup', description: 'Get instructions for installing the full Rosetta MCP server locally with PyRosetta support for scoring, optimization, introspection, and protocol execution.', inputSchema: { type: 'object', properties: {} }, annotations: { readOnlyHint: true, openWorldHint: false } },
    ];
}

// HTTP server
const server = http.createServer(async (req, res) => {
    const ip = req.headers['x-forwarded-for'] || req.socket.remoteAddress || 'unknown';
    const start = Date.now();

    // CORS headers
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET, POST, DELETE, OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type, Accept, Mcp-Session-Id');

    if (req.method === 'OPTIONS') {
        res.writeHead(204);
        res.end();
        return;
    }

    const url = req.url.split('?')[0];

    // Health check
    if (url === '/health' && req.method === 'GET') {
        const body = JSON.stringify({ status: 'ok', version: serverVersion, tools: REMOTE_TOOLS.size });
        res.writeHead(200, { 'Content-Type': 'application/json' });
        res.end(body);
        return;
    }

    // Server card for Smithery
    if (url === '/.well-known/mcp/server-card.json' && req.method === 'GET') {
        const body = JSON.stringify(getServerCard(), null, 2);
        res.writeHead(200, { 'Content-Type': 'application/json' });
        res.end(body);
        return;
    }

    // MCP endpoint
    if (url === '/mcp') {
        // Rate limit check
        if (!checkRateLimit(ip)) {
            res.writeHead(429, { 'Content-Type': 'application/json' });
            res.end(JSON.stringify({ error: 'Rate limit exceeded. Max 1000 requests per hour.' }));
            console.error(`[http] 429 rate-limited ip=${ip}`);
            return;
        }

        if (req.method === 'POST') {
            let body = '';
            req.on('data', chunk => { body += chunk; });
            req.on('end', async () => {
                try {
                    const message = JSON.parse(body);
                    const response = await handleMCPRequest(message);

                    if (response === null) {
                        // Notification -- no response needed
                        res.writeHead(202);
                        res.end();
                    } else {
                        const responseBody = JSON.stringify(response);
                        res.writeHead(200, { 'Content-Type': 'application/json' });
                        res.end(responseBody);
                    }

                    const elapsed = Date.now() - start;
                    const toolName = message.method === 'tools/call' ? (message.params && message.params.name) : message.method;
                    console.error(`[http] ${req.method} /mcp ip=${ip} tool=${toolName} ${elapsed}ms`);

                } catch (e) {
                    const errorResponse = JSON.stringify({
                        jsonrpc: '2.0', id: null,
                        error: { code: -32700, message: 'Parse error: ' + e.message }
                    });
                    res.writeHead(400, { 'Content-Type': 'application/json' });
                    res.end(errorResponse);
                    console.error(`[http] 400 parse-error ip=${ip} err=${e.message}`);
                }
            });
        } else if (req.method === 'GET') {
            // Streamable HTTP spec: GET opens SSE stream. We don't need it.
            res.writeHead(405, { 'Content-Type': 'application/json' });
            res.end(JSON.stringify({ error: 'Method not allowed. Use POST for MCP requests.' }));
        } else if (req.method === 'DELETE') {
            // Session termination -- stateless, not supported
            res.writeHead(405, { 'Content-Type': 'application/json' });
            res.end(JSON.stringify({ error: 'Session termination not supported (stateless server).' }));
        } else {
            res.writeHead(405);
            res.end();
        }
        return;
    }

    // Landing page for browser visitors
    if (url === '/' && req.method === 'GET') {
        res.writeHead(200, { 'Content-Type': 'text/html' });
        res.end(`<!DOCTYPE html>
<html><head><title>Rosetta MCP Server | Ariel J. Ben-Sasson</title>
<style>body{font-family:'Inter',-apple-system,BlinkMacSystemFont,sans-serif;max-width:700px;margin:0 auto;padding:60px 20px;color:#2D2A26;line-height:1.6;background-color:#FAF6F1;min-height:100vh}
body::before{content:'';position:fixed;top:0;left:0;width:100%;height:100%;pointer-events:none;z-index:9999;opacity:0.22;background-image:url("data:image/svg+xml,%3Csvg viewBox='0 0 256 256' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='noise'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.7' numOctaves='3' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23noise)'/%3E%3C/svg%3E")}
h1{color:#2D2A26}h2{color:#6B6560}a{color:#B8895A}a:hover{color:#A07848}code{background:rgba(232,213,196,0.4);padding:2px 6px;border-radius:4px;font-size:0.9em;color:#2D2A26}
.tools{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin:16px 0}.tool{background:linear-gradient(145deg,rgba(240,235,227,0.85),rgba(250,246,241,0.9),rgba(232,213,196,0.7));padding:8px 12px;border-radius:6px;font-size:0.9em;border:1px solid rgba(212,165,116,0.2)}
.tool.bridge{background:linear-gradient(145deg,rgba(212,165,116,0.25),rgba(232,213,196,0.4));border:1px dashed #B8895A}
.badge{display:inline-block;background:linear-gradient(135deg,#B8895A,#C9A06C);color:white;padding:2px 8px;border-radius:12px;font-size:0.8em}
.local-section{position:relative;margin:32px 0}
.local-header{display:flex;align-items:center;gap:12px;margin-bottom:12px}
.local-label{background:rgba(212,165,116,0.15);border:1px solid rgba(212,165,116,0.3);padding:4px 12px;border-radius:6px;font-size:0.85em;color:#6B6560}
.local-tools{display:grid;grid-template-columns:1fr 1fr 1fr;gap:6px}.local-tool{background:rgba(240,235,227,0.5);padding:6px 10px;border-radius:5px;font-size:0.8em;color:#6B6560;border:1px dashed rgba(212,165,116,0.25)}
.arrow-wrapper{position:relative;width:100%;height:90px;margin:-4px 0 -20px 0}
.arrow-wrapper svg{width:100%;height:100%}
.byline{font-size:0.95em;color:#6B6560;margin-top:-4px}</style></head>
<body>
<h1>Rosetta MCP Server <span class="badge">v${serverVersion}</span></h1>
<p class="byline">by <strong>Ariel J. Ben-Sasson</strong> (<a href="https://www.molcore.bio">molCore</a>) &mdash; for the RosettaCommons community &amp; any bio-engineers out there</p>
<p>Give your AI coding agent expert-level knowledge of protein modeling and design. Ask naive questions, get working code, translate and validate across Rosetta, PyRosetta, and Biotite.</p>
<p>This is a <a href="https://modelcontextprotocol.io">Model Context Protocol</a> server providing ${REMOTE_TOOLS.size} remote tools and 9 additional local tools.</p>

<h2>Connect your AI agent</h2>
<p><strong>Claude.ai:</strong> Settings &rarr; Connectors &rarr; Add &rarr; <code>https://mcp.molcore.bio/mcp</code></p>
<p><strong>Full local install</strong> (with scoring &amp; protocols): <code>npm install -g rosetta-mcp-server</code></p>

<h2>Available remotely</h2>
<div class="tools">
<div class="tool"><strong>get_rosetta_info</strong></div>
<div class="tool"><strong>get_rosetta_help</strong></div>
<div class="tool"><strong>validate_xml</strong></div>
<div class="tool"><strong>xml_to_pyrosetta</strong></div>
<div class="tool"><strong>rosetta_to_biotite</strong></div>
<div class="tool"><strong>biotite_to_rosetta</strong></div>
<div class="tool"><strong>translate_rosetta_script_to_biotite</strong></div>
<div class="tool"><strong>search_rosetta_web_docs</strong></div>
<div class="tool"><strong>get_rosetta_web_doc</strong></div>
<div class="tool bridge"><strong>get_local_setup</strong></div>
</div>

<div class="arrow-wrapper">
<svg viewBox="0 0 700 90" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M 530 2 C 520 20, 400 35, 250 38 S 100 42, 75 55 C 65 62, 60 70, 58 82" stroke="#B8895A" stroke-width="2" stroke-dasharray="5 4" fill="none" stroke-linecap="round"/>
<path d="M 52 75 L 58 86 L 64 75" stroke="#B8895A" stroke-width="2" fill="none" stroke-linecap="round" stroke-linejoin="round"/>
<text x="280" y="30" font-family="Inter,sans-serif" font-size="11" fill="#6B6560" font-style="italic" text-anchor="middle">sets up locally</text>
</svg>
</div>

<div class="local-section">
<div class="local-header">
<h2 style="margin:0">Requires local install</h2>
<span class="local-label">npm install -g rosetta-mcp-server</span>
</div>
<p style="font-size:0.9em;color:#6B6560;margin:4px 0 12px">These tools need PyRosetta or the Rosetta binary on your machine.</p>
<div class="local-tools">
<div class="local-tool">pyrosetta_score</div>
<div class="local-tool">pyrosetta_introspect</div>
<div class="local-tool">run_rosetta_scripts</div>
<div class="local-tool">rosetta_scripts_schema</div>
<div class="local-tool">check_pyrosetta</div>
<div class="local-tool">install_pyrosetta_installer</div>
<div class="local-tool">find_rosetta_scripts</div>
<div class="local-tool">get_cached_docs</div>
<div class="local-tool">python_env_info</div>
</div>
</div>

<h2>Links</h2>
<ul>
<li><a href="https://github.com/Arielbs/rosetta-mcp-server">GitHub</a></li>
<li><a href="https://www.npmjs.com/package/rosetta-mcp-server">npm package</a></li>
<li><a href="/health">Server status</a></li>
<li><a href="https://www.molcore.bio">molCore</a></li>
</ul>
</body></html>`);
        return;
    }

    // Unknown path
    res.writeHead(404, { 'Content-Type': 'application/json' });
    res.end(JSON.stringify({ error: 'Not found. Use POST /mcp for MCP requests, GET /health for status.' }));
});

server.listen(PORT, () => {
    console.error(`[rosetta-mcp-http] Listening on port ${PORT}`);
    console.error(`[rosetta-mcp-http] MCP endpoint: http://localhost:${PORT}/mcp`);
    console.error(`[rosetta-mcp-http] Health: http://localhost:${PORT}/health`);
    console.error(`[rosetta-mcp-http] Server card: http://localhost:${PORT}/.well-known/mcp/server-card.json`);
    console.error(`[rosetta-mcp-http] Remote tools: ${REMOTE_TOOLS.size}`);
    console.error(`[rosetta-mcp-http] Rate limit: ${RATE_LIMIT} req/hour per IP`);
});

process.on('SIGINT', () => {
    console.error('[rosetta-mcp-http] Shutting down...');
    server.close();
    process.exit(0);
});
