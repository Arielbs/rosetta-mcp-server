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
