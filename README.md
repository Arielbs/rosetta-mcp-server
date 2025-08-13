# Rosetta MCP Server

A Model Context Protocol (MCP) server that provides comprehensive access to Rosetta/PyRosetta functions, properties, and capabilities.

## Features

- **Complete Rosetta Function Database**: Access to all score functions, movers, filters, residue selectors, and task operations
- **XML Protocol Validation**: Validate your Rosetta XML protocol files
- **Interactive Help System**: Get detailed help on any Rosetta topic
- **Path Detection**: Automatically finds your Rosetta installation
- **PyRosetta Integration**: Detects if PyRosetta is available and provides enhanced functionality
- **Global Installation**: Available system-wide on your computer

## Installation

### Prerequisites

- Python 3.7+
- Node.js 14+
- Rosetta installation (optional but recommended)
- PyRosetta (optional but recommended)

### Global Installation

1. **Clone or download this repository:**
   ```bash
   cd ~/rosetta_mcp_server
   ```

2. **Install globally:**
   ```bash
   npm install -g .
   ```

3. **Verify installation:**
   ```bash
   rosetta-mcp-server --help
   ```

### Manual Installation

1. **Make the Python script executable:**
   ```bash
   chmod +x rosetta_mcp_server.py
   ```

2. **Test the Python server:**
   ```bash
   python3 rosetta_mcp_server.py info
   ```

## Usage

### Command Line Interface

```bash
# Get comprehensive Rosetta information
rosetta-mcp-server info

# Get help on specific topics
rosetta-mcp-server help score_functions
rosetta-mcp-server help movers
rosetta-mcp-server help filters

# Validate XML protocol files
rosetta-mcp-server validate path/to/protocol.xml
```

### MCP Integration

Add to your `~/.cursor/mcp.json`:

```json
{
  "mcpServers": {
    "rosetta": {
      "command": "rosetta-mcp-server",
      "args": []
    }
  }
}
```

## Available Tools

### 1. get_rosetta_info
Get comprehensive information about your Rosetta setup:
- Rosetta installation path
- PyRosetta availability
- Available score functions
- Common movers and filters
- Residue selectors and task operations
- Command line parameters

### 2. get_rosetta_help
Get detailed help on specific topics:
- `score_functions`: Score function explanations
- `movers`: Mover descriptions and usage
- `filters`: Filter explanations and parameters
- `xml`: XML protocol structure
- `parameters`: Command line parameter details

### 3. validate_xml
Validate your Rosetta XML protocol files:
- Syntax validation
- Structure checking
- Error reporting

### 4. list_functions
Get categorized lists of available functions:
- Score functions
- Movers
- Filters
- Residue selectors
- Task operations
- Parameters

## Configuration

The server automatically detects:
- Rosetta installation paths
- PyRosetta availability
- Available XML protocols
- System-specific configurations

## Examples

### Get All Available Rosetta Functions
```bash
rosetta-mcp-server info | jq '.common_movers'
```

### Validate Your Protocol
```bash
rosetta-mcp-server validate inputs/xml/your_protocol.xml
```

### Get Help on Score Functions
```bash
rosetta-mcp-server help score_functions
```

## Troubleshooting

### Common Issues

1. **Python not found**: Ensure `python3` is in your PATH
2. **Permission denied**: Use `sudo npm install -g .` for global installation
3. **Rosetta path not found**: The server will work with default values

### Testing

```bash
# Test Python server
python3 rosetta_mcp_server.py info

# Test Node.js wrapper
node rosetta_mcp_wrapper.js

# Test global installation
rosetta-mcp-server info
```

## Development

### Project Structure
```
rosetta_mcp_server/
├── rosetta_mcp_server.py    # Python server with Rosetta knowledge
├── rosetta_mcp_wrapper.js   # Node.js MCP wrapper
├── package.json             # NPM package configuration
└── README.md               # This file
```

### Adding New Functions

1. **Add to Python server** (`rosetta_mcp_server.py`)
2. **Add to Node.js wrapper** (`rosetta_mcp_wrapper.js`)
3. **Update MCP tools list**
4. **Test and document**

## License

MIT License - see LICENSE file for details.

## Support

For issues or questions:
1. Check the troubleshooting section
2. Test individual components
3. Verify your Rosetta installation
4. Check MCP configuration in Cursor

## Contributing

Contributions welcome! Please:
1. Test your changes thoroughly
2. Update documentation
3. Follow existing code style
4. Add appropriate error handling
