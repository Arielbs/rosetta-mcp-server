#!/usr/bin/env python3
"""
Simple test script for the Rosetta MCP Server
"""

import json
import subprocess
import sys

def test_mcp_server():
    """Test basic MCP server functionality"""
    
    print("🧪 Testing Rosetta MCP Server...")
    print("=" * 40)
    
    # Test 1: Initialize
    print("\n1️⃣ Testing initialization...")
    init_msg = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "1"}
        }
    }
    
    try:
        result = subprocess.run(
            ['rosetta-mcp-server'],
            input=json.dumps(init_msg) + '\n',
            text=True,
            capture_output=True,
            timeout=10
        )
        
        if result.returncode == 0:
            print("✅ Initialization successful")
            print(f"   Response: {result.stdout.strip()}")
        else:
            print("❌ Initialization failed")
            print(f"   Error: {result.stderr.strip()}")
            
    except subprocess.TimeoutExpired:
        print("⏰ Initialization timed out")
    except Exception as e:
        print(f"💥 Initialization error: {e}")
    
    # Test 2: List tools
    print("\n2️⃣ Testing tools list...")
    tools_msg = {
        "jsonrpc": "2.0",
        "id": 2,
        "method": "tools/list",
        "params": {}
    }
    
    try:
        result = subprocess.run(
            ['rosetta-mcp-server'],
            input=json.dumps(tools_msg) + '\n',
            text=True,
            capture_output=True,
            timeout=10
        )
        
        if result.returncode == 0:
            print("✅ Tools list successful")
            # Parse and show PyRosetta tools
            try:
                response = json.loads(result.stdout.strip())
                if 'result' in response and 'tools' in response['result']:
                    pyrosetta_tools = [t for t in response['result']['tools'] if 'pyrosetta' in t['name']]
                    print(f"   Found {len(pyrosetta_tools)} PyRosetta tools:")
                    for tool in pyrosetta_tools:
                        print(f"     - {tool['name']}: {tool['description']}")
                else:
                    print("   Response structure unexpected")
            except json.JSONDecodeError:
                print("   Could not parse JSON response")
        else:
            print("❌ Tools list failed")
            print(f"   Error: {result.stderr.strip()}")
            
    except subprocess.TimeoutExpired:
        print("⏰ Tools list timed out")
    except Exception as e:
        print(f"💥 Tools list error: {e}")
    
    # Test 3: PyRosetta info
    print("\n3️⃣ Testing PyRosetta info...")
    info_msg = {
        "jsonrpc": "2.0",
        "id": 3,
        "method": "tools/call",
        "params": {
            "name": "check_pyrosetta",
            "arguments": {}
        }
    }
    
    try:
        result = subprocess.run(
            ['rosetta-mcp-server'],
            input=json.dumps(info_msg) + '\n',
            text=True,
            capture_output=True,
            timeout=10
        )
        
        if result.returncode == 0:
            print("✅ PyRosetta check successful")
            try:
                response = json.loads(result.stdout.strip())
                if 'result' in response:
                    if response['result'].get('available'):
                        print(f"   ✅ PyRosetta is available (v{response['result'].get('version', 'unknown')})")
                    else:
                        print(f"   ❌ PyRosetta not available: {response['result'].get('error', 'unknown error')}")
                else:
                    print("   Response structure unexpected")
            except json.JSONDecodeError:
                print("   Could not parse JSON response")
        else:
            print("❌ PyRosetta check failed")
            print(f"   Error: {result.stderr.strip()}")
            
    except subprocess.TimeoutExpired:
        print("⏰ PyRosetta check timed out")
    except Exception as e:
        print(f"💥 PyRosetta check error: {e}")
    
    print("\n🎯 Test completed!")

if __name__ == "__main__":
    test_mcp_server()
