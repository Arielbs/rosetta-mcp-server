#!/usr/bin/env python3
"""
Test script for XML to PyRosetta translation
"""

import json
import subprocess
import sys

def test_xml_translation():
    """Test the XML to PyRosetta translation functionality"""
    
    print("🧪 Testing XML to PyRosetta Translation")
    print("=" * 50)
    
    # Test cases with properly formatted XML
    test_cases = [
        {
            "name": "Simple FastRelax",
            "xml": "<ROSETTASCRIPTS><MOVERS><FastRelax/></MOVERS></ROSETTASCRIPTS>"
        },
        {
            "name": "FastRelax with score function",
            "xml": '<ROSETTASCRIPTS><SCOREFXNS><ScoreFunction name="ref2015" weights="ref2015"/></SCOREFXNS><MOVERS><FastRelax name="relax" scorefxn="ref2015"/></MOVERS></ROSETTASCRIPTS>'
        },
        {
            "name": "Multiple components",
            "xml": '<ROSETTASCRIPTS><MOVERS><FastRelax name="relax"/><MinMover name="min"/></MOVERS><FILTERS><ScoreType score_type="total_score" threshold="0.0"/></FILTERS></ROSETTASCRIPTS>'
        },
        {
            "name": "Chain selector",
            "xml": '<ROSETTASCRIPTS><RESIDUE_SELECTORS><Chain chains="A"/></RESIDUE_SELECTORS></ROSETTASCRIPTS>'
        }
    ]
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n{i}️⃣ Testing: {test_case['name']}")
        print(f"XML: {test_case['xml']}")
        
        # Create the MCP request
        request = {
            "jsonrpc": "2.0",
            "id": i,
            "method": "tools/call",
            "params": {
                "name": "xml_to_pyrosetta",
                "arguments": {
                    "xml_content": test_case['xml'],
                    "include_comments": True
                }
            }
        }
        
        try:
            # Send request to rosetta-mcp-server
            result = subprocess.run(
                ['rosetta-mcp-server'],
                input=json.dumps(request),
                text=True,
                capture_output=True,
                timeout=30
            )
            
            if result.returncode == 0:
                try:
                    response = json.loads(result.stdout.strip())
                    if 'result' in response:
                        result_data = response['result']
                        if result_data.get('success'):
                            print("✅ Success!")
                            print(f"Components found: {result_data.get('components_found', {})}")
                            print(f"Python code length: {len(result_data.get('python_code', ''))} characters")
                            if result_data.get('python_code'):
                                print("Generated Python code preview:")
                                lines = result_data['python_code'].split('\n')[:10]
                                for line in lines:
                                    print(f"  {line}")
                                if len(result_data['python_code'].split('\n')) > 10:
                                    print("  ... (truncated)")
                        else:
                            print("❌ Failed:")
                            print(f"Error: {result_data.get('error', 'Unknown error')}")
                    else:
                        print("❌ No result in response")
                        print(f"Response: {response}")
                except json.JSONDecodeError:
                    print("❌ Failed to parse JSON response")
                    print(f"Raw output: {result.stdout}")
            else:
                print("❌ Command failed")
                print(f"Error: {result.stderr}")
                
        except subprocess.TimeoutExpired:
            print("❌ Request timed out")
        except FileNotFoundError:
            print("❌ rosetta-mcp-server not found. Make sure it's installed and in your PATH")
            break
        except Exception as e:
            print(f"❌ Unexpected error: {e}")
    
    print("\n✅ XML translation tests completed!")

if __name__ == "__main__":
    test_xml_translation()
