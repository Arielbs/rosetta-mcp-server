#!/bin/bash

echo "🧪 Testing XML to PyRosetta Translation"
echo "========================================"

# Test 1: Simple FastRelax
echo -e "\n1️⃣ Testing simple FastRelax XML:"
echo '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"xml_to_pyrosetta","arguments":{"xml_content":"<ROSETTASCRIPTS><MOVERS><FastRelax/></MOVERS></ROSETTASCRIPTS>"}}}' | rosetta-mcp-server

# Test 2: FastRelax with score function (using single quotes for XML attributes)
echo -e "\n2️⃣ Testing FastRelax with score function:"
echo '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"xml_to_pyrosetta","arguments":{"xml_content":"<ROSETTASCRIPTS><SCOREFXNS><ScoreFunction name=\"ref2015\" weights=\"ref2015\"/></SCOREFXNS><MOVERS><FastRelax name=\"relax\" scorefxn=\"ref2015\"/></MOVERS></ROSETTASCRIPTS>"}}}' | rosetta-mcp-server

# Test 3: Chain selector
echo -e "\n3️⃣ Testing Chain selector:"
echo '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"xml_to_pyrosetta","arguments":{"xml_content":"<ROSETTASCRIPTS><RESIDUE_SELECTORS><Chain chains=\"A\"/></RESIDUE_SELECTORS></ROSETTASCRIPTS>"}}}' | rosetta-mcp-server

# Test 4: Without comments
echo -e "\n4️⃣ Testing without comments:"
echo '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"xml_to_pyrosetta","arguments":{"xml_content":"<ROSETTASCRIPTS><MOVERS><FastRelax/></MOVERS></ROSETTASCRIPTS>","include_comments":false}}}' | rosetta-mcp-server

echo -e "\n✅ XML translation tests completed!"
