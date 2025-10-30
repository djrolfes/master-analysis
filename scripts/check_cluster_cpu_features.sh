#!/bin/bash
# Check CPU features across all SLURM nodes to find which support AVX2, AVX512, etc.

echo "Checking CPU features on cluster nodes..."
echo "========================================="
echo ""

# Get list of nodes
NODES=$(sinfo -N -h -o "%N" | sort -u)

echo "Node Name       | CPU Model                              | AVX | AVX2 | AVX512"
echo "----------------|----------------------------------------|-----|------|--------"

for node in $NODES; do
    # Run lscpu on each node and extract relevant info
    result=$(srun --nodelist=$node -w $node timeout 10 lscpu 2>/dev/null || echo "TIMEOUT")
    
    if [[ "$result" == "TIMEOUT" ]]; then
        printf "%-15s | %-38s | %-3s | %-4s | %-6s\n" "$node" "UNAVAILABLE" "?" "?" "?"
        continue
    fi
    
    model=$(echo "$result" | grep "Model name:" | sed 's/Model name: *//' | cut -c1-38)
    flags=$(echo "$result" | grep "Flags:" | sed 's/Flags: *//')
    
    # Check for instruction sets
    has_avx="No"
    has_avx2="No"
    has_avx512="No"
    
    if echo "$flags" | grep -q " avx "; then has_avx="Yes"; fi
    if echo "$flags" | grep -q " avx2 "; then has_avx2="Yes"; fi
    if echo "$flags" | grep -q "avx512"; then has_avx512="Yes"; fi
    
    printf "%-15s | %-38s | %-3s | %-4s | %-6s\n" "$node" "$model" "$has_avx" "$has_avx2" "$has_avx512"
done

echo ""
echo "Recommendation for package installation:"
echo "- Nodes WITH AVX2: Can use pre-compiled R packages (faster)"
echo "- Nodes WITHOUT AVX2: Must compile from source with -march=native"
echo ""
echo "To request a node WITH AVX2 for installation:"
echo "  srun --nodelist=<node_with_avx2> ./scripts/check_and_install_hadron.sh"
