import subprocess

tasks = [
    'distance_matrix',
    'inorder',
    'ladderize',
    'levelorder',
    'load_tree',
    'memory',
    'mrca',
    'postorder',
    'preorder',
    'rootdistorder',
    'total_branch_length'
]

# Variables for dynamic command
X = 100  # Example value for X
tool = 'treeswift_newick'  # Example value for tool

for task in tasks:
    # Build the command
    command = ['python3', 'scripts/our_time.py', 'data/tree_n%d.tre.gz' % X, tool, task]

    # Print the command for debugging
    print(f"Running command: {' '.join(command)}")  # This will print the full command to run

    # Run the command using subprocess
    try:
        measurement = subprocess.check_output(command, stderr=subprocess.PIPE).decode().strip()
        print(f"Measurement: {measurement}")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        print(f"Standard Output: {e.output.decode()}")
        print(f"Standard Error: {e.stderr.decode()}")
