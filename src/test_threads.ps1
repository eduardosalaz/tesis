# Define the output file name
$outputFileName = "output.txt"

# Get all the files in the specified directory
$files = Get-ChildItem -Path .\instances\1250_155_62\

# For each file in the directory
foreach ($file in $files) {
    # Run the Julia command with thread counts from 1 to 8
    for ($i = 1; $i -le 8; $i++) {
        $command = "julia -t $i .\src\heuristics\grasp_parallel.jl $($file.FullName) 30"
        Write-Host "Executing: $command"
        
        # Capture the output
        $output = Invoke-Expression $command
        
        # Format and write the output
        $formattedOutput = "$($file.Name)`t${i}`t${output}"
        Add-Content -Path $outputFileName -Value $formattedOutput
    }
}
