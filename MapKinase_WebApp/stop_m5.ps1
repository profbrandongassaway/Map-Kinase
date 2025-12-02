<#
stop_m5.ps1

Stop the MapKinase m5 webapp (uvicorn) that was started with `start_m5.ps1` or manually.
The script attempts to locate running python processes whose command line contains 'uvicorn' or 'm5_webapp'
and stops them. Optionally filter by port.

Usage:
# Stop any uvicorn/m5 processes found
PS> .\stop_m5.ps1

# Stop only processes whose command line mentions port 8003
PS> .\stop_m5.ps1 -Port 8003

# Force stop (default behavior uses -Force when stopping)
PS> .\stop_m5.ps1 -Force
#>

param(
    [int]$Port = 0,
    [switch]$Force
)

Write-Host "Searching for python/uvicorn processes (m5_webapp)..."

# Retrieve processes with a command line containing keywords
try {
    $procs = Get-CimInstance Win32_Process | Where-Object {
        ($_.CommandLine -and ($_.CommandLine -match 'uvicorn' -or $_.CommandLine -match 'm5_webapp' -or $_.CommandLine -match 'app.starlette_app'))
    }
} catch {
    Write-Host "Could not query processes: $_"
    exit 1
}

if ($Port -ne 0) {
    Write-Host "Filtering for processes mentioning port $Port"
    $procs = $procs | Where-Object { $_.CommandLine -match "$Port" }
}

if (-not $procs -or $procs.Count -eq 0) {
    Write-Host "No matching m5/uvicorn python processes found."
    exit 0
}

foreach ($p in $procs) {
    $pid = $p.ProcessId
    $cmd = $p.CommandLine
    Write-Host "Found PID=$pid -> $cmd"
    try {
        if ($Force) {
            Stop-Process -Id $pid -Force -ErrorAction Stop
        } else {
            Stop-Process -Id $pid -ErrorAction Stop
        }
        Write-Host "Stopped PID $pid"
    } catch {
        Write-Host "Failed to stop PID $pid: $_"
    }
}

Write-Host "Done."