<#
start_m5.ps1

Start the MapKinase m5 webapp (uvicorn) with the correct environment variables.

Usage examples:
# Start in foreground (logs printed to this console)
PS> .\start_m5.ps1

# Start in the background (new PowerShell window)
PS> .\start_m5.ps1 -Background

# Override file paths or port
PS> .\start_m5.ps1 -ProteinFile 'C:\path\to\BCp2_ProtMaphsaanno.txt' -PtmFile 'C:\path\to\BCp2_PhosMaphsanno.txt' -Port 8003 -Background
#>

Param(
    [string]$ProteinFile = "C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt",
    [string]$PtmFile = "C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt",
    [int]$Port = 8003,
    [switch]$Background
)

# Change to script directory (assumes the script lives in MapKinase_WebApp)
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
Set-Location $scriptDir

# Locate python: prefer system python, otherwise fallback to known Python 3.12 path
$pythonCmd = $null
try {
    $pythonCmd = (Get-Command python -ErrorAction SilentlyContinue).Source
} catch { $pythonCmd = $null }
if (-not $pythonCmd) {
    $fallback = 'C:\Users\clayt\AppData\Local\Programs\Python\Python312\python.exe'
    if (Test-Path $fallback) { $pythonCmd = $fallback } else {
        Write-Error "python executable not found. Please install Python or adjust the script to point to your python.exe."
        exit 1
    }
}

$uvicornArgs = "-m uvicorn m5_webapp:app.starlette_app --host 127.0.0.1 --port $Port --log-level info"

if (-not $Background) {
    Write-Host "Starting uvicorn in foreground (logs will appear here) with port $Port"
    $env:M5_PROTEIN_FILE = $ProteinFile
    $env:M5_PTM_FILE = $PtmFile
    Write-Host "M5_PROTEIN_FILE=$env:M5_PROTEIN_FILE"
    Write-Host "M5_PTM_FILE=$env:M5_PTM_FILE"
    Write-Host "Running: $pythonCmd $uvicornArgs"
    & $pythonCmd $uvicornArgs
} else {
    Write-Host "Starting uvicorn in background (a new PowerShell window will open)."
    # Build a command string that sets env vars in the new PowerShell and runs uvicorn
    $pwCmd = "`$env:M5_PROTEIN_FILE = '$ProteinFile'; `$env:M5_PTM_FILE = '$PtmFile'; Set-Location -Path '$scriptDir'; & '$pythonCmd' $uvicornArgs"
    Start-Process -FilePath 'powershell' -ArgumentList '-NoProfile','-NoExit','-Command',$pwCmd -WorkingDirectory $scriptDir -WindowStyle Normal -PassThru | Out-Null
    Write-Host "Background uvicorn launched (check the new PowerShell window for logs)."
}
