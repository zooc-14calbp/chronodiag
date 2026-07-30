[CmdletBinding()]
param(
    [string]$ProjectRoot = "D:\Dropbox\papers_2026\03f_eps",
    [string]$PythonExe,
    [switch]$SkipPlots,
    [switch]$SkipFullValidation
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Invoke-Native {
    param(
        [Parameter(Mandatory = $true)]
        [string]$FilePath,
        [Parameter(Mandatory = $false)]
        [string[]]$Arguments = @(),
        [Parameter(Mandatory = $false)]
        [string]$FailureMessage = "Native command failed."
    )
    $oldPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & $FilePath @Arguments
        $code = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $oldPreference
    }
    if ($code -ne 0) {
        throw "$FailureMessage Exit code: $code"
    }
}

$RepoRoot = Split-Path -Parent $PSScriptRoot
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path

if ([string]::IsNullOrWhiteSpace($PythonExe)) {
    $LocalAppData = [Environment]::GetFolderPath("LocalApplicationData")
    $DedicatedPython = Join-Path $LocalAppData "chronodiag\venv\Scripts\python.exe"
    if (Test-Path -LiteralPath $DedicatedPython -PathType Leaf) {
        $PythonExe = $DedicatedPython
    }
    else {
        $PythonCommand = Get-Command python -ErrorAction SilentlyContinue
        if ($null -eq $PythonCommand) {
            throw "Python was not found."
        }
        $PythonExe = $PythonCommand.Source
    }
}
$PythonExe = (Resolve-Path -LiteralPath $PythonExe).Path

$Rscript = Get-Command Rscript -ErrorAction SilentlyContinue
if ($null -eq $Rscript) {
    throw "Rscript was not found on PATH."
}

Write-Host "Repository:   $RepoRoot"
Write-Host "Project root: $ProjectRoot"
Write-Host "Python:       $PythonExe"
Write-Host "Rscript:      $($Rscript.Source)"
Write-Host ""

$OldPythonPath = $env:PYTHONPATH
$env:PYTHONPATH = Join-Path $RepoRoot "python\src"
try {
    Write-Host "[1/6] Audit value-blind complete cohorts"
    Invoke-Native `
        -FilePath $PythonExe `
        -Arguments @(
            (Join-Path $RepoRoot "scripts\audit_itrdb_empirical_cohorts_v01.py"),
            "--project-root", $ProjectRoot,
            "--repo-root", $RepoRoot
        ) `
        -FailureMessage "ITRDB empirical cohort audit failed."

    Write-Host ""
    Write-Host "[2/6] Detrend cores and freeze tree-level data"
    Invoke-Native `
        -FilePath $Rscript.Source `
        -Arguments @(
            (Join-Path $RepoRoot "scripts\prepare_itrdb_empirical_stage_b_v01.R"),
            "--project-root", $ProjectRoot,
            "--repo-root", $RepoRoot
        ) `
        -FailureMessage "ITRDB empirical preprocessing failed."

    Write-Host ""
    Write-Host "[3/6] Generate explicit empirical memberships"
    Invoke-Native `
        -FilePath $PythonExe `
        -Arguments @(
            (Join-Path $RepoRoot "scripts\generate_itrdb_empirical_memberships_v01.py"),
            "--repo-root", $RepoRoot
        ) `
        -FailureMessage "ITRDB empirical membership generation failed."

    Write-Host ""
    Write-Host "[4/6] Run composition-aware empirical diagnostics"
    Invoke-Native `
        -FilePath $PythonExe `
        -Arguments @(
            (Join-Path $RepoRoot "scripts\run_itrdb_empirical_stage_b_v01.py"),
            "--repo-root", $RepoRoot
        ) `
        -FailureMessage "ITRDB empirical Stage B analysis failed."

    Write-Host ""
    Write-Host "[5/6] Validate deterministic empirical outputs"
    Invoke-Native `
        -FilePath $PythonExe `
        -Arguments @(
            (Join-Path $RepoRoot "scripts\validate_itrdb_empirical_stage_b_v01.py"),
            "--repo-root", $RepoRoot
        ) `
        -FailureMessage "ITRDB empirical Stage B validation failed."

    Write-Host ""
    Write-Host "[6/6] Build empirical figures"
    $CanPlot = $false
    if (-not $SkipPlots) {
        $oldPreference = $ErrorActionPreference
        try {
            $ErrorActionPreference = "Continue"
            & $PythonExe -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('matplotlib') else 2)" *> $null
            $CanPlot = ($LASTEXITCODE -eq 0)
        }
        finally {
            $ErrorActionPreference = $oldPreference
        }
    }
    if ($SkipPlots) {
        Write-Warning "Empirical figure generation was skipped by request."
    }
    elseif (-not $CanPlot) {
        Write-Warning "matplotlib is not installed; empirical tables were completed but figure generation was skipped. Install the Python 'figures' extra and rerun the plot script."
    }
    else {
        Invoke-Native `
            -FilePath $PythonExe `
            -Arguments @(
                (Join-Path $RepoRoot "scripts\plot_itrdb_empirical_stage_b_v01.py"),
                "--repo-root", $RepoRoot
            ) `
            -FailureMessage "ITRDB empirical figure generation failed."
    }
}
finally {
    $env:PYTHONPATH = $OldPythonPath
}

if (-not $SkipFullValidation) {
    Write-Host ""
    Write-Host "Running the complete repository validation..."
    $ValidationArguments = @(
        "-NoProfile",
        "-ExecutionPolicy", "Bypass",
        "-File", (Join-Path $RepoRoot "scripts\validate_local.ps1"),
        "-ProjectRoot", $ProjectRoot,
        "-PythonExe", $PythonExe
    )
    Invoke-Native `
        -FilePath "powershell.exe" `
        -Arguments $ValidationArguments `
        -FailureMessage "The empirical workflow completed, but full repository validation failed."
}

Write-Host ""
Write-Host "ITRDB empirical Stage B workflow completed successfully."
Write-Host "Results:"
Write-Host "  $(Join-Path $RepoRoot 'analysis\itrdb_empirical_stage_b_v01\results\empirical_key_results_v01.csv')"
Write-Host "Figures:"
Write-Host "  $(Join-Path $RepoRoot 'analysis\itrdb_empirical_stage_b_v01\results\figures')"
