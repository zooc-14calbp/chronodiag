[CmdletBinding()]
param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PythonExe = "python",
    [string]$RscriptExe = "Rscript"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot = (Resolve-Path -LiteralPath $RepoRoot).Path
$PythonExe = (Get-Command $PythonExe -ErrorAction Stop).Source
$RscriptExe = (Get-Command $RscriptExe -ErrorAction Stop).Source
$AnalysisRoot = Join-Path $RepoRoot "analysis\yoshino_empirical_stage_c_v01"
$RawPath = Join-Path $RepoRoot "data\yoshino\raw\nr_all.rwl"

function Invoke-Native {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [Parameter(Mandatory = $false)][string[]]$Arguments = @(),
        [Parameter(Mandatory = $true)][string]$FailureMessage
    )
    $old = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & $FilePath @Arguments
        $code = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $old
    }
    if ($code -ne 0) {
        throw "$FailureMessage Exit code: $code"
    }
}

Write-Host "Repository: $RepoRoot"
Write-Host "Python:     $PythonExe"
Write-Host "Rscript:    $RscriptExe"
Write-Host ""

Write-Host "[1/6] Audit Yoshino specimen-mean source and value-blind cohort"
Invoke-Native -FilePath $PythonExe -Arguments @(
    (Join-Path $RepoRoot "scripts\audit_yoshino_empirical_stage_c_v01.py"),
    "--repo-root", $RepoRoot,
    "--raw-path", $RawPath,
    "--output-dir", $AnalysisRoot
) -FailureMessage "Yoshino source/cohort audit failed."

Write-Host ""
Write-Host "[2/6] Detrend frozen specimen-mean cohort"
Invoke-Native -FilePath $RscriptExe -Arguments @(
    (Join-Path $RepoRoot "scripts\prepare_yoshino_empirical_stage_c_v01.R"),
    "--repo-root", $RepoRoot,
    "--raw-path", $RawPath,
    "--output-dir", (Join-Path $AnalysisRoot "preprocessed")
) -FailureMessage "Yoshino preprocessing failed."

Write-Host ""
Write-Host "[3/6] Generate explicit Yoshino memberships"
Invoke-Native -FilePath $PythonExe -Arguments @(
    (Join-Path $RepoRoot "scripts\generate_yoshino_empirical_memberships_v01.py"),
    "--repo-root", $RepoRoot,
    "--preprocessed-dir", (Join-Path $AnalysisRoot "preprocessed"),
    "--output-dir", (Join-Path $AnalysisRoot "memberships")
) -FailureMessage "Yoshino membership generation failed."

Write-Host ""
Write-Host "[4/6] Run Yoshino composition-aware diagnostics"
Invoke-Native -FilePath $PythonExe -Arguments @(
    (Join-Path $RepoRoot "scripts\run_yoshino_empirical_stage_c_v01.py"),
    "--repo-root", $RepoRoot,
    "--preprocessed-dir", (Join-Path $AnalysisRoot "preprocessed"),
    "--membership-dir", (Join-Path $AnalysisRoot "memberships"),
    "--output-dir", (Join-Path $AnalysisRoot "results")
) -FailureMessage "Yoshino diagnostics failed."

Write-Host ""
Write-Host "[5/6] Validate deterministic Yoshino outputs"
Invoke-Native -FilePath $PythonExe -Arguments @(
    (Join-Path $RepoRoot "scripts\validate_yoshino_empirical_stage_c_v01.py"),
    "--repo-root", $RepoRoot
) -FailureMessage "Yoshino deterministic validation failed."

Write-Host ""
Write-Host "[6/6] Build Yoshino figures"
$old = $ErrorActionPreference
try {
    $ErrorActionPreference = "Continue"
    & $PythonExe -c "import matplotlib" *> $null
    $hasMatplotlib = ($LASTEXITCODE -eq 0)
}
finally {
    $ErrorActionPreference = $old
}
if ($hasMatplotlib) {
    Invoke-Native -FilePath $PythonExe -Arguments @(
        (Join-Path $RepoRoot "scripts\plot_yoshino_empirical_stage_c_v01.py"),
        "--repo-root", $RepoRoot
    ) -FailureMessage "Yoshino figure generation failed."
}
else {
    Write-Warning "matplotlib is absent; Yoshino tables completed but figures were skipped."
}

Write-Host ""
Write-Host "Running complete repository validation..."
Invoke-Native -FilePath "powershell.exe" -Arguments @(
    "-NoProfile", "-ExecutionPolicy", "Bypass",
    "-File", (Join-Path $RepoRoot "scripts\validate_local.ps1"),
    "-SkipYoshinoStageC"
) -FailureMessage "Repository validation failed after Yoshino Stage C."

Write-Host ""
Write-Host "Yoshino empirical Stage C workflow completed successfully."
Write-Host "Key results:"
Write-Host "  $(Join-Path $AnalysisRoot 'results\yoshino_key_results_v01.csv')"
