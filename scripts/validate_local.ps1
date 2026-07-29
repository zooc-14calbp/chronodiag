[CmdletBinding()]
param(
    [string]$ProjectRoot,
    [string]$PythonExe,
    [switch]$SkipRCheck
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = (Resolve-Path (Join-Path $RepoRoot "..\..")).Path
} else {
    $ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
}

if ([string]::IsNullOrWhiteSpace($PythonExe)) {
    $LocalAppData = [Environment]::GetFolderPath("LocalApplicationData")
    $DedicatedPython = if ([string]::IsNullOrWhiteSpace($LocalAppData)) {
        $null
    } else {
        Join-Path $LocalAppData "chronodiag\venv\Scripts\python.exe"
    }

    if (-not [string]::IsNullOrWhiteSpace($env:CHRONODIAG_PYTHON) -and
        (Test-Path -LiteralPath $env:CHRONODIAG_PYTHON -PathType Leaf)) {
        $PythonExe = $env:CHRONODIAG_PYTHON
    }
    elseif (-not [string]::IsNullOrWhiteSpace($DedicatedPython) -and
            (Test-Path -LiteralPath $DedicatedPython -PathType Leaf)) {
        $PythonExe = $DedicatedPython
    }
    else {
        $PythonCommand = Get-Command python -ErrorAction SilentlyContinue
        if ($null -eq $PythonCommand) {
            throw "Python was not found. Run setup_chronodiag_python_env_v04.ps1."
        }
        $PythonExe = $PythonCommand.Source
    }
}

$PythonExe = (Resolve-Path -LiteralPath $PythonExe).Path

Write-Host "Repository:   $RepoRoot"
Write-Host "Project root: $ProjectRoot"
Write-Host "Python:       $PythonExe"
Write-Host ""

$OldPythonPath = $env:PYTHONPATH
$env:PYTHONPATH = Join-Path $RepoRoot "python\src"
try {
    Write-Host "[1/3] Python dependency, unit, and parity tests"

    & $PythonExe -c "import importlib.util, sys; required=('numpy','pandas','pytest'); missing=[m for m in required if importlib.util.find_spec(m) is None]; print('Missing Python modules: ' + ', '.join(missing)) if missing else None; sys.exit(2 if missing else 0)"
    if ($LASTEXITCODE -ne 0) {
        throw "Required Python dependencies are missing. Run setup_chronodiag_python_env_v04.ps1."
    }

    & $PythonExe -m pytest (Join-Path $RepoRoot "python\tests") -q
    if ($LASTEXITCODE -ne 0) { throw "Python tests failed." }

    Write-Host ""
    Write-Host "[2/3] Hash and full ITRDB audit validation"
    & $PythonExe (Join-Path $RepoRoot "scripts\validate_pre_github.py") `
        --project-root $ProjectRoot
    if ($LASTEXITCODE -ne 0) { throw "ITRDB validation failed." }
}
finally {
    $env:PYTHONPATH = $OldPythonPath
}

Write-Host ""
Write-Host "[3/3] R validation"
if ($SkipRCheck) {
    Write-Warning "R validation was skipped by request."
} else {
    $Rscript = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($null -eq $Rscript) {
        Write-Warning "Rscript was not found. R validation remains for GitHub Actions."
    } else {
        Push-Location $RepoRoot
        try {
            & $Rscript.Source "r/parity_tests/test_shared_metrics.R"
            if ($LASTEXITCODE -ne 0) { throw "R shared-metric parity test failed." }

            & $Rscript.Source "r/parity_tests/test_shared_composition.R"
            if ($LASTEXITCODE -ne 0) { throw "R shared-composition parity test failed." }

            & $Rscript.Source -e "quit(status=if(requireNamespace('dplR', quietly=TRUE)) 0 else 2)"
            if ($LASTEXITCODE -eq 0) {
                & $Rscript.Source "r/parity_tests/test_shared_io.R"
                if ($LASTEXITCODE -ne 0) { throw "R shared-I/O parity test failed." }
            } else {
                Write-Warning "dplR is not installed; the RWL parity test was deferred to CI."
            }

            & $Rscript.Source -e "quit(status=if(all(vapply(c('roxygen2','rcmdcheck','testthat'), requireNamespace, logical(1), quietly=TRUE))) 0 else 2)"
            if ($LASTEXITCODE -eq 0) {
                & $Rscript.Source -e "roxygen2::roxygenise('r/chronodiag')"
                if ($LASTEXITCODE -ne 0) { throw "R documentation generation failed." }
                & $Rscript.Source -e "rcmdcheck::rcmdcheck('r/chronodiag', args='--no-manual', build_args='--no-manual', error_on='warning', check_dir='validation_output/rcheck')"
                if ($LASTEXITCODE -ne 0) { throw "R CMD check failed." }
            } else {
                Write-Warning "roxygen2, rcmdcheck, or testthat is missing; full R CMD check was deferred to CI."
            }
        } finally {
            Pop-Location
        }
    }
}

Write-Host ""
Write-Host "Local pre-GitHub validation completed."
