[CmdletBinding()]
param(
    [string]$ProjectRoot,
    [string]$GitHubOwner = "zooc-14calbp",
    [string]$RepositoryName = "chronodiag",
    [string]$RemoteUrl,
    [switch]$SkipValidation,
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

if (-not $SkipValidation) {
    $validateArgs = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", (Join-Path $PSScriptRoot "validate_local.ps1"), "-ProjectRoot", $ProjectRoot)
    if ($SkipRCheck) { $validateArgs += "-SkipRCheck" }
    & powershell @validateArgs
    if ($LASTEXITCODE -ne 0) { throw "Validation failed; Git initialization was not performed." }
}

$Git = Get-Command git -ErrorAction SilentlyContinue
if ($null -eq $Git) { throw "git was not found on PATH." }

Push-Location $RepoRoot
try {
    if (-not (Test-Path -LiteralPath ".git")) {
        & $Git.Source init -b main
        if ($LASTEXITCODE -ne 0) { throw "git init failed." }
    }

    $name = (& $Git.Source config --get user.name 2>$null)
    if ([string]::IsNullOrWhiteSpace($name)) {
        & $Git.Source config user.name "Hitoshi Yonenobu"
    }
    $email = (& $Git.Source config --get user.email 2>$null)
    if ([string]::IsNullOrWhiteSpace($email)) {
        & $Git.Source config user.email "yn@naruto-u.ac.jp"
    }

    & $Git.Source add --all
    if ($LASTEXITCODE -ne 0) { throw "git add failed." }

    & $Git.Source diff --cached --quiet
    if ($LASTEXITCODE -ne 0) {
        & $Git.Source commit -m "Initialize chronodiag dual-language validation baseline"
        if ($LASTEXITCODE -ne 0) { throw "git commit failed." }
    } else {
        Write-Host "No uncommitted changes were found."
    }

    $Gh = Get-Command gh -ErrorAction SilentlyContinue
    if ($null -ne $Gh -and [string]::IsNullOrWhiteSpace($RemoteUrl)) {
        & $Gh.Source auth status
        if ($LASTEXITCODE -eq 0) {
            & $Gh.Source repo create "$GitHubOwner/$RepositoryName" --private --source . --remote origin --push
            if ($LASTEXITCODE -ne 0) { throw "gh repo create/push failed." }
            Write-Host "Private GitHub repository created and pushed: $GitHubOwner/$RepositoryName"
            exit 0
        }
    }

    if (-not [string]::IsNullOrWhiteSpace($RemoteUrl)) {
        $origin = (& $Git.Source remote get-url origin 2>$null)
        if ([string]::IsNullOrWhiteSpace($origin)) {
            & $Git.Source remote add origin $RemoteUrl
        } elseif ($origin -ne $RemoteUrl) {
            throw "An origin remote already exists and differs from -RemoteUrl: $origin"
        }
        & $Git.Source push -u origin main
        if ($LASTEXITCODE -ne 0) { throw "git push failed." }
        Write-Host "Pushed to $RemoteUrl"
        exit 0
    }

    Write-Host ""
    Write-Host "Local Git repository and initial commit are ready."
    Write-Host "No authenticated GitHub CLI or -RemoteUrl was available, so no remote write was attempted."
    Write-Host "Create an empty private repository named '$RepositoryName', then rerun with:"
    Write-Host "  -RemoteUrl <the repository clone URL>"
} finally {
    Pop-Location
}
