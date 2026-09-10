New-Item -ItemType Directory -Force -Path "results" | Out-Null

$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logFile = "results\pytest_$timestamp.log"

python -m pytest -v 2>&1 | Tee-Object -FilePath $logFile

if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}